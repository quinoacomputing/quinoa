//******************************************************************************
/*!
  \file      src/Walker/Distributor.C
  \author    J. Bakosi
  \date      Mon 20 Jul 2015 07:48:36 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Distributor drives the time integration of differential equations
  \details   Distributor drives the time integration of differential equations.
    The implementation uses the Charm++ runtime system and is fully asynchronous,
    overlapping computation, communication as well as I/O. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Walker/distributor.ci.
*/
//******************************************************************************

#include <list>
#include <string>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <cmath>
#include <cstddef>
#include <sys/sysmacros.h>

#include <boost/format.hpp>
#include <boost/optional.hpp>

#include "Print.h"
#include "Tags.h"
#include "StatCtr.h"
#include "Exception.h"
#include "ParticleProperties.h"
#include "LoadDistributor.h"
#include "Distributor.h"
#include "Integrator.h"
#include "DiffEqStack.h"
#include "TxtStatWriter.h"
#include "PDFUtil.h"
#include "PDFWriter.h"
#include "Options/PDFFile.h"
#include "Options/PDFPolicy.h"
#include "Walker/InputDeck/InputDeck.h"
#include "walker.decl.h"

extern CProxy_Main mainProxy;

using walker::Distributor;

Distributor::Distributor( const ctr::CmdLine& cmdline ) :
  m_print( cmdline.get< tag::verbose >() ? std::cout : std::clog ),
  m_output( false, false ),
  m_it( 0 ),
  m_nchare( 0 ),
  m_t( 0.0 ),
  m_nameOrdinary( g_inputdeck.momentNames( tk::ctr::ordinary ) ),
  m_nameCentral( g_inputdeck.momentNames( tk::ctr::central ) ),
  m_ordinary( m_nameOrdinary.size(), 0.0 ),
  m_central( m_nameCentral.size(), 0.0 )
//******************************************************************************
// Constructor
//! \param[in] cmdline Data structure storing data from the command-line parser
//! \author J. Bakosi
//******************************************************************************
{
  // Compute load distribution given total work (= number of particles) and
  // user-specified virtualization
  uint64_t chunksize, remainder;
  m_nchare = tk::linearLoadDistributor(
               g_inputdeck.get< tag::cmd, tag::virtualization >(),
               g_inputdeck.get< tag::discr, tag::npar >(),
               CkNumPes(),
               chunksize,
               remainder );

  // Compute total number of particles distributed over all workers
  // Note that this number will not necessarily be the same as given by the
  // user, coming from g_inputdeck.get< tag::discr, tag::npar >(), since each
  // Charm++ chare array element constructor takes this chunksize argument,
  // which equals the number of particles the array element (worker) will work
  // on.
  m_npar = static_cast< tk::real >( m_nchare * chunksize );

  // Print out info on what will be done and how
  info( chunksize );

  // Start timer measuring total integration time
  m_timer.emplace_back();

  // Compute size of initial time step
  m_dt = computedt();

  // Construct and initialize map of statistical moments
  for (const auto& product : g_inputdeck.get< tag::stat >())
    m_moments[ product ] = 0.0;

  // Activate SDAG-wait for estimation of ordinary statistics
  wait4ord();
  // Activate SDAG-wait for estimation of central moments
  wait4cen();
  // Activate SDAG-wait for estimation of PDFs at select times
  if ( !(m_it % g_inputdeck.get< tag::interval, tag::pdf >()) ) wait4pdf();

  // Create statistics merger chare group collecting chare contributions
  CProxy_Collector collproxy = CProxy_Collector::ckNew( thisProxy );

  // Fire up asynchronous differential equation integrators
  m_intproxy =
    CProxy_Integrator::ckNew( thisProxy, collproxy, chunksize,
                              static_cast<int>( m_nchare ) );
}

void
Distributor::info( uint64_t chunksize ) const
//******************************************************************************
//  Print information at startup
//! \param[in] chunksize Chunk size, see Base/LoadDistribution.h
//! \author J. Bakosi
//******************************************************************************
{
  // Print out info data layout
  m_print.list( "Particle properties data layout policy (CMake: LAYOUT)",
                std::list< std::string >{ tk::ParProps().major() } );

  // Re-create differential equations stack for output
  DiffEqStack stack;

  // Print out information on factory
  m_print.eqlist( "Registered differential equations", stack.factory(),
                  stack.ntypes() );
  m_print.endpart();

  // Print out information on problem
  m_print.part( "Problem" );

  // Print out info on problem title
  if ( !g_inputdeck.get< tag::title >().empty() )
    m_print.title( g_inputdeck.get< tag::title >() );

  // Print out info on settings of selected differential equations
  m_print.diffeqs( "Differential equations integrated", stack.info() );

  // Print out info on RNGs selected
  // ...

  // Print I/O filenames
  m_print.section( "Output filenames" );
  if (!g_inputdeck.get< tag::stat >().empty())
    m_print.item( "Statistics", g_inputdeck.get< tag::cmd, tag::io, tag::stat >() );
  if (!g_inputdeck.get< tag::pdf >().empty())
    m_print.item( "PDF", g_inputdeck.get< tag::cmd, tag::io, tag::pdf >() );

  // Print discretization parameters
  m_print.section( "Discretization parameters" );
  m_print.item( "Number of time steps",
                g_inputdeck.get< tag::discr, tag::nstep >() );
  m_print.item( "Terminate time",
                g_inputdeck.get< tag::discr, tag::term >() );
  m_print.item( "Initial time step size",
                g_inputdeck.get< tag::discr, tag::dt >() );

  // Print output intervals
  m_print.section( "Output intervals" );
  m_print.item( "TTY", g_inputdeck.get< tag::interval, tag::tty>() );
  if (!g_inputdeck.get< tag::stat >().empty())
    m_print.item( "Statistics", g_inputdeck.get< tag::interval, tag::stat >() );
  if (!g_inputdeck.get< tag::pdf >().empty())
    m_print.item( "PDF", g_inputdeck.get< tag::interval, tag::pdf >() );

  // Print out statistics estimated
  m_print.statistics( "Statistical moments and distributions" );

  // Print out info on load distirubtion
  m_print.section( "Load distribution" );
  m_print.item( "Virtualization [0.0...1.0]",
                g_inputdeck.get< tag::cmd, tag::virtualization >() );
  m_print.item( "Number of processing elements", CkNumPes() );
  m_print.item( "Number of work units", m_nchare );
  m_print.item( "User load (# of particles)",
                g_inputdeck.get< tag::discr, tag::npar >() );
  m_print.item( "Chunksize (load per work unit)", chunksize );
  m_print.item( "Actual load (# of particles)",
                std::to_string( m_nchare * chunksize ) +
                " (=" +
                std::to_string( m_nchare ) + "*" +
                std::to_string( chunksize ) + ")" );

  // Print out time integration header
  if (g_inputdeck.get< tag::discr, tag::nstep >()) {
    header();
    tk::TxtStatWriter sw( !m_nameOrdinary.empty() || !m_nameCentral.empty() ?
                          g_inputdeck.get< tag::cmd, tag::io, tag::stat >() :
                          std::string(),
                          g_inputdeck.get< tag::flformat, tag::stat >(),
                          g_inputdeck.get< tag::prec, tag::stat >() );
    sw.header( m_nameOrdinary, m_nameCentral );
  }
}

tk::real
Distributor::computedt()
//******************************************************************************
// Compute size of next time step
//! \return Size of dt for the next time step
//! \author  J. Bakosi
//******************************************************************************
{
  // Simply return a constant user-defined dt for now
  return g_inputdeck.get< tag::discr, tag::dt >();
}

void
Distributor::init() const
//******************************************************************************
// Reduction target indicating that all integrators finished initialization
//! \details Upon calling this Charm++ reduction target, we simply put in a time
//!   stamp measuring setting the initial conditions.
//! \author J. Bakosi
//******************************************************************************
{
  mainProxy.timestamp( "Initial conditions", m_timer[0].dsec() );
}

void
Distributor::estimateOrd( tk::real* ord, std::size_t n )
//******************************************************************************
// Estimate ordinary moments
//! \param[in] ord Ordinary moments (sum) collected over all chares
//! \param[in] n Number of ordinary moments in array ord
//! \author J. Bakosi
//******************************************************************************
{
  Assert( n == m_ordinary.size(),
          "Number of ordinary moments contributed not equal to expected" );

  // Add contribution from PE to total sums, i.e., u[i] += v[i] for all i
  for (std::size_t i=0; i<m_ordinary.size(); ++i) m_ordinary[i] += ord[i];

  // Finish computing moments, i.e., divide sums by the number of samples
  for (auto& m : m_ordinary) m /= m_npar;

  // Activate SDAG trigger signaling that ordinary moments have been estimated
  estimateOrdDone();
}

void
Distributor::estimateCen( tk::real* cen, std::size_t n )
//******************************************************************************
// Estimate ordinary moments
//! \param[in] cen Central moments (sum) collected over all chares
//! \param[in] n Number of central moments in array cen
//! \author J. Bakosi
//******************************************************************************
{
  Assert( n == m_central.size(),
          "Number of central moments contributed not equal to expected" );

  // Add contribution from PE to total sums, i.e., u[i] += v[i] for all i
  for (std::size_t i=0; i<m_central.size(); ++i) m_central[i] += cen[i];

  // Finish computing moments, i.e., divide sums by the number of samples
  for (auto& m : m_central) m /= m_npar;

  // Activate SDAG trigger signaling that central moments have been estimated
  estimateCenDone();
}

void
Distributor::estimateOrdPDF( CkReductionMsg* msg )
//******************************************************************************
// Estimate ordinary PDFs
//! \param[in] Serialized tuple of vectors of uni-, bi-, and tri-variate PDFs
//! \author J. Bakosi
//******************************************************************************
{
  // Deserialize and merge PDFs
  std::tie( m_ordupdf, m_ordbpdf, m_ordtpdf ) = tk::merge( msg );

  // Activate SDAG trigger signaling that ordinary PDFs have been estimated
  estimateOrdPDFDone();
}

void
Distributor::estimateCenPDF( CkReductionMsg* msg )
//******************************************************************************
// Estimate central PDFs
//! \param[in] Serialized tuple of vectors of uni-, bi-, and tri-variate PDFs
//! \author J. Bakosi
//******************************************************************************
{
  // Deserialize and merge PDFs
  std::tie( m_cenupdf, m_cenbpdf, m_centpdf ) = tk::merge( msg );

  // Activate SDAG trigger signaling that central PDFs have been estimated
  estimateCenPDFDone();
}

void
Distributor::outStat()
//******************************************************************************
// Output statistics to file
//! \author J. Bakosi
//******************************************************************************
{
  // Append statistics file at selected times
  if (!(m_it % g_inputdeck.get< tag::interval, tag::stat >())) {
    tk::TxtStatWriter sw( !m_nameOrdinary.empty() || !m_nameCentral.empty() ?
                          g_inputdeck.get< tag::cmd, tag::io, tag::stat >() :
                          std::string(),
                          g_inputdeck.get< tag::flformat, tag::stat >(),
                          g_inputdeck.get< tag::prec, tag::stat >(),
                          std::ios_base::app );
    if (sw.stat( m_it, m_t, m_ordinary, m_central ))
      m_output.get< tag::stat >() = true;
  }
}

void
Distributor::outPDF()
//******************************************************************************
// Output PDFs to file
//! \author J. Bakosi
//******************************************************************************
{
  // Output PDFs at selected times
  if ( !(m_it % g_inputdeck.get< tag::interval, tag::pdf >()) ) {
    auto n = outUniPDF();              // Output univariate PDFs to file(s)
    n += outBiPDF();                   // Output bivariate PDFs to file(s)
    n += outTriPDF();                  // Output trivariate PDFs to file(s)
    if (n) m_output.get< tag::pdf >() = true; // Signal that PDFs were written
  }
}

void
Distributor::writeUniPDF( const tk::UniPDF& p, int& cnt )
//******************************************************************************
// Write univariate PDF to file
//! \param[in] p Univariate PDF to output
//! \param[inout] cnt Count up number of PDFs written
//! \author J. Bakosi
//******************************************************************************
{
  // Get PDF metadata
  const auto info =
    tk::ctr::pdfInfo< 1 >( g_inputdeck.get< tag::discr, tag::binsize >(),
                           g_inputdeck.get< tag::cmd, tag::io, tag::pdfnames >(),
                           g_inputdeck.get< tag::discr, tag::extent >(),
                           g_inputdeck.get< tag::pdf >(),
                           cnt++ );

  // Construct PDF file name: base name + '_' + pdf name
  std::string filename =
    g_inputdeck.get< tag::cmd, tag::io, tag::pdf >() + '_' + info.name;

  // Augment PDF filename by time stamp if PDF output file policy is multiple
  if (g_inputdeck.get< tag::selected, tag::pdfpolicy >() ==
      tk::ctr::PDFPolicyType::MULTIPLE)
    filename += '_' + std::to_string( m_t );

  // Augment PDF filename by '.txt' extension
  filename += ".txt";

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfw( filename,
                      g_inputdeck.get< tag::flformat, tag::pdf >(),
                      g_inputdeck.get< tag::prec, tag::pdf >() );

  // Output PDF
  pdfw.writeTxt( p, info );
}

void
Distributor::writeBiPDF( const tk::BiPDF& p, int& cnt )
//******************************************************************************
// Write bivariate PDF to file
//! \param[in] p Bivariate PDF to output
//! \param[inout] cnt Count up number of PDFs written
//! \author J. Bakosi
//******************************************************************************
{
  // Get PDF metadata
  const auto info =
    tk::ctr::pdfInfo< 2 >( g_inputdeck.get< tag::discr, tag::binsize >(),
                           g_inputdeck.get< tag::cmd, tag::io, tag::pdfnames >(),
                           g_inputdeck.get< tag::discr, tag::extent >(),
                           g_inputdeck.get< tag::pdf >(),
                           cnt++ );

  // Construct PDF file name: base name + '_' + pdf name
  std::string filename =
    g_inputdeck.get< tag::cmd, tag::io, tag::pdf >() + '_' + info.name;

  // Augment PDF filename by time stamp if PDF output file policy is multiple
  if (g_inputdeck.get< tag::selected, tag::pdfpolicy >() ==
      tk::ctr::PDFPolicyType::MULTIPLE)
    filename += '_' + std::to_string( m_t );

  const auto& filetype = g_inputdeck.get< tag::selected, tag::pdffiletype >();

  // Augment PDF filename by the appropriate extension
  if (filetype == tk::ctr::PDFFileType::TXT)
    filename += ".txt";
  else if (filetype == tk::ctr::PDFFileType::GMSHTXT ||
           filetype == tk::ctr::PDFFileType::GMSHBIN )
    filename += ".gmsh";
  else if (filetype == tk::ctr::PDFFileType::EXODUSII)
    filename += ".exo";
  else Throw( "Unkown PDF file type attempting to output bivariate PDF" );

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfw( filename,
                      g_inputdeck.get< tag::flformat, tag::pdf >(),
                      g_inputdeck.get< tag::prec, tag::pdf >() );

  // Output PDF
  if (filetype == tk::ctr::PDFFileType::TXT)
    pdfw.writeTxt( p, info );
  else if (filetype == tk::ctr::PDFFileType::GMSHTXT)
    pdfw.writeGmshTxt( p, info,
                       g_inputdeck.get< tag::selected, tag::pdfctr >() );
  else if (filetype == tk::ctr::PDFFileType::GMSHBIN)
    pdfw.writeGmshBin( p, info,
                       g_inputdeck.get< tag::selected, tag::pdfctr >() );
  else if (filetype == tk::ctr::PDFFileType::EXODUSII)
    pdfw.writeExodusII( p, info,
                        m_it,
                        g_inputdeck.get< tag::selected, tag::pdfctr >() );
}

void
Distributor::writeTriPDF( const tk::TriPDF& p, int& cnt )
//******************************************************************************
// Write trivariate PDF to file
//! \param[in] p Trivariate PDF to output
//! \param[inout] cnt Count up number of PDFs written
//! \author J. Bakosi
//******************************************************************************
{
  // Get PDF metadata
  const auto info =
    tk::ctr::pdfInfo< 3 >( g_inputdeck.get< tag::discr, tag::binsize >(),
                           g_inputdeck.get< tag::cmd, tag::io, tag::pdfnames >(),
                           g_inputdeck.get< tag::discr, tag::extent >(),
                           g_inputdeck.get< tag::pdf >(),
                           cnt++ );

  // Construct PDF file name: base name + '_' + pdf name
  std::string filename =
    g_inputdeck.get< tag::cmd, tag::io, tag::pdf >() + '_' + info.name;

  // Augment PDF filename by time stamp if PDF output file policy is multiple
  if (g_inputdeck.get< tag::selected, tag::pdfpolicy >() ==
      tk::ctr::PDFPolicyType::MULTIPLE)
    filename += '_' + std::to_string( m_t );

  const auto& filetype = g_inputdeck.get< tag::selected, tag::pdffiletype >();

  // Augment PDF filename by the appropriate extension
  if (filetype == tk::ctr::PDFFileType::TXT)
    filename += ".txt";
  else if (filetype == tk::ctr::PDFFileType::GMSHTXT ||
           filetype == tk::ctr::PDFFileType::GMSHBIN )
    filename += ".gmsh";
  else if (filetype == tk::ctr::PDFFileType::EXODUSII)
    filename += ".exo";
  else Throw( "Unkown PDF file type attempting to output trivariate PDF" );

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfw( filename,
                      g_inputdeck.get< tag::flformat, tag::pdf >(),
                      g_inputdeck.get< tag::prec, tag::pdf >() );

  // Output PDF
  if (filetype == tk::ctr::PDFFileType::TXT)
    pdfw.writeTxt( p, info );
  else if (filetype == tk::ctr::PDFFileType::GMSHTXT)
     pdfw.writeGmshTxt( p, info,
                        g_inputdeck.get< tag::selected, tag::pdfctr >() );
  else if (filetype == tk::ctr::PDFFileType::GMSHBIN)
     pdfw.writeGmshBin( p, info,
                        g_inputdeck.get< tag::selected, tag::pdfctr >() );
  else if (filetype == tk::ctr::PDFFileType::EXODUSII)
    pdfw.writeExodusII( p, info,
                        m_it,
                        g_inputdeck.get< tag::selected, tag::pdfctr >() );
}

int
Distributor::outUniPDF()
//******************************************************************************
// Output all requested univariate PDFs to file(s)
//! \return Number of PDFs written
//! \author J. Bakosi
//******************************************************************************
{
  int cnt = 0;
  for (const auto& p : m_ordupdf) writeUniPDF( p, cnt );
  for (const auto& p : m_cenupdf) writeUniPDF( p, cnt );
  return cnt;
}

int
Distributor::outBiPDF()
//******************************************************************************
// Output all requested bivariate PDFs to file(s)
//! \return Number of PDFs written
//! \author J. Bakosi
//******************************************************************************
{
  int cnt = 0;
  for (const auto& p : m_ordbpdf) writeBiPDF( p, cnt );
  for (const auto& p : m_cenbpdf) writeBiPDF( p, cnt );
  return cnt;
}

int
Distributor::outTriPDF()
//******************************************************************************
// Output all requested trivariate PDFs to file(s)
//! \return Number of PDFs written
//! \author J. Bakosi
//******************************************************************************
{
  int cnt = 0;
  for (const auto& p : m_ordtpdf) writeTriPDF( p, cnt );
  for (const auto& p : m_centpdf) writeTriPDF( p, cnt );
  return cnt;
}

void
Distributor::evaluateTime()
//******************************************************************************
// Evaluate time step, compute new time step size, decide if it is time to quit
//! \author J. Bakosi
//******************************************************************************
{
  // Increase number of iterations taken
  ++m_it;

  // Compute size of next time step
  m_dt = computedt();

  // Advance physical time
  m_t += m_dt;

  // Get physical time at which to terminate
  const auto term = g_inputdeck.get< tag::discr, tag::term >();

  // Truncate the size of last time step
  if (m_t > term) m_t = term;

  // Echo one-liner info on time step
  report();

  // Finish if either max iterations or max time reached 
  if ( std::fabs(m_t - term) > std::numeric_limits< tk::real >::epsilon() &&
       m_it < g_inputdeck.get< tag::discr, tag::nstep >() ) {

    if (g_inputdeck.stat()) {
      // Update map of statistical moments
      std::size_t ord = 0;
      std::size_t cen = 0;
      for (const auto& product : g_inputdeck.get< tag::stat >())
        if (tk::ctr::ordinary( product ))
          m_moments[ product ] = m_ordinary[ ord++ ];
        else
          m_moments[ product ] = m_central[ cen++ ];

      // Zero statistics counters and accumulators
      std::fill( begin(m_ordinary), end(m_ordinary), 0.0 );
      std::fill( begin(m_central), end(m_central), 0.0 );

      // Re-activate SDAG-wait for estimation of ordinary statistics for next step
      wait4ord();
      // Re-activate SDAG-wait for estimation of central moments for next step
      wait4cen();

      // Selectively re-activate SDAG-wait for estimation of PDFs for next step
      if ( !(m_it % g_inputdeck.get< tag::interval, tag::pdf >()) ) {
        // Zero PDF counters and accumulators
        for (auto& p : m_ordupdf) p.zero();
        for (auto& p : m_ordbpdf) p.zero();
        for (auto& p : m_ordtpdf) p.zero();
        for (auto& p : m_cenupdf) p.zero();
        for (auto& p : m_cenbpdf) p.zero();
        for (auto& p : m_centpdf) p.zero();
        wait4pdf();
      }
    }

    // Continue with next time step with all integrators
    m_intproxy.advance( m_dt, m_it, m_moments );

  } else {

    // Normal finish, print out reason
    const auto term = g_inputdeck.get< tag::discr, tag::term >();
    const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
    m_print.endsubsection();
    if (m_it >= g_inputdeck.get< tag::discr, tag::nstep >())
       m_print.note( "Normal finish, maximum number of iterations reached: " +
                     std::to_string( nstep ) );
     else 
       m_print.note( "Normal finish, maximum time reached: " +
                     std::to_string( term ) );
    // Quit
    mainProxy.finalize();

  }
}

void
Distributor::nostat()
//******************************************************************************
//  Charm++ reduction target enabling shortcutting sync points if no stats
//! \details This reduction target is called if there are no statistics nor PDFs
//!   to be estimated and thus some synchronization points can be skipped. Upon
//!   this call we simply finish up the time step as usual.
//! \author J. Bakosi
//******************************************************************************
{
  evaluateTime();
}

void
Distributor::header() const
//******************************************************************************
// Print out time integration header
//! \author J. Bakosi
//******************************************************************************
{
  m_print.inthead( "Time integration", "Differential equations testbed",
    "Legend: it - iteration count\n"
    "         t - time\n"
    "        dt - time step size\n"
    "       ETE - estimated time elapsed (h:m:s)\n"
    "       ETA - estimated time for accomplishment (h:m:s)\n"
    "       out - output-saved flags (S: statistics, P: PDFs)\n",
    "\n      it             t            dt        ETE        ETA   out\n"
      " ---------------------------------------------------------------\n" );
}

void
Distributor::report()
//******************************************************************************
// Print out one-liner report on time step
//! \author J. Bakosi
//******************************************************************************
{
  if (!(m_it % g_inputdeck.get< tag::interval, tag::tty >())) {

    // estimated time elapsed and for accomplishment
    tk::Timer::Watch ete, eta;
    m_timer[0].eta( g_inputdeck.get< tag::discr, tag::term >(), m_t,
                    g_inputdeck.get< tag::discr, tag::nstep >(), m_it,
                    ete, eta );

    // Output one-liner
    m_print << std::setfill(' ') << std::setw(8) << m_it << "  "
            << std::scientific << std::setprecision(6)
            << std::setw(12) << m_t << "  "
            << m_dt << "  "
            << std::setfill('0')
            << std::setw(3) << ete.hrs.count() << ":"
            << std::setw(2) << ete.min.count() << ":"
            << std::setw(2) << ete.sec.count() << "  "
            << std::setw(3) << eta.hrs.count() << ":"
            << std::setw(2) << eta.min.count() << ":"
            << std::setw(2) << eta.sec.count() << "  ";

    // Augment one-liner with output indicators
    if (m_output.get< tag::stat >()) m_print << 'S';
    if (m_output.get< tag::pdf >()) m_print << 'P';

    m_print << '\n';
  }

  // Reset output indicators
  m_output.get< tag::stat >() = false;
  m_output.get< tag::pdf >() = false;
}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "distributor.def.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
