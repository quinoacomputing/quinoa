// *****************************************************************************
/*!
  \file      src/Walker/Distributor.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Distributor drives the time integration of differential equations
  \details   Distributor drives the time integration of differential equations.
    The implementation uses the Charm++ runtime system and is fully asynchronous,
    overlapping computation, communication as well as I/O. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Walker/distributor.ci.
*/
// *****************************************************************************

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

#include "NoWarning/format.h"
#include <boost/optional.hpp>

#include "Macro.h"
#include "Print.h"
#include "Tags.h"
#include "StatCtr.h"
#include "Exception.h"
#include "Particles.h"
#include "LoadDistributor.h"
#include "Distributor.h"
#include "Integrator.h"
#include "DiffEqStack.h"
#include "TxtStatWriter.h"
#include "PDFReducer.h"
#include "PDFWriter.h"
#include "Options/PDFFile.h"
#include "Options/PDFPolicy.h"
#include "Walker/InputDeck/InputDeck.h"
#include "NoWarning/walker.decl.h"

extern CProxy_Main mainProxy;

using walker::Distributor;

Distributor::Distributor( const ctr::CmdLine& cmdline ) :
  __dep(),
  m_print( cmdline.get< tag::verbose >() ? std::cout : std::clog ),
  m_output( false, false ),
  m_it( 0 ),
  m_npar( 0 ),
  m_t( 0.0 ),
  m_dt( computedt() ),
  m_intproxy(),
  m_timer(),
  m_nameOrdinary( g_inputdeck.momentNames( tk::ctr::ordinary ) ),
  m_nameCentral( g_inputdeck.momentNames( tk::ctr::central ) ),
  m_ordinary( m_nameOrdinary.size(), 0.0 ),
  m_central( m_nameCentral.size(), 0.0 ),
  m_ordupdf(),
  m_ordbpdf(),
  m_ordtpdf(),
  m_cenupdf(),
  m_cenbpdf(),
  m_centpdf(),
  m_tables(),
  m_moments()
// *****************************************************************************
// Constructor
//! \param[in] cmdline Data structure storing data from the command-line parser
//! \author J. Bakosi
// *****************************************************************************
{
  // Compute load distribution given total work (= number of particles) and
  // user-specified virtualization
  uint64_t chunksize, remainder;
  auto nchare = tk::linearLoadDistributor(
                  g_inputdeck.get< tag::cmd, tag::virtualization >(),
                  g_inputdeck.get< tag::discr, tag::npar >(),
                  CkNumPes(),
                  chunksize,
                  remainder );

  // Compute total number of particles distributed over all workers. Note that
  // this number will not necessarily be the same as given by the user, coming
  // from g_inputdeck.get< tag::discr, tag::npar >(), since each Charm++ chare
  // array element constructor takes this chunksize argument, which equals the
  // number of particles the array element (worker) will work on.
  m_npar = static_cast< tk::real >( nchare * chunksize );

  // Print out info on what will be done and how
  info( chunksize, nchare );

  // Output header for statistics output file
  tk::TxtStatWriter sw( !m_nameOrdinary.empty() || !m_nameCentral.empty() ?
                        g_inputdeck.get< tag::cmd, tag::io, tag::stat >() :
                        std::string(),
                        g_inputdeck.get< tag::flformat, tag::stat >(),
                        g_inputdeck.get< tag::prec, tag::stat >() );
  sw.header( m_nameOrdinary, m_nameCentral, m_tables.first );

  // Print out time integration header
  m_print.endsubsection();
  m_print.diag( "Starting time stepping ..." );
  header();

  // Start timer measuring total integration time
  m_timer.emplace_back();

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
                              static_cast<int>( nchare ) );
}

void
Distributor::info( uint64_t chunksize, std::size_t nchare )
// *****************************************************************************
//  Print information at startup
//! \param[in] chunksize Chunk size, see Base/LoadDistribution.h
//! \param[in] nchare Total number of Charem++ Integrator chares doing work
//! \author J. Bakosi
// *****************************************************************************
{
  m_print.part( "Factory" );

  // Print out info data layout
  m_print.list( "Particle properties data layout (CMake: PARTICLE_DATA_LAYOUT)",
                std::list< std::string >{ tk::Particles::layout() } );

  // Re-create differential equations stack for output
  DiffEqStack stack;

  // Print out information on factory
  m_print.eqlist( "Registered differential equations",
                  stack.factory(), stack.ntypes() );
  m_print.endpart();

  // Instantiate tables to sample and output to statistics file
  m_tables = stack.tables();

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
  m_print.item( "Number of work units", nchare );
  m_print.item( "User load (# of particles)",
                g_inputdeck.get< tag::discr, tag::npar >() );
  m_print.item( "Chunksize (load per work unit)", chunksize );
  m_print.item( "Actual load (# of particles)",
                std::to_string( nchare * chunksize ) +
                " (=" +
                std::to_string( nchare ) + "*" +
                std::to_string( chunksize ) + ")" );
}

tk::real
Distributor::computedt()
// *****************************************************************************
// Compute size of next time step
//! \return Size of dt for the next time step
//! \author  J. Bakosi
// *****************************************************************************
{
  // Simply return a constant user-defined dt for now
  return g_inputdeck.get< tag::discr, tag::dt >();
}

void
Distributor::estimateOrd( tk::real* ord, int n )
// *****************************************************************************
// Estimate ordinary moments
//! \param[in] ord Ordinary moments (sum) collected over all chares
//! \param[in] n Number of ordinary moments in array ord
//! \author J. Bakosi
// *****************************************************************************
{
  #ifdef NDEBUG
  IGNORE(n);
  #endif

  Assert( static_cast<std::size_t>(n) == m_ordinary.size(),
          "Number of ordinary moments contributed not equal to expected" );

  // Add contribution from PE to total sums, i.e., u[i] += v[i] for all i
  for (std::size_t i=0; i<m_ordinary.size(); ++i) m_ordinary[i] += ord[i];

  // Finish computing moments, i.e., divide sums by the number of samples
  for (auto& m : m_ordinary) m /= m_npar;

  // Activate SDAG trigger signaling that ordinary moments have been estimated
  estimateOrdDone();
}

void
Distributor::estimateCen( tk::real* cen, int n )
// *****************************************************************************
// Estimate ordinary moments
//! \param[in] cen Central moments (sum) collected over all chares
//! \param[in] n Number of central moments in array cen
//! \author J. Bakosi
// *****************************************************************************
{
  #ifdef NDEBUG
  IGNORE(n);
  #endif

  Assert( static_cast<std::size_t>(n) == m_central.size(),
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
// *****************************************************************************
// Estimate ordinary PDFs
//! \param[in] msg Serialized vectors of uni-, bi-, and tri-variate PDFs
//! \author J. Bakosi
// *****************************************************************************
{
  // Deserialize final PDFs
  PUP::fromMem creator( msg->getData() );
  creator | m_ordupdf;
  creator | m_ordbpdf;
  creator | m_ordtpdf;

  delete msg;

  // Activate SDAG trigger signaling that ordinary PDFs have been estimated
  estimateOrdPDFDone();
}

void
Distributor::estimateCenPDF( CkReductionMsg* msg )
// *****************************************************************************
// Estimate central PDFs
//! \param[in] Serialized vectors of uni-, bi-, and tri-variate PDFs
//! \author J. Bakosi
// *****************************************************************************
{
  // Deserialize final PDFs
  PUP::fromMem creator( msg->getData() );
  creator | m_cenupdf;
  creator | m_cenbpdf;
  creator | m_centpdf;

  delete msg;

  // Activate SDAG trigger signaling that central PDFs have been estimated
  estimateCenPDFDone();
}

void
Distributor::outStat()
// *****************************************************************************
// Output statistics to file
//! \author J. Bakosi
// *****************************************************************************
{
  // lambda to sample tables to write to statistics file
  auto extra = [this]() -> std::vector< tk::real > {
    std::vector< tk::real > x;
    for (const auto& t : this->m_tables.second)
      x.push_back( tk::sample(m_t,t) );
    return x;
  };

  // Append statistics file at selected times
  if (!((m_it+1) % g_inputdeck.get< tag::interval, tag::stat >())) {
    tk::TxtStatWriter sw( !m_nameOrdinary.empty() || !m_nameCentral.empty() ?
                          g_inputdeck.get< tag::cmd, tag::io, tag::stat >() :
                          std::string(),
                          g_inputdeck.get< tag::flformat, tag::stat >(),
                          g_inputdeck.get< tag::prec, tag::stat >(),
                          std::ios_base::app );
    if (sw.stat( m_it, m_t, m_ordinary, m_central, extra() ))
      m_output.get< tag::stat >() = true;
  }
}

void
Distributor::outPDF()
// *****************************************************************************
// Output PDFs to file
//! \author J. Bakosi
// *****************************************************************************
{
  // Output PDFs at selected times
  if ( !((m_it+1) % g_inputdeck.get< tag::interval, tag::pdf >()) ) {
    outUniPDF();                        // Output univariate PDFs to file(s)
    outBiPDF();                         // Output bivariate PDFs to file(s)
    outTriPDF();                        // Output trivariate PDFs to file(s)
    m_output.get< tag::pdf >() = true;  // Signal that PDFs were written
  }
}

void
Distributor::writeUniPDF( const tk::UniPDF& p,
                          tk::ctr::Moment m,
                          std::size_t idx )
// *****************************************************************************
// Write univariate PDF to file
//! \param[in] p Univariate PDF to output
//! \param[in] m ORDINARY or CENTRAL PDF we are writing
//! \param[in] idx Index of the PDF of all ordinary or central PDFs requested
//! \author J. Bakosi
// *****************************************************************************
{
  // Get PDF metadata
  const auto nfo =
    tk::ctr::pdfInfo< 1 >( g_inputdeck.get< tag::discr, tag::binsize >(),
                           g_inputdeck.get< tag::cmd, tag::io, tag::pdfnames >(),
                           g_inputdeck.get< tag::discr, tag::extent >(),
                           g_inputdeck.get< tag::pdf >(),
                           m,
                           idx );

  // Construct PDF file name: base name + '_' + pdf name
  std::string filename =
    g_inputdeck.get< tag::cmd, tag::io, tag::pdf >() + '_' + nfo.name;

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
  pdfw.writeTxt( p, nfo );
}

void
Distributor::writeBiPDF( const tk::BiPDF& p,
                         tk::ctr::Moment m,
                         std::size_t idx )
// *****************************************************************************
// Write bivariate PDF to file
//! \param[in] p Bivariate PDF to output
//! \param[in] m ORDINARY or CENTRAL PDF we are writing
//! \param[in] idx Index of the PDF of all ordinary or central PDFs requested
//! \author J. Bakosi
// *****************************************************************************
{
  // Get PDF metadata
  const auto nfo =
    tk::ctr::pdfInfo< 2 >( g_inputdeck.get< tag::discr, tag::binsize >(),
                           g_inputdeck.get< tag::cmd, tag::io, tag::pdfnames >(),
                           g_inputdeck.get< tag::discr, tag::extent >(),
                           g_inputdeck.get< tag::pdf >(),
                           m,
                           idx );

  // Construct PDF file name: base name + '_' + pdf name
  std::string filename =
    g_inputdeck.get< tag::cmd, tag::io, tag::pdf >() + '_' + nfo.name;

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
    pdfw.writeTxt( p, nfo );
  else if (filetype == tk::ctr::PDFFileType::GMSHTXT)
    pdfw.writeGmshTxt( p, nfo,
                       g_inputdeck.get< tag::selected, tag::pdfctr >() );
  else if (filetype == tk::ctr::PDFFileType::GMSHBIN)
    pdfw.writeGmshBin( p, nfo,
                       g_inputdeck.get< tag::selected, tag::pdfctr >() );
  else if (filetype == tk::ctr::PDFFileType::EXODUSII)
    pdfw.writeExodusII( p, nfo,
                        m_it,
                        g_inputdeck.get< tag::selected, tag::pdfctr >() );
}

void
Distributor::writeTriPDF( const tk::TriPDF& p,
                          tk::ctr::Moment m,
                          std::size_t idx )
// *****************************************************************************
// Write trivariate PDF to file
//! \param[in] p Trivariate PDF to output
//! \param[in] m ORDINARY or CENTRAL PDF we are writing
//! \param[in] idx Index of the PDF of all ordinary or central PDFs requested
//! \author J. Bakosi
// *****************************************************************************
{
  // Get PDF metadata
  const auto nfo =
    tk::ctr::pdfInfo< 3 >( g_inputdeck.get< tag::discr, tag::binsize >(),
                           g_inputdeck.get< tag::cmd, tag::io, tag::pdfnames >(),
                           g_inputdeck.get< tag::discr, tag::extent >(),
                           g_inputdeck.get< tag::pdf >(),
                           m,
                           idx );

  // Construct PDF file name: base name + '_' + pdf name
  std::string filename =
    g_inputdeck.get< tag::cmd, tag::io, tag::pdf >() + '_' + nfo.name;

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
    pdfw.writeTxt( p, nfo );
  else if (filetype == tk::ctr::PDFFileType::GMSHTXT)
     pdfw.writeGmshTxt( p, nfo,
                        g_inputdeck.get< tag::selected, tag::pdfctr >() );
  else if (filetype == tk::ctr::PDFFileType::GMSHBIN)
     pdfw.writeGmshBin( p, nfo,
                        g_inputdeck.get< tag::selected, tag::pdfctr >() );
  else if (filetype == tk::ctr::PDFFileType::EXODUSII)
    pdfw.writeExodusII( p, nfo,
                        m_it,
                        g_inputdeck.get< tag::selected, tag::pdfctr >() );
}

void
Distributor::outUniPDF()
// *****************************************************************************
// Output all requested univariate PDFs to file(s)
//! \author J. Bakosi
// *****************************************************************************
{
  std::size_t idx = 0;
  for (const auto& p : m_ordupdf)
    writeUniPDF( p, tk::ctr::Moment::ORDINARY, idx++ );
  idx = 0;
  for (const auto& p : m_cenupdf)
    writeUniPDF( p, tk::ctr::Moment::CENTRAL, idx++ );
}

void
Distributor::outBiPDF()
// *****************************************************************************
// Output all requested bivariate PDFs to file(s)
//! \return Number of PDFs written
//! \author J. Bakosi
// *****************************************************************************
{
  std::size_t idx = 0;
  for (const auto& p : m_ordbpdf)
    writeBiPDF( p, tk::ctr::Moment::ORDINARY, idx++ );
  idx = 0;
  for (const auto& p : m_cenbpdf) {
    writeBiPDF( p, tk::ctr::Moment::CENTRAL, idx++ );
  }
}

void
Distributor::outTriPDF()
// *****************************************************************************
// Output all requested trivariate PDFs to file(s)
//! \return Number of PDFs written
//! \author J. Bakosi
// *****************************************************************************
{
  std::size_t idx = 0;
  for (const auto& p : m_ordtpdf) {
    writeTriPDF( p, tk::ctr::Moment::ORDINARY, idx++ );
  }
  idx = 0;
  for (const auto& p : m_centpdf) {
    writeTriPDF( p, tk::ctr::Moment::CENTRAL, idx++ );
  }
}

void
Distributor::evaluateTime()
// *****************************************************************************
// Evaluate time step, compute new time step size, decide if it is time to quit
//! \author J. Bakosi
// *****************************************************************************
{
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();

  // Increase number of iterations taken
  ++m_it;
  // Advance physical time
  m_t += m_dt;
  // Truncate the size of last time step
  if (m_t > term) m_t = term;
  // Compute size of next time step
  m_dt = computedt();
  // Echo one-liner info on time step
  report();

  // Finish if either max iterations or max time reached 
  if ( std::fabs(m_t-term) > eps && m_it < nstep ) {

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

      // Re-activate SDAG-wait for estimation of ordinary stats for next step
      wait4ord();
      // Re-activate SDAG-wait for estimation of central moments for next step
      wait4cen();

      // Selectively re-activate SDAG-wait for estimation of PDFs for next step
      if ( !(m_it % g_inputdeck.get< tag::interval, tag::pdf >()) ) wait4pdf();
    }

    // Continue with next time step with all integrators
    m_intproxy.advance( m_dt, m_t, m_it, m_moments );

  } else finish();
}

void
Distributor::finish()
// *****************************************************************************
// Normal finish of time stepping
//! \author J. Bakosi
// *****************************************************************************
{
  // Print out reason for stopping
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

void
Distributor::nostat()
// *****************************************************************************
//  Charm++ reduction target enabling shortcutting sync points if no stats
//! \details This reduction target is called if there are no statistics nor PDFs
//!   to be estimated and thus some synchronization points can be skipped. Upon
//!   this call we simply finish up the time step as usual.
//! \author J. Bakosi
// *****************************************************************************
{
  evaluateTime();
}

void
Distributor::header() const
// *****************************************************************************
// Print out time integration header
//! \author J. Bakosi
// *****************************************************************************
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
// *****************************************************************************
// Print out one-liner report on time step
//! \author J. Bakosi
// *****************************************************************************
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

    // Reset output indicators
    m_output.get< tag::stat >() = false;
    m_output.get< tag::pdf >() = false;

    m_print << std::endl;
  }
}

#include "NoWarning/distributor.def.h"
