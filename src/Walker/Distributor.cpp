// *****************************************************************************
/*!
  \file      src/Walker/Distributor.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

#include "NoWarning/format.hpp"

#include "Macro.hpp"
#include "Print.hpp"
#include "Tags.hpp"
#include "StatCtr.hpp"
#include "Exception.hpp"
#include "Particles.hpp"
#include "LoadDistributor.hpp"
#include "Distributor.hpp"
#include "Integrator.hpp"
#include "DiffEqStack.hpp"
#include "TxtStatWriter.hpp"
#include "PDFReducer.hpp"
#include "PDFWriter.hpp"
#include "Options/PDFFile.hpp"
#include "Options/PDFPolicy.hpp"
#include "Walker/InputDeck/InputDeck.hpp"
#include "NoWarning/walker.decl.h"

extern CProxy_Main mainProxy;

using walker::Distributor;

Distributor::Distributor() :
  m_output( { false, false } ),
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
// *****************************************************************************
{
  // Get command line object reference
  const auto& cmd = g_inputdeck.get< tag::cmd >();

  // Compute load distribution given total work (= number of particles) and
  // user-specified virtualization
  uint64_t chunksize = 0, remainder = 0;
  auto nchare = tk::linearLoadDistributor(
                  cmd.get< tag::virtualization >(),
                  g_inputdeck.get< tag::discr, tag::npar >(),
                  CkNumPes(),
                  chunksize,
                  remainder );
  Assert( chunksize != 0, "Chunksize must not be zero" );

  // Compute total number of particles distributed over all workers. Note that
  // this number will not necessarily be the same as given by the user, coming
  // from g_inputdeck.get< tag::discr, tag::npar >(), since each Charm++ chare
  // array element constructor takes this chunksize argument, which equals the
  // number of particles the array element (worker) will work on.
  m_npar = static_cast< tk::real >( nchare * chunksize );

  auto print = printer();

  // Print out info on what will be done and how
  info( print, chunksize, nchare );

  // Output header for statistics output file
  tk::TxtStatWriter sw( !m_nameOrdinary.empty() || !m_nameCentral.empty() ?
                        cmd.get< tag::io, tag::stat >() :
                        std::string(),
                        g_inputdeck.get< tag::flformat, tag::stat >(),
                        g_inputdeck.get< tag::prec, tag::stat >() );
  sw.header( m_nameOrdinary, m_nameCentral, m_tables.first );

  // Print out time integration header
  print.endsubsection();
  print.diag( "Starting time stepping ..." );
  header( print );

  // Start timer measuring total integration time
  m_timer.emplace_back();

  // Construct and initialize map of statistical moments
  for (const auto& product : g_inputdeck.get< tag::stat >())
    m_moments[ product ] = 0.0;

  // Activate SDAG-wait for estimation of ordinary statistics
  thisProxy.wait4ord();
  // Activate SDAG-wait for estimation of PDFs at select times
  thisProxy.wait4pdf();

  // Create statistics merger chare group collecting chare contributions
  CProxy_Collector collproxy = CProxy_Collector::ckNew( thisProxy );

  // Create partcle writer Charm++ chare nodegroup
  tk::CProxy_ParticleWriter particlewriter =
    tk::CProxy_ParticleWriter::ckNew( cmd.get< tag::io, tag::particles >() );

  // Fire up asynchronous differential equation integrators
  m_intproxy =
    CProxy_Integrator::ckNew( thisProxy, collproxy, particlewriter, chunksize,
                              static_cast<int>( nchare ) );
}

void
Distributor::info( const WalkerPrint& print,
                   uint64_t chunksize,
                   std::size_t nchare )
// *****************************************************************************
//  Print information at startup
//! \param[in] print Pretty printer object to use for printing
//! \param[in] chunksize Chunk size, see Base/LoadDistribution.h
//! \param[in] nchare Total number of Charem++ Integrator chares doing work
// *****************************************************************************
{
  // Get command line object reference
  const auto& cmd = g_inputdeck.get< tag::cmd >();

  print.part( "Factory" );

  // Print out info data layout
  print.list( "Particle properties data layout (CMake: PARTICLE_DATA_LAYOUT)",
              std::list< std::string >{ tk::Particles::layout() } );

  // Re-create differential equations stack for output
  DiffEqStack stack;

  // Print out information on factory
  print.eqlist( "Registered differential equations",
                stack.factory(), stack.ntypes() );
  print.endpart();

  // Instantiate tables to sample and output to statistics file
  m_tables = stack.tables();

  // Print out information on problem
  print.part( "Problem" );

  // Print out info on problem title
  if ( !g_inputdeck.get< tag::title >().empty() )
    print.title( g_inputdeck.get< tag::title >() );

  // Print out info on settings of selected differential equations
  print.diffeqs( "Differential equations integrated", stack.info() );

  // Print out info on RNGs selected
  // ...

  // Print I/O filenames
  print.section( "Output filenames" );
  if (!g_inputdeck.get< tag::stat >().empty())
    print.item( "Statistics", cmd.get< tag::io, tag::stat >() );
  if (!g_inputdeck.get< tag::pdf >().empty())
    print.item( "PDF", cmd.get< tag::io, tag::pdf >() );
  if (!g_inputdeck.get< tag::param, tag::position, tag::depvar >().empty())
    print.item( "Particle positions", cmd.get< tag::io, tag::particles >() );

  // Print discretization parameters
  print.section( "Discretization parameters" );
  print.item( "Number of time steps",
              g_inputdeck.get< tag::discr, tag::nstep >() );
  print.item( "Terminate time",
              g_inputdeck.get< tag::discr, tag::term >() );
  print.item( "Initial time step size",
              g_inputdeck.get< tag::discr, tag::dt >() );

  // Print output intervals
  print.section( "Output intervals" );
  const auto& interval = g_inputdeck.get< tag::output, tag::iter >();
  print.item( "TTY", interval.get< tag::tty >() );
  if (!g_inputdeck.get< tag::stat >().empty())
    print.item( "Statistics", interval.get< tag::stat >() );
  if (!g_inputdeck.get< tag::pdf >().empty())
    print.item( "PDF", interval.get< tag::pdf >() );
  if (!g_inputdeck.get< tag::param, tag::position, tag::depvar >().empty())
    print.item( "Particles", interval.get< tag::particles >() );

  // Print out statistics estimated
  print.statistics( "Statistical moments and distributions" );

  // Print out info on load distirubtion
  print.section( "Load distribution" );
  print.item( "Virtualization [0.0...1.0]",
              g_inputdeck.get< tag::cmd, tag::virtualization >() );
  print.item( "Number of work units", nchare );
  print.item( "User load (# of particles)",
              g_inputdeck.get< tag::discr, tag::npar >() );
  print.item( "Chunksize (load per work unit)", chunksize );
  print.item( "Actual load (# of particles)",
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
// *****************************************************************************
{
  // Simply return a constant user-defined dt for now
  return g_inputdeck.get< tag::discr, tag::dt >();
}

void
Distributor::estimateOrd( tk::real* ord, [[maybe_unused]] int n )
// *****************************************************************************
// Estimate ordinary moments
//! \param[in] ord Ordinary moments (sum) collected over all chares
//! \param[in] n Number of ordinary moments in array ord
// *****************************************************************************
{
  Assert( static_cast<std::size_t>(n) == m_ordinary.size(),
          "Number of ordinary moments contributed not equal to expected" );

  // Add contribution from PE to total sums, i.e., u[i] += v[i] for all i
  for (std::size_t i=0; i<m_ordinary.size(); ++i) m_ordinary[i] += ord[i];

  // Finish computing moments, i.e., divide sums by the number of samples
  // cppcheck-suppress useStlAlgorithm
  for (auto& m : m_ordinary) m /= m_npar;

  // Activate SDAG trigger signaling that ordinary moments have been estimated
  estimateOrdDone();
}

void
Distributor::estimateCen( tk::real* cen, [[maybe_unused]] int n )
// *****************************************************************************
// Estimate ordinary moments
//! \param[in] cen Central moments (sum) collected over all chares
//! \param[in] n Number of central moments in array cen
// *****************************************************************************
{
  Assert( static_cast<std::size_t>(n) == m_central.size(),
          "Number of central moments contributed not equal to expected" );

  // Add contribution from PE to total sums, i.e., u[i] += v[i] for all i
  for (std::size_t i=0; i<m_central.size(); ++i) m_central[i] += cen[i];

  // Finish computing moments, i.e., divide sums by the number of samples
  // cppcheck-suppress useStlAlgorithm
  for (auto& m : m_central) m /= m_npar;

  // Activate SDAG trigger signaling that central moments have been estimated
  estimateCenDone();
}

void
Distributor::estimateOrdPDF( CkReductionMsg* msg )
// *****************************************************************************
// Estimate ordinary PDFs
//! \param[in] msg Serialized vectors of uni-, bi-, and tri-variate PDFs
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
//! \param[in] msg Serialized vectors of uni-, bi-, and tri-variate PDFs
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
// *****************************************************************************
{
  // lambda to sample tables to write to statistics file
  auto extra = [this]() -> std::vector< tk::real > {
    std::vector< tk::real > x( m_tables.second.size() );
    std::size_t j = 0;
    for (const auto& t : m_tables.second) x[ j++ ] = tk::sample<1>(m_t,t)[0];
    return x;
  };

  // Append statistics file at selected times
  if (!((m_it+1) % g_inputdeck.get< tag::output, tag::iter, tag::stat >())) {
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
// *****************************************************************************
{
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto pdffreq = g_inputdeck.get< tag::output, tag::iter, tag::pdf >();

  // output PDFs at t=0 (regardless of whether it was requested), or at
  // selected times, or in the last time step (regardless of whether it was
  // requested
  if ( m_it == 0 ||
       !((m_it+1) % pdffreq) ||
       (std::fabs(m_t+m_dt-term) < eps || (m_it+1) >= nstep) )
  {
    // Generate iteration count and physical time for PDF output. In the first
    // iteration, the particles are NOT advanced, see Integration::advance(),
    // and we write it=0 and time=0.0 into the PDF files. For the rest of the
    // iterations we write the iteration count and the physical time
    // corresponding to the iteration just completed.
    auto it = m_it == 0 ? m_it : m_it + 1;
    auto t = m_it == 0 ? m_t : m_t + m_dt;

    outUniPDF( it, t );                 // Output univariate PDFs to file(s)
    outBiPDF( it, t );                  // Output bivariate PDFs to file(s)
    outTriPDF( it, t );                 // Output trivariate PDFs to file(s)
    m_output.get< tag::pdf >() = true;  // Signal that PDFs were written
  }
}

void
Distributor::writeUniPDF( std::uint64_t it,
                          tk::real t,
                          const tk::UniPDF& p,
                          tk::ctr::Moment m,
                          std::size_t idx )
// *****************************************************************************
// Write univariate PDF to file
//! \param[in] it Iteration count to write in output file
//! \param[in] t Physical time to write in output file
//! \param[in] p Univariate PDF to output
//! \param[in] m ORDINARY or CENTRAL PDF we are writing
//! \param[in] idx Index of the PDF of all ordinary or central PDFs requested
// *****************************************************************************
{
  // Get PDF metadata
  const auto nfo =
    tk::ctr::pdfInfo< 1 >( g_inputdeck.get< tag::discr, tag::binsize >(),
                           g_inputdeck.get< tag::cmd, tag::io, tag::pdfnames >(),
                           g_inputdeck.get< tag::discr, tag::extent >(),
                           g_inputdeck.get< tag::pdf >(),
                           m,
                           idx,
                           it,
                           t );

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
Distributor::writeBiPDF( std::uint64_t it,
                         tk::real t,
                         const tk::BiPDF& p,
                         tk::ctr::Moment m,
                         std::size_t idx )
// *****************************************************************************
// Write bivariate PDF to file
//! \param[in] it Iteration count to write in output file
//! \param[in] t Physical time to write in output file
//! \param[in] p Bivariate PDF to output
//! \param[in] m ORDINARY or CENTRAL PDF we are writing
//! \param[in] idx Index of the PDF of all ordinary or central PDFs requested
// *****************************************************************************
{
  // Get PDF metadata
  const auto nfo =
    tk::ctr::pdfInfo< 2 >( g_inputdeck.get< tag::discr, tag::binsize >(),
                           g_inputdeck.get< tag::cmd, tag::io, tag::pdfnames >(),
                           g_inputdeck.get< tag::discr, tag::extent >(),
                           g_inputdeck.get< tag::pdf >(),
                           m,
                           idx,
                           it,
                           t );

  // Construct PDF file name: base name + '_' + pdf name
  std::string filename =
    g_inputdeck.get< tag::cmd, tag::io, tag::pdf >() + '_' + nfo.name;

  // Augment PDF filename by time stamp if PDF output file policy is multiple
  if (g_inputdeck.get< tag::selected, tag::pdfpolicy >() ==
      tk::ctr::PDFPolicyType::MULTIPLE)
    filename += '_' + std::to_string( m_t );

  const auto& filetype = g_inputdeck.get< tag::selected, tag::filetype >();

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
                        g_inputdeck.get< tag::selected, tag::pdfctr >() );
}

void
Distributor::writeTriPDF( std::uint64_t it,
                          tk::real t,
                          const tk::TriPDF& p,
                          tk::ctr::Moment m,
                          std::size_t idx )
// *****************************************************************************
// Write trivariate PDF to file
//! \param[in] it Iteration count to write in output file
//! \param[in] t Physical time to write in output file
//! \param[in] p Trivariate PDF to output
//! \param[in] m ORDINARY or CENTRAL PDF we are writing
//! \param[in] idx Index of the PDF of all ordinary or central PDFs requested
// *****************************************************************************
{
  // Get PDF metadata
  const auto nfo =
    tk::ctr::pdfInfo< 3 >( g_inputdeck.get< tag::discr, tag::binsize >(),
                           g_inputdeck.get< tag::cmd, tag::io, tag::pdfnames >(),
                           g_inputdeck.get< tag::discr, tag::extent >(),
                           g_inputdeck.get< tag::pdf >(),
                           m,
                           idx,
                           it,
                           t );

  // Construct PDF file name: base name + '_' + pdf name
  std::string filename =
    g_inputdeck.get< tag::cmd, tag::io, tag::pdf >() + '_' + nfo.name;

  // Augment PDF filename by time stamp if PDF output file policy is multiple
  if (g_inputdeck.get< tag::selected, tag::pdfpolicy >() ==
      tk::ctr::PDFPolicyType::MULTIPLE)
    filename += '_' + std::to_string( m_t );

  const auto& filetype = g_inputdeck.get< tag::selected, tag::filetype >();

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
                        g_inputdeck.get< tag::selected, tag::pdfctr >() );
}

void
Distributor::outUniPDF( std::uint64_t it, tk::real t )
// *****************************************************************************
// Output all requested univariate PDFs to file(s)
//! \param[in] it Iteration count to write in output file
//! \param[in] t Physical time to write in output file
// *****************************************************************************
{
  std::size_t idx = 0;
  for (const auto& p : m_ordupdf)
    writeUniPDF( it, t, p, tk::ctr::Moment::ORDINARY, idx++ );
  idx = 0;
  for (const auto& p : m_cenupdf)
    writeUniPDF( it, t, p, tk::ctr::Moment::CENTRAL, idx++ );
}

void
Distributor::outBiPDF( std::uint64_t it, tk::real t )
// *****************************************************************************
// Output all requested bivariate PDFs to file(s)
//! \param[in] it Iteration count to write in output file
//! \param[in] t Physical time to write in output file
//! \return Number of PDFs written
// *****************************************************************************
{
  std::size_t idx = 0;
  for (const auto& p : m_ordbpdf)
    writeBiPDF( it, t, p, tk::ctr::Moment::ORDINARY, idx++ );
  idx = 0;
  for (const auto& p : m_cenbpdf) {
    writeBiPDF( it, t, p, tk::ctr::Moment::CENTRAL, idx++ );
  }
}

void
Distributor::outTriPDF( std::uint64_t it, tk::real t )
// *****************************************************************************
// Output all requested trivariate PDFs to file(s)
//! \param[in] it Iteration count to write in output file
//! \param[in] t Physical time to write in output file
//! \return Number of PDFs written
// *****************************************************************************
{
  std::size_t idx = 0;
  for (const auto& p : m_ordtpdf) {
    writeTriPDF( it, t, p, tk::ctr::Moment::ORDINARY, idx++ );
  }
  idx = 0;
  for (const auto& p : m_centpdf) {
    writeTriPDF( it, t, p, tk::ctr::Moment::CENTRAL, idx++ );
  }
}

void
Distributor::evaluateTime()
// *****************************************************************************
// Evaluate time step, compute new time step size, decide if it is time to quit
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
      thisProxy.wait4ord();
      // Re-activate SDAG-wait for estimation of PDFs for next step
      thisProxy.wait4pdf();
    }

    // Continue with next time step with all integrators
    m_intproxy.advance( m_dt, m_t, m_it, m_moments );

  } else finish();
}

void
Distributor::finish()
// *****************************************************************************
// Normal finish of time stepping
// *****************************************************************************
{
  // Print out reason for stopping
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();

  auto print = printer();

  print.endsubsection();
  if (m_it >= g_inputdeck.get< tag::discr, tag::nstep >())
     print.note( "Normal finish, maximum number of iterations reached: " +
                 std::to_string( nstep ) );
   else 
     print.note( "Normal finish, maximum time reached: " +
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
// *****************************************************************************
{
  evaluateTime();
}

void
Distributor::header( const WalkerPrint& print ) const
// *****************************************************************************
// Print out time integration header
//! \param[in] print Pretty printer object to use for printing
// *****************************************************************************
{
  print.inthead( "Time integration", "Differential equations testbed",
    "Legend: it - iteration count\n"
    "         t - time\n"
    "        dt - time step size\n"
    "       ETE - estimated time elapsed (h:m:s)\n"
    "       ETA - estimated time for accomplishment (h:m:s)\n"
    "       out - status flags, legend:\n"
    "             s - statistics output\n"
    "             p - PDFs output\n"
    "             x - particle positions output\n",
    "\n      it             t            dt        ETE        ETA   out\n"
      " ---------------------------------------------------------------\n" );
}

void
Distributor::report()
// *****************************************************************************
// Print out one-liner report on time step
// *****************************************************************************
{
  if (!(m_it % g_inputdeck.get< tag::output, tag::iter, tag::tty >())) {

  const auto parfreq =
    g_inputdeck.get< tag::output, tag::iter, tag::particles >();
  const auto poseq =
    !g_inputdeck.get< tag::param, tag::position, tag::depvar >().empty();

    // estimated time elapsed and for accomplishment
    tk::Timer::Watch ete, eta;
    m_timer[0].eta( g_inputdeck.get< tag::discr, tag::term >(), m_t,
                    g_inputdeck.get< tag::discr, tag::nstep >(), m_it,
                    ete, eta );

    auto print = printer();

    // Output one-liner
    print << std::setfill(' ') << std::setw(8) << m_it << "  "
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
    if (m_output.get< tag::stat >()) print << 's';
    if (m_output.get< tag::pdf >()) print << 'p';
    if (poseq && !(m_it % parfreq)) print << 'x';

    // Reset output indicators
    m_output.get< tag::stat >() = false;
    m_output.get< tag::pdf >() = false;

    print << std::endl;
  }
}

#include "NoWarning/distributor.def.h"
