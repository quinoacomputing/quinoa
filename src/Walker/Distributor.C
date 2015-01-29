//******************************************************************************
/*!
  \file      src/Walker/Distributor.C
  \author    J. Bakosi
  \date      Thu 29 Jan 2015 09:20:50 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Distributor drives the time integration of differential equations
  \details   Distributor drives the time integration of differential equations.
    The implementation uses the Charm++ runtime system and is fully asynchronous,
    overlapping computation, communication as well I/O. The algorithm utilizes
    the structured dagger (SDAG) Charm++ functionality. The high-level overview
    of the algorithm structure and how it interfaces with Charm++ is discussed
    in the Charm++ interface file src/Walker/distributor.ci.
*/
//******************************************************************************

#include <Distributor.h>
#include <Integrator.h>
#include <DiffEqStack.h>
#include <TxtStatWriter.h>
#include <PDFWriter.h>
#include <walker.decl.h>
#include <flip_map.h>

extern CProxy_Main mainProxy;

using walker::Distributor;

Distributor::Distributor( const ctr::CmdLine& cmdline ) :
  m_print( cmdline.get< tag::verbose >() ? std::cout : std::clog ),
  m_count( 0, 0, 0, 0, 0, 0 ),
  m_output( false, false ),
  m_it( 0 ),
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
  // Compute load distribution given total work and specified virtualization
  uint64_t chunksize, remainder;
  computeLoadDistribution( chunksize, remainder );

  // Print out info on what will be done and how
  info( chunksize, remainder );

  // Start timer measuring total integration time
  m_timer.emplace_back();

  // Compute size of initial time step
  const auto dt = computedt();

  // Activate SDAG-wait for estimation of ordinary statistics
  wait4ord();
  // Activate SDAG-wait for estimation of central moments
  wait4cen();
  // Activate SDAG-wait for estimation of PDFs at select times
  if ( !(m_it % g_inputdeck.get< tag::interval, tag::pdf >()) ) wait4pdf();

  // Fire up asynchronous differential equation integrators
  for (uint64_t i = 1; i < m_count.get< tag::chare >(); ++i)
    m_proxy.push_back( CProxyInt::ckNew( thisProxy, chunksize, dt, m_it ) );
  m_proxy.push_back(CProxyInt::ckNew(thisProxy, chunksize+remainder, dt, m_it));
}

void
Distributor::computeLoadDistribution( uint64_t& chunksize, uint64_t& remainder )
//******************************************************************************
//  Compute load distribution for given total work and virtualization
//! \param[inout] chunksize Chunk size, see detailed description
//! \param[inout] remainder Remainder, see detailed description
//! \details Compute load distibution (number of chares and chunksize) based on
//!   total work (total number of particles) and virtualization
//!
//!   The virtualization parameter, specified by the user, is a real number
//!   between 0.0 and 1.0, inclusive, which controls the degree of
//!   virtualization or over-decomposition. Independent of the value of
//!   virtualization the work is approximately evenly distributed among the
//!   available processing elements. For zero virtualization (no
//!   over-decomposition), the work is simply decomposed into total_work/numPEs,
//!   which yields the smallest number of Charm++ chares and the largest chunks
//!   of work units. The other extreme is unity virtualization, which decomposes
//!   the total work into the smallest size work units possible, yielding the
//!   largest number of Charm++ chares. Obviously, the optimum will be between
//!   0.0 and 1.0, depending on the problem.
//!
//!   The formula implemented uses the simplest (linear) relationship between
//!   the virtualization parameter and the number of work units with the
//!   extremes described above. The formula is given by
//!
//!   chunksize = (1 - n) * v + n;
//!
//!   where
//!
//!    - n = npar/npes
//!
//!    - npar = number of particles, representing the total work
//!
//!    - npes = number of hardware processing elements
//!
//! \author J. Bakosi
//******************************************************************************
{
  // Get virtualization parameter
  const auto v = g_inputdeck.get< tag::cmd, tag::virtualization >();
  Assert( v > -std::numeric_limits< tk::real >::epsilon() &&
          v < 1.0+std::numeric_limits< tk::real >::epsilon(),
          "Virtualization parameter must be between [0.0...1.0]" );

  // Get total number of particles (represents total work)
  const auto npar = g_inputdeck.get< tag::discr, tag::npar >();

  // Query number of processing elements
  const auto npe = CkNumPes();

  // Compute minimum number of work units
  const auto n = npar/npe;

  // Compute work unit size based on the linear formula above
  chunksize = (1.0-n)*v + n;

  // Compute number of work units with size computed ignoring remainder
  m_count.get< tag::chare >() = npar / chunksize;

  // Compute remainder of work if the above number of units were to be created
  remainder = npar - m_count.get< tag::chare >() * chunksize;

  // Redistribute remainder among the work units for a more equal distribution
  chunksize += remainder / m_count.get< tag::chare >();

  // Compute new remainder (after redistribution of the previous remainder)
  remainder = npar - m_count.get< tag::chare >() * chunksize;
}

void
Distributor::info( uint64_t chunksize, uint64_t remainder ) const
//******************************************************************************
//  Print information at startup
//! \param[in] chunksize Chunk size, see computeLoadDistribution()
//! \param[in] remainder Remainder, see computeLoadDistribution()
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
    m_print.section( "Title", g_inputdeck.get< tag::title >() );

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
  m_print.item( "Load (number of particles)",
                g_inputdeck.get< tag::discr, tag::npar >() );
  m_print.item( "Number of processing elements", CkNumPes() );
  m_print.item( "Number of work units",
                std::to_string( m_count.get< tag::chare >() ) + " (" +
                std::to_string( m_count.get< tag::chare >()-1 ) + "*" +
                std::to_string( chunksize ) + "+" +
                std::to_string( chunksize+remainder ) + ")" );

  // Print out time integration header
  if (g_inputdeck.get< tag::discr, tag::nstep >()) {
    header();
    tk::TxtStatWriter sw( !m_nameOrdinary.empty() || !m_nameCentral.empty() ?
                          g_inputdeck.get< tag::cmd, tag::io, tag::stat >() :
                          std::string() );
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
Distributor::init()
//******************************************************************************
// Wait for all integrators to finish initialization
//! \author  J. Bakosi
//******************************************************************************
{
  // Increase number of integrators completing initialization
  ++m_count.get< tag::init >();

  // Wait for all integrators completing initialization
  if (m_count.get< tag::init >() == m_count.get< tag::chare >())
    mainProxy.timestamp( "Initial conditions", m_timer[0].dsec() );
}

void
Distributor::estimateOrd( const std::vector< tk::real >& ord )
//******************************************************************************
// Wait for all integrators to finish accumulation of the ordinary moments
//! \param[in] ord Partially accumulated ordinary moments contributed by caller
//! \author J. Bakosi
//******************************************************************************
{
  // Increase number of integrators completing the accumulation of the ordinary
  // moments
  ++m_count.get< tag::ordinary >();

  // Add contribution from PE to total sums, i.e., u[i] += v[i] for all i
  for (std::size_t i=0; i<m_ordinary.size(); ++i) m_ordinary[i] += ord[i];

  // Wait for all integrators completing accumulation of ordinary moments
  if (m_count.get< tag::ordinary >() == m_count.get< tag::chare >()) {
    // Finish computing moments, i.e., divide sums by the number of samples
    for (auto& m : m_ordinary) m /= g_inputdeck.get< tag::discr, tag::npar >();
    // Activate SDAG trigger signaling that ordinary moments have been estimated
    estimateOrdDone();
  }
}

void
Distributor::estimateCen( const std::vector< tk::real >& cen )
//******************************************************************************
// Wait for all integrators to finish accumulation of the central moments
//! \param[in] cen Partially accumulated central moments contributed by caller
//! \author J. Bakosi
//******************************************************************************
{
  // Increase number of integrators completing the accumulation of the central
  // moments
  ++m_count.get< tag::central >();

  // Add contribution from PE to total sums, i.e., u[i] += v[i] for all i
  for (std::size_t i=0; i<m_central.size(); ++i) m_central[i] += cen[i];

  // Wait for all integrators completing accumulation of central moments
  if (m_count.get< tag::central >() == m_count.get< tag::chare >()) {
    // Finish computing moments, i.e., divide sums by the number of samples
    for (auto& m : m_central) m /= g_inputdeck.get< tag::discr, tag::npar >();
    // Activate SDAG trigger signaling that central moments have been estimated
    estimateCenDone();
  }
}

void
Distributor::estimateOrdPDF( const std::vector< tk::UniPDF >& updf,
                             const std::vector< tk::BiPDF >& bpdf,
                             const std::vector< tk::TriPDF >& tpdf )
//******************************************************************************
// Wait for all integrators to finish accumulation of ordinary PDFs
//! \param[in] updf Partially accumulated univariate PDFs contributed by caller
//! \param[in] bpdf Partially accumulated bivariate PDFs contributed by caller
//! \param[in] tpdf Partially accumulated trivariate PDFs contributed by caller
//! \author J. Bakosi
//******************************************************************************
{
  // Increase number of integrators completing the accumulation of ordinary PDFs
  ++m_count.get< tag::ordpdf >();

  // Add contribution from PE to total sums
  std::size_t i = 0;
  m_ordupdf.resize( updf.size() );
  for (const auto& p : updf) m_ordupdf[i++].addPDF( p );

  i = 0;
  m_ordbpdf.resize( bpdf.size() );
  for (const auto& p : bpdf) m_ordbpdf[i++].addPDF( p );

  i = 0;
  m_ordtpdf.resize( tpdf.size() );
  for (const auto& p : tpdf) m_ordtpdf[i++].addPDF( p );

  // Wait for all integrators completing accumulation of ordinary PDFs
  if (m_count.get< tag::ordpdf >() == m_count.get< tag::chare >()) {
    // Activate SDAG trigger signaling that ordinary PDFs have been estimated
    estimateOrdPDFDone();
  }
}

void
Distributor::estimateCenPDF( const std::vector< tk::UniPDF >& updf,
                             const std::vector< tk::BiPDF >& bpdf,
                             const std::vector< tk::TriPDF >& tpdf )
//******************************************************************************
// Wait for all integrators to finish accumulation of central PDFs
//! \param[in] updf Partially accumulated univariate PDFs contributed by caller
//! \param[in] bpdf Partially accumulated bivariate PDFs contributed by caller
//! \param[in] tpdf Partially accumulated trivariate PDFs contributed by caller
//! \author J. Bakosi
//******************************************************************************
{
  // Increase number of integrators completing the accumulation of PDFs
  ++m_count.get< tag::cenpdf >();

  // Add contribution from PE to total sums
  std::size_t i = 0;
  m_cenupdf.resize( updf.size() );
  for (const auto& p : updf) m_cenupdf[i++].addPDF( p );

  i = 0;
  m_cenbpdf.resize( bpdf.size() );
  for (const auto& p : bpdf) m_cenbpdf[i++].addPDF( p );

  i = 0;
  m_centpdf.resize( tpdf.size() );
  for (const auto& p : tpdf) m_centpdf[i++].addPDF( p );

  // Wait for all integrators completing accumulation of PDFs
  if (m_count.get< tag::cenpdf >() == m_count.get< tag::chare >()) {
    // Activate SDAG trigger signaling that central PDFs have been estimated
    estimateCenPDFDone();
  }
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
                          std::string(), std::ios_base::app );
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
    std::size_t n = outUniPDF();       // Output univariate PDFs to file(s)
    n += outBiPDF();                   // Output bivariate PDFs to file(s)
    n += outTriPDF();                  // Output trivariate PDFs to file(s)
    if (n) m_output.get< tag::pdf >() = true; // Signal that PDFs were written
  }
}

void
Distributor::writeUniPDF( const tk::UniPDF& p, std::size_t& cnt )
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
                      g_inputdeck.get< tag::selected, tag::float_format >(),
                      g_inputdeck.get< tag::discr, tag::precision >() );

  // Output PDF
  pdfw.writeTxt( p, info );
}

void
Distributor::writeBiPDF( const tk::BiPDF& p, std::size_t& cnt )
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
                      g_inputdeck.get< tag::selected, tag::float_format >(),
                      g_inputdeck.get< tag::discr, tag::precision >() );

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
    pdfw.writeExodusII( p, info, static_cast< int >( m_it ),
                        g_inputdeck.get< tag::selected, tag::pdfctr >() );
}

void
Distributor::writeTriPDF( const tk::TriPDF& p, std::size_t& cnt )
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
                      g_inputdeck.get< tag::selected, tag::float_format >(),
                      g_inputdeck.get< tag::discr, tag::precision >() );

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
    pdfw.writeExodusII( p, info, static_cast< int >( m_it ),
                        g_inputdeck.get< tag::selected, tag::pdfctr >() );
}

std::size_t
Distributor::outUniPDF()
//******************************************************************************
// Output all requested univariate PDFs to file(s)
//! \return Number of PDFs written
//! \author J. Bakosi
//******************************************************************************
{
  std::size_t cnt = 0;
  for (const auto& p : m_ordupdf) writeUniPDF( p, cnt );
  for (const auto& p : m_cenupdf) writeUniPDF( p, cnt );
  return cnt;
}

std::size_t
Distributor::outBiPDF()
//******************************************************************************
// Output all requested bivariate PDFs to file(s)
//! \return Number of PDFs written
//! \author J. Bakosi
//******************************************************************************
{
  std::size_t cnt = 0;
  for (const auto& p : m_ordbpdf) writeBiPDF( p, cnt );
  for (const auto& p : m_cenbpdf) writeBiPDF( p, cnt );
  return cnt;
}

std::size_t
Distributor::outTriPDF()
//******************************************************************************
// Output all requested trivariate PDFs to file(s)
//! \return Number of PDFs written
//! \author J. Bakosi
//******************************************************************************
{
  std::size_t cnt = 0;
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
  const auto dt = computedt();

  // Advance physical time
  m_t += dt;

  // Get physical time at which to terminate
  const auto term = g_inputdeck.get< tag::discr, tag::term >();

  // Truncate the size of last time step
  if (m_t > term) m_t = term;

  // Echo one-liner info on time step
  report();

  // Finish if either max iterations or max time reached 
  if (std::fabs(m_t - term) > std::numeric_limits< tk::real >::epsilon() &&
    m_it < g_inputdeck.get< tag::discr, tag::nstep >()) {

    // Zero statistics counters and accumulators
    m_count.get< tag::ordinary >() = m_count.get< tag::central >() = 0;
    std::fill( begin(m_ordinary), end(m_ordinary), 0.0 );
    std::fill( begin(m_central), end(m_central), 0.0 );

    // Re-activate SDAG-wait for estimation of ordinary statistics for next step
    wait4ord();
    // Re-activate SDAG-wait for estimation of central moments for next step
    wait4cen();

    // Selectively re-activate SDAG-wait for estimation of PDFs for next step
    if ( !(m_it % g_inputdeck.get< tag::interval, tag::pdf >()) ) {
      // Zero PDF counters and accumulators
      m_count.get< tag::ordpdf >() = m_count.get< tag::cenpdf >() = 0;
      for (auto& p : m_ordupdf) p.zero();
      for (auto& p : m_ordbpdf) p.zero();
      for (auto& p : m_ordtpdf) p.zero();
      for (auto& p : m_cenupdf) p.zero();
      for (auto& p : m_cenbpdf) p.zero();
      for (auto& p : m_centpdf) p.zero();
      wait4pdf();
    }

    // Continue with next time step with all integrators
    for (auto& p : m_proxy) p.advance( dt, m_it );

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
            << g_inputdeck.get< tag::discr, tag::dt >() << "  "
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

    m_print << '\n';
  }
}

#include <distributor.def.h>
