//******************************************************************************
/*!
  \file      src/Integrator/Distributor.C
  \author    J. Bakosi
  \date      Mon 13 Oct 2014 10:18:18 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Distributor drives the time integration of differential equations
  \details   Distributor drives the time integration of differential equations
*/
//******************************************************************************

#include <Distributor.h>
#include <Integrator.h>
#include <DiffEqStack.h>
#include <TxtStatWriter.h>
#include <PDFWriter.h>
#include <quinoa.decl.h>
#include <flip_map.h>

extern CProxy_Main mainProxy;

using quinoa::Distributor;

Distributor::Distributor( const ctr::CmdLine& cmdline ) :
  m_print( cmdline.get< tk::tag::verbose >() ? std::cout : std::clog ),
  m_count( 0, 0, 0, 0, 0 ),
  m_output( false, false ),
  m_it( 0 ),
  m_t( 0.0 ),
  m_nameOrdinary( g_inputdeck.momentNames( ctr::ordinary ) ),
  m_nameCentral( g_inputdeck.momentNames( ctr::central ) ),
  m_ordinary( m_nameOrdinary.size(), 0.0 ),
  m_central( m_nameCentral.size(), 0.0 ),
  m_updf( g_inputdeck.npdf< 1 >() ),
  m_bpdf( g_inputdeck.npdf< 2 >() ),
  m_tpdf( g_inputdeck.npdf< 3 >() )
//******************************************************************************
// Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Print out info data layout
  m_print.list( "Particle properties data layout policy (CMake: LAYOUT)",
                std::list< std::string >{ ParProps().major() } );

  // Re-create differential equations stack for output
  DiffEqStack stack;

  // Print out information on factory
  m_print.eqlist( "Registered differential equations", stack.factory(),
                  stack.ntypes() );
  m_print.endpart();

  // Compute load distribution given total work and specified virtualization
  uint64_t chunksize, remainder;
  computeLoadDistribution( chunksize, remainder );

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
  m_print.item( "Input", g_inputdeck.get< tag::cmd, tag::io, tag::input >() );
  m_print.item( "Output", g_inputdeck.get< tag::cmd, tag::io, tag::output >() );
  m_print.item( "Glob", g_inputdeck.get< tag::cmd, tag::io, tag::glob >() );
  m_print.item( "Statistics", g_inputdeck.get< tag::cmd, tag::io, tag::stat >() );
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
  m_print.item( "Dump", g_inputdeck.get< tag::interval, tag::dump>() );
  m_print.item( "Glob", g_inputdeck.get< tag::interval, tag::glob >() );
  m_print.item( "Statistics", g_inputdeck.get< tag::interval, tag::stat >() );
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
    TxtStatWriter sw( !m_nameOrdinary.empty() || !m_nameCentral.empty() ?
                      g_inputdeck.get< tag::cmd, tag::io, tag::stat >() :
                      std::string() );
    sw.header( g_inputdeck.plotOrdinary(), m_nameOrdinary, m_nameCentral );
  }

  // Start timer measuring total integration time
  m_timer.emplace_back();

  // Compute size of initial time step
  const auto dt = computedt();

  // Fire up integrators
  for (uint64_t i = 1; i < m_count.get< tag::chare >(); ++i)
    m_proxy.push_back( CProxyInt::ckNew( thisProxy, chunksize, dt ) );
  m_proxy.push_back( CProxyInt::ckNew( thisProxy, chunksize+remainder, dt ) );
}

void
Distributor::computeLoadDistribution( uint64_t& chunksize, uint64_t& remainder )
//******************************************************************************
//  Compute load distribution for given total work and virtualization
//! \author J. Bakosi
//******************************************************************************
{
  // Compute load distibution (number of chares and chunksize) based on total
  // work (total number of particles) and virtualization

  // The virtualization parameter, specified by the user, is a real number
  // between 0.0 and 1.0, inclusive, which controls the degree of virtualization
  // or over-decomposition. Independent of the value of virtualization the work
  // is approximately evenly distributed among the available processing
  // elements. For zero virtualization (no over-decomposition), the work is
  // simply decomposed into total_work/numPEs, which yields the smallest number
  // of Charm++ chares and the largest chunks of work units. The other extreme
  // is unity virtualization, which decomposes the total work into the smallest
  // size work units possible, yielding the largest number of Charm++ chares.
  // Obviously, the optimum will be between 0.0 and 1.0, depending on the
  // problem.
  //
  // The formula below uses the simplest (linear) relationship between the
  // virtualization parameter and the number of work units with the extremes
  // described above. The formula is given by
  //
  // chunksize = (1 - n) * v + n;
  //
  // where
  //         n = npar/npes
  //      npar = number of particles, representing the total work
  //      npes = number of hardware processing elements

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

tk::real
Distributor::computedt()
//******************************************************************************
// Compute size of next time step
//! \author  J. Bakosi
//******************************************************************************
{
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
  if (m_count.get< tag::init >() == m_count.get< tag::chare >()) {
    mainProxy.timestamp( "Initial conditions", m_timer[0].dsec() );
  }
}

void
Distributor::estimateOrd( const std::vector< tk::real >& ord )
//******************************************************************************
// Wait for all integrators to finish accumulation of the ordinary moments
//! \author  J. Bakosi
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
    // Continue with accumulation for central moments with all integrators
    for (auto& p : m_proxy) p.accumulateCen( m_ordinary );
  }
}

void
Distributor::estimateCen( const std::vector< tk::real >& cen )
//******************************************************************************
// Wait for all integrators to finish accumulation of the central moments
//! \author  J. Bakosi
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

    // Append statistics file at selected times
    if (!(m_it % g_inputdeck.get< tag::interval, tag::stat >())) {
      TxtStatWriter sw( !m_nameOrdinary.empty() || !m_nameCentral.empty() ?
                        g_inputdeck.get< tag::cmd, tag::io, tag::stat >() :
                        std::string(), std::ios_base::app );
      sw.stat( m_it, m_t, m_ordinary, m_central, g_inputdeck.plotOrdinary() );
      m_output.get< tag::stat >() = true;
    }

    // Zero accumulator counters and total-sums for next time step
    m_count.get< tag::ordinary >() = m_count.get< tag::central >() = 0;
    std::fill( begin(m_ordinary), end(m_ordinary), 0.0 );    
    std::fill( begin(m_central), end(m_central), 0.0 );    

    // Continue with accumulation for PDFs with all integrators
    for (auto& p : m_proxy) p.accumulatePDF();
  }
}

void
Distributor::estimatePDF( const std::vector< UniPDF >& oupdf,
                          const std::vector< BiPDF >& obpdf,
                          const std::vector< TriPDF >& otpdf,
                          const std::vector< UniPDF >& cupdf,
                          const std::vector< BiPDF >& cbpdf,
                          const std::vector< TriPDF >& ctpdf )
//******************************************************************************
// Wait for all integrators to finish accumulation of PDFs
//! \author  J. Bakosi
//******************************************************************************
{
  // Increase number of integrators completing the accumulation of PDFs
  ++m_count.get< tag::pdf >();

  // Add contribution from PE to total sums
  std::size_t i=0;
  for (const auto& p : oupdf) m_updf[i++].addPDF( p );
  for (const auto& p : cupdf) m_updf[i++].addPDF( p );
  i = 0;
  for (const auto& p : obpdf) m_bpdf[i++].addPDF( p );
  for (const auto& p : cbpdf) m_bpdf[i++].addPDF( p );
  i = 0;
  for (const auto& p : otpdf) m_tpdf[i++].addPDF( p );
  for (const auto& p : ctpdf) m_tpdf[i++].addPDF( p );

  // Wait for all integrators completing accumulation of PDFs
  if (m_count.get< tag::pdf >() == m_count.get< tag::chare >()) {

    // Output PDFs at selected times
    if ( !(m_it % g_inputdeck.get< tag::interval, tag::pdf >()) &&
         (!m_updf.empty() || !m_bpdf.empty() || !m_tpdf.empty()) )
    {
      outUniPDF();                       // Output univariate PDFs to file(s)
      outBiPDF();                        // Output bivariate PDFs to file(s)
      outTriPDF();                       // Output trivariate PDFs to file(s)
      m_output.get< tag::pdf >() = true; // Signal that PDFs were written
    }

    // Zero accumulator counter and PDFs for next time step
    m_count.get< tag::pdf >() = 0;
    for (auto& p : m_updf) p.zero();
    for (auto& p : m_bpdf) p.zero();
    for (auto& p : m_tpdf) p.zero();

    // Decide if it is time to quit
    evaluateTime();
  }
}

void
Distributor::outUniPDF()
//******************************************************************************
// Output univariate PDFs to file(s)
//! \author  J. Bakosi
//******************************************************************************
{
  std::size_t i = 0;
  for (const auto& p : m_updf) {

    // Get PDF metadata
    const auto info = g_inputdeck.pdf< 1 >( i++ );

    // Construct PDF file name: base name + '_' + pdf name
    std::string filename =
      g_inputdeck.get< tag::cmd, tag::io, tag::pdf >() + '_' + info.name;

    // Augment PDF filename by time stamp if PDF output file policy is multiple
    if (g_inputdeck.get< tag::selected, tag::pdfpolicy >() ==
        ctr::PDFPolicyType::MULTIPLE)
      filename += '_' + std::to_string( m_t );

    // Augment PDF filename by '.txt' extension
    filename += ".txt";

    // Create new PDF file (overwrite if exists)
    PDFWriter pdfw( filename,
                    g_inputdeck.get< tag::selected, tag::float_format >(),
                    g_inputdeck.get< tag::discr, tag::precision >() );

    // Output PDF
    pdfw.writeTxt( p, info );
  }
}

void
Distributor::outBiPDF()
//******************************************************************************
// Output bivariate PDFs to file(s)
//! \author  J. Bakosi
//******************************************************************************
{
  std::size_t i = 0;
  for (const auto& p : m_bpdf) {

    // Get PDF metadata
    const auto info = g_inputdeck.pdf< 2 >( i++ );

    // Construct PDF file name: base name + '_' + pdf name
    std::string filename =
      g_inputdeck.get< tag::cmd, tag::io, tag::pdf >() + '_' + info.name;

    // Augment PDF filename by time stamp if PDF output file policy is multiple
    if (g_inputdeck.get< tag::selected, tag::pdfpolicy >() ==
        ctr::PDFPolicyType::MULTIPLE)
      filename += '_' + std::to_string( m_t );

    const auto& filetype = g_inputdeck.get< tag::selected, tag::pdffiletype >();

    // Augment PDF filename by the appropriate extension
    if (filetype == ctr::PDFFileType::TXT)
      filename += ".txt";
    else if (filetype == ctr::PDFFileType::GMSHTXT ||
             filetype == ctr::PDFFileType::GMSHBIN )
      filename += ".gmsh";
    else if (filetype == ctr::PDFFileType::EXODUSII)
      filename += ".exo";
    else Throw( "Unkown PDF file type attempting to output bivariate PDF" );

    // Create new PDF file (overwrite if exists)
    PDFWriter pdfw( filename,
                    g_inputdeck.get< tag::selected, tag::float_format >(),
                    g_inputdeck.get< tag::discr, tag::precision >() );

    // Output PDF
    if (filetype == ctr::PDFFileType::TXT)
      pdfw.writeTxt( p, info );
    else if (filetype == ctr::PDFFileType::GMSHTXT)
      pdfw.writeGmshTxt( p, info,
                         g_inputdeck.get< tag::selected, tag::pdfctr >() );
    else if (filetype == ctr::PDFFileType::GMSHBIN)
      pdfw.writeGmshBin( p, info,
                         g_inputdeck.get< tag::selected, tag::pdfctr >() );
    else if (filetype == ctr::PDFFileType::EXODUSII)
      pdfw.writeExodusII( p, info, static_cast< int >( m_it ),
                          g_inputdeck.get< tag::selected, tag::pdfctr >() );
  }
}

void
Distributor::outTriPDF()
//******************************************************************************
// Output trivariate PDFs to file(s)
//! \author  J. Bakosi
//******************************************************************************
{
  std::size_t i = 0;
  for (const auto& p : m_tpdf) {

    // Get PDF metadata
    const auto info = g_inputdeck.pdf< 3 >( i++ );

    // Construct PDF file name: base name + '_' + pdf name
    std::string filename =
      g_inputdeck.get< tag::cmd, tag::io, tag::pdf >() + '_' + info.name;

    // Augment PDF filename by time stamp if PDF output file policy is multiple
    if (g_inputdeck.get< tag::selected, tag::pdfpolicy >() ==
        ctr::PDFPolicyType::MULTIPLE)
      filename += '_' + std::to_string( m_t );

    const auto& filetype = g_inputdeck.get< tag::selected, tag::pdffiletype >();

    // Augment PDF filename by the appropriate extension
    if (filetype == ctr::PDFFileType::TXT)
      filename += ".txt";
    else if (filetype == ctr::PDFFileType::GMSHTXT ||
             filetype == ctr::PDFFileType::GMSHBIN )
      filename += ".gmsh";
    else if (filetype == ctr::PDFFileType::EXODUSII)
      filename += ".exo";
    else Throw( "Unkown PDF file type attempting to output trivariate PDF" );

    // Create new PDF file (overwrite if exists)
    PDFWriter pdfw( filename,
                    g_inputdeck.get< tag::selected, tag::float_format >(),
                    g_inputdeck.get< tag::discr, tag::precision >() );

    // Output PDF
    if (filetype == ctr::PDFFileType::TXT)
      pdfw.writeTxt( p, info );
    else if (filetype == ctr::PDFFileType::GMSHTXT)
       pdfw.writeGmshTxt( p, info,
                          g_inputdeck.get< tag::selected, tag::pdfctr >() );
    else if (filetype == ctr::PDFFileType::GMSHBIN)
       pdfw.writeGmshBin( p, info,
                          g_inputdeck.get< tag::selected, tag::pdfctr >() );
    else if (filetype == ctr::PDFFileType::EXODUSII)
      pdfw.writeExodusII( p, info, static_cast< int >( m_it ),
                          g_inputdeck.get< tag::selected, tag::pdfctr >() );
  }
}

void
Distributor::evaluateTime()
//******************************************************************************
// Evaluate time, decide if it is time to quit
//! \author  J. Bakosi
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
    // Continue with next time step with all integrators
    for (auto& p : m_proxy) p.advance( dt );
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
//! \author  J. Bakosi
//******************************************************************************
{
  m_print.inthead( "Time integration",
    ctr::MonteCarlo().name(g_inputdeck.get< tag::selected, tag::montecarlo >()),
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
//! \author  J. Bakosi
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
