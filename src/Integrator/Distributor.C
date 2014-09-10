//******************************************************************************
/*!
  \file      src/Integrator/Distributor.C
  \author    J. Bakosi
  \date      Tue 09 Sep 2014 03:51:21 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Distributor drives the time integration of differential equations
  \details   Distributor drives the time integration of differential equations
*/
//******************************************************************************

#include <Distributor.h>
#include <Integrator.h>
#include <DiffEqStack.h>
#include <quinoa.decl.h>

extern CProxy_Main mainProxy;

using quinoa::Distributor;

Distributor::Distributor( const ctr::CmdLine& cmdline ) :
  m_print( cmdline.get< tk::tag::verbose >() ? std::cout : std::clog ),
  m_ninit( 0 ),
  m_nAccOrd( 0 ),
  m_nAccCen( 0 ),
  m_numchares( 0 ),
  m_it( 0 ),
  m_t( 0.0 ),
  m_plotOrdinary( g_inputdeck.plotOrdinary() ),
  m_nameOrdinary( g_inputdeck.momentNames( ctr::ordinary ) ),
  m_nameCentral( g_inputdeck.momentNames( ctr::central ) ),
  m_ordinary( m_nameOrdinary.size(), 0.0 ),
  m_central( m_nameCentral.size(), 0.0 ),
  m_statWriter( g_inputdeck.get< tag::cmd, tag::io, tag::stat >() ),
  m_output( false, false )
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
  m_print.statistics( "Statistics and probabilities" );

  // Print out info on load distirubtion
  m_print.section( "Load distribution" );
  m_print.item( "Virtualization [0.0...1.0]",
                g_inputdeck.get< tag::cmd, tag::virtualization >() );
  m_print.item( "Load (number of particles)",
                g_inputdeck.get< tag::discr, tag::npar >() );
  m_print.item( "Number of processing elements", CkNumPes() );
  m_print.item( "Number of work units",
                std::to_string( m_numchares ) + " (" +
                std::to_string( m_numchares-1 ) + "*" +
                std::to_string( chunksize ) + "+" +
                std::to_string( chunksize+remainder ) + ")" );

  // Print out time integration header
  if (g_inputdeck.get< tag::discr, tag::nstep >()) {
    header();
    m_statWriter.header( m_plotOrdinary, m_nameOrdinary, m_nameCentral );
  }

  // Start timer measuring total integration time
  m_timer.emplace_back();

  // Compute size of initial time step
  const auto dt = computedt();

  // Fire up integrators
  for (auto i=decltype(m_numchares){1}; i<m_numchares; ++i)
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
  m_numchares = npar/chunksize;

  // Compute remainder of work if the above number of units were to be created
  remainder = npar - m_numchares*chunksize;

  // Redistribute remainder among the work units for a more equal distribution
  chunksize += remainder/m_numchares;

  // Compute new remainder (after redistribution of the previous remainder)
  remainder = npar - m_numchares*chunksize;
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
  ++m_ninit;

  // Wait for all integrators completing initialization
  if ( m_ninit == m_numchares ) {
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
  ++m_nAccOrd;

  // Add contribution from PE to total sums, i.e., u[i] += v[i] for all i
  for (std::size_t i=0; i<m_ordinary.size(); ++i) m_ordinary[i] += ord[i];

  // Wait for all integrators completing accumulation of ordinary moments
  if ( m_nAccOrd == m_numchares ) {
    // Finish computing moments, i.e., divide sums by the number of samples
    for (auto& m : m_ordinary) m /= g_inputdeck.get< tag::discr, tag::npar >();
    // Continue with accumulation for central moments with all integrators
    for (auto& p : m_proxy) p.estimateCen( m_ordinary );
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
  ++m_nAccCen;

  // Add contribution from PE to total sums, i.e., u[i] += v[i] for all i
  for (std::size_t i=0; i<m_central.size(); ++i) m_central[i] += cen[i];

  // Wait for all integrators completing accumulation of central moments
  if ( m_nAccCen == m_numchares ) {

    // Finish computing moments, i.e., divide sums by the number of samples
    for (auto& m : m_central) m /= g_inputdeck.get< tag::discr, tag::npar >();

     // Append statistics file at selected times
    if (!(m_it % g_inputdeck.get< tag::interval, tag::stat >())) {
      m_statWriter.writeStat( m_it, m_t, m_ordinary, m_central, m_plotOrdinary );
      m_output.get< tag::stat >() = true;
    }

    // Zero accumulator counters and total-sums for next time step
    m_nAccOrd = m_nAccCen = 0;
    std::fill( begin(m_ordinary), end(m_ordinary), 0.0 );    
    std::fill( begin(m_central), end(m_central), 0.0 );    

    // Decide if it is time to quit
    evaluateTime();
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
