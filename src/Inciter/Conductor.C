//******************************************************************************
/*!
  \file      src/Inciter/Conductor.C
  \author    J. Bakosi
  \date      Mon 23 Nov 2015 02:28:55 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Conductor drives the time integration of a PDE
  \details   Conductor drives the time integration of a PDE
    The implementation uses the Charm++ runtime system and is fully asynchronous,
    overlapping computation, communication as well as I/O. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/conductor.ci.
*/
//******************************************************************************

#include <string>
#include <iostream>
#include <cstddef>

#include "Conductor.h"
#include "Spawner.h"
#include "ContainerUtil.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "inciter.decl.h"

extern CProxy_Main mainProxy;

using inciter::Conductor;

Conductor::Conductor( std::size_t npoin, uint64_t nchare ) :
  m_print( g_inputdeck.get<tag::cmd,tag::verbose>() ? std::cout : std::clog ),
  m_nchare( static_cast< int >( nchare ) ),
  m_eval( 0 ),
  m_init( 0 ),
  m_it( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_dt( computedt() ),
  m_stage( 0 ),
  m_arrTimestampCnt( 0 ),
  m_grpTimestampCnt( 0 ),
  m_arrPerfstatCnt( 0 ),
  m_grpPerfstatCnt( 0 )
//******************************************************************************
//  Constructor
//! \param[in] npoin Total number of points in computational mesh
//! \param[in] nchare Total number of chares
//! \author J. Bakosi
//******************************************************************************
{
  // Print out information on problem
  m_print.part( "Problem" );

  // Print out info on problem title
  if ( !g_inputdeck.get< tag::title >().empty() )
    m_print.title( g_inputdeck.get< tag::title >() );

  // Print discretization parameters
  m_print.section( "Discretization parameters" );
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto dt = g_inputdeck.get< tag::discr, tag::dt >();
  m_print.item( "Number of time steps", nstep );
  m_print.item( "Start time", t0 );
  m_print.item( "Terminate time", term );
  m_print.item( "Initial time step size", dt );
  // If the desired max number of time steps is larger than zero, and the
  // termination time is larger than the initial time, and the initial time step
  // is smaller than the duration of the time to be simulated, we have work to
  // do, otherwise, finish right away
  if ( nstep != 0 && term > t0 && dt < term-t0 ) {

   // Print I/O filenames
    m_print.section( "Output filenames" );
    m_print.item( "Field",
      g_inputdeck.get< tag::cmd, tag::io, tag::output >() + "_<PEid>" );

    // Print output intervals
    m_print.section( "Output intervals" );
    m_print.item( "TTY", g_inputdeck.get< tag::interval, tag::tty>() );
    m_print.item( "Field", g_inputdeck.get< tag::interval, tag::field >() );
    m_print.endsubsection();

    // Create linear system merger chare group collecting chare contributions to
    // a linear system
    m_linsysmerger = LinSysMergerProxy::ckNew( thisProxy, npoin );

    m_print.diagstart( "Creating workers ..." );

    m_timer[ TimerTag::CREATE ].zero();

    // Create chare group spawning asynchronous performers
    m_spawner = SpawnerProxy::ckNew( m_nchare, thisProxy );
    for (int p=0; p<CkNumPes(); ++p) m_spawner[ p ].create( m_linsysmerger );

  } else finish();      // stop if no time stepping requested
}

void
Conductor::timestamp(
  const std::vector< std::pair< std::string, tk::real > >& stamp,
  std::map< std::string, std::vector< tk::real > >& map,
  int& counter,
  int max )
//******************************************************************************
// Collect and compute averages over time stamps contributed by chares
//! \param[in] stamp Vector of time stamps contributed
//! \param[inout] map Map in which to sore labels and vector of times
//! \param[inout] counter Counter to use
//! \param[in] max Max number of contributions to expect
//! \author J. Bakosi
//******************************************************************************
{
  // Collect time stamps with the same name. These come from different chares
  // performing their portion of some work.
  for (const auto& s : stamp) map[ s.first ].push_back( s.second );

  // If all expected have arrived, get ready for the next round
  if (++counter == max) counter = 0;
}

void
Conductor::perfstat( const std::vector< std::pair< std::string, tk::real > >& p,
                     std::map< std::string, std::vector< tk::real > >& map,
                     int& counter,
                     int max )
//******************************************************************************
// Collect and compute performance statistics contributed by chares
//! \param[in] p Vector of time stamps contributed
//! \param[inout] map Map in which to sore labels and vector of times
//! \param[inout] counter Counter to use
//! \param[in] max Max number of contributions to expect
//! \author J. Bakosi
//******************************************************************************
{
  // Collect performance statistics with the same name. These come from
  // different chares performing their portion of some work.
  for (const auto& s : p) map[ s.first ].push_back( s.second );

  // If all expected have arrived, get ready for the next round
  if (++counter == max) counter = 0;
}

void
Conductor::finalReport()
//******************************************************************************
// Send collected timer and performance data to host
//! \author J. Bakosi
//******************************************************************************
{
  // Compute average over chare array and send to main proxy for printing
  mainProxy.timestamp(
    tk::average( m_arrTimestamp,
                 ", avg of " + std::to_string(m_nchare) + " chares" ) );

  // Compute average over chare group and send to main proxy for printing
  mainProxy.timestamp(
    tk::average( m_grpTimestamp,
                 ", avg of " + std::to_string(CkNumPes()) + " PEs" ) );

  // Lambda for sending performance statistics to host
  auto perf = []( std::map< std::string, std::vector< tk::real > >& map,
                  int max,
                  const std::string& of )
  {
    const auto nc = std::to_string( max );
    // Compute average over chares and send to main proxy for printing
    const auto avg = tk::average( map, ", avg of " + nc + " " + of );
    mainProxy.perfstat( avg );
    // Compute variance over chares and send to main proxy for printing
    mainProxy.perfstat(
      tk::variance( map, avg, ", var of " + nc + " " + of ) );
  };

  // Send performance statistics of chare array to host
  perf( m_arrPerfstat, m_nchare, "chares" );
  perf( m_grpPerfstat, CkNumPes(), "PEs" );

  // Quit
  mainProxy.finalize();
}

void
Conductor::created()
//******************************************************************************
// Reduction target indicating that all Spawner chare groups have finished
// creating their Charm++ Performer worker chare array elements (initializing
// their mesh element ids they will work on)
//! \details Once this is done on all PE, indicated by entering this
//!   function, we issue a broadcast to all Spawners (on PEs) to issue their
//!   workers to start setup.
//! \author J. Bakosi
//******************************************************************************
{
  mainProxy.timestamp( "Create " + std::to_string(m_nchare) + " workers on " +
                       std::to_string(CkNumPes()) + " PEs",
                       tk::query(m_timer, TimerTag::CREATE) );
  m_print.diagend( "done" );
  m_print.diagstart( "Setting up workers ..." );
  m_timer[ TimerTag::SETUP ].zero();
  m_spawner.setup();
}

void
Conductor::rowcomplete()
//******************************************************************************
// Reduction target indicating that all linear system merger branches have done
// their part of storing and exporting global row ids
//! \details This function is a Charm++ reduction target that is called when
//!   all linear system merger branches have done their part of storing and
//!   exporting global row ids. Once this is done, we issue a broadcast to
//!   all Spawners and thus implicitly all Performer chares to continue with
//!   the initialization step.
//! \author J. Bakosi
//******************************************************************************
{
  mainProxy.timestamp( "Setup " + std::to_string(m_nchare) + " workers on " +
                       std::to_string(CkNumPes()) + " PEs",
                       tk::query(m_timer, TimerTag::SETUP) );
  m_linsysmerger.rowsreceived();
  m_print.diagend( "done" );
  m_print.diagstart( "Initializing workers ..." );
  m_timer[ TimerTag::INITIALIZE ].zero();
  m_spawner.init( m_dt );
}

void
Conductor::initcomplete()
//******************************************************************************
//  Reduction target indicating that all Performer chares have finished their
//  initialization step and have already continued with start time stepping
//! \author J. Bakosi
//******************************************************************************
{
  if ( ++m_init == CkNumPes() ) {
    mainProxy.timestamp( "Initialize " + std::to_string(m_nchare) + " workers "
                         "on " + std::to_string(CkNumPes()) + " PEs",
                         tk::query(m_timer, TimerTag::INITIALIZE) );
    m_init = 0; // get ready for next time
    m_print.diagend( "done" );
    m_print.diagstart( "Starting time stepping ...\n" );
    header();   // print out time integration header
  }
}

void
Conductor::evaluateTime()
//******************************************************************************
//  Evaluate time step: decide if it is time to quit
//! \author J. Bakosi
//******************************************************************************
{
  if ( ++m_eval == CkNumPes() ) {
    m_eval = 0; // get ready for next time
    // Get physical time at which to terminate
    const auto term = g_inputdeck.get< tag::discr, tag::term >();
    const auto eps = std::numeric_limits< tk::real >::epsilon();
    const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
    // if not final stage of time step or if neither max iterations nor max time
    // reached, will continue by telling all linear system merger group elements
    // to prepare for a new rhs, otherwise, finish
    if ( m_stage < 1 || (std::fabs(m_t-term) > eps && m_it < nstep) )
      m_linsysmerger.enable_wait4rhs();
    else
      finish();
  }
}

void
Conductor::advance()
//******************************************************************************
//  Reduction target indicating that all Performer chares have finished their
//  initialization step
//! \author J. Bakosi
//******************************************************************************
{
  // Update stage in multi-stage time stepping
  if (m_stage < 1) {    // if not final stage, continue with next stage

    m_spawner.advance( ++m_stage, m_dt, m_it, m_t );

  } else {              // if final stage, evaluate time

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

    // Continue with next time step (at stage 0) with all integrators
    m_stage = 0;
    m_spawner.advance( m_stage, m_dt, m_it, m_t );

  }
}

tk::real
Conductor::computedt()
//******************************************************************************
// Compute size of next time step
//! \return Size of dt for the next time step
//! \author  J. Bakosi
//******************************************************************************
{
  // Simply return the constant user-defined initial dt for now
  return g_inputdeck.get< tag::discr, tag::dt >();
}

void
Conductor::finish()
//******************************************************************************
// Normal finish of time stepping
//! \author J. Bakosi
//******************************************************************************
{
  // Print out reason for stopping
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  m_print.endsubsection();
  if (m_it >= g_inputdeck.get< tag::discr, tag::nstep >())
     m_print.note( "Normal finish, maximum number of iterations reached: " +
                   std::to_string( nstep ) );
   else
     m_print.note( "Normal finish, maximum time reached: " +
                   std::to_string( g_inputdeck.get<tag::discr,tag::term>() ) );

  // Send timer and performance data to main proxy and quit
  finalReport();
}

void
Conductor::header()
//******************************************************************************
// Print out time integration header
//! \author J. Bakosi
//******************************************************************************
{
  m_print.inthead( "Time integration", "Unstructured-mesh PDE solver testbed",
    "Legend: it - iteration count\n"
    "         t - time\n"
    "        dt - time step size\n"
    "       ETE - estimated time elapsed (h:m:s)\n"
    "       ETA - estimated time for accomplishment (h:m:s)\n"
    "       out - output-saved flags (F: field)\n",
    "\n      it             t            dt        ETE        ETA   out\n"
      " ---------------------------------------------------------------\n" );
  m_timer[ TimerTag::TIMESTEP ].zero();
}

void
Conductor::report()
//******************************************************************************
// Print out one-liner report on time step
//! \author J. Bakosi
//******************************************************************************
{
  if (!(m_it % g_inputdeck.get< tag::interval, tag::tty >())) {

    // estimate time elapsed and time for accomplishment
    tk::Timer::Watch ete, eta;
    const auto& timer = tk::cref_find( m_timer, TimerTag::TIMESTEP );
    timer.eta( g_inputdeck.get< tag::discr, tag::term >() -
                 g_inputdeck.get< tag::discr, tag::t0 >(),
               m_t -  g_inputdeck.get< tag::discr, tag::t0 >(),
               g_inputdeck.get< tag::discr, tag::nstep >(),
               m_it,
               ete,
               eta );

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
    if (!(m_it % g_inputdeck.get< tag::interval, tag::field >()))
      m_print << 'F';

    m_print << '\n';
  }
}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "conductor.def.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
