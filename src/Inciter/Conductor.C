//******************************************************************************
/*!
  \file      src/Inciter/Conductor.C
  \author    J. Bakosi
  \date      Wed 21 Oct 2015 07:58:24 AM MDT
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
#include "ContainerUtil.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "inciter.decl.h"

extern CProxy_Main mainProxy;

namespace inciter {

extern std::size_t g_npoin;
extern std::vector< std::map< std::size_t, std::vector<std::size_t> > > g_pcomm;
extern std::vector< std::vector< std::size_t > > g_point;

} // inciter::

using inciter::Conductor;

Conductor::Conductor() :
  m_print( g_inputdeck.get<tag::cmd,tag::verbose>() ? std::cout : std::clog ),
  m_timer( 1 ), // start a timer
  m_nchare( static_cast< int >( g_pcomm.empty() ? 1 : g_pcomm.size() ) ),
  m_it( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_dt( computedt() ),
  m_stage( 0 ),
  m_arrTimestampCnt( 0 ),
  m_grpTimestampCnt( 0 ),
  m_arrPerfstatCnt( 0 ),
  m_grpPerfstatCnt( 0 )
//******************************************************************************
// Constructor
//! \author J. Bakosi
//******************************************************************************
{
  // Print out info and time stepping header
  info();

  // Create linear system merger chare group collecting chare contributions
  m_lsmproxy = LinSysMergerProxy::ckNew( thisProxy, g_npoin );

  // Charm++ map object for custom initial placement of chare array elements
  auto map = tk::CProxy_UnsMeshMap::ckNew( g_npoin, g_point );
  //auto map = tk::CProxy_LinearMap::ckNew( m_nchare );
  //auto map = CProxy_RRMap::ckNew();      // round robin
  CkArrayOptions opts( m_nchare );
  opts.setMap( map );

  // Fire up array of asynchronous performers
  m_perfproxy = PerformerProxy::ckNew( thisProxy, m_lsmproxy, opts );
  //m_perfproxy = PerformerProxy::ckNew( thisProxy, m_lsmproxy, m_nchare );
  m_perfproxy.doneInserting();
}

void
Conductor::info() const
//******************************************************************************
//  Print information at startup
//! \author J. Bakosi
//******************************************************************************
{
  // Print out information on problem
  m_print.part( "Problem" );

  // Print out info on problem title
  if ( !g_inputdeck.get< tag::title >().empty() )
    m_print.title( g_inputdeck.get< tag::title >() );

  // Print I/O filenames
  m_print.section( "Output filenames" );
  m_print.item( "Field",
              g_inputdeck.get< tag::cmd, tag::io, tag::output >() + "_<PEid>" );

  // Print discretization parameters
  m_print.section( "Discretization parameters" );
  m_print.item( "Number of time steps",
                g_inputdeck.get< tag::discr, tag::nstep >() );
  m_print.item( "Start time",
                g_inputdeck.get< tag::discr, tag::t0 >() );
  m_print.item( "Terminate time",
                g_inputdeck.get< tag::discr, tag::term >() );
  m_print.item( "Initial time step size",
                g_inputdeck.get< tag::discr, tag::dt >() );

  // Print output intervals
  m_print.section( "Output intervals" );
  m_print.item( "TTY", g_inputdeck.get< tag::interval, tag::tty>() );
  m_print.item( "Field", g_inputdeck.get< tag::interval, tag::field >() );

  // Print out time integration header
  if (g_inputdeck.get< tag::discr, tag::nstep >()) header();
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
Conductor::evaluateTime()
//******************************************************************************
// Evaluate time step, compute new time step size, decide if it is time to quit
//! \author J. Bakosi
//******************************************************************************
{
  // Update stage in multi-stage time stepping
  if (m_stage < 1) {    // if not final stage, continue with next stage

    m_perfproxy.advance( ++m_stage, m_dt, m_it, m_t );

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

    // Finish if either max iterations or max time reached
    if ( std::fabs(m_t - term) > std::numeric_limits< tk::real >::epsilon() &&
         m_it < g_inputdeck.get< tag::discr, tag::nstep >() ) {

      // Continue with next time step (at stage 0) with all integrators
      m_stage = 0;
      m_perfproxy.advance( m_stage, m_dt, m_it, m_t );

    } else {

      finish();

    }
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
Conductor::header() const
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
}

void
Conductor::report()
//******************************************************************************
// Print out one-liner report on time step
//! \author J. Bakosi
//******************************************************************************
{
  if (!(m_it % g_inputdeck.get< tag::interval, tag::tty >())) {

    // estimated time elapsed and for accomplishment
    tk::Timer::Watch ete, eta;
    m_timer[0].eta( g_inputdeck.get< tag::discr, tag::term >() -
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
