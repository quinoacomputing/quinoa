//******************************************************************************
/*!
  \file      src/Inciter/Conductor.C
  \author    J. Bakosi
  \date      Wed 20 May 2015 01:51:36 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Conductor drives the time integration of the Euler equations
  \details   Conductor drives the time integration of the Euler equations.
    The implementation uses the Charm++ runtime system and is fully asynchronous,
    overlapping computation, communication as well as I/O. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/conductor.ci.
*/
//******************************************************************************

#include <Conductor.h>
#include <inciter.decl.h>
#include <LinearMap.h>
#include <UnsMeshMap.h>

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
  m_arrTimestampCnt( 0 ),
  m_grpTimestampCnt( 0 ),
  m_arrPerfstatCnt( 0 ),
  m_grpPerfstatCnt( 0 )
//******************************************************************************
// Constructor
//! \author J. Bakosi
//******************************************************************************
{
  // Create linear system merger chare group collecting chear contributions
  m_lsmproxy = LinSysMergerProxy::ckNew( thisProxy, g_npoin );

  // Charm++ map object for custom initial placement of chare array elements
  auto map = tk::CProxy_UnsMeshMap::ckNew( g_npoin, g_point );
  //auto map = tk::CProxy_LinearMap::ckNew( m_nchare );
  //auto map = CProxy_RRMap::ckNew();      // round robin
  CkArrayOptions opts( m_nchare );
  opts.setMap( map );

  // Fire up array of asynchronous performers
  m_perfproxy = PerfProxy::ckNew( thisProxy, m_lsmproxy, opts );
  //m_perfproxy = PerfProxy::ckNew( thisProxy, m_lsmproxy, m_nchare );
  m_perfproxy.doneInserting();
}

void
Conductor::init() const
//******************************************************************************
// Reduction target indicating that all members of LinSysMerger have finished
// their portion of initializing the linear system distributed across all PEs
//! \author J. Bakosi
//******************************************************************************
{
  mainProxy.timestamp( "Initialize distributed linear system",
                       m_timer[0].dsec() );
  mainProxy.finalize();
}

void
Conductor::arrTimestamp(
  const std::vector< std::pair< std::string, tk::real > >& stamp )
//******************************************************************************
// Forward a vector of time stamps from (Performer) chare array to the main
// proxy
//! \param[in] stamp Vector of time stamps contributed
//! \author J. Bakosi
//******************************************************************************
{
  timestamp( stamp, m_arrTimestamp, m_arrTimestampCnt, m_nchare, "chares" );
  m_lsmproxy.trigger_timestamp_complete();
}

void
Conductor::grpTimestamp(
  const std::vector< std::pair< std::string, tk::real > >& stamp )
//******************************************************************************
// Forward a vector of time stamps from (LinSysMerger) chare group branches to
// the main proxy
//! \param[in] stamp Vector of time stamps contributed
//! \author J. Bakosi
//******************************************************************************
{
  timestamp( stamp, m_grpTimestamp, m_grpTimestampCnt, CkNumPes(), "PEs" );
  m_lsmproxy.trigger_perfstat_complete();
}

void
Conductor::arrPerfstat(
  const std::vector< std::pair< std::string, tk::real > >& p )
//******************************************************************************
// Forward a vector of performance statistics from (Performer) chare array
// elements to the main proxy
//! \param[in] p Vector of performance statistics contributed
//! \author J. Bakosi
//******************************************************************************
{
  perfstat( p, m_arrPerfstat, m_arrPerfstatCnt, m_nchare, "chares" );
}

void
Conductor::grpPerfstat(
  const std::vector< std::pair< std::string, tk::real > >& p )
//******************************************************************************
// Forward a vector of performance statistics from (LinSysMerger) chare group
// branches to the main proxy
//! \param[in] p Vector of performance statistics contributed
//! \author J. Bakosi
//******************************************************************************
{
  perfstat( p, m_grpPerfstat, m_grpPerfstatCnt, CkNumPes(), "PEs" );
}

void
Conductor::timestamp(
  const std::vector< std::pair< std::string, tk::real > >& stamp,
  std::map< std::string, std::vector< tk::real > >& map,
  int& counter,
  int max,
  const std::string& of )
//******************************************************************************
// Collect and compute averages over time stamps contributed by chares
//! \param[in] stamp Vector of time stamps contributed
//! \param[inout] map Map in which to sore labels and vector of times
//! \param[inout] counter Counter to use
//! \param[in] max Max number of contributions to expect
//! \param[in] of Of what the averages computed over
//! \author J. Bakosi
//******************************************************************************
{
  // Collect time stamps with the same name. These come from different chares
  // performing their portion of some work.
  for (const auto& s : stamp) map[ s.first ].push_back( s.second );

  if (++counter == max) {       // if all expected have arrived
    // Compute average over chares and send to main proxy for printing
    mainProxy.timestamp(
      tk::average( map, ", avg of " + std::to_string(max) + " " + of ) );
    // Prepare for next round
    counter = 0;
    map.clear();
  }
}

void
Conductor::perfstat( const std::vector< std::pair< std::string, tk::real > >& p,
                     std::map< std::string, std::vector< tk::real > >& map,
                     int& counter,
                     int max,
                     const std::string& of )
//******************************************************************************
// Collect and compute performance statistics contributed by chares
//! \param[in] p Vector of time stamps contributed
//! \param[inout] map Map in which to sore labels and vector of times
//! \param[inout] counter Counter to use
//! \param[in] max Max number of contributions to expect
//! \param[in] of Of what the statistics computed over
//! \author J. Bakosi
//******************************************************************************
{
  // Collect performance statistics with the same name. These come from
  // different chares performing their portion of some work.
  for (const auto& s : p) map[ s.first ].push_back( s.second );

  if (++counter == max) {       // if all expected have arrived
    const auto nc = std::to_string( max );
    // Compute average over chares and send to main proxy for printing
    const auto avg = tk::average( map, ", avg of " + nc + " " + of );
    mainProxy.perfstat( avg );
    // Compute variance over chares and send to main proxy for printing
    mainProxy.perfstat(
      tk::variance( map, avg, ", var of " + nc + " " + of ) );
    // Prepare for next round
    counter = 0;
    map.clear();
  }
}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <conductor.def.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
