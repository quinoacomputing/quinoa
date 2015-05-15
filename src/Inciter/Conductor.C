//******************************************************************************
/*!
  \file      src/Inciter/Conductor.C
  \author    J. Bakosi
  \date      Fri 15 May 2015 12:35:18 PM MDT
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

extern CProxy_Main mainProxy;

namespace inciter {

extern std::size_t g_npoin;
extern std::vector< std::map< std::size_t, std::vector<std::size_t> > > g_pcomm;

} // inciter::

using inciter::Conductor;

Conductor::Conductor() :
  m_print( g_inputdeck.get<tag::cmd,tag::verbose>() ? std::cout : std::clog ),
  m_timer( 1 ), // start timer
  m_nchare( static_cast< int >( g_pcomm.empty() ? 1 : g_pcomm.size() ) ),
  m_charecnt( 0 )
//******************************************************************************
// Constructor
//! \author J. Bakosi
//******************************************************************************
{
  // Create linear system merger chare group collecting chear contributions
  m_lsmproxy = LinSysMergerProxy::ckNew( thisProxy, g_npoin );

  // Fire up array of asynchronous performers
  m_perfproxy = PerfProxy::ckNew( thisProxy, m_lsmproxy, m_nchare );
}

void
Conductor::init() const
//******************************************************************************
// Reduction target indicating that all members of LinSysMerger have finished
// their portion of initializing the linear system distributed across all PEs
//! \author J. Bakosi
//******************************************************************************
{
  mainProxy.timestamp(
    "Initialize Charm++ chare arrays and distributed linear system",
    m_timer[0].dsec() );
  mainProxy.finalize();
}

void
Conductor::timestamp(
  const std::vector< std::pair< std::string, tk::real > >& stamp )
//******************************************************************************
// Forward a vector of time stamps to the main proxy
//! \param[in] stamp Vector of time stamps contributed
//! \author J. Bakosi
//******************************************************************************
{
  // Collect time stamps with the same name. These come from different chares
  // performing their portion of some work and we will display the average only.
  for (const auto& s : stamp) m_avg[ s.first ].push_back( s.second );

  // If all Performer chares have contributed, compute averages and forward
  if (++m_charecnt == m_nchare) {

    // Compute average times for time stamps with the same name
    std::vector< std::pair< std::string, tk::real > > s;
    for (const auto& t : m_avg) {
      tk::real sum = 0.0;
      for (auto v : t.second) sum += v;
      s.emplace_back(
        t.first + ", avg of " + std::to_string(m_nchare) + " chares",
        sum/t.second.size() );
    }

    // Send to main proxy for printing at the end
    mainProxy.timestamp( s );

    // Prepare for next round
    m_charecnt = 0;
    m_avg.clear();
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
