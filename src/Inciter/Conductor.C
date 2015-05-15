//******************************************************************************
/*!
  \file      src/Inciter/Conductor.C
  \author    J. Bakosi
  \date      Wed 13 May 2015 06:12:42 AM MDT
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
  m_timer( 1 )  // start timer
//******************************************************************************
// Constructor
//! \author J. Bakosi
//******************************************************************************
{
  // Get number of chares
  auto nchare = static_cast< int >( g_pcomm.empty() ? 1 : g_pcomm.size() );

  // Create linear system merger chare group collecting chear contributions
  m_lsmproxy = LinSysMergerProxy::ckNew( thisProxy, g_npoin );

  // Fire up array of asynchronous performers
  m_perfproxy = PerfProxy::ckNew( thisProxy, m_lsmproxy, nchare );
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

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <conductor.def.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
