//******************************************************************************
/*!
  \file      src/Inciter/Conductor.C
  \author    J. Bakosi
  \date      Fri 10 Apr 2015 05:29:47 PM MDT
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

extern
  std::vector< std::map< std::size_t, std::vector< std::size_t > > > g_comm;

} // inciter::

using inciter::Conductor;

Conductor::Conductor() :
  m_print( g_inputdeck.get<tag::cmd,tag::verbose>() ? std::cout : std::clog ),
  m_count( 0 ),
  m_timer( 1 )  // start timer
//******************************************************************************
// Constructor
//! \author J. Bakosi
//******************************************************************************
{
  // Get number of chares
  auto nchare = static_cast< int >( g_comm.empty() ? 1 : g_comm.size() );

  // Create linear system merger chare group collecting chear contributions
  m_lsmproxy = LinSysMergerProxy::ckNew( thisProxy );

  // Fire up array of asynchronous performers
  m_perfproxy = PerfProxy::ckNew( thisProxy, m_lsmproxy, nchare );
}

void
Conductor::init()
//******************************************************************************
// Wait for all members of LinSysMerger to create their portion of the linear
// system distributed across all PEs
//! \author J. Bakosi
//******************************************************************************
{
// Increase number of performers completing initialization
  ++m_count.get< tag::init >();

  // Wait for all performers completing initialization
  if (m_count.get< tag::init >() == CkNumPes()) {
    mainProxy.timestamp(
      "Initialize Charm++ chare arrays and distributed linear system",
      m_timer[0].dsec() );
    mainProxy.finalize();
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
