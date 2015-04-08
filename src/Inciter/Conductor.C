//******************************************************************************
/*!
  \file      src/Inciter/Conductor.C
  \author    J. Bakosi
  \date      Wed 08 Apr 2015 08:32:33 AM MDT
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
#include <Performer.h>
#include <inciter.decl.h>

extern CProxy_Main mainProxy;

using inciter::Conductor;

Conductor::Conductor() :
  m_print( g_inputdeck.get< tag::cmd, tag::verbose >() ? std::cout : std::clog )
//******************************************************************************
// Constructor
//! \author J. Bakosi
//******************************************************************************
{
  // Start timer measuring total integration time
  m_timer.emplace_back();

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
