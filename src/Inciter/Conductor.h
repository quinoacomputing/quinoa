//******************************************************************************
/*!
  \file      src/Inciter/Conductor.h
  \author    J. Bakosi
  \date      Wed 08 Apr 2015 08:32:00 AM MDT
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
#ifndef Conductor_h
#define Conductor_h

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <conductor.decl.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <InciterPrint.h>

namespace inciter {

//! Conductor drives the time integration of the Euler equations
class Conductor : public CBase_Conductor {

  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  //Conductor_SDAG_CODE

  public:
    //! Constructor
    explicit Conductor();

    //! Finish initialization
    void init();

  private:
    using CProxyPerf = CProxy_Performer< CProxy_Conductor >;

    //! Pretty printer
    InciterPrint m_print;

    //! Counters of performer chares completing a function
    tk::tuple::tagged_tuple< tag::init,  int,
                             tag::chare, int > m_count;

    std::vector< tk::Timer > m_timer;           //!< Timers
};

} // inciter::

#endif // Conductor_h
