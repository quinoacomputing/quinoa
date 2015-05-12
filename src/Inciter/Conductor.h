//******************************************************************************
/*!
  \file      src/Inciter/Conductor.h
  \author    J. Bakosi
  \date      Tue 12 May 2015 09:27:56 AM MDT
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
#include <Performer.h>

namespace inciter {

//! Conductor drives the time integration of the Euler equations
class Conductor : public CBase_Conductor {

  public:
    //! Constructor
    explicit Conductor();

    //! \brief Reduction target indicating that all members of LinSysMerger have
    //!   finished their portion of initializing the linear system distributed
    //!   across all PEs
    void init() const;

  private:
    using PerfProxy = CProxy_Performer;
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Conductor >;

    //! Pretty printer
    InciterPrint m_print;

    //! Timers
    std::vector< tk::Timer > m_timer;

    //! Charm++ proxy to array of Performer chares
    PerfProxy m_perfproxy;

    //! Charm++ proxy to group of linear system merger chares
    LinSysMergerProxy m_lsmproxy;
};

} // inciter::

#endif // Conductor_h
