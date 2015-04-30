//******************************************************************************
/*!
  \file      src/Inciter/Conductor.h
  \author    J. Bakosi
  \date      Wed 29 Apr 2015 11:49:15 AM MDT
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

  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  //Conductor_SDAG_CODE

  public:
    //! Constructor
    explicit Conductor();

    //! \brief Reduction target indicating that all Performer chares have
    //!   registered with the linear system merger, LinSysMerger
    //! \details This function is a Charm++ reduction target that is called when
    //!   all Performer chares have registered with their their local branch of
    //!   the linear system merger group, LinSysMerger. Once this is done, we
    //!   issue a broadcast to all Performer chares to initialize their portion
    //!   of the linear system.
    void registered() { m_perfproxy.initLinearSystem(); }

    //! \brief Reduction target indicating that all Performer chares have
    //!   submitted their contribution toward initialization of the linear
    //!   system distributed across all PEs
    //! \details This function is a Charm++ reduction target that is called when
    //!   all Performer chares have submitted their contribution toward
    //!   initialization of the linear system distributed across all PEs by
    //!   submitting their contribution to their local branch of the linear
    //!   system merger group, LinSysMerger. Once this is done, we issue a
    //!   broadcast to all members of LinSysMerger (one per PE) to create the
    //!   linear system.
    //! \author J. Bakosi
    void linsysinit() { m_lsmproxy.createLinearSystem(); }

    //! \brief Wait for all members of LinSysMerger to create their portion of
    //!   the linear system distrubted across all PEs
    void init();

  private:
    using PerfProxy = CProxy_Performer;
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Conductor >;

    //! Pretty printer
    InciterPrint m_print;

    //! Counters of performer chares completing a function
    tk::tuple::tagged_tuple< tag::init, int > m_count;

    //! Timers
    std::vector< tk::Timer > m_timer;

    //! Charm++ proxy to array of Performer chares
    PerfProxy m_perfproxy;

    //! Charm++ proxy to group of linear system merger chares
    LinSysMergerProxy m_lsmproxy;
};

} // inciter::

#endif // Conductor_h
