//******************************************************************************
/*!
  \file      src/Inciter/Conductor.h
  \author    J. Bakosi
  \date      Fri 15 May 2015 03:01:35 PM MDT
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

    //! Forward a vector of time stamps to the main proxy
    void timestamp(
      const std::vector< std::pair< std::string, tk::real > >& stamp );

  private:
    using PerfProxy = CProxy_Performer;
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Conductor >;

    InciterPrint m_print;               //!< Pretty printer
    std::vector< tk::Timer > m_timer;   //!< Timers
    int m_nchare;                       //!< Number of performer chares
    int m_charecnt;                     //!< Chare counter
    PerfProxy m_perfproxy;              //!< Performer chare array
    LinSysMergerProxy m_lsmproxy;       //!< Linear system merger chare group
    std::map< std::string, std::vector< tk::real > > m_avg; //!< Avg time stamps
};

} // inciter::

#endif // Conductor_h
