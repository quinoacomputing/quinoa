//******************************************************************************
/*!
  \file      src/Inciter/Conductor.h
  \author    J. Bakosi
  \date      Sun 14 Jun 2015 09:26:12 PM MDT
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

#include <map>
#include <vector>
#include <iosfwd>
#include <utility>

#include "Timer.h"
#include "Types.h"
#include "InciterPrint.h"
#include "Performer.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "conductor.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

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

    //! \brief Forward a vector of time stamps from (Performer) chare array to
    //!   the main proxy
    void arrTimestamp(
      const std::vector< std::pair< std::string, tk::real > >& stamp );

    //! \brief Forward a vector of time stamps from (LinSysMerger) chare group
    //!   branches to the main proxy
    void grpTimestamp(
      const std::vector< std::pair< std::string, tk::real > >& stamp );

    //! \brief Forward a vector of performance statistics from (Performer) chare
    //!   array elements to the main proxy
    void arrPerfstat(
      const std::vector< std::pair< std::string, tk::real > >& p );

    //! \brief Forward a vector of performance statistics from (LinSysMerger)
    //!   chare group branches to the main proxy
    void grpPerfstat(
      const std::vector< std::pair< std::string, tk::real > >& p );

  private:
    using PerformerProxy = CProxy_Performer;
    using LinSysMergerProxy =
      tk::CProxy_LinSysMerger< CProxy_Conductor, CProxy_Performer >;

    InciterPrint m_print;               //!< Pretty printer
    std::vector< tk::Timer > m_timer;   //!< Timers
    int m_nchare;                       //!< Number of performer chares
    int m_arrTimestampCnt;              //!< Time stamp chare array counter
    int m_grpTimestampCnt;              //!< Time stamp chare group counter
    int m_arrPerfstatCnt;               //!< Perfstat chare array counter
    int m_grpPerfstatCnt;               //!< Perfstat chare group counter
    PerformerProxy m_perfproxy;         //!< Performer chare array
    LinSysMergerProxy m_lsmproxy;       //!< Linear system merger chare group
    //! Time stamps merged from chare array elements
    std::map< std::string, std::vector< tk::real > > m_arrTimestamp;
    //! Time stamps merged from chare group elements
    std::map< std::string, std::vector< tk::real > > m_grpTimestamp;
    //! Performance statistics merged from chare array elements
    std::map< std::string, std::vector< tk::real > > m_arrPerfstat;
    //! Performance statistics merged from chare group elements
    std::map< std::string, std::vector< tk::real > > m_grpPerfstat;

    //! Collect and compute averages of time stamps contributed by chares
    void timestamp(
      const std::vector< std::pair< std::string, tk::real > >& stamp,
      std::map< std::string, std::vector< tk::real > >& map,
      int& counter,
      int max,
      const std::string& of );

    //! Collect and compute performance statistics contributed by chares
    void perfstat( const std::vector< std::pair< std::string, tk::real > >& p,
                   std::map< std::string, std::vector< tk::real > >& map,
                   int& counter,
                   int max,
                   const std::string& of );
};

} // inciter::

#endif // Conductor_h
