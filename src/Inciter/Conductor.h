//******************************************************************************
/*!
  \file      src/Inciter/Conductor.h
  \author    J. Bakosi
  \date      Wed 20 Jan 2016 05:55:18 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Conductor drives the time integration of a PDE
  \details   Conductor drives the time integration of a PDE
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
#include <unordered_map>
#include <unordered_set>

#include "Timer.h"
#include "Types.h"
#include "InciterPrint.h"
#include "Partitioner.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "linsysmerger.decl.h"
#include "conductor.decl.h"
#include "performer.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace inciter {

//! Conductor drives the time integration of a PDE
class Conductor : public CBase_Conductor {

  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  Conductor_SDAG_CODE

  public:
    //! Constructor
    explicit Conductor();

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished reading their part of the computational mesh graph and
    //!   we are ready to compute the computational load
    void load( uint64_t nelem );

    //! Reduction target collecting global mesh node IDs from PEs
    void nodes( CkReductionMsg* msg );

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished setting up the necessary data structures for
    //!   partitioning the computational mesh and we are ready for partitioning
    void partition();

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished distributing the mesh node IDs after partitioning and
    //!   we are ready to start reordering mesh node IDs
    void flatten();

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished preparing their chunk of the mesh connectivity and ready
    //!   for a new order
    void flattened() { trigger_flatten_complete(); }

    //! \brief Reduction target estimating the average communication cost of
    //!   merging the linear system
    void aveCost( tk::real c );

    //! \brief Reduction target estimating the standard deviation of the
    //!   communication cost of merging the linear system
    void stdCost( tk::real c );

    //! Reduction target indicating that all chare groups are ready for workers
    void setup() { m_performer.setup(); }

    //! \brief Reduction target indicating that all linear system merger
    //!   branches have done their part of storing and exporting global row ids
    void rowcomplete();

    //! \brief Reduction target indicating that all Performer chares have
    //!   finished their initialization step
    void initcomplete();

    //! \brief Reduction target indicating that all Performer chares have
    //!   finished a time step and it is time to decide whether to continue
    void evaluateTime();

    //! \brief Reduction target indicating that all ...
    void advance();

    //! Normal finish of time stepping
    void finish();

  private:
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Conductor,
                                                       CProxy_Performer >;
    using PerformerProxy = CProxy_Performer;
    using PartitionerProxy = CProxy_Partitioner< CProxy_Conductor,
                                                 CProxy_Performer,
                                                 LinSysMergerProxy >;

    InciterPrint m_print;               //!< Pretty printer
    int m_nchare;                       //!< Number of performer chares
    uint64_t m_it;                      //!< Iteration count
    tk::real m_t;                       //!< Physical time
    tk::real m_dt;                      //!< Physical time step size
    uint8_t m_stage;                    //!< Stage in multi-stage time stepping
    LinSysMergerProxy m_linsysmerger;   //!< Linear system merger group proxy
    PerformerProxy m_performer;         //!< Performer chare array proxy
    PartitionerProxy m_partitioner;     //!< Partitioner group proxy
    //! Average communication cost of merging the linear system
    tk::real m_avcost;
    //! \brief Communication map for all PEs used for node reordering
    //! \details This map, for all PEs, associates the list of global mesh point
    //!   indices to fellow PE IDs which a PE will receive new node ID numbers
    //!   during reordering. Only data that will be received from PEs with a
    //!   lower index are stored.
    std::vector<
      std::unordered_map< int, std::set< std::size_t > > > m_communication;
    //! Start IDs for each PE for reordering nodes
    std::vector< std::size_t > m_start;
    //! Timer tags
    enum class TimerTag { TIMESTEP };
    //! Timers
    std::map< TimerTag, tk::Timer > m_timer;

    //! Reorder mesh node IDs owned on each PE
    void reorder();

    //! Compute size of next time step
    tk::real computedt();

    //! Print out time integration header
    void header();

    //! Print out one-liner report on time step
    void report();
};

} // inciter::

#endif // Conductor_h
