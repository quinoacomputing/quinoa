// *****************************************************************************
/*!
  \file      src/Inciter/Conductor.h
  \author    J. Bakosi
  \date      Mon 11 Jul 2016 12:47:46 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Conductor drives the time integration of a PDE
  \details   Conductor drives the time integration of a PDE
    The implementation uses the Charm++ runtime system and is fully asynchronous,
    overlapping computation, communication as well as I/O. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/conductor.ci.
*/
// *****************************************************************************
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
#include "VectorReducer.h"

#include "NoWarning/performer.decl.h"

namespace inciter {

//! Conductor drives the time integration of a PDE
class Conductor : public CBase_Conductor {

  public:
    //! Constructor
    explicit Conductor();

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished reading their part of the computational mesh graph and
    //!   we are ready to compute the computational load
    void load( uint64_t nelem );

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished setting up the necessary data structures for
    //!   partitioning the computational mesh and we are ready for partitioning
    void partition();

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished flattening its global mesh node IDs and they are ready
    //!   for computing the communication maps required for node ID reordering
    void flattened() { m_partitioner.gather(); }

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

    //! Reduction target initiating verification of the boundary conditions set
    void verifybc( CkReductionMsg* msg );

    //! Reduction target as a 2nd (final) of the verification of BCs
    void doverifybc( CkReductionMsg* msg );

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
    //! Timer tags
    enum class TimerTag { TIMESTEP };
    //! Timers
    std::map< TimerTag, tk::Timer > m_timer;
    //! \brief Aggregate 'old' (as in file) node ID list at which LinSysMerger
    //!   sets boundary conditions, see also Partitioner.h
    std::vector< std::size_t > m_linsysbc;

    //! Compute size of next time step
    tk::real computedt();

    //! Print out time integration header
    void header();

    //! Print out one-liner report on time step
    void report();
};

} // inciter::

#endif // Conductor_h
