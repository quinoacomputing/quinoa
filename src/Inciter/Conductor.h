// *****************************************************************************
/*!
  \file      src/Inciter/Conductor.h
  \author    J. Bakosi
  \date      Fri 29 Jul 2016 02:35:43 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Conductor drives the time integration of systems of systems of PDEs
  \details   Conductor drives the time integration of systems of systems of
    PDEs.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/conductor.ci.

    #### Call graph ####
    The following is a directed acyclic graph (DAG) that outlines the
    asynchronous algorithm implemented in this class The detailed discussion of
    the algorithm is given in the Charm++ interface file conductor.ci, which
    also repeats the graph below using ASCII graphics. On the DAG orange
    fills denote global synchronization points, orange frames with white fill
    are partial synchronization points that overlap with other tasks, and dashed
    lines are potential shortcuts that allow jumping over some of the task-graph
    under some circumstances or optional code paths (taken, e.g., only in DEBUG
    mode). See the detailed discussion in conductor.ci.
    \dot
    digraph "Conductor SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      Diag [ label="Diag"
              tooltip="chares contribute diagnostics"
              URL="\ref inciter::Performer::diagnostics"];
      Eval [ label="Eval"
              tooltip="evaluate time at the end of the time step"
              URL="\ref inciter::Conductor::evaluateTime"];
      Rep [ label="Rep"
              tooltip="output one-liner report"
              URL="\ref inciter::Conductor::report"];
      Diag -> Rep [ style="dashed" ];
      Eval -> Rep [ style="solid" ];
    }
    \enddot
    \include Inciter/conductor.ci
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
#include "NoWarning/particlewriter.decl.h"

namespace inciter {

//! Conductor drives the time integration of a PDE
class Conductor : public CBase_Conductor {

  public:
    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wunused-parameter"
      #pragma clang diagnostic ignored "-Wdeprecated-declarations"
    #elif defined(__GNUC__)
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-parameter"
      #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    #elif defined(__INTEL_COMPILER)
      #pragma warning( push )
      #pragma warning( disable: 1478 )
    #endif
    // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
    // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
    Conductor_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(__GNUC__)
    #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

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

    //! \brief Reduction target indicating that all workers have registered
    //!    with their particle writer branches
    void regcomplete();

    //! Reduction target initiating verification of the boundary conditions set
    void verifybc( CkReductionMsg* msg );

    //! Reduction target as a 2nd (final) of the verification of BCs
    void doverifybc( CkReductionMsg* msg );

    //! \brief Reduction target indicating that all Performer chares have
    //!   finished their initialization step
    void initcomplete();

    //! \brief Reduction target optionally collecting diagnostics, e.g.,
    //!   residuals, from all Performer chares
    void diagnostics( tk::real* d, std::size_t n );

    //! \brief Reduction target indicating that Performer chares contribute no
    //!    diagnostics and we ready to output the one-liner report
    void diagcomplete() { trigger_diag_complete(); }

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
    using TrackerProxy = CProxy_Tracker< PerformerProxy >;
    using ParticleWriterProxy = tk::CProxy_ParticleWriter;
    using PartitionerProxy = CProxy_Partitioner< CProxy_Conductor,
                                                 CProxy_Performer,
                                                 LinSysMergerProxy,
                                                 TrackerProxy,
                                                 ParticleWriterProxy >;

    InciterPrint m_print;                //!< Pretty printer
    int m_nchare;                        //!< Number of performer chares
    uint64_t m_it;                       //!< Iteration count
    tk::real m_t;                        //!< Physical time
    tk::real m_dt;                       //!< Physical time step size
    uint8_t m_stage;                     //!< Stage in multi-stage time stepping
    LinSysMergerProxy m_linsysmerger;    //!< Linear system merger group proxy
    PerformerProxy m_performer;          //!< Performer chare array proxy
    TrackerProxy m_tracker;              //!< Tracker chare array proxy
    ParticleWriterProxy m_particlewriter;//!< Particle writer group proxy
    PartitionerProxy m_partitioner;      //!< Partitioner group proxy
    //! Average communication cost of merging the linear system
    tk::real m_avcost;
    //! Total number of mesh nodes
    std::size_t m_npoin;
    //! Timer tags
    enum class TimerTag { TIMESTEP };
    //! Timers
    std::map< TimerTag, tk::Timer > m_timer;
    //! \brief Aggregate 'old' (as in file) node ID list at which LinSysMerger
    //!   sets boundary conditions, see also Partitioner.h
    std::vector< std::size_t > m_linsysbc;
    //! Diagnostics
    std::vector< tk::real > m_diag;

    //! Compute size of next time step
    tk::real computedt();

    //! Print out time integration header
    void header();

    //! Print out one-liner report on time step
    void report();
};

} // inciter::

#endif // Conductor_h
