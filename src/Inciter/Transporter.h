// *****************************************************************************
/*!
  \file      src/Inciter/Transporter.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Transporter drives the time integration of transport equations
  \details   Transporter drives the time integration of transport equations.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/transporter.ci.

    #### Call graph ####
    The following is a directed acyclic graph (DAG) that outlines the
    asynchronous algorithm implemented in this class The detailed discussion of
    the algorithm is given in the Charm++ interface file transporter.ci, which
    also repeats the graph below using ASCII graphics. On the DAG orange
    fills denote global synchronization points that contain or eventually lead
    to global reductions. Dashed lines are potential shortcuts that allow
    jumping over some of the task-graph under some circumstances or optional
    code paths (taken, e.g., only in DEBUG mode). See the detailed discussion in
    transporter.ci.
    \dot
    digraph "Transporter SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      Diag [ label="Diag"
              tooltip="chares contribute diagnostics"
              URL="\ref inciter::Carrier::diagnostics"];
      Out [ label="Out"
              tooltip="particles output to file"
              URL="\ref inciter::Carrier::doWriteParticles"];
      Eval [ label="Eval"
              tooltip="evaluate time at the end of the time step"
              URL="\ref inciter::Transporter::evaluateTime"];
      Diag -> Eval [ style="solid" ];
      Out -> Eval [ style="solid" ];
    }
    \enddot
    \include Inciter/transporter.ci
*/
// *****************************************************************************
#ifndef Transporter_h
#define Transporter_h

#include <map>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "Timer.h"
#include "Types.h"
#include "InciterPrint.h"
#include "Partitioner.h"
#include "VectorReducer.h"
#include "ParticleWriter.h"
#include "Progress.h"

#include "NoWarning/carrier.decl.h"

namespace inciter {

//! Transporter drives the time integration of transport equations
class Transporter : public CBase_Transporter {

  public:
    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wunused-parameter"
      #pragma clang diagnostic ignored "-Wdeprecated-declarations"
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-parameter"
      #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    #elif defined(__INTEL_COMPILER)
      #pragma warning( push )
      #pragma warning( disable: 1478 )
    #endif
    // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
    // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
    Transporter_SDAG_CODE
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #elif defined(STRICT_GNUC)
      #pragma GCC diagnostic pop
    #elif defined(__INTEL_COMPILER)
      #pragma warning( pop )
    #endif

    //! Constructor
    explicit Transporter();

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished reading their part of the computational mesh graph and
    //!   we are ready to compute the computational load
    void load( uint64_t nelem );

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished setting up the necessary data structures for
    //!   partitioning the computational mesh and we are ready for partitioning
    void partition();

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished distributing its global mesh node IDs and they are ready
    //!   for preparing (flattening) their owned mesh node IDs for reordering
    void distributed();

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
    void setup();

    //! Non-reduction target for receiving progress report on reading mesh graph
    void pegraph() { m_progGraph.inc<0>(); }

    //! Non-reduction target for receiving progress report on partitioning mesh
    void pepartitioned() { m_progPart.inc<0>(); }
    //! Non-reduction target for receiving progress report on distributing mesh
    void pedistributed() { m_progPart.inc<1>(); }

    //! Non-reduction target for receiving progress report on flattening mesh
    void peflattened() { m_progReorder.inc<0>(); }
    //! Non-reduction target for receiving progress report on node ID mask
    void pemask() { m_progReorder.inc<1>(); }
    //! Non-reduction target for receiving progress report on reordering mesh
    void pereordered() { m_progReorder.inc<2>(); }
    //! Non-reduction target for receiving progress report on computing bounds
    void pebounds() { m_progReorder.inc<3>(); }

    //! Non-reduction target for receiving progress report on computing row IDs
    void perowcomplete() { m_progSetup.inc<0>(); }
    //! Non-reduction target for receiving progress report on matching BCs
    void chbcmatched() { m_progSetup.inc<1>(); }
    //! Non-reduction target for receiving progress report on computing BCs
    void pebccomplete() { m_progSetup.inc<2>(); }

    //! Non-reduction target for receiving progress report on the initial guess
    void pesolcomplete() { m_progInit.inc<0>(); }
    //! Non-reduction target for receiving progress report on outputing ICs
    void chic() { m_progInit.inc<1>(); }
    //! Non-reduction target for receiving progress report on computing the LHS
    void chlhs() { m_progInit.inc<2>(); }

    //! Non-reduction target for receiving progress report on computing the RHS
    void chrhs() { m_progStep.inc<0>(); }
    //! Non-reduction target for receiving progress report on solving the system
    void pesolve() { m_progStep.inc<1>(); }
    //! Non-reduction target for receiving progress report on limiting
    void chlim() { m_progStep.inc<2>(); }
    //! Non-reduction target for receiving progress report on tracking particles
    void chtrack() { m_progStep.inc<3>(); }

    //! \brief Reduction target indicating that all linear system merger
    //!   branches have done their part of storing and exporting global row ids
    void rowcomplete();

    //! \brief Reduction target indicating that all Carriers have finished
    //!   computing/receiving their part of the nodal volumes
    void volcomplete();

    //! \brief Reduction target yielding a single minimum time step size across
    //!   all workers
    void dt( tk::real* d, std::size_t n );

    //! Reduction target initiating verification of the boundary conditions set
    void verifybc( CkReductionMsg* msg );

    //! Reduction target as a 2nd (final) of the verification of BCs
    void doverifybc( CkReductionMsg* msg );

    //! \brief Reduction target indicating that all Carrier chares have
    //!   finished their initialization step
    void initcomplete();

    //! Reduction target indicating the particle communication is complete
    void parcomcomplete() { m_carrier.out(); }

    //! \brief Reduction target indicating that all workers have sent their
    //!   number of particles to be output
    void nparcomplete() { m_carrier.doWriteParticles(); }

    //! \brief Reduction target optionally collecting diagnostics, e.g.,
    //!   residuals, from all Carrier chares
    void diagnostics( tk::real* d, std::size_t n );

    //! \brief Reduction target indicating that Carrier chares contribute no
    //!    diagnostics and we ready to output the one-liner report
    void diagcomplete() { diag_complete(); }

    //! \brief Reduction target indicating that all particles writers have
    //!   finished outputing particles to file
    //! \details This function is a Charm++ reduction target that is called when
    //!   all carrier chares have finished communicating particles
    void outcomplete() { out_complete(); }

    //! \brief Reduction target indicating that the linear system mergers are
    //!   ready for the next time step
    void computedt() { m_carrier.dt(); }

    //! Normal finish of time stepping
    void finish();

    //! \brief Reduction target outputing diagnostics
    void verified() { m_print.diag( "AEC verified" ); }

  private:
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Transporter,
                                                       CProxy_Carrier,
                                                       AuxSolverLumpMassDiff >;
    using CarrierProxy = CProxy_Carrier;
    using ParticleWriterProxy = tk::CProxy_ParticleWriter< CProxy_Transporter >;
    using PartitionerProxy = CProxy_Partitioner< CProxy_Transporter,
                                                 CarrierProxy,
                                                 LinSysMergerProxy,
                                                 ParticleWriterProxy >;

    InciterPrint m_print;                //!< Pretty printer
    int m_nchare;                        //!< Number of carrier chares
    uint64_t m_it;                       //!< Iteration count
    tk::real m_t;                        //!< Physical time
    tk::real m_dt;                       //!< Physical time step size
    uint8_t m_stage;                     //!< Stage in multi-stage time stepping
    LinSysMergerProxy m_linsysmerger;    //!< Linear system merger group proxy
    CarrierProxy m_carrier;              //!< Carrier chare array proxy
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
    // Progress object for task "Partitioning and distributing mesh"
    tk::Progress< 2 > m_progPart;
    // Progress object for task "Creating partitioners and reading mesh graph"
    tk::Progress< 1 > m_progGraph;
    // Progress object for task "Reordering mesh"
    tk::Progress< 4 > m_progReorder;
    // Progress object for task "Computing row IDs, querying BCs, ..."
    tk::Progress< 3 > m_progSetup;
    // Progress object for task "Setting and output ICs, ..."
    tk::Progress< 3 > m_progInit;
    // Progress object for sub-tasks of a time step
    tk::Progress< 4 > m_progStep;

    //! Print out time integration header
    void header();

    //! Evaluate time step and output one-liner report
    void evaluateTime();
};

} // inciter::

#endif // Transporter_h
