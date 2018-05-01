// *****************************************************************************
/*!
  \file      src/Inciter/Transporter.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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
      Load [ label="Load"
              tooltip="Load is computed"
              URL="\ref inciter::Transporter::load"];
      Centroid [ label="Centroid"
              tooltip="Cell centroids have been computed"
              URL="\ref inciter::Transporter::centroid"];
      Refined [ label="Refined"
              tooltip="Mesh has been optionally initially refined"
              URL="\ref inciter::Transporter::refined"];
      Part [ label="Part"
              tooltip="Partition mesh"
              URL="\ref inciter::Partitioner::partition"];
      Load -> Part [ style="solid" ];
      Centroid -> Part [ style="solid" ];
      Refined -> Part [ style="solid" ];
      MinStat [ label="MinStat"
              tooltip="chares contribute to minimum mesh cell statistics"
              URL="\ref inciter::Discretization::stat"];
      MaxStat [ label="MaxStat"
              tooltip="chares contribute to maximum mesh cell statistics"
              URL="\ref inciter::Discretization::stat"];
      SumStat [ label="SumStat"
              tooltip="chares contribute to sum mesh cell statistics"
              URL="\ref inciter::Discretization::stat"];
      PDFStat [ label="PDFStat"
              tooltip="chares contribute to PDF mesh cell statistics"
              URL="\ref inciter::Discretization::stat"];
      Com [ label="Com"
              tooltip="chares contributed to mesh cell statistics"
              URL="\ref tk::Solver::charecom::"];
      Stat [ label="Stat"
              tooltip="all chares have contributed to mesh cell statistics"
              URL="\ref inciter::Discretization::stat"];
      Setup [ label="Setup"
              tooltip="start computing row IDs, querying BCs, outputing mesh"
              URL="\ref inciter::Transporter::setup"];
      MinStat -> Stat [ style="solid" ];
      MaxStat -> Stat [ style="solid" ];
      SumStat -> Stat [ style="solid" ];
      PDFStat -> Stat [ style="solid" ];
      Com -> Stat [ style="solid" ];
      Stat -> Setup [ style="solid" ];
    }
    \enddot
    \include Inciter/transporter.ci

    #### Call graph with the interactions of Transporter, Solver, and MatCG ####
    The following a DAG documentaion of the itneractions _among_ the classes
    Transporter (single chare), Solver (chare group), and MatCG (chare array).
    The prefixes p_, s_, t_ c_, respectively, denote Partitioner, Solver,
    Transporter, and MatCG, which help identifying which class' member function
    the label is associated with. (These prefixes only show up in the source of
    "dot", used to generate the visual graph. Note that MatMatCG can be thought
    of as a child that inherits from class Discretization. Both MatCG and
    Discretization are Charm++ chare arrays whose elements are abound together
    when migrated and via class Scheme they are used in a runtime polymorphic
    fashion. This means when the prefix c is used in the DAG below, the member
    function might be in the (base) Discretization class instead of in MatCG.
    Note that the graph below is partial as it only partially contains how this
    hooks into Partitioner.
    \dot
    digraph "Transporter-Solver-MatCG SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      p_dcreate [ label="Partitioner::create"
              tooltip="Partitioners create Discretization (base) workers"
              URL="\ref inciter::Partitioner::createWorkers"];
      p_ccreate [ label="Partitioner::createWorkers"
              tooltip="Partitioners create MatCG (child) workers"
              URL="\ref inciter::Partitioner::createWorkers"];
      s_nchare [ label="Solver::nchare"
              tooltip="set number of worker chares expected to contribute"
              URL="\ref tk::Solver::nchare"];
      s_bounds [ label="Solver::bounds"
              tooltip="receive lower and upper global node ids"
              URL="\ref tk::Solver::bounds"];
      t_coord [ label="Transporter::coord"
              tooltip="read mesh coordinates"
              URL="\ref inciter::Transporter::coord"];
      t_stat [ label="Transporter::stat"
              tooltip="output mesh statistics"
              URL="\ref inciter::Transporter::stat"];
      c_setup [ label="MatCG::setup"
              tooltip="workers start setting up PE communications, output mesh"
              URL="\ref inciter::MatCG::setup"];
      s_sol [ label="Solver::sol"
              tooltip="assemble unknown/solution vector"
              URL="\ref tk::Solver:charesol"];
      s_lhs [ label="Solver::lhs"
              tooltip="assemble LHS"
              URL="\ref tk::Solver:charelhs"];
      c_dt [ label="MatCG::dt"
              tooltip="workers compute size of next time step"
              URL="\ref inciter::MatCG::dt"];
      c_advance [ label="MatCG::advance"
              tooltip="advance to next time step"
              URL="\ref inciter::MatCG::advance"];
      c_rhs [ label="MatCG::rhs"
              tooltip="workers compute new RHS"
              URL="\ref inciter::MatCG::rhs"];
      s_rhs [ label="Solver::rhs"
              tooltip="solver assemble new RHS"
              URL="\ref tk::Solver::rhs"];
      s_solve [ label="Solver::solve"
              tooltip="solve linear system"
              URL="\ref tk::Solver::solve"];
      s_com [ label="Solver::charecom"
              tooltip="setup PE communication maps"
              URL="\ref tk::Solver::charecom"];
      c_update [ label="MatCG::update"
              tooltip="workers update their field data to new solution"
              URL="\ref inciter::MatCG::updateLow"];
      t_diag [ label="Transporter::diag"
              tooltip="diagnostics have been output"
              URL="\ref inciter::Transporter::diagnostics"];
      s_next [ label="Solver::next"
              tooltip="prepare for next time step"
              URL="\ref tk::Solver::next"];
      t_comfinal [ label="Transporter::comfinal"
              tooltip="communication maps are complete on Trarnsporter"
              URL="\ref tk::Solver::comfinal"];
      c_vol [ label="MatCG::vol"
              tooltip="compute nodal mesh volumes"
              URL="\ref tk::MatCG::vol"];
      t_vol [ label="Transporter::vol"
              tooltip="nodal mesh volumes complete, start computing total volume"
              URL="\ref inicter::Transporter::vol"];
      t_start [ label="Transporter::start"
              tooltip="start time stepping"
              URL="\ref inicter::Transporter::start"];
      c_totalvol [ label="MatCG::totalvol"
              tooltip="compute total mesh volume"
              URL="\ref tk::MatCG::totalvol"];
      t_totalvol [ label="Transporter::totalvol"
              tooltip="total mesh volume computed, start with mesh stats"
              URL="\ref inicter::Transporter::totalvol"];
      c_minstat [ label="MatCG::stat(min)"
              tooltip="compute mesh statistics finding minima"
              URL="\ref inciter::MatCG::stat"];
      c_maxstat [ label="MatCG::stat(max)"
              tooltip="compute mesh statistics finding maxima"
              URL="\ref inciter::MatCG::stat"];
      c_sumstat [ label="MatCG::stat(sum)"
              tooltip="compute mesh statistics finding sums"
              URL="\ref inciter::MatCG::stat"];
      c_pdfstat [ label="MatCG::stat(pdf)"
              tooltip="compute mesh statistics computing PDFs"
              URL="\ref inciter::MatCG::stat"];
      t_minstat [ label="Transporter::minstat"
              tooltip="compute mesh statistics finding global minima"
              URL="\ref inciter::Transporter::minstat"];
      t_maxstat [ label="Transporter::maxstat"
              tooltip="compute mesh statistics finding global maxima"
              URL="\ref inciter::Transporter::maxstat"];
      t_sumstat [ label="Transporter::sumstat"
              tooltip="compute mesh statistics finding global sums"
              URL="\ref inciter::Transporter::sumstat"];
      t_pdfstat [ label="Transporter::pdfstat"
              tooltip="compute mesh statistics computing global PDFs"
              URL="\ref inciter::Transporter::pdfstat"];
      p_dcreate -> s_bounds [ style="solid" ];
      s_nchare -> t_coord [ style="solid" ];
      s_bounds -> t_coord [ style="solid" ];
      t_coord -> c_vol [ style="solid" ];
      c_vol -> t_vol [ style="solid" ];
      t_vol -> c_totalvol [ style="solid" ];
      p_ccreate -> s_com [ style="solid" ];
      s_com -> t_comfinal [ style="solid" ];
      c_totalvol -> t_totalvol [ style="solid" ];
      t_totalvol -> p_ccreate [ style="solid" ];
      t_totalvol -> c_minstat [ style="solid" ];
      t_totalvol -> c_maxstat [ style="solid" ];
      t_totalvol -> c_sumstat [ style="solid" ];
      t_totalvol -> c_pdfstat [ style="solid" ];
      c_minstat -> t_minstat [ style="solid" ];
      c_maxstat -> t_maxstat [ style="solid" ];
      c_sumstat -> t_sumstat [ style="solid" ];
      c_pdfstat -> t_pdfstat [ style="solid" ];
      t_minstat -> t_stat [ style="solid" ];
      t_maxstat -> t_stat [ style="solid" ];
      t_sumstat -> t_stat [ style="solid" ];
      t_pdfstat -> t_stat [ style="solid" ];
      t_comfinal-> t_stat [ style="solid" ];
      t_stat -> c_setup [ style="solid" ];
      c_setup -> s_sol [ style="solid" ];
      c_setup -> s_lhs [ style="solid" ];
      c_setup -> t_start [ style="solid" ];
      t_start -> c_dt [ style="solid" ];
      c_dt -> c_advance [ style="solid" ];
      c_advance -> c_rhs [ style="solid" ];
      c_rhs -> s_rhs [ style="solid" ];
      s_sol -> s_solve [ style="solid" ];
      s_lhs -> s_solve [ style="solid" ];
      s_rhs -> s_solve [ style="solid" ];
      c_update -> t_diag [ style="solid" ];
      s_solve -> c_update [ style="solid" ];
      t_diag -> s_next [ style="solid" ];
      s_next -> c_dt [ style="solid" ];
    }
    \enddot

    #### Call graph documenting the interactions of Transporter and DG ####
    The following a DAG documentaion of the itneractions _among_ the classes
    Transporter (single chare) and DG (chare array). The prefixes p_, t_ d_,
    respectively, denote Partitioner, Transporter, and DG, which help
    identifying which class' member function the label is associated with.
    (These prefixes only show up in the source of "dot", used to generate the
    visual graph. Note that DG can be thought of as a child that inherits from
    class Discretization. Both DG and Discretization are Charm++ chare arrays
    whose elements are abound together when migrated and via class Scheme they
    are used in a runtime polymorphic fashion. This means when the prefix c is
    used in the DAG below, the member function might be in the (base)
    Discretization class instead of in DG. Note that the graph below is partial
    as it only partially contains how this hooks into Partitioner.
    \dot
    digraph "Transporter-DG SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      p_dcreate [ label="Partitioner::create"
              tooltip="Partitioners create Discretization (base) workers"
              URL="\ref inciter::Partitioner::createWorkers"];
      p_ccreate [ label="Partitioner::createWorkers"
              tooltip="Partitioners create DG (child) workers"
              URL="\ref inciter::Partitioner::createWorkers"];
      t_coord [ label="Transporter::coord"
              tooltip="read mesh coordinates"
              URL="\ref inciter::Transporter::coord"];
      t_stat [ label="Transporter::stat"
              tooltip="output mesh statistics"
              URL="\ref inciter::Transporter::stat"];
      d_setup [ label="DG::setup"
              tooltip="workers start setting up PE communications, output mesh"
              URL="\ref inciter::DG::setup"];
      t_comfinal [ label="Transporter::comfinal"
              tooltip="communications maps final"
              URL="\ref tk::Transporter:comfinal"];
      d_dt [ label="DG::dt"
              tooltip="workers compute size of next time step"
              URL="\ref inciter::DG::dt"];
      d_advance [ label="DG::advance"
              tooltip="advance to next time step"
              URL="\ref inciter::DG::advance"];
      t_diag [ label="Transporter::diag"
              tooltip="diagnostics have been output"
              URL="\ref inciter::Transporter::diagnostics"];
      d_vol [ label="DG::vol"
              tooltip="compute nodal mesh volumes"
              URL="\ref tk::DG::vol"];
      t_vol [ label="Transporter::vol"
              tooltip="nodal mesh volumes complete, start computing total volume"
              URL="\ref inicter::Transporter::vol"];
      d_totalvol [ label="DG::totalvol"
              tooltip="compute total mesh volume"
              URL="\ref tk::DG::totalvol"];
      t_totalvol [ label="Transporter::totalvol"
              tooltip="total mesh volume computed, start with mesh stats"
              URL="\ref inicter::Transporter::totalvol"];
      d_minstat [ label="DG::stat(min)"
              tooltip="compute mesh statistics finding minima"
              URL="\ref inciter::DG::stat"];
      d_maxstat [ label="DG::stat(max)"
              tooltip="compute mesh statistics finding maxima"
              URL="\ref inciter::DG::stat"];
      d_sumstat [ label="DG::stat(sum)"
              tooltip="compute mesh statistics finding sums"
              URL="\ref inciter::DG::stat"];
      d_pdfstat [ label="DG::stat(pdf)"
              tooltip="compute mesh statistics computing PDFs"
              URL="\ref inciter::DG::stat"];
      t_minstat [ label="Transporter::minstat"
              tooltip="compute mesh statistics finding global minima"
              URL="\ref inciter::Transporter::minstat"];
      t_maxstat [ label="Transporter::maxstat"
              tooltip="compute mesh statistics finding global maxima"
              URL="\ref inciter::Transporter::maxstat"];
      t_sumstat [ label="Transporter::sumstat"
              tooltip="compute mesh statistics finding global sums"
              URL="\ref inciter::Transporter::sumstat"];
      t_pdfstat [ label="Transporter::pdfstat"
              tooltip="compute mesh statistics computing global PDFs"
              URL="\ref inciter::Transporter::pdfstat"];
      t_start [ label="Transporter::start"
              tooltip="start time stepping"
              URL="\ref inicter::Transporter::start"];
      p_dcreate -> t_coord [ style="solid" ];
      t_coord -> d_vol [ style="solid" ];
      d_vol -> t_vol [ style="solid" ];
      t_vol -> d_totalvol [ style="solid" ];
      d_totalvol -> t_totalvol [ style="solid" ];
      t_totalvol -> p_ccreate [ style="solid" ];
      t_totalvol -> d_minstat [ style="solid" ];
      t_totalvol -> d_maxstat [ style="solid" ];
      t_totalvol -> d_sumstat [ style="solid" ];
      t_totalvol -> d_pdfstat [ style="solid" ];
      d_minstat -> t_minstat [ style="solid" ];
      d_maxstat -> t_maxstat [ style="solid" ];
      d_sumstat -> t_sumstat [ style="solid" ];
      d_pdfstat -> t_pdfstat [ style="solid" ];
      t_minstat -> t_stat [ style="solid" ];
      t_maxstat -> t_stat [ style="solid" ];
      t_sumstat -> t_stat [ style="solid" ];
      t_pdfstat -> t_stat [ style="solid" ];
      t_stat -> d_setup [ style="solid" ];
      p_ccreate -> t_comfinal [ style="solid" ];
      t_comfinal -> t_stat [ style="solid" ];
      d_setup -> t_start [ style="solid" ];
      t_start -> d_dt [ style="solid" ];
      d_dt -> d_advance [ style="solid" ];
      d_advance -> t_diag [ style="solid" ];
      t_diag -> d_dt [ style="solid" ];
    }
    \enddot

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
#include "Progress.h"
#include "Scheme.h"
#include "BoundaryConditions.h"

namespace inciter {

//! Transporter drives the time integration of transport equations
class Transporter : public CBase_Transporter {

  private:
    //! Indices for progress report on mesh read (and prep for partitioning)
    enum ProgMesh{ READ=0, REFINE, CENTROID };
    //! Indices for progress report on mesh partitioning
    enum ProgPart{ PART=0, DIST };
    //! Indices for progress report on mesh reordering
    enum ProgReord{ FLAT=0, GATHER, QUERY, MASK, REORD, BOUND };

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

    //! Reduction target indicating that the mesh has been read from file
    void load( uint64_t nelem );

    //! \brief Reduction target indicating that optional initial mesh refinement
    //!   has been completed on all PEs
    void refined();

    //! Reduction target indicating that centroids have been computed all PEs
    void centroid();

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished distributing its global mesh node IDs and they are ready
    //!   for preparing (flattening) their owned mesh node IDs for reordering
    void distributed();

    //! \brief Reduction target indicating that all Partitioner chare groups
    //!   have finished flattening its global mesh node IDs and they are ready
    //!   for computing the communication maps required for node ID reordering
    void flattened();

    //! Reduction target estimating the average communication cost among all PEs
    void aveCost( tk::real c );

    //! \brief Reduction target estimating the standard deviation of the
    //!   communication cost among all PEs
    void stdCost( tk::real c );

    //! \brief Reduction target indicating that all chare groups are ready for
    //!   workers to read their mesh coordinates
    void coord();

    //! Non-reduction target for receiving progress report on reading mesh
    void peread() { m_progMesh.inc< READ >(); }
    //! Non-reduction target for receiving progress report on mesh refinement
    void perefined() { m_progMesh.inc< REFINE >(); }
    //! Non-reduction target for receiving progress report on mesh centroids
    void pecentroid() { m_progMesh.inc< CENTROID >(); }

    //! Non-reduction target for receiving progress report on partitioning mesh
    void pepartitioned() { m_progPart.inc< PART >(); }
    //! Non-reduction target for receiving progress report on distributing mesh
    void pedistributed() { m_progPart.inc< DIST >(); }

    //! Non-reduction target for receiving progress report on flattening mesh
    void peflattened() { m_progReorder.inc< FLAT >(); }
    //! Non-reduction target for receiving progress report on node ID gather
    void pegather() { m_progReorder.inc< GATHER >(); }
    //! Non-reduction target for receiving progress report on node ID query
    void pequery() { m_progReorder.inc< QUERY >(); }
    //! Non-reduction target for receiving progress report on node ID mask
    void pemask() { m_progReorder.inc< MASK >(); }
    //! Non-reduction target for receiving progress report on reordering mesh
    void pereordered() { m_progReorder.inc< REORD >(); }
    //! Non-reduction target for receiving progress report on computing bounds
    void pebounds() { m_progReorder.inc< BOUND >(); }

    //! \brief Reduction target indicating that the communication has been
    //!    established among PEs
    void comfinal();

    //! Reduction target summing total mesh volume
    void totalvol( tk::real v );

    //! \brief Reduction target indicating that all workers have finished
    //!   computing/receiving their part of the nodal volumes
    void vol();

    //! \brief Reduction target yielding the minimum mesh statistics across
    //!   all workers
    void minstat( tk::real* d, std::size_t n );

    //! \brief Reduction target yielding the maximum mesh statistics across
    //!   all workers
    void maxstat( tk::real* d, std::size_t n );

    //! \brief Reduction target yielding the sum of mesh statistics across
    //!   all workers
    void sumstat( tk::real* d, std::size_t n );

    //! \brief Reduction target yielding PDF of mesh statistics across all
    //!    workers
    void pdfstat( CkReductionMsg* msg );

    //! \brief Reduction target optionally collecting diagnostics, e.g.,
    //!   residuals, from all  worker chares
    void diagnostics( CkReductionMsg* msg );

    //! Start time stepping
    void start();

    //! Reset linear solver for next time step
    void next() { m_solver.next(); }

    //! Normal finish of time stepping
    void finish();

  private:
    InciterPrint m_print;                //!< Pretty printer
    int m_nchare;                        //!< Number of worker chares
    uint64_t m_nelem;                    //!< Total number of mesh elements
    uint64_t m_chunksize;                //!< Number of elements per PE
    uint64_t m_remainder;                //!< Number elements added to last PE
    tk::CProxy_Solver m_solver;          //!< Linear system solver group proxy
    CProxy_BoundaryConditions m_bc;      //!< Boundary conditions group proxy
    Scheme m_scheme;                     //!< Discretization scheme
    CProxy_Partitioner m_partitioner;    //!< Partitioner group proxy
    //! Average communication cost of merging the linear system
    tk::real m_avcost;
     //! Total mesh volume
    tk::real m_V;
    //! Total number of mesh nodes
    std::size_t m_npoin;
    //! Minimum mesh statistics
    std::array< tk::real, 2 > m_minstat;
    //! Maximum mesh statistics
    std::array< tk::real, 2 > m_maxstat;
    //! Average mesh statistics
    std::array< tk::real, 2 > m_avgstat;
    //! Timer tags
    enum class TimerTag { MESH_PREP=0 };
    //! Timers
    std::map< TimerTag, tk::Timer > m_timer;
    //! \brief Aggregate 'old' (as in file) node ID list at which Solver
    //!   sets boundary conditions, see also Partitioner.h
    std::vector< std::size_t > m_linsysbc;
    //! \brief Progress object for task "Creating partitioners, reading, and
    //!    optionally refining mesh"
    tk::Progress< 3 > m_progMesh;
    // Progress object for task "Partitioning and distributing mesh"
    tk::Progress< 2 > m_progPart;
    // Progress object for task "Reordering mesh"
    tk::Progress< 6 > m_progReorder;

    //! Create linear solver group
    void createSolver();

    //! Create mesh partitioner and boundary condition object group
    void createPartitioner();

    //! Start partitioning the mesh
    void partition();

    //! Configure and write diagnostics file header
    void diagHeader();

    //! Echo diagnostics on mesh statistics
    void stat();

    //! Query variable names for all equation systems to be integrated
    //! \param[in] eq Equation system whose variable names to query
    //! \param[in,out] var Vector of strings to which we append the variable
    //!   names for this equation. We append as many strings as many scalar
    //!   variables are in the equation system given by eq.
    template< class Eq >
    void varnames( const Eq& eq, std::vector< std::string >& var ) {
      auto o = eq.names();
      var.insert( end(var), begin(o), end(o) );
    }
};

} // inciter::

#endif // Transporter_h
