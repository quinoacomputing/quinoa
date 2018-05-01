// *****************************************************************************
/*!
  \file      src/LinSys/Solver.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ chare linear system merger group to solve a linear system
  \details   Charm++ chare linear system merger group used to collect and
    assemble the left hand side matrix (lhs), the right hand side (rhs) vector,
    and the solution (unknown) vector from individual worker
    chares. Beside collection and assembly, the system is also solved. The
    solution is outsourced to hypre, an MPI-only library. Once the solution is
    available, the individual worker chares are updated with the new solution.

    This class assembles and solves two linear systems, whose rhs vectors may
    change during time stepping. One of the two linear systems is called a
    high-order and the other one is the low-order linear system.

    Characteristics of the low-order linear system: (1) the left hand side
    matrix is diagonal, (2) the right hand side vector is a combination of the
    high-order right hand side vector and another vector, assembled separately
    (and overlapped with that of the high-order system). This dual system
    solution is done here as required by the flux-corrected transport algorithm
    used for transport equations in inciter.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/LinSys/solver.ci.

    #### Call graph ####
    The following is a directed acyclic graph (DAG) that outlines the
    asynchronous algorithm implemented in this class The detailed discussion of
    the algorithm is given in the Charm++ interface file solver.ci. On the
    DAG orange fills denote global synchronization points that contain or
    eventually lead to global reductions. Dashed lines are potential shortcuts
    that allow jumping over some of the task-graph under some circumstances or
    optional code paths (taken, e.g., only in DEBUG mode). See the detailed
    discussion in solver.ci.
    \dot
    digraph "Solver SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      ChCom [ label="ChCom"
              tooltip="chares contribute their global row IDs"
              URL="\ref tk::Solver::charecom"];
      ChBC [ label="ChBC"
             tooltip="chares contribute their global node IDs at which
                      they can set boundary conditions"
             URL="\ref tk::Solver::charebc"];
      ComComplete   [ label="ComComplete"  style="solid"
              tooltip="all linear system solver branches have done heir part of
              storing and exporting global row ids"
              URL="\ref inciter::Transporter::comcomplete"];
      Setup [ label="Setup"
              tooltip="Workers start setting and outputing ICs, computing
                       initial dt, computing LHS"
              URL="\ref inciter::MatCG::setup"];
      dt [ label="dt"
           tooltip="Worker chares compute their minimum time step size"
            style="solid"];
      ComFinal  [ label="ComFinal"
              tooltip="start converting row IDs to hypre format"
              URL="\ref tk::Solver::comfinal"];
      t_start [ label="Transporter::start"
              tooltip="start time stepping"
              URL="\ref inicter::Transporter::start"];
      LhsBC [ label="LhsBC"
              tooltip="set boundary conditions on the left-hand side matrix"
              URL="\ref tk::Solver::lhsbc"];
      RhsBC [ label="RhsBC"
              tooltip="set boundary conditions on the right-hand side vector"
              URL="\ref tk::Solver::rhsbc"];
      ChSol [ label="ChSol"
              tooltip="chares contribute their solution vector nonzeros"
              URL="\ref tk::Solver::charesol"];
      ChLhs [ label="ChLhs"
              tooltip="chares contribute their left hand side matrix nonzeros"
              URL="\ref tk::Solver::charelhs"];
      ChLowLhs [ label="ChLowLhs"
              tooltip="chares contribute their low-order left hand side matrix"
              URL="\ref tk::Solver::charelow"];
      ChRhs [ label="ChRhs"
              tooltip="chares contribute their right hand side vector nonzeros"
              URL="\ref tk::Solver::charesol"];
      ChLowRhs [ label="ChLowRhs"
              tooltip="chares contribute their contribution to the low-order
                       right hand side vector"
              URL="\ref tk::Solver::charelowrhs"];
      HypreRow [ label="HypreRow"
              tooltip="convert global row ID vector to hypre format"
              URL="\ref tk::Solver::hyprerow"];
      HypreSol [ label="HypreSol"
              tooltip="convert solution vector to hypre format"
              URL="\ref tk::Solver::hypresol"];
      HypreLhs [ label="HypreLhs"
              tooltip="convert left hand side matrix to hypre format"
              URL="\ref tk::Solver::hyprelhs"];
      HypreRhs [ label="HypreRhs"
              tooltip="convert right hand side vector to hypre format"
              URL="\ref tk::Solver::hyprerhs"];
      Sol [ label="Sol"
              tooltip="fill/set solution vector"
              URL="\ref tk::Solver::sol"];
      Lhs [ label="Lhs"
              tooltip="fill/set left hand side matrix"
              URL="\ref tk::Solver::lhs"];
      Rhs [ label="Rhs"
              tooltip="fill/set right hand side vector"
              URL="\ref tk::Solver::rhs"];
      Solve [ label="Solve" tooltip="solve linear system"
              URL="\ref tk::Solver::solve"];
      LowSolve [ label="LowSolve" tooltip="solve low-order linear system"
              URL="\ref tk::Solver::lowsolve"];
      Upd [ label="Upd" tooltip="update solution"
                 style="solid"
                URL="\ref tk::Solver::updateSol"];
      LowUpd [ label="LowUpd" tooltip="update low-order solution"
               style="solid"
               URL="\ref tk::Solver::updateLowol"];
      ChCom -> ComComplete [ style="solid" ];
      ComComplete -> Setup [ style="solid" ];
      ComFinal -> HypreRow [ style="solid" ];
      ChLhs -> LhsBC [ style="solid" ];
      ChRhs -> RhsBC [ style="solid" ];
      LhsBC -> HypreLhs [ style="solid" ];
      RhsBC -> HypreRhs [ style="solid" ];
      RhsBC -> LowSolve [ style="solid" ];
      ChLowRhs -> LowSolve [ style="solid" ];
      ChLowLhs -> LowSolve [ style="solid" ];
      Setup -> ComFinal [ style="solid" ];
      Setup -> ChSol [ style="solid" ];
      Setup -> ChLhs [ style="solid" ];
      Setup -> ChLowLhs [ style="solid" ];
      Setup -> t_start [ style="solid" ];
      t_start -> dt [ style="solid" ];
      dt -> ChRhs [ style="solid" ];
      dt -> ChLowRhs [ style="solid" ];
      dt -> ChBC [ style="solid" ];
      ChBC -> LhsBC [ style="solid" ];
      ChBC -> RhsBC [ style="solid" ];
      HypreRow -> Sol [ style="solid" ];
      HypreRow -> Lhs [ style="solid" ];
      HypreRow -> Rhs [ style="solid" ];
      ChSol -> HypreSol -> Sol -> Solve [ style="solid" ];
      HypreLhs -> Lhs -> Solve [ style="solid" ];
      HypreRhs -> Rhs -> Solve [ style="solid" ];
      Solve -> Upd [ style="solid" ];
      LowSolve -> LowUpd [ style="solid" ];
    }
    \enddot
    \include LinSys/solver.ci
*/
// *****************************************************************************
#ifndef Solver_h
#define Solver_h

#include <vector>
#include <map>
#include <unordered_map>

#include "Types.h"
#include "Tags.h"
#include "Fields.h"
#include "TaggedTuple.h"
#include "HypreMatrix.h"
#include "HypreVector.h"
#include "HypreSolver.h"

#include "NoWarning/solver.decl.h"
#include "NoWarning/matcg.decl.h"

namespace tk {

//! Solver shadow class used to fire off a reduction different from Solver
//! \details Solver shadow class constructor used to fire off a reduction
//!   different from Solver to avoid the runtime error "mis-matched client
//!   callbacks in reduction messages"
class SolverShadow : public CBase_SolverShadow { public: SolverShadow(); };

//! Linear system merger and solver Charm++ chare group class
//! \details Instantiations of Solver comprise a processor aware Charm++
//!   chare group. When instantiated, a new object is created on each PE and not
//!   more (as opposed to individual chares or chare array object elements). The
//!   group's elements are used to collect information from all chare objects
//!   that happen to be on a given PE. See also the Charm++ interface file
//!   solver.ci.
class Solver : public CBase_Solver {

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-parameter"
    #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    #pragma GCC diagnostic ignored "-Wunused-parameter"
  #elif defined(__INTEL_COMPILER)
    #pragma warning( push )
    #pragma warning( disable: 1478 )
  #endif
  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  Solver_SDAG_CODE
  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #elif defined(__INTEL_COMPILER)
    #pragma warning( pop )
  #endif

  public:
    //! Constructor
    Solver( CProxy_SolverShadow sh,
            const std::vector< CkCallback >& cb,
            std::size_t n,
            bool /*feedback*/ );

    //! Configure Charm++ reduction types for concatenating BC nodelists
    static void registerReducers();

    //! Receive lower and upper global node IDs all PEs will operate on
    void bounds( int p, std::size_t lower, std::size_t upper );

    //! Prepare for next step
    void next();

    //! Set number of worker chares expected to contribute on my PE
    void nchare( int n );

    //! Chares contribute their global row ids for establishing communications
    void charecom( const inciter::CProxy_MatCG& worker,
                   int fromch,
                   const std::vector< std::size_t >& row );

    //! Signal the runtime system that the workers have been created
    void created();

    //! Receive global row ids from fellow group branches
    void addrow( int fromch, int frompe, const std::set< std::size_t >& row );

    //! Acknowledge received row ids
    void recrow();

    //! Chares contribute their solution nonzero values
    void charesol( int fromch,
                   const std::vector< std::size_t >& gid,
                   const Fields& solution );

    //! Receive solution vector nonzeros from fellow group branches
    void addsol( int fromch,
                 const std::map< std::size_t,
                                 std::vector< tk::real > >& solution );

    //! Chares contribute their matrix nonzero values
    void charelhs( int fromch,
                   const std::vector< std::size_t >& gid,
                   const std::pair< std::vector< std::size_t >,
                                    std::vector< std::size_t > >& psup,
                   const tk::Fields& lhsd,
                   const tk::Fields& lhso );

    //! Receive matrix nonzeros from fellow group branches
    void addlhs( int fromch,
                 const std::map< std::size_t,
                                 std::map< std::size_t,
                                           std::vector< tk::real > > >& l );

    //! Chares contribute their rhs nonzero values
    void charerhs( int fromch,
                   const std::vector< std::size_t >& gid,
                   const Fields& r );

    //! Receive+add right-hand side vector nonzeros from fellow group branches
    void addrhs( int fromch,
                 const std::map< std::size_t, std::vector< tk::real > >& r );

    //! Chares contribute to the rhs of the low-order linear system
    void charelowrhs( int fromch,
                      const std::vector< std::size_t >& gid,
                      const Fields& lowrhs );

    //! Receive+add low-order rhs vector nonzeros from fellow group branches
    void addlowrhs( int fromch,
                    const std::map< std::size_t,
                                    std::vector< tk::real > >& lowrhs );

    //! Chares contribute to the lhs of the low-order linear system
    void charelowlhs( int fromch,
                      const std::vector< std::size_t >& gid,
                      const Fields& lowlhs );

    //! \brief Receive and add lhs vector to the low-order system from fellow
    //!   group branches
    void addlowlhs( int fromch,
                    const std::map< std::size_t,
                                    std::vector< tk::real > >& lowlhs );

    //! All communications have been establised among PEs
    void comfinal();

    //! Chares query Dirichlet boundary conditions
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    const std::unordered_map< std::size_t,
            std::vector< std::pair< bool, tk::real > > >&
      dirbc() { return m_bc; }

    //! \brief Chares contribute their global row ids and associated Dirichlet
    //!   boundary condition values at which they set BCs
    void charebc( const std::unordered_map< std::size_t,
                          std::vector< std::pair< bool, tk::real > > >& bc );

    //! Reduction target collecting the final aggregated BC node list map
    void addbc( CkReductionMsg* msg );

    //! \brief Chares contribute their numerical and analytical solutions
    //!   nonzero values for computing diagnostics
    void charediag( int fromch,
                    uint64_t it,
                    tk::real t,
                    tk::real dt,
                    const std::vector< std::size_t >& gid,
                    const Fields& u,
                    const Fields& a,
                    const std::vector< tk::real >& v );

    //! \brief Receive numerical and analytical solution vector nonzeros from
    //!   fellow group branches for computing diagnostics
    void adddiag( int fromch,
                  std::map< std::size_t,
                    std::vector< std::vector< tk::real > > >& solution );

    //! Return processing element for global mesh row id
    int pe( std::size_t gid );

  private:
    CProxy_SolverShadow m_shadow;
    //! Charm++ reduction callbacks associated to compile-time tags
    tk::tuple::tagged_tuple<
        tag::com,   CkCallback
      , tag::coord, CkCallback
      , tag::diag,  CkCallback
    > m_cb;
    std::size_t m_ncomp;       //!< Number of scalar components per unknown
    std::size_t m_nchare;      //!< Number of chares contributing to my PE
    std::size_t m_ncomm;       //!< Number of chares finished commaps on my PE
    std::size_t m_nperow;      //!< Number of fellow PEs to send row ids to
    std::size_t m_nchbc;       //!< Number of chares we received bcs from
    std::size_t m_lower;       //!< Lower index of the global rows on my PE
    std::size_t m_upper;       //!< Upper index of the global rows on my PE
    uint64_t m_it;             //!< Iteration count (original in Discretization)
    tk::real m_t;              //!< Physical time (original in Discretization)
    tk::real m_dt;             //!< Time step size (original in Discretization)
    bool m_feedback;           //!< Whether to send sub-task feedback to host
    //! Ids of workers on my PE
    std::vector< int > m_myworker;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the global row ids
    std::map< int, std::vector< std::size_t > > m_rowimport;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the solution/unknown vector
    std::map< int, std::vector< std::size_t > > m_solimport;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the left-hand side matrix
    std::map< int, std::vector< std::size_t > > m_lhsimport;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the righ-hand side vector
    std::map< int, std::vector< std::size_t > > m_rhsimport;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the low-order rhs vector
    std::map< int, std::vector< std::size_t > > m_lowrhsimport;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the low-order lhs vector
    std::map< int, std::vector< std::size_t > > m_lowlhsimport;
    //! Part of global row indices owned by my PE
    std::set< std::size_t > m_row;
    //! \brief Part of unknown/solution vector owned by my PE
    //! \details Vector of values (for each scalar equation solved) associated
    //!   to global mesh point row IDs
    std::map< std::size_t, std::vector< tk::real > > m_sol;
    //! \brief Part of left-hand side matrix owned by my PE
    //! \details Nonzero values (for each scalar equation solved) associated to
    //!   global mesh point row and column IDs.
    std::map< std::size_t,
              std::map< std::size_t, std::vector< tk::real > > > m_lhs;
    //! \brief Part of right-hand side vector owned by my PE
    //! \details Vector of values (for each scalar equation solved) associated
    //!   to global mesh point row ids
    std::map< std::size_t, std::vector< tk::real > > m_rhs;
    //! \brief Part of low-order right-hand side vector owned by my PE
    //! \details Vector of values (for each scalar equation solved) associated
    //!   to global mesh point row ids. This vector collects the rhs
    //!   terms to be combined with the rhs to produce the low-order rhs.
    std::map< std::size_t, std::vector< tk::real > > m_lowrhs;
    //! \brief Part of the low-order system left-hand side vector owned by my PE
    //! \details Vector of values (for each scalar equation solved) associated
    //!   to global mesh point row ids. This vector collects the nonzero values
    //!   of the low-order system lhs "matrix" solution.
    std::map< std::size_t, std::vector< tk::real > > m_lowlhs;
    tk::hypre::HypreVector m_x; //!< Hypre vector to store the solution/unknowns
    tk::hypre::HypreMatrix m_A; //!< Hypre matrix to store the left-hand side
    tk::hypre::HypreVector m_b; //!< Hypre vector to store the right-hand side
    //! Hypre solver
    tk::hypre::HypreSolver m_solver;
    //! Row indices for my PE
    std::vector< int > m_hypreRows;
    //! Number of matrix columns/rows on my PE
    std::vector< int > m_hypreNcols;
    //! Matrix column indices for rows on my PE
    std::vector< int > m_hypreCols;
    //! Matrix nonzero values for my PE
    std::vector< tk::real > m_hypreMat;
    //! RHS vector nonzero values for my PE
    std::vector< tk::real > m_hypreRhs;
    //! Solution vector nonzero values for my PE
    std::vector< tk::real > m_hypreSol;
    //! Global->local row id map for sending back solution vector parts
    std::map< std::size_t, std::size_t > m_lid;
    //! \brief PEs associated to lower and upper global row indices
    //! \details These are the divisions at which the linear system is divided
    //!   along PE boundaries.
    std::map< std::pair< std::size_t, std::size_t >, int > m_div;
    //! \brief PEs associated to global mesh point indices
    //! \details This is used to cache the PE associated to mesh nodes
    //!   communicated, so that a quicker-than-linear-cost search can be used to
    //!   find the PE for a communicated mesh node after the node is in the
    //!   cache.
    std::map< std::size_t, int > m_pe;
    //! \brief Values (for each scalar equation solved) of Dirichlet boundary
    //!   conditions assigned to global node IDs we set
    //! \details The map key is the global mesh node/row ID, the value is a
    //!   vector of pairs in which the bool (.first) indicates whether the
    //!   boundary condition value (.second) is set at the given node. The size
    //!   of the vectors is the number of PDEs integrated times the number of
    //!   scalar components in all PDEs. This BC map stores all row IDs at which
    //!   Dirichlet boundary conditions are prescribed, i.e., including boundary
    //!   conditions set across all PEs, not just the ones need to be set on
    //!   this PE.
    std::unordered_map< std::size_t,
                        std::vector< std::pair< bool, tk::real > > > m_bc;
    //! Matrix non-zero coefficients for Dirichlet boundary conditions
    std::unordered_map< std::size_t,
                        std::map< std::size_t, std::vector<tk::real> > > m_bca;
    //! Worker proxy
    inciter::CProxy_MatCG m_worker;

    //! Check if we have done our part in storing and exporting global row ids
    bool comcomplete() const;

    //! Check if our portion of the solution vector values is complete
    bool solcomplete() const { return m_solimport == m_rowimport; }

    //! Check if our portion of the matrix values is complete
    bool lhscomplete() const { return m_lhsimport == m_rowimport; }

    //! Check if our portion of the right-hand side vector values is complete
    bool rhscomplete() const { return m_rhsimport == m_rowimport; }

    //! Check if our portion of the low-order rhs vector values is complete
    bool lowrhscomplete() const { return m_lowrhsimport == m_rowimport; }

    //! Check if our portion of the low-order lhs vector values is complete
    bool lowlhscomplete() const { return m_lowlhsimport == m_rowimport; }

    //! Build Hypre data for our portion of the global row ids
    void hyprerow();

    //! Set boundary conditions on the left-hand side matrix
    void lhsbc();

    //! Set Dirichlet boundary conditions on the right-hand side vector
    void rhsbc();

    //! Build Hypre data for our portion of the solution vector
    void hypresol();

    //! Build Hypre data for our portion of the matrix
    void hyprelhs();

    //! Build Hypre data for our portion of the right-hand side vector
    void hyprerhs();

    //! Set our portion of values of the distributed solution vector
    void sol();

    //! Set our portion of values of the distributed matrix
    void lhs();

    //! Set our portion of values of the distributed right-hand side vector
    void rhs();

    //! Assemble distributed solution vector
    void assemblesol();

    //! Assemble distributed matrix
    void assemblelhs();

    //! Assemble distributed right-hand side vector
    void assemblerhs();

    //! Update solution vector in our PE's workers
    void updateSol();

    //! Solve hyigh-order linear system
    void solve();

    //! Update low-order solution vector in our PE's workers
    void updateLowSol();

    //! Solve low-order linear system
    void lowsolve();

    //! Update diagnostics vector
    void updateDiag( std::size_t row,
                     std::vector< tk::real >&& u,
                     std::vector< tk::real >&& a,
                     tk::real v );

    //! Compute diagnostics (residuals) and contribute them back to host
    void diagnostics();
};

} // tk::

#endif // Solver_h
