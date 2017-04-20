// *****************************************************************************
/*!
  \file      src/LinSys/LinSysMerger.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Charm++ chare linear system merger group to solve a linear system
  \details   Charm++ chare linear system merger group used to collect and
    assemble the left hand side matrix (lhs), the right hand side (rhs) vector,
    and the solution (unknown) vector from individual worker (e.g., Carrier)
    chares. Beside collection and assembly, the system is also solved. The
    solution is outsourced to hypre, an MPI-only library. Once the solution is
    available, the individual worker chares are updated with the new solution.

    In the basic configuration this class assembles and solves a single linear
    system, whose rhs vector may change during time stepping. In a more advanced
    configuration, an additional, auxiliary, linear system can also be solved.
    The basic configuration is configured by instantiating this class using the
    AuxSolver = AuxSolverNull template argument, while the one that solves the
    auxiliary system is configured by an AuxSolver class that provides a number
    of member functions. For possible auxiliary solver classes, see, e.g.,
    src/Inciter/AuxSolver.h.

    When LinSysMerger is configured to solve an auxiliary solution beside its
    primary solution, the class behind the AuxSolver template argument, must
    have static member functions only, which may provide the necessary
    functionality to collect, assemble, solve, and update an auxiliary linear
    system. The requirements on the auxiliary linear system are: (1) the left
    hand side matrix must be diagonal, (2) the right hand side vector is
    (optionally) a combination of the primary right hand side vector and another
    vector, assembled separately (and overlapped with the primary system). This
    enables configuring the auxiliary system to be, e.g., the low order solution
    as required by the flux-corrected transport algorithm used for transport
    equations in inciter. However, it also enables entirely removing the
    auxiliary solution from LinSysMerger, via inciter::AuxSolverNull. Tthe
    Charm++ SDAG logic does not change depending on whether an auxiliary
    solution is performed or not. The composition is done at compile time,
    depending on how LinSysMerger is instantiated.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/LinSys/linsysmerger.ci. Note
    that the SDAG logic is the same regardless whether an auxiliary solution is
    performed.

    #### Call graph ####
    The following is a directed acyclic graph (DAG) that outlines the
    asynchronous algorithm implemented in this class The detailed discussion of
    the algorithm is given in the Charm++ interface file linsysmerger.ci. On the
    DAG orange fills denote global synchronization points that contain or
    eventually lead to global reductions. Dashed lines are potential shortcuts
    that allow jumping over some of the task-graph under some circumstances or
    optional code paths (taken, e.g., only in DEBUG mode). See the detailed
    discussion in linsysmrger.ci.
    \dot
    digraph "LinSysMerger SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      ChRow [ label="ChRow"
              tooltip="chares contribute their global row IDs"
              URL="\ref tk::LinSysMerger::charerow"];
      ChBC [ label="ChBC"
             tooltip="chares contribute their global node IDs at which
                      they can set boundary conditions"
             URL="\ref tk::LinSysMerger::charebc"];
      RowComplete   [ label="RowComplete" color="#e6851c" style="filled"
              tooltip="all linear system merger branches have done heir part of
              storing and exporting global row ids"
              URL="\ref inciter::Transporter::rowcomplete"];
      Init [ label="Init"
              tooltip="Worker start setting and outputing ICs, computing
                       initial dt, computing LHS"];
      dt [ label="dt"
           tooltip="Worker chares compute their minimum time step size"
           color="#e6851c" style="filled"];
      Ver  [ label="Ver"
              tooltip="start optional verifications, query BCs, and converting
                       row IDs to hypre format"
              URL="\ref tk::LinSysMerger::rowsreceived"];
      LhsBC [ label="LhsBC"
              tooltip="set boundary conditions on the left-hand side matrix"
              URL="\ref tk::LinSysMerger::lhsbc"];
      RhsBC [ label="RhsBC"
              tooltip="set boundary conditions on the right-hand side vector"
              URL="\ref tk::LinSysMerger::rhsbc"];
      ChSol [ label="ChSol"
              tooltip="chares contribute their solution vector nonzeros"
              URL="\ref tk::LinSysMerger::charesol"];
      ChLhs [ label="ChLhs"
              tooltip="chares contribute their left hand side matrix nonzeros"
              URL="\ref tk::LinSysMerger::charelhs"];
      ChAuxLhs [ label="ChAuxLhs"
              tooltip="chares contribute their auxiliary left hand side matrix"
              URL="\ref tk::LinSysMerger::chareaux"];
      ChRhs [ label="ChRhs"
              tooltip="chares contribute their right hand side vector nonzeros"
              URL="\ref tk::LinSysMerger::charesol"];
      ChAuxRhs [ label="ChAuxRhs"
              tooltip="chares contribute their contribution to the auxiliary
                       right hand side vector"
              URL="\ref tk::LinSysMerger::chareauxrhs"];
      HypreRow [ label="HypreRow"
              tooltip="convert global row ID vector to hypre format"
              URL="\ref tk::LinSysMerger::hyprerow"];
      HypreSol [ label="HypreSol"
              tooltip="convert solution vector to hypre format"
              URL="\ref tk::LinSysMerger::hypresol"];
      HypreLhs [ label="HypreLhs"
              tooltip="convert left hand side matrix to hypre format"
              URL="\ref tk::LinSysMerger::hyprelhs"];
      HypreRhs [ label="HypreRhs"
              tooltip="convert right hand side vector to hypre format"
              URL="\ref tk::LinSysMerger::hyprerhs"];
      FillSol [ label="FillSol"
              tooltip="fill/set solution vector"
              URL="\ref tk::LinSysMerger::sol"];
      FillLhs [ label="FillLhs"
              tooltip="fill/set left hand side matrix"
              URL="\ref tk::LinSysMerger::lhs"];
      FillRhs [ label="FillRhs"
              tooltip="fill/set right hand side vector"
              URL="\ref tk::LinSysMerger::rhs"];
      AsmSol [ label="AsmSol"
              tooltip="assemble solution vector"
              URL="\ref tk::LinSysMerger::assemblesol"];
      AsmLhs [ label="AsmLhs"
              tooltip="assemble left hand side matrix"
              URL="\ref tk::LinSysMerger::assemblelhs"];
      AsmRhs [ label="AsmRhs"
              tooltip="assemble right hand side vector"
              URL="\ref tk::LinSysMerger::assemblerhs"];
      Solve [ label="Solve" tooltip="solve linear system"
              URL="\ref tk::LinSysMerger::solve"];
      AuxSolve [ label="AuxSolve" tooltip="solve auxiliary linear system"
              URL="\ref tk::LinSysMerger::auxsolve"];
      Upd [ label="Upd" tooltip="update solution"
                color="#e6851c" style="filled"
                URL="\ref tk::LinSysMerger::updateSol"];
      AuxUpd [ label="AuxUpd" tooltip="update auxiliary solution"
               color="#e6851c"style="filled"
               URL="\ref tk::LinSysMerger::updateAuxSol"];
      ChRow -> RowComplete [ style="solid" ];
      ChBC -> RowComplete Ver [ style="solid" ];
      RowComplete -> Init [ style="solid" ];
      RowComplete -> Ver [ style="solid" ];
      Ver -> HypreRow [ style="solid" ];
      ChLhs -> LhsBC [ style="solid" ];
      ChRhs -> RhsBC [ style="solid" ];
      LhsBC -> HypreLhs [ style="solid" ];
      RhsBC -> HypreRhs [ style="solid" ];
      RhsBC -> AuxSolve [ style="solid" ];
      ChAuxRhs -> AuxSolve [ style="solid" ];
      ChAuxLhs -> AuxSolve [ style="solid" ];
      Init -> ChSol [ style="solid" ];
      Init -> ChLhs [ style="solid" ];
      Init -> ChAuxLhs [ style="solid" ];
      Init -> dt [ style="solid" ];
      dt -> ChRhs [ style="solid" ];
      dt -> ChAuxRhs [ style="solid" ];
      HypreRow -> FillSol [ style="solid" ];
      HypreRow -> FillLhs [ style="solid" ];
      HypreRow -> FillRhs [ style="solid" ];
      ChSol -> HypreSol -> FillSol -> AsmSol -> Solve [ style="solid" ];
      HypreLhs -> FillLhs -> AsmLhs -> Solve [ style="solid" ];
      HypreRhs -> FillRhs -> AsmRhs -> Solve [ style="solid" ];
      Solve -> Upd [ style="solid" ];
      AuxSolve -> AuxUpd [ style="solid" ];
    }
    \enddot
    \include LinSys/linsysmerger.ci
*/
// *****************************************************************************
#ifndef LinSysMerger_h
#define LinSysMerger_h

#include <vector>
#include <map>
#include <unordered_map>
#include <utility>
#include <numeric>
#include <algorithm>
#include <iosfwd>
#include <cstddef>

#include "Types.h"
#include "Exception.h"
#include "ContainerUtil.h"
#include "PUPUtil.h"
#include "Fields.h"
#include "HypreMatrix.h"
#include "HypreVector.h"
#include "HypreSolver.h"
#include "VectorReducer.h"
#include "HashMapReducer.h"
#include "AuxSolver.h"

#include "NoWarning/linsysmerger.decl.h"
#include "NoWarning/transporter.decl.h"

namespace tk {

extern CkReduction::reducerType BCVectorMerger;
extern CkReduction::reducerType BCMapMerger;
extern CkReduction::reducerType BCValMerger;

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wundefined-func-template"
#endif

//! Linear system merger Charm++ chare group class
//! \details Instantiations of LinSysMerger comprise a processor aware Charm++
//!   chare group. When instantiated, a new object is created on each PE and not
//!   more (as opposed to individual chares or chare array object elements). The
//!   group's elements are used to collect information from all chare objects
//!   that happen to be on a given PE. See also the Charm++ interface file
//!   linsysmerger.ci.
//! \author J. Bakosi
template< class HostProxy, class WorkerProxy, class AuxSolver  >
class LinSysMerger : public CBase_LinSysMerger< HostProxy,
                                                WorkerProxy,
                                                AuxSolver > {

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
  LinSysMerger_SDAG_CODE
  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #elif defined(__INTEL_COMPILER)
    #pragma warning( pop )
  #endif

  private:
    using Group = CBase_LinSysMerger< HostProxy, WorkerProxy, AuxSolver >;
    using GroupIdx = CkIndex_LinSysMerger< HostProxy, WorkerProxy, AuxSolver >;

  public:
    //! Constructor
    //! \param[in] host Charm++ host proxy
    //! \param[in] worker Charm++ worker proxy
    //! \param[in] s Mesh node IDs mapped to side set ids
    //! \param[in] n Total number of scalar components in the linear system
    //! \param[in] feedback Whether to send sub-task feedback to host    
    LinSysMerger( const HostProxy& host,
                  const WorkerProxy& worker,
                  const std::map< int, std::vector< std::size_t > >& s,
                  std::size_t n,
                  bool feedback ) :
      __dep(),
      m_host( host ),
      m_worker( worker ),
      m_side( s ),
      m_ncomp( n ),
      m_nchare( 0 ),
      m_nperow( 0 ),
      m_nchbc( 0 ),
      m_lower( 0 ),
      m_upper( 0 ),
      m_feedback( feedback ),
      m_myworker(),
      m_rowimport(),
      m_solimport(),
      m_lhsimport(),
      m_rhsimport(),
      m_auxrhsimport(),
      m_auxlhsimport(),
      m_diagimport(),
      m_row(),
      m_sol(),
      m_lhs(),
      m_rhs(),
      m_auxrhs(),
      m_auxlhs(),
      m_diag(),
      m_x(),
      m_A(),
      m_b(),
      m_solver(),
      m_hypreRows(),
      m_hypreNcols(),
      m_hypreCols(),
      m_hypreMat(),
      m_hypreRhs(),
      m_hypreSol(),
      m_lid(),
      m_div(),
      m_pe(),
      m_bc()
    {
      // Activate SDAG waits
      wait4row();
      wait4lhsbc();
      wait4rhsbc();
      wait4sol();
      wait4lhs();
      wait4rhs();
      wait4hypresol();
      wait4hyprelhs();
      wait4hyprerhs();
      wait4fillsol();
      wait4filllhs();
      wait4fillrhs();
      wait4asm();
      wait4aux();
      wait4solve();
      wait4auxsolve();
    }

    //! \brief Configure Charm++ reduction types for concatenating BC nodelists
    //! \details Since this is a [nodeinit] routine, see linsysmerger.ci, the
    //!   Charm++ runtime system executes the routine exactly once on every
    //!   logical node early on in the Charm++ init sequence. Must be static as
    //!   it is called without an object. See also: Section "Initializations at
    //!   Program Startup" at in the Charm++ manual
    //!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
    static void registerBCMerger() {
      BCVectorMerger = CkReduction::addReducer( tk::mergeVector );
      BCMapMerger = CkReduction::addReducer(
                      tk::mergeHashMap< std::size_t,
                        std::vector< std::pair< bool, tk::real > > > );
      BCValMerger = CkReduction::addReducer(
                      tk::mergeHashMap< std::size_t,
                        std::vector< std::pair< bool, tk::real > > > );
    }

    //! Receive lower and upper global node IDs all PEs will operate on
    //! \param[in] p PE whose bounds being received
    //! \param[in] lower Lower index of the global rows on my PE
    //! \param[in] upper Upper index of the global rows on my PE
    void bounds( int p, std::size_t lower, std::size_t upper ) {
      Assert( lower < upper, "Lower bound must be lower than the upper bound: "
              "(" + std::to_string(lower) + "..." +  std::to_string(upper) +
              ") sent by PE " + std::to_string(p) );
      // Store our bounds
      if (p == CkMyPe()) {
        m_lower = lower;
        m_upper = upper;
      }
      // Store inverse of PE-division map stored on all PE
      m_div[ {lower,upper} ] = p;
      // If we have all PEs' bounds, signal the runtime system to continue
      if (m_div.size() == static_cast<std::size_t>(CkNumPes())) {
        // Create my PE's lhs matrix distributed across all PEs
        m_A.create( m_lower*m_ncomp, m_upper*m_ncomp );
        // Create my PE's rhs and unknown vectors distributed across all PEs
        m_b.create( m_lower*m_ncomp, m_upper*m_ncomp );
        m_x.create( m_lower*m_ncomp, m_upper*m_ncomp );
        // Create linear solver
        m_solver.create();
        // Signal back to host that setup of workers can start
        signal2host_setup( m_host );
      }
    }

    //! Re-enable SDAG waits for rebuilding the right-hand side vector only
    void enable_wait4rhs() {
      wait4rhs();
      wait4rhsbc();
      wait4hyprerhs();
      wait4fillrhs();
      wait4asm();
      wait4aux();
      wait4solve();
      wait4auxsolve();
      m_rhsimport.clear();
      m_auxrhsimport.clear();
      m_diagimport.clear();
      m_rhs.clear();
      m_auxrhs.clear();
      m_hypreRhs.clear();
      m_diag.clear();
      auxlhs_complete();
      hyprerow_complete();
      asmsol_complete();
      asmlhs_complete();
      signal2host_computedt( m_host );
    }

    //! Chares register on my PE
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void checkin() { ++m_nchare; }

    //! Chares contribute their global row ids
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] row Global mesh point (row) indices contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charerow( int fromch, const std::vector< std::size_t >& row ) {
      // Collect global ids of workers on my PE
      m_myworker.push_back( fromch );
      // Store rows owned and pack those to be exported, also build import map
      // used to test for completion
      std::map< int, std::set< std::size_t > > exp;
      for (auto gid : row) {
        if (gid >= m_lower && gid < m_upper) {  // if own
          m_rowimport[ fromch ].push_back( gid );
          m_row.insert( gid );
        } else exp[ pe(gid) ].insert( gid );
      }
      // Export non-owned parts to fellow branches that own them
      m_nperow += exp.size();
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        Group::thisProxy[ tope ].addrow( fromch, CkMyPe(), p.second );
      }
      checkifrowcomplete();
    }
    //! Receive global row ids from fellow group branches
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] frompe PE contribution coming from
    //! \param[in] row Global mesh point (row) indices received
    void addrow( int fromch, int frompe, const std::set< std::size_t >& row ) {
      for (auto r : row) {
        m_rowimport[ fromch ].push_back( r );
        m_row.insert( r );
      }
      Group::thisProxy[ frompe ].recrow();
    }
    //! Acknowledge received row ids
    void recrow() {
      --m_nperow;
      checkifrowcomplete();
    }

    //! Chares contribute their solution nonzero values
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] gid Global row indices of the vector contributed
    //! \param[in] solution Portion of the unknown/solution vector contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charesol( int fromch,
                   const std::vector< std::size_t >& gid,
                   const Fields& solution )
    {
      Assert( gid.size() == solution.nunk(),
              "Size of solution and row ID vectors must equal" );
      // Store solution vector nonzero values owned and pack those to be
      // exported, also build import map used to test for completion
      std::map< int, std::map< std::size_t, std::vector< tk::real > > > exp;
      for (std::size_t i=0; i<gid.size(); ++i)
        if (gid[i] >= m_lower && gid[i] < m_upper) {    // if own
          m_solimport[ fromch ].push_back( gid[i] );
          m_sol[ gid[i] ] = solution[i];
        } else {
          exp[ pe(gid[i]) ][ gid[i] ] = solution[i];
        }
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        Group::thisProxy[ tope ].addsol( fromch, p.second );
      }
      checkifsolcomplete();
    }
    //! Receive solution vector nonzeros from fellow group branches
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] solution Portion of the unknown/solution vector contributed,
    //!   containing global row indices and values for all components
    void addsol( int fromch,
                 const std::map< std::size_t,
                                 std::vector< tk::real > >& solution )
    {
      for (const auto& r : solution) {
        m_solimport[ fromch ].push_back( r.first );
        m_sol[ r.first ] = r.second;
      }
      checkifsolcomplete();
    }

    //! Chares contribute their matrix nonzero values
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] gid Global row indices of the matrix chunk contributed
    //! \param[in] psup Points surrounding points using local indices. See also
    //!   tk::genPsup().
    //! \param[in] lhsd Portion of the left-hand side matrix contributed,
    //!   containing non-zero values (for all scalar components of the equations
    //!   solved) as a sparse matrix diagonal
    //! \param[in] lhso Portion of the left-hand side matrix contributed,
    //!   containing non-zero values (for all scalar components of the equations
    //!   solved) as a sparse matrix off-diagonal entries in compressed row
    //!   storage format
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charelhs( int fromch,
                   const std::vector< std::size_t >& gid,
                   const std::pair< std::vector< std::size_t >,
                                    std::vector< std::size_t > >& psup,
                   const tk::Fields& lhsd,
                   const tk::Fields& lhso )
    {
      Assert( psup.second.size()-1 == gid.size(),
              "Number of mesh points and number of global IDs unequal" );
      Assert( psup.second.size()-1 == lhsd.nunk(),
              "Number of mesh points and number of diagonals unequal" );
      Assert( psup.first.size() == lhso.nunk(),
              "Number of off-diagonals and their number of indices unequal" );
      // Store matrix nonzero values owned and pack those to be exported, also
      // build import map used to test for completion
      std::map< int, std::map< std::size_t,
                               std::map< std::size_t,
                                         std::vector< tk::real > > > > exp;
      for (std::size_t i=0; i<gid.size(); ++i)
        if (gid[i] >= m_lower && gid[i] < m_upper) {  // if own
          m_lhsimport[ fromch ].push_back( gid[i] );
          auto& row = m_lhs[ gid[i] ];
          row[ gid[i] ] += lhsd[i];
          for (auto j=psup.second[i]+1; j<=psup.second[i+1]; ++j)
            row[ gid[ psup.first[j] ] ] += lhso[j];
        } else {
          auto& row = exp[ pe(gid[i]) ][ gid[i] ];
          row[ gid[i] ] = lhsd[i];
          for (auto j=psup.second[i]+1; j<=psup.second[i+1]; ++j)
            row[ gid[ psup.first[j] ] ] = lhso[j];
        }
      // Export non-owned matrix rows values to fellow branches that own them
      for (const auto& p : exp)
        Group::thisProxy[ p.first ].addlhs( fromch, p.second );
      if (lhscomplete()) lhs_complete();
    }
    //! Receive matrix nonzeros from fellow group branches
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] l Portion of the left-hand side matrix contributed,
    //!   containing global row and column indices and non-zero values
    void addlhs( int fromch,
                 const std::map< std::size_t,
                                 std::map< std::size_t,
                                           std::vector< tk::real > > >& l )
    {
      for (const auto& r : l) {
        m_lhsimport[ fromch ].push_back( r.first );
        auto& row = m_lhs[ r.first ];
        for (const auto& c : r.second) row[ c.first ] += c.second;
      }
      if (lhscomplete()) lhs_complete();
    }

    //! Chares contribute their rhs nonzero values
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] gid Global row indices of the vector contributed
    //! \param[in] r Portion of the right-hand side vector contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charerhs( int fromch,
                   const std::vector< std::size_t >& gid,
                   const Fields& r )
    {
      Assert( gid.size() == r.nunk(),
              "Size of right-hand side and row ID vectors must equal" );
      // Store+add vector of nonzero values owned and pack those to be exported
      std::map< int, std::map< std::size_t, std::vector< tk::real > > > exp;
      for (std::size_t i=0; i<gid.size(); ++i)
        if (gid[i] >= m_lower && gid[i] < m_upper) {  // if own
          m_rhsimport[ fromch ].push_back( gid[i] );
          m_rhs[ gid[i] ] += r[i];
        } else
          exp[ pe(gid[i]) ][ gid[i] ] = r[i];
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        Group::thisProxy[ tope ].addrhs( fromch, p.second );
      }
      if (rhscomplete()) rhs_complete();
    }
    //! Receive+add right-hand side vector nonzeros from fellow group branches
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] r Portion of the right-hand side vector contributed,
    //!   containing global row indices and values
    void addrhs( int fromch,
                 const std::map< std::size_t, std::vector< tk::real > >& r ) {
      for (const auto& l : r) {
        m_rhsimport[ fromch ].push_back( l.first );
        m_rhs[ l.first ] += l.second;
      }
      if (rhscomplete()) rhs_complete();
    }

    //! Chares contribute to the rhs of the auxiliary linear system
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] gid Global row indices of the vector contributed
    //! \param[in] auxrhs Portion of the auxiliary rhs vector contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void chareauxrhs( int fromch,
                      const std::vector< std::size_t >& gid,
                      const Fields& auxrhs )
    {
      AuxSolver::chareauxrhs( Group::thisProxy, this, m_lower, m_upper, fromch,
                              gid, auxrhs, m_auxrhsimport, m_auxrhs );
      if (auxrhscomplete()) auxrhs_complete();
    }
    //! Receive+add auxiliary rhs vector nonzeros from fellow group branches
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] auxrhs Portion of the auxiliary rhs vector contributed,
    //!   containing global row indices and values
    void addauxrhs( int fromch, const std::map< std::size_t,
                                   std::vector< tk::real > >& auxrhs )
    {
      AuxSolver::addauxrhs( fromch, auxrhs, m_auxrhsimport, m_auxrhs );
      if (auxrhscomplete()) auxrhs_complete();
    }

    //! Chares contribute to the lhs of the auxiliary linear system
    //! \param[in] fromch Charm chare array the contribution coming from
    //! \param[in] gid Global row indices of the vector contributed
    //! \param[in] auxlhs Portion of the auxiliary lhs vector contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void chareauxlhs( int fromch,
                      const std::vector< std::size_t >& gid,
                      const Fields& auxlhs )
    {
      AuxSolver::chareauxlhs( Group::thisProxy, this, m_lower, m_upper, fromch,
                              gid, auxlhs, m_auxlhsimport, m_auxlhs );
      if (auxlhscomplete()) auxlhs_complete();
    }
    //! \brief Receive and add lhs vector to the auxiliary system from fellow
    //!   group branches
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] auxlhs Portion of the lhs vector contributed to the auxiliary
    //!   linear system, containing global row indices and values
    void addauxlhs( int fromch, const std::map< std::size_t,
                                  std::vector< tk::real > >& auxlhs )
    {
      AuxSolver::addauxlhs( fromch, auxlhs, m_auxlhsimport, m_auxlhs );
      if (auxlhscomplete()) auxlhs_complete();
    }

    //! Assert that all global row indices have been received on my PE
    //! \details The assert consists of three necessary conditions, which
    //!   together comprise the sufficient condition that all global row indices
    //!   have been received owned by this PE.
    void rowsreceived() {
      Assert( // 1. have heard from every chare on my PE
              m_myworker.size() == m_nchare &&
              // 2. number of rows equals that of the expected on my PE
              m_row.size() == m_upper-m_lower &&
              // 3. all fellow PEs have received my row ids contribution
              m_nperow == 0,
              // if any of the above is unsatisfied, the row ids are incomplete
              "Row ids are incomplete on PE " + std::to_string(CkMyPe()) );
      // now that the global row ids are complete, build Hypre data from it
      hyprerow();
    }

    //! Chares query side set info
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    const std::map< int, std::vector< std::size_t > >& side() { return m_side; }

    //! Chares query Dirichlet boundary conditions
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    const std::unordered_map< std::size_t,
            std::vector< std::pair< bool, tk::real > > >&
      dirbc() { return m_bc; }

    //! \brief Chares contribute their global row ids and associated Dirichlet
    //!   boundary condition values at which they set BCs
    //! \param[in] bc Vector of pairs of bool and BC value (BC vector)
    //!   associated to global node IDs at which the boundary condition is set.
    //!   Here the bool indicates whether the BC value is set at the given node
    //!   by the user. The size of the vectors is the number of PDEs integrated
    //!   times the number of scalar components in all PDEs.
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charebc( const std::unordered_map< std::size_t,
                          std::vector< std::pair< bool, tk::real > > >& bc )
    {
      // Associate BC vectors to mesh nodes owned
      if (!bc.empty())  // only if chare has anything to offer
        for (const auto& n : bc) {
          Assert( n.second.size() == m_ncomp, "The total number of scalar "
          "components does not equal that of set in the BC vector." );
          m_bc[ n.first ] = n.second;
        }
      // Forward all BC vectors received to fellow branches
      if (++m_nchbc == m_nchare) {
        auto stream = tk::serialize( m_bc );
        Group::contribute( stream.first, stream.second.get(), BCMapMerger,
                      CkCallback(GroupIdx::addbc(nullptr),Group::thisProxy) );
      }
    }
    // Reduction target collecting the final aggregated BC node list map
    void addbc( CkReductionMsg* msg ) {
      PUP::fromMem creator( msg->getData() );
      creator | m_bc;
      delete msg;
      if (m_feedback) m_host.pebccomplete();    // send progress report to host
      bc_complete();
    }

    //! \brief Chares contribute their numerical and analytical solutions
    //!   nonzero values for computing diagnostics
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] gid Global row indices of the vector contributed
    //! \param[in] num Portion of the numerical unknown/solution vector
    //! \param[in] an Portion of the analytical solution vector
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charediag( int fromch,
                    const std::vector< std::size_t >& gid,
                    const Fields& num,
                    const Fields& an )
    {
      Assert( gid.size() == num.nunk(),
              "Size of numerical solution and row ID vectors must equal" );
      // Store numerical and analytical solution vector nonzero values owned and
      // pack those to be exported, also build import map used to test for
      // completion
      std::map< int, std::map< std::size_t,
                       std::vector< std::vector< tk::real > > > > exp;
      for (std::size_t i=0; i<gid.size(); ++i)
        if (gid[i] >= m_lower && gid[i] < m_upper) {    // if own
          m_diagimport[ fromch ].push_back( gid[i] );
          m_diag[ gid[i] ] = {{ num[i], an[i] }};
        } else
          exp[ pe(gid[i]) ][ gid[i] ] = {{ num[i], an[i] }};
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        Group::thisProxy[ tope ].adddiag( fromch, p.second );
      }
      if (diagcomplete()) diagnostics();
    }
    //! \brief Receive numerical and analytical solution vector nonzeros from
    //!   fellow group branches for computing diagnostics
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] solution Portion of the solution vector contributed,
    //!   containing global row indices and values for all components of the
    //!   numerical and analytical solutions (if available)
    void adddiag( int fromch,
                  const std::map< std::size_t,
                          std::vector< std::vector< tk::real > > >& solution )
    {
      for (const auto& r : solution) {
        m_diagimport[ fromch ].push_back( r.first );
        m_diag[ r.first ] = r.second;
      }
      if (diagcomplete()) diagnostics();
    }

    //! Check if our portion of the right-hand side vector values is complete
    //! \return True if all parts of the right-hand side vector have been
    //!   received
    bool rhscomplete() const { return m_rhsimport == m_rowimport; }
    //! Check if our portion of the auxiliary rhs vector values is complete
    //! \return True if all parts of the auxiliary rhs vector have been received
    bool auxrhscomplete() const { return m_auxrhsimport == m_rowimport; }
    //! Check if our portion of the auxiliary lhs vector values is complete
    //! \return True if all parts of the auxiliary lhs vector have been received
    bool auxlhscomplete() const { return m_auxlhsimport == m_rowimport; }

    //! Return processing element for global mesh row id
    //! \param[in] gid Global mesh point (matrix or vector row) id
    //! \details First we attempt to the point index in the cache. If that
    //!   fails, we resort to a linear search across the division map. Once the
    //!   PE is found, we store it in the cache, so next time the search is
    //!   quicker. This procedure must find the PE for the id.
    //! \return PE that owns global row id
    int pe( std::size_t gid ) {
      int p = -1;
      auto it = m_pe.find( gid );
      if (it != end(m_pe))
        p = it->second;
      else
        for (const auto& d : m_div)
          if (gid >= d.first.first && gid < d.first.second)
            p = m_pe[ gid ] = d.second;
      Assert( p >= 0, "PE not found for node id " + std::to_string(gid) );
      Assert( p < CkNumPes(), "Assigning to nonexistent PE" );
      return p;
    }

  private:
    HostProxy m_host;           //!< Host proxy
    WorkerProxy m_worker;       //!< Worker proxy
    //! Global (as in file) mesh node IDs mapped to side set ids of the mesh
    std::map< int, std::vector< std::size_t > > m_side;
    std::size_t m_ncomp;        //!< Number of scalar components per unknown
    std::size_t m_nchare;       //!< Number of chares contributing to my PE
    std::size_t m_nperow;       //!< Number of fellow PEs to send row ids to
    std::size_t m_nchbc;        //!< Number of chares we received bcs from
    std::size_t m_lower;        //!< Lower index of the global rows on my PE
    std::size_t m_upper;        //!< Upper index of the global rows on my PE
    bool m_feedback;            //!< Whether to send sub-task feedback to host
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
    //!   id during the communication of the auxiliary rhs vector
    std::map< int, std::vector< std::size_t > > m_auxrhsimport;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the auxiliary lhs vector
    std::map< int, std::vector< std::size_t > > m_auxlhsimport;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the solution/unknown vector for
    //!   computing diagnostics, e.g., residuals
    std::map< int, std::vector< std::size_t > > m_diagimport;
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
    //! \brief Part of auxiliary right-hand side vector owned by my PE
    //! \details Vector of values (for each scalar equation solved) associated
    //!   to global mesh point row ids. This vector collects the rhs
    //!   terms to be combined with the rhs to produce the auxiliary rhs.
    std::map< std::size_t, std::vector< tk::real > > m_auxrhs;
    //! \brief Part of the auxiliary system left-hand side vector owned by my PE
    //! \details Vector of values (for each scalar equation solved) associated
    //!   to global mesh point row ids. This vector collects the nonzero values
    //!   of the auxiliary system lhs "matrix" solution.
    std::map< std::size_t, std::vector< tk::real > > m_auxlhs;
    //! \brief Part of solution vector owned by my PE
    //! \details Vector of values (for each scalar equation solved) associated
    //!   to global mesh point row IDs. The map-value is a pair of vectors,
    //!   where the first one is the numerical and the second one is the
    //!   analytical solution. If the analytical solution for a PDE is not
    //!   defined, it is the initial condition.
    //! \see charediag(), adddiag(), e.g., inciter::Carrier::diagnostics()
    std::map< std::size_t, std::vector< std::vector< tk::real > > > m_diag;
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

    //! Check if we have done our part in storing and exporting global row ids
    //! \details This does not mean the global row ids on our PE is complete
    //!   (which is tested by an assert in rowsreceived), only that we have done
    //!   our part of receiving contributions from chare array groups storing
    //!   the parts that we own and have sent the parts we do not own to fellow
    //!   PEs, i.e., we have nothing else to export. Only when all other fellow
    //!   branches have received all contributions are the row ids complete on
    //!   all PEs. This latter condition can only be tested after the global
    //!   reduction initiated by signal2host_row_complete, which is called when
    //!   all fellow branches have returned true from rowcomplete.
    //! \see rowsreceived()
    //! \return True if we have done our part storing and exporting row ids
    bool rowcomplete() const {
      return // have heard from every chare on my PE
             m_myworker.size() == m_nchare &&
             // all fellow PEs have received my row ids contribution
             m_nperow == 0;
    }
    //! Check if our portion of the solution vector values is complete
    //! \return True if all parts of the unknown/solution vector have been
    //!   received
    bool solcomplete() const { return m_solimport == m_rowimport; }
    //! Check if our portion of the matrix values is complete
    //! \return True if all parts of the left-hand side matrix have been
    //!   received
    bool lhscomplete() const { return m_lhsimport == m_rowimport; }
    //! \brief Check if our portion of the solution vector values (for
    //!   diagnostics) is complete
    //! \return True if all parts of the unknown/solution vector have been
    //!   received
    bool diagcomplete() const { return m_diagimport == m_rowimport; }

    //! Check if contributions to global row IDs are complete
    //! \details If so, send progress report to host that this sub-task is done,
    //!   and tell the runtime system that this is complete.
    void checkifrowcomplete() {
      if (rowcomplete()) {
        if (m_feedback) m_host.perowcomplete();
        row_complete();
      };
    }
    //! Check if contributions to unknown/solution vector are complete
    //! \details If so, send progress report to host that this sub-task is done,
    //!   and tell the runtime system that this is complete.
    void checkifsolcomplete() {
      if (solcomplete()) {
        if (m_feedback) m_host.pesolcomplete();
        sol_complete();
      }
    }

    //! Build Hypre data for our portion of the global row ids
    //! \note Hypre only likes one-based indexing. Zero-based row indexing fails
    //!   to update the vector with HYPRE_IJVectorGetValues().
    void hyprerow() {
      for (auto r : m_row) {
        std::vector< int > h( m_ncomp );
        std::iota( begin(h), end(h), r*m_ncomp+1 );
        m_hypreRows.insert( end(m_hypreRows), begin(h), end(h) );
      }
      hyprerow_complete(); hyprerow_complete(); hyprerow_complete();
    }

    //! Set boundary conditions on the left-hand side matrix
    void lhsbc() {
      Assert( lhscomplete(),
              "Nonzero values of distributed matrix on PE " +
              std::to_string( CkMyPe() ) + " is incomplete: cannot set BCs" );

      for (auto& r : m_lhs) {
        auto it = m_bc.find( r.first );
        if (it != end(m_bc)) {
          auto& diag = tk::ref_find( r.second, r.first );
          for (std::size_t i=0; i<m_ncomp; ++i)
            if (it->second[i].first) {
              // zero columns in BC row
              for (auto& col : r.second) col.second[i] = 0.0;
              // put 1.0 in diagonal of BC row
              diag[i] = 1.0;
            }
        }
      }
      lhsbc_complete();
    }

    //! Set Dirichlet boundary conditions on the right-hand side vector
    //! \details Since we solve for the solution increment, this amounts to
    //!    enforcing zero rhs (no solution increment) at BC nodes
    void rhsbc() {
      Assert( rhscomplete(),
              "Values of distributed right-hand-side vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete: cannot set BCs" );
      for (const auto& n : m_bc) {
        if (n.first >= m_lower && n.first < m_upper) {
          auto& r = tk::ref_find( m_rhs, n.first );
          for (std::size_t i=0; i<m_ncomp; ++i)
            if (n.second[i].first)
              r[i] = 0.0;
        }
      }
      rhsbc_complete(); rhsbc_complete();
    }

    //! Build Hypre data for our portion of the solution vector
    void hypresol() {
      Assert( solcomplete(),
              "Values of distributed solution vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete" );
      std::size_t i = 0;
      for (const auto& r : m_sol) {
        m_lid[ r.first ] = i++;
        m_hypreSol.insert( end(m_hypreSol), begin(r.second), end(r.second) );
      }
      hypresol_complete();
    }
    //! Build Hypre data for our portion of the matrix
    //! \note Hypre only likes one-based indexing. Zero-based row indexing fails
    //!   to update the vector with HYPRE_IJVectorGetValues().
    void hyprelhs() {
      Assert( lhscomplete(),
              "Nonzero values of distributed matrix on PE " +
              std::to_string( CkMyPe() ) + " is incomplete: cannot convert" );
      for (const auto& r : m_lhs) {
        for (std::size_t i=0; i<m_ncomp; ++i) {
          m_hypreNcols.push_back( static_cast< int >( r.second.size() ) );
          for (const auto& c : r.second) {
            m_hypreCols.push_back( static_cast< int >( c.first*m_ncomp+i+1 ) );
            m_hypreMat.push_back( c.second[i] );
          }
        }
      }
      hyprelhs_complete();
    }
    //! Build Hypre data for our portion of the right-hand side vector
    void hyprerhs() {
      Assert( rhscomplete(),
              "Values of distributed right-hand-side vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete: cannot convert" );
      for (const auto& r : m_rhs)
        m_hypreRhs.insert( end(m_hypreRhs), begin(r.second), end(r.second) );
      hyprerhs_complete();
    }

    //! Set our portion of values of the distributed solution vector
    void sol() {
      Assert( m_hypreSol.size() == m_hypreRows.size(), "Solution vector values "
              "incomplete on PE " + std::to_string(CkMyPe()) );
      // Set our portion of the vector values
      m_x.set( static_cast< int >( (m_upper - m_lower)*m_ncomp ),
               m_hypreRows.data(),
               m_hypreSol.data() );
      fillsol_complete();
    }
    //! Set our portion of values of the distributed matrix
    void lhs() {
      Assert( m_hypreMat.size() == m_hypreCols.size(), "Matrix values "
              "incomplete on PE " + std::to_string(CkMyPe()) );
      // Set our portion of the matrix values
      m_A.set( static_cast< int >( (m_upper - m_lower)*m_ncomp ),
               m_hypreNcols.data(),
               m_hypreRows.data(),
               m_hypreCols.data(),
               m_hypreMat.data() );
      filllhs_complete();
    }
    //! Set our portion of values of the distributed right-hand side vector
    void rhs() {
      Assert( m_hypreRhs.size() == m_hypreRows.size(), "RHS vector values "
              "incomplete on PE " + std::to_string(CkMyPe()) );
      // Set our portion of the vector values
      m_b.set( static_cast< int >( (m_upper - m_lower)*m_ncomp ),
               m_hypreRows.data(),
               m_hypreRhs.data() );
      fillrhs_complete();
    }

    //! Assemble distributed solution vector
    void assemblesol() {
      m_x.assemble();
      asmsol_complete();
    }
    //! Assemble distributed matrix
    void assemblelhs() {
      m_A.assemble();
      asmlhs_complete();
    }
    //! Assemble distributed right-hand side vector
    void assemblerhs() {
      m_b.assemble();
      asmrhs_complete();
    }

    //! Update solution vector in our PE's workers
    void updateSol() {
      // Get solution vector values for our PE
      m_x.get( static_cast< int >( (m_upper - m_lower)*m_ncomp ),
               m_hypreRows.data(),
               m_hypreSol.data() );
      // Group solution vector by workers and send each the parts back to
      // workers that own them
      for (const auto& w : m_solimport) {
        std::vector< std::size_t > gid;
        std::vector< tk::real > solution;
        for (auto r : w.second) {
          const auto it = m_sol.find( r );
          if (it != end(m_sol)) {
            gid.push_back( it->first );
            auto i = tk::cref_find( m_lid, it->first );
            using diff_type = typename decltype(m_hypreSol)::difference_type;
            auto b = static_cast< diff_type >( i*m_ncomp );
            auto e = static_cast< diff_type >( (i+1)*m_ncomp );
            solution.insert( end(solution),
                             std::next( begin(m_hypreSol), b ),
                             std::next( begin(m_hypreSol), e ) );
          } else
            Throw( "Can't find global row id " + std::to_string(r) +
                   " to export in solution vector" );
        }
        m_worker[ w.first ].updateSol( gid, solution );
      }
    }

    //! Solve linear system
    void solve() {
      m_solver.solve( m_A, m_b, m_x );
      // send progress report to host
      if (m_feedback) m_host.pesolve();
      solve_complete();
    }

    //! Update auxiliary solution vector in our PE's workers
    void updateAuxSol() {
      AuxSolver::update( m_worker, m_solimport, m_auxrhs );
    }

    //! Solve auxiliary linear system
    void auxsolve() {
      // Set boundary conditions on the auxiliary right hand side vector
      for (const auto& n : m_bc)
        if (n.first >= m_lower && n.first < m_upper) {
          auto& r = tk::ref_find( m_auxrhs, n.first );
          for (std::size_t i=0; i<m_ncomp; ++i)
            if (n.second[i].first)
              r[i] = 0.0;
        }
      // Solve auxiliary system
      AuxSolver::solve( this, m_ncomp, m_rhs, m_auxlhs, m_auxrhs );
      auxsolve_complete();
    }

    //! Compute diagnostics (residuals) and contribute them back to host
    //! \details Diagnostics: L1 norm for all components
    void diagnostics() {
      Assert( diagcomplete(),
              "Values of distributed solution vector (for diagnostics) on PE " +
              std::to_string( CkMyPe() ) + " is incomplete" );
      std::vector< tk::real > diag( m_ncomp*2, 0.0 );
      for (const auto& s : m_diag) {
        Assert( s.second.size() == 2, "Size of diagnostics vector must be 2" );
        // Compute L1 norm of the numerical solution
        for (std::size_t c=0; c<m_ncomp; ++c)
          diag[c] += std::abs( s.second[0][c] );
        // Compute L1 norm of the numerical - analytical solution
        for (std::size_t c=0; c<m_ncomp; ++c)
          diag[m_ncomp+c] += std::abs( s.second[0][c] - s.second[1][c] );
      }
      // Contribute to diagnostics across all PEs
      signal2host_diag( m_host, diag );
    }

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wdocumentation"
    #endif
    /** @name Host signal calls
      * \brief These functions signal back to the host via a global reduction
      *   originating from each PE branch
      * \details Singal calls contribute to a reduction on all branches (PEs)
      *   of LinSysMerger to the host, e.g., inciter::CProxy_Transporter, given
      *   by the template argument HostProxy. The signal functions are overloads
      *   on the specialization, e.g., inciter::CProxy_Transporter, of the
      *   LinSysMerger template. They create Charm++ reduction targets via
      *   creating a callback that invokes the typed reduction client, where
      *   host is the proxy on which the reduction target method, given by the
      *   string followed by "redn_wrapper_", e.g., rowcomplete(), is called
      *   upon completion of the reduction.
      *
      *   Note that we do not use Charm++'s CkReductionTarget macro here,
      *   but instead explicitly generate the code that that macro would
      *   generate. To explain why, here is Charm++'s CkReductionTarget macro's
      *   definition, given in ckreduction.h:
      *   \code{.cpp}
      *      #define CkReductionTarget(me, method) \
      *        CkIndex_##me::redn_wrapper_##method(NULL)
      *   \endcode
      *   This macro takes arguments 'me' (a class name) and 'method' a member
      *   function of class 'me' and generates the call
      *   'CkIndex_<class>::redn_wrapper_<method>(NULL)'. With the overloads the
      *   signal2* functions generate, we do the above macro's job for
      *   LinSysMerger specialized by HostProxy, hard-coded here, as well its
      *   reduction target. This is required since
      *    * Charm++'s CkReductionTarget macro's preprocessing happens earlier
      *      than type resolution and the string of the template argument would
      *      be substituted instead of the type specialized (which is not what
      *      we want here), and
      *    * the template argument class, e.g, CProxy_Transporter, is in a
      *      namespace different than that of LinSysMerger. When a new class is
      *      used to specialize LinSysMerger, the compiler will alert that a new
      *      overload needs to be defined.
      *
      * \note This simplifies client-code, e.g., inciter::Transporter, which now
      *   requires no explicit book-keeping with counters, etc. Also a reduction
      *   (instead of a direct call to the host) better utilizes the
      *   communication network as computational nodes can send their aggregated
      *   contribution to other nodes on a network instead of all chares sending
      *   their (smaller) contributions to the same host, (hopefully)
      *   implemented using a tree among the PEs.
      * \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html,
      *   Sections "Processor-Aware Chare Collections" and "Chare Arrays".
      * */
    ///@{
    //! \brief Signal back to host that the initialization of the row indices of
    //!   the linear system is complete
    void signal2host_row_complete( const inciter::CProxy_Transporter& host ) {
      using inciter::CkIndex_Transporter;
      Group::contribute(
        CkCallback( CkIndex_Transporter::redn_wrapper_rowcomplete(NULL), host ) );
    }
    //! \brief Signal back to host that enabling the SDAG waits for assembling
    //!    the right-hand side is complete and ready for a new advance in time
    void signal2host_computedt( const inciter::CProxy_Transporter& host ) {
      using inciter::CkIndex_Transporter;
      Group::contribute(
       CkCallback( CkIndex_Transporter::redn_wrapper_computedt(NULL), host ) );
    }
    //! \brief Signal back to host that receiving the inverse PE-division map is
    //!  complete and we are ready for Prformers to start their setup.
    void signal2host_setup( const inciter::CProxy_Transporter& host ) {
      using inciter::CkIndex_Transporter;
      Group::contribute(
       CkCallback( CkIndex_Transporter::redn_wrapper_setup(NULL), host ) );
    }
    //! Contribute diagnostics back to host
    void signal2host_diag( const inciter::CProxy_Transporter& host,
                           const std::vector< tk::real >& diag ) {
      using inciter::CkIndex_Transporter;
      Group::contribute( static_cast< int >( diag.size() * sizeof(tk::real) ),
                         diag.data(), CkReduction::sum_double,
        CkCallback( CkReductionTarget( Transporter, diagnostics), host ) );
    }
    ///@}
    #if defined(__clang__)
      #pragma GCC diagnostic pop
    #endif
};

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

} // tk::

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wreorder"
  #pragma clang diagnostic ignored "-Wunused-variable"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
  #pragma GCC diagnostic ignored "-Wreorder"
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wunused-parameter"
  #pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#define CK_TEMPLATES_ONLY
#include "linsysmerger.def.h"
#undef CK_TEMPLATES_ONLY

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // LinSysMerger_h
