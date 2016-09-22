// *****************************************************************************
/*!
  \file      src/LinSys/LinSysMerger.h
  \author    J. Bakosi
  \date      Mon 19 Sep 2016 02:33:51 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Charm++ chare linear system merger group to solve a linear system
  \details   Charm++ chare linear system merger group used to collect and
    assemble the left hand side matrix, the right hand side vector, and the
    solution (unknown) vector from individual worker (e.g., Carrier) chares. The
    solution is outsourced to hypre, an MPI-only library.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/LinSys/linsysmerger.ci.

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
      ChBCs [ label="ChBCs"
              tooltip="chares contribute their global node IDs at which
                       they can set boundary conditions"
              URL="\ref tk::LinSysMerger::charebc"];
      RowComplete   [ label="RowComplete" color="#e6851c" style="filled"
              tooltip="all linear system merger branches have done heir part of
              storing and exporting global row ids"
              URL="\ref inciter::Transporter::rowcomplete"];
      Init [ label="Init"
              tooltip="Carrier chares start initialization (setting ICs, sending
              unknown/solution vectors for assembly to LinSysMerger, start
              computing LHS, RHS, sending both for assembly)"
              URL="\ref inciter::Carrier::init"];
      Ver  [ label="Ver"
              tooltip="start optional verifications, query BCs, and converting
                       row IDs to hypre format"
              URL="\ref tk::LinSysMerger::rowsreceived"];
      VerRow [ label="VerRow"
              tooltip="optional verification ensuring consistent row IDs"
              URL="\ref tk::LinSysMerger::rowsreceived"];
      VerBCs [ label="VerBCs" color="#e6851c" style="filled"
              tooltip="optional verification ensuring consistent BC node IDs"
              URL="\ref tk::LinSysMerger::verifybc"];
      QueryBCVal [ label="QueryBCVal"
              tooltip="query boundary condition values from Carrier chares"
              URL="\ref tk::LinSysMerger::querybcval"];
      LhsBC [ label="LhsBC"
              tooltip="set boundary conditions on the left-hand side matrix"
              URL="\ref tk::LinSysMerger::lhsbc"];
      RhsBC [ label="RhsBC"
              tooltip="set boundary conditions on the right-hand side vector"
              URL="\ref tk::LinSysMerger::rhsbc"];
      LumpBC [ label="LumpBC"
              tooltip="set boundary conditions on the lumped mass left-hand side
              matrix" URL="\ref tk::LinSysMerger::lumpbc"];
      ChSol [ label="ChSol"
              tooltip="chares contribute their solution vector nonzeros"
              URL="\ref tk::LinSysMerger::charesol"];
      ChLhs [ label="ChLhs"
              tooltip="chares contribute their left hand side matrix nonzeros"
              URL="\ref tk::LinSysMerger::charelhs"];
      ChLump [ label="ChLump"
              tooltip="chares contribute their lumped mass left hand side matrix"
              URL="\ref tk::LinSysMerger::charelump"];
      ChRhs [ label="ChRhs"
              tooltip="chares contribute their right hand side vector nonzeros"
              URL="\ref tk::LinSysMerger::charesol"];
      ChDiff [ label="ChDiff"
              tooltip="chares contribute their mass diffusion contribution to
                       the right hand side vector"
              URL="\ref tk::LinSysMerger::charediff"];
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
              tooltip="fill/set lefth hand side matrix"
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
      Solve [ label="Solve" tooltip="solve high order linear system"
              URL="\ref tk::LinSysMerger::solve"];
      LowSolve [ label="LowSolve" tooltip="solve low order linear system"
              URL="\ref tk::LinSysMerger::lowsolve"];
      HighUpd [ label="HighUpd" tooltip="update high order solution"
                color="#e6851c" style="filled"
                URL="\ref tk::LinSysMerger::updateHighSol"];
      LowUpd [ label="LowUpd" tooltip="update low order solution"
               color="#e6851c"style="filled"
               URL="\ref tk::LinSysMerger::updateLowSol"];
      ChRow -> RowComplete [ style="solid" ];
      ChBCs -> RowComplete Ver [ style="solid" ];
      RowComplete -> Init [ style="solid" ];
      RowComplete -> Ver [ style="solid" ];
      Ver -> HypreRow [ style="solid" ];
      Ver -> VerRow [ style="dashed" ];
      Ver -> VerBCs [ style="dashed" ];
      Ver -> QueryBCVal [ style="solid" ];
      VerBCs -> HighUpd [ style="dashed" ];
      VerBCs -> LowUpd [ style="dashed" ];
      QueryBCVal -> LhsBC [ style="solid" ];
      QueryBCVal -> RhsBC [ style="solid" ];
      QueryBCVal -> LumpBC [ style="solid" ];
      ChLhs -> LhsBC [ style="solid" ];
      ChRhs -> RhsBC [ style="solid" ];
      ChLump -> LumpBC [ style="solid" ];
      LhsBC -> HypreLhs [ style="solid" ];
      RhsBC -> HypreRhs [ style="solid" ];
      LumpBC -> LowSolve [ style="solid" ];
      ChRhs -> LowSolve [ style="solid" ];
      ChDiff -> LowSolve [ style="solid" ];
      Init -> ChSol [ style="solid" ];
      Init -> ChLhs [ style="solid" ];
      Init -> ChLump [ style="solid" ];
      Init -> ChRhs [ style="solid" ];
      Init -> ChDiff [ style="solid" ];
      HypreRow -> FillSol [ style="solid" ];
      HypreRow -> FillLhs [ style="solid" ];
      HypreRow -> FillRhs [ style="solid" ];
      ChSol -> HypreSol -> FillSol -> AsmSol -> Solve [ style="solid" ];
      HypreLhs -> FillLhs -> AsmLhs -> Solve [ style="solid" ];
      HypreRhs -> FillRhs -> AsmRhs -> Solve [ style="solid" ];
      Solve -> HighUpd [ style="solid" ];
      LowSolve -> LowUpd [ style="solid" ];
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
template< class HostProxy, class WorkerProxy  >
class LinSysMerger : public CBase_LinSysMerger< HostProxy, WorkerProxy > {

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-parameter"
    #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #elif defined(__GNUC__)
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
  #elif defined(__GNUC__)
    #pragma GCC diagnostic pop
  #elif defined(__INTEL_COMPILER)
    #pragma warning( pop )
  #endif

  private:
    using Group = CBase_LinSysMerger< HostProxy, WorkerProxy >;

  public:
    //! Constructor
    //! \param[in] host Charm++ host proxy
    //! \param[in] worker Charm++ worker proxy
    //! \param[in] ncomp Total number of scalar components in the linear system
    LinSysMerger( const HostProxy& host,
                  const WorkerProxy& worker,
                  const std::map< int, std::vector< std::size_t > >& side,
                  std::size_t ncomp ) :
      __dep(),
      m_host( host ),
      m_worker( worker ),
      m_side( side ),
      m_ncomp( ncomp ),
      m_nchare( 0 ),
      m_nperow( 0 ),
      m_nchbc( 0 ),
      m_nchbcval( 0 ),
      m_nchdiag( 0 ),
      m_lower( 0 ),
      m_upper( 0 ),
      m_myworker(),
      m_rowimport(),
      m_solimport(),
      m_lhsimport(),
      m_rhsimport(),
      m_diffimport(),
      m_lumpimport(),
      m_row(),
      m_sol(),
      m_lhs(),
      m_rhs(),
      m_diff(),
      m_lump(),
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
      m_bc(),
      m_oldbc(),
      m_bcval()
    {
      // Activate SDAG waits
      wait4row();
      wait4lhsbc();
      wait4rhsbc();
      wait4lumpbc();
      wait4sol();
      wait4lhs();
      wait4rhs();
      wait4loworder();
      wait4hypresol();
      wait4hyprelhs();
      wait4hyprerhs();
      wait4fillsol();
      wait4filllhs();
      wait4fillrhs();
      wait4asm();
      wait4solve();
      wait4lowsolve();
      #ifdef NDEBUG     // skip verification of BCs in RELEASE mode
      ver_complete(); ver_complete();
      #endif
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
                      tk::mergeHashMap< int, std::vector< std::size_t > > );
      BCValMerger = CkReduction::addReducer(
                      tk::mergeHashMap< std::size_t,
                        std::vector< std::pair< bool, tk::real > > > );
    }

    //! Receive lower and upper global node IDs all PEs will operate on
    //! \param[in] pe PE whose bounds being received
    //! \param[in] lower Lower index of the global rows on my PE
    //! \param[in] upper Upper index of the global rows on my PE
    void bounds( int pe, std::size_t lower, std::size_t upper ) {
      // Store our bounds
      if (pe == CkMyPe()) {
        m_lower = lower;
        m_upper = upper;
      }
      // Store inverse of PE-division map stored on all PE
      m_div[ {lower,upper} ] = pe;
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
      wait4loworder();
      wait4solve();
      wait4lowsolve();
      m_rhsimport.clear();
      m_diffimport.clear();
      m_rhs.clear();
      m_diff.clear();
      m_hypreRhs.clear();
      m_bcval.clear();
      hyprerow_complete();
      asmsol_complete();
      asmlhs_complete();
      ver_complete(); ver_complete();
      lumpbc_complete();
      signal2host_advance( m_host );
      querybcval();
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
      if (rowcomplete()) row_complete();
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
      if (rowcomplete()) row_complete();
    }

    //! Chares contribute their solution nonzero values
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] gid Global row indices of the vector contributed
    //! \param[in] sol Portion of the unknown/solution vector contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charesol( int fromch,
                   const std::vector< std::size_t >& gid,
                   const Fields& sol )
    {
      Assert( gid.size() == sol.nunk(),
              "Size of solution and row ID vectors must equal" );
      // Store solution vector nonzero values owned and pack those to be
      // exported, also build import map used to test for completion
      std::map< int, std::map< std::size_t, std::vector< tk::real > > > exp;
      for (std::size_t i=0; i<gid.size(); ++i)
        if (gid[i] >= m_lower && gid[i] < m_upper) {    // if own
          m_solimport[ fromch ].push_back( gid[i] );
          m_sol[ gid[i] ] = sol[i];
        } else {
          exp[ pe(gid[i]) ][ gid[i] ] = sol[i];
        }
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        Group::thisProxy[ tope ].addsol( fromch, p.second );
      }
      if (solcomplete()) sol_complete();
    }
    //! Receive solution vector nonzeros from fellow group branches
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] sol Portion of the unknown/solution vector contributed,
    //!   containing global row indices and values for all components
    void addsol( int fromch,
                 const std::map< std::size_t, std::vector< tk::real > >& sol )
    {
      for (const auto& r : sol) {
        m_solimport[ fromch ].push_back( r.first );
        m_sol[ r.first ] = r.second;
      }
      if (solcomplete()) sol_complete();
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
    //! \param[in] lhs Portion of the left-hand side matrix contributed,
    //!   containing global row and column indices and non-zero values
    void addlhs( int fromch,
                 const std::map< std::size_t,
                                 std::map< std::size_t,
                                           std::vector< tk::real > > >& lhs )
    {
      for (const auto& r : lhs) {
        m_lhsimport[ fromch ].push_back( r.first );
        auto& row = m_lhs[ r.first ];
        for (const auto& c : r.second) row[ c.first ] += c.second;
      }
      if (lhscomplete()) lhs_complete();
    }

    //! Chares contribute their rhs nonzero values
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] gid Global row indices of the vector contributed
    //! \param[in] rhs Portion of the right-hand side vector contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charerhs( int fromch,
                   const std::vector< std::size_t >& gid,
                   const Fields& rhs )
    {
      Assert( gid.size() == rhs.nunk(),
              "Size of right-hand side and row ID vectors must equal" );
      // Store+add vector of nonzero values owned and pack those to be exported
      std::map< int, std::map< std::size_t, std::vector< tk::real > > > exp;
      for (std::size_t i=0; i<gid.size(); ++i)
        if (gid[i] >= m_lower && gid[i] < m_upper) {  // if own
          m_rhsimport[ fromch ].push_back( gid[i] );
          m_rhs[ gid[i] ] += rhs[i];
        } else
          exp[ pe(gid[i]) ][ gid[i] ] = rhs[i];
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        Group::thisProxy[ tope ].addrhs( fromch, p.second );
      }
      if (rhscomplete()) { rhs_complete(); rhs_complete(); }
    }
    //! Receive+add right-hand side vector nonzeros from fellow group branches
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] rhs Portion of the right-hand side vector contributed,
    //!   containing global row indices and values
    void addrhs( int fromch,
                 const std::map< std::size_t, std::vector< tk::real > >& rhs ) {
      for (const auto& r : rhs) {
        m_rhsimport[ fromch ].push_back( r.first );
        m_rhs[ r.first ] += r.second;
      }
      if (rhscomplete()) { rhs_complete(); rhs_complete(); }
    }

    //! Chares contribute their mass diffusion rhs to low order system
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] gid Global row indices of the vector contributed
    //! \param[in] diff Portion of the mass diffusion rhs vector contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charediff( int fromch,
                    const std::vector< std::size_t >& gid,
                    const Fields& diff )
    {
      Assert( gid.size() == diff.nunk(),
              "Size of mass diffusion rhs and row ID vectors must equal" );
      // Store+add vector of nonzero values owned and pack those to be exported
      std::map< int, std::map< std::size_t, std::vector< tk::real > > > exp;
      for (std::size_t i=0; i<gid.size(); ++i)
        if (gid[i] >= m_lower && gid[i] < m_upper) {  // if own
          m_diffimport[ fromch ].push_back( gid[i] );
          m_diff[ gid[i] ] += diff[i];
        } else
          exp[ pe(gid[i]) ][ gid[i] ] = diff[i];
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        Group::thisProxy[ tope ].adddiff( fromch, p.second );
      }
      if (diffcomplete()) diff_complete();
    }
    //! Receive+add massdiffusion rhs vector nonzeros from fellow group branches
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] diff Portion of the right-hand side vector contributed,
    //!   containing global row indices and values
    void adddiff( int fromch, const std::map< std::size_t,
                                std::vector< tk::real > >& diff )
    {
      for (const auto& r : diff) {
        m_diffimport[ fromch ].push_back( r.first );
        m_diff[ r.first ] += r.second;
      }
      if (diffcomplete()) diff_complete();
    }

    //! Chares contribute their lumped mass lhs to low order system
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] gid Global row indices of the vector contributed
    //! \param[in] mass Portion of the lumped mass lhs vector contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charelump( int fromch,
                    const std::vector< std::size_t >& gid,
                    const Fields& mass )
    {
      Assert( gid.size() == mass.nunk(),
              "Size of mass diffusion rhs and row ID vectors must equal" );
      // Store+add vector of nonzero values owned and pack those to be exported
      std::map< int, std::map< std::size_t, std::vector< tk::real > > > exp;
      for (std::size_t i=0; i<gid.size(); ++i)
        if (gid[i] >= m_lower && gid[i] < m_upper) {  // if own
          m_lumpimport[ fromch ].push_back( gid[i] );
          m_lump[ gid[i] ] += mass[i];
        } else
          exp[ pe(gid[i]) ][ gid[i] ] = mass[i];
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        Group::thisProxy[ tope ].addlump( fromch, p.second );
      }
      if (lumpcomplete()) lump_complete();
    }
    //! Receive+add lumped mass lhs vector from fellow group branches
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] mass Portion of the lumped mass lhs vector contributed,
    //!   containing global row indices and values
    void addlump( int fromch, const std::map< std::size_t,
                                std::vector< tk::real > >& mass )
    {
      for (const auto& r : mass) {
        m_lumpimport[ fromch ].push_back( r.first );
        m_lump[ r.first ] += r.second;
      }
      if (lumpcomplete()) lump_complete();
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
      // verify if the BCs to be set match the user's BCs
      #ifndef NDEBUG
      verifybc();
      #endif
      // start collecting boundary condition values
      querybcval();
      // now that the global row ids are complete, build Hypre data from it
      hyprerow();
    }

    //! Chares query side set info
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    const std::map< int, std::vector< std::size_t > >& side() { return m_side; }

    //! Chares offer their global row ids at which they can set BCs
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] bc List of old and new global mesh point (row) indices mapped
    //!   to sides
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charebc( int fromch, const std::vector< std::size_t >& bc ) {
      // Store nodes owned
      if (!bc.empty()) {        // only of chare has anything to offer
        auto& b = m_bc[ fromch ];
        b.insert( end(b), begin(bc), end(bc) );
      }
      // Forward all bcs received to fellow branches
      if (++m_nchbc == m_nchare) {
        auto stream = tk::serialize( m_bc );
        CkCallback cb( CkIndex_LinSysMerger< HostProxy, WorkerProxy >::
                         addbc(nullptr), Group::thisProxy );
        Group::contribute( stream.first, stream.second.get(), BCMapMerger, cb );
      }
    }
    // Reduction target collecting the final aggregated BC node list map
    void addbc( CkReductionMsg* msg ) {
      PUP::fromMem creator( msg->getData() );
      creator | m_bc;
      delete msg;
      bc_complete();
    }

    //! Chares return old global node IDs
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] oldids Vector of old (as in file) global node ID
    //! \details This is the second step to verify if the BCs we will set
    //!   matches the user's BCs. Here we receive the old node IDs. Once all of
    //    them are received, we flatten them and send them to the host, which
    //    does the verification.
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void oldID( int fromch, const std::vector< std::size_t >& oldids ) {
       m_oldbc[ fromch ] = oldids;
       // if we have heard from every chare, flatten, unique, and send to host
       if (m_oldbc.size() == m_bc.size()) {
         std::vector< std::size_t > b;
         for (const auto& c : m_oldbc)
           b.insert( end(b), c.second.cbegin(), c.second.cend() );
         signal2host_verifybc( m_host, b );
       }
    }

    //! Chares contribute their BC values
    //! \param[in] bcv Vector of pairs of bool and BC value associated to global
    //!   node IDs at which the boundary condition is set. Here the bool
    //!   indicates whether the BC value is set at the given node by the user.
    //!   The size of the vectors is the number of PDEs integrated times the
    //!   number of scalar components in all PDEs.
    void charebcval( const std::unordered_map< std::size_t,
                       std::vector< std::pair< bool, tk::real > > >& bcv )
    {
      for (auto& n : bcv) {
        Assert( n.second.size() == m_ncomp, "The total number of scalar "
          "components does not equal that of set in the BC data structure." );
        m_bcval[ n.first ] = n.second;
      }
      if (--m_nchbcval == 0) {
        bcval_complete(); bcval_complete(); bcval_complete();
      }
    }

    //! \brief Receive BC node/row IDs and values from all other PEs and zero
    //!   the distributed matrix columns at which BCs are set
    void bcval( CkReductionMsg* msg ) {
      std::unordered_map< std::size_t,
                          std::vector< std::pair< bool, tk::real > > > bcval;
      PUP::fromMem creator( msg->getData() );
      creator | bcval;
      delete msg;
      for (const auto& n : bcval) {
        auto& b = n.second;
        for (std::size_t i=0; i<m_ncomp; ++i)
          if (b[i].first)       // zero our portion of the column
            for (auto& r : m_lhs) {
              auto it = r.second.find( n.first );
              if (it != end(r.second) && r.first != n.first)
                it->second[ i ] = 0.0;
            }
      }
      lhsbc_complete();
    }

    //! Reduction target indicating that verification of the BCs are complete
    void vercomplete() { ver_complete(); ver_complete(); }

    //! Compute diagnostics (residuals) and contribute them back to host
    //! \details Diagnostics: L1 norm for all components
    void diagnostics() {
      std::vector< tk::real > diag( m_ncomp, 0.0 );
      for (std::size_t i=0; i<m_hypreSol.size()/m_ncomp; ++i)
        for (std::size_t c=0; c<m_ncomp; ++c)
          diag[c] += std::abs( m_hypreSol[i*m_ncomp+c] );
      // if we have heard from every chare, signal back to host
      if (++m_nchdiag == m_nchare) {
        signal2host_diag( m_host, diag );
        m_nchdiag = 0;
      }
    }

  private:
    HostProxy m_host;           //!< Host proxy
    WorkerProxy m_worker;       //!< Worker proxy
    //! Global (as in file) row id lists mapped to all side set ids of the mesh
    std::map< int, std::vector< std::size_t > > m_side;
    std::size_t m_ncomp;        //!< Number of scalar components per unknown
    std::size_t m_nchare;       //!< Number of chares contributing to my PE
    std::size_t m_nperow;       //!< Number of fellow PEs to send row ids to
    std::size_t m_nchbc;        //!< Number of chares we received bcs from
    std::size_t m_nchbcval;     //!< Number of chares we received bc values from
    std::size_t m_nchdiag;      //!< Number of chares we received diags from
    std::size_t m_lower;        //!< Lower index of the global rows on my PE
    std::size_t m_upper;        //!< Upper index of the global rows on my PE
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
    //!   id during the communication of the mass diffusion rhs vector
    std::map< int, std::vector< std::size_t > > m_diffimport;
    //! \brief Import map associating a list of global row ids to a worker chare
    //!   id during the communication of the lumped mass lhs vector
    std::map< int, std::vector< std::size_t > > m_lumpimport;
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
    //! \brief Part of mass diffusion right-hand side vector owned by my PE
    //! \details Vector of values (for each scalar equation solved) associated
    //!   to global mesh point row ids. This vector collects the mass diffusion
    //!   terms to be added to the right-hand side for the low order solution.
    std::map< std::size_t, std::vector< tk::real > > m_diff;
    //! \brief Part of lumped mass left-hand side vector owned by my PE
    //! \details Vector of values (for each scalar equation solved) associated
    //!   to global mesh point row ids. This vector collects the nonzero values
    //!   of the lumped mass matrix for the low order solution.
    std::map< std::size_t, std::vector< tk::real > > m_lump;
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
    //! \brief Map associating lists of new global row ids to chare ids at which
    //!   we set boundary conditions on
    //! \details 'New' as in producing contiguous-row-id linear system
    //!   contributions, see also Partitioner.h.
    std::unordered_map< int, std::vector< std::size_t > > m_bc;
    //! Flat list of old global row ids at boundary conditions are set
    //! \details 'Old' as in file, see also Partitioner.h.
    std::unordered_map< int, std::vector< std::size_t > > m_oldbc;
    //! \brief Values (for each scalar equation solved) of Dirichlet boundary
    //!   conditions assigned to global node IDs we set
    //! \details The map key is the global mesh node/row ID, the value is a
    //!   vector of pairs in which the bool (.first) indicates whether the
    //!   boundary condition value (.second) is set at the given node. The size
    //!   of the vectors is the number of PDEs integrated times the number of
    //!   scalar components in all PDEs.
    std::unordered_map< std::size_t,
                        std::vector< std::pair< bool, tk::real > > > m_bcval;

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
    //! Check if our portion of the right-hand side vector values is complete
    //! \return True if all parts of the right-hand side vector have been
    //!   received
    bool rhscomplete() const { return m_rhsimport == m_rowimport; }
    //! Check if our portion of the mass diffusion rhs vector values is complete
    //! \return True if all parts of the mass diffusion rhs vector have been
    //!   received
    bool diffcomplete() const { return m_diffimport == m_rowimport; }
    //! Check if our portion of the lumped mass lhs vector values is complete
    //! \return True if all parts of the lumped mass lhs vector have been
    //!   received
    bool lumpcomplete() const { return m_lumpimport == m_rowimport; }

    //! Verify if the BCs we will set match the user's BCs
    //! \details This is only the first step: we send out queries to workers
    //!   which then return with the old node IDs.
    void verifybc() {
      // Make BC node lists unqiue
      for (auto& s : m_bc) tk::unique( s.second );
      if (m_bc.empty()) { // if no BCs, we are done (nothing to verify)
        ver_complete(); ver_complete();
      } else // query old (as in file) node IDs from chares that offered them
        for (const auto& c : m_bc)
          m_worker[ c.first ].oldID( CkMyPe(), c.second );
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

    //! Collect boundary condition values from workers
    void querybcval() {
      // Create flat node lists at which we set BCs associated to chares
      std::unordered_map< int, std::vector< std::size_t > > b;
      for (const auto& c : m_bc) {
        auto& q = b[ c.first ];
        for (auto n : c.second)
          if (n >= m_lower && n < m_upper)  // if own
            q.push_back( n );
      }
      // Query BC values from chares
      if (b.empty()) {
        bcval_complete(); bcval_complete(); bcval_complete();
      } else {
        m_nchbcval = b.size();
        for (const auto& c : b) m_worker[ c.first ].bcval( CkMyPe(), c.second );
      }
    }

    //! Set boundary conditions on the left-hand side matrix
    //! \details Here we only zero the row and put 1.0 in the diagonal, zeroing
    //!   the column requires communication among all fellow PEs and is done at
    //!   the receiving side in bcval().
    void lhsbc() {
      Assert( lhscomplete(),
              "Nonzero values of distributed matrix on PE " +
              std::to_string( CkMyPe() ) + " is incomplete: cannot set BCs" );
      for (const auto& n : m_bcval) {
        auto& row = tk::ref_find( m_lhs, n.first );
        auto& b = n.second;
        for (std::size_t i=0; i<m_ncomp; ++i)
          if (b[i].first) { // zero row and put 1.0 in diagonal
            for (auto& col : row) col.second[i] = 0.0;
            tk::ref_find( row, n.first )[ i ] = 1.0;
          }
      }
      // Export BC node IDs and values to all other PEs so that all of the
      // distributed matrix columns can be zeroed on all PEs where BCs are set
      auto stream = tk::serialize( m_bcval );
      CkCallback c( CkIndex_LinSysMerger< HostProxy, WorkerProxy >::
                      bcval(nullptr), Group::thisProxy );
      Group::contribute( stream.first, stream.second.get(), BCValMerger, c );
    }

    //! Set boundary conditions on the right-hand side vector
    void rhsbc() {
      Assert( rhscomplete(),
              "Values of distributed right-hand-side vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete: cannot set BCs" );
      for (const auto& n : m_bcval) {
        auto& row = tk::ref_find( m_rhs, n.first );
        auto& b = n.second;
        for (std::size_t i=0; i<m_ncomp; ++i)
          if (b[i].first) row[i] = b[i].second;  // put in BC value
      }
      rhsbc_complete();
    }

    //! Set boundary conditions on the lumped mass left-hand side matrix
    //! \details We put 1.0 in the diagonal and we are done.
    void lumpbc() {
      Assert( lumpcomplete(),
              "Values of distributed lumped mass lhs vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete: cannot set BCs" );
      for (const auto& n : m_bcval) {
        auto& row = tk::ref_find( m_lump, n.first );
        auto& b = n.second;
        for (std::size_t i=0; i<m_ncomp; ++i)
          if (b[i].first) row[i] = 1.0; // put 1.0 in diagonal
      }
      lumpbc_complete();
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

    //! Update high order solution vector in our PE's workers
    void updateHighSol() {
      // Get solution vector values for our PE
      m_x.get( static_cast< int >( (m_upper - m_lower)*m_ncomp ),
               m_hypreRows.data(),
               m_hypreSol.data() );
      // Group solution vector by workers and send each the parts back to
      // workers that own them
      for (const auto& w : m_solimport) {
        std::vector< std::size_t > gid;
        std::vector< tk::real > sol;
        for (auto r : w.second) {
          const auto it = m_sol.find( r );
          if (it != end(m_sol)) {
            gid.push_back( it->first );
            auto i = tk::cref_find( m_lid, it->first );
            using diff_type = typename decltype(m_hypreSol)::difference_type;
            auto b = static_cast< diff_type >( i*m_ncomp );
            auto e = static_cast< diff_type >( (i+1)*m_ncomp );
            sol.insert( end(sol),
                        std::next( begin(m_hypreSol), b ),
                        std::next( begin(m_hypreSol), e ) );
          } else
            Throw( "Can't find global row id " + std::to_string(r) +
                   " to export in high order solution vector" );
        }
        m_worker[ w.first ].updateHighSol( gid, sol );
      }
    }

    //! Solve high order linear system
    void solve() {
      m_solver.solve( m_A, m_b, m_x );
      solve_complete();
    }

    //! Test if all keys of two maps are equal
    //! \param[in] a 1st map to compare
    //! \param[in] b 2nd map to compare
    //! \return True if the maps have the same size and all keys (and only the
    //!   keys) of the two maps are equal
    //! \note It is an error to call this function with unequal-size maps,
    //!   triggering an exception in DEBUG mode.
    //! \note Operator != is used to compare the map keys.
    template< typename Key >
    bool keyEqual( const std::map< Key, std::vector< tk::real > >& a,
                   const std::map< Key, std::vector< tk::real > >& b ) const
    {
      Assert( a.size() == b.size(), "Size mismatch comparing maps" );
      auto ia = a.cbegin();
      auto ib = b.cbegin();
      while (ia != a.cend()) {
        if (ia->first != ib->first) return false;
        ++ia;
        ++ib;
      }
      return true;
    }

    //! Update low order solution vector in our PE's workers
    void updateLowSol() {
      // Group solution vector by workers and send each the parts back to
      // workers that own them
      for (const auto& w : m_solimport) {
        std::vector< std::size_t > gid;
        std::vector< tk::real > sol;
        for (auto r : w.second) {
          const auto it = m_diff.find( r );
          if (it != end(m_diff)) {
            gid.push_back( it->first );
            sol.insert( end(sol), begin(it->second), end(it->second) );
          } else
            Throw( "Can't find global row id " + std::to_string(r) +
                   " to export in low order solution vector" );
        }
        m_worker[ w.first ].updateLowSol( gid, sol );
      }
    }

    //! Solve low order linear system
    void lowsolve() {
      Assert( rhscomplete(),
              "Values of distributed right-hand-side vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete: cannot solve low "
              "order system" );
      Assert( diffcomplete(),
              "Values of distributed mass diffusion rhs vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete: cannot solve low "
              "order system" );
      Assert( lumpcomplete(),
              "Values of distributed lumped mass lhs vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete: cannot solve low "
              "order system" );
      Assert( keyEqual( m_rhs, m_diff ), "Row IDs of rhs and mass diffusion "
              "rhs vector unequal on PE " + std::to_string( CkMyPe() ) + ": "
              "cannot solve low order system" );
      Assert( keyEqual( m_rhs, m_lump ), "Row IDs of rhs and lumped mass lhs "
              "vector unequal on PE " + std::to_string( CkMyPe() ) + ": cannot "
              "solve low order system" );
      auto ir = m_rhs.cbegin();
      auto id = m_diff.begin();
      auto im = m_lump.cbegin();
      while (ir != m_rhs.cend()) {
        const auto& r = ir->second;
        const auto& m = im->second;
        auto& d = id->second;
        Assert( r.size()==m_ncomp && m.size()==m_ncomp && d.size()==m_ncomp,
                "Wrong number of components in solving the low order system" );
        for (std::size_t i=0; i<m_ncomp; ++i) d[i] = (r[i]+d[i])/m[i];
        ++ir; ++id; ++im;
      }
      lowsolve_complete();
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
    //! Send nodelists of BCs we set to host to verify aggregated BCs
    void signal2host_verifybc( const inciter::CProxy_Transporter& host,
                               const std::vector< std::size_t >& bc ) {
      using inciter::CkIndex_Transporter;
      auto stream = tk::serialize( bc );
      CkCallback c( CkIndex_Transporter::verifybc(nullptr), host );
      Group::contribute( stream.first, stream.second.get(), BCVectorMerger, c );
    }
    //! \brief Signal back to host that enabling the SDAG waits for assembling
    //!    the right-hand side is complete and ready for a new advance in time
    void signal2host_advance( const inciter::CProxy_Transporter& host ) {
      using inciter::CkIndex_Transporter;
      Group::contribute(
       CkCallback( CkIndex_Transporter::redn_wrapper_advance(NULL), host ) );
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
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
  #pragma GCC diagnostic ignored "-Wreorder"
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Weffc++"
  #pragma GCC diagnostic ignored "-Wunused-parameter"
  #pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#define CK_TEMPLATES_ONLY
#include "linsysmerger.def.h"
#undef CK_TEMPLATES_ONLY

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // LinSysMerger_h
