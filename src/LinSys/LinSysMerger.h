// *****************************************************************************
/*!
  \file      src/LinSys/LinSysMerger.h
  \author    J. Bakosi
  \date      Sun 15 May 2016 08:12:13 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Charm++ chare linear system merger group to solve a linear system
  \details   Charm++ chare linear system merger group used to collect and
    assemble the left hand side matrix, the right hand side vector, and the
    solution (unknown) vector from individual worker (Performer) chares. The
    solution is outsourced to hypre, an MPI-only library.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/LinSys/linsysmerger.ci.

    #### Call graph ####
    The following is a directed acyclic graph (DAG) that outlines the
    asynchronous algorithm implemented in this class The detailed discussion of
    the algorithm is given in the Charm++ interface file linsysmerger.ci,
    which also repeats the graph below using ASCII graphics. On the DAG orange
    fills denote global synchronization points, orange frames with white fill
    are partial synchronization points that overlap with other tasks, and dashed
    lines are potential shortcuts that allow jumping over some of the task-graph
    under some circumstances. See the detailed discussion in linsysmrger.ci.
    \dot
    digraph "LinSysMerger SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      ChSol [ label="ChSol"
              tooltip="chares contribute their solution vector nonzeros"
              URL="\ref tk::LinSysMerger::charesol"];
      ChLhs [ label="ChLhs"
              tooltip="chares contribute their left hand side matrix nonzeros"
              URL="\ref tk::LinSysMerger::charelhs"];
      ChRhs [ label="ChRhs"
              tooltip="chares contribute their right hand side vector nonzeros"
              URL="\ref tk::LinSysMerger::charesol"];
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
      Sol [ label="Sol" tooltip="solve linear system" color="#e6851c"
            URL="\ref tk::LinSysMerger::solve"];
      Upd [ label="Upd" tooltip="update solution" color="#e6851c"
            URL="\ref tk::LinSysMerger::updateSolution"];
      ChSol -> HypreSol -> FillSol -> AsmSol -> Sol [ style="solid" ];
      ChLhs -> HypreLhs -> FillLhs -> AsmLhs -> Sol [ style="solid" ];
      ChRhs -> HypreRhs -> FillRhs -> AsmRhs -> Sol [ style="solid" ];
      Sol -> Upd [ style="solid" ];
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
#include <iosfwd>
#include <cstddef>

#include "Types.h"
#include "Exception.h"
#include "ContainerUtil.h"
#include "MeshNodes.h"
#include "HypreMatrix.h"
#include "HypreVector.h"
#include "HypreSolver.h"

#include "NoWarning/conductor.decl.h"

namespace tk {

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
                  std::size_t ncomp ) :
      __dep(),
      m_host( host ),
      m_worker( worker ),
      m_ncomp( ncomp ),
      m_nchare( 0 ),
      m_nperow( 0 ),
      m_lower( 0 ),
      m_upper( 0 ),
      m_myworker(),
      m_rowimport(),
      m_solimport(),
      m_lhsimport(),
      m_rhsimport(),
      m_row(),
      m_sol(),
      m_lhs(),
      m_rhs(),
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
      m_pe()
    {
      // Activate SDAG waits for assembling lhs, rhs, and solution
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
      wait4hyprerhs();
      wait4fillrhs();
      wait4asm();
      m_rhsimport.clear();
      m_rhs.clear();
      m_hypreRhs.clear();
      trigger_asmsol_complete();
      trigger_asmlhs_complete();
      signal2host_advance( m_host );
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
        Group::thisProxy[ tope ].
          addrow( fromch, CkMyPe(), p.second );
      }
      if (rowcomplete()) signal2host_row_complete( m_host );
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
      if (rowcomplete()) signal2host_row_complete( m_host );
    }

    //! Chares contribute their solution nonzero values
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] gid Global row indices of the vector contributed
    //! \param[in] sol Portion of the unknown/solution vector contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charesol( int fromch,
                   const std::vector< std::size_t >& gid,
                   const MeshNodes& sol )
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
      if (solcomplete()) trigger_sol_complete();
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
      if (solcomplete()) trigger_sol_complete();
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
                   const tk::MeshNodes& lhsd,
                   const tk::MeshNodes& lhso )
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
      if (lhscomplete()) trigger_lhs_complete();
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
      if (lhscomplete()) trigger_lhs_complete();
    }

    //! Chares contribute their rhs nonzero values
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] gid Global row indices of the vector contributed
    //! \param[in] rhs Portion of the right-hand side vector contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charerhs( int fromch,
                   const std::vector< std::size_t >& gid,
                   const MeshNodes& rhs )
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
      if (rhscomplete()) trigger_rhs_complete();
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
      if (rhscomplete()) trigger_rhs_complete();
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
              // if any of the above is not satisfied, the row ids are
              // incomplete
              "Row ids are incomplete on PE " + std::to_string( CkMyPe() ) );
      // now that the global row ids are complete, build Hypre data from it
      hyprerow();
    }

  private:
    HostProxy m_host;           //!< Host proxy
    WorkerProxy m_worker;       //!< Worker proxy
    std::size_t m_ncomp;        //!< Number of scalar components per unknown
    std::size_t m_nchare;       //!< Number of chares contributing to my PE
    std::size_t m_nperow;       //!< Number of fellow PEs to send row ids to
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

    //! Build Hypre data for our portion of the global row ids
    //! \note Hypre only likes one-based indexing. Zero-based row indexing fails
    //!   to update the vector with HYPRE_IJVectorGetValues().
    void hyprerow() {
      for (auto r : m_row) {
        decltype(m_hypreRows) h( m_ncomp );
        std::iota( begin(h), end(h), r*m_ncomp+1 );
        m_hypreRows.insert( end(m_hypreRows), begin(h), end(h) );
      }
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
      trigger_hypresol_complete();
    }
    //! Build Hypre data for our portion of the matrix
    //! \note Hypre only likes one-based indexing. Zero-based row indexing fails
    //!   to update the vector with HYPRE_IJVectorGetValues().
    void hyprelhs() {
      Assert( lhscomplete(),
              "Nonzero values of distributed matrix on PE " +
              std::to_string( CkMyPe() ) + " is incomplete" );
      Assert( m_lhs.size() == m_hypreRows.size() / m_ncomp,
              "Left-hand side matrix incomplete on PE " +
              std::to_string(CkMyPe()) );
      for (const auto& r : m_lhs) {
        for (decltype(m_ncomp) i=0; i<m_ncomp; ++i) {
          m_hypreNcols.push_back( static_cast< int >( r.second.size() ) );
          for (const auto& c : r.second) {
            m_hypreCols.push_back( static_cast< int >( c.first*m_ncomp+i+1 ) );
            m_hypreMat.push_back( c.second[i] );
          }
        }
      }
      trigger_hyprelhs_complete();
    }
    //! Build Hypre data for our portion of the right-hand side vector
    void hyprerhs() {
      Assert( rhscomplete(),
              "Values of distributed right-hand-side vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete" );
      for (const auto& r : m_rhs)
        m_hypreRhs.insert( end(m_hypreRhs), begin(r.second), end(r.second) );
      trigger_hyprerhs_complete();
    }

    //! Set our portion of values of the distributed solution vector
    void sol() {
      Assert( m_hypreSol.size() == m_hypreRows.size(), "Solution vector values "
              "incomplete on PE " + std::to_string(CkMyPe()) );
      // Set our portion of the vector values
      m_x.set( static_cast< int >( (m_upper - m_lower)*m_ncomp ),
               m_hypreRows.data(),
               m_hypreSol.data() );
      trigger_fillsol_complete();
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
      trigger_filllhs_complete();
    }
    //! Set our portion of values of the distributed right-hand side vector
    void rhs() {
      Assert( m_hypreRhs.size() == m_hypreRows.size(), "RHS vector values "
              "incomplete on PE " + std::to_string(CkMyPe()) );
      // Set our portion of the vector values
      m_b.set( static_cast< int >( (m_upper - m_lower)*m_ncomp ),
               m_hypreRows.data(),
               m_hypreRhs.data() );
      trigger_fillrhs_complete();
    }

    //! Assemble distributed solution vector
    void assemblesol() {
      m_x.assemble();
      trigger_asmsol_complete();
    }
    //! Assemble distributed matrix
    void assemblelhs() {
      m_A.assemble();
      trigger_asmlhs_complete();
    }
    //! Assemble distributed right-hand side vector
    void assemblerhs() {
      m_b.assemble();
      trigger_asmrhs_complete();
    }

    //! Update solution vector in our PE's performers
    void updateSolution() {
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
                   " to export in solution vector" );
        }
        m_worker[ w.first ].updateSolution( gid, sol );
      }
    }

    //! Solve linear system
    void solve() {
      m_solver.solve( m_A, m_b, m_x );
      updateSolution();
    }

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wdocumentation"
    #endif
    /** @name Host signal calls
      * \brief These functions signal back to the host via a global reduction
      *   originating from each PE branch
      * \details Singal calls contribute to a reduction on all branches (PEs)
      *   of LinSysMerger to the host, e.g., inciter::CProxy_Conductor, given by
      *   the template argument HostProxy. The signal functions are overloads on
      *   the specialization, e.g., inciter::CProxy_Conductor, of the
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
      *    * the template argument class, e.g, CProxy_Conductor, is in a
      *      namespace different than that of LinSysMerger. When a new class is
      *      used to specialize LinSysMerger, the compiler will alert that a new
      *      overload needs to be defined.
      *
      * \note This simplifies client-code, e.g., inciter::Conductor, which now
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
    void signal2host_row_complete( const inciter::CProxy_Conductor& host ) {
      using inciter::CkIndex_Conductor;
      Group::contribute(
        CkCallback( CkIndex_Conductor::redn_wrapper_rowcomplete(NULL), host ) );
    }
    //! \brief Signal back to host that enabling the SDAG waits for assembling
    //!    the right-hand side is complete and ready for a new advance in time
    void signal2host_advance( const inciter::CProxy_Conductor& host ) {
      using inciter::CkIndex_Conductor;
      Group::contribute(
       CkCallback( CkIndex_Conductor::redn_wrapper_advance(NULL), host ) );
    }
    //! \brief Signal back to host that receiving the inverse PE-division map is
    //!  complete and we are ready for Prformers to start their setup.
    void signal2host_setup( const inciter::CProxy_Conductor& host ) {
      using inciter::CkIndex_Conductor;
      Group::contribute(
       CkCallback( CkIndex_Conductor::redn_wrapper_setup(NULL), host ) );
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
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
  #pragma GCC diagnostic ignored "-Wreorder"
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Weffc++"
  #pragma GCC diagnostic ignored "-Wunused-parameter"
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
