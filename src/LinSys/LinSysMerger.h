//******************************************************************************
/*!
  \file      src/LinSys/LinSysMerger.h
  \author    J. Bakosi
  \date      Mon 28 Dec 2015 12:34:50 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Linear system merger
  \details   Linear system merger.
*/
//******************************************************************************
#ifndef LinSysMerger_h
#define LinSysMerger_h

#include <vector>
#include <map>
#include <unordered_map>
#include <utility>
#include <iosfwd>
#include <cstddef>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "linsysmerger.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include "Types.h"
#include "Timer.h"
#include "Exception.h"
#include "ContainerUtil.h"
#include "HypreMatrix.h"
#include "HypreVector.h"
#include "HypreSolver.h"
#include "Performer.h"

namespace tk {

//! Linear system merger Charm++ chare group class
//! \details Instantiations of LinSysMerger comprise a processor aware Charm++
//!   chare group. When instantiated, a new object is created on each PE and not
//!   more (as opposed to individual chares or chare array object elements). The
//!   group's elements are used to collect information from all chare objects
//!   that happen to be on a given PE. See also the Charm++ interface file
//!   linsysmerger.ci. The class is templated on the host as well as the work
//!   proxy so that the same code (parameterized by the template arguments) can
//!   be generated for interacting with different types of Charm++ proxies.
//! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html
//! \author J. Bakosi
template< class HostProxy, class WorkerProxy  >
class LinSysMerger : public CBase_LinSysMerger< HostProxy, WorkerProxy > {

  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  LinSysMerger_SDAG_CODE

  private:
    using Group = CBase_LinSysMerger< HostProxy, WorkerProxy >;

  public:
    //! Constructor
    //! \param[in] host Charm++ host proxy
    //! \param[in] div Lower and upper global row indices associated to a PE.
    //!   These are the divisions at which the linear system is divided at along
    //!   PE boundaries.
    LinSysMerger(
      HostProxy& host,
      const std::map< int, std::pair< std::size_t, std::size_t > >& div )
    :
      m_host( host ),
      m_lower( tk::cref_find(div,CkMyPe()).first ),
      m_upper( tk::cref_find(div,CkMyPe()).second ),
      m_nchare( 0 ),
      m_nperow( 0 )
    {
      // Invert and store PE-division map
      for (const auto& p : div) m_div[ p.second ] = p.first;
      // Create distributed linear system
      create();
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
      signal2host_wait4rhs_complete( m_host );
    }

    //! Chares register on my PE
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void checkin() { ++m_nchare; }

    //! Chares contribute their global row ids
    //! \param[in] worker Worker proxy contribution coming from
    //! \param[in] chgid Charm chare global array index contribution coming from
    //! \param[in] chlid Charm chare local array index contribution coming from
    //! \param[in] row Global mesh point (row) indices contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charerow( WorkerProxy& worker,
                   int chgid,
                   int chlid,
                   const std::vector< std::size_t >& row )
    {
      // Store worker proxy associated to chare id
      m_worker[ chgid ] = { worker, chlid };
      // Collect global ids of workers on my PE
      m_myworker.push_back( chgid );
      // Store rows owned and pack those to be exported, also build import map
      // used to test for completion
      std::map< int, std::set< std::size_t > > exp;
      for (auto gid : row) {
        if (gid >= m_lower && gid < m_upper) {  // if own
          m_rowimport[ chgid ].push_back( gid );
          m_row.insert( gid );
        } else exp[ pe(gid) ].insert( gid );
      }
      // Export non-owned parts to fellow branches that own them
      m_nperow += exp.size();
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        Group::thisProxy[ tope ].
          addrow( worker, chgid, chlid, CkMyPe(), p.second );
      }
      if (rowcomplete()) signal2host_row_complete( m_host );
    }
    //! Receive global row ids from fellow group branches
    //! \param[in] worker Worker proxy contribution coming from
    //! \param[in] chgid Charm chare global array index contribution coming from
    //! \param[in] chlid Charm chare local array index contribution coming from
    //! \param[in] frompe PE contribution coming from
    //! \param[in] row Global mesh point (row) indices received
    void addrow( WorkerProxy worker,
                 int chgid,
                 int chlid,
                 int frompe,
                 const std::set< std::size_t >& row )
    {
      // Store worker proxy associated to chare id
      m_worker[ chgid ] = { worker, chlid };
      for (auto r : row) {
        m_rowimport[ chgid ].push_back( r );
        m_row.insert( r );
      }
      Group::thisProxy[ frompe ].recrow();
    }
    //! Acknowledge received row ids
    void recrow() {
      --m_nperow;
      if (rowcomplete()) signal2host_row_complete( m_host );
    };

    //! Chares contribute their solution nonzero values
    //! \param[in] chgid Charm chare global array index contribution coming from
    //! \param[in] sol Portion of the unknown/solution vector contributed,
    //!   containing global row indices and values
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charesol( int chgid, const std::map< std::size_t, tk::real >& sol ) {
      m_timer[ TimerTag::SOL ]; // start measuring merging of solution vector
      // Store solution vector nonzero values owned and pack those to be
      // exported, also build import map used to test for completion
      std::map< int, std::map< std::size_t, tk::real > > exp;
      for (const auto& r : sol) {
        auto gid = r.first;
        if (gid >= m_lower && gid < m_upper) {  // if own
          m_solimport[ chgid ].push_back( gid );
          m_sol[gid] = r.second;
        } else {
          exp[ pe(gid) ][ gid ] = r.second;
        }
      }
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        Group::thisProxy[ tope ].addsol( chgid, p.second );
      }
      if (solcomplete()) trigger_sol_complete();
    }
    //! Receive solution vector nonzeros from fellow group branches
    //! \param[in] chgid Chare id contribution coming from
    //! \param[in] sol Portion of the unknown/solution vector contributed,
    //!   containing global row indices and values
    void addsol( int chgid, const std::map< std::size_t, tk::real >& sol ) {
      for (const auto& r : sol) {
        m_solimport[ chgid ].push_back( r.first );
        m_sol[ r.first ] = r.second;
      }
      if (solcomplete()) trigger_sol_complete();
    }

    //! Chares contribute their matrix nonzero values
    //! \param[in] chgid Charm chare global array index contribution coming from
    //! \param[in] lhs Portion of the left-hand side matrix contributed,
    //!   containing global row and column indices and non-zero values
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charelhs( int chgid,
                   const std::map< std::size_t,
                                   std::map< std::size_t, tk::real > >& lhs )
    {
      m_timer[ TimerTag::LHS ]; // start measuring merging of lhs matrix
      // Store matrix nonzero values owned and pack those to be exported, also
      // build import map used to test for completion
      std::map< int, std::map< std::size_t,
                               std::map< std::size_t, tk::real > > > exp;
      for (const auto& r : lhs) {
        auto gid = r.first;
        if (gid >= m_lower && gid < m_upper) {  // if own
          m_lhsimport[ chgid ].push_back( gid );
          auto& row = m_lhs[gid];
          for (const auto& c : r.second) row[ c.first ] += c.second;
        } else {
          exp[ pe(gid) ][ gid ] = r.second;
        }
      }
      // Export non-owned matrix rows values to fellow branches that own them
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        Group::thisProxy[ tope ].addlhs( chgid, p.second );
      }
      if (lhscomplete()) trigger_lhs_complete();
    }
    //! Receive matrix nonzeros from fellow group branches
    //! \param[in] chgid Charm chare global array index contribution coming from
    //! \param[in] lhs Portion of the left-hand side matrix contributed,
    //!   containing global row and column indices and non-zero values
    void addlhs( int chgid,
                 const std::map< std::size_t,
                                 std::map< std::size_t, tk::real > >& lhs ) {
      for (const auto& r : lhs) {
        m_lhsimport[ chgid ].push_back( r.first );
        auto& row = m_lhs[ r.first ];
        for (const auto& c : r.second) row[ c.first ] += c.second;
      }
      if (lhscomplete()) trigger_lhs_complete();
    }

    //! Chares contribute their rhs nonzero values
    //! \param[in] chgid Charm chare global array index contribution coming from
    //! \param[in] rhs Portion of the right-hand side vector contributed,
    //!   containing global row indices and values
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charerhs( int chgid, const std::map< std::size_t, tk::real >& rhs ) {
      m_timer[ TimerTag::RHS ]; // start measuring merging of rhs vector
      // Store vector nonzero values owned and pack those to be exported
      std::map< int, std::map< std::size_t, tk::real > > exp;
      for (const auto& r : rhs) {
        auto gid = r.first;
        if (gid >= m_lower && gid < m_upper) {  // if own
          m_rhsimport[ chgid ].push_back( gid );
          m_rhs[gid] += r.second;
        } else {
          exp[ pe(gid) ][ gid ] = r.second;
        }
      }
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        Group::thisProxy[ tope ].addrhs( chgid, p.second );
      }
      if (rhscomplete()) trigger_rhs_complete();
    }
    //! Receive right-hand side vector nonzeros from fellow group branches
    //! \param[in] chgid Charm chare global array index contribution coming from
    //! \param[in] rhs Portion of the right-hand side vector contributed,
    //!   containing global row indices and values
    void addrhs( int chgid, const std::map< std::size_t, tk::real >& rhs ) {
      for (const auto& r : rhs) {
        m_rhsimport[ chgid ].push_back( r.first );
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
    std::size_t m_lower;        //!< Lower index of the global rows on my PE
    std::size_t m_upper;        //!< Upper index of the global rows on my PE
    std::size_t m_nchare;       //!< Number of chares contributing to my PE
    std::size_t m_nperow;       //!< Number of fellow PEs to send row ids to
    //! \brief Pair of worker proxy and local chare id associated to chare
    //!   global ids we receive contributions from
    std::map< int, std::pair< WorkerProxy, int > > m_worker;
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
    //! \brief Part of unknown/solution vector owned by my PE: global mesh point
    //!   row ids and values
    std::map< std::size_t, tk::real > m_sol;
    //! \brief Part of left-hand side matrix owned by my PE: global mesh point
    //!   row and column ids, and nonzero value
    std::map< std::size_t, std::map< std::size_t, tk::real > > m_lhs;
    //! \brief Part of right-hand side vector owned by my PE: global mesh point
    //!   row ids and values
    std::map< std::size_t, tk::real > m_rhs;
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
    //! Time stamps
    std::vector< std::pair< std::string, tk::real > > m_timestamp;
    //! Timer labels
    enum class TimerTag { LHS, RHS, SOL };
    //! Timers
    std::map< TimerTag, tk::Timer > m_timer;
    //! Performance statistics
    std::vector< std::pair< std::string, tk::real > > m_perfstat;
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

    //! \brief Create linear system, i.e., left-hand side matrix, vector of
    //!   unknowns, right-hand side vector, solver perform their initialization
    void create() {
      tk::Timer t;
      // Create my PE's lhs matrix distributed across all PEs
      m_A.create( m_lower, m_upper );
      // Create my PE's rhs and unknown vectors distributed across all PEs
      m_b.create( m_lower, m_upper );
      m_x.create( m_lower, m_upper );
      // Create linear solver
      m_solver.create();
      m_timestamp.emplace_back( "Create distributed linear system", t.dsec() );
    }

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
      for (auto r : m_row) m_hypreRows.push_back( static_cast< int >( r+1 ) );
    }

    //! Build Hypre data for our portion of the solution vector
    void hypresol() {
      tk::Timer t;
      Assert( solcomplete(),
              "Values of distributed solution vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete" );
      std::size_t i = 0;
      for (const auto& r : m_sol) {
        m_lid[ r.first ] = i++;
        m_hypreSol.push_back( r.second );
      }
      m_timestamp.emplace_back( "Build Hypre data for solution vector",
                                t.dsec() );
      trigger_hypresol_complete();
    }
    //! Build Hypre data for our portion of the matrix
    //! \note Hypre only likes one-based indexing. Zero-based row indexing fails
    //!   to update the vector with HYPRE_IJVectorGetValues().
    void hyprelhs() {
      tk::Timer t;
      Assert( lhscomplete(),
              "Nonzero values of distributed matrix on PE " +
              std::to_string( CkMyPe() ) + " is incomplete" );
      Assert( m_lhs.size() == m_hypreRows.size(),
              "Left-hand side matrix incomplete on " +
              std::to_string(CkMyPe()) );
      for (auto& r : m_lhs) {
        m_hypreNcols.push_back( static_cast< int >( r.second.size() ) );
        for (const auto& c : r.second) {
           m_hypreCols.push_back( static_cast< int >( c.first+1 ) );
           m_hypreMat.push_back( c.second );
        }
      }
      m_timestamp.emplace_back( "Build Hypre data for lhs matrix", t.dsec() );
      trigger_hyprelhs_complete();
    }
    //! Build Hypre data for our portion of the right-hand side vector
    void hyprerhs() {
      tk::Timer t;
      Assert( rhscomplete(),
              "Values of distributed right-hand-side vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete" );
      for (const auto& r : m_rhs) m_hypreRhs.push_back( r.second );
      m_timestamp.emplace_back( "Build Hypre data for rhs vector", t.dsec() );
      trigger_hyprerhs_complete();
    }

    //! Set our portion of values of the distributed solution vector
    void sol() {
      tk::Timer t;
      Assert( m_hypreSol.size() == m_hypreRows.size(),
              "Solution vector values incomplete on " +
              std::to_string(CkMyPe()) );
      // Set our portion of the vector values
      m_x.set( static_cast< int >( m_upper - m_lower ),
               m_hypreRows.data(),
               m_hypreSol.data() );
      m_timestamp.emplace_back( "Fill solution vector", t.dsec() );
      trigger_fillsol_complete();
    }
    //! Set our portion of values of the distributed matrix
    void lhs() {
      tk::Timer t;
      Assert( m_hypreMat.size() == m_hypreCols.size(),
              "Matrix values incomplete on " + std::to_string(CkMyPe()) );
      // Set our portion of the matrix values
      m_A.set( static_cast< int >( m_upper - m_lower ),
               m_hypreNcols.data(),
               m_hypreRows.data(),
               m_hypreCols.data(),
               m_hypreMat.data() );
      m_timestamp.emplace_back( "Fill left-hand side matrix", t.dsec() );
      trigger_filllhs_complete();
    }
    //! Set our portion of values of the distributed right-hand side vector
    void rhs() {
      tk::Timer t;
      Assert( m_hypreRhs.size() == m_hypreRows.size(),
              "RHS vector values incomplete on " + std::to_string(CkMyPe()) );
      // Set our portion of the vector values
      m_b.set( static_cast< int >( m_upper - m_lower  ),
               m_hypreRows.data(),
               m_hypreRhs.data() );
      m_timestamp.emplace_back( "Fill right-hand side vector", t.dsec() );
      trigger_fillrhs_complete();
    }

    //! Assemble distributed solution vector
    void assemblesol() {
      tk::Timer t;
      m_x.assemble();
      m_timestamp.emplace_back( "Assemble solution vector", t.dsec() );
      trigger_asmsol_complete();
    }
    //! Assemble distributed matrix
    void assemblelhs() {
      tk::Timer t;
      m_A.assemble();
      m_timestamp.emplace_back( "Assemble left-hand side matrix", t.dsec() );
      trigger_asmlhs_complete();
    }
    //! Assemble distributed right-hand side vector
    void assemblerhs() {
      tk::Timer t;
      m_b.assemble();
      m_timestamp.emplace_back( "Assemble right-hand side vector", t.dsec() );
      trigger_asmrhs_complete();
    }

    //! Update solution vector in our PE's performers
    void updateSolution() {
      // Get solution vector values for our PE
      m_x.get( static_cast< int >( m_upper - m_lower ),
               m_hypreRows.data(),
               m_hypreSol.data() );

      // Group solution vector by workers and send each the parts back to
      // workers that own them
      for (const auto& w : m_solimport) {
        std::map< std::size_t, tk::real > sol;
        for (auto r : w.second) {
          const auto it = m_sol.find( r );
          if (it != end(m_sol))
            sol.emplace( it->first, m_hypreSol[tk::val_find(m_lid,it->first)] );
          else
            Throw( "Can't find global row id " + std::to_string(r) +
                   " to export in solution vector" );
        }
        // Find pair of local chare id and proxy of worker chare to update
        const auto& ch = tk::ref_find( m_worker, w.first );
        ch.first[ ch.second ].updateSolution( sol );
      }
    }

    //! Solve linear system
    void solve() {
      tk::Timer t;
      //m_A.print( "hypre_mat" );
      //m_b.print( "hypre_b" );
      //m_x.print( "hypre_x" );
      m_solver.solve( m_A, m_b, m_x );
      m_timestamp.emplace_back( "Solve linear system", t.dsec() );
      //m_x.print( "hypre_sol" );
      updateSolution();
    }

    //! Send timers and performance statistics to host for collection
    void sendTimers() {
       m_host.grpTimestamp( m_timestamp );
       m_host.grpPerfstat( m_perfstat );
       m_timestamp.clear();
       m_perfstat.clear();
       m_timer.clear();
    }

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
        CkCallback( CkIndex_Conductor::redn_wrapper_rowcomplete(NULL), host )
      );
    }
    ///@}
    ///@{
    //! \brief Signal back to host that enabling the SDAG waits for assembling
    //!    the right-hand side is complete and ready for a new advance in time
    void signal2host_wait4rhs_complete( const inciter::CProxy_Conductor& host ) {
      using inciter::CkIndex_Conductor;
      Group::contribute(
        CkCallback( CkIndex_Conductor::redn_wrapper_advance(NULL), host )
      );
    }
    ///@}
};

} // tk::

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#define CK_TEMPLATES_ONLY
#include "linsysmerger.def.h"
#undef CK_TEMPLATES_ONLY

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // LinSysMerger_h
