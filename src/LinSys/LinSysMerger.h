//******************************************************************************
/*!
  \file      src/LinSys/LinSysMerger.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 03:23:09 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Linear system merger
  \details   Linear system merger.
*/
//******************************************************************************
#ifndef LinSysMerger_h
#define LinSysMerger_h

#include <vector>
#include <map>
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
    //! \param[in] npoin Total number of mesh points
    LinSysMerger( HostProxy& host, std::size_t npoin ) :
      m_host( host ),
      m_chunksize( npoin / static_cast<std::size_t>(CkNumPes()) ),
      m_lower( static_cast<std::size_t>(CkMyPe()) * m_chunksize ),
      m_upper( m_lower + m_chunksize ),
      m_ownpts( 0 ),
      m_compts( 0 )
    {
      auto remainder = npoin % static_cast<std::size_t>(CkNumPes());
      if (remainder && CkMyPe() == CkNumPes()-1) m_upper += remainder;
      //std::cout << CkMyPe() << ": [" << m_lower << "..." << m_upper << ")\n";
      // Create distributed linear system
      create();
      // Activate SDAG waits
      wait4rows();
      wait4lhs();
      wait4rhs();
      wait4sol();
      wait4hyprelhs();
      wait4hyprerhs();
      wait4hypresol();
      wait4filllhs();
      wait4fillrhs();
      wait4fillsol();
      wait4asm();
      wait4stat();
    }

    //! \brief Create linear system, i.e., left-hand side matrix, vector of
    //!   unknowns and right-hand side vector and perform their initialization
    void create() {
      tk::Timer t;
      // Create my PE's lhs matrix distributed across all PEs
      m_A.create( m_lower, m_upper );
      // Create my PE's rhs and unknown vectors distributed across all PEs
      m_b.create( m_lower, m_upper );
      m_x.create( m_lower, m_upper );
      m_timestamp.emplace_back( "Create distributed linear system", t.dsec() );
    }

    //! Chares contribute their global row ids
    //! \param[in] worker Worker proxy contribution coming from
    //! \param[in] id Charm chare array index contribution coming from
    //! \param[in] row Global mesh point (row) indices contributed
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charerow( WorkerProxy& worker,
                   int id,
                   const std::vector< std::size_t >& row )
    {
      // Store worker proxy
      m_worker = worker;
      // Collect ids of workers on my PE
      m_myworker.push_back( id );
      // Store rows owned and pack those to be exported
      std::map< std::size_t, std::set< std::size_t > > exp;
      for (auto gid : row) {
        if (gid >= m_lower && gid < m_upper)    // if own
          m_rows.push_back( static_cast<int>(gid) );
        else
          exp[ pe(gid) ].insert( gid );
      }
      tk::unique( m_rows );
      // Export non-owned global rows to fellow branches that own them
      for (const auto& p : exp)
        Group::thisProxy[ static_cast<int>(p.first) ].addrow( p.second );
      // If our portion is complete, we are done
      if (rowcomplete()) trigger_row_complete();
    }

    //! Receive global row ids from fellow group branches
    //! \param[in] row Global mesh point (row) indices received
    void addrow( const std::set< std::size_t >& row ) {
      for (const auto& r : row) m_rows.push_back( static_cast<int>(r) );
      tk::unique( m_rows );
      if (rowcomplete()) trigger_row_complete();
    }

    //! Chares contribute their matrix nonzero values
    //! \param[in] lhs Portion of the left-hand side matrix contributed,
    //!   containing global row and column indices and non-zero values
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charelhs( const std::map< std::size_t,
                                   std::map< std::size_t, tk::real > >& lhs )
    {
      m_timer[ TimerTag::LHS ]; // start measuring merging of lhs
      // Store matrix nonzero values owned and pack those to be exported
      std::map< std::size_t,
                std::map< std::size_t,
                          std::map< std::size_t, tk::real > > > exp;
      for (const auto& r : lhs) {
        auto rowsize = r.second.size() * sizeof( decltype(exp)::value_type );
        auto gid = r.first;
        if (gid >= m_lower && gid < m_upper) {  // if own
          m_lhs[gid] = r.second;
          m_ownpts += rowsize;
        } else {
          exp[ pe(gid) ][ gid ] = r.second;
          m_compts += rowsize;
        }
      }
      // Export non-owned matrix rows values to fellow branches that own them
      for (const auto& p : exp)
        Group::thisProxy[ static_cast<int>(p.first) ].addlhs( p.second );
      // If our portion is complete, we are done
      if (lhscomplete()) trigger_lhs_complete();
    }

    //! Receive matrix nonzeros from fellow group branches
    //! \param[in] lhs Portion of the left-hand side matrix contributed,
    //!   containing global row and column indices and non-zero values
    void addlhs( const std::map< std::size_t,
                                 std::map< std::size_t, tk::real > >& lhs ) {
      for (const auto& r : lhs) m_lhs[ r.first ] = r.second;
      if (lhscomplete()) trigger_lhs_complete();
    }

    //! Chares contribute their rhs nonzero values
    //! \param[in] rhs Portion of the right-hand side vector contributed,
    //!   containing global row indices and values
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charerhs( const std::map< std::size_t, tk::real >& rhs ) {
      m_timer[ TimerTag::RHS ]; // start measuring merging of rhs
      // Store vector nonzero values owned and pack those to be exported
      std::map< std::size_t, std::map< std::size_t, tk::real > > exp;
      for (const auto& r : rhs) {
        auto gid = r.first;
        if (gid >= m_lower && gid < m_upper)    // if own
          m_rhs[gid] = r.second;
        else
          exp[ pe(gid) ][ gid ] = r.second;
      }
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp)
        Group::thisProxy[ static_cast<int>(p.first) ].addrhs( p.second );
      // If our portion is complete, we are done
      if (rhscomplete()) trigger_rhs_complete();
    }

    //! Receive RHS vector nonzeros from fellow group branches
    //! \param[in] rhs Portion of the right-hand side vector contributed,
    //!   containing global row indices and values
    void addrhs( const std::map< std::size_t, tk::real >& rhs ) {
      for (const auto& r : rhs) m_rhs[ r.first ] = r.second;
      if (rhscomplete()) trigger_rhs_complete();
    }

    //! Chares contribute their solution nonzero values
    //! \param[in] id Charm chare array index contribution coming from
    //! \param[in] sol Portion of the unknown/solution vector contributed,
    //!   containing global row indices and values
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charesol( int id, const std::map< std::size_t, tk::real >& sol ) {
      m_timer[ TimerTag::SOL ]; // start measuring merging of rhs
      // Store vector nonzero values owned and pack those to be exported
      std::map< std::size_t, std::map< std::size_t, tk::real > > exp;
      for (const auto& r : sol) {
        auto gid = r.first;
        if (gid >= m_lower && gid < m_upper) {  // if own
          m_comm[id].push_back( gid );
          m_sol[gid] = r.second;
        } else {
          exp[ pe(gid) ][ gid ] = r.second;
        }
      }
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp)
        Group::thisProxy[ static_cast<int>(p.first) ].addsol( id, p.second );
      // If our portion is complete, we are done
      if (solcomplete()) trigger_sol_complete();
    }

    //! Receive solution vector nonzeros from fellow group branches
    //! \param[in] id Charm chare array index contribution coming from
    //! \param[in] sol Portion of the unknown/solution vector contributed,
    //!   containing global row indices and values
    void addsol( int id, const std::map< std::size_t, tk::real >& sol ) {
      for (const auto& r : sol) {
        m_comm[id].push_back( r.first );
        m_sol[ r.first ] = r.second;
      }
      if (solcomplete()) trigger_sol_complete();
    }

  private:
    HostProxy m_host;           //!< Host proxy
    WorkerProxy m_worker;       //!< Worker proxy
    std::size_t m_chunksize;    //!< Number of rows the first npe-1 PE own
    std::size_t m_lower;        //!< Lower index of the global rows for my PE
    std::size_t m_upper;        //!< Upper index of the global rows for my PE
    std::vector< int > m_myworker; //!< Ids of workers on my PE
    //! \brief Export map associating the list of global row (mesh point) ids
    //!   owned to the chare array element index a contribution came from
    std::map< int, std::vector< std::size_t > > m_comm;
    tk::hypre::HypreMatrix m_A; //!< Hypre matrix to store the lhs
    tk::hypre::HypreVector m_b; //!< Hypre vector to store the rhs
    tk::hypre::HypreVector m_x; //!< Hypre vector to store the unknowns
    //! Sparse matrix: global mesh point row and column ids, and nonzero value
    std::map< std::size_t, std::map< std::size_t, tk::real > > m_lhs;
    //! Right-hand side vector: global mesh point row ids and values
    std::map< std::size_t, tk::real > m_rhs;
    //! Unknown/solution vector: global mesh point row ids and values
    std::map< std::size_t, tk::real > m_sol;
    std::vector< int > m_rows;  //!< Row indices for my PE
    std::vector< int > m_ncols; //!< Number of matrix columns/rows for my PE
    std::vector< int > m_cols;  //!< Matrix column indices for rows for my PE
    std::vector< tk::real > m_hypreMat; //!< Matrix nonzero values for my PE
    std::vector< tk::real > m_hypreRhs; //!< RHS vector nonzero values for my PE
    std::vector< tk::real > m_hypreSol; //!< Sol vector nonzero values for my PE
    std::size_t m_ownpts;       //!< Size (in bytes) of owned matrix nonzeros
    std::size_t m_compts;       //!< size (in bytes) of communicated nonzeros
    //! Time stamps
    std::vector< std::pair< std::string, tk::real > > m_timestamp;
    enum class TimerTag { LHS, RHS, SOL };      //!< Timer labels
    std::map< TimerTag, tk::Timer > m_timer;    //!< Timers
    //! Performance statistics
    std::vector< std::pair< std::string, tk::real > > m_perfstat;

    //! Return processing element for global mesh row id
    //! \param[in] gid Global mesh point (matrix or vector row) id
    //! \return PE that owns global row id
    std::size_t pe( std::size_t gid ) {
      auto pe = gid / m_chunksize;
      if (pe == CkNumPes()) --pe;
      return pe;
    }

    //! Check if our portion of the global row ids is complete
    //! \return True if all owned rows have been received
    bool rowcomplete() const { return m_rows.size() == m_upper-m_lower; }

    //! Check if our portion of the matrix values is complete
    //! \return True if all parts of the left-hand side matrix have been
    //!   received
    bool lhscomplete() const {
      return m_lhs.size() == m_upper-m_lower &&
             m_lhs.cbegin()->first == m_lower &&
             (--m_lhs.cend())->first == m_upper-1;
    }

    //! Check if our portion of the right-hand side vector values is complete
    //! \return True if all parts of the right-hand side vector have been
    //!   received
    bool rhscomplete() const {
      return m_rhs.size() == m_upper-m_lower &&
             m_rhs.cbegin()->first == m_lower &&
             (--m_rhs.cend())->first == m_upper-1;
    }

    //! Check if our portion of the solution vector values is complete
    //! \return True if all parts of the unknown/solution vector have been
    //!   received
    bool solcomplete() const {
      return m_sol.size() == m_upper-m_lower &&
             m_sol.cbegin()->first == m_lower &&
             (--m_sol.cend())->first == m_upper-1;
    }

    //! Instruct performers on our PE to start building the linear system
    void buildSystem() {
      for (auto w : m_myworker) m_worker[ w ].buildSystem();
    }

    //! Update solution vector in our PE's performers
    void updateSolution() {
      for (const auto& w : m_comm) {
        std::map< std::size_t, tk::real > sol;
        for (auto r : w.second) {
          const auto it = m_sol.find( r );
          if (it != end(m_sol))
            sol.emplace( it->first, -it->second );
          else
            Throw( "Can't find global row id " + std::to_string(r) +
                   " to export in solution vector" );
        }
        m_worker[ w.first ].updateSolution( sol );
      }
    }

    //! Build Hypre data for our portion of the matrix
    void hyprelhs() {
      tk::Timer t;
      Assert( lhscomplete(),
              "Nonzero values of distributed matrix on PE " +
              std::to_string( CkMyPe() ) + " is incomplete" );
      for (const auto& r : m_lhs) {
        //m_rows.push_back( static_cast< int >( r.first ) );
        m_ncols.push_back( static_cast< int >( r.second.size() ) );
        for (const auto& c : r.second) {
           m_cols.push_back( static_cast< int >( c.first ) );
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

    //! Build Hypre data for our portion of the solution vector
    void hypresol() {
      tk::Timer t;
      Assert( solcomplete(),
              "Values of distributed solution vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete" );
      for (const auto& r : m_sol) m_hypreSol.push_back( r.second );
      m_timestamp.emplace_back( "Build Hypre data for solution vector", t.dsec() );
      trigger_hypresol_complete();
    }

    //! Set our portion of values of the distributed matrix
    void lhs() {
      tk::Timer t;
      Assert( m_hypreMat.size() == m_cols.size(),
              "Matrix values incomplete on " + std::to_string(CkMyPe()) );
      // Set our portion of the matrix values
      m_A.set( static_cast< int >( m_upper - m_lower ),
               m_ncols.data(),
               m_rows.data(),
               m_cols.data(),
               m_hypreMat.data() );
      m_timestamp.emplace_back( "Fill left-hand side matrix", t.dsec() );
      trigger_filllhs_complete();
    }

    //! Set our portion of values of the distributed right-hand side vector
    void rhs() {
      tk::Timer t;
      Assert( m_hypreRhs.size() == m_rows.size(),
              "RHS vector values incomplete on " + std::to_string(CkMyPe()) );
      // Set our portion of the vector values
      m_b.set( static_cast< int >( m_upper - m_lower ),
               m_rows.data(),
               m_hypreSol.data() );
      m_timestamp.emplace_back( "Fill right-hand side vector", t.dsec() );
      trigger_fillrhs_complete();
    }

    //! Set our portion of values of the distributed solution vector
    void sol() {
      tk::Timer t;
      Assert( m_hypreSol.size() == m_rows.size(),
              "Solution vector values incomplete on " + std::to_string(CkMyPe()) );
      // Set our portion of the vector values
      m_x.set( static_cast< int >( m_upper - m_lower ),
               m_rows.data(),
               m_hypreRhs.data() );
      m_timestamp.emplace_back( "Fill solution vector", t.dsec() );
      trigger_fillsol_complete();
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

    //! Assemble distributed solution vector
    void assemblesol() {
      tk::Timer t;
      m_x.assemble();
      m_timestamp.emplace_back( "Assemble solution vector", t.dsec() );
      trigger_asmsol_complete();
    }

    //! Signal back to host that the initialization of the matrix is complete
    //! \details This function contributes to a reduction on all branches (PEs)
    //!   of LinSysMerger to the host, inciter::CProxy_Conductor, given by a
    //!   template argument. This is an overload on the specialization,
    //!   inciter::CProxy_Conductor, of the LinSysMerger template. It creates a
    //!   Charm++ reduction target via creating a callback that invokes the
    //!   typed reduction client, where host is the proxy on which the
    //!   reduction target method, init(), is called upon completion of the
    //!   reduction. Note that we do not use Charm++'s CkReductionTarget macro,
    //!   but explicitly generate the code that the macro would generate. To
    //!   explain why here is Charm++'s CkReductionTarget macro's definition,
    //!   defined in ckreduction.h:
    //!   \code{.cpp}
    //!      #define CkReductionTarget(me, method) \
    //!        CkIndex_##me::redn_wrapper_##method(NULL)
    //!   \endcode
    //!   which takes arguments 'me' (a class name) and 'method' a member
    //!   function of class 'me' and generates the call
    //!   'CkIndex_<class>::redn_wrapper_<method>(NULL)'. With this overload to
    //!   contributeTo() we do the above macro's job for LinSysMerger
    //!   specialized on class inciter::CProxy_Conductor and its init()
    //!   reduction target. This is required since (1) Charm++'s
    //!   CkReductionTarget macro's preprocessing happens earlier than type
    //!   resolution and the string of the template argument would be
    //!   substituted instead of the type specialized (which not what we want
    //!   here), and (2) the template argument class, CProxy_Conductor, is in a
    //!   namespace different than that of LinSysMerger. When a new class is
    //!   used to specialize LinSysMerger, the compiler will alert that a new
    //!   overload needs to be defined.
    //! \note This simplifies client-code, e.g., Conductor, which now requires
    //!   no explicit book-keeping with counters, etc. Also a reduction (instead
    //!   of a direct call to the host) better utilizes the communication
    //!   network as computational nodes can send their aggregated contribution
    //!   to other nodes on a network instead of all chares sending their
    //!   (smaller) contributions to the same host.
    //! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html,
    //!   Sections "Processor-Aware Chare Collections" and "Chare Arrays".
    void init_complete( const inciter::CProxy_Conductor& host ) {
      using inciter::CkIndex_Conductor;
      Group::contribute(
        CkCallback( CkIndex_Conductor::redn_wrapper_init(NULL), host )
      );
    }
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
