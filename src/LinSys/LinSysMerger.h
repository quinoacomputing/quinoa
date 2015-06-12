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
#include "HypreMatrix.h"
#include "HypreVector.h"

namespace tk {

//! Linear system merger Charm++ chare group class
//! \details Instantiations of LinSysMerger comprise a processor aware Charm++
//!   chare group. When instantiated, a new object is created on each PE and not
//!   more (as opposed to individual chares or chare array object elements). The
//!   group's elements are used to collect information from all chare objects
//!   that happen to be on a given PE. See also the Charm++ interface file
//!   linsysmerger.ci. The class is templated so that the same code
//!   (parameterized by the template arguments) can be generated for interacting
//!   with different types of Charm++ proxies.
//! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html
//! \author J. Bakosi
template< class HostProxy  >
class LinSysMerger : public CBase_LinSysMerger< HostProxy > {

  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  LinSysMerger_SDAG_CODE

  private:
   using Group = CBase_LinSysMerger< HostProxy >;

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
      wait4lhs();
      wait4rhs();
      wait4hypremat();
      wait4hyprevec();
      wait4filllhs();
      wait4fillrhs();
      wait4asm();
      wait4stat();
    }

    void create() {
      tk::Timer t;
      // Create my PE's lhs matrix distributed across all PEs
      m_A.create( m_lower, m_upper );
      // Create my PE's rhs and unknown vectors distributed across all PEs
      m_b.create( m_lower, m_upper );
      m_x.create( m_lower, m_upper );
      m_timestamp.emplace_back( "Create distributed linear system", t.dsec() );
    }

    //! Chares contribute their matrix nonzero values
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
          auto pe = gid / m_chunksize;
          if (pe == CkNumPes()) --pe;
          exp[pe][gid] = r.second;
          m_compts += rowsize;
        }
      }

//       std::cout << CkMyPe() << ": nrows=" << m_lhs.size() << ", ";
//       for (const auto& r : m_lhs) {
//         std::cout << "(" << r.first << ":" << r.second.size() << ") ";
//         for (const auto& c : r.second) std::cout << c.first << " ";
//       }
//       for (const auto& p : exp) {
//         std::cout << "e:" << p.first << " ";
//         for (const auto& r : p.second) {
//           std::cout << "(" << r.first << ":" << r.second.size() << ") ";
//           for (const auto& c : r.second) std::cout << c.first << " ";
//         }
//       }
//       std::cout << std::endl;

      // Export non-owned matrix rows values to fellow branches that own them
      for (const auto& p : exp)
        Group::thisProxy[ static_cast<int>(p.first) ].addlhs( p.second );
      // If our portion is complete, we are done
      if (lhscomplete()) trigger_lhs_complete();
    }

    //! Receive matrix nonzeros from fellow group branches
    void addlhs( const std::map< std::size_t,
                                 std::map< std::size_t, tk::real > >& lhs ) {
//       std::cout << "import on " << CkMyPe() << ": nrows=" << lhs.size() << ", ";
//       for (const auto& r : lhs) {
//         std::cout << "(" << r.first << ":" << r.second.size() << ") ";
//         for (const auto& c : r.second) std::cout << c.first << " ";
//       }
//       std::cout << std::endl;
      for (const auto& r : lhs) m_lhs[ r.first ] = r.second;
      if (lhscomplete()) trigger_lhs_complete();
    }

    //! Chares contribute their rhs nonzero values
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charerhs( const std::map< std::size_t, tk::real >& rhs ) {
      m_timer[ TimerTag::RHS ]; // start measuring merging of rhs
      // Store vector nonzero values owned and pack those to be exported
      std::map< std::size_t, std::map< std::size_t, tk::real > > exp;
      for (const auto& r : rhs) {
        auto gid = r.first;
        if (gid >= m_lower && gid < m_upper) {  // if own
          m_rhs[gid] = r.second;
        } else {
          auto pe = gid / m_chunksize;
          if (pe == CkNumPes()) --pe;
          exp[pe][gid] = r.second;
        }
      }
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp)
        Group::thisProxy[ static_cast<int>(p.first) ].addrhs( p.second );
      // If our portion is complete, we are done
      if (rhscomplete()) trigger_rhs_complete();
    }

    //! Receive vector nonzeros from fellow group branches
    void addrhs( const std::map< std::size_t, tk::real >& rhs ) {
      for (const auto& r : rhs) m_rhs[ r.first ] = r.second;
      if (rhscomplete()) trigger_rhs_complete();
    }

  private:
    HostProxy m_host;           //!< Host proxy
    std::size_t m_chunksize;    //!< Number of rows the first npe-1 PE own
    std::size_t m_lower;        //!< Lower index of the global rows for my PE
    std::size_t m_upper;        //!< Upper index of the global rows for my PE
    tk::hypre::HypreMatrix m_A; //!< Hypre matrix to store the lhs
    tk::hypre::HypreVector m_b; //!< Hypre vector to store the rhs
    tk::hypre::HypreVector m_x; //!< Hypre vector to store the unknowns
    //! Sparse matrix: global mesh point row and column ids, and nonzero value
    std::map< std::size_t, std::map< std::size_t, tk::real > > m_lhs;
    //! Right-hand side vector: global mesh point row ids and values
    std::map< std::size_t, tk::real > m_rhs;
    std::vector< int > m_rows;  //!< Row indices for my PE
    std::vector< int > m_ncols; //!< Number of matrix columns/rows for my PE
    std::vector< int > m_cols;  //!< Matrix column indices for rows for my PE
    std::vector< tk::real > m_mat;  //!< Matrix nonzero values for my PE
    std::vector< tk::real > m_vec;  //!< Vector nonzero values for my PE
    std::size_t m_ownpts;       //!< Size (in bytes) of owned matrix nonzeros
    std::size_t m_compts;       //!< size (in bytes) of communicated nonzeros
    //! Time stamps
    std::vector< std::pair< std::string, tk::real > > m_timestamp;
    enum class TimerTag { LHS, RHS };           //!< Timer labels
    std::map< TimerTag, tk::Timer > m_timer;    //!< Timers
    //! Performance statistics
    std::vector< std::pair< std::string, tk::real > > m_perfstat;

    //! Check if our portion of the matrix values is complete
    bool lhscomplete() const {
      return m_lhs.size() == m_upper-m_lower &&
             m_lhs.cbegin()->first == m_lower &&
             (--m_lhs.cend())->first == m_upper-1;
    }

    //! Check if our portion of the vector values is complete
    bool rhscomplete() const {
      return m_rhs.size() == m_upper-m_lower &&
             m_rhs.cbegin()->first == m_lower &&
             (--m_rhs.cend())->first == m_upper-1;
    }

    //! Build Hypre data for our portion of the matrix
    void hyprelhs() {
      tk::Timer t;
      Assert( lhscomplete(),
              "Nonzero values of distributed matrix on PE " +
              std::to_string( CkMyPe() ) + " is incomplete" );
      for (const auto& r : m_lhs) {
        m_rows.push_back( static_cast< int >( r.first ) );
        m_ncols.push_back( static_cast< int >( r.second.size() ) );
        for (const auto& c : r.second) {
           m_cols.push_back( static_cast< int >( c.first ) );
           m_mat.push_back( c.second );
         }
      }

//       std::cout << CkMyPe() << ": ";
//       for (const auto& r : m_lhs) {
//         std::cout << "(" << r.first << ") ";
//         for (const auto& c : r.second)
//           std::cout << c.first << " ";
//       }
//       std::cout << '\n';

      m_timestamp.emplace_back( "Build Hypre data for lhs matrix", t.dsec() );
      trigger_hyprelhs_complete();
    }

    //! Build Hypre data for our portion of the vector
    void hyprerhs() {
      tk::Timer t;
      Assert( rhscomplete(),
              "Values of distributed vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete" );
      for (const auto& r : m_rhs) m_vec.push_back( r.second );
      m_timestamp.emplace_back( "Build Hypre data for rhs vector", t.dsec() );
      trigger_hyprerhs_complete();
    }

    //! Set our portion of values of the distributed matrix
    void lhs() {
      tk::Timer t;
      Assert( m_mat.size() == m_cols.size(),
              "Matrix values incomplete on " + std::to_string(CkMyPe()) );
      // Set our portion of the matrix values
      m_A.set( static_cast< int >( m_upper - m_lower ),
               m_ncols.data(),
               m_rows.data(),
               m_cols.data(),
               m_mat.data() );
      m_timestamp.emplace_back( "Fill lhs matrix", t.dsec() );
      trigger_filllhs_complete();
    }

    //! Set our portion of values of the distributed vector
    void rhs() {
      tk::Timer t;
      Assert( m_vec.size() == m_rows.size(),
              "Vector values incomplete on " + std::to_string(CkMyPe()) );
      // Set our portion of the vector values
      m_b.set( static_cast< int >( m_upper - m_lower ),
               m_rows.data(),
               m_vec.data() );
      m_timestamp.emplace_back( "Fill rhs vector", t.dsec() );
      trigger_fillrhs_complete();
    }

    //! Assemble distributed matrix
    void assemblelhs() {
      tk::Timer t;
      m_A.assemble();
      m_timestamp.emplace_back( "Assemble lhs matrix", t.dsec() );
      trigger_asmlhs_complete();
    }

    //! Assemble distributed vector
    void assemblerhs() {
      tk::Timer t;
      m_b.assemble();
      m_timestamp.emplace_back( "Assemble rhs vector", t.dsec() );
      trigger_asmrhs_complete();
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
