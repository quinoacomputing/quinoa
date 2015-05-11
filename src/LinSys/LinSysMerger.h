//******************************************************************************
/*!
  \file      src/LinSys/LinSysMerger.h
  \author    J. Bakosi
  \date      Mon 11 May 2015 03:07:16 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Linear system merger
  \details   Linear system merger.
*/
//******************************************************************************
#ifndef LinSysMerger_h
#define LinSysMerger_h

#include <iostream>     // NOT REALLY NEEDED
#include <numeric>
#include <limits>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <linsysmerger.decl.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <Types.h>
#include <HypreMatrix.h>
#include <Exception.h>
#include <ContainerUtil.h>

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
      m_upper( m_lower + m_chunksize )
    {
      auto remainder = npoin % static_cast<std::size_t>(CkNumPes());
      if (remainder && CkMyPe() == CkNumPes()-1) m_upper += remainder;
      std::cout << CkMyPe() << ": [" << m_lower << "..." << m_upper << ")\n";
      // Create my PE's part of the lhs matrix distributed across all PEs
      m_A.create( m_lower, m_upper );
      // Start SDAG waits
      wait4nz();
      wait4fill();
      wait4asm();
    }

    //! Chares contribute their matrix nonzero structure
    //! \param[in] point Global mesh point ids sent by chares on our PE
    //! \param[in] psup Points surrounding points, see tk::genPsup()
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void charenz( const std::vector< std::size_t >& point,
                  const std::pair< std::vector< std::size_t >,
                                   std::vector< std::size_t > >& psup )
    {
      Assert( point.size() == psup.second.size()-1,
              "Number of owned points must equal in global id vector and "
              "derived data, points surrounding points, sent by chare." );
      // Lambda to add all column indices to a row
      auto storenz = [ &psup ]( std::size_t p, std::vector< std::size_t >& row )
      {
        for (auto i=psup.second[p]+1; i<=psup.second[p+1]; ++i)
          row.push_back( psup.first[i] );
      };
      // Store matrix nonzero locations owned and pack those to be exported
      std::map< std::size_t,
                std::map< std::size_t, std::vector< std::size_t > > > exp;
      for (std::size_t p=0; p<psup.second.size()-1; ++p) {
        auto gid = point[ p ];
        if (gid >= m_lower && gid < m_upper)    // if own
          storenz( p, m_psup[gid] );
        else {
          auto pe = gid / m_chunksize;
          if (pe == CkNumPes()) --pe;
          storenz( p, exp[pe][gid] );
        }
      }
      // If all chares contributed and our portion is complete
      if (complete()) nz_complete();
      // Export non-owned matrix rows to fellow branches that own them
      for (const auto& p : exp)
        Group::thisProxy[ static_cast<int>(p.first) ].add( p.second );
    }

    //! Receive mesh point ids from fellow group branches
    //! \param[in] transfer Mesh point ids received
    void add( const std::map< std::size_t, std::vector< std::size_t > >& psup )
    {
      // Store imported matrix nonzero structure contributed from a chare
      for (const auto& p : psup) m_psup[ p.first ] = p.second;
      if (complete()) nz_complete();
    }

  private:
    HostProxy m_host;           //!< Host proxy
    std::size_t m_chunksize;    //!< Number of rows the first npe-1 PE own
    std::size_t m_lower;        //!< Lower index of the global rows for my PE
    std::size_t m_upper;        //!< Upper index of the global rows for my PE
    tk::hypre::HypreMatrix m_A; //!< Hypre matrix
    std::map< std::size_t, std::vector< std::size_t > > m_psup; //! Nonzeros

    //! Check if our portion of the matrix non-zero structure is complete
    //! \return True if our portion of the distributed matrix is complete
    //! \details Since we known what rows we own, there is no need to explicitly
    //!   send and receive information rows imported. This function is used to
    //!   test whether we have received all of the non-zero structure that we
    //!   own in the distributed matrix.
    bool complete() {
      if ( m_psup.size() == m_upper-m_lower &&
            m_psup.begin()->first == m_lower &&
          (--m_psup.end())->first == m_upper-1 )
        return true;
      else
        return false;
    }

    // Create and set non-zero matrix values
    void fillMatrix() {
      Assert( !m_psup.empty(),
              "Distributed matrix empty on PE " + std::to_string( CkMyPe() ) );
      Assert( complete(),
              "Lower and upper row indices of distributed matrix on PE " +
              std::to_string( CkMyPe() ) + " is incomplete" );
      // Build Hypre data for my PE
      std::vector< int > rows;    // Row indices from my PE
      std::vector< int > ncols;   // Number of matrix columns/row from my PE
      std::vector< int > cols;    // Matrix column indices for rows from my PE
      for (const auto& p : m_psup) {
        rows.push_back( static_cast< int >( p.first ) );
        ncols.push_back( static_cast< int >( p.second.size()+1 ) );
        cols.push_back( static_cast< int >( p.first ) );
        for (auto c : p.second) cols.push_back( static_cast< int >( c ) );
      }
      std::vector< tk::real > vals;
      for (auto v : cols) vals.push_back( 1 );
      // Set our portion of the matrix values
      m_A.set( static_cast< int >( m_upper - m_lower ),
               ncols.data(),
               rows.data(),
               cols.data(),
               vals.data() );
      // Activate SDAG trigger signaling that our matrix part has been filled
      filled();
    }

    //! Assemble distributed matrix
    void assemble() {
      m_A.assemble();
      assembled();
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
      Group::contribute(
        CkCallback( inciter::CkIndex_Conductor::redn_wrapper_init(NULL), host )
      );
    }
};

} // tk::

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#define CK_TEMPLATES_ONLY
#include <linsysmerger.def.h>
#undef CK_TEMPLATES_ONLY

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // LinSysMerger_h
