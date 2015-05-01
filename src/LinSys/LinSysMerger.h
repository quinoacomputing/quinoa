//******************************************************************************
/*!
  \file      src/LinSys/LinSysMerger.h
  \author    J. Bakosi
  \date      Thu 30 Apr 2015 10:54:25 PM MDT
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
    LinSysMerger( HostProxy& host ) :
      m_host( host ),
      m_lower( std::numeric_limits< std::size_t >::max() ),
      m_upper( 0 ),
      m_nCharePes( 0 ),
      m_nPePts( 0 ) {}

    //! Register chare on my PE to expect contribution from
    //! \param[in] id Charm++ chare array index the contribution coming from
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void checkin( int id )
    { m_chare.insert( static_cast< std::size_t >( id ) ); }

    //! \brief Function called by chares on my PE to contribute their part of
    //!   the structure of the global distributed matrix and vector
    //! \param[in] workerarray Charm++ proxy of the array of worker chares.
    //!    Also templated on the worker proxy type so the same code can be
    //!    generated for different types of workers and the linear system merger
    //!    is functionally separate from workers that the merger integrates.
    //! \param[in] lower Lower global index (inclusive) of the merged system
    //! \param[in] upper Upper global index (exclusive) of the merged system
    //! \param[in] exp Chare-export map associating chares to mesh points sent
    //! \param[in] point Global mesh point ids sent by chares on our PE
    //! \param[in] psup Points surrounding points, see tk::genPsup()
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    template< class WorkProxy >
    void structure( WorkProxy& workerarray,
                    std::size_t lower,
                    std::size_t upper,
                    const std::map< std::size_t,
                                    std::vector< std::size_t > >& exp,
                    const std::vector< std::size_t >& point,
                    const std::pair< std::vector< std::size_t >,
                                     std::vector< std::size_t > >& psup )
    {
      // Add global mesh point ids sent by chare on our PE to local store
      m_point.reserve( m_point.size() + point.size() );
      m_point.insert( end(m_point), begin(point), end(point) );
      std::cout << "PE " << CkMyPe() << " chares: ";
      for (auto c : m_chare) std::cout << c << " ";
      std::cout << " ci: [" << lower << "..."
                << upper << ") " << "npoin: " << point.size() << "\n";
      Assert( point.size() == psup.second.size()-1,
              "Number of owned points must equal in global id vector and "
              "derived data, points surrounding points, sent by chare." );
      // Find lower and upper row indices for my PE
      if (lower < m_lower) m_lower = lower;
      if (upper > m_upper) m_upper = upper;
      // Activate SDAG-wait for finishing initialization
      wait4init();
      // Build chare-export map for my PE. This is similar to the export map of
      // chares (see e.g., inciter::Performer::m_export), but only contains
      // export information needed between PEs. The chare-export map associates
      // chares not owned by my PE to a vector of point ids sent. To build the
      // map we loop through chares the caller of this function (a worker chare)
      // exports to and find out if a chare that it exports to is owned by
      // another PE; if so, store its point ids sent. This is only the first
      // step in constructing the export map we need, since ultimately we need
      // PEs associated to exported point ids. Thus we also construct another
      // map associating chare ids (we need to send points to) to PEs which will
      // be used to address the chare-export map in a PE-export map fashion.
      for (const auto& c : exp)
        if (m_chare.find( c.first ) == m_chare.end())
          for (auto p : c.second) m_chExport[ c.first ].push_back( p );
      // This is separated from the above loop so that counting the number of
      // worker chares expected whose PEs we need to export to are final
      // before we start querying the chares' PEs. If these async calls were
      // fired up in the above loop, together with counting up the number,
      // there could be a chance that workPe() prematurely stops waiting for
      // some chares.
      if (m_chExport.empty()) {  // no PEs to export to, our work is done here
        pointsExported(); 
      } else
        for (const auto& c : exp)
          if (m_chare.find( c.first ) == m_chare.end())
            workerarray[ static_cast<int>(c.first) ].pe( CkMyPe() );

      for (std::size_t p=0; p<psup.second.size()-1; ++p) {
        //auto lid = m_lower + p;
        auto gid = point[ p ];
        m[ gid ].push_back( gid );
        for (auto i=psup.second[p]+1; i<=psup.second[p+1]; ++i)
          m[ gid ].push_back( psup.first[i] );
      }
    }

    //! Receive and associate worker's chare id to its PE we need to export to
    //! \param[in] ch Chare worker array index PE coming from
    //! \param[in] pe PE the chare happens to be on
    void workPe( int ch, int pe ) {
      ++m_nCharePes;
      // Construct PE-export map from chare-export map. That is, create a new
      // map that associates mesh points sent to PEs based on the map
      // associating mesh points sent to chares. This will allow fewer and
      // larger messages when my PE's contributions will be sent to fellow group
      // branches, since we will send a single vector of mesh points to a PE
      // instead of several smaller vectors per chare contributing.
      for (const auto& c : m_chExport)
        if (c.first == ch)
          for (auto p : c.second)
            m_peExport[ static_cast< std::size_t >( pe ) ].push_back( p );
      // If all worker chares whose PEs we need to export to have contributed
      // their PE, make PE-export vectors unique and activate SDAG trigger.
      if (m_nCharePes == m_chExport.size()) {
        // Make PE-export map of vectors contain only unique point ids
        for (auto& c : m_peExport) tk::unique( c.second );
        // Activate SDAG trigger signaling that all worker chares whose PEs we
        // need to export have contributed their PE
//         std::cout << "peEx on " << CkMyPe() << ": ";
//         for (const auto& b : m_peExport) {
//           std::cout << b.first << "( ";
//           for (auto p : b.second) std::cout << p << " ";
//           std::cout << ") ";
//         }
//         std::cout << std::endl;
        // Send mesh point ids to fellow group branches my PE contributes to
        exportPoints();
      }
    }

    void createLinearSystem() {
//       if (!m_chExport.empty()) {
//         std::cout << "chEx on " << CkMyPe() << " c: ";
//         for (const auto& l : m_chExport) {
//           std::cout << "(" << l.first << ") ";
//           for (auto c : l.second) std::cout << c << " ";
//         }
//         std::cout << std::endl;
//       }
      // Create my PE's part of the lhs matrix distributed across all PEs
      m_A.create( m_lower, m_upper );
      // Activate SDAG trigger signaling that our part of the linear system has
      // been created
      systemCreated();
     }

    //! Count number of PEs we exported mesh points to (signaling back 'roger')
    void rgr() {
      ++m_nPePts;
      // If all PEs we export mesh points to have successfully received their
      // points, activate SDAG trigger signaling that this task is completed
      if (m_nPePts == m_peExport.size()) {
//         std::sort( begin(m_point), end(m_point) );
//         std::cout << "final gids on PE " << CkMyPe() << " (" << m_point.size()
//                   << ") : ";
//         for (auto p : m_point) std::cout << p << " ";
//         std::cout << std::endl;
        pointsExported();
      }
    }

    //! Receive mesh point ids from fellow group branches
    //! \param[in] caller Caller PE index we send back 'roger' to
    //! \param[in] transfer Mesh point ids received
    void add( int caller, const std::vector< std::size_t >& transfer )
    {
//       std::cout << "points import on PE " << CkMyPe() << " from " << caller
//                 << ": ";
//       for (auto p : transfer) std::cout << p << " ";
//       std::cout << std::endl;
//       // Add global mesh point ids sent by fellow group branch to local store
//       m_point.reserve( m_point.size() + transfer.size() );
//       m_point.insert( end(m_point), begin(transfer), end(transfer) );
//       std::cout << "now gids on PE " << CkMyPe() << ": ";
//       for (auto p : m_point) std::cout << p << " ";
//       std::cout << std::endl;
      // Send Roger, i.e., "I have received all of the last transmission" back
      Group::thisProxy[ caller ].rgr();
    }

    //! Assemble matrix
    void assemble() {
      m_A.assemble();
      // Activate SDAG trigger signaling that the matrix has been assembled
      assembled();
    }

  private:
    HostProxy m_host;           //!< Host proxy

    std::size_t m_lower;        //!< Lower index of the global rows for my PE
    std::size_t m_upper;        //!< Upper index of the global rows for my PE

    //! Global point ids of owned and received mesh points
    std::vector< std::size_t > m_point;

    //! Number of worker chares contributed whose PEs we need to export to
    int m_nCharePes;
    //! Number of fellow group branches that have received points from us
    int m_nPePts;

    //! Unique chare IDs on my PE
    std::set< std::size_t > m_chare;

    //! Export map associating chares not owned by my PE to points sent
    std::map< std::size_t, std::vector< std::size_t > > m_chExport;
    //! Export map associating PEs to points sent
    std::map< std::size_t, std::vector< std::size_t > > m_peExport;

    tk::hypre::HypreMatrix m_A; //!< Hypre matrix

    std::map< std::size_t, std::vector< std::size_t > > m;

    //! Send mesh point ids to fellow group branches my PE contributes to
    void exportPoints() {
      // Bool array to store PEs we do not export to
      std::vector< bool > pes( static_cast<std::size_t>(CkNumPes()), true );
      //pes[ static_cast< std::size_t >( CkMyPe() ) ] = false; // no expt to self
      // Send vector of mesh points to all PEs we export to
      for (const auto& c : m_peExport) {
        Group::thisProxy[ static_cast<int>(c.first) ].add( CkMyPe(), c.second );
        pes[ c.first ] = false;
      }
      // Activate SDAG trigger on PEs we don't export to signaling that "points
      // have been exported"
      for (std::size_t p=0; p<pes.size(); ++p)
        if (pes[p]) {
          //std::cout << CkMyPe() << ": does not send to " << p << std::endl;
          Group::thisProxy[ static_cast<int>(p) ].pointsExported();
        }
    }

    // Set non-zero matrix values for my PE
    void fillMatrix() {
      std::cout << "pi: [" << m_lower << "..." << m_upper << ")\n";
      std::cout << "m on PE " << CkMyPe() << ": ";
      for (const auto& l : m) {
        std::cout << "(" << l.first << ") ";
        for (auto c : l.second) std::cout << c << " ";
      }
      std::cout << std::endl;
      for (const auto& l : m)
        Assert( l.first >= m_lower || l.first < m_upper,
                "Matrix contributions collected from chares on PE " +
                std::to_string( CkMyPe() ) + " are out of bounds" );

      // Build Hypre data for my PE
      std::vector< int > rows( m_upper - m_lower, 0 ); // Row indices from my PE
      std::iota( begin(rows), end(rows), static_cast<int>(m_lower) );
      //for (auto i=m_lower; i<m_upper; ++i)
      //  rows.push_back( static_cast< int >( i ) );
      std::vector< int > cols;    // Matrix column indices for rows from my PE
      std::vector< int > ncols;   // Number of matrix columns/row from my PE
      for (auto& l : m) {
        tk::unique( l.second );
        ncols.push_back( static_cast< int >( l.second.size() ) );
        for (auto c : l.second) cols.push_back( static_cast< int >( c ) );
      }
      int i = 0;
      std::vector< tk::real > vals;
      for (auto v : cols) vals.push_back( static_cast< tk::real >( i++ ) );
      // Set our portion of the matrix values
      m_A.set( static_cast< int >( m_upper-m_lower),
               ncols.data(),
               rows.data(),
               cols.data(),
               vals.data() );
      // Activate SDAG trigger signaling that our matrix part has been filled
      filled();
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
