//******************************************************************************
/*!
  \file      src/Inciter/Partitioner.h
  \author    J. Bakosi
  \date      Wed 02 Dec 2015 09:42:45 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Charm++ chare partitioner group used to perform mesh partitioning
  \details   Charm++ chare partitioner group used to parform mesh partitioning.
*/
//******************************************************************************
#ifndef Partitioner_h
#define Partitioner_h

#include <unordered_map>
#include <numeric>

#include "ExodusIIMeshReader.h"
#include "ContainerUtil.h"
#include "ZoltanInterOp.h"
#include "Inciter/InputDeck/InputDeck.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "partitioner.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::unordered_map< int, std::vector< std::size_t > > g_element;

//! Partitioner Charm++ chare group class
//! \details Instantiations of Partitioner comprise a processor aware Charm++
//!   chare group. When instantiated, a new object is created on each PE and not
//!   more (as opposed to individual chares or chare array object elements). See
//!   also the Charm++ interface file partitioner.ci. The class is templated so
//!   that the same code (parameterized by the template arguments) can be
//!   generated for interacting with different types of Charm++ proxies.
//! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html
//! \author J. Bakosi
template< class HostProxy >
class Partitioner : public CBase_Partitioner< HostProxy > {

  private:
    using Group = CBase_Partitioner< HostProxy >;

  public:
    //! Constructor
    //! \param[in] hostproxy Host Charm++ proxy we are being called from
    Partitioner( HostProxy& host ) : m_host( host ), m_npe( 0 ) {
      // Create ExodusII mesh reader
      tk::ExodusIIMeshReader
        er( g_inputdeck.get< tag::cmd, tag::io, tag::input >() );
      // Read our chunk of the mesh graph from file
      readGraph( er );
      // If a geometric partitioner is selected, compute element centroid
      // coordinates
      const auto alg = g_inputdeck.get< tag::selected, tag::partitioner >();
      if ( alg == tk::ctr::PartitioningAlgorithmType::RCB ||
           alg == tk::ctr::PartitioningAlgorithmType::RIB )
        computeCentroids( er );
      else
        signal2host_setup_complete( m_host );
    }

    //! Partition the computational mesh
    //! \param[in] nchare Number of parts the mesh will be partitioned into
    void partition( int nchare ) {
      m_nchare = nchare;
      const auto alg = g_inputdeck.get< tag::selected, tag::partitioner >();
      const auto che = tk::zoltan::geomPartMesh( alg,
                                                 m_centroid,
                                                 m_gelemid,
                                                 m_tetinpoel.size()/4,
                                                 nchare );
      Assert( che.size() == m_tetinpoel.size()/4, "Size of ownership array "
              "does not equal the number of mesh graph elements" );
      // Construct global mesh element ids for each chare
      m_partelem = elemOwner( che );
      g_element = m_partelem;
      distribute();
    }

    //! Reorder global mesh node IDs
    void reorder() {
      int chunksize, mynchare;
      std::tie( chunksize, mynchare ) = chareDistribution();
      Assert(m_elem.size() == mynchare, "Global mesh element ids associated "
             "to chares on PE " + std::to_string( CkMyPe() ) + " is incomplete");

//       for (const auto& c : m_elem)
//         std::cout << CkMyPe() << ": " << c.first << ", " << c.second.size() << std::endl;

      signal2host_reordering_complete( m_host );
    }

  private:
    //! Host proxy
    HostProxy m_host;
    //! Number of fellow PEs to send elem IDs to
    std::size_t m_npe;
    //! Tetrtahedron element connectivity of our chunk of the mesh
    std::vector< std::size_t > m_tetinpoel;
    //! Global element IDs we read (our chunk of the mesh)
    std::vector< std::size_t > m_gelemid;
    //! Element centroid coordinates of our chunk of the mesh
    std::array< std::vector< tk::real >, 3 > m_centroid;
    //! Total number of chares across all PEs
    int m_nchare;
    //! \brief Global mesh element ids associated to chare IDs on this PE
    //!   resulting from partitioning
    std::unordered_map< int, std::vector< std::size_t > > m_partelem;
    //! Global mesh element ids associated to chares owned by this PE
    std::unordered_map< int, std::vector< std::size_t > > m_elem;

    //! Read mesh graph from file, a chunk by each PE group
    void readGraph( tk::ExodusIIMeshReader& er )
    {
      // Get number of mesh points and number of tetrahedron elements in file
      er.readElemBlockIDs();
      auto nel = er.nel( tk::ExoElemType::TET );
      // Read our chunk of tetrahedron element connectivity from file and also
      // store the list of global element indices for our chunk of the mesh
      auto chunk = nel / CkNumPes();
      auto from = CkMyPe() * chunk;
      auto till = from+chunk;
      if (CkMyPe() == CkNumPes()-1) till += nel % CkNumPes();
      std::array< std::size_t, 2 > ext = { { static_cast<std::size_t>(from),
                                             static_cast<std::size_t>(till) } };
      er.readElements( ext, tk::ExoElemType::TET, m_tetinpoel );
      m_gelemid.resize( static_cast<std::size_t>(till-from) );
      std::iota( begin(m_gelemid), end(m_gelemid), from );
      signal2host_graph_complete( m_host, m_tetinpoel.size()/4 );
    }

    // Compute element centroid coordinates
    void computeCentroids( tk::ExodusIIMeshReader& er )
    {
      auto gid = m_tetinpoel;
      tk::unique( gid );
      // Read node coordinates of our chunk of the mesh elements from file
      auto coord = er.readNodes( tk::extents(gid) );
      // Make room for element centroid coordinates
      auto& cx = m_centroid[0];
      auto& cy = m_centroid[1];
      auto& cz = m_centroid[2];
      auto num = m_tetinpoel.size()/4;
      cx.resize( num );
      cy.resize( num );
      cz.resize( num );
      // Compute element centroids for our chunk of the mesh elements
      for (std::size_t e=0; e<num; ++e) {
        const auto& a = tk::cref_find( coord, m_tetinpoel[e*4+0] );
        const auto& b = tk::cref_find( coord, m_tetinpoel[e*4+1] );
        const auto& c = tk::cref_find( coord, m_tetinpoel[e*4+2] );
        const auto& d = tk::cref_find( coord, m_tetinpoel[e*4+3] );
        cx[e] = (a[0] + b[0] + c[0] + d[0]) / 4.0;
        cy[e] = (a[1] + b[1] + c[1] + d[1]) / 4.0;
        cz[e] = (a[2] + b[2] + c[2] + d[2]) / 4.0;
      }
      signal2host_setup_complete( m_host );
    }

    //! Construct global mesh element ids for each chare
    //! \param[in] che Chares of elements: array of chare ownership IDs mapping
    //!   graph elements to Charm++ chares. Size: number of elements in the
    //!   chunk of the mesh graph on this PE.
    //! \return Vector of global mesh element ids owned by each chare
    //! \note The chare IDs, as keys in the map constructed here, are simply the
    //!   chare IDs returned by the partitioner assigning mesh elements to these
    //!   chares. It does not mean these chare IDs are owned on this PE.
    //!   Distributing the chare IDs and their associated global element IDs to
    //!   their owner PEs are done by the member function distribute().
    std::unordered_map< int, std::vector< std::size_t > >
    elemOwner( const std::vector< std::size_t >& che ) const {
      Assert( che.size() == m_gelemid.size(), "The size of the global element "
              "index and the chare element arrays must equal" );

      std::unordered_map< int, std::vector< std::size_t > > element;

      for (std::size_t e=0; e<che.size(); ++e)
        element[ static_cast<int>(che[e]) ].push_back( m_gelemid[e] );

      Assert( !element.empty(),
              "No elements assigned to chares on one of the MPI ranks" );

      // This check should always be done, as it can result from incorrect user
      // input compared to the mesh size and not due to programmer error.
      for(const auto& c : element)
        ErrChk( !c.second.empty(),
                R"(Overdecomposition of the mesh is too large compared to the
                   number of work units computed based on the degree of
                   virtualization desired. As a result, there would be at least
                   one work unit with no mesh elements to work on, i.e., nothing
                   to do. Solution 1: decrease the virtualization to a lower
                   value using the command-line argument '-u'. Solution 2:
                   decrease the number processing elements (PEs) using the
                   charmrun command-line argument '+pN' where N is the number of
                   PEs, which implicitly increases the size (and thus decreases
                   the number) of work units.)" );

      return element;
    } 

    //! Distribute global mesh element IDs to their owner PEs
    //! \details Chare ids are distributed to PEs in a linear continguous order
    //!   with the last PE taking the remainder if the number of PEs is not
    //!   divisible by the number chares. For example, if nchare=7 and npes=3,
    //!   the chare distribution is PE0: 0 1, PE1: 2 3, and PE2: 4 5 6.
    //!   As a result of this distribution, all PEs will have in their m_elem
    //!   map filled with the global mesh elements IDs associated to the Charm++
    //!   chare IDs each PE owns.
    void distribute() {
      int chunksize, mynchare;
      std::tie( chunksize, mynchare ) = chareDistribution();
      for (int c=0; c<mynchare; ++c) {
        auto cid = CkMyPe() * chunksize + c;  // compute owned chare ID
        const auto it = m_partelem.find(cid); // attempt to find its elements
        if (it != end(m_partelem)) {          // if found
          m_elem.insert( *it );               // move over owned key-value pairs
          m_partelem.erase( it );             // remove chare ID and elements
        }
      }
      // Construct export map associating those map entries (mesh element
      // indices associated to chare IDs) owned by chares we do not own
      std::unordered_map< int,
        std::unordered_map< int, std::vector< std::size_t > > > exp;
      for (const auto& c : m_partelem) exp[ pe(c.first) ].insert( c );
      // Export chare IDs and element IDs we do not own to fellow PEs
      m_npe = exp.size();
      for (const auto& p : exp)
        Group::thisProxy[ p.first ].add( CkMyPe(), p.second );
      if (recvaliens()) signal2host_distribution_complete( m_host );
    }

  public:
    //! Receive mesh element indices associated to chares we own
    //! \param[in] element Mesh element indices associated to chare IDs
    void add( int frompe,
             const std::unordered_map< int, std::vector< std::size_t > >& elem )
    {
      for (const auto& c : elem) {
        Assert( pe(c.first) == CkMyPe(), "PE " + std::to_string(CkMyPe()) +
                " received a chareid-elemidx-vector pair whose chare it does"
                " not own" );
        m_elem[ c.first ] = c.second;
      }
      Group::thisProxy[ frompe ].recv();
    }

    //! Acknowledge received element IDs
    void recv() {
      --m_npe;
      if (recvaliens()) signal2host_distribution_complete( m_host );
    }

  private:

    std::pair< int, int >
    chareDistribution() const {
      auto chunksize = m_nchare / CkNumPes();
      auto mynchare = chunksize;
      if (CkMyPe() == CkNumPes()-1) mynchare += m_nchare % CkNumPes();
      return { chunksize, mynchare };
    }

    //! Test if all fellow PEs have received my element ids contributions
    bool recvaliens() const { return m_npe == 0; }

    //! Return processing element for chare id
    //! \param[in] id Chare id
    //! \return PE that creates the chare
    int pe( int id ) const {
      auto pe = id / (m_nchare / CkNumPes());
      if (pe == CkNumPes()) --pe;
      return pe;
    }

    //! Signal back to host that we have done our part of reading the mesh graph
    void signal2host_graph_complete( const CProxy_Conductor& host,
                                     uint64_t nelem )
    {
      Group::contribute(sizeof(uint64_t), &nelem, CkReduction::sum_int,
                        CkCallback( CkReductionTarget(Conductor,load), host ));
    }
    //! Signal back to host that we are ready for partitioning the mesh
    void signal2host_setup_complete( const CProxy_Conductor& host ) {
      Group::contribute(
        CkCallback(CkIndex_Conductor::redn_wrapper_partition(NULL), m_host ));
    }
    //! \brief Signal back to host that we have done our part of distributing
    //!   mesh element IDs after partitioning
    void signal2host_distribution_complete( const CProxy_Conductor& host ) {
      Group::contribute(
        CkCallback(CkIndex_Conductor::redn_wrapper_reorder(NULL), m_host ));
    }
    //! Signal back to host that we have done our part of reordering the mesh
    void signal2host_reordering_complete( const CProxy_Conductor& host ) {
      Group::contribute(
        CkCallback(CkIndex_Conductor::redn_wrapper_prepared(NULL), m_host ));
    }
};

} // inciter::

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#define CK_TEMPLATES_ONLY
#include "partitioner.def.h"
#undef CK_TEMPLATES_ONLY

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // Partitioner_h
