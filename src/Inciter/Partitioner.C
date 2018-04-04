// *****************************************************************************
/*!
  \file      src/Inciter/Partitioner.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ chare partitioner group used to perform mesh partitioning
  \details   Charm++ chare partitioner group used to parform mesh partitioning.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation, communication as well as I/O. The
    algorithm utilizes the structured dagger (SDAG) Charm++ functionality. The
    high-level overview of the algorithm structure and how it interfaces with
    Charm++ is discussed in the Charm++ interface file
    src/Inciter/partitioner.ci. See also src/Inciter/Partitioner.h for the
    asynchronous call graph.
*/
// *****************************************************************************

#include "Partitioner.h"
#include "DerivedData.h"
#include "Reorder.h"
#include "Inciter/Options/Scheme.h"
#include "AMR/mesh_adapter.h"
#include "MeshReader.h"
#include "Around.h"
#include "ExodusIIMeshWriter.h"
#include "UnsMesh.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::Partitioner;

Partitioner::Partitioner(
  const std::vector< CkCallback >& cb,
  const CProxy_Transporter& host,
  const tk::CProxy_Solver& solver,
  const CProxy_BoundaryConditions& bc,
  const Scheme& scheme,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel ) :
  m_cb( cb[0], cb[1], cb[2], cb[3], cb[4], cb[5], cb[6], cb[7] ),
  m_host( host ),
  m_solver( solver ),
  m_bc( bc ),
  m_scheme( scheme ),
  m_npe( 0 ),
  m_reqNodes(),
  m_reqEdges(),
  m_start( 0 ),
  m_noffset( 0 ),
  m_nquery( 0 ),
  m_nmask( 0 ),
  m_tetinpoel(),
  m_gelemid(),
  m_coord(),
  m_centroid(),
  m_nchare( 0 ),
  m_lower( 0 ),
  m_upper( 0 ),
  m_ncomm(),
  m_ecomm(),
  m_ncommunication(),
  m_ecommunication(),
  m_nodeset(),
  m_edgeset(),
  m_linnodes(),
  m_linedges(),
  m_chinpoel(),
  m_chfilenodes(),
  m_chedgenodes(),
  m_cost( 0.0 ),
  m_bnodechares(),
  m_edgechares(),
  m_msum(),
  m_msumed(),
  m_bface( bface ),
  m_triinpoel( triinpoel )
// *****************************************************************************
//  Constructor
//! \param[in] cb Charm++ callbacks
//! \param[in] host Host Charm++ proxy we are being called from
//! \param[in] solver Linear system solver proxy
//! \param[in] bc Boundary conditions group proxy
//! \param[in] scheme Discretization scheme
//! \param[in] bface Face lists mapped to side set ids
//! \param[in] triinpoel Interconnectivity of points and boundary-face
// *****************************************************************************
{
  // Create mesh reader
  MeshReader mr( g_inputdeck.get< tag::cmd, tag::io, tag::input >(),
                 static_cast< std::size_t >( CkNumPes() ),
                 static_cast< std::size_t >( CkMyPe() ) );

  // Read our chunk of the mesh graph from file
  mr.readGraph(  m_gelemid, m_tetinpoel );

  // Send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.peread();

  // Sum number of elements across all PEs (will define total load)
  uint64_t nelem = m_gelemid.size();
  contribute( sizeof(uint64_t), &nelem, CkReduction::sum_int,
              m_cb.get< tag::load >() );

  // Compute local from global mesh data
  auto el = tk::global2local( m_tetinpoel );
  const auto& inp = std::get< 0 >( el );        // Local mesh connectivity
  const auto& gid = std::get< 1 >( el );        // Local->global node IDs
  const auto& lid = std::get< 2 >( el );        // Global->local node IDs

  // Read our chunk of the mesh node coordinates from file
  m_coord = mr.readCoords( gid );

  // Optionally refine mesh if requested
  refine( inp );

  // Compute cell centroids if a geometric partitioner is selected
  computeCentroids( lid );
}

void
Partitioner::partition( int nchare )
// *****************************************************************************
//  Partition the computational mesh
//! \param[in] nchare Number of parts the mesh will be partitioned into
// *****************************************************************************
{
  m_nchare = nchare;
  const auto alg = g_inputdeck.get< tag::selected, tag::partitioner >();
  const auto che = tk::zoltan::geomPartMesh( alg,
                                             m_centroid,
                                             m_gelemid,
                                             m_tetinpoel.size()/4,
                                             nchare );

  // send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pepartitioned();

  Assert( che.size() == m_gelemid.size(), "Size of ownership array does "
          "not equal the number of mesh graph elements" );

  // Construct global mesh node ids for each chare and distribute
  distribute( chareNodes(che) );

  // Free storage of element connectivity, element centroids, and element
  // IDs as they are no longer needed after the mesh partitioning.
  tk::destroy( m_gelemid );
  tk::destroy( m_centroid );
}

void
Partitioner::offset( int p, std::size_t u )
// *****************************************************************************
//  Receive number of uniquely assigned global mesh node IDs from lower PEs
//! \param[in] p PE ID
//! \param[in] u Number of mesh node IDs PE p will assign IDs to
//! \details This function computes the offset each PE will need to start
//!   assigning its new node IDs from (for those nodes that are not assigned
//!   new IDs by any PEs with lower indices). The offset for a PE is the
//!   offset for the previous PE plus the number of node IDs the previous PE
//!   (uniquely) assigns new IDs for minus the number of node IDs the
//!   previous PE receives from others (lower PEs). This is computed here in
//!   a parallel/distributed fashion by each PE sending its number of node
//!   IDs (that it uniquely assigns) to all PEs. Note that each PE would
//!   only need to send this information to higher PEs, but instead this
//!   function is called in a broadcast fashion, because that is more
//!   efficient than individual calls to only the higher PEs. Therefore when
//!   computing the offsets, we only count the lower PEs. When this is done,
//!   we have the precise communication map as well as the start offset on
//!   all PEs and so we can start the distributed global mesh node ID
//!   reordering.
// *****************************************************************************
{
  if (p < CkMyPe()) m_start += u;
  if (++m_noffset == static_cast<std::size_t>(CkNumPes())) reorder();
}


void
Partitioner::request( int p, const std::unordered_set< std::size_t >& nd )
// *****************************************************************************
//  Request new global node IDs for old node IDs
//! \param[in] p PE request coming from and to which we send new IDs to
//! \param[in] nd Set of old node IDs whose new IDs are requested
// *****************************************************************************
{
  // Queue up requesting PE and node IDs
  m_reqNodes.push_back( { p, nd } );
  // Trigger SDAG wait signaling that node IDs have been requested from us
  nodes_requested_complete();
}

void
Partitioner::request( int p, const tk::UnsMesh::Edges& ed )
// *****************************************************************************
//  Request new global node IDs for edges
//! \param[in] p PE request coming from and to which we send new IDs to
//! \param[in] ed Set of edges whose new IDs are requested
// *****************************************************************************
{
  // Queue up requesting PE and node IDs
  m_reqEdges.push_back( { p, ed } );
  // Trigger SDAG wait signaling that node IDs have been requested from us
  nodes_requested_complete();
}

void
Partitioner::neworder(const std::unordered_map< std::size_t, std::size_t >& nd)
// *****************************************************************************
//  Receive new (reordered) global node IDs
//! \param[in] nd Map associating new to old node IDs
// *****************************************************************************
{
  // Signal to the runtime system that we have participated in reordering
  participated_complete();
  // Store new node IDs associated to old ones
  for (const auto& p : nd) m_linnodes[ p.first ] = p.second;
  // If all our nodes have new IDs assigned, signal that to the runtime
  if (m_linnodes.size() == m_nodeset.size()) nodesreorder_complete();
}

void
Partitioner::neworder( const tk::UnsMesh::EdgeNodes& ed )
// *****************************************************************************
//  Receive new global node IDs associated to edge-nodes
//! \param[in] ed Map associating node IDs to edges
// *****************************************************************************
{
  // Signal to the runtime system that we have participated in reordering
  participated_complete();
  // Store node IDs associated to edge
  for (const auto& e : ed) m_linedges[ e.first ] = e.second;
  // If all our edges have new IDs assigned, signal that to the runtime
  if (m_linedges.size() == m_edgeset.size()) edgesreorder_complete();
}

void
Partitioner::add( int frompe,
  const std::unordered_map< int, std::vector< std::size_t > >& n )
// *****************************************************************************
//! Receive mesh node IDs associated to chares we own
//! \param[in] n Mesh node indices associated to chare IDs
//! \param[in] frompe PE call coming from
// *****************************************************************************
{
  for (const auto& c : n) {
    Assert( pe(c.first) == CkMyPe(), "PE " + std::to_string(CkMyPe()) +
            " received a chareid-nodeidx-vector pair whose chare it does"
            " not own" );
    auto& inpoel = m_chinpoel[ c.first ];
    inpoel.insert( end(inpoel), begin(c.second), end(c.second) );
  }

  thisProxy[ frompe ].recv();
}

void
Partitioner::recv()
// *****************************************************************************
//  Acknowledge received node IDs
// *****************************************************************************
{
  if (--m_npe == 0) contribute( m_cb.get< tag::distributed >() );
}

void
Partitioner::flatten()
// *****************************************************************************
//  Prepare owned mesh node IDs for reordering
//! \details The 'flatten' is used here as a concatenation of a data
//!   structure that stores date categorized by chares owned on this PE. The
//!   result of the flattening is thus a simpler data structure that is no
//!   longer categorized by (or associated to) chares.
// *****************************************************************************
{
  // Make sure we are not fed garbage
  Assert( m_chinpoel.size() ==
            static_cast< std::size_t >( chareDistribution()[1] ),
          "Global mesh nodes ids associated to chares on PE " +
          std::to_string( CkMyPe() ) + " is incomplete" );

  // Collect chare IDs we own associated to edges
  for (const auto& c : m_chedgenodes)
    for (const auto& e : c.second)
      m_edgechares[ e.first ].push_back( c.first );
  // Make chare IDs (associated to edges) unique
  for (auto& c : m_edgechares) tk::unique( c.second );

  // Flatten node IDs of elements our chares operate on
  for (const auto& c : m_chinpoel)
    for (auto i : c.second)
      m_nodeset.insert( i );

  // Flatten edges of elements our chares operate on
  for (const auto& c : m_chedgenodes)
    for (const auto& e : c.second)
      m_edgeset.insert( e.first );

  // Find chare-boundary nodes of all chares on this PE
  for (const auto& c : m_chinpoel) {    // for all chare connectivities
    // generate local ids and connectivity from global connectivity
    auto el = tk::global2local( c.second );
    const auto& inpoel = std::get< 0 >( el );   // local connectivity
    const auto& gid = std::get< 1 >( el );      // local->global node ids
    auto esup = tk::genEsup( inpoel, 4 );       // elements surrounding points
    auto esuel = tk::genEsuelTet( inpoel, esup ); // elems surrounding elements
    // collect chare boundary nodes
    for (std::size_t e=0; e<esuel.size()/4; ++e) {
      auto mark = e*4;
      for (std::size_t f=0; f<4; ++f) {
        if (esuel[mark+f] == -1) {
          auto A = gid[ inpoel[ mark+tk::lpofa[f][0] ] ];
          auto B = gid[ inpoel[ mark+tk::lpofa[f][1] ] ];
          auto C = gid[ inpoel[ mark+tk::lpofa[f][2] ] ];
          m_bnodechares[ A ].push_back( c.first );
          m_bnodechares[ B ].push_back( c.first );
          m_bnodechares[ C ].push_back( c.first );
        }
      }
    }
  }
  // Make node chare lists unique
  for (auto& n : m_bnodechares) tk::unique( n.second );

  // send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.peflattened();

  // Signal host that we are ready for computing the communication map,
  // required for parallel distributed global mesh node reordering
  contribute( m_cb.get< tag::flattened >() );
}

void
Partitioner::lower( std::size_t low )
// *****************************************************************************
//  Receive lower bound of node IDs our PE operates on after reordering
//! \param[in] low Lower bound of node IDs assigned to us
// *****************************************************************************
{
  m_lower = low;
  lower_complete();
}

void
Partitioner::stdCost( tk::real av )
// *****************************************************************************
//  Compute the variance of the communication cost
//! \param[in] av Average of the communication cost
//! \details Computing the standard deviation is done via computing and
//!   summing up the variances on each PE and asynchronously reducing the
//!   sum to our host.
// *****************************************************************************
{
  tk::real var = (m_cost-av)*(m_cost-av);
  contribute( sizeof(tk::real), &var, CkReduction::sum_double,
              m_cb.get< tag::stdcost >() );
}

tk::real
Partitioner::cost( std::size_t l, std::size_t u )
// *****************************************************************************
//  Compute communication cost on our PE
//! \param[in] l Lower global node ID this PE works on
//! \param[in] u Upper global node ID this PE works on
//! \return Communication cost for our PE
//! \details The cost is a real number between 0 and 1, defined as the
//!   number of mesh points we do not own, i.e., need to send to some other
//!   PE, divided by the total number of points we contribute to. The lower
//!   the better.
// *****************************************************************************
{
  std::size_t ownpts = 0, compts = 0;
  for (auto p : m_nodeset) if (p >= l && p < u) ++ownpts; else ++compts;

  // Free storage of unique global node IDs chares on our PE will contribute to
  // as it is no longer needed after computing the communication cost.
  tk::destroy( m_nodeset );

  return static_cast<tk::real>(compts) / static_cast<tk::real>(ownpts + compts);
}

void
Partitioner::gather()
// *****************************************************************************
//  Start gathering global node IDs this PE will need to receive (instead of
//  assign) during reordering
// *****************************************************************************
{
  // Exctract boundary chare nodes to vector
  std::vector< std::size_t > bnodes( m_bnodechares.size() );
  std::size_t i = 0;
  for (const auto& n : m_bnodechares) bnodes[ i++ ] = n.first;

  // Send chare boundary nodes to all PEs in a broadcast fashion
  thisProxy.query( CkMyPe(), bnodes );

  // send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pegather();
}

void
Partitioner::query( int p, const std::vector< std::size_t >& bnodes )
// *****************************************************************************
//  Query our global node IDs by other PEs so they know if they are to receive
//  IDs for those from during reordering
//! \param[in] p Querying PE
//! \param[in] nodes List of global mesh node IDs to query
//! \details Note that every PE calls this function in a broadcast fashion,
//!   including our own. However, to compute the correct result, this would
//!   only be necessary for PEs whose ID is higher than ours. However, the
//!   broadcast (calling everyone) is more efficient. This also results in a
//!   simpler logic, because every PE goes through this single call path.
//!   The returned mask is simply a boolean array signaling if the node ID
//!   is found (owned).
// *****************************************************************************
{
  std::unordered_map< std::size_t, std::vector< int > > cn;
  for (auto j : bnodes) {
    const auto it = m_nodeset.find( j );
    if (it != end(m_nodeset)) {
      const auto& c = tk::cref_find( m_bnodechares, j );
      auto& chares = cn[j];
      chares.insert( end(chares), begin(c), end(c) );
    }
  }

  // when we have heard from all PEs, send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() &&
       ++m_nquery == static_cast<std::size_t>(CkNumPes()) )
    m_host.pequery();

  thisProxy[ p ].mask( CkMyPe(), cn );
}

void
Partitioner::mask(
  int p,
  const std::unordered_map< std::size_t, std::vector< int > >& cn )
// *****************************************************************************
//  Receive mask of to-be-received global mesh node IDs
//! \param[in] p The PE uniquely assigns the node IDs marked listed in ch
//! \param[in] cn Vector containing the set of potentially multiple chare
//!   IDs that we own (i.e., contribute to) for all of our node IDs.
//! \details Note that every PE will call this function, since query() was
//!   called in a broadcast fashion and query() answers to every PE once.
//!   This is more efficient than calling only the PEs from which we would
//!   have to receive results from. Thus the incoming results are only
//!   interesting from PEs with lower IDs than ours.
// *****************************************************************************
{
  // Store the old global mesh node IDs associated to chare IDs bordering the
  // mesh chunk held by and associated to chare IDs we own. This loop computes
  // m_msum, a symmetric chare-node communication map, that associates (in its
  // inner map) a unique set of global node IDs to areceiving chare ID, both
  // associated (in its outer map) to a sending chare ID.

  // the mesh chunk held by and associated to chare IDs we own
  for (const auto& h : cn) {
    const auto& chares = tk::ref_find( m_bnodechares, h.first );
    for (auto c : chares) {           // surrounded chares
      auto& sch = m_msum[c];
      for (auto s : h.second)         // surrounding chares
        if (s != c) sch[ s ].insert( h.first );
    }
  }

  // Associate global mesh node IDs to lower PEs we will need to receive from
  // during node reordering. The choice of associated container is std::map,
  // which is ordered (vs. unordered, hash-map). This is required by the
  // following operation that makes the mesh node IDs unique in the
  // communication map. (We are called in an unordered fashion, so we need to
  // collect from all PEs and then we need to make the node IDs unique, keeping
  // only the lowest PEs a node ID is associated with.) Note that m_ncomm is an
  // asymmetric PE-node communication map, associating a set of unique global
  // node IDs to PE IDs from which we (this PE) will need to receive newly
  // assigned node IDs during global mesh node reordering. This map is
  // asymmetric, beacuse of the agreement of the reordering that if a mesh node
  // is shared by multiple PEs, the PE with the lowest ID gets to assign a new
  // ID to it and all others must receive it instead of assigning it.

  if (p < CkMyPe()) {
    auto& id = m_ncomm[ p ];
    for (const auto& h : cn) id.insert( h.first );
  }

  if (++m_nmask == static_cast<std::size_t>(CkNumPes())) {
    // Make sure we have received all we need
    Assert( m_ncomm.size() == static_cast<std::size_t>(CkMyPe()),
            "Communication map size on PE " +
            std::to_string(CkMyPe()) + " must equal " +
            std::to_string(CkMyPe()) );
    // Fill new hash-map, keeping only unique node IDs obtained from the
    // lowest possible PEs
    for (auto c=m_ncomm.cbegin(); c!=m_ncomm.cend(); ++c) {
      auto& n = m_ncommunication[ c->first ];
      for (auto j : c->second)
        if (std::none_of( m_ncomm.cbegin(), c,
             [ j ]( const typename decltype(m_ncomm)::value_type& s )
             { return s.second.find(j) != end(s.second); } )) {
          n.insert(j);
        }
      if (n.empty()) m_ncommunication.erase( c->first );
    }
    // Free storage of temporary communication map used to receive global
    // mesh node IDs as it is no longer needed once the final communication
    // map is generated.
    tk::destroy( m_ncomm );
    // Count up total number of nodes and (nodes associated to edges) we
    // will need receive during reordering
    std::size_t nrecv = 0, erecv = 0;
    for (const auto& u : m_ncommunication) nrecv += u.second.size();

    // send progress report to host
    if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pemask();

    // Compute number of mesh node IDs we will assign IDs to
    auto nuniq = m_nodeset.size() - nrecv + m_edgeset.size() - erecv;

    // Start computing PE offsets for node reordering
    thisProxy.offset( CkMyPe(), nuniq );
  }
}

void
Partitioner::computeCentroids(
  const std::unordered_map< std::size_t, std::size_t >& lid )
// *****************************************************************************
//  Compute element centroid coordinates
//! \param[in] lid Map from global to local node IDs for our chunk of the mesh
// *****************************************************************************
{
  const auto alg = g_inputdeck.get< tag::selected, tag::partitioner >();

  if ( tk::ctr::PartitioningAlgorithm().geometric(alg) ) {
    const auto& x = std::get< 0 >( m_coord );
    const auto& y = std::get< 1 >( m_coord );
    const auto& z = std::get< 2 >( m_coord );
  
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
      auto A = tk::cref_find( lid, m_tetinpoel[e*4+0] );
      auto B = tk::cref_find( lid, m_tetinpoel[e*4+1] );
      auto C = tk::cref_find( lid, m_tetinpoel[e*4+2] );
      auto D = tk::cref_find( lid, m_tetinpoel[e*4+3] );
      cx[e] = (x[A] + x[B] + x[C] + x[D]) / 4.0;
      cy[e] = (y[A] + y[B] + y[C] + y[D]) / 4.0;
      cz[e] = (z[A] + z[B] + z[C] + z[D]) / 4.0;
    }
  }

  contribute( m_cb.get< tag::centroid >() );
}

std::unordered_map< int, std::vector< std::size_t > >
Partitioner::chareNodes( const std::vector< std::size_t >& che ) const
// *****************************************************************************
//  Construct global mesh node ids for each chare
//! \param[in] che Chares of elements: array of chare ownership IDs mapping
//!   graph elements to Charm++ chares. Size: number of elements in the
//!   chunk of the mesh graph on this PE.
//! \return Vector of global mesh node ids connecting elements owned by each
//!   chare on this PE
//! \note The chare IDs, as keys in the map constructed here, are simply the
//!   chare IDs returned by the partitioner assigning mesh elements to these
//!   chares. It does not mean that these chare IDs are owned on this PE.
// *****************************************************************************
{
  Assert( che.size() == m_gelemid.size(), "The size of the global element "
          "index and the chare element arrays must equal" );
  Assert( che.size() == m_tetinpoel.size()/4, "The size of the mesh "
          "connectivity / 4 and the chare element arrays must equal" );

  // Categorize global mesh node ids of elements by chares
  std::unordered_map< int, std::vector< std::size_t > > nodes;
  for (std::size_t e=0; e<che.size(); ++e) {
    auto& c = nodes[ static_cast<int>(che[e]) ];
    for (std::size_t n=0; n<4; ++n) c.push_back( m_tetinpoel[e*4+n] );
  }

  // Make sure all PEs have chares assigned
  Assert( !nodes.empty(), "No nodes have been assigned to chares on PE " +
          std::to_string(CkMyPe()) );

  // This check should always be done, hence ErrChk and not Assert, as it
  // can result from particular pathological combinations of (1) too large
  // degree of virtualization, (2) too many PEs, and/or (3) too small of a
  // mesh and not due to programmer error.
  for(const auto& c : nodes)
    ErrChk( !c.second.empty(),
            "Overdecomposition of the mesh is too large compared to the "
            "number of work units computed based on the degree of "
            "virtualization desired. As a result, there would be at least "
            "one work unit with no mesh elements to work on, i.e., nothing "
            "to do. Solution 1: decrease the virtualization to a lower "
            "value using the command-line argument '-u'. Solution 2: "
            "decrease the number processing elements (PEs) using the "
            "charmrun command-line argument '+pN' where N is the number of "
            "PEs, which implicitly increases the size (and thus decreases "
            "the number) of work units.)" );

  return nodes;
}

void
Partitioner::distribute(
  std::unordered_map< int, std::vector< std::size_t > >&& n )
// *****************************************************************************
//  Distribute global mesh node IDs to their owner PEs
//! \param[in] n Global mesh node IDs connecting elements associated to
//!   chare IDs on this PE resulting from partitioning the mesh elements.
//!   Note that this data is moved in.
// *****************************************************************************
{
  auto dist = chareDistribution();

  for (int c=0; c<dist[1]; ++c) {
    auto chid = CkMyPe() * dist[0] + c;   // compute owned chare ID
    const auto it = n.find( chid );       // attempt to find its nodes
    if (it != end(n)) {                   // if found
      m_chinpoel.insert( *it );           // move over owned key-value pairs
      n.erase( it );                      // remove chare ID and nodes
    }
    Assert( n.find(chid) == end(n), "Not all owned node IDs stored" );
  }

  // Construct export map associating those map entries (mesh node indices
  // associated to chare IDs) owned by chares we do not own. Outer key: PE
  // to export to, inner key: chare ID, value: vector of global node IDs
  std::unordered_map< int,
    std::unordered_map< int, std::vector< std::size_t > > > exp;
  for (auto&& c : n) exp[ pe(c.first) ].insert( std::move(c) );

  // Export chare IDs and node IDs we do not own to fellow PEs
  m_npe = exp.size();
  for (const auto& p : exp)
    thisProxy[ p.first ].add( CkMyPe(), p.second );

  // send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pedistributed();

  if (m_npe == 0) contribute( m_cb.get< tag::distributed >() );
}

std::array< int, 2 >
Partitioner::chareDistribution() const
// *****************************************************************************
//  Compute chare distribution
//! \return Chunksize, i.e., number of chares per all PEs except the last
//!   one, and the number of chares for my PE
//! \details Chare ids are distributed to PEs in a linear continguous order
//!   with the last PE taking the remainder if the number of PEs is not
//!   divisible by the number chares. For example, if nchare=7 and npe=3,
//!   the chare distribution is PE0: 0 1, PE1: 2 3, and PE2: 4 5 6. As a
//!   result of this distribution, all PEs will have their chare-categorized
//!   element connectivity filled with the global mesh node IDs associated
//!   to the Charm++ chare IDs each PE owns.
// *****************************************************************************
{
  auto chunksize = m_nchare / CkNumPes();
  auto mynchare = chunksize;
  if (CkMyPe() == CkNumPes()-1) mynchare += m_nchare % CkNumPes();
  return {{ chunksize, mynchare }};
}

void
Partitioner::reorder()
// *****************************************************************************
//  Reorder global mesh node IDs
// *****************************************************************************
{
  // Activate SDAG waits for having requests arrive from other PEs for some
  // of our node IDs; and for computing/receiving lower and upper bounds of
  // global node IDs our PE operates on after reordering
  wait4prep();
  wait4bounds();
  wait4reorder();

  // In serial signal to the runtime system that we have participated in
  // reordering. This is required here because this is only triggered if
  // communication is required during mesh node reordering. See also
  // particioner.ci.
  if (CkNumPes() == 1) participated_complete();

  // Send out request for new global node IDs for nodes we do not reorder
  for (const auto& c : m_ncommunication)
    thisProxy[ c.first ].request( CkMyPe(), c.second );

  // Send out request for new global node IDs for edges we do not reorder
  for (const auto& e : m_ecommunication)
    thisProxy[ e.first ].request( CkMyPe(), e.second );

  // Lambda to decide if node ID is being assigned a new ID by us
  auto ownnode = [ this ]( std::size_t p ) {
    using Set = typename std::remove_reference<
                  decltype(m_ncommunication) >::type::value_type;
    return !std::any_of( m_ncommunication.cbegin(), m_ncommunication.cend(),
                         [&](const Set& s)
                         { return s.second.find(p) != s.second.cend(); } );
  };

  // Lambda to decide if edge-node ID is being assigned a new ID by us
  auto ownedge = [ this ]( const tk::UnsMesh::Edge& e ) {
    using Set = typename std::remove_reference<
                  decltype(m_ecommunication) >::type::value_type;
    return !std::any_of( m_ecommunication.cbegin(), m_ecommunication.cend(),
                         [&](const Set& s)
                         { return s.second.find(e) != s.second.cend(); } );
  };

  // Reorder our chunk of the mesh node IDs by looping through all of our
  // node IDs (resulting from reading our chunk of the mesh cells). We test
  // if we are to assign a new ID to a node ID, and if so, we assign new ID,
  // i.e., reorder, by constructing a map associating new to old IDs. We
  // also count up the reordered nodes, which also serves as the new node
  // id.
  for (auto p : m_nodeset)
    if (ownnode(p))
      m_linnodes[ p ] = m_start++;

  // Reorder our chunk of the mesh edges by looping through all of our edges
  // (resulting from initial uniform refinement of our chunk of the mesh
  // cells). We test if we are to assign a new ID to an edge, and if so, we
  // assign new ID, i.e., reorder, by constructing a map associating new
  // node IDs to edges. We also count up the reordered edge-nodes, which
  // also serves as the new node id.
  for (const auto& e : m_edgeset)
    if (ownedge(e)) m_linedges[ e ] = m_start++;

  // Trigger SDAG wait indicating that reordering own node IDs are complete
  reorderowned_complete();

  // If all our nodes have new IDs assigned, signal that to the runtime
  if (m_linnodes.size() == m_nodeset.size()) nodesreorder_complete();

  // If all our edges have new IDs assigned, signal that to the runtime
  if (m_linedges.size() == m_edgeset.size()) edgesreorder_complete();
}

int
Partitioner::pe( int id ) const
// *****************************************************************************
//  Return processing element for chare id
//! \param[in] id Chare id
//! \return PE that creates the chare
//! \details This is computed based on a simple contiguous linear
//!   distribution of chare ids to PEs.
// *****************************************************************************
{
  Assert( m_nchare > 0, "Number of chares must be a positive number" );
  auto p = id / (m_nchare / CkNumPes());
  if (p >= CkNumPes()) p = CkNumPes()-1;
  Assert( p < CkNumPes(), "Assigning to nonexistent PE" );
  return p;
}

void
Partitioner::prepare()
// *****************************************************************************
//  Associate new node IDs to old ones and return them to the requestor(s)
// *****************************************************************************
{
  // Signal to the runtime system that we have participated in reordering
  participated_complete();

  // Find and return new node IDs to sender
  for (const auto& r : m_reqNodes) {
    std::unordered_map< std::size_t, std::size_t > n;
    for (auto p : r.second) n[ p ] = tk::cref_find( m_linnodes, p );
    thisProxy[ r.first ].neworder( n );
    tk::destroy( n );
  }

  tk::destroy( m_reqNodes ); // Clear queue of requests just fulfilled

  // Find and return new node IDs associated to edges to sender
  for (const auto& r : m_reqEdges) {
    tk::UnsMesh::EdgeNodes n;
    for (const auto& e : r.second) n[ e ] = tk::cref_find( m_linedges, e );
    thisProxy[ r.first ].neworder( n );
    tk::destroy( n );
  }

  tk::destroy( m_reqEdges ); // Clear queue of requests just fulfilled

  wait4prep();      // Re-enable SDAG wait for preparing new node requests

  // Re-enable trigger signaling that reordering of owned node IDs are
  // complete right away
  reorderowned_complete();
}

void
Partitioner::generate_compact_inpoel()
// *****************************************************************************
//  Generate compact mesh connectivity
// *****************************************************************************
{
  // Concatenate mesh connectivities of our chares
  tk::destroy(m_tetinpoel);

  std::size_t nn = 0;

  // M_chinpoel contains an array of size chars.
  // Each slot in the array contains "connectivity"
  // This is a vector of size_t. This is tets?

  // This loop counts up the total number of elements in all lsots of
  // m_chinpoel
  for (const auto &c : m_chinpoel)
  {
      nn += c.second.size();
  }

  // Resize global array to avoid resize later
  m_tetinpoel.resize(nn);

  // Flatten out m_chinpoel into m_tetinpoel
  size_t n = 0;
  for (const auto &c : m_chinpoel) {
      for (auto i : c.second) {
          m_tetinpoel[n++] = i;
      }
  }
}

// void
// Partitioner::refine()
// // *****************************************************************************
// //  Uniformly refine our mesh replacing each tetrahedron with 8 new ones
// // *****************************************************************************
// {
//   generate_compact_inpoel();
// 
//   // Create AMR object
//   AMR::mesh_adapter_t* mesh_adapter = new AMR::mesh_adapter_t();
// 
//   // Generate unique edges (nodes connected to nodes)?
// 
//   // shift to zero-based node IDs
//   // Find min/max to shift
//   auto minmax = std::minmax_element(begin(m_tetinpoel), end(m_tetinpoel));
// 
//   std::array<std::size_t, 2> ext{{*minmax.first, *minmax.second}};
//   auto nnode = ext[1] - ext[0] + 1;
// 
//   mesh_adapter->init(m_tetinpoel, nnode);
// 
//   // Do uniform refinement
//   mesh_adapter->uniform_refinement();
// 
//   // TODO: Do we need to replicate the shift?
// 
//   /*
//   // Perform shift
//   for (auto &i : m_tetinpoel) {
//   i -= ext[0];
//   }
// 
//   // TODO: Is there a reason we want to not enforce 0 based no ids
//   // always? (i.e not put them back?)
// 
//   // Generate elements surrounding points..? If we have tets, why do we
//   // need to do this? Can't we just generate edges?
//   auto esup = tk::genEsup(m_tetinpoel, 4);
// 
//   for (auto &i : m_tetinpoel) i += ext[0];  // shift back node IDs
// 
//   std::unordered_map<std::size_t, std::unordered_set<std::size_t> > star;
// 
//   for (std::size_t j = 0; j < nnode; ++j) {
//   for (std::size_t i = esup.second[j] + 1; i <= esup.second[j + 1]; ++i) {
//   for (std::size_t n = 0; n < 4; ++n) {
//   auto p = ext[0] + j;
//   auto q = m_tetinpoel[esup.first[i] * 4 + n];
//   if (p < q) star[p].insert(q);
//   if (p > q) star[q].insert(p);
//   }
//   }
//   }
// 
//   // Starting node ID (on all PEs) while assigning new edge-nodes
//   nnode = tk::ExodusIIMeshReader(g_inputdeck.get<tag::cmd, tag::io,
//   tag::input>()).readHeader();
// 
//   // Add new edge-nodes
//   tk::UnsMesh::EdgeNodes edgenodes;
//   for (const auto &s : star) {
//   for (auto q : s.second) {
//   edgenodes[ {{ s.first, q }} ] = nnode++;
//   }
//   }
//   */
// 
//   tk::destroy(m_tetinpoel);
// 
//   // Generate maps associating new node IDs (as in producing
//   // contiguous-row-id linear system contributions)to edges (a pair of old
//   // node IDs) in tk::UnsMesh::EdgeNodes maps, associated to and categorized
//   // by chares. Note that the new edge-node IDs assigned here will be
//   // overwritten with globally unique node IDs after reordering.
//   for (const auto& conn : m_chinpoel) {
//       auto& en = m_chedgenodes[ conn.first ];
//       for (std::size_t e=0; e<conn.second.size()/4; ++e) {
//           const auto A = conn.second[e*4+0];
//           const auto B = conn.second[e*4+1];
//           const auto C = conn.second[e*4+2];
//           const auto D = conn.second[e*4+3];
//           // Look up the added node IDs based on old ids {A,B}
//           // (vector)
// 
//           /*
//              const auto AB = tk::cref_find( edgenodes, {{ A,B }} );
//              const auto AC = tk::cref_find( edgenodes, {{ A,C }} );
//              const auto AD = tk::cref_find( edgenodes, {{ A,D }} );
//              const auto BC = tk::cref_find( edgenodes, {{ B,C }} );
//              const auto BD = tk::cref_find( edgenodes, {{ B,D }} );
//              const auto CD = tk::cref_find( edgenodes, {{ C,D }} );
//              */
// 
//           // TODO: We should likely check the return values here
//           const int AB = mesh_adapter->node_connectivity.find(A, B);
//           const int AC = mesh_adapter->node_connectivity.find(A, C);
//           const int AD = mesh_adapter->node_connectivity.find(A, D);
//           const int BC = mesh_adapter->node_connectivity.find(B, C);
//           const int BD = mesh_adapter->node_connectivity.find(B, D);
//           const int CD = mesh_adapter->node_connectivity.find(C, D);
// 
//           en[ {{A,B}} ] = static_cast<size_t>(AB);
//           en[ {{A,C}} ] = static_cast<size_t>(AC);
//           en[ {{A,D}} ] = static_cast<size_t>(AD);
//           en[ {{B,C}} ] = static_cast<size_t>(BC);
//           en[ {{B,D}} ] = static_cast<size_t>(BD);
//           en[ {{C,D}} ] = static_cast<size_t>(CD);
//       }
//   }
// 
//   generate_compact_inpoel();
// 
//   delete mesh_adapter;
// 
//   // TODO: This only needs to have set en? Which is m_chedgenodes.
// }


void
Partitioner::refine( const std::vector< std::size_t >& inpoel )
// *****************************************************************************
// Optionally refine mesh if requested by user
// *****************************************************************************
{
  const auto ir = g_inputdeck.get< tag::amr, tag::init >();
  if (!ir.empty()) {

    auto orig_inpoel = inpoel;
    auto orig_coord = m_coord;

    AMR::mesh_adapter_t refiner( orig_inpoel );

    for (std::size_t level=0; level<2; ++level) {

      refiner.uniform_refinement();
      auto refined_inpoel = refiner.tet_store.get_active_inpoel();

//       std::cout << CkMyPe() << ":c: ";
//       for (auto p : orig_inpoel) std::cout << p << ' ';
//       std::cout << '\n';
//       std::cout << CkMyPe() << ":rc: ";
//       for (auto p : refined_inpoel) std::cout << p << ' ';
//       std::cout << '\n';

      // find number of nodes in old mesh
      auto minmax = std::minmax_element( begin(orig_inpoel), end(orig_inpoel) );
      Assert( *minmax.first == 0, "node ids should start from zero" );
      auto npoin = *minmax.second + 1;

      // generate edges surrounding points in old mesh
      auto esup = tk::genEsup( orig_inpoel, 4 );
      auto psup = tk::genPsup( orig_inpoel, 4, esup );

      auto refined_coord = orig_coord;
      auto& x = refined_coord[0];
      auto& y = refined_coord[1];
      auto& z = refined_coord[2];

      auto refined_minmax =
        std::minmax_element( begin(refined_inpoel), end(refined_inpoel) );
      Assert( *refined_minmax.first == 0, "node ids should start from zero" );
      auto refined_npoin = *refined_minmax.second + 1;

      x.resize( refined_npoin );
      y.resize( refined_npoin );
      z.resize( refined_npoin );

      // generate coordinates for newly added nodes
      //std::cout << CkMyPe() << ": npoin:" << npoin << " -> " << refined_npoin << ": ";
      for (std::size_t p=0; p<npoin; ++p)
        for (auto q : tk::Around(psup,p)) {
          auto e = refiner.node_connectivity.find( p, q );
          if (e > 0) {
            auto E = static_cast< std::size_t >( e );
            Assert( E < x.size(), "Indexing out of refined_coord" );
            x[E] = (x[p]+x[q])/2.0;
            y[E] = (y[p]+y[q])/2.0;
            z[E] = (z[p]+z[q])/2.0;
          } else std::cout << level << ": " << e << ':' << p << "." << q << ' ';
        }
      //std::cout << '\n';

//       // reorder new mesh
//       const auto refined_psup =
//         tk::genPsup( refined_inpoel, 4, tk::genEsup( refined_inpoel, 4 ) );
//       auto map = tk::renumber( refined_psup );
//       tk::remap( refined_inpoel, map );
//       tk::remap( x, map );
//       tk::remap( y, map );
//       tk::remap( z, map );

      tk::UnsMesh refmesh( refined_inpoel, x, y, z );
      tk::ExodusIIMeshWriter mr( "refmesh." + std::to_string(level) + ".exo",
                                 tk::ExoWriter::CREATE );
      mr.writeMesh( refmesh );

      orig_inpoel = refined_inpoel;
      orig_coord = refined_coord;

    } // level
  } // if enable uniform

  // send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.perefined();

  contribute( m_cb.get< tag::refined >() );
}


void
Partitioner::reordered()
// *****************************************************************************
//  Compute final result of reordering
//! \details This member function is called when both those node IDs that we
//!   assign a new ordering to as well as those assigned new IDs by other
//!   PEs have been reordered (and we contribute to) and we are ready (on
//!   this PE) to compute our final result of the reordering.
// *****************************************************************************
{
  // Free memory used by communication maps used to store nodes and
  // edge-nodes and associated PEs during reordering.
  tk::destroy( m_ncommunication );
  tk::destroy( m_ecommunication );

  // Free memory used by maps associating a list of chare IDs to old (as in
  // file) global mesh node IDs and to edges as no longer needed.
  tk::destroy( m_bnodechares );
  tk::destroy( m_edgechares );

  // Construct maps associating old node IDs (as in file) to new node IDs
  // (as in producing contiguous-row-id linear system contributions)
  // associated to chare IDs (outer key).
  for (const auto& c : m_chinpoel) {
    auto& nodes = m_chfilenodes[ c.first ];
    for (auto p : c.second)
      nodes[ tk::cref_find(m_linnodes,p) ] = p;
  }

  // Update node IDs of edges, i.e., the map values
  for (auto& c : m_chedgenodes)
    for (auto& e : c.second)
       e.second = tk::ref_find( m_linedges, e.first );

  if (!m_chedgenodes.empty()) {

    // Update chare-categorized element connectivities with new nodes and
    // newly added edge-nodes during initial unifor mesh refinement
    decltype(m_chinpoel) refinpoel;
    for (const auto& chi : m_chinpoel) {
      auto& ri = refinpoel[ chi.first ];
      const auto& edgenodes = tk::cref_find( m_chedgenodes, chi.first );
      for (std::size_t e=0; e<chi.second.size()/4; ++e) {
        auto A = chi.second[e*4+0];
        auto B = chi.second[e*4+1];
        auto C = chi.second[e*4+2];
        auto D = chi.second[e*4+3];
        const auto nA = tk::cref_find( m_linnodes, A );
        const auto nB = tk::cref_find( m_linnodes, B );
        const auto nC = tk::cref_find( m_linnodes, C );
        const auto nD = tk::cref_find( m_linnodes, D );
        const auto AB = tk::cref_find( edgenodes, {{ A,B }} );
        const auto AC = tk::cref_find( edgenodes, {{ A,C }} );
        const auto AD = tk::cref_find( edgenodes, {{ A,D }} );
        const auto BC = tk::cref_find( edgenodes, {{ B,C }} );
        const auto BD = tk::cref_find( edgenodes, {{ B,D }} );
        const auto CD = tk::cref_find( edgenodes, {{ C,D }} );
        // update connectivity of our mesh chunk
        std::vector< std::size_t > newelems{{ nA, AB, AC, AD,
                                              nB, BC, AB, BD,
                                              nC, AC, BC, CD,
                                              nD, AD, CD, BD,
                                              BC, CD, AC, BD,
                                              AB, BD, AC, AD,
                                              AB, BC, AC, BD,
                                              AC, BD, CD, AD }};
        ri.insert( end(ri), begin(newelems), end(newelems) );
      }
    }
    m_chinpoel = std::move( refinpoel );

    // Update chare-categorized mesh nodes surrounding our mesh chunk with
    // the reordered node IDs
    decltype(m_msum) newmsum;
    for (const auto& c : m_msum) {
      auto& m = newmsum[ c.first ];
      for (const auto& e : c.second) {
        auto& s = m[ e.first ];
        for (auto n : e.second)
          s.insert( tk::cref_find( m_linnodes, n ) );
      }
    }
    m_msum = std::move( newmsum );

    // Add newly added edge-nodes to chare-categorized mesh nodes
    // surrounding our mesh chunk
    for (const auto& c : m_msumed) {
      auto& sur = tk::cref_find( m_msum, c.first );
      const auto& edgenodes = tk::cref_find( m_chedgenodes, c.first );
      for (const auto& e : c.second) {
        auto& s = tk::ref_find( sur, e.first );
        for (const auto& ed : e.second)
          s.insert( tk::cref_find( edgenodes, ed ) );
      }
    }
    tk::destroy( m_msumed );

    // Update node IDs of edges, i.e., the map keys
    for (auto& c : m_chedgenodes) {
      tk::UnsMesh::EdgeNodes edgenodes;
      for (auto& e : c.second)
         edgenodes[ {{ tk::ref_find(m_linnodes,e.first[0]),
                       tk::ref_find(m_linnodes,e.first[1]) }} ] = e.second;
      c.second = std::move( edgenodes );
    }

  } else {

    // Update chare-categorized elem connectivities with the reordered node IDs
    for (auto& c : m_chinpoel)
      for (auto& p : c.second)
         p = tk::cref_find( m_linnodes, p );

    // Update chare-categorized mesh chare-nodes comm map with the reordered
    // node IDs
    for (auto& c : m_msum)
      for (auto& s : c.second) {
        decltype(s.second) n;
        for (auto p : s.second) {
          n.insert( tk::cref_find( m_linnodes, p ) );
        }
        s.second = std::move( n );
      }

  }

  // Update unique global node IDs chares on our PE will contribute to with
  // the reordered node IDs
  tk::destroy( m_nodeset );
  for (const auto& c : m_chinpoel)
    for (auto i : c.second)
      m_nodeset.insert( i );

  // send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pereordered();

  // Compute lower and upper bounds of reordered node IDs our PE operates on
  bounds();
}

void
Partitioner::bounds()
// *****************************************************************************
// Compute lower and upper bounds of reordered node IDs our PE operates on
// \details This function computes the global row IDs at which the linear
//   system will have a PE boundary. We simply find the largest node ID
//   assigned on each PE by the reordering and use that as the upper global
//   row index. Note that while this rarely results in equal number of rows
//   assigned to PEs, potentially resulting in some load-imbalance, it
//   yields a pretty good division reducing communication costs during the
//   assembly of the linear system, which is more important than a slight
//   (FLOP) load imbalance. Since the upper index for PE 1 is the same as
//   the lower index for PE 2, etc., we find the upper indices and then the
//   lower indices for all PEs are communicated.
// *****************************************************************************
{
  m_upper = 0;

  using P1 = std::pair< const std::size_t, std::size_t >;
  for (const auto& c : m_chfilenodes) {
    auto x = std::max_element( begin(c.second), end(c.second),
             [](const P1& a, const P1& b){ return a.first < b.first; } );
    if (x->first > m_upper) m_upper = x->first;
  }

  using P2 = std::pair< const tk::UnsMesh::Edge, std::size_t >;
  for (const auto& c : m_chedgenodes) {
    auto x = std::max_element( begin(c.second), end(c.second),
             [](const P2& a, const P2& b){ return a.second < b.second; } );
    if (x->second > m_upper) m_upper = x->second;
  }

  // The bounds are the dividers (global mesh point indices) at which the
  // linear system assembly is divided among PEs. However, Hypre and thus
  // Solver expect exclusive upper indices, so we increase the last one by
  // one here. Note that the cost calculation, Partitioner::cost() also
  // expects exclusive upper indices.
  if (CkMyPe() == CkNumPes()-1) ++m_upper;

  // Tell the runtime system that the upper bound has been computed
  upper_complete();

  // Set lower index for PE 0 as 0
  if (CkMyPe() == 0) lower(0);

  // All PEs except the last one send their upper indices as the lower index for
  // PE+1
  if (CkMyPe() < CkNumPes()-1)
    thisProxy[ CkMyPe()+1 ].lower( m_upper );
}

void
Partitioner::create()
// *****************************************************************************
// Create chare array elements on this PE and assign the global mesh element IDs
// they will operate on
// *****************************************************************************
{
  // send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pebounds();

  // Initiate asynchronous reduction across all Partitioner objects computing
  // the average communication cost of merging the linear system
  m_cost = cost( m_lower, m_upper );
  contribute( sizeof(tk::real), &m_cost, CkReduction::sum_double,
              m_cb.get< tag::avecost >() );

  // Create worker chare array elements
  createDiscWorkers();

  // Broadcast our bounds of global node IDs to all matrix solvers
  const auto scheme = g_inputdeck.get< tag::selected, tag::scheme >();
  if (scheme == ctr::SchemeType::MatCG || scheme == ctr::SchemeType::DiagCG)
    m_solver.bounds( CkMyPe(), m_lower, m_upper );
  else // if no MatCG, no matrix solver, continue
    contribute( m_cb.get< tag::coord >() );
}

void
Partitioner::createDiscWorkers()
// *****************************************************************************
//  Create Discretization chare array elements on this PE
//! \details We create chare array elements by calling the insert() member
//!   function, which allows specifying the PE on which the array element is
//!   created. and we send each chare array element the chunk of mesh it will
//!   operate on.
// *****************************************************************************
{
  auto dist = chareDistribution();

  for (int c=0; c<dist[1]; ++c) {
    // Compute chare ID
    auto cid = CkMyPe() * dist[0] + c;
    // Guard those searches that operate on empty containers in serial
    typename decltype(m_msum)::mapped_type msum;
    if (!m_msum.empty()) msum = tk::cref_find( m_msum, cid );
    typename decltype(m_chedgenodes)::mapped_type edno;
    if (!m_chedgenodes.empty()) edno = tk::cref_find( m_chedgenodes, cid );
    // Create worker array element
    m_scheme.discInsert( cid, m_host, m_bc,
      tk::cref_find(m_chinpoel,cid), msum, tk::cref_find(m_chfilenodes,cid),
      edno, m_nchare, CkMyPe() );
  }

  // Free storage for unique global mesh nodes chares on our PE will
  // contribute to in a linear system as no longer needed.
  tk::destroy( m_nodeset );
  // Free storage for unique global mesh edges whose nodes chares on our
  // PE will contribute to in a linear system as no longer needed.
  tk::destroy( m_edgeset );
  // Free map storing new node IDs associated to edges categorized by chares
  // owned as no linger needed after creating workers.
  tk::destroy( m_chedgenodes );
  // Free storage of global mesh node IDs associated to chare IDs bordering
  // the mesh chunk held by and associated to chare IDs we own as it is no
  // longer needed after creating the workers.
  tk::destroy( m_msum );
}

void
Partitioner::createWorkers()
// *****************************************************************************
//  Create worker chare array elements on this PE
//! \details We create chare array elements by calling the insert() member
//!   function, which allows specifying the PE on which the array element is
//!   created. and we send each chare array element the chunk of mesh it will
//!   operate on.
// *****************************************************************************
{
  auto dist = chareDistribution();

  // Prepare data to pass to and put in a request to create worker chares
  for (int c=0; c<dist[1]; ++c) {

    // Compute chare ID (linear distribution across PEs)
    auto cid = CkMyPe() * dist[0] + c;

    // Find mesh connectivity for this chare (with new global ids)
    auto chinpoel = tk::cref_find( m_chinpoel, cid );

    // Generate map associating new(value) to file(key) node ids for this chare
    decltype(m_linnodes) newnodes;
    const auto& chfilenodes = tk::cref_find( m_chfilenodes, cid );
    for (auto i : chinpoel) newnodes[ tk::cref_find(chfilenodes,i) ] = i;

    // Generate face set for this chare for potentially faster searches
    tk::UnsMesh::FaceSet faceset;
    for (std::size_t e=0; e<chinpoel.size()/4; ++e) { // for all tets in chare
      auto mark = e*4;
      for (std::size_t f=0; f<4; ++f) // for all tet faces
        faceset.insert( {{{ chinpoel[ mark + tk::lpofa[f][0] ],
                            chinpoel[ mark + tk::lpofa[f][1] ],
                            chinpoel[ mark + tk::lpofa[f][2] ] }}} );
    }
  
    // Generate input face data for class FaceData
    std::vector< std::size_t > chtriinpoel;
    std::unordered_map< int, std::vector< std::size_t > > chbface;
    std::size_t cnt = 0;

    for (const auto& ss : m_bface)  // for all phsyical boundaries (sidesets)
      for (auto f : ss.second) {    // for all faces on this physical boundary
        auto f1 = newnodes.find( m_triinpoel[f*3+0] );
        auto f2 = newnodes.find( m_triinpoel[f*3+1] );
        auto f3 = newnodes.find( m_triinpoel[f*3+2] );
        // if all 3 nodes of the physical boundary face are on this chare
        if (f1 != end(newnodes) && f2 != end(newnodes) && f3 != end(newnodes)) {
          std::array< std::size_t, 3 > t{{f1->second, f2->second, f3->second}};
          // if this boundary face is on this chare
          if (faceset.find(t) != end(faceset)) {
            // store face connectivity with new (global) node ids of this chare
            for (std::size_t i=0; i<3; ++i) chtriinpoel.push_back( t[i] );
            // store physical boundary face id associated to sideset id
            chbface[ ss.first ].push_back( cnt++ );
          }
        }
      }

    // Face data class
    FaceData fd( chinpoel, chbface, chtriinpoel );

    // Make sure (bound) base is already created and accessible
    Assert( m_scheme.get()[cid].ckLocal() != nullptr, "About to pass nullptr" );

    // Create worker array element
    m_scheme.insert( cid, m_scheme.get(), m_solver, fd, CkMyPe() );
  }

  tk::destroy( m_bface );
  tk::destroy( m_triinpoel );
}

#include "NoWarning/partitioner.def.h"
