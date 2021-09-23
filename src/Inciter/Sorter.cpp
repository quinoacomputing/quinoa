// *****************************************************************************
/*!
  \file      src/Inciter/Sorter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh sorter for global distributed mesh reordering
  \see       Sorter.h for more info.
*/
// *****************************************************************************

#include <vector>
#include <algorithm>

#include "Sorter.hpp"
#include "Reorder.hpp"
#include "DerivedData.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::Sorter;

Sorter::Sorter( std::size_t meshid,
                const CProxy_Transporter& transporter,
                const tk::CProxy_MeshWriter& meshwriter,
                const tk::SorterCallback& cbs,
                const std::vector< Scheme >& scheme,
                CkCallback reorderRefiner,
                const std::vector< std::size_t >& ginpoel,
                const tk::UnsMesh::CoordMap& coordmap,
                const tk::UnsMesh::Chunk& el,
                const std::map< int, std::vector< std::size_t > >& bface,
                const std::vector< std::size_t >& triinpoel,
                const std::map< int, std::vector< std::size_t > >& bnode,
                int nchare ) :
  m_meshid( meshid ),
  m_host( transporter ),
  m_meshwriter( meshwriter ),
  m_cbs( cbs ),
  m_scheme( scheme ),
  m_reorderRefiner( reorderRefiner ),
  m_ginpoel( ginpoel ),
  m_coordmap( coordmap ),
  m_el( el ),
  m_nbnd( 0 ),
  m_bface( bface ),
  m_triinpoel( triinpoel ),
  m_bnode( bnode ),
  m_nchare( nchare ),
  m_nodeset( begin(ginpoel), end(ginpoel) ),
  m_noffset( 0 ),
  m_nodech(),
  m_chnode(),
  m_edgech(),
  m_chedge(),
  m_msum(),
  m_reordcomm(),
  m_start( 0 ),
  m_newnodes(),
  m_newcoordmap(),
  m_reqnodes(),
  m_lower( 0 ),
  m_upper( 0 )
// *****************************************************************************
//  Constructor: prepare owned mesh node IDs for reordering
//! \param[in] meshid Mesh ID
//! \param[in] transporter Transporter (host) Charm++ proxy
//! \param[in] meshwriter Mesh writer Charm++ proxy
//! \param[in] cbs Charm++ callbacks for Sorter
//! \param[in] scheme Discretization schemes (one per mesh)
//! \param[in] reorderRefiner Callback to use to send reordered mesh to Refiner
//! \param[in] ginpoel Mesh connectivity (this chare) using global node IDs
//! \param[in] coordmap Mesh node coordinates (this chare) for global node IDs
//! \param[in] bface Face lists mapped to side set ids
//! \param[in] triinpoel Interconnectivity of points and boundary-faces
//! \param[in] bnode Node ids mapped to side set ids
//! \param[in] nchare Total number of Charm++ worker chares
// *****************************************************************************
{
  // Ensure boundary face ids will not index out of face connectivity
  Assert( std::all_of( begin(m_bface), end(m_bface),
            [&](const auto& s)
            { return std::all_of( begin(s.second), end(s.second),
                       [&](auto f){ return f*3+2 < m_triinpoel.size(); } ); } ),
          "Boundary face data structures inconsistent" );
}

void
Sorter::setup( std::size_t npoin )
// *****************************************************************************
// Setup chare mesh boundary node communication map
//! \param[in] npoin Total number of mesh points in mesh. Note that the number
//!   of mesh points does not have to be exactly the total number of points in
//!   the mesh. It can be a larger number, but not less. This is only used here
//!   to assign nodes to workers that will assign ids to mesh nodes during node
//!   reordering.
// *****************************************************************************
{
  // Compute the number of nodes (chunksize) a chare will build a node
  // communication map for. We compute two values of chunksize: one for when
  // the global node ids are abounded between [0...npoin-1], inclusive, and
  // another one for when the global node ids are assigned by a hash algorithm
  // during initial mesh refinement. In the latter case, the maximum
  // representable value of a std::size_t is assumed to be the large global node
  // id and is used to compute the chunksize. To compute the bin id, we attempt
  // to use the first chunksize first: if it gives a chare id that is
  // (strictly) lower than the number of chares, that's good. If not, we compute
  // the bin id based on the second chunksize, which almost always will give a
  // bin id strictly lower than the number of chares, except if the global node
  // id assigned by the hash algorithm in Refiner hits the maximum
  // representable number in std::size_t. If that is the case, we just assign
  // that node to the last chare.
  auto N = static_cast< std::size_t >( m_nchare );
  std::array< std::size_t, 2 > chunksize{{
     npoin / N, std::numeric_limits< std::size_t >::max() / N }};

  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();

  // Find chare-boundary nodes and edges of our mesh chunk. This algorithm
  // collects the global mesh node ids and edges on the chare boundary. A node
  // is on a chare boundary if it belongs to a face of a tetrahedron that has
  // no neighbor tet at a face. The edge is on the chare boundary if its first
  // edge-end point is on a chare boundary. The nodes are categorized to bins
  // that will be sent to different chares to build point-to-point
  // communication maps across all chares. The binning is determined by the
  // global node id divided by the chunksizes. See discussion above on how we
  // use two chunksizes for global node ids assigned by the hash algorithm in
  // Refiner (if initial mesh refinement has been done).
  tk::CommMaps chbnd;
  auto el = tk::global2local( m_ginpoel );      // generate local mesh data
  const auto& inpoel = std::get< 0 >( el );     // local connectivity
  auto esup = tk::genEsup( inpoel, 4 );         // elements surrounding points
  auto esuel = tk::genEsuelTet( inpoel, esup ); // elems surrounding elements
  for (std::size_t e=0; e<esuel.size()/4; ++e) {
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f)
      if (esuel[mark+f] == -1)
        for (std::size_t n=0; n<3; ++n) {
          auto g = m_ginpoel[ mark+tk::lpofa[f][n] ];
          auto bin = g / chunksize[0];
          if (bin >= N) bin = g / chunksize[1];
          if (bin >= N) bin = N - 1;
          Assert( bin < N, "Will index out of number of chares" );
          auto& b = chbnd[ static_cast< int >( bin ) ];
          b.get< tag::node >().insert( g );
          if (scheme == ctr::SchemeType::ALECG) {
            auto h = m_ginpoel[ mark + tk::lpofa[ f ][ tk::lpoet[n][1] ] ];
            b.get< tag::edge >().insert( { std::min(g,h), std::max(g,h) } );
          }
        }
  }

  // Send boundary data in bins to chares that will compute communication maps
  // for the data in the bin. These bins form a distributed table.  Note that
  // we only send data to those chares that have data to work on. The receiving
  // sides do not know in advance if they receive messages or not.  Completion
  // is detected by having the receiver respond back and counting the responses
  // on the sender side, i.e., this chare.
  m_nbnd = chbnd.size();
  if (m_nbnd == 0)
    contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
                m_cbs.get< tag::queried >() );
  else
    for (const auto& [ targetchare, bnd ] : chbnd)
      thisProxy[ targetchare ].query( thisIndex, bnd );
}

void
Sorter::query( int fromch, const tk::AllCommMaps& bnd )
// *****************************************************************************
// Incoming query for a list of mesh nodes for which this chare compiles node
// communication maps
//! \param[in] fromch Sender chare ID
//! \param[in] bnd Chare-boundary data from another chare
// *****************************************************************************
{
  // Store incoming nodes in node->chare and its inverse, chare->node, maps
  const auto& nodes = bnd.get< tag::node >();
  for (auto n : nodes) m_nodech[ n ].push_back( fromch );
  m_chnode[ fromch ].insert( begin(nodes), end(nodes) );

  // Store incoming edges in edge->chare and its inverse, chare->edge, maps
  const auto& edges = bnd.get< tag::edge >();
  for (const auto& e : edges) m_edgech[ e ].push_back( fromch );
  m_chedge[ fromch ].insert( begin(edges), end(edges) );

  // Report back to chare message received from
  thisProxy[ fromch ].recvquery();
}

void
Sorter::recvquery()
// *****************************************************************************
// Receive receipt of boundary node lists to query
// *****************************************************************************
{
  if (--m_nbnd == 0)
    contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
                m_cbs.get< tag::queried >() );
}

void
Sorter::response()
// *****************************************************************************
//  Respond to boundary node list queries
// *****************************************************************************
{
  std::unordered_map< int, tk::CommMaps > exp;

  // Compute node communication map to be sent back to chares
  for (const auto& [ neighborchare, bndnodes ] : m_chnode) {
    auto& nc = exp[ neighborchare ];
    for (auto n : bndnodes)
      for (auto d : tk::cref_find(m_nodech,n))
        if (d != neighborchare)
          nc[d].get< tag::node >().insert( n );
  }

  // Compute edge communication map to be sent back to chares
  for (const auto& [ neighborchare, bndedges ] : m_chedge) {
    auto& ec = exp[ neighborchare ];
    for (const auto& e : bndedges)
      for (auto d : tk::cref_find(m_edgech,e))
        if (d != neighborchare)
          ec[d].get< tag::edge >().insert( e );
  }

  // Send communication maps to chares that issued a query to us. Communication
  // maps were computed above for those chares that queried this map from us.
  // This data form a distributed table and we only work on a chunk of it. Note
  // that we only send data back to those chares that have queried us. The
  // receiving sides do not know in advance if the receive messages or not.
  // Completion is detected by having the receiver respond back and counting
  // the responses on the sender side, i.e., this chare.
  m_nbnd = exp.size();
  if (m_nbnd == 0)
    contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
                m_cbs.get< tag::responded >() );
  else
    for (const auto& [ targetchare, maps ] : exp)
      thisProxy[ targetchare ].bnd( thisIndex, maps );
}

void
Sorter::bnd( int fromch, const tk::CommMaps& msum )
// *****************************************************************************
// Receive boundary node communication maps for our mesh chunk
//! \param[in] fromch Sender chare ID
//! \param[in] msum Communication map(s) assembled by chare fromch
// *****************************************************************************
{
  for (const auto& [ neighborchare, maps ] : msum) {
    auto& m = m_msum[ neighborchare ];
    const auto& nodemap = maps.get< tag::node >();
    m.get< tag::node >().insert( begin(nodemap), end(nodemap) );
    const auto& edgemap = maps.get< tag::edge >();
    m.get< tag::edge >().insert( begin(edgemap), end(edgemap) );
  }

  // Report back to chare message received from
  thisProxy[ fromch ].recvbnd();
}

void
Sorter::recvbnd()
// *****************************************************************************
// Receive receipt of boundary node communication map
// *****************************************************************************
{
  if (--m_nbnd == 0)
    contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
                m_cbs.get< tag::responded >() );
}

void
Sorter::start()
// *****************************************************************************
//  Start reordering (if enabled)
// *****************************************************************************
{
  // Keep only those edges in edge comm map whose both end-points are in the
  // node comm map
  for (auto& [ neighborchare, maps ] : m_msum) {
    const auto& nodes = maps.get< tag::node >();
    tk::EdgeSet edges;
    for (const auto& e : maps.get< tag::edge >())
      if (nodes.find(e[0]) != end(nodes) && nodes.find(e[1]) != end(nodes))
        edges.insert( e );
    maps.get< tag::edge >() = std::move(edges);
  }

  if (g_inputdeck.get< tag::cmd, tag::feedback >()) m_host.chcomm();

  tk::destroy( m_nodech );
  tk::destroy( m_chnode );

  if (g_inputdeck.get< tag::discr, tag::pelocal_reorder >())
    mask();   // continue with mesh node reordering if requested (or required)
  else
    createDiscWorkers();  // skip mesh node reordering
}

void
Sorter::mask()
// *****************************************************************************
//  Start preparing for mesh node reordering in parallel
// *****************************************************************************
{
  // Compute asymmetric communcation map that will be used for reordering. This
  // communication map is asymmetric because it associates global mesh node IDs
  // to chares only with lower IDs than thisIndex. That is because this chare
  // will need to receive new (reorderd) node IDs only from chares with lower
  // IDs than thisIndex during node reordering. Since it only stores data for
  // lower chare IDs, it is asymmetric. Note that because of this algorithm the
  // type of m_msum is an ordered map, because of the std::none_of() algorithm
  // needs to look at ALL chares this chare potentially communicates nodes with
  // that have lower chare IDs that thisIndex. Since the map is ordered, it can
  // walk through from the beginning of m_msum until the outer loop variable c,
  // which is the chare ID the outer loop works on in a given cycle.
  for (auto c=m_msum.cbegin(); c!=m_msum.cend(); ++c)
    if (thisIndex > c->first) {
      auto& n = m_reordcomm[ c->first ];
      for (auto j : c->second.get< tag::node >())
        if (std::none_of( m_msum.cbegin(), c,
             [j]( const auto& s ) {
               const auto& nodemap = s.second.template get< tag::node >();
               return nodemap.find(j) != end(nodemap); } ))
        {
          n.insert(j);
        }
      if (n.empty()) m_reordcomm.erase( c->first );
    }

  // Count up total number of nodes this chare will need to receive
  auto nrecv = tk::sumvalsize( m_reordcomm );

  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.chmask();

  // Compute number of mesh node IDs we will assign IDs to
  auto nuniq = m_nodeset.size() - nrecv;

  // Start computing offsets for node reordering
  thisProxy.offset( thisIndex, nuniq );
}

void
Sorter::offset( int c, std::size_t u )
// *****************************************************************************
//  Receive number of uniquely assigned global mesh node IDs from chares with
//  lower IDs than thisIndex
//! \param[in] c Chare ID
//! \param[in] u Number of mesh node IDs chare c will assign IDs to
//! \details This function computes the offset each chare will need to start
//!   assigning its new node IDs from. The offset for a chare is the
//!   offset for the previous chare plus the number of node IDs the previous
//!   chare (uniquely) assigns new IDs for minus the number of node IDs the
//!   previous chare receives from others (lower chares). This is computed here
//!   in a parallel/distributed fashion by each chare sending its number of node
//!   IDs (that it uniquely assigns) to all chares. Note that each chare would
//!   only need to send this information to chares with higher IDs, but instead
//!   this function is called in a broadcast fashion, because that is more
//!   efficient than individual calls to only chares with higher IDs. Therefore
//!   when computing the offsets, we only count the lower chares. When this is
//!   done, we have the precise asymmetric communication map as well as the
//!   start offset on all chares and so we can start the distributed global mesh
//!   node ID reordering.
// *****************************************************************************
{
  if (c < thisIndex) m_start += u;
  if (++m_noffset == m_nchare) reorder();
}

void
Sorter::reorder()
// *****************************************************************************
//  Reorder global mesh node IDs
// *****************************************************************************
{
  // Activate SDAG waits for arriving requests from other chares requesting new
  // node IDs for node IDs we assign new IDs to during reordering; and for
  // computing/receiving lower and upper bounds of global node IDs our chare's
  // linear system will operate on after reordering.
  thisProxy[ thisIndex ].wait4prep();

  // Send out request for new global node IDs for nodes we do not reorder
  for (const auto& [ targetchare, nodes ] : m_reordcomm)
    thisProxy[ targetchare ].request( thisIndex, nodes );

  // Lambda to decide if node is assigned a new ID by this chare. If node is not
  // found in the asymmetric communication map, it is owned, i.e., this chare
  // assigns its new id.
  auto ownnode = [ this ]( std::size_t p ) {
    return std::all_of( m_reordcomm.cbegin(), m_reordcomm.cend(),
                        [&](const auto& s)
                        { return s.second.find(p) == s.second.cend(); } );
  };

  // Reorder our chunk of the mesh node IDs. Looping through all of our node
  // IDs, we test if we are to assign a new ID to a node ID, and if so, we
  // assign a new ID, i.e., reorder, by constructing a map associating new to
  // old IDs (m_newnodes). We also count up the reordered nodes, which serves as
  // the new node id. We also store the node coordinates associated to the new
  // node ID.
  for (auto p : m_nodeset)
    if (ownnode(p)) {
      m_newnodes[ p ] = m_start;        // assign new node ID (reorder)
      m_newcoordmap.emplace( m_start, tk::cref_find(m_coordmap,p) );
      ++m_start;
    }

  // Trigger SDAG wait indicating that reordering our node IDs are complete
  reorderowned_complete();

  // If all our nodes have new IDs assigned, reordering complete on this chare
  if (m_newnodes.size() == m_nodeset.size()) finish();
}

void
Sorter::request( int c, const std::unordered_set< std::size_t >& nd )
// *****************************************************************************
//  Request new global node IDs for old node IDs
//! \param[in] c Chare request coming from and to which we send new IDs to
//! \param[in] nd Set of old node IDs whose new IDs are requested
// *****************************************************************************
{
  // Queue up requesting chare and node IDs
  m_reqnodes.push_back( { c, nd } );
  // Trigger SDAG wait signaling that node IDs have been requested from us
  nodes_requested_complete();
}

void
Sorter::prepare()
// *****************************************************************************
//  Find new node IDs for old ones and return them to the requestor(s)
// *****************************************************************************
{
  // Find and return new node IDs to sender
  for (const auto& [ requestorchare, nodes ] : m_reqnodes) {
    std::unordered_map< std::size_t,
      std::tuple< std::size_t, tk::UnsMesh::Coord > > n;
    for (auto p : nodes) {
      auto newid = tk::cref_find( m_newnodes, p );
      n.emplace( p,
        std::make_tuple( newid, tk::cref_find(m_newcoordmap,newid) ) );
    }
    thisProxy[ requestorchare ].neworder( n );
  }

  tk::destroy( m_reqnodes ); // Clear queue of requests just fulfilled

  // Re-enable SDAG wait for preparing new node requests
  thisProxy[ thisIndex ].wait4prep();

  // Re-enable trigger signaling that reordering of owned node IDs are
  // complete right away
  reorderowned_complete();
}

void
Sorter::neworder( const std::unordered_map< std::size_t,
                        std::tuple< std::size_t, tk::UnsMesh::Coord > >& nodes )
// *****************************************************************************
//  Receive new (reordered) global node IDs
//! \param[in] nodes Map associating new to old node IDs
// *****************************************************************************
{
  // Store new node IDs associated to old ones, and node coordinates associated
  // to new node IDs.
  for (const auto& [ oldid, newnodes ] : nodes) {
    auto newid = std::get< 0 >( newnodes );
    m_newnodes[ oldid ] = newid;
    m_newcoordmap.emplace( newid, std::get< 1 >( newnodes ) );
  }

  // If all our nodes have new IDs assigned, reorder complete on this PE
  if (m_newnodes.size() == m_nodeset.size()) finish();
}

void
Sorter::finish()
// *****************************************************************************
//  Compute final result of reordering
//! \details Reordering is now complete on this chare. We now remap all mesh
//!   data to reflect the new ordering.
// *****************************************************************************
{
  // Update elem connectivity with the reordered node IDs
  tk::remap( m_ginpoel, m_newnodes );

  // Update node coordinate map with the reordered IDs
  m_coordmap = m_newcoordmap;

  // Update mesh chunk data structure held in our state with new node order
  m_el = tk::global2local( m_ginpoel );

  // Update symmetric chare-node communication map with the reordered IDs
  for (auto& [ neighborchare, maps ] : m_msum) {

    tk::NodeSet n;
    for (auto p : maps.get< tag::node >())
      n.insert( tk::cref_find( m_newnodes, p ) );
    maps.get< tag::node >() = std::move( n );

    tk::EdgeSet e;
    for (const auto& ed : maps.get< tag::edge >()) {
      e.insert( { tk::cref_find(m_newnodes,ed[0]),
                  tk::cref_find(m_newnodes,ed[1]) } );
    }
    maps.get< tag::edge >() = std::move( e );

  }

  // Update boundary face-node connectivity with the reordered node IDs
  tk::remap( m_triinpoel, m_newnodes );

  // Update boundary node lists with the reordered node IDs
  for (auto& [ setid, nodes ] : m_bnode) tk::remap( nodes, m_newnodes );

  // Update mesh in Refiner after reordering
  m_reorderRefiner.send();

  // Progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.chreordered();

  createDiscWorkers();
}

void
Sorter::mesh( std::vector< std::size_t >& ginpoel,
              tk::UnsMesh::CoordMap& coordmap,
              std::vector< std::size_t >& triinpoel,
              std::map< int, std::vector< std::size_t > >& bnode )
// *****************************************************************************
// Update mesh data we hold for whoever calls this function
//! \param[in,out] ginpoel Mesh connectivity using global IDs
//! \param[in,out] coordmap Map of mesh node coordinates
//! \param[in,out] triinpoel Boundary face-node connectivity
//! \param[in] bnode Node lists of side sets
// *****************************************************************************
{
  ginpoel = m_ginpoel;
  coordmap = m_coordmap;
  triinpoel = m_triinpoel;
  bnode = m_bnode;
}

void
Sorter::createDiscWorkers()
// *****************************************************************************
//  Create Discretization chare array elements on this PE
//! \details We create chare array elements by calling the insert() member
//!   function, which allows specifying the PE on which the array element is
//!   created. and we send each chare array element the chunk of mesh it will
//!   operate on.
// *****************************************************************************
{
  std::vector< CProxy_Discretization > disc;
  for (auto& d : m_scheme) disc.push_back( d.disc() );

  // Create worker array element using Charm++ dynamic chare array element
  // insertion: last arg: PE chare is created on. See also Charm++ manual, Sec.
  // "Dynamic Insertion".

  m_scheme[m_meshid].disc()[ thisIndex ].insert( m_meshid, disc,
    m_scheme[m_meshid].fct(), m_scheme[m_meshid].ale(),
    m_scheme[m_meshid].conjugategradients(), m_host,
    m_meshwriter, m_coordmap, m_el, m_msum, m_nchare );

  contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
              m_cbs.get< tag::discinserted >() );
}

void
Sorter::createWorkers()
// *****************************************************************************
//  Create worker chare array element
// *****************************************************************************
{
  // Make sure (bound) base is already created and accessible
  Assert( m_scheme[m_meshid].disc()[thisIndex].ckLocal() != nullptr,
          "About to pass nullptr" );

  // Create worker array element using Charm++ dynamic chare array element
  // insertion: 1st arg: chare id, other args: Discretization's child ctor args.
  // See also Charm++ manual, Sec. "Dynamic Insertion".

  m_scheme[m_meshid].insert( thisIndex, m_scheme[m_meshid].disc(), m_bface,
                             m_bnode, m_triinpoel );

  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.chcreated();

  contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
              m_cbs.get< tag::workinserted >() );

  // Free up some memory
  tk::destroy( m_ginpoel );
  tk::destroy( m_coordmap );
  tk::destroy( m_bface );
  tk::destroy( m_triinpoel );
  tk::destroy( m_bnode );
  tk::destroy( m_nodeset );
  tk::destroy( m_nodech );
  tk::destroy( m_chnode );
  tk::destroy( m_msum );
  tk::destroy( m_reordcomm );
  tk::destroy( m_newnodes );
  tk::destroy( m_reqnodes );
}

#include "NoWarning/sorter.def.h"
