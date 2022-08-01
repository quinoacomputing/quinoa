// *****************************************************************************
/*!
  \file      src/Inciter/Partitioner.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare partitioner nodegroup used to perform mesh
             partitioning
  \details   Charm++ chare partitioner nodegroup used to perform mesh read and
             partitioning, one worker per compute node.
*/
// *****************************************************************************

#include <numeric>

#include "Partitioner.hpp"
#include "DerivedData.hpp"
#include "Reorder.hpp"
#include "MeshReader.hpp"
#include "CGPDE.hpp"
#include "DGPDE.hpp"
#include "Inciter/Options/Scheme.hpp"
#include "UnsMesh.hpp"
#include "ContainerUtil.hpp"
#include "Callback.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< CGPDE > g_cgpde;
extern std::vector< DGPDE > g_dgpde;

} // inciter::

using inciter::Partitioner;

Partitioner::Partitioner(
  std::size_t meshid,
  const std::string& filename,
  const tk::PartitionerCallback& cbp,
  const tk::RefinerCallback& cbr,
  const tk::SorterCallback& cbs,
  const CProxy_Transporter& host,
  const CProxy_Refiner& refiner,
  const CProxy_Sorter& sorter,
  const tk::CProxy_MeshWriter& meshwriter,
  const std::vector< Scheme >& scheme,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& faces,
  const std::map< int, std::vector< std::size_t > >& bnode ) :
  m_meshid( meshid ),
  m_cbp( cbp ),
  m_cbr( cbr ),
  m_cbs( cbs ),
  m_host( host ),
  m_refiner( refiner ),
  m_sorter( sorter ),
  m_meshwriter( meshwriter ),
  m_scheme( scheme ),
  m_ginpoel(),
  m_coord(),
  m_inpoel(),
  m_lid(),
  m_elemBlockId(),
  m_ndist( 0 ),
  m_nchare( 0 ),
  m_nface(),
  m_chinpoel(),
  m_chcoordmap(),
  m_chbface(),
  m_chtriinpoel(),
  m_chbnode(),
  m_chelemblockid(),
  m_bface( bface ),
  m_bnode( bnode )
// *****************************************************************************
//  Constructor
//! \param[in] meshid Mesh ID
//! \param[in] filename Input mesh filename to read from
//! \param[in] cbp Charm++ callbacks for Partitioner
//! \param[in] cbr Charm++ callbacks for Refiner
//! \param[in] cbs Charm++ callbacks for Sorter
//! \param[in] host Host Charm++ proxy we are being called from
//! \param[in] refiner Mesh refiner proxy
//! \param[in] sorter Mesh reordering (sorter) proxy
//! \param[in] meshwriter Mesh writer proxy
//! \param[in] scheme Discretization scheme
//! \param[in] bface File-internal elem ids of side sets (whole mesh)
//! \param[in] faces Elem-relative face ids of side sets (whole mesh)
//! \param[in] bnode Node lists of side sets (whole mesh)
// *****************************************************************************
{
  // Create mesh reader
  tk::MeshReader mr( filename );

  // Read this compute node's chunk of the mesh (graph and coords) from file
  std::vector< std::size_t > triinpoel;
  mr.readMeshPart( m_ginpoel, m_inpoel, triinpoel, m_lid, m_coord,
                   m_elemBlockId, CkNumNodes(), CkMyNode() );

  // Compute triangle connectivity for side sets, reduce boundary face for side
  // sets to this compute node only and to compute-node-local face ids
  m_triinpoel = mr.triinpoel( m_bface, faces, m_ginpoel, triinpoel );

  // Reduce boundary node lists (global ids) for side sets to this compute node
  // only
  ownBndNodes( m_lid, m_bnode );

  // Sum number of cells across distributed mesh
  std::vector< std::size_t > meshdata{ meshid, m_ginpoel.size()/4 };
  contribute( meshdata, CkReduction::sum_ulong, m_cbp.get< tag::load >() );
}

void
Partitioner::ownBndNodes(
  const std::unordered_map< std::size_t, std::size_t >& lid,
  std::map< int, std::vector< std::size_t > >& bnode )
// *****************************************************************************
// Keep only those nodes for side sets that reside on this compute node
//! \param[in] lid Global->local node IDs of elements of this compute node's
//!   mesh chunk
//! \param[in,out] bnode Global ids of nodes for side sets for whole mesh
//! \details This function overwrites the input boundary node lists map with the
//!    nodes that reside on the caller compute node.
// *****************************************************************************
{
  std::map< int, std::vector< std::size_t > > bnode_own;

  for (const auto& [ setid, nodes ] : bnode) {
    auto& b = bnode_own[ setid ];
    for (auto n : nodes) {
      auto i = lid.find( n );
      if (i != end(lid)) b.push_back( n );
    }
    if (b.empty()) bnode_own.erase( setid );
  }

  bnode = std::move(bnode_own);
}

void
Partitioner::partition( int nchare )
// *****************************************************************************
//  Partition the computational mesh into a number of chares
//! \param[in] nchare Number of parts the mesh will be partitioned into
//! \details This function calls the mesh partitioner to partition the mesh. The
//!   number of partitions equals the number nchare argument which must be no
//!   lower than the number of compute nodes.
// *****************************************************************************
{
  Assert( nchare >= CkNumNodes(), "Number of chares must not be lower than the "
                                  "number of compute nodes" );

  // Generate element IDs for Zoltan
  std::vector< long > gelemid( m_ginpoel.size()/4 );
  std::iota( begin(gelemid), end(gelemid), 0 );

  m_nchare = nchare;
  const auto alg = g_inputdeck.get< tag::selected, tag::partitioner >();
  const auto che = tk::zoltan::geomPartMesh( alg,
                                             centroids( m_inpoel, m_coord ),
                                             gelemid,
                                             nchare );

  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pepartitioned();

  contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
              m_cbp.get< tag::partitioned >() );

  Assert( che.size() == gelemid.size(), "Size of ownership array (chare ID "
          "of elements) after mesh partitioning does not equal the number of "
          "mesh graph elements" );

  // Categorize mesh elements (given by their gobal node IDs) by target chare
  // and distribute to their compute nodes based on mesh partitioning.
  distribute( categorize( che ) );
}

void
Partitioner::addMesh(
  int fromnode,
  const std::unordered_map< int,        // chare id
          std::tuple<
            std::vector< std::size_t >, // tet connectivity
            tk::UnsMesh::CoordMap,      // node coords
            std::unordered_map< int, std::vector< std::size_t > >, // bface conn
            std::unordered_map< int, std::vector< std::size_t > >, // bnodes
            std::vector< int >          // elem-blocks
          > >& chmesh )
// *****************************************************************************
//  Receive mesh associated to chares we own after refinement
//! \param[in] fromnode Compute node call coming from
//! \param[in] chmesh Map associating mesh connectivities to global node ids
//!   and node coordinates for mesh chunks we are assigned by the partitioner
// *****************************************************************************
{
  // Store mesh connectivity and global node coordinates categorized by chares.
  // The send side also writes to the data written here, so concat.
  for (const auto& [ chareid, chunk ] : chmesh) {
    Assert( node(chareid) == CkMyNode(), "Compute node "
            + std::to_string(CkMyNode()) +
            " received a mesh whose chare it does not own" );
    // Store domain element (tetrahedron) connectivity
    const auto& inpoel = std::get< 0 >( chunk );
    auto& inp = m_chinpoel[ chareid ];  // will store tetrahedron connectivity
    inp.insert( end(inp), begin(inpoel), end(inpoel) );
    // Store own mesh block id and associated tet ids
    const auto& elblock = std::get<4>( chunk );
    auto& cheb = m_chelemblockid[ chareid ];
    cheb.insert( end(cheb), begin(elblock), end(elblock) );
    // Store mesh node coordinates associated to global node IDs
    const auto& coord = std::get< 1 >( chunk );
    Assert( tk::uniquecopy(inpoel).size() == coord.size(), "Size mismatch" );
    auto& chcm = m_chcoordmap[ chareid ];     // will store node coordinates
    chcm.insert( begin(coord), end(coord) );  // concatenate node coords
    // Store boundary side set id + face ids + face connectivities
    const auto& bconn = std::get< 2 >( chunk );
    auto& bface = m_chbface[ chareid ];  // for side set id + boundary face ids
    auto& t = m_chtriinpoel[ chareid ];  // for boundary face connectivity
    auto& f = m_nface[ chareid ];        // use counter for chare
    for (const auto& [ setid, faceids ] : bconn) {
      auto& b = bface[ setid ];
      for (std::size_t i=0; i<faceids.size()/3; ++i) {
        b.push_back( f++ );
        t.push_back( faceids[i*3+0] );
        t.push_back( faceids[i*3+1] );
        t.push_back( faceids[i*3+2] );
      }
    }
    // Store boundary side set id + node lists
    const auto& bnode = std::get< 3 >( chunk );
    auto& nodes = m_chbnode[ chareid ];  // for side set id + boundary nodes
    for (const auto& [ setid, bnodes ] : bnode) {
      auto& b = nodes[ setid ];
      b.insert( end(b), begin(bnodes), end(bnodes) );
    }
  }

  thisProxy[ fromnode ].recvMesh();
}

int
Partitioner::node( int id ) const
// *****************************************************************************
//  Return nodegroup id for chare id
//! \param[in] id Chare id
//! \return Nodegroup that creates the chare
//! \details This is computed based on a simple contiguous linear
//!   distribution of chare ids to compute nodes.
// *****************************************************************************
{
  Assert( m_nchare > 0, "Number of chares must be a positive number" );
  auto p = id / (m_nchare / CkNumNodes());
  if (p >= CkNumNodes()) p = CkNumNodes()-1;
  Assert( p < CkNumNodes(), "Assigning to nonexistent node" );
  return p;
}

void
Partitioner::recvMesh()
// *****************************************************************************
//  Acknowledge received mesh chunk and its nodes after mesh refinement
// *****************************************************************************
{
  if (--m_ndist == 0) {
    if (g_inputdeck.get< tag::cmd, tag::feedback >()) m_host.pedistributed();
    contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
                m_cbp.get< tag::distributed >() );
  }
}

void
Partitioner::refine()
// *****************************************************************************
// Optionally start refining the mesh
// *****************************************************************************
{
  auto dist = distribution( m_nchare );

  std::size_t error = 0;
  if (m_chinpoel.size() < static_cast<std::size_t>(dist[1])) {

    error = 1;

  } else {

    for (int c=0; c<dist[1]; ++c) {
      // compute chare ID
      auto cid = CkMyNode() * dist[0] + c;

      Assert(tk::cref_find(m_chinpoel,cid).size()/4 ==
        tk::cref_find(m_chelemblockid,cid).size(),
        "incorrect number of elements in elem-block id vector");
      // create refiner Charm++ chare array element using dynamic insertion
      m_refiner[ cid ].insert( m_meshid,
                               m_host,
                               m_sorter,
                               m_meshwriter,
                               m_scheme,
                               m_cbr,
                               m_cbs,
                               tk::cref_find(m_chinpoel,cid),
                               tk::cref_find(m_chcoordmap,cid),
                               tk::cref_find(m_chbface,cid),
                               tk::cref_find(m_chtriinpoel,cid),
                               tk::cref_find(m_chbnode,cid),
                               tk::cref_find(m_chelemblockid,cid),
                               m_nchare );
    }

  }

  tk::destroy( m_ginpoel );
  tk::destroy( m_coord );
  tk::destroy( m_inpoel );
  tk::destroy( m_lid );
  tk::destroy( m_elemBlockId );
  tk::destroy( m_nface );
  tk::destroy( m_nodech );
  tk::destroy( m_linnodes );
  tk::destroy( m_chinpoel );
  tk::destroy( m_chcoordmap );
  tk::destroy( m_chbface );
  tk::destroy( m_chtriinpoel );
  tk::destroy( m_chbnode );
  tk::destroy( m_chelemblockid );
  tk::destroy( m_bnodechares );
  tk::destroy( m_bface );
  tk::destroy( m_triinpoel );
  tk::destroy( m_bnode );

  std::vector< std::size_t > meshdata{ m_meshid, error };
  contribute( meshdata, CkReduction::max_ulong, m_cbp.get<tag::refinserted>() );
}

std::array< std::vector< tk::real >, 3 >
Partitioner::centroids( const std::vector< std::size_t >& inpoel,
                        const tk::UnsMesh::Coords& coord )
// *****************************************************************************
//  Compute element centroid coordinates
//! \param[in] inpoel Mesh connectivity with local ids
//! \param[in] coord Node coordinates
//! \return Centroids for all cells on this compute node
// *****************************************************************************
{
  Assert( tk::uniquecopy(inpoel).size() == coord[0].size(), "Size mismatch" );

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Make room for element centroid coordinates
  std::array< std::vector< tk::real >, 3 > cent;
  auto& cx = cent[0];
  auto& cy = cent[1];
  auto& cz = cent[2];
  auto num = inpoel.size()/4;
  cx.resize( num );
  cy.resize( num );
  cz.resize( num );

  // Compute element centroids for mesh passed in
  for (std::size_t e=0; e<num; ++e) {
    auto A = inpoel[e*4+0];
    auto B = inpoel[e*4+1];
    auto C = inpoel[e*4+2];
    auto D = inpoel[e*4+3];
    cx[e] = (x[A] + x[B] + x[C] + x[D]) / 4.0;
    cy[e] = (y[A] + y[B] + y[C] + y[D]) / 4.0;
    cz[e] = (z[A] + z[B] + z[C] + z[D]) / 4.0;
  }

  return cent;
}

std::unordered_map< int, Partitioner::MeshData >
Partitioner::categorize( const std::vector< std::size_t >& target ) const
// *****************************************************************************
// Categorize mesh data by target
//! \param[in] target Target chares of mesh elements, size: number of
//!   elements in the chunk of the mesh graph on this compute node.
//! \return Vector of global mesh node ids connecting elements owned by each
//!   target chare.
// *****************************************************************************
{
  Assert( target.size() == m_ginpoel.size()/4, "Size mismatch");

  using Face = tk::UnsMesh::Face;

  // Build hash map associating side set id to boundary faces
  std::unordered_map< Face, int,
                      tk::UnsMesh::Hash<3>, tk::UnsMesh::Eq<3> > faceside;
  for (const auto& [ setid, faceids ] : m_bface)
    for (auto f : faceids)
      faceside[ {{ m_triinpoel[f*3+0],
                   m_triinpoel[f*3+1],
                   m_triinpoel[f*3+2] }} ] = setid;

  // Build hash map associating side set ids to boundary nodes
  std::unordered_map< std::size_t, std::unordered_set< int > > nodeside;
  for (const auto& [ setid, nodes ] : m_bnode)
    for (auto n : nodes)
      nodeside[ n ].insert( setid );

  // Categorize mesh data (tets, node coordinates, and boundary data) by target
  // chare based on which chare the partitioner assigned elements (tets) to
  std::unordered_map< int, MeshData > chmesh;
  for (std::size_t e=0; e<target.size(); ++e) {
    // Construct a tetrahedron with global node ids
    tk::UnsMesh::Tet t{{ m_ginpoel[e*4+0], m_ginpoel[e*4+1],
                         m_ginpoel[e*4+2], m_ginpoel[e*4+3] }};
    // Categorize tetrahedron (domain element) connectivity
    auto& mesh = chmesh[ static_cast<int>(target[e]) ];
    auto& inpoel = std::get< 0 >( mesh );
    inpoel.insert( end(inpoel), begin(t), end(t) );
    // Categorize mesh block ids and associated local tet ids (local tet ids
    // basically correspond to the inpoel entries updated in line above).
    bool added(false);
    auto& elblock = std::get< 3 >( mesh );
    for (const auto& [blid, elset] : m_elemBlockId) {
      if (elset.find(e) != elset.end()) {
        elblock.push_back(blid);
        added = true;
      }
    }
    Assert(added, "Tet element not found in any block after partitioning");
    // Categorize boundary face connectivity
    auto& bconn = std::get< 1 >( mesh );
    std::array<Face,4> face{{ {{t[0],t[2],t[1]}}, {{t[0],t[1],t[3]}},
                              {{t[0],t[3],t[2]}}, {{t[1],t[2],t[3]}} }};
    for (const auto& f : face) {
      auto it = faceside.find( f );
      if (it != end(faceside)) {
        auto& s = bconn[ it->second ];
        s.insert( end(s), begin(f), end(f) );
      }
    }
    // Categorize boundary node lists
    auto& bnode = std::get< 2 >( mesh );
    for (const auto& n : t) {
      auto it = nodeside.find( n );
      if (it != end(nodeside))
        for (auto s : it->second)
          bnode[ s ].push_back( n );
    }
  }

  // Make boundary node lists unique per side set
  for (auto& c : chmesh)
    for (auto& n : std::get<2>(c.second))
       tk::unique( n.second );

  // Make sure all compute nodes have target chares assigned
  Assert( !chmesh.empty(), "No elements have been assigned to a chare" );

  // This check should always be done, hence ErrChk and not Assert, as it
  // can result from particular pathological combinations of (1) too large
  // degree of virtualization, (2) too many compute nodes, and/or (3) too small
  // of a mesh and not due to programmer error.
  for(const auto& c : chmesh)
    ErrChk( !std::get<0>(c.second).empty(),
            "Overdecomposition of the mesh is too large compared to the "
            "number of work units computed based on the degree of "
            "virtualization desired. As a result, there would be at least "
            "one work unit with no mesh elements to work on, i.e., nothing "
            "to do. Solution 1: decrease the virtualization to a lower "
            "value using the command-line argument '-u'. Solution 2: "
            "decrease the number processing elements (PEs and/or compute "
            "nodes) using the charmrun command-line argument '+pN' where N is "
            "the number of PEs (or in SMP-mode in combination with +ppn to "
            "reduce the number of compute nodes), which implicitly increases "
            "the size (and thus decreases the number) of work units.)" );

  return chmesh;
}

tk::UnsMesh::CoordMap
Partitioner::coordmap( const std::vector< std::size_t >& inpoel )
// *****************************************************************************
// Extract coordinates associated to global nodes of a mesh chunk
//! \param[in] inpoel Mesh connectivity
//! \return Map storing the coordinates of unique nodes associated to global
//!    node IDs in mesh given by inpoel
// *****************************************************************************
{
  Assert( inpoel.size() % 4 == 0, "Incomplete mesh connectivity" );

  tk::UnsMesh::CoordMap map;

  for (auto g : tk::uniquecopy(inpoel)) {
     auto i = tk::cref_find( m_lid, g );
     auto& c = map[g];
     c[0] = m_coord[0][i];
     c[1] = m_coord[1][i];
     c[2] = m_coord[2][i];
  }

  Assert( tk::uniquecopy(inpoel).size() == map.size(), "Size mismatch" );

  return map;
}

void
Partitioner::distribute( std::unordered_map< int, MeshData >&& mesh )
// *****************************************************************************
// Distribute mesh to target compute nodes after mesh partitioning
//! \param[in] mesh Mesh data categorized by target by target chares
// *****************************************************************************
{
  auto dist = distribution( m_nchare );

  // Extract mesh data whose chares are on ("owned by") this compute node
  for (int c=0; c<dist[1]; ++c) {
    auto chid = CkMyNode() * dist[0] + c; // compute owned chare ID
    const auto it = mesh.find( chid );    // attempt to find its mesh data
    if (it != end(mesh)) {                // if found
      // Store own tetrahedron connectivity
      const auto& inpoel = std::get<0>( it->second );
      auto& inp = m_chinpoel[ chid ];     // will store own mesh connectivity
      inp.insert( end(inp), begin(inpoel), end(inpoel) );
      // Store own mesh block id and associated tet ids
      const auto& elblock = std::get<3>( it->second );
      auto& cheb = m_chelemblockid[ chid ];
      cheb.insert( end(cheb), begin(elblock), end(elblock) );
      // Store own node coordinates
      auto& chcm = m_chcoordmap[ chid ];  // will store own node coordinates
      auto cm = coordmap( inpoel );       // extract node coordinates 
      chcm.insert( begin(cm), end(cm) );  // concatenate node coords
      // Store own boundary face connectivity
      const auto& bconn = std::get<1>( it->second );
      auto& bface = m_chbface[ chid ];    // will store own boundary faces
      auto& t = m_chtriinpoel[ chid ];    // wil store own boundary face conn
      auto& f = m_nface[ chid ];          // use counter for chare
      for (const auto& [ setid, faceids ] : bconn) {
        auto& b = bface[ setid ];
        for (std::size_t i=0; i<faceids.size()/3; ++i) {
          b.push_back( f++ );
          t.push_back( faceids[i*3+0] );
          t.push_back( faceids[i*3+1] );
          t.push_back( faceids[i*3+2] );
        }
      }
      // Store own boundary node lists
      const auto& bnode = std::get<2>( it->second );
      auto& nodes = m_chbnode[ chid ];    // will store own boundary nodes
      for (const auto& [ setid, nodeids ] : bnode) {
        auto& b = nodes[ setid ];
        b.insert( end(b), begin(nodeids), end(nodeids) );
      }
      // Remove chare ID and mesh data
      mesh.erase( it );
    }
    Assert( mesh.find(chid) == end(mesh), "Not all owned mesh data stored" );
  }

  // Construct export map (associating mesh connectivities with global node
  // indices and node coordinates) for mesh chunks associated to chare IDs
  // owned by chares we do not own.
  std::unordered_map< int,                     // target compute node
    std::unordered_map< int,                   // chare ID
      std::tuple<
        // (domain-element) tetrahedron connectivity
        std::vector< std::size_t >,
        // (domain) node IDs & coordinates
        tk::UnsMesh::CoordMap,
        // boundary side set + face connectivity
        std::unordered_map< int, std::vector< std::size_t > >,
        // boundary side set + node list
        std::unordered_map< int, std::vector< std::size_t > >,
        // Mesh block ids + local tet ids
        std::vector< int >
      > > > exp;

  for (const auto& c : mesh)
    exp[ node(c.first) ][ c.first ] =
      std::make_tuple( std::get<0>(c.second),
                       coordmap(std::get<0>(c.second)),
                       std::get<1>(c.second),
                       std::get<2>(c.second),
                       std::get<3>(c.second) );

  // Export chare IDs and mesh we do not own to fellow compute nodes
  if (exp.empty()) {
    if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pedistributed();
    contribute( sizeof(std::size_t), &m_meshid, CkReduction::nop,
                m_cbp.get< tag::distributed >() );
  } else {
     m_ndist += exp.size();
     for (const auto& [ targetnode, chunk ] : exp)
       thisProxy[ targetnode ].addMesh( CkMyNode(), chunk );
  }
}

std::array< int, 2 >
Partitioner::distribution( int npart ) const
// *****************************************************************************
//  Compute chare (partition) distribution
//! \param[in] npart Total number of chares (partitions) to distribute
//! \return Chunksize, i.e., number of chares per all compute nodes except the
//!   last one, and the number of chares for this compute node.
//! \details Chare ids are distributed to compute nodes in a linear continguous
//!   order with the last compute node taking the remainder if the number of
//!   compute nodes is not divisible by the number chares. For example, if
//!   nchare=7 and nnode=3, the chare distribution is node0: 0 1, node1: 2 3,
//!   and node2: 4 5 6. As a result of this distribution, all compute nodes will
//!   have their chare-categorized element connectivity filled with the global
//!   mesh node IDs associated to the Charm++ chare IDs each compute node owns.
// *****************************************************************************
{
  auto chunksize = npart / CkNumNodes();
  auto mynchare = chunksize;
  if (CkMyNode() == CkNumNodes()-1) mynchare += npart % CkNumNodes();
  return {{ chunksize, mynchare }};
}

#include "NoWarning/partitioner.def.h"
