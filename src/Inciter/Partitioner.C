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

#include <algorithm>

#include "Partitioner.h"
#include "DerivedData.h"
#include "Reorder.h"
#include "MeshReader.h"
#include "Around.h"
#include "ExodusIIMeshWriter.h"
#include "CGPDE.h"
#include "DGPDE.h"
#include "AMR/Error.h"
#include "Inciter/Options/Scheme.h"
#include "Inciter/Options/AMRInitial.h"
#include "UnsMesh.h"
#include "ContainerUtil.h"
#include "Callback.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< CGPDE > g_cgpde;
extern std::vector< DGPDE > g_dgpde;

} // inciter::

using inciter::Partitioner;

Partitioner::Partitioner(
  const tk::PartitionerCallback& cbp,
  const tk::RefinerCallback& cbr,
  const tk::SorterCallback& cbs,
  const CProxy_Transporter& host,
  const tk::CProxy_Solver& solver,
  const CProxy_Refiner& refiner,
  const CProxy_Sorter& sorter,
  const Scheme& scheme,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel ) :
  m_cbp( cbp ),
  m_cbr( cbr ),
  m_cbs( cbs ),
  m_host( host ),
  m_solver( solver ),
  m_refiner( refiner ),
  m_sorter( sorter ),
  m_scheme( scheme ),
  m_el(),
  m_reqNodes(),
  m_start( 0 ),
  m_ndist( 0 ),
  m_noffset( 0 ),
  m_nquery( 0 ),
  m_nmask( 0 ),
  m_ginpoel(),
  m_coord(),
  m_coordmap(),
  m_nchare( 0 ),
  m_chinpoel(),
  m_chcoordmap(),
  m_bnodechares(),
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
  MeshReader mr( g_inputdeck.get< tag::cmd, tag::io, tag::input >() );

  // Read our chunk of the mesh graph from file
  mr.readGraph( m_ginpoel, CkNumPes(), CkMyPe() );

  // Compute local data from global mesh connectivity (m_inpoel, m_gid, m_lid)
  m_el = tk::global2local( m_ginpoel );

  // Read our chunk of the mesh node coordinates from file
  m_coord = mr.readCoords( m_gid );

  // Compute number of cells across whole problem
  auto nelem = m_ginpoel.size()/4;
  contribute( sizeof(uint64_t), &nelem, CkReduction::sum_int,
              m_cbp.get< tag::load >() );
}

void
Partitioner::partition( int nchare )
// *****************************************************************************
//  Partition the computational mesh into a number of chares
//! \param[in] nchare Number of parts the mesh will be partitioned into
//! \details This function calls the mesh partitioner to partition the mesh. The
//!   number of partitions equals the number nchare argument which must be no
//!   lower than the number of PEs.
// *****************************************************************************
{
  Assert( nchare >= CkNumPes(),
          "Number of chares must not be lower than the number of PEs" );

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

  Assert( che.size() == gelemid.size(), "Size of ownership array (chare ID "
          "of elements) after mesh partitioning does not equal the number of "
          "mesh graph elements" );

  m_chcoordmap.clear();
  m_chinpoel.clear();

  // Categorize mesh elements (given by their gobal node IDs) by target chare
  // and distribute to their PEs based on mesh partitioning.
  distribute( categorize( che, m_ginpoel ) );
}

void
Partitioner::addMesh( int frompe,
                      const std::unordered_map< int,
                              std::tuple< std::vector< std::size_t >,
                                          tk::UnsMesh::CoordMap > >& chmesh )
// *****************************************************************************
//  Receive mesh associated to chares we own after refinement
//! \param[in] frompe PE call coming from
//! \param[in] cmesh Map associating mesh connectivities to global node ids
//!   and node coordinates for mesh chunks we are assigned by the partitioner
// *****************************************************************************
{
  // Store mesh connectivity and global node coordinates categorized by chares
  for (const auto& c : chmesh) {
    Assert( pe(c.first) == CkMyPe(), "PE " + std::to_string(CkMyPe()) +
            " received a mesh whose chare it does not own" );
    // The send side also writes to this so append/concat
    auto& inpoel = m_chinpoel[ c.first ];
    const auto& mesh = std::get< 0 >( c.second );
    inpoel.insert( end(inpoel), begin(mesh), end(mesh) );
    // Store coordinates associated to global node IDs. The send side also
    // writes to this, so concat.
    const auto& coord = std::get< 1 >( c.second );
    Assert( tk::uniquecopy(mesh).size() == coord.size(), "Size mismatch" );
    auto& chcm = m_chcoordmap[ c.first ];  // store node coordinates per chare
    chcm.insert( begin(coord), end(coord) );  // concatenate node coords
  }

  thisProxy[ frompe ].recvMesh();
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
Partitioner::recvMesh()
// *****************************************************************************
//  Acknowledge received mesh chunk and its nodes after mesh refinement
// *****************************************************************************
{
  if (--m_ndist == 0) {
    if (g_inputdeck.get< tag::cmd, tag::feedback >()) m_host.pedistributed();
    contribute( m_cbp.get< tag::distributed >() );
  }
}

void
Partitioner::refine()
// *****************************************************************************
// Optionally start refining the mesh
// *****************************************************************************
{
  auto dist = distribution( m_nchare );

  //m_cbr.get< tag::update >() =
  //  CkCallback( CkIndex_Partitioner::updateBoundaryMesh(nullptr), thisProxy );

  for (int c=0; c<dist[1]; ++c) {
    // compute chare ID
    auto cid = CkMyPe() * dist[0] + c;
    // create refiner Charm++ chare array element using dynamic insertion
    m_refiner[ cid ].insert( m_host,
                             m_sorter,
                             m_solver,
                             m_scheme,
                             m_cbr,
                             m_cbs,
                             tk::cref_find(m_chinpoel,cid),
                             tk::cref_find(m_chcoordmap,cid),
                             m_bface,
                             m_triinpoel,
                             m_nchare,
                             CkMyPe() );
  }

  contribute( m_cbp.get< tag::refinserted >() );
}

std::array< std::vector< tk::real >, 3 >
Partitioner::centroids( const std::vector< std::size_t >& inpoel,
                        const tk::UnsMesh::Coords& coord )
// *****************************************************************************
//  Compute element centroid coordinates
//! \param[in] inpoel Mesh connectivity with local ids
//! \param[ib] coord Node coordinates
//! \return Centroids for all cells on this PE
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

std::unordered_map< int, std::vector< std::size_t > >
Partitioner::categorize( const std::vector< std::size_t >& target,
                         const std::vector< std::size_t >& inpoel ) const
// *****************************************************************************
// Categorize mesh elements (given by their gobal node IDs) by target
//! \param[in] target Targets (chares or PEs) of mesh elements, size: number of
//!   elements in the chunk of the mesh graph on this PE.
//! \param[in] inpoel Mesh connectivity to distribute elemets from
//! \return Vector of global mesh node ids connecting elements owned by each
//!   target (chare or PE)
// *****************************************************************************
{
  Assert( target.size() == inpoel.size()/4, "Size mismatch");

  // Categorize global mesh node ids of elements by target
  std::unordered_map< int, std::vector< std::size_t > > nodes;
  for (std::size_t e=0; e<target.size(); ++e) {
    auto& c = nodes[ static_cast<int>(target[e]) ];
    for (std::size_t n=0; n<4; ++n) c.push_back( inpoel[e*4+n] );
  }

  // Make sure all PEs have targets assigned
  Assert( !nodes.empty(), "No nodes have been assigned to chares on PE " );

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
Partitioner::distribute(
 std::unordered_map< int, std::vector< std::size_t > >&& elems )
// *****************************************************************************
// Distribute mesh to target PEs after mesh partitioning
//! \param[in] elems Mesh cells (with global node IDs) categorized by target
//!   chares
// *****************************************************************************
{
  auto dist = distribution( m_nchare );

  // Extract those mesh connectivities whose chares live on this PE
  for (int c=0; c<dist[1]; ++c) {
    auto chid = CkMyPe() * dist[0] + c;   // compute owned chare ID
    const auto it = elems.find( chid );   // attempt to find its nodes
    if (it != end(elems)) {               // if found
      auto& inp = m_chinpoel[ chid ];     // store own mesh
      inp.insert( end(inp), begin(it->second), end(it->second) );
      auto& chcm = m_chcoordmap[ chid ];         // store own node coordinates
      auto cm = coordmap( it->second );   // extract node coordinates 
      chcm.insert( begin(cm), end(cm) );  // concatenate node coords
      elems.erase( it );                  // remove chare ID and nodes
    }
    Assert( elems.find(chid) == end(elems), "Not all owned node IDs stored" );
  }

  // Construct export map associating mesh connectivities with global node
  // indices and node coordinates for mesh chunks associated to chare IDs (inner
  // key) owned by chares we do not own.
  std::unordered_map< int,                              // target PE
    std::unordered_map< int,                            // chare ID
      std::tuple< std::vector< std::size_t >,           // mesh connectivity
                  tk::UnsMesh::CoordMap > > > exp;      // node ID & coords
  for (const auto& c : elems)
    exp[ pe(c.first) ][ c.first ] =
      std::make_tuple( c.second, coordmap(c.second) );

  // Export chare IDs and mesh we do not own to fellow PEs
  if (exp.empty()) {
    if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pedistributed();
    contribute( m_cbp.get< tag::distributed >() );
  } else {
     m_ndist = exp.size();
     for (const auto& p : exp)
       thisProxy[ p.first ].addMesh( CkMyPe(), p.second );
  }
}

std::array< int, 2 >
Partitioner::distribution( int npart ) const
// *****************************************************************************
//  Compute chare (partition) distribution
//! \param[in] npart Total number of chares (partitions) to distribute
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
  auto chunksize = npart / CkNumPes();
  auto mynchare = chunksize;
  if (CkMyPe() == CkNumPes()-1) mynchare += npart % CkNumPes();
  return {{ chunksize, mynchare }};
}

#include "NoWarning/partitioner.def.h"
