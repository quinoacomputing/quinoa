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
#include "AMR/mesh_adapter.h"
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

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< CGPDE > g_cgpde;
extern std::vector< DGPDE > g_dgpde;

} // inciter::

using inciter::Partitioner;

Partitioner::Partitioner(
  const std::vector< CkCallback >& cb,
  const CProxy_Transporter& host,
  const tk::CProxy_Solver& solver,
  const Scheme& scheme,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel ) :
  m_cb( cb[0], cb[1], cb[2], cb[3], cb[4], cb[5], cb[6], cb[7] ),
  m_host( host ),
  m_solver( solver ),
  m_scheme( scheme ),
  m_npeDist( 0 ),
  m_npe( 0 ),
  m_ndist( 0 ),
  m_nedge( 0 ),
  m_nref( 0 ),
  m_extra( 1 ),
  m_pe(),
  m_initref(),
  m_el(),
  m_edgenode(),
  m_bndEdges(),
  m_edgenodePe(),
  m_reqNodes(),
  m_start( 0 ),
  m_noffset( 0 ),
  m_nquery( 0 ),
  m_nmask( 0 ),
  m_ginpoel(),
  m_coord(),
  m_coordmap(),
  m_nchare( 0 ),
  m_lower( 0 ),
  m_upper( 0 ),
  m_ncomm(),
  m_ncommunication(),
  m_nodeset(),
  m_nodech(),
  m_linnodes(),
  m_chinpoel(),
  m_chcoordmap(),
  m_chfilenodes(),
  m_cost( 0.0 ),
  m_bnodechares(),
  m_msum(),
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

  // Send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.peread();

  // Store initial mesh refinement type list in reverse order
  m_initref = g_inputdeck.get< tag::amr, tag::init >();
  std::reverse( begin(m_initref), end(m_initref) );

  if ( !g_inputdeck.get< tag::amr, tag::init >().empty() ||
       !g_inputdeck.get< tag::amr, tag::edge >().empty() )
    partref();          // if initial mesh refinement configured, partition
  else
    finishref();        // if not, continue
}

void
Partitioner::partref()
// *****************************************************************************
//  Partition the mesh to NumPes partitions (before an initial refinement step)
//! \details This function calls the mesh partitioner to partition (or
//!   re-partition) the (current) mesh as a first step for an initial mesh
//!   refinement step. The number of partitions always equals the numher of PEs.
// *****************************************************************************
{
  // Generate element IDs for Zoltan
  std::vector< long > gelemid( m_inpoel.size()/4 );
  std::iota( begin(gelemid), end(gelemid), 0 );

{auto initref = g_inputdeck.get< tag::amr, tag::init >();
auto level = initref.size() - m_initref.size();
tk::UnsMesh rm( m_inpoel, m_coord );
tk::ExodusIIMeshWriter mwr( "initref.partref." + std::to_string(level) + '.' +
                             std::to_string(CkMyPe()),
                            tk::ExoWriter::CREATE );

mwr.writeMesh( rm );}

  // Partition the mesh using Zoltan to number of PEs parts
  const auto alg = g_inputdeck.get< tag::selected, tag::partitioner >();
  const auto pel = tk::zoltan::geomPartMesh( alg,
                                             centroids( m_inpoel, m_coord ),
                                             gelemid,
                                             CkNumPes() );

  Assert( pel.size() == gelemid.size(), "Size of ownership array (PE of "
          "elements) after mesh partitioning does not equal the number of mesh "
          "graph elements" );

  // Prepare for a step of initial mesh refinement
  m_extra = 1;
  m_bndEdges.clear();
  m_pe.clear();
  m_edgenodePe.clear();
  m_edgenode.clear();
  m_coordmap.clear();
  auto g = m_ginpoel;
  m_ginpoel.clear();

  // Categorize mesh elements (given by their gobal node IDs) by target PE and
  // distribute to their PEs based on mesh partitioning.
  distributePe( categorize( pel, g ) );
}

void
Partitioner::partchare( int nchare )
// *****************************************************************************
//  Partition the computational mesh into a number of chares
//! \param[in] nchare Number of parts the mesh will be partitioned into
//! \details This function calls the mesh partitioner to partition the (current)
//!   mesh after potentially a number of initial mesh refinement steps. The
//!   number of partitions equals the number nchare argument which could be
//!   larger but not smaller than the number of PEs.
// *****************************************************************************
{
  Assert( nchare >= CkNumPes(),
          "Number of chares must not be lower than the number of PEs" );

  // Generate element IDs for Zoltan
  std::vector< long > gelemid( m_ginpoel.size()/4 );
  std::iota( begin(gelemid), end(gelemid), 0 );

tk::UnsMesh rm( m_inpoel, m_coord );
tk::ExodusIIMeshWriter mwr( "initref.partchare." + std::to_string(CkMyPe()),
                            tk::ExoWriter::CREATE );

mwr.writeMesh( rm );

  m_nchare = nchare;
  const auto alg = g_inputdeck.get< tag::selected, tag::partitioner >();
  const auto che = tk::zoltan::geomPartMesh( alg,
                                             centroids( m_inpoel, m_coord ),
                                             gelemid,
                                             nchare );

  // send progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pepartitioned();

  Assert( che.size() == gelemid.size(), "Size of ownership array (chare ID "
          "of elements) after mesh partitioning does not equal the number of "
          "mesh graph elements" );

  m_coordmap.clear();
  m_chinpoel.clear();

  // Categorize mesh elements (given by their gobal node IDs) by target chare
  // and distribute to their PEs based on mesh partitioning.
  distributeCh( categorize( che, m_ginpoel ) );
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
Partitioner::neworder( const std::unordered_map< std::size_t,
                        std::tuple< std::size_t, tk::UnsMesh::Coord > >& nodes )
// *****************************************************************************
//  Receive new (reordered) global node IDs
//! \param[in] nd Map associating new to old node IDs
// *****************************************************************************
{
  // Signal to the runtime system that we have participated in reordering
  participated_complete();

  // Store new node IDs associated to old ones, and node coordinates associated
  // to new node IDs. Since multiple chares can contribute to a single node, we
  // store such shared node coordinates for all chares that contribute.
  for (const auto& p : nodes) {
    auto id = std::get< 0 >( p.second );
    auto coord = std::get< 1 >( p.second );
    m_linnodes[ p.first ] = id;
    for (auto c : tk::cref_find(m_nodech,p.first))
      m_chcoordmap[ c ].emplace( id, coord );
  }

  // If all our nodes have new IDs assigned, reorder complete on this PE
  if (m_linnodes.size() == m_nodeset.size()) reordered();
}

void
Partitioner::addChMesh( int frompe,
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
    m_coordmap.insert( begin(coord), end(coord) );
  }

  thisProxy[ frompe ].recvChMesh();
}

void
Partitioner::recvChMesh()
// *****************************************************************************
//  Acknowledge received mesh chunk and its nodes after mesh refinement
// *****************************************************************************
{
  if (--m_npe == 0) {
    contribute( m_cb.get< tag::distributed >() );
    if ( g_inputdeck.get< tag::cmd, tag::feedback >() )
      m_host.pedistributed();
  }
}

void
Partitioner::addPeMesh( int frompe,
                        const std::vector< std::size_t >& inpoel,
                        const tk::UnsMesh::CoordMap& cm )
// *****************************************************************************
//! Receive mesh elements and their node coordinates after partitioning
//! \param[in] frompe PE call coming from
//! \param[in] inpoel Mesh connectivity with global node IDs
//! \param[in] cm Coordinates associated to global node IDs
// *****************************************************************************
{
  Assert( tk::uniquecopy(inpoel).size() == cm.size(), "Size mismatch" );

  // Store mesh connectivity. The send side also writes to this, so concat.
  m_ginpoel.insert( end(m_ginpoel), begin(inpoel), end(inpoel) );

  // Store coordinates associated to global node IDs. The send side also writes
  // to this, so concat.
  m_coordmap.insert( begin(cm), end(cm) );

  // Acknowledge receipt of mesh
  thisProxy[ frompe ].recvPeMesh();
}

void
Partitioner::recvPeMesh()
// *****************************************************************************
//  Acknowledge received mesh chunk and its nodes
// *****************************************************************************
{
  if (++m_ndist == m_npeDist) contribute( m_cb.get< tag::refdistributed >() );
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
  std::size_t l = 0;
  for (const auto& i : m_chinpoel) {      // for all meshes passed in
    // generate local ids and connectivity from global connectivity
    auto el = tk::global2local( i.second );
    const auto& inpoel = std::get< 0 >( el );   // local connectivity
    const auto& lid = std::get< 2 >( el );      // global->local node ids
    auto np = tk::npoin( inpoel );
    tk::UnsMesh::Coords coord;
    coord[0].resize( np );
    coord[1].resize( np );
    coord[2].resize( np );
    for (const auto& c : m_coordmap) {
      auto p = lid.find( c.first );
      if (p != lid.end()) {
        Assert( p->second < np, "Indexing out of coord vector" );
        coord[0][p->second] = c.second[0];
        coord[1][p->second] = c.second[1];
        coord[2][p->second] = c.second[2];
      }
    }
    tk::UnsMesh refmesh( inpoel, coord );
    tk::ExodusIIMeshWriter mw( "initref.flatten." + std::to_string(l++) + '.' +
                               std::to_string(CkMyPe()),
                               tk::ExoWriter::CREATE );
    mw.writeMesh( refmesh );
  }

  // Make sure all cells of all chare meshes have non-negative Jacobians
  Assert( positiveJacobians( m_chinpoel, m_coordmap ),
          "Chare-partitioned mesh cell Jacobian non-positive" );

  // Make sure we are not fed garbage
  Assert( m_chinpoel.size() ==
            static_cast< std::size_t >( distribution(m_nchare)[1] ),
          "Global mesh nodes ids associated to chares on PE " +
          std::to_string( CkMyPe() ) + " is incomplete" );

  // Flatten node IDs of elements our chares operate on and associate global
  // node ID to chare IDs.
  for (const auto& c : m_chinpoel)
    for (auto i : c.second) {
      m_nodeset.insert( i );
      m_nodech[ i ].push_back( c.first );
    }

  m_chcoordmap.clear();

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


bool
Partitioner::positiveJacobians(
  const std::unordered_map< int, std::vector< std::size_t > > chinpoel,
  const tk::UnsMesh::CoordMap& cm )
// *****************************************************************************
// Test for positivity of the Jacobian for all cells in multiple meshes
//! \param[in] chinpoel Connectivities of multiple meshes assigned to their
//!   chare IDs
//! \param[in] cm Coordinates associated to global node IDs of of all meshes
//!   passed in
//! \return True if Jacobians of all cells of all meshes are positive
// *****************************************************************************
{
  for (const auto& i : chinpoel) {      // for all meshes passed in
    // generate local ids and connectivity from global connectivity
    auto el = tk::global2local( i.second );
    const auto& inpoel = std::get< 0 >( el );   // local connectivity
    const auto& lid = std::get< 2 >( el );      // global->local node ids
    auto np = tk::npoin( inpoel );
    tk::UnsMesh::Coords coord;
    coord[0].resize( np );
    coord[1].resize( np );
    coord[2].resize( np );
    for (const auto& c : cm) {
      auto p = lid.find( c.first );
      if (p != lid.end()) {
        Assert( p->second < np, "Indexing out of coord vector" );
        coord[0][p->second] = c.second[0];
        coord[1][p->second] = c.second[1];
        coord[2][p->second] = c.second[2];
      }
    }
    if (!tk::positiveJacobians( inpoel, coord )) return false;
  }

  return true;
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
Partitioner::mask( int p, const std::unordered_map< std::size_t,
                                  std::vector< int > >& cn )
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
    // Count up total number of nodes we will need receive during reordering
    std::size_t nrecv = 0;
    for (const auto& u : m_ncommunication) nrecv += u.second.size();

    // send progress report to host
    if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pemask();

    // Compute number of mesh node IDs we will assign IDs to
    auto nuniq = m_nodeset.size() - nrecv;

    // Start computing PE offsets for node reordering
    thisProxy.offset( CkMyPe(), nuniq );
  }
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
Partitioner::distributePe(
  std::unordered_map< int, std::vector< std::size_t > >&& elems )
// *****************************************************************************
// Distribute mesh to their PEs during initial mesh refinement
//! \param[in] elems Mesh cells (with global node IDs) categorized by target PEs
// *****************************************************************************
{
  // Store own mesh connectivity
  auto i = elems.find( CkMyPe() );
  if (i != end(elems)) {
    // Store our mesh chunk. The receive side also writes to this, so concat.
    m_ginpoel.insert( end(m_ginpoel), begin(i->second), end(i->second) );
    // store coordinates associated to global nodes of our mesh chunk
    auto cm = coordmap( i->second );
    // the receive side also writes to this, so concatenate
    m_coordmap.insert( begin(cm), end(cm) );
    // remove our mesh chunk from list (the rest will be exported)
    elems.erase( i );
  }

  m_nedge = 0;

  // Export connectivities to other PEs
  if (elems.empty())
    contribute( m_cb.get< tag::refdistributed >() );
  else {
    m_ndist = 0;
    m_npeDist = elems.size();
    for (const auto& p : elems)
      thisProxy[ p.first ].addPeMesh( CkMyPe(), p.second, coordmap(p.second) );
  }
}

void
Partitioner::distributeCh(
 std::unordered_map< int, std::vector< std::size_t > >&& elems )
// *****************************************************************************
// Distribute mesh chunk to their PEs after initial mesh refinement
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
      auto cm = coordmap( it->second );   // extract node coordinates 
      m_coordmap.insert( begin(cm), end(cm) );  // concatenate node coords
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
    contribute( m_cb.get< tag::distributed >() );
    // send progress report to host
    if ( g_inputdeck.get< tag::cmd, tag::feedback >() )
      m_host.pedistributed();
  } else {
     m_npe = exp.size();
     for (const auto& p : exp)
       thisProxy[ p.first ].addChMesh( CkMyPe(), p.second );
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

void
Partitioner::reorder()
// *****************************************************************************
//  Reorder global mesh node IDs
// *****************************************************************************
{
  // Activate SDAG waits for having requests arrive from other PEs for some
  // of our node IDs; and for computing/receiving lower and upper bounds of
  // global node IDs our PE operates on after reordering
  thisProxy[ CkMyPe() ].wait4prep();
  thisProxy[ CkMyPe() ].wait4bounds();

  // In serial signal to the runtime system that we have participated in
  // reordering. This is required here because this is only triggered if
  // communication is required during mesh node reordering. See also
  // particioner.ci.
  if (CkNumPes() == 1) participated_complete();

  // Send out request for new global node IDs for nodes we do not reorder
  for (const auto& c : m_ncommunication)
    thisProxy[ c.first ].request( CkMyPe(), c.second );

  // Lambda to decide if node ID is being assigned a new ID by us
  auto ownnode = [ this ]( std::size_t p ) {
    using Set = typename std::remove_reference<
                  decltype(m_ncommunication) >::type::value_type;
    return !std::any_of( m_ncommunication.cbegin(), m_ncommunication.cend(),
                         [&](const Set& s)
                         { return s.second.find(p) != s.second.cend(); } );
  };

  // Reorder our chunk of the mesh node IDs by looping through all of our node
  // IDs. We test if this PE is to assign a new ID to a node ID, and if so, we
  // assign a new ID, i.e., reorder, by constructing a map associating new to
  // old IDs. We also count up the reordered nodes, which also serves as the new
  // node id. Also, we store the coordinates associated to the new node ID for
  // each chare on this PE. Since multiple chares can contribute to a single
  // node, we store such shared node coordinates for all chares that contribute.
  for (auto p : m_nodeset)
    if (ownnode(p)) {
      m_linnodes[ p ] = m_start;
      auto coord = tk::cref_find( m_coordmap, p );
      for (auto c : tk::cref_find(m_nodech,p))
        m_chcoordmap[ c ].emplace( m_start, coord );
      ++m_start;
    }

  // Trigger SDAG wait indicating that reordering own node IDs are complete
  reorderowned_complete();

  // If all our nodes have new IDs assigned, reordering complete on this PE
  if (m_linnodes.size() == m_nodeset.size()) reordered();
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
    std::unordered_map< std::size_t,
      std::tuple< std::size_t, tk::UnsMesh::Coord > > n;
    for (auto p : r.second)
      n.emplace( p, std::make_tuple( tk::cref_find(m_linnodes,p),
                                     tk::cref_find(m_coordmap,p) ) );
    thisProxy[ r.first ].neworder( n );
    tk::destroy( n );
  }

  tk::destroy( m_reqNodes ); // Clear queue of requests just fulfilled

  // Re-enable SDAG wait for preparing new node requests
  thisProxy[ CkMyPe() ].wait4prep();

  // Re-enable trigger signaling that reordering of owned node IDs are
  // complete right away
  reorderowned_complete();
}

void
Partitioner::bndEdges()
// *****************************************************************************
// Generate boundary edges and send them to all PEs
//! \details This step happens when the mesh chunk on this PE has been
//!   distributed after partitioning during an initial mesh refinement step. At
//!   this point we have a contiguous chunk of the mesh on this PE as
//!   determined by the partitioner. The next step is to extract the edges on
//!   the boundary only. The boundary edges (shared by multiple PEs) will be
//!   agreed on a refinement that yields a conforming mesh across PE boundaries.
// *****************************************************************************
{
  Assert( !m_ginpoel.empty(), "No elements are assigned to PE" );

  // Compute local data from global mesh connectivity (m_inpoel, m_gid, m_lid)
  m_el = tk::global2local( m_ginpoel );

  // Generate boundary edges of our mesh chunk
  tk::UnsMesh::EdgeSet bnded;
  auto esup = tk::genEsup( m_inpoel, 4 );         // elements surrounding points
  auto esuel = tk::genEsuelTet( m_inpoel, esup ); // elems surrounding elements
  for (std::size_t e=0; e<esuel.size()/4; ++e) {
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f) {
      if (esuel[mark+f] == -1) {
        auto A = m_gid[ m_inpoel[ mark+tk::lpofa[f][0] ] ];
        auto B = m_gid[ m_inpoel[ mark+tk::lpofa[f][1] ] ];
        auto C = m_gid[ m_inpoel[ mark+tk::lpofa[f][2] ] ];
        bnded.insert( {{{A,B}}} );
        bnded.insert( {{{B,C}}} );
        bnded.insert( {{{C,A}}} );
        Assert( m_lid.find( A ) != end(m_lid), "Local node ID not found" );
        Assert( m_lid.find( B ) != end(m_lid), "Local node ID not found" );
        Assert( m_lid.find( C ) != end(m_lid), "Local node ID not found" );
      }
    }
  }

  // Export boundary edges to all PEs
  thisProxy.addBndEdges( CkMyPe(), bnded );
}

void
Partitioner::addBndEdges( int frompe, const tk::UnsMesh::EdgeSet& ed )
// *****************************************************************************
//! Receive boundary edges from all PEs (including this one)
//! \param[in] frompe PE call coming from
//! \param[in] ed Edges on frompe's boundary (with global node IDs)
// *****************************************************************************
{
  // Store incoming boundary edges
  m_bndEdges[ frompe ].insert( begin(ed), end(ed) );

  if (++m_nedge == static_cast<std::size_t>(CkNumPes())) {
    // Compute unique set of PEs that share at least a single edge with this PE
    const auto& ownedges = tk::cref_find( m_bndEdges, CkMyPe() );
    for (const auto& p : m_bndEdges)    // for all PEs
      if (p.first != CkMyPe())          // for all PEs other than this PE
        for (const auto& e : p.second)  // for all boundary edges
          if (ownedges.find(e) != end(ownedges))
            m_pe.insert( p.first );     // if edge is shared, store its PE

    refine();
  }
}

void
Partitioner::refine()
// *****************************************************************************
//  Do a single step of initial mesh refinement based on user-input
//! \details This is a single step in a potentially multiple-entry list of
//!   initial adaptive mesh refinement steps. Distribution of the PE-boundary
//!   edges has preceded this step, so that boundary edges (shared by multiple
//!   PEs) can agree on a refinement that yields a conforming mesh across PE
//!   boundaries.
// *****************************************************************************
{
  // Convert node coordinates associated to global node IDs to a flat vector
  auto npoin = m_coordmap.size();
  Assert( m_gid.size() == npoin, "Size mismatch" );
  m_coord[0].resize( npoin );
  m_coord[1].resize( npoin );
  m_coord[2].resize( npoin );
  for (const auto& c : m_coordmap) {
    auto i = tk::cref_find( m_lid, c.first );
    Assert( i < npoin, "Indexing out of coordinate map" );
    m_coord[0][i] = c.second[0];
    m_coord[1][i] = c.second[1];
    m_coord[2][i] = c.second[2];
  }

  // Query user input for initial mesh refinement type list
  auto initref = g_inputdeck.get< tag::amr, tag::init >();
  // Determine which level this is
  auto level = initref.size() - m_initref.size();

  // Output mesh before this initial refinement step
  tk::UnsMesh refmesh( m_inpoel, m_coord );
  tk::ExodusIIMeshWriter mw( "initref.b." + std::to_string(level) + '.' +
                               std::to_string(CkMyPe()),
                             tk::ExoWriter::CREATE );
  mw.writeMesh( refmesh );

  for (const auto& e : tk::cref_find(m_bndEdges,CkMyPe())) {
    IGNORE(e);
    Assert( m_lid.find( e[0] ) != end( m_lid ) &&
            m_lid.find( e[1] ) != end( m_lid ),
            "Boundary edge not found before refinement" );
  }

  // Refine mesh based on next initial refinement type
  if (!m_initref.empty()) {
    auto r = m_initref.back();    // consume (reversed) list from back
    if (r == ctr::AMRInitialType::UNIFORM)
      uniformRefine();
    else if (r == ctr::AMRInitialType::INITIAL_CONDITIONS)
      errorRefine();
    else Throw( "Initial AMR type not implemented" );
  }

  // Additionally refine mesh based on user explicitly tagging edges
  userRefine();

  for (const auto& e : tk::cref_find(m_bndEdges,CkMyPe())) {
    IGNORE(e);
    Assert( m_lid.find( e[0] ) != end( m_lid ) &&
            m_lid.find( e[1] ) != end( m_lid ),
            "Boundary edge not found after refinement" );
  }

  // Ensure valid mesh after refinement
  Assert( tk::positiveJacobians( m_inpoel, m_coord ),
          "Refined mesh cell Jacobian non-positive" );

  // Export added nodes on our mesh chunk boundary to other PEs
  if (m_pe.empty())
    contribute( sizeof(std::size_t), &m_extra, CkReduction::max_int,
                m_cb.get< tag::matched >() );
  else {
    m_nref = 0;
    for (auto p : m_pe) {       // for all PEs we share at least an edge with
      // For all boundary edges of PE p, find out if we have added a new
      // node to it, and if so, export parents->(newid,coords) to p.
      tk::UnsMesh::EdgeNodeCoord exp;
      for (const auto& e : tk::cref_find(m_bndEdges,p)) {
        auto i = m_edgenode.find(e);
        if (i != end(m_edgenode)) exp[ e ] = i->second;
      }
      thisProxy[ p ].addRefBndEdges( CkMyPe(), exp );
    }
  }
}

void
Partitioner::addRefBndEdges( int frompe, const tk::UnsMesh::EdgeNodeCoord& ed )
// *****************************************************************************
//! Receive newly added mesh node IDs on our PE boundary
//! \param[in] frompe PE call coming from
//! \param[in] ed Newly added node IDs associated to parent nodes on PE boundary
//! \details Receive newly added global node IDs and coordinates associated to
//!   global parent IDs of edges on our mesh chunk boundary.
// *****************************************************************************
{
  // Save/augment buffer of edge-node (IDs and coords) categorized by sender PE
  m_edgenodePe[ frompe ].insert( begin(ed), end(ed) );
  // Acknowledge receipt of PE-boundary edges to sender
  thisProxy[ frompe ].recvRefBndEdges();
}

void
Partitioner::recvRefBndEdges()
// *****************************************************************************
//  Acknowledge received newly added node IDs to edges shared among multiple PEs
// *****************************************************************************
{
  // When we have heard from all PEs we share at least a single edge with,
  // contribute the number of extra edges that this mesh refinement step has
  // found that were not refined by this PE but were refined by other PEs this
  // PE shares the edge with. A global maximum will then be computed on the
  // number of extra edges appearing in Transporter::matched() and that is then
  // used to decide if a new correction step is needed. If this is called for
  // the first time in a given initial mesh refinement step, i.e., not after a
  // correction step, m_extra = 1 on all PEs, so a correction step is assumed
  // to be required.
  if (++m_nref == m_pe.size()) {
    contribute( sizeof(std::size_t), &m_extra, CkReduction::max_int,
                m_cb.get< tag::matched >() );
  }
}

void
Partitioner::correctref()
// *****************************************************************************
//  Correct refinement to arrive at a conforming mesh across PE boundaries
//! \details This function is called repeatedly until there is not a a single
//!    edge that needs correction for the whole distributed problem to arrive at
//!    a conforming mesh across PE boundaries during this initial mesh
//!    refinement step.
// *****************************************************************************
{
  // Storage for edges that still need a new node to yield a conforming mesh
  tk::UnsMesh::EdgeSet extra;

  // Ensure that the same global node ID has been assigned by all PEs and that
  // the new nodes have the same coordinates generated by potentially multiple
  // PEs sharing the refined edge. This is done by searching for all edges that
  // we share with and refined by other PEs: (1) If the incoming edge is found
  // among our refined ones, we ensure the newly assigned global IDs equal
  // (independently assigned by multiple PEs) and also that the new node
  // coordinates equal to machine precision. (2) If the incoming edge is not
  // found among our refined ones, we need to correct the mesh to make it
  // conforming since the edge has been refined by the remote PE. We collect
  // these extra edges, and run a correction refinement, whose result then
  // needs to be communicated again as the new refinement step may introduce
  // new edges that other PEs did not refine but are shared.
  for (const auto& p : m_edgenodePe)        // for all PEs we share edges with
    for (const auto& e : p.second) {        // for all refined edges on p.first
      auto i = m_edgenode.find( e.first );  // find refined edge given parents
      if (i != end(m_edgenode)) {           // found same added node on edge
        // locally assigned added node ID and coordinates: i->second
        // remotely assigned added node ID and coordinates: e.second
        // Ensure global IDs of newly added nodes are the same
        Assert( std::get< 0 >( i->second ) == std::get< 0 >( e.second ),
                "Remotely and locally assigned global ids mismatch" );
        // Ensure coordinates are the same
        Assert( std::abs( std::get<1>(i->second) - std::get<1>(e.second) ) <
                  std::numeric_limits<tk::real>::epsilon() &&
                std::abs( std::get<2>(i->second) - std::get<2>(e.second) ) <
                  std::numeric_limits<tk::real>::epsilon() &&
                std::abs( std::get<3>(i->second) - std::get<3>(e.second) ) <
                  std::numeric_limits<tk::real>::epsilon(),
                "Remote and local added node coordinates mismatch" );
      } else {  // remote PE added node on edge but we did not
        // Make sure we know about this boundary-PE edge (that we did not refine)
        Assert( m_bndEdges.find( CkMyPe() )->second.find( e.first ) !=
                m_bndEdges.find( CkMyPe() )->second.end(),
                "Local node IDs of boundary edge not found" );
        // Save edge (given by parent global node IDs) to which the remote PE
        // has added a new node but we did not. Will need to correct the mesh so
        // it conforms across PEs.
        extra.insert( {{ { tk::cref_find( m_lid, e.first[0] ),
                           tk::cref_find( m_lid, e.first[1] ) } }} );
      }
    }

  // Store number of extra edges on this PE which this PE did not add but was
  // refined by another PE, so now we need to tag and refine them and propagate
  // reconnection of neighbor cells to arrive at conforming mesh across PE
  // boundaries.
  m_extra = extra.size();

  // Refine mesh triggered by nodes added on PE-boundary edges by other PEs
  // PEs
  correctRefine( extra );

  // Since refining edges that we originally did not but other PEs did may
  // result in refining new edges that may be shared along PE boundaries, we
  // now need to communicate these edges and potentially repeat the correction
  // step. This must happen until all PEs that share edges can agree that there
  // are no more edges to correct. Only then this refinement step can be
  // considered complete.
  if (m_pe.empty())
    contribute( sizeof(std::size_t), &m_extra, CkReduction::max_int,
                m_cb.get< tag::matched >() );
  else {
    m_nref = 0;
    for (auto p : m_pe) {       // for all PEs we share at least an edge with
      // For all boundary edges of PE p, find out if we have added a new
      // node to it, and if so, export parents->(newid,coords) to p.
      tk::UnsMesh::EdgeNodeCoord exp;
      for (const auto& e : tk::cref_find(m_bndEdges,p)) {
        auto i = m_edgenode.find(e);
        if (i != end(m_edgenode)) exp[ e ] = i->second;
      }
      thisProxy[ p ].addRefBndEdges( CkMyPe(), exp );
    }
  }
}

void
Partitioner::nextref()
// *****************************************************************************
// Decide wether to continue with another step of initial mesh refinement
//! \details At this point the mesh has been refined and all PEs have received
//!   a map associating the global IDs and the coordinates of a node added to
//!   an edge during initial mesh refinement associated to all other PEs the
//!   edges are shared with. Now the mesh is corrected so that it conforms
//!   across PE-boundaries by tagging those edges for refinement that have been
//!   refined by at least a PE. This concludes this initial mesh refinement
//!   step, and we continue if there are more steps configured by the user.
// *****************************************************************************
{
  // Remove initial mesh refinement step from list
  if (!m_initref.empty()) m_initref.pop_back();

  if (!m_initref.empty())       // Continue to next initial refinement step
    partref();
  else {                        // Finish list of initial mesh refinement steps
    // Output final mesh after initial mesh refinement
    tk::UnsMesh refmesh( m_inpoel, m_coord );
    tk::ExodusIIMeshWriter mw( "initref.final." + std::to_string(CkMyPe()),
                               tk::ExoWriter::CREATE );
    mw.writeMesh( refmesh );
    // Finish initial mesh refinement
    finishref();
  }
}

void
Partitioner::finishref()
// *****************************************************************************
// Finish initial mesh refinement
//! \details This function is called as after initial mesh refinement has
//!   finished. If initial mesh reifnement was not configured by the user, this
//!   is the point where we continue after the constructor, by computing the
//!   total number of elements across the whole problem.
// *****************************************************************************
{
  // Compute final number of cells across whole problem
  auto nelem = m_ginpoel.size()/4;
  contribute( sizeof(uint64_t), &nelem, CkReduction::sum_int,
              m_cb.get< tag::refined >() );
}

void
Partitioner::uniformRefine()
// *****************************************************************************
// Do uniform mesh refinement
// *****************************************************************************
{
  // Instantiate mesh refiner
  AMR::mesh_adapter_t refiner( m_inpoel );

  // Do uniform refinement
  refiner.uniform_refinement();

  // Update mesh coordinates and connectivity
  updateMesh( refiner );
}

void
Partitioner::errorRefine()
// *****************************************************************************
// Do error-based mesh refinement
// *****************************************************************************
{
  // Instantiate mesh refiner
  AMR::mesh_adapter_t refiner( m_inpoel );

  // Find number of nodes in old mesh
  auto npoin = tk::npoin( m_inpoel );
  // Generate edges surrounding points in old mesh
  auto esup = tk::genEsup( m_inpoel, 4 );
  auto psup = tk::genPsup( m_inpoel, 4, esup );

  // Evaluate initial conditions at mesh nodes
  auto u = nodeinit( npoin, esup );

  // Get the indices (in the system of systems) of refinement variables and the
  // error indicator configured
  const auto& refidx = g_inputdeck.get< tag::amr, tag::id >();
  auto errtype = g_inputdeck.get< tag::amr, tag::error >();

  // Compute errors in ICs and define refinement criteria for edges
  std::vector< edge_t > edge;
  std::vector< real_t > crit;
  AMR::Error error;
  for (std::size_t p=0; p<npoin; ++p)   // for all mesh nodes on this PE
    for (auto q : tk::Around(psup,p)) { // for all nodes surrounding p
       tk::real cmax = 0.0;
       edge_t e(p,q);
       for (auto i : refidx) {          // for all refinement variables
         auto c = error.scalar( u, e, i, m_coord, m_inpoel, esup, errtype );
         if (c > cmax) cmax = c;        // find max error at edge
       }
       if (cmax > 0.0) {         // if nonzero error, will pass edge to refiner
         edge.push_back( e );
         crit.push_back( cmax );
       }
     }

  Assert( edge.size() == crit.size(), "Size mismatch" );

  // Do error-based refinement
  refiner.error_refinement( edge, crit );

  // Update mesh coordinates and connectivity
  updateMesh( refiner );
}

void
Partitioner::userRefine()
// *****************************************************************************
// Do mesh refinement based on user explicitly tagging edges
// *****************************************************************************
{
  // Instantiate mesh refiner
  AMR::mesh_adapter_t refiner( m_inpoel );

  // Find number of nodes in old mesh
  auto npoin = tk::npoin( m_inpoel );
  // Generate edges surrounding points in old mesh
  auto esup = tk::genEsup( m_inpoel, 4 );
  auto psup = tk::genPsup( m_inpoel, 4, esup );

  // Get user-defined node-pairs (edges) to tag for refinement
  const auto& edgenodelist = g_inputdeck.get< tag::amr, tag::edge >();
  tk::UnsMesh::EdgeSet edgeset;
  for (std::size_t i=0; i<edgenodelist.size()/2; ++i)
    edgeset.insert( {{ edgenodelist[i*2+0], edgenodelist[i*2+1] }} );

  // Compute errors in ICs and define refinement criteria for edges
  std::vector< edge_t > edge;
  std::vector< real_t > crit;
  for (std::size_t p=0; p<npoin; ++p)        // for all mesh nodes on this PE
    for (auto q : tk::Around(psup,p)) {      // for all nodes surrounding p
      tk::UnsMesh::Edge e{{p,q}};
      if (edgeset.find(e) != end(edgeset)) { // tag edge if on user's list
        edge.push_back( edge_t(e[0],e[1]) );
        crit.push_back( 1.0 );
      }
    }

  Assert( edge.size() == crit.size(), "Size mismatch" );

  // Do error-based refinement
  refiner.error_refinement( edge, crit );

  // Update mesh coordinates and connectivity
  updateMesh( refiner );
}

tk::Fields
Partitioner::nodeinit( std::size_t npoin,
                       const std::pair< std::vector< std::size_t >,
                          std::vector< std::size_t > >& esup )
// *****************************************************************************
// Evaluate initial conditions (IC) at mesh nodes
//! \param[in] npoin Number points in mesh (on this PE)
//! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
//! \return Initial conditions (evaluated at t0) at nodes
// *****************************************************************************
{
  auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();
  auto nprop = g_inputdeck.get< tag::component >().nprop();

  // Will store nodal ICs
  tk::Fields u( m_coord[0].size(), nprop );

  // Evaluate ICs differently depending on nodal or cell-centered discretization
  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
  if (scheme == ctr::SchemeType::MatCG || scheme == ctr::SchemeType::DiagCG) {

    // Node-centered: evaluate ICs for all scalar components integrated
    for (const auto& eq : g_cgpde) eq.initialize( m_coord, u, t0 );

  } else if (scheme == ctr::SchemeType::DG) {

    // Initialize cell-centered unknowns
    tk::Fields ue( m_inpoel.size()/4, nprop );
    auto geoElem = tk::genGeoElemTet( m_inpoel, m_coord );
    for (const auto& eq : g_dgpde) eq.initialize( geoElem, ue, t0 );

    // Transfer initial conditions from cells to nodes
    for (std::size_t p=0; p<npoin; ++p) {       // for all mesh nodes on this PE
      std::vector< tk::real > up( nprop, 0.0 );
      tk::real vol = 0.0;
      for (auto e : tk::Around(esup,p)) {       // for all cells around node p
        // compute nodal volume: every element contributes their volume / 4
        vol += geoElem(e,0,0) / 4.0;
        // sum cell value to node weighed by cell volume / 4
        for (std::size_t c=0; c<nprop; ++c)
          up[c] += ue[e][c] * geoElem(e,0,0) / 4.0;
      }
      // store nodal value
      for (std::size_t c=0; c<nprop; ++c) u(p,c,0) = up[c] / vol;
    }

  } else Throw( "Nodal initialization not handled for discretization scheme" );

  Assert( u.nunk() == m_coord[0].size(), "Size mismatch" );
  Assert( u.nprop() == nprop, "Size mismatch" );

  return u;
}

void
Partitioner::correctRefine( const tk::UnsMesh::EdgeSet& extra )
// *****************************************************************************
// Do mesh refinement correcting PE-boundary edges
//! \param[in] extra Unique set of edges that need a new node on PE boundaries
// *****************************************************************************
{
  if (!extra.empty()) {
    // Instantiate mesh refiner
    AMR::mesh_adapter_t refiner( m_inpoel );
  
    // Generate list of edges that need to be corrected
    std::vector< edge_t > edge;
    for (const auto& e : extra) edge.push_back( edge_t(e[0],e[1]) );
    std::vector< real_t > crit( edge.size(), 1.0 );
  
    // Do refinement including edges that need to be corrected
    refiner.error_refinement( edge, crit );
  
    // Update mesh coordinates and connectivity
    updateMesh( refiner );
  }
}

void
Partitioner::updateMesh( AMR::mesh_adapter_t& refiner )
// *****************************************************************************
// Update mesh after refinement
//! \param[in] refiner Mesh refiner (AMR) object
// *****************************************************************************
{
  // Get refined mesh connectivity
  const auto& refinpoel = refiner.tet_store.get_active_inpoel();
  Assert( refinpoel.size()%4 == 0, "Inconsistent refined mesh connectivity" );

  // Generate unique node lists of old and refined mesh using local ids
  std::unordered_set< std::size_t > old( m_inpoel.cbegin(), m_inpoel.cend() );
  std::unordered_set< std::size_t > ref( refinpoel.cbegin(), refinpoel.cend() );

  updateVolumeMesh( refiner, old, ref );

  updateBoundaryMesh( refiner, old, ref );

  // Update mesh connectivity with local node IDs
  m_inpoel = refinpoel;

  // Update mesh connectivity with new global node ids
  m_ginpoel = m_inpoel;
  Assert( tk::uniquecopy(m_ginpoel).size() == m_coord[0].size(),
          "Size mismatch" );
  for (auto& i : m_ginpoel) i = m_gid[i];
}

void
Partitioner::updateVolumeMesh( AMR::mesh_adapter_t& refiner,
                               const std::unordered_set< std::size_t >& old,
                               const std::unordered_set< std::size_t >& ref )
// *****************************************************************************
//  Update volume mesh after mesh refinement
//! \param[in] refiner Mesh refiner (AMR) object
//! \param[in] old Unique nodes of the old (unrefined) mesh using local ids
//! \param[in] ref Unique nodes of the refined mesh using local ids
// *****************************************************************************
{
  auto& x = m_coord[0];
  auto& y = m_coord[1];
  auto& z = m_coord[2];

  // Resize node coordinate storage to accommodate refined mesh nodes
  auto npoin = ref.size();
  x.resize( npoin );
  y.resize( npoin );
  z.resize( npoin );
  m_gid.resize( npoin, std::numeric_limits< std::size_t >::max() );

  // Generate coordinates and ids to newly added nodes after refinement step
  for (auto r : ref) {               // for all unique nodes of the refined mesh
    if (old.find(r) == end(old)) {   // if node is newly added (in this step)
      // get (local) parent ids of newly added node
      auto p = refiner.node_connectivity.get( r );
      Assert( old.find(p[0]) != end(old) && old.find(p[1]) != end(old),
              "Parent(s) not in old mesh" );
      Assert( r >= old.size(), "Attempting to overwrite node with added one" );
      // generate coordinates for newly added node
      x[r] = (x[p[0]] + x[p[1]])/2.0;
      y[r] = (y[p[0]] + y[p[1]])/2.0;
      z[r] = (z[p[0]] + z[p[1]])/2.0;
      decltype(p) gp{{ m_gid[p[0]], m_gid[p[1]] }}; // global parent ids
      // generate new global ID for newly added node
      auto g = tk::UnsMesh::EdgeHash()( gp );
      // ensure newly generated node id has not yet been used
      Assert( g >= old.size(), "Hashed id overwriting old id" );
      Assert( m_coordmap.find(g) == end(m_coordmap),
              "Hash collision: ID already exist" );
      // assign new global ids to local->global and to global->local maps
      m_gid[r] = g;
      Assert( m_lid.find(g) == end(m_lid),
              "Overwriting entry global->local node ID map" );
      m_lid[g] = r;
      // assign new coordinates to new global node id
      Assert( m_coordmap.find(g) == end(m_coordmap),
              "Overwriting entry coordmap" );
      m_coordmap.insert( {g, {{x[r], y[r], z[r]}}} );
      // assign new coordinates and new global node id to global parent id pair
      m_edgenode[ gp ] = std::make_tuple( g, x[r], y[r], z[r] );
    }
  }

  Assert( m_gid.size() == m_lid.size(), "Size mismatch" );

  Assert( std::none_of( begin(m_gid), end(m_gid), [](std::size_t i){
            return i == std::numeric_limits< std::size_t >::max(); } ),
          "Not all local->global node IDs have been assigned" );
}

void
Partitioner::updateBoundaryMesh( AMR::mesh_adapter_t& refiner,
                                 const std::unordered_set< std::size_t >& old,
                                 const std::unordered_set< std::size_t >& ref )
// *****************************************************************************
// Update boundary data structures after mesh refinement
//! \param[in] refiner Mesh refiner (AMR) object
//! \param[in] old Unique nodes of the old (unrefined) mesh using local ids
//! \param[in] ref Unique nodes of the refined mesh using local ids
// *****************************************************************************
{
  IGNORE(ref);  // to avoid compiler warning when asserts are optimized away

  using Face = tk::UnsMesh::Face;

  // Will associate a pair of side set id and adjacent tet id to a boundary
  // triangle face
  using BndFaces = std::unordered_map< Face,
                                       std::pair< int, std::size_t >,
                                       tk::UnsMesh::FaceHasher,
                                       tk::UnsMesh::FaceEq >;

  // Generate data structure that associates the pair of side set id and
  // adjacent tet id to a boundary triangle face for all boundary faces. After
  // this loop we will have all tets adjacent to boundary faces, where the
  // "boundary" includes all physical boundaries (where the user may assign
  // boundary conditions via side sets from the mesh file) as well as
  // PE-boundaries due to domain decomposition.
  BndFaces bnd;
  auto oldesuel = tk::genEsuelTet( m_inpoel, tk::genEsup(m_inpoel,4) );
  for (std::size_t e=0; e<oldesuel.size()/4; ++e) {     // for tets on this PE
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f) {
      if (oldesuel[mark+f] == -1) {  // if a face does not have an adjacent  tet
        Face t{{ m_gid[ m_inpoel[ mark+tk::lpofa[f][0] ] ],
                 m_gid[ m_inpoel[ mark+tk::lpofa[f][1] ] ],
                 m_gid[ m_inpoel[ mark+tk::lpofa[f][2] ] ] }};
        // associate tet id adjacent to boundary face to boundary face, at this
        // point we fill in -1 for the side set id, since we don't know it yet
        bnd[t] = {-1,e};
      }
    }
  }

  // Assign side set ids to faces on the physical boundary only
  for (const auto& ss : m_bface)  // for all phsyical boundaries (sidesets)
    for (auto f : ss.second) {    // for all faces on this physical boundary
      Face t{{ m_triinpoel[f*3+0], m_triinpoel[f*3+1], m_triinpoel[f*3+2] }};
      auto cf = bnd.find( t );
      if (cf != end(bnd)) cf->second.first = ss.first;
    }

  // Remove PE-boundary faces (to which the above loop did not assign a set id)
  auto bndcopy = bnd;
  for (const auto& f : bndcopy) if (f.second.first == -1) bnd.erase( f.first );
  tk::destroy( bndcopy );

  // Now in bnd we have a pair of side set ids of all tet ids that are adjacent
  // to all physical boundary faces that are associated to a side set in the
  // mesh file and only along faces on the physical boundary and not along faces
  // on the PE-boundary.

  // storage for boundary faces associated to side-set IDs of the refined mesh
  decltype(m_bface) bface;              // will become m_bface
  // storage for boundary faces-node connectivity of the refined mesh
  decltype(m_triinpoel) triinpoel;      // will become m_triinpoel
  // face id counter
  std::size_t facecnt = 0;

  // Lambda to associate a boundary face and connectivity to a side set.
  // Parameter 's' is the list of faces (ids) to add new face to. Parameter 'f'
  // is the triangle face connecttivity to add.
  auto addBndFace = [ &facecnt, &triinpoel ]
                    ( std::vector< std::size_t >& s, const Face& f )
  {
    s.push_back( facecnt++ );
    triinpoel.push_back( f[0] );
    triinpoel.push_back( f[1] );
    triinpoel.push_back( f[2] );
  };

  // Lambda to find the nodes of the parent face of a child face. Parameter
  // 'face' denotes the node ids of the child face whose parent face we are
  // looking for. This search may find 3 or 4 parent nodes, depending on whether
  // the child shares or does not share a face with the parent tet,
  // respectively.
  auto parentFace = [ &old, &refiner ]( const Face& face ){
    std::unordered_set< std::size_t > s;// will store nodes of parent face
    for (auto n : face) {               // for all 3 nodes of the face
      if (old.find(n) != end(old))      // if child node found in old mesh,
        s.insert( n );                  // that node is also in the parent face
      else {                            // if child node is a newly added one
        auto p = refiner.node_connectivity.get( n );  // find its parent nodes
        s.insert( begin(p), end(p) );                 // and store both uniquely
      }
    }
    // Ensure all parent nodes are part of the old (unrefined) mesh
    Assert( std::all_of( begin(s), end(s), [ &old ]( std::size_t j ){
                           return old.find(j) != end(old); } ),
            "Old mesh nodes not in old mesh" );
    // Return unique set of nodes of the parent face of child face
    return s;
  };

  // Generate boundary face data structures after refinement step
  for (const auto& f : bnd) {
    // construct face of old mesh boundary face using local ids
    Face oldface{{ tk::cref_find( m_lid, f.first[0] ),
                   tk::cref_find( m_lid, f.first[1] ),
                   tk::cref_find( m_lid, f.first[2] ) }};
    // will associate to side set id of old (unrefined) mesh boundary face
    auto& side = bface[ f.second.first ];
    // query number of children of boundary tet adjacent to boundary face
    auto nc = refiner.tet_store.data( f.second.second ).num_children;
    if (nc == 0) {      // if boundary tet is not refined, add its boundary face
      addBndFace( side, f.first );
    } else {            // if boundary tet is refined
      const auto& tets = refiner.tet_store.tets;
      for (decltype(nc) i=0; i<nc; ++i ) {      // for all child tets
        // get child tet id
        auto childtet = refiner.tet_store.get_child_id( f.second.second, i );
        auto ct = tets.find( childtet );
        Assert( ct != end(tets), "Child tet not found" );
        // ensure all nodes of child tet are in refined mesh
        Assert( ref.find(ct->second[0]) != end(ref) &&
                ref.find(ct->second[1]) != end(ref) &&
                ref.find(ct->second[2]) != end(ref) &&
                ref.find(ct->second[3]) != end(ref),
                "Boundary child tet node id not found in refined mesh" );
        // get nodes of child tet
        auto A = ct->second[0];
        auto B = ct->second[1];
        auto C = ct->second[2];
        auto D = ct->second[3];
        // form all 4 faces of child tet
        std::array<Face, 4> face{{{{A,C,B}}, {{A,B,D}}, {{A,D,C}}, {{B,C,D}}}};
        for (const auto& rf : face) {   // for all faces of child tet
          // find nodes of the parent face of child face
          auto parfac = parentFace( rf );
          // if the child shares a face with its parent and all 3 nodes of the
          // parent face of the child's face are the same, the child face is on
          // the same boundary as the parent face
          if ( parfac.size() == 3 &&
               parfac.find(oldface[0]) != end(parfac) &&
               parfac.find(oldface[1]) != end(parfac) &&
               parfac.find(oldface[2]) != end(parfac) )
          {
            addBndFace( side, {m_gid[rf[0]], m_gid[rf[1]], m_gid[rf[2]]} );
          }
        }
      }
    }
  }

  // Update boundary face data structures
  m_bface = std::move(bface);
  m_triinpoel = std::move(triinpoel);
}

void
Partitioner::reordered()
// *****************************************************************************
//  Compute final result of reordering
//! \details At this point the node coordinates on all PEs have been updated to
//!   be consistent with the new ordering. We continue by updating other data,
//!   such as mesh connectivities of chares.
// *****************************************************************************
{
  tk::destroy( m_bnodechares );
  tk::destroy( m_ncommunication );

  // Construct maps associating old node IDs (as in file) to new node IDs
  // (as in producing contiguous-row-id linear system contributions)
  // associated to chare IDs (outer key).
  for (const auto& c : m_chinpoel) {
    auto& nodes = m_chfilenodes[ c.first ];
    for (auto p : c.second)
      nodes[ tk::cref_find(m_linnodes,p) ] = p;
  }

  // Update chare-categorized elem connectivities with the reordered node IDs
  for (auto& c : m_chinpoel)
    for (auto& p : c.second)
       p = tk::cref_find( m_linnodes, p );

  // Update chare-categorized chare-mesh-nodes comm map with the reordered IDs
  for (auto& c : m_msum)
    for (auto& s : c.second) {
      decltype(s.second) n;
      for (auto p : s.second) {
        n.insert( tk::cref_find( m_linnodes, p ) );
      }
      s.second = std::move( n );
    }

  // Update unique global node IDs chares on our PE will contribute to with
  // the reordered node IDs
  m_nodeset.clear();
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
  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
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
  auto dist = distribution( m_nchare );

  for (int c=0; c<dist[1]; ++c) {
    // Compute chare ID
    auto cid = CkMyPe() * dist[0] + c;
    // Guard those searches that operate on empty containers in serial
    typename decltype(m_msum)::mapped_type msum;
    if (!m_msum.empty()) msum = tk::cref_find( m_msum, cid );
    // Create worker array element using Charm++ dynamic chare array element
    // insertion: 1st arg: chare id, last arg: PE chare is created on, middle
    // args: Discretization ctor args. See also Charm++ manual, Sec. "Dynamic
    // Insertion".
    m_scheme.discInsert( cid, m_host, tk::cref_find(m_chinpoel,cid),
      tk::cref_find(m_chcoordmap,cid), msum, m_nchare, CkMyPe() );
  }

  // Free storage for unique global mesh nodes chares on our PE will
  // contribute to in a linear system as no longer needed.
  tk::destroy( m_nodeset );
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
  auto dist = distribution( m_nchare );

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
  
    // Generate boundary face and node ids (after mesh node reordering)
    std::vector< std::size_t > chtriinpoel;
    std::unordered_map< int, std::vector< std::size_t > > chbface;
    std::map< int, std::vector< std::size_t > > chbnode;
    std::size_t cnt = 0;

    for (const auto& ss : m_bface)  // for all phsyical boundaries (sidesets)
      for (auto f : ss.second) {    // for all faces on this physical boundary
        // attempt to find face nodes on this chare
        auto f1 = newnodes.find( m_triinpoel[f*3+0] );
        auto f2 = newnodes.find( m_triinpoel[f*3+1] );
        auto f3 = newnodes.find( m_triinpoel[f*3+2] );
        // if face node on this chare, assign its new (global) id to side set
        auto& n = chbnode[ ss.first ];
        if (f1 != end(newnodes)) n.push_back( f1->second );
        if (f2 != end(newnodes)) n.push_back( f2->second );
        if (f3 != end(newnodes)) n.push_back( f3->second );
        // if all 3 nodes of the physical boundary face are on this chare
        if (f1 != end(newnodes) && f2 != end(newnodes) && f3 != end(newnodes)) {
          // Create face with new node ids (after mesh node reordering)
          std::array< std::size_t, 3 > t{{f1->second, f2->second, f3->second}};
          // if this boundary face is on this chare
          if (faceset.find(t) != end(faceset)) {
            // store face connectivity with new (global) node ids of this chare
            chtriinpoel.insert( end(chtriinpoel), begin(t), end(t) );
            // generate/store physical boundary face id associated to sideset id
            chbface[ ss.first ].push_back( cnt++ );
          }
        }
      }

    // Make boundary node IDs unique for each physical boundary (side set)
    for (auto& s : chbnode) tk::unique( s.second );

    // Face data class
    FaceData fd( chinpoel, chbface, chbnode, chtriinpoel );

    // Make sure (bound) base is already created and accessible
    Assert( m_scheme.get()[cid].ckLocal() != nullptr, "About to pass nullptr" );

    // Create worker array element using Charm++ dynamic chare array element
    // insertion: 1st arg: chare id, last arg: PE chare is created on, middle
    // args: Discretization's child ctor args. See also Charm++ manual, Sec.
    // "Dynamic Insertion".
    m_scheme.insert( cid, m_scheme.get(), m_solver, fd, CkMyPe() );
  }

  tk::destroy( m_bface );
  tk::destroy( m_triinpoel );
}

#include "NoWarning/partitioner.def.h"
