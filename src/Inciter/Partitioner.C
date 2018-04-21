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
#include "AMR/Error.h"
#include "Inciter/Options/Scheme.h"
#include "Inciter/Options/AMRInitial.h"
#include "UnsMesh.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< CGPDE > g_cgpde;

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
  m_npeDist( 0 ),
  m_npe( 0 ),
  m_ndist( 0 ),
  m_nedge( 0 ),
  m_nref( 0 ),
  m_pe(),
  m_initref(),
  m_el(),
  m_edgenode(),
  m_bndEdges(),
  m_edgeNodeCoord(),
  m_reqNodes(),
  m_start( 0 ),
  m_noffset( 0 ),
  m_nquery( 0 ),
  m_nmask( 0 ),
  m_ginpoel(),
  m_rinpoel(),
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

  if (!g_inputdeck.get< tag::amr, tag::init >().empty())
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
  // Move mesh connectivity to under a new name. This also clears the original
  // one. The new one (rinpoel) is treated during communication as a source, and
  // the old one (ginpoel) will serve as the destination for the new one.
  m_rinpoel = std::move( m_ginpoel );

  // Generate element IDs for Zoltan
  std::vector< long > gelemid( m_rinpoel.size()/4 );
  std::iota( begin(gelemid), end(gelemid), 0 );

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
  m_bndEdges.clear();
  m_pe.clear();
  m_edgeNodeCoord.clear();
  m_edgenode.clear();
  m_coordmap.clear();

  // Categorize mesh elements (given by their gobal node IDs) by target PE and
  // distribute to their PEs based on mesh partitioning.
  distributePE( categorize( pel, m_rinpoel ) );
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
//! \param[in] mesh Map associating mesh connectivities with global node indices
//!   and node coordinates for mesh chunks we are assigned by the partitioner
// *****************************************************************************
{
  // Store mesh connectivity and global node coordinates categorized by chares
  for (const auto& c : chmesh) {
    Assert( pe(c.first) == CkMyPe(), "PE " + std::to_string(CkMyPe()) +
            " received a mesh whose chare it does not own" );
    // The send also writers to this so append/concat
    auto& inpoel = m_chinpoel[ c.first ];
    const auto& mesh = std::get< 0 >( c.second );
    inpoel.insert( end(inpoel), begin(mesh), end(mesh) );
    // Store coordinates associated to global node IDs. The send side also
    // writes to this, so concat.
    const auto& coord = std::get< 1 >( c.second );
    Assert( tk::cunique(mesh).size() == coord.size(), "Size mismatch" );
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
  Assert( tk::cunique(inpoel).size() == cm.size(), "Size mismatch" );

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
  Assert( tk::cunique(inpoel).size() == coord[0].size(), "Size mismatch" );

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
// *****************************************************************************
{
  tk::UnsMesh::CoordMap map;

  for (auto g : tk::cunique(inpoel)) {
     auto i = tk::cref_find( m_lid, g );
     auto& c = map[g];
     c[0] = m_coord[0][i];
     c[1] = m_coord[1][i];
     c[2] = m_coord[2][i];
  }

  Assert( tk::cunique(inpoel).size() == map.size(), "Size mismatch" );

  return map;
}

void
Partitioner::distributePE(
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
  } else Throw( "No elements are assigned to PE" );

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
      m_chinpoel.insert( *it );           // move over owned key-value pairs
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
  wait4prep();
  wait4bounds();

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

  wait4prep();      // Re-enable SDAG wait for preparing new node requests

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
        for (const auto& e : p.second)  // for all boundary edge
          if (ownedges.find(e) != end(ownedges))
            m_pe.insert( p.first );

    m_nref = 0;
    refine();
  }
}

void
Partitioner::refine()
// *****************************************************************************
//  Optionally refine mesh
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
    m_coord[0][i] = c.second[0];
    m_coord[1][i] = c.second[1];
    m_coord[2][i] = c.second[2];
  }

  // Query user input for initial mesh refinement type list
  auto initref = g_inputdeck.get< tag::amr, tag::init >();
  Assert( !initref.empty(), "No initial mesh refinement steps configured" );
  // Determine which level this is
  auto level = initref.size() - m_initref.size();

  // Output mesh before this initial refinement step
  tk::UnsMesh refmesh( m_inpoel, m_coord );
  tk::ExodusIIMeshWriter mw( "initref." + std::to_string(level) + '.' +
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
  //std::cout << CkMyPe() << " ref start " << level << '\n';
  auto r = m_initref.back();    // consume (reversed) list from back
  if (r == ctr::AMRInitialType::UNIFORM)
    uniformRefine();
  else if (r == ctr::AMRInitialType::INITIAL_CONDITIONS)
    errorRefine();
  else Throw( "Initial AMR type not implemented" );
  std::cout << CkMyPe() << " ref finish " << level << '\n';

  for (const auto& e : tk::cref_find(m_bndEdges,CkMyPe())) {
    IGNORE(e);
    Assert( m_lid.find( e[0] ) != end( m_lid ) &&
            m_lid.find( e[1] ) != end( m_lid ),
            "Boundary edge not found after refinement" );
  }

  // Ensure valid mesh after refinement
  Assert( tk::positiveJacobians( m_inpoel, m_coord ),
          "Refined mesh Jacobian non-positive" );

  // Export added nodes on our mesh chunk boundary to other PEs
  if (m_pe.empty())
    contribute( m_cb.get< tag::matched >() );
  else
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
  // Save to buffer categorized by sender PE
  m_edgeNodeCoord[ frompe ] = ed;
  // Acknowledge receipt of PE-boundary edges to sender
  thisProxy[ frompe ].recvRefBndEdges();
}

void
Partitioner::recvRefBndEdges()
// *****************************************************************************
//  Acknowledge received newly added node IDs to edges shared among multiple PEs
// *****************************************************************************
{
  if (++m_nref == m_pe.size()) contribute( m_cb.get< tag::matched >() );
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
  // Storage for edges that still need a new node to yield a conforming mesh
  tk::UnsMesh::EdgeSet extra;

  // Ensure that the same global node ID has been assigned by all PEs and that
  // the new nodes have the same coordinates generated by potentially multiple
  // PEs sharing the refined edge. This is done by searching for all edges that
  // we share with and refined by other PEs: (1) If the incoming edge is found
  // among our refined ones, we ensure the newly assigned global IDs equal
  // (independently assigned PEs) and also that the new node coordinates equal
  // to machine precision. (2) If the incoming edge is not found among our
  // refined ones, we need to correct the mesh to make it conforming since the
  // edge has been refined by the remote PE. We collect these extra edges, and
  // run a correction refinement.
  for (const auto& p : m_edgeNodeCoord)     // for all PEs we share edges with
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
        Assert( m_bndEdges.find( CkMyPe() )->second.find( e.first ) !=
                m_bndEdges.find( CkMyPe() )->second.end(),
                "Local node IDs of boundary edge not found" );
        // Save edge (given by parent global node IDs) to which the remote PE
        // has added a new node but we did not. Will need to correct mesh so it
        // conforms across PEs.
        extra.insert( {{ { tk::cref_find( m_lid, e.first[0] ),
                           tk::cref_find( m_lid, e.first[1] ) } }} );
      }
    }

  // Correct PE-boundary edges
  correctRefine( extra );

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

  auto& x = m_coord[0];

  // Set initial conditions for all PDEs (use CG for now)
  auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();
  tk::Fields u( x.size(), g_inputdeck.get< tag::component >().nprop() );

  for (const auto& eq : g_cgpde) eq.initialize( m_coord, u, t0 );

//   // hard-code jump at x=0.5 for now
//   for (std::size_t i=0; i<u.nunk(); ++i)
//     if (x[i] > 0.5) u(i,0,0) = 0.0; else u(i,0,0) = 1.0;

  // Find number of nodes in old mesh
  auto npoin = tk::npoin( m_inpoel );
  // Generate edges surrounding points in old mesh
  auto esup = tk::genEsup( m_inpoel, 4 );
  auto psup = tk::genPsup( m_inpoel, 4, esup );

  std::vector< edge_t > edge;
  std::vector< real_t > crit;

  // Compute errors in initial condition and define refinement critera for edges
  AMR::Error error;
  for (std::size_t p=0; p<npoin; ++p)
    for (auto q : tk::Around(psup,p)) {
       edge_t e(p,q);
       edge.push_back( e );
       auto c = error.scalar( u, e, 0, m_coord, m_inpoel, esup,
                              g_inputdeck.get< tag::amr, tag::error >() );
       crit.push_back( c );
     }

  // Do error-based refinement
  refiner.error_refinement( edge, crit );

  // Update mesh coordinates and connectivity
  updateMesh( refiner );
}

void
Partitioner::correctRefine( const tk::UnsMesh::EdgeSet& extra )
// *****************************************************************************
// Do mesh refinement correcting PE-boundary edges
//! \param[in] Unique set of edges that need a new edge
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
// *****************************************************************************
{
  // Get refined mesh connectivity
  const auto& refinpoel = refiner.tet_store.get_active_inpoel();
  Assert( refinpoel.size()%4 == 0, "Inconsistent refined mesh connectivity" );

  std::unordered_set< std::size_t > old( m_inpoel.cbegin(), m_inpoel.cend() );
  std::unordered_set< std::size_t > ref( refinpoel.cbegin(), refinpoel.cend() );

//   std::size_t npoin = 0;
//   for (auto r : rid) if (oid.find(r) != end(oid)) ++npoin;

  auto npoin = ref.size();

  auto& x = m_coord[0];
  auto& y = m_coord[1];
  auto& z = m_coord[2];

  x.resize( npoin );
  y.resize( npoin );
  z.resize( npoin );
  m_gid.resize( npoin );

//   tk::UnsMesh::EdgeNodes prem;  // parents of removed child
//   for (auto o : old)
//     if (ref.find(o) == end(ref))
//       prem[ refiner.node_connectivity.get(o) ] = o;

// std::cout << "old: ";
// for (auto p : old) std::cout << p << ' ';
// std::cout << '\n';
// std::cout << "ref: ";
// for (auto p : ref) std::cout << p << ' ';
// std::cout << '\n';

// std::cout << "prem: ";
// for (const auto& p : prem)
//   std::cout << p.first[0] << ',' << p.first[1] << ':' << p.second << ' ';
// std::cout << '\n';

  for (auto r : ref) {
    if (old.find(r) == end(old)) {
      auto p = refiner.node_connectivity.get( r );
      Assert( old.find(p[0]) != end(old) && old.find(p[1]) != end(old),
              "Parent(s) not in old mesh" );
      Assert( r >= old.size(), "Attempting to overwrite node with added one" );
      x[r] = (x[p[0]] + x[p[1]])/2.0;
      y[r] = (y[p[0]] + y[p[1]])/2.0;
      z[r] = (z[p[0]] + z[p[1]])/2.0;
      decltype(p) gp{{ m_gid[p[0]], m_gid[p[1]] }}; // global parent ids
      //auto g = ( std::hash< std::size_t >()( gp[0] ) ^
      //           std::hash< std::size_t >()( gp[1] ) ) + 10000;
      //auto g = tk::UnsMesh::EdgeHash()( gp ) + 10000;
      //std::cout << CkMyPe() << " hash: " << gp[0] << ',' << gp[1] << ": " << tk::UnsMesh::EdgeHash()( gp ) << '\n';
      auto g = tk::UnsMesh::EdgeHash()( gp );
      Assert( g >= old.size(), "Hashed id overwriting old id" );
      m_gid[r] = g;
      m_lid[g] = r;
      m_coordmap.insert( {g, {{x[r], y[r], z[r]}}} );
      m_edgenode[ gp ] = std::make_tuple( g, x[r], y[r], z[r] );
    }
//     auto i = prem.find(p);
//     if (i != end(prem)) {
// std::cout << "switch: " << r << ", " << i->second << '\n';
//       x[r] = x[ i->second ];
//       y[r] = y[ i->second ];
//       z[r] = z[ i->second ];
//     }
  }

  // Update mesh connectivity with local node IDs
  m_inpoel = refinpoel;

  // Update mesh connectivity with new global node ids
  m_ginpoel = m_inpoel;
  Assert( tk::cunique(m_ginpoel).size() == m_coord[0].size(), "Size mismatch" );
  for (auto& i : m_ginpoel) i = m_gid[i];
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
  auto dist = distribution( m_nchare );

  for (int c=0; c<dist[1]; ++c) {
    // Compute chare ID
    auto cid = CkMyPe() * dist[0] + c;
    // Guard those searches that operate on empty containers in serial
    typename decltype(m_msum)::mapped_type msum;
    if (!m_msum.empty()) msum = tk::cref_find( m_msum, cid );
    // Create worker array element
    m_scheme.discInsert( cid, m_host, m_bc, tk::cref_find(m_chinpoel,cid),
      tk::cref_find(m_chcoordmap,cid), msum, tk::cref_find(m_chfilenodes,cid),
      m_nchare, CkMyPe() );
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
