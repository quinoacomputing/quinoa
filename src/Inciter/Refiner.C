// *****************************************************************************
/*!
  \file      src/Inciter/Refiner.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Mesh refiner for interfacing the mesh refinement library
  \see       Refiner.h for more info.
*/
// *****************************************************************************

#include <vector>
#include <algorithm>

#include "Refiner.h"
#include "Reorder.h"
#include "AMR/mesh_adapter.h"
#include "AMR/Error.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "CGPDE.h"
#include "DGPDE.h"
#include "DerivedData.h"
#include "UnsMesh.h"
#include "Around.h"
#include "ExodusIIMeshWriter.h"
#include "HashMapReducer.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< CGPDE > g_cgpde;
extern std::vector< DGPDE > g_dgpde;

static CkReduction::reducerType BndEdgeMerger;

} // inciter::

using inciter::Refiner;

Refiner::Refiner( const CProxy_Transporter& transporter,
                  const CProxy_Sorter& sorter,
                  const tk::CProxy_Solver& solver,
                  const Scheme& scheme,                
                  const tk::RefinerCallback& cbr,
                  const tk::SorterCallback& cbs,
                  const std::vector< std::size_t >& ginpoel,
                  const tk::UnsMesh::CoordMap& coordmap,
                  const std::map< int, std::vector< std::size_t > >& bface,
                  const std::vector< std::size_t >& triinpoel,
                  int nchare ) :
  m_host( transporter ),
  m_sorter( sorter ),
  m_solver( solver ),
  m_scheme( scheme ),
  m_cbr( cbr ),
  m_cbs( cbs ),
  m_ginpoel( ginpoel ),
  m_coordmap( coordmap ),
  m_bface( bface ),
  m_triinpoel( triinpoel ),
  m_nchare( nchare ),
  m_el( tk::global2local( ginpoel ) ),     // fills m_inpoel, m_gid, m_lid
  m_initref( g_inputdeck.get< tag::amr, tag::init >() ),
  m_refiner( m_inpoel ),
  m_nedge( 0 ),
  m_nref( 0 ),
  m_extra( 1 ),
  m_ch(),
  m_edgenode(),
  m_edgenodeCh(),
  m_bndEdges()
// *****************************************************************************
//  Constructor
//! \param[in] bface Face lists mapped to side set ids
//! \param[in] triinpoel Interconnectivity of points and boundary-face
// *****************************************************************************
{
  Assert( !m_ginpoel.empty(), "No elements assigned to refiner chare" );

  // Reverse initial mesh refinement type list (will pop from back)
  std::reverse( begin(m_initref), end(m_initref) );

  // If initial mesh refinement is configured, start initial mesh refinement.
  // See also tk::grm::check_amr_errors in Control/Inciter/InputDeck/Ggrammar.h.
  if (g_inputdeck.get< tag::amr, tag::initamr >())
    start();
  else
    finish();
}

void
Refiner::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types
//! \details Since this is a [nodeinit] routine, the runtime system executes the
//!   routine exactly once on every logical node early on in the Charm++ init
//!   sequence. Must be static as it is called without an object. See also:
//!   Section "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  BndEdgeMerger = CkReduction::addReducer(
                    tk::mergeHashMap< decltype(m_bndEdges)::key_type,
                                      decltype(m_bndEdges)::mapped_type > );
}

void
Refiner::start()
// *****************************************************************************
// Prepare for next step of mesh refinement
// *****************************************************************************
{
  Assert( (!g_inputdeck.get< tag::amr, tag::init >().empty()) ||
          (!g_inputdeck.get< tag::amr, tag::init >().empty()),
          "Neither initial mesh refinement type list nor user-defined "
          "edge list given" );

  m_extra = 1;  // assume at least a single step of correction is needed
  m_bndEdges.clear();
  m_ch.clear();
  m_edgenodeCh.clear();
  m_edgenode.clear();

  // Generate boundary edges
  bndEdges();
}

void
Refiner::bndEdges()
// *****************************************************************************
// Generate boundary edges and send them to all chares
//! \details Extract edges on the boundary only. The boundary edges (shared by
//!   multiple chares) will be agreed on a refinement that yields a conforming
//!   mesh across chares boundaries.
// *****************************************************************************
{
  // (Re-)compute local from global mesh data (m_inpoel, m_gid, m_lid)
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

  // Aggregate boundary edges across all refiner chares
  decltype(m_bndEdges) bnd{{ thisIndex, std::move(bnded) }};
  auto stream = tk::serialize( bnd );
  contribute( stream.first, stream.second.get(), BndEdgeMerger,
    CkCallback(CkIndex_Refiner::addBndEdges(nullptr),thisProxy) );
}

void
Refiner::addBndEdges( CkReductionMsg* msg )
// *****************************************************************************
//! Receive boundary edges from all refiner chares (including this one)
//! \param[in] msg Charm++ message containing the aggregated map of bnd edges
// *****************************************************************************
{
   PUP::fromMem creator( msg->getData() );
   creator | m_bndEdges;
   delete msg;

   // Compute unique set of chares that share at least a single edge with us
   const auto& ownedges = tk::cref_find( m_bndEdges, thisIndex );
   for (const auto& p : m_bndEdges)    // for all chares
     if (p.first != thisIndex)         // for all chares other than this one
       for (const auto& e : p.second)  // for all boundary edges
         if (ownedges.find(e) != end(ownedges))
           m_ch.insert( p.first );     // if edge is shared, store its chare id

   refine();
}

void
Refiner::refine()
// *****************************************************************************
//  Do a single step of initial mesh refinement based on initial mesh ref list
//! \details This is a single step in a potentially multiple-entry list of
//!   initial adaptive mesh refinement steps. Distribution of the chare-boundary
//!   edges has preceded this step, so that boundary edges (shared by multiple
//!   chares) can agree on a refinement that yields a conforming mesh across
//!   chare boundaries.
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

  for (const auto& e : tk::cref_find(m_bndEdges,thisIndex)) {
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

  for (const auto& e : tk::cref_find(m_bndEdges,thisIndex)) {
    IGNORE(e);
    Assert( m_lid.find( e[0] ) != end( m_lid ) &&
            m_lid.find( e[1] ) != end( m_lid ),
            "Boundary edge not found after refinement" );
  }

  // Ensure valid mesh after refinement
  Assert( tk::positiveJacobians( m_inpoel, m_coord ),
          "Refined mesh cell Jacobian non-positive" );

  // Communicate extra edges
  comExtra();
}


void
Refiner::comExtra()
// *****************************************************************************
// Communicate refined edges after a refinement step
// *****************************************************************************
{
  // Export added nodes on our mesh chunk boundary to other chares (in serial)
  if (m_ch.empty())
    contribute( sizeof(std::size_t), &m_extra, CkReduction::max_ulong,
                m_cbr.get< tag::matched >() );
  else {
    m_nref = 0;
    for (auto c : m_ch) {       // for all chars we share at least an edge with
      // For all boundary edges of chare c, find out if we have added a new
      // node to it, and if so, export parents->(newid,coords) to c.
      tk::UnsMesh::EdgeNodeCoord exp;
      for (const auto& e : tk::cref_find(m_bndEdges,c)) {
        auto i = m_edgenode.find(e);
        if (i != end(m_edgenode)) exp[ e ] = i->second;
      }
      thisProxy[ c ].addRefBndEdges( thisIndex, exp );
    }
  }
}

void
Refiner::addRefBndEdges( int fromch, const tk::UnsMesh::EdgeNodeCoord& ed )
// *****************************************************************************
//! Receive newly added mesh node IDs on our chare boundary
//! \param[in] fromch Chare call coming from
//! \param[in] ed Newly added node IDs associated to parent nodes on chare
//!   boundary
//! \details Receive newly added global node IDs and coordinates associated to
//!   global parent IDs of edges on our mesh chunk boundary.
// *****************************************************************************
{
  // Save/augment buffer of edge-node (IDs and coords) categorized by sender
  // chare
  m_edgenodeCh[ fromch ].insert( begin(ed), end(ed) );
  // Acknowledge receipt of chare-boundary edges to sender
  thisProxy[ fromch ].recvRefBndEdges();
}

void
Refiner::recvRefBndEdges()
// *****************************************************************************
//  Acknowledge received newly added nodes shared with other chares
// *****************************************************************************
{
  // When we have heard from all chares we share at least a single edge with,
  // contribute the number of extra edges that this mesh refinement step has
  // found that were not refined by this chare but were refined by other chares
  // this chare shares the edge with. A global maximum is then be computed on
  // the umber of extra edges appearing in Transporter::matched() which is then
  // used to decide if a new correction step is needed. If this is called for
  // the first time in a given initial mesh refinement step, i.e., not after a
  // correction step, m_extra=1 on all chares, so a correction step is assumed
  // to be required.
  if (++m_nref == m_ch.size()) {
    contribute( sizeof(std::size_t), &m_extra, CkReduction::max_ulong,
                m_cbr.get< tag::matched >() );
  }
}

void
Refiner::correctref()
// *****************************************************************************
//  Correct refinement to arrive at conforming mesh across chare boundaries
//! \details This function is called repeatedly until there is not a a single
//!    edge that needs correction for the whole distributed problem to arrive at
//!    a conforming mesh across chare boundaries during this initial mesh
//!    refinement step.
// *****************************************************************************
{
  // Storage for edges that still need a new node to yield a conforming mesh
  tk::UnsMesh::EdgeSet extra;

  // Ensure that the same global node ID has been assigned by all chares and
  // that the new nodes have the same coordinates generated by potentially
  // multiple chares sharing the refined edge. This is done by searching for all
  // edges that we share with and refined by other chares: (1) If the incoming
  // edge is found among our refined ones, we ensure the newly assigned global
  // IDs equal (independently assigned by multiple chares, this is also a test
  // on the quality of the hash algorithm) and also verify that the new node
  // coordinates equal to machine precision. (2) If the incoming edge is not
  // found among our refined ones, we need to correct the mesh to make it
  // conforming since the edge has been refined by the remote chare. We collect
  // these extra edges, and run a correction refinement, whose result then needs
  // to be communicated again as the new refinement step may introduce new edges
  // that other chares did not yet refine but are shared.
  for (const auto& c : m_edgenodeCh)       // for all chares we share edges with
    for (const auto& e : c.second) {       // for all refined edges on c.first
      auto i = m_edgenode.find( e.first ); // find refined edge given parents
      if (i != end(m_edgenode)) {          // found same added node on edge
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
      } else {  // remote chare added node on edge but we did not
        // Make sure we know about this chare-boundary edge (we did not refine)
        Assert( m_bndEdges.find( thisIndex )->second.find( e.first ) !=
                m_bndEdges.find( thisIndex )->second.end(),
                "Local node IDs of boundary edge not found" );
        // Save edge (given by parent global node IDs) to which the remote chare
        // has added a new node but we did not. Will need to correct the mesh so
        // it conforms across chare boundaries.
        extra.insert( {{ { tk::cref_find( m_lid, e.first[0] ),
                           tk::cref_find( m_lid, e.first[1] ) } }} );
      }
    }

  // Store number of extra edges on this chare which this chare did not add but
  // was refined by another chare, so now we need to tag and refine them and
  // propagate reconnection of neighbor cells to arrive at conforming mesh
  // across chare boundaries.
  m_extra = extra.size();

  // Refine mesh triggered by nodes added on chare-boundary by other chares
  correctRefine( extra );

  // Communicate extra edges. Since refining edges that we did not but other
  // chares did may result in refining new edges that may also be shared along
  // chare boundaries, we now need to communicate these edges and potentially
  // repeat the correction step. This must happen until all chares that share
  // edges can agree that there are no more edges to correct (in which case the
  // above loop finds no extra edges). Only then can this refinement step be
  // considered complete.
  comExtra();
}

void
Refiner::nextref()
// *****************************************************************************
// Decide wether to continue with another step of mesh refinement
//! \details This concludes this mesh refinement step, and we continue if there
//!   are more steps configured by the user.
// *****************************************************************************
{
  // Remove initial mesh refinement step from list
  if (!m_initref.empty()) m_initref.pop_back();

  if (!m_initref.empty())       // Continue to next initial refinement step
    start();
  else {                        // Finish mesh refinement
    // Output final mesh after initial mesh refinement
    tk::UnsMesh refmesh( m_inpoel, m_coord );
    tk::ExodusIIMeshWriter mw( "initref.final." + std::to_string(thisIndex),
                               tk::ExoWriter::CREATE );
    mw.writeMesh( refmesh );
    finish();
  }
}

void
Refiner::finish()
// *****************************************************************************
// Finish initial mesh refinement
//! \details This function is called as after initial mesh refinement has
//!   finished. If initial mesh reifnement was not configured by the user, this
//!   is the point where we continue after the constructor, by computing the
//!   total number of elements across the whole problem.
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::feedback >()) m_host.chrefined();

  // create sorter Charm++ chare array elements using dynamic insertion
  m_sorter[ thisIndex ].insert( m_host, m_solver, m_cbs, m_scheme, m_ginpoel,
    m_coordmap, m_bface, m_triinpoel, m_nchare, CkMyPe() );
 
  // Compute final number of cells across whole problem
  std::vector< std::uint64_t > mesh{ m_ginpoel.size()/4, m_coord[0].size() };
  contribute( mesh, CkReduction::sum_ulong, m_cbr.get< tag::refined >() );
}

void
Refiner::uniformRefine()
// *****************************************************************************
// Do uniform mesh refinement
// *****************************************************************************
{
  // Do uniform refinement
  m_refiner.uniform_refinement();

  // Update mesh coordinates and connectivity
  updateMesh();
}

void
Refiner::errorRefine()
// *****************************************************************************
// Do error-based mesh refinement
// *****************************************************************************
{
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
  for (std::size_t p=0; p<npoin; ++p)   // for all mesh nodes on this chare
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
  m_refiner.error_refinement( edge, crit );

  // Update mesh coordinates and connectivity
  updateMesh();
}

void
Refiner::userRefine()
// *****************************************************************************
// Do mesh refinement based on user explicitly tagging edges
// *****************************************************************************
{
  // Find number of nodes in old mesh
  auto npoin = tk::npoin( m_inpoel );
  // Generate edges surrounding points in old mesh
  auto esup = tk::genEsup( m_inpoel, 4 );
  auto psup = tk::genPsup( m_inpoel, 4, esup );

  // Get user-defined node-pairs (edges) to tag for refinement
  const auto& edgenodelist = g_inputdeck.get< tag::amr, tag::edge >();
  tk::UnsMesh::EdgeSet edgeset;
  for (std::size_t i=0; i<edgenodelist.size()/2; ++i)
    edgeset.insert( {{ {edgenodelist[i*2+0], edgenodelist[i*2+1]} }} );

  // Compute errors in ICs and define refinement criteria for edges
  std::vector< edge_t > edge;
  std::vector< real_t > crit;
  for (std::size_t p=0; p<npoin; ++p)        // for all mesh nodes on this chare
    for (auto q : tk::Around(psup,p)) {      // for all nodes surrounding p
      tk::UnsMesh::Edge e{{p,q}};
      if (edgeset.find(e) != end(edgeset)) { // tag edge if on user's list
        edge.push_back( edge_t(e[0],e[1]) );
        crit.push_back( 1.0 );
      }
    }

  Assert( edge.size() == crit.size(), "Size mismatch" );

  // Do error-based refinement
  m_refiner.error_refinement( edge, crit );

  // Update mesh coordinates and connectivity
  updateMesh();
}

tk::Fields
Refiner::nodeinit( std::size_t npoin,
                       const std::pair< std::vector< std::size_t >,
                          std::vector< std::size_t > >& esup )
// *****************************************************************************
// Evaluate initial conditions (IC) at mesh nodes
//! \param[in] npoin Number points in mesh (on this chare)
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
    for (std::size_t p=0; p<npoin; ++p) {    // for all mesh nodes on this chare
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
Refiner::correctRefine( const tk::UnsMesh::EdgeSet& extra )
// *****************************************************************************
// Do mesh refinement correcting chare-boundary edges
//! \param[in] extra Unique edges that need a new node on chare boundaries
// *****************************************************************************
{
  if (!extra.empty()) {
    // Generate list of edges that need to be corrected
    std::vector< edge_t > edge;
    for (const auto& e : extra) edge.push_back( edge_t(e[0],e[1]) );
    std::vector< real_t > crit( edge.size(), 1.0 );
  
    // Do refinement including edges that need to be corrected
    m_refiner.error_refinement( edge, crit );
  
    // Update mesh coordinates and connectivity
    updateMesh();
  }
}

void
Refiner::updateMesh()
// *****************************************************************************
// Update mesh after refinement
// *****************************************************************************
{
  // Get refined mesh connectivity
  const auto& refinpoel = m_refiner.tet_store.get_active_inpoel();
  Assert( refinpoel.size()%4 == 0, "Inconsistent refined mesh connectivity" );

  // Generate unique node lists of old and refined mesh using local ids
  std::unordered_set< std::size_t > old( m_inpoel.cbegin(), m_inpoel.cend() );
  std::unordered_set< std::size_t > ref( refinpoel.cbegin(), refinpoel.cend() );

  updateVolumeMesh( old, ref );

  updateBoundaryMesh( old, ref );

  // Update mesh connectivity with local node IDs
  m_inpoel = refinpoel;

  // Update mesh connectivity with new global node ids
  m_ginpoel = m_inpoel;
  Assert( tk::uniquecopy(m_ginpoel).size() == m_coord[0].size(),
          "Size mismatch" );
  for (auto& i : m_ginpoel) i = m_gid[i];
}

void
Refiner::updateVolumeMesh( const std::unordered_set< std::size_t >& old,
                               const std::unordered_set< std::size_t >& ref )
// *****************************************************************************
//  Update volume mesh after mesh refinement
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
      auto p = m_refiner.node_connectivity.get( r );
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
Refiner::updateBoundaryMesh( const std::unordered_set< std::size_t >& old,
                                 const std::unordered_set< std::size_t >& ref )
// *****************************************************************************
// Update boundary data structures after mesh refinement
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
  // chare-boundaries due to domain decomposition.
  BndFaces bnd;
  auto oldesuel = tk::genEsuelTet( m_inpoel, tk::genEsup(m_inpoel,4) );
  for (std::size_t e=0; e<oldesuel.size()/4; ++e) {    // for tets on this chare
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

  // Remove chare-boundary faces (to which the above loop did not assign set id)
  auto bndcopy = bnd;
  for (const auto& f : bndcopy) if (f.second.first == -1) bnd.erase( f.first );
  tk::destroy( bndcopy );

  // Now in bnd we have a pair of side set ids of all tet ids that are adjacent
  // to all physical boundary faces that are associated to a side set in the
  // mesh file and only along faces on the physical boundary and not along faces
  // on the chare-boundary.

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
  auto parentFace = [ &old, this ]( const Face& face ){
    std::unordered_set< std::size_t > s;// will store nodes of parent face
    for (auto n : face) {               // for all 3 nodes of the face
      if (old.find(n) != end(old))      // if child node found in old mesh,
        s.insert( n );                  // that node is also in the parent face
      else {                            // if child node is a newly added one
        // find its parent nodes and store both uniquely
        auto p = this->m_refiner.node_connectivity.get( n );
        s.insert( begin(p), end(p) );
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
    auto nc = m_refiner.tet_store.data( f.second.second ).num_children;
    if (nc == 0) {      // if boundary tet is not refined, add its boundary face
      addBndFace( side, f.first );
    } else {            // if boundary tet is refined
      const auto& tets = m_refiner.tet_store.tets;
      for (decltype(nc) i=0; i<nc; ++i ) {      // for all child tets
        // get child tet id
        auto childtet = m_refiner.tet_store.get_child_id( f.second.second, i );
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
            addBndFace( side, {{m_gid[rf[0]], m_gid[rf[1]], m_gid[rf[2]]}} );
          }
        }
      }
    }
  }

  // Update boundary face data structures
  m_bface = std::move(bface);
  m_triinpoel = std::move(triinpoel);
}

#include "NoWarning/refiner.def.h"
