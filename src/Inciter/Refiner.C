// *****************************************************************************
/*!
  \file      src/Inciter/Refiner.C
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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
#include "Centering.h"
#include "Around.h"
#include "Sorter.h"
#include "HashMapReducer.h"
#include "Discretization.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< CGPDE > g_cgpde;
extern std::vector< DGPDE > g_dgpde;

static CkReduction::reducerType BndEdgeMerger;

} // inciter::

using inciter::Refiner;

Refiner::Refiner( const CProxy_Transporter& transporter,
                  const CProxy_Sorter& sorter,
                  const tk::CProxy_MeshWriter& meshwriter,
                  const Scheme& scheme,
                  const tk::RefinerCallback& cbr,
                  const tk::SorterCallback& cbs,
                  const std::vector< std::size_t >& ginpoel,
                  const tk::UnsMesh::CoordMap& coordmap,
                  const std::map< int, std::vector< std::size_t > >& bface,
                  const std::vector< std::size_t >& triinpoel,
                  const std::map< int, std::vector< std::size_t > >& bnode,
                  int nchare ) :
  m_host( transporter ),
  m_sorter( sorter ),
  m_meshwriter( meshwriter ),
  m_scheme( scheme ),
  m_cbr( cbr ),
  m_cbs( cbs ),
  m_ginpoel( ginpoel ),
  m_el( tk::global2local( ginpoel ) ),     // fills m_inpoel, m_gid, m_lid
  m_coordmap( coordmap ),
  m_coord( flatcoord(coordmap) ),
  m_bface( bface ),
  m_bnode( bnode ),
  m_triinpoel( triinpoel ),
  m_nchare( nchare ),
  m_initial( true ),
  m_initref( g_inputdeck.get< tag::amr, tag::init >() ),
  m_ninitref( g_inputdeck.get< tag::amr, tag::init >().size() ),
  m_refiner( m_inpoel ),
  m_nref( 0 ),
  m_extra( 0 ),
  m_ch(),
  m_edgedata(),
  m_edgedataCh(),
  m_bndEdges(),
  m_msumset(),
  m_oldTetIdMap(),
  m_oldTets(),
  m_addedNodes(),
  m_addedTets(),
  m_prevnTets( 0 ),
  m_grand( begin(m_inpoel), end(m_inpoel) )
// *****************************************************************************
//  Constructor
//! \param[in] transporter Transporter (host) proxy
//! \param[in] sorter Mesh reordering (sorter) proxy
//! \param[in] meshwriter Mesh writer proxy
//! \param[in] scheme Discretization scheme
//! \param[in] cbr Charm++ callbacks for Refiner
//! \param[in] cbs Charm++ callbacks for Sorter
//! \param[in] ginpoel Mesh connectivity (this chare) using global node IDs
//! \param[in] coordmap Mesh node coordinates (this chare) for global node IDs
//! \param[in] bface File-internal elem ids of side sets
//! \param[in] bnode Node lists of side sets
//! \param[in] triinpoel Triangle face connectivity with global IDs
//! \param[in] nchare Total number of refiner chares (chare array elements)
// *****************************************************************************
{
  Assert( !m_ginpoel.empty(), "No elements assigned to refiner chare" );
  Assert( tk::positiveJacobians( m_inpoel, m_coord ),
          "Input mesh to Refiner Jacobian non-positive" );
  Assert( !tk::leakyPartition(
            tk::genEsuelTet( m_inpoel, tk::genEsup(m_inpoel,4) ),
            m_inpoel, m_coord ),
          "Input mesh to Refiner leaky" );
  Assert( tk::conforming( m_inpoel, m_coord ),
          "Input mesh to Refiner not conforming" );

  // Reverse initial mesh refinement type list (will pop from back)
  std::reverse( begin(m_initref), end(m_initref) );

  // If initial mesh refinement is configured, start initial mesh refinement.
  // See also tk::grm::check_amr_errors in Control/Inciter/InputDeck/Ggrammar.h.
  if (g_inputdeck.get< tag::amr, tag::t0ref >() && m_ninitref > 0)
    t0ref();
  else
    endt0ref();
}

void
Refiner::sendProxy()
// *****************************************************************************
// Send Refiner proxy to Discretization objects
//! \details This should be called when bound Discretization chare array
//!   elements have already been created.
// *****************************************************************************
{
  // Make sure (bound) Discretization chare is already created and accessible
  Assert( m_scheme.get()[thisIndex].ckLocal() != nullptr,
          "About to dereference nullptr" );

  // Pass Refiner Charm++ chare proxy to fellow (bound) Discretization object
  m_scheme.get()[thisIndex].ckLocal()->setRefiner( thisProxy );
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
Refiner::reorder()
// *****************************************************************************
// Query Sorter and update local mesh with the reordered one
// *****************************************************************************
{
  m_sorter[thisIndex].ckLocal()->mesh( m_ginpoel, m_coordmap, m_triinpoel );

  // Update local mesh data based on data just received from Sorter
  m_el = tk::global2local( m_ginpoel );     // fills m_inpoel, m_gid, m_lid
  m_coord = flatcoord( m_coordmap );

  // WARNING: This re-creates the AMR lib which is probably not what we
  // ultimately want, beacuse this deletes its history recorded during initial
  // (t<0) refinement. However, this appears to correctly update the local mesh
  // based on the reordered one (from Sorter) at least when t0ref is off.
  m_refiner = AMR::mesh_adapter_t( m_inpoel );
}

tk::UnsMesh::Coords
Refiner::flatcoord( const tk::UnsMesh::CoordMap& coordmap )
// *****************************************************************************
// Generate flat coordinate data from coordinate map
//! \param[in] coordmap Coordinates associated to global node IDs of mesh chunk
//! \return Flat coordinate data
// *****************************************************************************
{
  tk::UnsMesh::Coords coord;

  // Convert node coordinates associated to global node IDs to a flat vector
  auto npoin = coordmap.size();
  Assert( m_gid.size() == npoin, "Size mismatch" );
  coord[0].resize( npoin );
  coord[1].resize( npoin );
  coord[2].resize( npoin );
  for (const auto& c : coordmap) {
    auto i = tk::cref_find( m_lid, c.first );
    Assert( i < npoin, "Indexing out of coordinate map" );
    coord[0][i] = c.second[0];
    coord[1][i] = c.second[1];
    coord[2][i] = c.second[2];
  }

  return coord;
}

void
Refiner::dtref( const std::map< int, std::vector< std::size_t > >& bface,
                const std::map< int, std::vector< std::size_t > >& bnode,
                const std::vector< std::size_t >& triinpoel )
// *****************************************************************************
// Start mesh refinement (during time stepping, t>0)
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] bnode Boundary-node lists mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
// *****************************************************************************
{
  m_initial = false;

  // Update boundary node lists
  m_bface = bface;
  m_bnode = bnode;
  m_triinpoel = triinpoel;

  start();
}

void
Refiner::t0ref()
// *****************************************************************************
// Output mesh to file before a new step mesh refinement
// *****************************************************************************
{
  Assert( m_ninitref > 0, "No initial mesh refinement steps configured" );
  // Output initial mesh to file
  auto l = m_ninitref - m_initref.size();  // num initref steps completed
  auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();
  if (l == 0)
    writeMesh( "t0ref", l, t0-1.0,
      CkCallback( CkIndex_Refiner::start(), thisProxy[thisIndex] ) );
  else
    start();
}

void
Refiner::start()
// *****************************************************************************
//  Start new step of initial mesh refinement
// *****************************************************************************
{
  m_extra = 0;
  m_bndEdges.clear();
  m_ch.clear();
  m_edgedataCh.clear();

  updateEdgeData();

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
  // Generate boundary edges of our mesh chunk
  tk::UnsMesh::EdgeSet bnded;
  auto esup = tk::genEsup( m_inpoel, 4 );         // elements surrounding points
  auto esuel = tk::genEsuelTet( m_inpoel, esup ); // elems surrounding elements
  for (std::size_t e=0; e<esuel.size()/4; ++e) {
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f) {
      if (esuel[mark+f] == -1) {
        auto A = m_ginpoel[ mark+tk::lpofa[f][0] ];
        auto B = m_ginpoel[ mark+tk::lpofa[f][1] ];
        auto C = m_ginpoel[ mark+tk::lpofa[f][2] ];
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

  contribute( m_cbr.get< tag::edges >() );
}

void
Refiner::refine()
// *****************************************************************************
//  Do a single step of mesh refinement
//! \details During initial (t<0) mesh refinement, this is a single step in a
//!   potentially multiple-entry list of initial adaptive mesh refinement steps.
//!   Distribution of the chare-boundary edges must have preceded this step, so
//!   that boundary edges (shared by multiple chares) can agree on a refinement
//!   that yields a conforming mesh across chare boundaries.
//!   During-timestepping (dtref) mesh refinement this is called once, as we
//!   only do a single step during time stepping.
// *****************************************************************************
{
  // Perform leak test on old mesh
  Assert( !tk::leakyPartition(
            tk::genEsuelTet( m_inpoel, tk::genEsup(m_inpoel,4) ),
            m_inpoel, m_coord ),
          "Mesh partition before refinement leaky" );

  for (const auto& e : tk::cref_find(m_bndEdges,thisIndex)) {
    IGNORE(e);
    Assert( m_lid.find( e[0] ) != end( m_lid ) &&
            m_lid.find( e[1] ) != end( m_lid ),
            "Boundary edge not found before refinement" );
  }

  if (m_initial) {      // if initial (before t=0) AMR

    // Refine mesh based on next initial refinement type
    if (!m_initref.empty()) {
      auto r = m_initref.back();    // consume (reversed) list from its back
      if (r == ctr::AMRInitialType::UNIFORM)
        uniformRefine();
      else if (r == ctr::AMRInitialType::INITIAL_CONDITIONS)
        errorRefine();
      else if (r == ctr::AMRInitialType::COORDINATES)
        coordRefine();
      else if (r == ctr::AMRInitialType::EDGELIST)
        edgelistRefine();
      else Throw( "Initial AMR type not implemented" );
    }

  } else {              // if AMR during time stepping (t>0)

    if (g_inputdeck.get< tag::amr, tag::dtref_uniform >())
      uniformRefine();
    else
      errorRefine();

  }

//   for (const auto& e : tk::cref_find(m_bndEdges,thisIndex)) {
//     IGNORE(e);
//     Assert( m_lid.find( e[0] ) != end( m_lid ) &&
//             m_lid.find( e[1] ) != end( m_lid ),
//             "Boundary edge not found after refinement" );
//   }

  // Communicate extra edges
  comExtra();
}

void
Refiner::comExtra()
// *****************************************************************************
// Communicate extra edges along chare boundaries
// *****************************************************************************
{
  // Export extra added nodes on our mesh chunk boundary to other chares
  if (m_ch.empty())
    correctref();
  else {
    for (auto c : m_ch) {       // for all chares we share at least an edge with
      // send refinement data of all our edges
      thisProxy[ c ].addRefBndEdges( thisIndex, m_edgedata, m_intermediates );
    }
  }
}

void
Refiner::addRefBndEdges(
  int fromch,
  const AMR::EdgeData& ed,
  const std::unordered_set< std::size_t >& intermediates )
// *****************************************************************************
//! Receive tagged edges on our chare boundary
//! \param[in] fromch Chare call coming from
//! \param[in] ed Tagged edges on chare boundary
//! \param[in] intermediates Intermediate nodes
// *****************************************************************************
{
  // Save/augment buffer of edge data for each sender chare
  m_edgedataCh[ fromch ].insert( begin(ed), end(ed) );

  // Add the intermediates to mesh refiner lib
  for (const auto i : intermediates) {
    auto l = m_lid.find( i ); // convert to local node ids
    //auto l1 = tk::cref_find( m_lid, r.first[0] );
    //auto it = m_edgedata.find( r.first );
    if (l != end(m_lid)) {
      m_refiner.tet_store.intermediate_list.insert( l->second );
    }
  }

  // Heard from every worker we share at least a single edge with
  if (++m_nref == m_ch.size()) {
    m_nref = 0;
    correctref();
  }
}

void
Refiner::correctref()
// *****************************************************************************
//  Correct extra edges to arrive at conforming mesh across chare boundaries
//! \details This function is called repeatedly until there is not a a single
//!    edge that needs correction for the whole distributed problem to arrive at
//!    a conforming mesh across chare boundaries during a mesh refinement step.
// *****************************************************************************
{
  m_refiner.lock_intermediates();

  auto unlocked = AMR::Edge_Lock_Case::unlocked;

  // Storage for edge data that need correction to yield a conforming mesh
  AMR::EdgeData extra;

  // loop through all edges shared with other chares
  for (const auto& c : m_edgedataCh)       // for all chares we share edges with
    for (const auto& r : c.second) {       // for all edges shared with c.first
      // find local data of remote edge { r.first[0] - r.first[1] }
      auto it = m_edgedata.find( r.first );

      if (it != end(m_edgedata))
      {
        const auto& local = it->second;
        const auto& remote = r.second;
        auto local_needs_refining = local.first;
        auto local_lock_case = local.second;
        auto remote_needs_refining = remote.first;
        auto remote_lock_case = remote.second;

        auto local_needs_refining_orig = local_needs_refining;
        auto local_lock_case_orig = local_lock_case;

        Assert( !(local_lock_case > unlocked && local_needs_refining),
                "Invalid local edge: locked & needs refining" );
        Assert( !(remote_lock_case > unlocked && remote_needs_refining),
                "Invalid remote edge: locked & needs refining" );

        // compute lock from local and remote locks as most restrictive
        local_lock_case = std::max( local_lock_case, remote_lock_case );

        if (local_lock_case > unlocked) {
            local_needs_refining = false;
        }

        if (local_lock_case == unlocked && remote_needs_refining)
        {
          local_needs_refining = true;
        }

        if (local_lock_case != local_lock_case_orig ||
            local_needs_refining != local_needs_refining_orig) {

           auto l1 = tk::cref_find( m_lid, r.first[0] );
           auto l2 = tk::cref_find( m_lid, r.first[1] );

           Assert( l1 != l2, "Edge end-points local ids are the same" );

           /*
           std::cout << thisIndex << " correcting "
             << r.first[0] << '-' << r.first[1] << ": l{"
             << local_needs_refining_orig << ',' << local_lock_case_orig
             << "} + r{" << remote_needs_refining_orig << ','
             << remote_lock_case_orig << "} -> {" << local_needs_refining
             << ',' << local_lock_case << "}\n";
           */

           extra[ {{ std::min(l1,l2), std::max(l1,l2) }} ] =
             { local_needs_refining, local_lock_case };
         }

      }
    }

  m_extra = extra.size();

  correctRefine( extra );

  // Aggregate number of extra edges that still need correction
  std::vector< std::size_t > m{ m_extra, m_edgedata.size() };
  contribute( m, CkReduction::sum_ulong, m_cbr.get< tag::matched >() );
}

void
Refiner::correctRefine( const AMR::EdgeData& extra )
// *****************************************************************************
// Do mesh refinement correcting chare-boundary edges
//! \param[in] extra Unique edges that need a new node on chare boundaries
// *****************************************************************************
{
  if (!extra.empty()) {
    // Do refinement including edges that need to be corrected
    m_refiner.mark_error_refinement_corr( extra );
    // Update our extra-edge store based on refiner
    updateEdgeData();
  }
}

void
Refiner::updateEdgeData()
// *****************************************************************************
// Query AMR lib and update our local store of edge data
// *****************************************************************************
{

  using Edge = tk::UnsMesh::Edge;
  const auto& ref_edges = m_refiner.tet_store.edge_store.edges;

  m_edgedata.clear();
  m_intermediates.clear();

  // This currently takes ALL edges from the AMR lib, i.e., on the whole
  // domain. We should eventually only collect edges here that are shared with
  // other chares.
  for (const auto& e : ref_edges) {
    const auto& ed = e.first.get_data();
    const auto ged = Edge{{ m_gid[ ed[0] ], m_gid[ ed[1] ] }};
    //std::cout << "Building to send " << ed[0] << "-" << ed[1] << " aka " << m_gid[ ed[0] ] << "-" << m_gid[ ed[1] ] << " val " << e.second.needs_refining << " and " << e.second.lock_case <<  std::endl;
    m_edgedata[ ged ] = { e.second.needs_refining, e.second.lock_case };
  }

  // Build intermediates to send. This currently takes ALL intermediates from
  // the AMR lib, i.e., on the whole domain. We should eventually only collect
  // edges here that are shared with other chares.
  for (const auto& i : m_refiner.tet_store.intermediate_list) {
     m_intermediates.insert( m_gid[i] );
  }
}

void
Refiner::writeMesh( const std::string& basefilename,
                    uint64_t itr,
                    tk::real t,
                    CkCallback c ) const
// *****************************************************************************
//  Output mesh to file(s)
//! \param[in] basefilename File name to append to
//! \param[in] itr Iteration count since a new mesh
//! \param[in] t "Physical time" to write to file. "Time" here is used to
//!   designate a new time step at which the mesh is saved.
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  // Prepare element fields with mesh refinement data
  std::vector< std::string > elemfieldnames{ "refinement level", "cell type" };
  auto& tet_store = m_refiner.tet_store;
  std::vector< std::vector< tk::real > > elemfields{
    tet_store.get_refinement_level_list(), tet_store.get_cell_type_list() };

  // Prepare solution field names: depvar + component id for all eqs
  auto nprop = g_inputdeck.get< tag::component >().nprop();
  auto solfieldnames = g_inputdeck.get<tag::component>().depvar( g_inputdeck );
  Assert( solfieldnames.size() == nprop, "Size mismatch" );

  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
  const auto centering = ctr::Scheme().centering( scheme );
  auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();

  std::vector< std::vector< tk::real > > nodefields;
  std::vector< std::string > nodefieldnames;

  // Prepare node or element fields for output to file
  if (centering == tk::Centering::NODE) {

    // Augment element field names with solution variable names + field ids
    nodefieldnames.insert( end(nodefieldnames),
                           begin(solfieldnames), end(solfieldnames) );

    // Evaluate initial conditions on current mesh at t0
    tk::Fields u( m_coord[0].size(), nprop );
    for (const auto& eq : g_cgpde) eq.initialize( m_coord, u, t0 );

    // Extract all scalar components from solution for output to file
    for (std::size_t i=0; i<nprop; ++i)
      nodefields.push_back( u.extract( i, 0 ) );

  } else if (centering == tk::Centering::ELEM) {

    // Augment element field names with solution variable names + field ids
    elemfieldnames.insert( end(elemfieldnames),
                           begin(solfieldnames), end(solfieldnames) );

    auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
    tk::Fields lhs( m_inpoel.size()/4, ndof*nprop );

    // Generate left hand side for DG initialize
    auto geoElem = tk::genGeoElemTet( m_inpoel, m_coord );
    for (const auto& eq : g_dgpde) eq.lhs( geoElem, lhs );

    // Evaluate initial conditions on current mesh at t0
    auto u = lhs;
    for (const auto& eq : g_dgpde)
      eq.initialize( lhs, m_inpoel, m_coord, u, t0, m_inpoel.size()/4 );

    // Extract all scalar components from solution for output to file
    for (std::size_t i=0; i<nprop; ++i)
      elemfields.push_back( u.extract( i, 0 ) );
  }

  // Output mesh
  m_meshwriter[ CkNodeFirst( CkMyNode() ) ].
    write( /*meshoutput = */ true, /*fieldoutput = */ true, itr, 1, t,
           thisIndex, basefilename, m_inpoel, m_coord, m_bface, m_bnode,
           tk::remap(m_triinpoel,m_lid), elemfieldnames, nodefieldnames,
           elemfields, nodefields, c );
}

void
Refiner::eval()
// *****************************************************************************
// Refine mesh and decide how to continue
//! \details First the mesh refiner object is called to perform a single step
//!   of mesh refinement. Then, if this function is called during a step
//!   (potentially multiple levels of) initial AMR, it evaluates whether to do
//!   another one. If it is called during time stepping, this concludes the
//!   single mesh refinement step and the new mesh is sent to the PDE worker
//!   (Discretization).
// *****************************************************************************
{
  // Save old tet id map before performing refinement
  m_oldTetIdMap = m_refiner.tet_store.get_active_id_mapping();

  // Save old tet ids before performing refinement
  m_prevnTets = m_oldTets.size();       // save number tets before refinement
  m_oldTets.clear();
  for (const auto& t : m_refiner.tet_store.tets) m_oldTets.insert( t.first );

  m_refiner.perform_refinement();

  updateMesh();

  if (m_initial) {      // if initial (before t=0) AMR
    auto l = m_ninitref - m_initref.size() + 1;  // num initref steps completed
    auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();
    // Generate times equally subdividing t0-1...t0 to ninitref steps
    auto t =
      t0 - 1.0 + static_cast<tk::real>(l)/static_cast<tk::real>(m_ninitref);
    auto itr = l;
    // Output mesh after refinement step
    writeMesh( "t0ref", itr, t,
               CkCallback( CkIndex_Refiner::next(), thisProxy[thisIndex] ) );
  } else next();
}

void
Refiner::next()
// *****************************************************************************
// Continue after finishing a refinement step
// *****************************************************************************
{
  if (m_initial) {      // if initial (before t=0) AMR

    // Remove initial mesh refinement step from list
    if (!m_initref.empty()) m_initref.pop_back();
    // Continue to next initial AMR step or finish
    if (!m_initref.empty()) t0ref(); else endt0ref();

  } else {              // if AMR during time stepping (t>0)

    // Augment node communication map with newly added nodes on chare-boundary
    for (const auto& c : m_edgedataCh) {
      auto& nodes = tk::ref_find( m_msumset, c.first );
      for (const auto& e : c.second) {
        // If parent nodes were part of the node communication map for chare
        if (nodes.find(e.first[0]) != end(nodes) &&
            nodes.find(e.first[1]) != end(nodes))
        {
          // Add new node if local id was generated for it
          auto n = tk::UnsMesh::Hash<2>()( e.first );
          if (m_lid.find(n) != end(m_lid)) nodes.insert( n );
        }
      }
    }

    // Send new mesh and solution back to PDE worker
    Assert( m_scheme.get()[thisIndex].ckLocal() != nullptr,
            "About to use nullptr" );
    auto e = tk::element< SchemeBase::ProxyElem >
                        ( m_scheme.getProxy(), thisIndex );
    std::unordered_map< int, std::vector< std::size_t > > msum;
    for (const auto& c : m_msumset) {
      auto& n = msum[ c.first ];
      n.insert( end(n), c.second.cbegin(), c.second.cend() );
    }

    boost::apply_visitor(
      ResizeAfterRefined( m_ginpoel, m_el, m_coord, m_addedNodes, m_addedTets,
        msum, m_bface, m_bnode, m_triinpoel ), e );

  }
}

void
Refiner::endt0ref()
// *****************************************************************************
// Finish initial mesh refinement
//! \details This function is called as after initial mesh refinement has
//!   finished. If initial mesh reifnement was not configured by the user, this
//!   is the point where we continue after the constructor, by computing the
//!   total number of elements across the whole problem.
// *****************************************************************************
{
  // create sorter Charm++ chare array elements using dynamic insertion
  m_sorter[ thisIndex ].insert( m_host, m_meshwriter, m_cbs, m_scheme,
    CkCallback( CkIndex_Refiner::reorder(), thisProxy[thisIndex] ),
    m_ginpoel, m_coordmap, m_bface, m_triinpoel, m_bnode, m_nchare );

  // Compute final number of cells across whole problem
  std::vector< std::size_t > meshsize{{ m_ginpoel.size()/4,
                                        m_coord[0].size() }};
  contribute( meshsize, CkReduction::sum_ulong, m_cbr.get< tag::refined >() );
}

void
Refiner::uniformRefine()
// *****************************************************************************
// Do uniform mesh refinement
// *****************************************************************************
{
  // Do uniform refinement
  m_refiner.mark_uniform_refinement();

  // Update our extra-edge store based on refiner
  updateEdgeData();

  // Set number of extra edges to be zero, skipping correction (if this is the
  // only step in this refinement step)
  m_extra = 0;
}

void
Refiner::errorRefine()
// *****************************************************************************
// Do error-based mesh refinement
// *****************************************************************************
{
  // Find number of nodes in old mesh
  auto npoin = tk::npoin_in_graph( m_inpoel );
  // Generate edges surrounding points in old mesh
  auto esup = tk::genEsup( m_inpoel, 4 );
  auto psup = tk::genPsup( m_inpoel, 4, esup );

  // Get solution whose error to evaluate
  tk::Fields u;
  if (m_initial) {      // if initial (before t=0) AMR

    // Evaluate initial conditions at mesh nodes
    u = nodeinit( npoin, esup );

  } else {              // if AMR during time stepping (t>0)

    // ...

  }

  // Get the indices (in the system of systems) of refinement variables and the
  // error indicator configured
  const auto& refidx = g_inputdeck.get< tag::amr, tag::id >();
  auto errtype = g_inputdeck.get< tag::amr, tag::error >();

  using AMR::edge_t;

  // Compute errors in ICs and define refinement criteria for edges
  std::vector< edge_t > edge;
  AMR::Error error;
  for (std::size_t p=0; p<npoin; ++p)   // for all mesh nodes on this chare
  {
    for (auto q : tk::Around(psup,p)) { // for all nodes surrounding p
       tk::real cmax = 0.0;
       edge_t e(p,q);
       for (auto i : refidx) {          // for all refinement variables
         auto c = error.scalar( u, e, i, m_coord, m_inpoel, esup, errtype );
         if (c > cmax) cmax = c;        // find max error at edge
       }
       if (cmax > 0.8) {         // if nonzero error, will pass edge to refiner
         edge.push_back( e );
       }
     }
  }

  // Do error-based refinement
  m_refiner.mark_error_refinement( edge );

  // Update our extra-edge store based on refiner
  updateEdgeData();

  // Set number of extra edges to a nonzero number, triggering correction
  m_extra = 1;
}

void
Refiner::edgelistRefine()
// *****************************************************************************
// Do mesh refinement based on user explicitly tagging edges
// *****************************************************************************
{
  // Get user-defined node-pairs (edges) to tag for refinement
  const auto& edgenodelist = g_inputdeck.get< tag::amr, tag::edge >();

  if (!edgenodelist.empty()) {  // if user explicitly tagged edges
    // Find number of nodes in old mesh
    auto npoin = tk::npoin_in_graph( m_inpoel );
    // Generate edges surrounding points in old mesh
    auto esup = tk::genEsup( m_inpoel, 4 );
    auto psup = tk::genPsup( m_inpoel, 4, esup );

    tk::UnsMesh::EdgeSet useredges;
    for (std::size_t i=0; i<edgenodelist.size()/2; ++i)
      useredges.insert( {{ {edgenodelist[i*2+0], edgenodelist[i*2+1]} }} );

    using AMR::edge_t;

    // Tag edges the user configured
    std::vector< edge_t > edge;
    //std::cout << thisIndex << ": ";
    for (std::size_t p=0; p<npoin; ++p)        // for all mesh nodes on this chare
      for (auto q : tk::Around(psup,p)) {      // for all nodes surrounding p
        tk::UnsMesh::Edge e{{ m_gid[p], m_gid[q] }};
        //std::cout << e[0] << ',' << e[1] << ' ';
        auto f = useredges.find(e);
        if (f != end(useredges)) { // tag edge if on user's list
          edge.push_back( edge_t(p,q) );
          useredges.erase( f );
        }
      }
    //std::cout << std::endl;

    //std::cout << thisIndex << ": " << edge.size() << std::endl;

    if (!useredges.empty()) {
      std::cout << "Edges tagged but not found on chare " << thisIndex << ": ";
      for (const auto& e : useredges) std::cout << e[0] << ',' << e[1] << ' ';
    }
    std::cout << std::endl;

    std::cout << thisIndex << " tet store: ";
    for (const auto& t : m_refiner.tet_store.tets)
      std::cout << t.second[0] << ',' << t.second[1] << ','
                << t.second[2] << ',' << t.second[3] << ' ';
    std::cout << std::endl;

    // Do error-based refinement
    m_refiner.mark_error_refinement( edge );

    // Update our extra-edge store based on refiner
    updateEdgeData();

    // Set number of extra edges to a nonzero number, triggering correction
    m_extra = 1;
  }
}

void
Refiner::coordRefine()
// *****************************************************************************
// Do mesh refinement based on tagging edges based on end-point coordinates
// *****************************************************************************
{
  // Get user-defined half-world coordinates
  auto xminus = g_inputdeck.get< tag::amr, tag::xminus >();
  auto xplus = g_inputdeck.get< tag::amr, tag::xplus >();
  auto yminus = g_inputdeck.get< tag::amr, tag::yminus >();
  auto yplus = g_inputdeck.get< tag::amr, tag::yplus >();
  auto zminus = g_inputdeck.get< tag::amr, tag::zminus >();
  auto zplus = g_inputdeck.get< tag::amr, tag::zplus >();

  // The default is the largest representable double
  auto rmax = std::numeric_limits< kw::amr_xminus::info::expect::type >::max();
  auto eps =
    std::numeric_limits< kw::amr_xminus::info::expect::type >::epsilon();

  // Decide if user has configured the half-world
  bool xm = std::abs(xminus - rmax) > eps ? true : false;
  bool xp = std::abs(xplus - rmax) > eps ? true : false;
  bool ym = std::abs(yminus - rmax) > eps ? true : false;
  bool yp = std::abs(yplus - rmax) > eps ? true : false;
  bool zm = std::abs(zminus - rmax) > eps ? true : false;
  bool zp = std::abs(zplus - rmax) > eps ? true : false;

  using AMR::edge_t;

  if (xm || xp || ym || yp || zm || zp) {       // if any half-world configured
    // Find number of nodes in old mesh
    auto npoin = tk::npoin_in_graph( m_inpoel );
    // Generate edges surrounding points in old mesh
    auto esup = tk::genEsup( m_inpoel, 4 );
    auto psup = tk::genPsup( m_inpoel, 4, esup );
    // Get access to node coordinates
    const auto& x = m_coord[0];
    const auto& y = m_coord[1];
    const auto& z = m_coord[2];
    // Compute edges to be tagged for refinement
    std::vector< edge_t > edge;
    for (std::size_t p=0; p<npoin; ++p)        // for all mesh nodes on this chare
      for (auto q : tk::Around(psup,p)) {      // for all nodes surrounding p
        tk::UnsMesh::Edge e{{p,q}};

        bool t = true;
        if (xm) { if (x[p]>xminus && x[q]>xminus) t = false; }
        if (xp) { if (x[p]<xplus && x[q]<xplus) t = false; }
        if (ym) { if (y[p]>yminus && y[q]>yminus) t = false; }
        if (yp) { if (y[p]<yplus && y[q]<yplus) t = false; }
        if (zm) { if (z[p]>zminus && z[q]>zminus) t = false; }
        if (zp) { if (z[p]<zplus && z[q]<zplus) t = false; }

        if (t) edge.push_back( edge_t(e[0],e[1]) );
      }

    // Do error-based refinement
    m_refiner.mark_error_refinement( edge );

    // Update our extra-edge store based on refiner
    updateEdgeData();

    // Set number of extra edges to a nonzero number, triggering correction
    m_extra = 1;
  }
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
  const auto centering = ctr::Scheme().centering( scheme );
  if (centering == tk::Centering::NODE) {

    // Node-centered: evaluate ICs for all scalar components integrated
    for (const auto& eq : g_cgpde) eq.initialize( m_coord, u, t0 );

  } else if (centering == tk::Centering::ELEM) {

    auto esuel = tk::genEsuelTet( m_inpoel, esup ); // elems surrounding elements
    // Initialize cell-based unknowns
    tk::Fields ue( m_inpoel.size()/4, nprop );
    auto lhs = ue;
    auto geoElem = tk::genGeoElemTet( m_inpoel, m_coord );
    for (const auto& eq : g_dgpde)
      eq.lhs( geoElem, lhs );
    for (const auto& eq : g_dgpde)
      eq.initialize( lhs, m_inpoel, m_coord, ue, t0, esuel.size()/4 );

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

  } else Throw( "Scheme centring not handled for nodal initialization" );

  Assert( u.nunk() == m_coord[0].size(), "Size mismatch" );
  Assert( u.nprop() == nprop, "Size mismatch" );

  return u;
}

void
Refiner::updateMesh()
// *****************************************************************************
// Update old mesh after refinement
// *****************************************************************************
{
  // Get refined mesh connectivity
  const auto& refinpoel = m_refiner.tet_store.get_active_inpoel();
  Assert( refinpoel.size()%4 == 0, "Inconsistent refined mesh connectivity" );
  Assert( tk::conforming( m_inpoel, m_coord ),
          "Mesh not conforming after refinement" );

  // Generate unique node lists of old and refined mesh using local ids
  std::unordered_set< std::size_t > old( begin(m_inpoel), end(m_inpoel) );
  std::unordered_set< std::size_t > ref( begin(refinpoel), end(refinpoel) );

  // Get nodal communication map from Discretization worker
  if (!m_initial) m_msumset = m_scheme.get()[thisIndex].ckLocal()->msumset();

  // Update mesh and solution after refinement
  newVolMesh( old, ref );
  newBndMesh( old, ref );

  // Save/update unique nodes of the current 'old' mesh that will become a
  // 'grandparent' mesh in the next refinement step with local ids. This mesh is
  // used in the next refinement step (if any) to find parent nodes of faces of
  // backtracker tets that will have undergone 2->8 or 4->8 refinement.
  m_grand = old;

  // Update mesh connectivity with local node IDs
  m_inpoel = m_refiner.tet_store.get_active_inpoel();

  // Update mesh connectivity with new global node ids
  m_ginpoel = m_inpoel;
  Assert( tk::uniquecopy(m_ginpoel).size() == m_coord[0].size(),
          "Size mismatch" );
  // cppcheck-suppress useStlAlgorithm
  for (auto& i : m_ginpoel) i = m_gid[i];

  // Update flat coordinates storage
  //m_coord = flatcoord( m_coordmap );

  // Ensure valid mesh after refinement
  Assert( tk::positiveJacobians( m_inpoel, m_coord ),
          "Refined mesh cell Jacobian non-positive" );

  // Perform leak test on new mesh
  Assert( !tk::leakyPartition(
            tk::genEsuelTet( m_inpoel, tk::genEsup(m_inpoel,4) ),
            m_inpoel, m_coord ),
          "Refined mesh partition leaky" );

  Assert( tk::conforming( m_inpoel, m_coord ),
          "Mesh not conforming after updating mesh after mesh refinement" );
}

void
Refiner::newVolMesh( const std::unordered_set< std::size_t >& old,
                     const std::unordered_set< std::size_t >& ref )
// *****************************************************************************
//  Compute new volume mesh after mesh refinement
//! \param[in] old Unique nodes of the old (unrefined) mesh using local ids
//! \param[in] ref Unique nodes of the refined mesh using local ids
// *****************************************************************************
{
  auto& x = m_coord[0];
  auto& y = m_coord[1];
  auto& z = m_coord[2];

  // Resize node coordinates, global ids, and added-nodes map
  auto npoin = ref.size();
  x.resize( npoin );
  y.resize( npoin );
  z.resize( npoin );
  m_gid.resize( npoin, std::numeric_limits< std::size_t >::max() );
  m_addedNodes.clear();

  // Generate coordinates and ids to newly added nodes after refinement step
  for (auto r : ref) {               // for all unique nodes of the refined mesh
    if (old.find(r) == end(old)) {   // if node is newly added (in this step)
      // get (local) parent ids of newly added node
      auto p = m_refiner.node_connectivity.get( r );
      Assert( old.find(p[0]) != end(old) && old.find(p[1]) != end(old),
              "Parent(s) not in old mesh" );
      Assert( r >= old.size(), "Attempting to overwrite node with added one" );
      // global parent ids
      decltype(p) gp{{ m_gid[p[0]], m_gid[p[1]] }};
      // generate new global ID for newly added node
      auto g = tk::UnsMesh::Hash<2>()( gp );

      // if node added by AMR lib has not yet been added to Refiner's new mesh
      if (m_coordmap.find(g) == end(m_coordmap)) {
        // ensure newly generated node id has not yet been used
        Assert( g >= old.size(), "Hashed id overwriting old id" );
        // generate coordinates for newly added node
        x[r] = (x[p[0]] + x[p[1]])/2.0;
        y[r] = (y[p[0]] + y[p[1]])/2.0;
        z[r] = (z[p[0]] + z[p[1]])/2.0;
        // store newly added node id and their parent ids (local ids)
        m_addedNodes[r] = p;
        // assign new global ids to local->global and to global->local maps
        m_gid[r] = g;
        Assert( m_lid.find(g) == end(m_lid),
                "Overwriting entry global->local node ID map" );
        m_lid[g] = r;
        // assign new coordinates to new global node id
        Assert( m_coordmap.find(g) == end(m_coordmap),
                "Overwriting entry already in coordmap" );
        m_coordmap.insert( {g, {{x[r], y[r], z[r]}}} );
      }
    }
  }

  Assert( m_gid.size() == m_lid.size(), "Size mismatch" );

  Assert( std::none_of( begin(m_gid), end(m_gid), [](std::size_t i){
            return i == std::numeric_limits< std::size_t >::max(); } ),
          "Not all local->global node IDs have been assigned" );

  // Extract new tet ids from AMR object after refinement step
  const auto& tet_store = m_refiner.tet_store;
  std::vector< std::size_t > newtets;
  for (const auto& t : tet_store.tets)
    if (m_oldTets.find(t.first) == end(m_oldTets))
      newtets.push_back( t.first );

  // Invert AMR's tet id map
  std::unordered_map< std::size_t, std::size_t > newids;
  std::size_t j = 0;
  for (auto t : m_refiner.tet_store.get_active_id_mapping())
    newids[t] = j++;

  // Generate child->parent tet id map after refinement step
  m_addedTets.clear();
  for (auto n : newtets) {
     auto parent = tet_store.data( n ).parent_id;
     Assert( parent < m_oldTets.size(),
             "Parent tet id will index out of old solution vector" );
     auto child = tk::cref_find( newids, n );
     Assert( child < m_oldTets.size() + newtets.size(),
             "New tet id will index out of new solution vector" );
     m_addedTets[ child ] = parent - m_prevnTets;
  }
}

Refiner::BndFaceData
Refiner::boundary()
// *****************************************************************************
//  Generate boundary data structures used to update refined boundary faces and
//  nodes of side sets
//! \return A tuple of boundary face data. See type definition.
//! \details The output of this function is used to regenerate physical boundary
//!   face and node data structures after refinement, see updateBndFaces() and
//!   updateBndNodes().
// *****************************************************************************
{
  using Face = tk::UnsMesh::Face;
  using Tet = tk::UnsMesh::Tet;

  // Generate the inverse of AMR's tet store
  std::unordered_map< Tet,
                      std::size_t,
                      tk::UnsMesh::Hash<4>,
                      tk::UnsMesh::Eq<4> > invtets;
  for (const auto& t : m_refiner.tet_store.tets) invtets[ t.second ] = t.first;

  // Generate data structure that associates the id of a tet adjacent to a
  // boundary triangle face for all (physical and chare) boundary faces.
  BndFaceTets facetets;
  auto oldesuel = tk::genEsuelTet( m_inpoel, tk::genEsup(m_inpoel,4) );
  for (std::size_t e=0; e<oldesuel.size()/4; ++e) {
    auto m = e*4;
    for (std::size_t f=0; f<4; ++f) {
      if (oldesuel[m+f] == -1) {  // if a face does not have an adjacent tet
        Face b{{ m_ginpoel[ m+tk::lpofa[f][0] ],
                 m_ginpoel[ m+tk::lpofa[f][1] ],
                 m_ginpoel[ m+tk::lpofa[f][2] ] }};
        Tet t{{ m_inpoel[m+0], m_inpoel[m+1], m_inpoel[m+2], m_inpoel[m+3] }};
        // associate tet id to adjacent (physical or chare) boundary face
        facetets[b] = tk::cref_find(invtets,t);
      }
    }
  }

  // Generate unique set of faces for each side set
  BndFaceSets facesets;
  for (const auto& s : m_bface) {  // for all phsyical boundaries (sidesets)
    auto& faces = facesets[ s.first ];
    for (auto f : s.second) {
      faces.insert(
        {{{ m_triinpoel[f*3+0], m_triinpoel[f*3+1], m_triinpoel[f*3+2] }}} );
    }
  }

  // Associate a unique set of side set ids to boundary faces of the parents of
  // backtrackers (that underwent a 2->8 or a 4->8 refinement). Note that this
  // will only store data for backtrackers and will associate side set ids to
  // faces of the grandparent mesh.
  const auto& tet_store = m_refiner.tet_store;
  BndFaceSides grandfacesides;
  for (const auto& f : facetets) {
    auto nc = tet_store.data( f.second ).children.size();
    if (nc > 0) {  // if tet is refined
      auto ss = faceSides( facesets, f.first );    // query side sets for face
      for (decltype(nc) i=0; i<nc; ++i) {  // for all child tets
        auto childtet = tet_store.get_child_id( f.second, i );
        auto parent = tet_store.get_parent_id( childtet );
        if (parent != f.second) {  // detect backtrackers: 2->8 or 4->8
          // construct face of old mesh boundary face using local ids
          Face face{{ tk::cref_find( m_lid, f.first[0] ),
                      tk::cref_find( m_lid, f.first[1] ),
                      tk::cref_find( m_lid, f.first[2] ) }};
          auto par = parents( m_grand, face );  // find parent nodes of face
          if (par.size() == 3) {
            std::vector< std::size_t > g( begin(par), end(par) );
            grandfacesides[ {{ m_gid[g[0]], m_gid[g[1]], m_gid[g[2]] }} ] = ss;
          }
        }
      }
    }
  }

  return { facetets, grandfacesides, facesets };
}

std::unordered_set< int >
Refiner::faceSides( const BndFaceSets& facesets, const tk::UnsMesh::Face& t )
// *****************************************************************************
//  Return a set of side set ids for a face, given by 3 node ids
//! \param[in] facesets Boundary face sets to search in
//! \param[in] t Triangle face, given by 3 node ids
//! \return A unique set of side set ids in which the face is found or an empty
//!    set if the face was not found.
// *****************************************************************************
{
  std::unordered_set< int > ss;

  for (const auto& s : facesets)  // for all phsyical boundaries
    if (s.second.find(t) != end(s.second))
      ss.insert( s.first );

  return ss;
}

void
Refiner::newBndMesh( const std::unordered_set< std::size_t >& old,
                     const std::unordered_set< std::size_t >& ref )
// *****************************************************************************
// Update boundary data structures after mesh refinement
//! \param[in] old Unique nodes of the old (unrefined) mesh using local ids
//! \param[in] ref Unique nodes of the refined mesh using local ids
// *****************************************************************************
{
  // Generate boundary face data structures used to regenerate boundary face
  // and node data after mesh refinement
  auto bnd = boundary();

  // Regerate boundary faces and nodes after mesh refinement
  updateBndFaces( old, ref, bnd );
  updateBndNodes( old, ref, bnd );
}

void
Refiner::updateBndFaces( const std::unordered_set< std::size_t >& old,
                         const std::unordered_set< std::size_t >& ref,
                         const BndFaceData& bnd )
// *****************************************************************************
// Regenerate boundary faces after mesh refinement step
//! \param[in] old Unique nodes of the old (unrefined) mesh using local ids
//! \param[in] ref Unique nodes of the refined mesh using local ids
//! \param[in] bnd Tuple of boundary face data. See type definition.
// *****************************************************************************
{
  IGNORE(ref);  // to avoid compiler warning when asserts are optimized away

  using Face = tk::UnsMesh::Face;

  // storage for boundary faces associated to side-set IDs of the refined mesh
  decltype(m_bface) bface;              // will become m_bface
  // storage for boundary faces-node connectivity of the refined mesh
  decltype(m_triinpoel) triinpoel;      // will become m_triinpoel
  // face id counter
  std::size_t facecnt = 0;
  // will collect unique faces added for each side set
  BndFaceSets bndfacesets;

  // See Refiner::boundary()
  const auto& facetets = std::get< 0 >( bnd );
  const auto& grandfacesides = std::get< 1 >( bnd );
  const auto& facesets = std::get< 2 >( bnd );

  // Lambda to associate a boundary face and connectivity to a side set.
  // Argument 's' is the list of faces (ids) to add the new face to. Argument
  // 'ss' is the side set id to which the face is added. Argument 'f' is the
  // triangle face connectivity to add.
  auto addBndFace = [&]( std::vector< std::size_t >& s, int ss, const Face& f )
  {
    // only add face if it has not yet been aded to this side set
    if (bndfacesets[ ss ].insert( f ).second) {
      s.push_back( facecnt++ );
      triinpoel.insert( end(triinpoel), begin(f), end(f) );
    }
  };

  // Lambda to evaluate and optionally add a child face to a side set if its
  // grandparent was on the side set. Argument 'par' is the unique set of nodes
  // of the face to search for in the grandparent mesh using loca ids. Argument
  // 'childface' is the 3 local node ids of the child's face to add. Argument
  // 'ss' is the side set id to which the face is being added. Argument faces
  // is the list of faces to add to.
  auto evalGrandFace = [&]( const std::unordered_set< std::size_t >& par,
                            const Face& childface,
                            int ss,
                            std::vector< std::size_t >& faces )
  {
    Assert( par.size() == 3, "Node list to search for grandparent must be a"
            " face with exactly 3 nodes" );
    // Construct a vector from the set of nodes passed in
    std::vector< std::size_t > g( begin(par), end(par) );
    // Construct a face from the set of nodes passed in using global ids
    Face gf{{ m_gid[g[0]], m_gid[g[1]], m_gid[g[2]] }};
    // Attempt to find parents of nodes of the face in the grandparent mesh
    auto p = grandfacesides.find( gf );
    // If the face was found among the grandparents and the grandparent was on
    // the same side side set, add child face to the side set.
    if (p != end(grandfacesides) && p->second.find(ss) != end(p->second)) {
      std::cout << 'g';
      addBndFace( faces, ss, {{ m_gid[ childface[0] ], m_gid[ childface[1] ],
                                m_gid[ childface[2] ] }} );
    }
  };

  std::cout << thisIndex << ':';
  // Regenerate boundary faces after refinement step
  const auto& tet_store = m_refiner.tet_store;
  for (const auto& f : facetets) {  // for all boundary faces in old mesh
    // construct face of old mesh boundary face using local ids
    Face oldFace{{ tk::cref_find( m_lid, f.first[0] ),
                   tk::cref_find( m_lid, f.first[1] ),
                   tk::cref_find( m_lid, f.first[2] ) }};
    // for all side sets of the face, match children's faces to side sets
    for (const auto& ss : faceSides(facesets,f.first)) {
      // will associate to side set id of old (unrefined) mesh boundary face
      auto& faces = bface[ ss ];
      // query number of children of boundary tet adjacent to boundary face
      auto nc = tet_store.data( f.second ).children.size();
      if (nc == 0) {    // if boundary tet is not refined, add its boundary face
        addBndFace( faces, ss, f.first );
      } else {          // if boundary tet is refined
        const auto& tets = tet_store.tets;
        for (decltype(nc) i=0; i<nc; ++i ) {      // for all child tets
          // get child tet id
          auto childtet = tet_store.get_child_id( f.second, i );
          auto ct = tets.find( childtet );
          Assert( ct != end(tets), "Child tet not found" );
          // ensure all nodes of child tet are in refined mesh
          Assert( std::all_of( begin(ct->second), end(ct->second),
                    [&]( std::size_t n ){ return ref.find(n) != end(ref); } ),
                  "Boundary child tet node id not found in refined mesh" );
          // get nodes of child tet
          auto A = ct->second[0];
          auto B = ct->second[1];
          auto C = ct->second[2];
          auto D = ct->second[3];
          // form all 4 faces of child tet
          std::array<Face,4> face{{{{A,C,B}}, {{A,B,D}}, {{A,D,C}}, {{B,C,D}}}};
          // search all faces of all children and match them to side sets
          for (const auto& rf : face) {   // for all faces of child tet
            // find nodes of the parent face of child face
            auto par = parents( old, rf );
            if (par.size() == 3) {
              if (std::all_of( begin(oldFace), end(oldFace),
                    [&]( std::size_t n ){ return par.find(n) != end(par); } ))
              { // If all 3 nodes of the old mesh boundary face are the same as
                // those of the child's parent face, the child's face is also
                // on the same side set as the parent's face.
                addBndFace(faces,ss,{{m_gid[rf[0]],m_gid[rf[1]],m_gid[rf[2]]}});
              } else {
                // If the number of the child's nodes' parent nodes are 3 but
                // they are not all the 3 nodes of the old mesh boundary face,
                // the child face nodes may still match a boundary face of the
                // parent of the old (i.e., the grandparent) mesh.
                evalGrandFace( par, rf, ss, faces );
              }
            } else if (par.size() == 4) {
              // If the number of the child's node's parent nodes are 4, the
              // child face nodes may still match a boundary face of the parent
              // of the old (i.e., the grandparent) mesh, but only if number of
              // the nodes of the grandparent face of the child's face is 3.
              auto pf = parents( m_grand, par );
              if (pf.size() == 3) evalGrandFace( pf, rf, ss, faces );
            }
          }
        }
      }
    }
  }
  std::cout << '\n';

  // Update boundary face data structures
  m_bface = std::move(bface);
  m_triinpoel = std::move(triinpoel);

std::cout << thisIndex << " update: " << tk::sumvalsize(m_bface) << '\n';
}

void
Refiner::updateBndNodes( const std::unordered_set< std::size_t >& old,
                         const std::unordered_set< std::size_t >& ref,
                         const BndFaceData& bnd )
// *****************************************************************************
// Update boundary nodes after mesh refinement
//! \param[in] old Unique nodes of the old (unrefined) mesh using local ids
//! \param[in] ref Unique nodes of the refined mesh using local ids
//! \return Two maps, both associating a pair of side set id and adjacent tet id
//!   (value) to a partition-boundary triangle face (key) given by three global
//!   node IDs. The first map contains face -> tet id for the old (unrefined
//!   mesh) and the second one contains face -> tet id for parents of tets that
//!   underwent 2->8 or 4->8 refinement before the old mesh, i.e., backtrackers.
// *****************************************************************************
{
  IGNORE(ref);  // to avoid compiler warning when asserts are optimized away

  using Edge = tk::UnsMesh::Edge;

  // Generate unique set of nodes for each side set
  std::unordered_map< int, std::unordered_set< std::size_t > > bnodeset;
  for (const auto& ss : m_bnode)  // for all phsyical boundaries (sidesets)
    bnodeset[ ss.first ].insert( begin(ss.second), end(ss.second) );

  // Lambda to find the nodes of the parent edge of a node
  auto parentEdge = [ &old, this ]( std::size_t c ) -> Edge {
    if (old.find(c) != end(old)) // if node is in old mesh, return doubled input
      return {{ c, c }};
    else
      return this->m_refiner.node_connectivity.get( c );
  };

  // Lambda to find a global node ID among the nodelists of side sets. Return
  // all side set ids in which the node is found (or an empty vector if the
  // node was not found).
  auto bndNode = [ &bnodeset ]( std::size_t p ){
    std::vector< int > ss;
    for (const auto& s : bnodeset)  // for all phsyical boundaries (sidesets)
      if (s.second.find(p) != end(s.second))
        ss.push_back( s.first );
    return ss;
  };

  // storage for boundary nodes associated to side-set IDs of the refined mesh
  decltype(m_bnode) bnode;              // will become m_node

  const auto& facetets = std::get< 0 >( bnd );

  // Regenerate boundary node lists after refinement step
  for (const auto& f : facetets) {  // for all boundary faces in old mesh
    // query number of children of boundary tet adjacent to boundary face
    auto nc = m_refiner.tet_store.data( f.second ).children.size();
    if (nc == 0) {  // if boundary tet is not refined, add its boundary node
      for (auto n : f.first) {
        auto ss = bndNode( n );
        for (auto s : ss) bnode[ s ].push_back( n );
      }
    } else {        // if boundary tet is refined
      const auto& tets = m_refiner.tet_store.tets;
      for (decltype(nc) i=0; i<nc; ++i ) {      // for all child tets
        // get child tet id
        auto childtet = m_refiner.tet_store.get_child_id( f.second, i );
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
        // form all 6 edges of child tet
        std::array<Edge,6> edge{{ {{A,B}}, {{B,C}}, {{A,C}},
                                  {{A,D}}, {{B,D}}, {{C,D}} }};
        for (const auto& re : edge)
          for (auto c : re) {
            auto p = parentEdge( c );
            auto ss1 = bndNode( m_gid[p[0]] );
            auto ss2 = bndNode( m_gid[p[1]] );
            for (auto s1 : ss1)
              for (auto s2 : ss2)
                if (s1 == s2)
                  bnode[ s1 ].push_back( m_gid[c] );
          }
      }
    }
  }

  // Make boundary node IDs unique for each physical boundary (side set)
  for (auto& s : bnode) tk::unique( s.second );

  // Update boundary node lists
  m_bnode = std::move(bnode);
}

#include "NoWarning/refiner.def.h"
