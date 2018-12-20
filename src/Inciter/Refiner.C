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
                  const tk::CProxy_Solver& solver,
                  const Scheme& scheme,
                  const tk::RefinerCallback& cbr,
                  const tk::SorterCallback& cbs,
                  const std::vector< std::size_t >& ginpoel,
                  const tk::UnsMesh::CoordMap& coordmap,
                  const std::map< int, std::vector< std::size_t > >& belem,
                  const std::vector< std::size_t >& triinpoel,
                  const std::map< int, std::vector< std::size_t > >& bnode,
                  int nchare ) :
  m_host( transporter ),
  m_sorter( sorter ),
  m_solver( solver ),
  m_scheme( scheme ),
  m_schemeproxy(),
  m_cbr( cbr ),
  m_cbs( cbs ),
  m_ginpoel( ginpoel ),
  m_el( tk::global2local( ginpoel ) ),     // fills m_inpoel, m_gid, m_lid
  m_coordmap( coordmap ),
  m_coord( flatcoord(coordmap) ),
  m_belem( belem ),
  m_triinpoel( triinpoel ),
  m_bnode( bnode ),
  m_nchare( nchare ),
  m_initial( true ),
  m_t( 0.0 ),
  m_initref( g_inputdeck.get< tag::amr, tag::init >() ),
  m_ninitref( g_inputdeck.get< tag::amr, tag::init >().size() ),
  m_refiner( m_inpoel ),
  m_nref( 0 ),
  m_extra( 0 ),
  m_ch(),
  m_edgedata(),
  m_edgedataCh(),
  m_bndEdges(),
  m_u(),
  m_msum()
// *****************************************************************************
//  Constructor
//! \param[in] transporter Transporter (host) proxy
//! \param[in] sorter Mesh reordering (sorter) proxy
//! \param[in] solver Linear system solver proxy
//! \param[in] scheme Discretization scheme
//! \param[in] cbr Charm++ callbacks for Refiner
//! \param[in] cbs Charm++ callbacks for Sorter
//! \param[in] ginpoel Mesh connectivity (this chare) using global node IDs
//! \param[in] coordmap Mesh node coordinates (this chare) for global node IDs
//! \param[in] belem File-internal elem ids of side sets (caller PE)
//! \param[in] triinpoel Triangle face connectivity with global IDs (caller PE)
//! \param[in] bnode Node lists of side sets (caller PE)
//! \param[in] nchare Total number of refiner chares (chare array elements)
// *****************************************************************************
{
  Assert( !m_ginpoel.empty(), "No elements assigned to refiner chare" );

  usesAtSync = true;    // Enable migration at AtSync

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
Refiner::dtref( tk::real t,
                const SchemeBase::Proxy& scheme,
                const std::map< int, std::vector< std::size_t > >& bnode )
// *****************************************************************************
// Start mesh refinement (during time stepping, t>0)
//! \param[in] t Physical time
//! \param[in] s Discretization scheme Charm++ proxy we interoperate with
//! \param[in] bnode Node lists of side sets
// *****************************************************************************
{
  m_initial = false;
  m_t = t;

  // Update boundary node lists
  m_bnode = bnode;

  // Store discretization scheme proxy
  m_schemeproxy = scheme;

  t0ref();
}

void
Refiner::t0ref()
// *****************************************************************************
// Start new step of initial mesh refinement (before t>0)
// *****************************************************************************
{
  if (m_initial) {
    Assert( m_ninitref > 0, "No initial mesh refinement steps configured" );
    // Output initial mesh to file
    auto l = m_ninitref - m_initref.size();  // num initref steps completed
    auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();
    if (l == 0) writeMesh( "t0ref", l, t0-1.0 );
  }

  m_extra = 0;
  m_bndEdges.clear();
  m_ch.clear();
  m_edgedataCh.clear();

  updateEdgeData();

  // Generate boundary edges
  bndEdges();
}

void
Refiner::writeMesh( const std::string& prefix, uint64_t it, tk::real t )
// *****************************************************************************
// Write mesh to file
//! \param[in] prefix Prefix to mesh filename
//! \param[in] it Iteration number counting number of meshes changed during
//!   multiple file output. (This can be the refinement level for initial (t<0)
//!   AMR, or ...
//! \param[in] t Physical time ...
// *****************************************************************************
{
  tk::ExodusIIMeshWriter
    mw( prefix + ".e-s"
        + '.' + std::to_string( it )       // create new file for new mesh
        + '.' + std::to_string( m_nchare )   // total number of workers
        + '.' + std::to_string( thisIndex ), // new file per worker
        tk::ExoWriter::CREATE );

  // Prepare boundary data for file output
  decltype(m_bnode) bnode;
  if (m_nchare == 1) {  // do not write boundary data in parallel
    // Convert boundary node lists to local ids for output
    bnode = m_bnode;
    for (auto& s : bnode) for (auto& p : s.second) p = tk::cref_find(m_lid,p);
  }

  // Output mesh
  tk::UnsMesh refmesh( m_inpoel, m_coord, bnode );
  mw.writeMesh( refmesh );

  // Output element-centered scalar fields on recent mesh refinement step
  mw.writeTimeStamp( 1, t );
  auto& tet_store = m_refiner.tet_store;
  mw.writeElemVarNames( { "refinement level", "cell type" } );
  mw.writeElemScalar( 1, 1, tet_store.get_refinement_level_list() );
  mw.writeElemScalar( 1, 2, tet_store.get_cell_type_list() );
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

   refine();
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

    //errorRefine();
    uniformRefine();

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
      thisProxy[ c ].addRefBndEdges( thisIndex, m_edgedata );
    }
  }
}

void
Refiner::addRefBndEdges( int fromch, const AMR::EdgeData& ed )
// *****************************************************************************
//! Receive tagged edges on our chare boundary
//! \param[in] fromch Chare call coming from
//! \param[in] ed Tagged edges on chare boundary
// *****************************************************************************
{
  // Save/augment buffer of edge data for each sender chare
  m_edgedataCh[ fromch ].insert( begin(ed), end(ed) );
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
  auto unlocked = AMR::Edge_Lock_Case::unlocked;

  // Storage for edge data that need correction to yield a conforming mesh
  AMR::EdgeData extra;

  // loop through all edges shared with other chares
  for (const auto& c : m_edgedataCh)       // for all chares we share edges with
    for (const auto& r : c.second) {       // for all edges shared with c.first
      // find local data of remote edge { r.first[0] - r.first[1] }
      auto it = m_edgedata.find( r.first );
      if (it != end(m_edgedata)) {

        const auto& local = it->second;
        const auto& remote = r.second;
        auto local_needs_refining = local.first;
        auto local_lock_case = local.second;
        auto remote_needs_refining = remote.first;
        auto remote_lock_case = remote.second;

        auto local_needs_refining_orig = local_needs_refining;
        auto local_lock_case_orig = local_lock_case;
        auto remote_needs_refining_orig = remote_needs_refining;
        auto remote_lock_case_orig = remote_lock_case;

        Assert( !(local_lock_case > unlocked && local_needs_refining),
                "Invalid local edge: locked & needs refining" );
        Assert( !(remote_lock_case > unlocked && remote_needs_refining),
                "Invalid remote edge: locked & needs refining" );

        if (remote_lock_case > unlocked) remote_lock_case = unlocked;

        // compute lock from local and remote locks as most restrictive
        local_lock_case = std::max( local_lock_case, remote_lock_case );

        if (local_lock_case > unlocked) local_needs_refining = false;

        if (local_lock_case == unlocked && remote_needs_refining)
          local_needs_refining = true;

        if (local_lock_case != local_lock_case_orig ||
            local_needs_refining != local_needs_refining_orig) {

           auto l1 = tk::cref_find( m_lid, r.first[0] );
           auto l2 = tk::cref_find( m_lid, r.first[1] );
           Assert( l1 != l2, "Edge end-points local ids are the same" );

           std::cout << thisIndex << " correcting "
             << r.first[0] << '-' << r.first[1] << ": l{"
             << local_needs_refining_orig << ',' << local_lock_case_orig
             << "} + r{" << remote_needs_refining_orig << ','
             << remote_lock_case_orig << "} -> {" << local_needs_refining
             << ',' << local_lock_case << "}\n";

           extra[ {{ std::min(l1,l2), std::max(l1,l2) }} ] =
             { local_needs_refining, local_lock_case };
         }

      }
    }

  m_extra = extra.size();

  correctRefine( extra );

  // Aggregate number of extra edges that still need correction
  contribute( sizeof(std::size_t), &m_extra, CkReduction::max_ulong,
              m_cbr.get< tag::matched >() );
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
    m_refiner.error_refinement_corr( extra );
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
  for (const auto& e : ref_edges) {
    const auto& ed = e.first.get_data();
    const auto ged = Edge{{ m_gid[ ed[0] ], m_gid[ ed[1] ] }};
    m_edgedata[ ged ] = { e.second.needs_refining, e.second.lock_case };
  }
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
  m_refiner.perform_refinement();
  updateMesh();

  AtSync();   // Migrate here if needed

  if (m_initial) {      // if initial (before t=0) AMR

    // Output mesh after recent step of initial mesh refinement
    auto l = m_ninitref - m_initref.size() + 1;  // num initref steps completed
    auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();
    // Generate times for file output equally subdividing t0-1...t0 to
    // m_ninitref steps
    tk::real t = t0 - 1.0 +
                 static_cast<tk::real>(l)/static_cast<tk::real>(m_ninitref);
    writeMesh( "t0ref", l, t );
    // Remove initial mesh refinement step from list
    if (!m_initref.empty()) m_initref.pop_back();
    // Continue to next initial AMR step or finish
    if (!m_initref.empty()) t0ref(); else endt0ref();

  } else {              // if AMR during time stepping (t>0)

    // Output mesh after recent step of mesh refinement during time stepping
    writeMesh( "dtref", 0, m_t );

//     // Augment node comm. map with newly added nodes on chare-boundary edges
//     for (const auto& c : m_edgedataCh) {
//       auto& nodes = tk::ref_find( m_msum, c.first );
//       for (const auto& n : c.second)
//         nodes.push_back( n.second.first );
//     }

    // Send new mesh and solution back to PDE worker
    Assert( m_scheme.get()[thisIndex].ckLocal() != nullptr,
            "About to use nullptr" );
    auto e = tk::element< SchemeBase::ProxyElem >( m_schemeproxy, thisIndex );
    boost::apply_visitor( Resize(m_el,m_coord,m_u,m_msum,m_bnode), e );

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
  m_sorter[ thisIndex ].insert( m_host, m_solver, m_cbs, m_scheme, m_ginpoel,
    m_coordmap, m_belem, m_triinpoel, m_bnode, m_nchare, CkMyPe() );

  // Compute final number of cells across whole problem
  std::size_t nelem = m_ginpoel.size()/4;
  contribute( sizeof(std::size_t), &nelem, CkReduction::sum_ulong,
              m_cbr.get< tag::refined >() );
}

void
Refiner::uniformRefine()
// *****************************************************************************
// Do uniform mesh refinement
// *****************************************************************************
{
  // Do uniform refinement
  m_refiner.uniform_refinement();

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
  auto npoin = tk::npoin( m_inpoel );
  // Generate edges surrounding points in old mesh
  auto esup = tk::genEsup( m_inpoel, 4 );
  auto psup = tk::genPsup( m_inpoel, 4, esup );

  // Get solution whose error to evaluate
  tk::Fields u;
  if (m_initial) {      // if initial (before t=0) AMR

    // Evaluate initial conditions at mesh nodes
    u = nodeinit( npoin, esup );

  } else {              // if AMR during time stepping (t>0)

    // Get old solution from worker (pointer to soln from bound array element)
    Assert( m_scheme.get()[thisIndex].ckLocal() != nullptr,
            "About to use nullptr" );
    auto e = tk::element< SchemeBase::ProxyElem >( m_schemeproxy, thisIndex );
    boost::apply_visitor( Solution(u), e );

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

  // Do error-based refinement
  m_refiner.error_refinement( edge );

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
    auto npoin = tk::npoin( m_inpoel );
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

    std::cout << thisIndex << ": " << edge.size() << std::endl;

    if (!useredges.empty()) {
      std::cout << "Edges tagged but not found on chare " << thisIndex << ": ";
      for (const auto& e : useredges) std::cout << e[0] << ',' << e[1] << ' ';
    }
    std::cout << std::endl;

    // Do error-based refinement
    m_refiner.error_refinement( edge );

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
    auto npoin = tk::npoin( m_inpoel );
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
    m_refiner.error_refinement( edge );

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
  if (centering == ctr::Centering::NODE) {

    // Node-centered: evaluate ICs for all scalar components integrated
    for (const auto& eq : g_cgpde) eq.initialize( m_coord, u, t0 );

  } else if (centering == ctr::Centering::ELEM) {

    // Initialize cell-based unknowns
    tk::Fields ue( m_inpoel.size()/4, nprop );
    auto lhs = ue;
    auto geoElem = tk::genGeoElemTet( m_inpoel, m_coord );
    for (const auto& eq : g_dgpde)
      eq.lhs( geoElem, lhs );
    for (const auto& eq : g_dgpde)
      eq.initialize( lhs, m_inpoel, m_coord, ue, t0 );

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

  // Generate unique node lists of old and refined mesh using local ids
  std::unordered_set< std::size_t > old( m_inpoel.cbegin(), m_inpoel.cend() );
  std::unordered_set< std::size_t > ref( refinpoel.cbegin(), refinpoel.cend() );

  // Get old solution from worker (pointer to soln from bound array element)
  if (!m_initial) {
    Assert( m_scheme.get()[thisIndex].ckLocal() != nullptr,
            "About to use nullptr" );
    auto e = tk::element< SchemeBase::ProxyElem >( m_schemeproxy, thisIndex );
    boost::apply_visitor( Solution(m_u), e );
    // Get nodal communication map from Discretization worker
    m_msum = m_scheme.get()[thisIndex].ckLocal()->Msum();
  }

  // Update mesh connectivity with local node IDs
  m_inpoel = m_refiner.tet_store.get_active_inpoel();

  // Update mesh and solution after refinement
  newVolMesh( old, ref );
  newBndMesh( old, ref );

  // Update mesh connectivity with new global node ids
  m_ginpoel = m_inpoel;
  Assert( tk::uniquecopy(m_ginpoel).size() == m_coord[0].size(),
          "Size mismatch" );
  for (auto& i : m_ginpoel) i = m_gid[i];

  // Update flat coordinates storage
  m_coord = flatcoord( m_coordmap );

  // Ensure valid mesh after refinement
  Assert( tk::positiveJacobians( m_inpoel, m_coord ),
          "Refined mesh cell Jacobian non-positive" );

  // Perform leak test on new mesh
  Assert( !tk::leakyPartition(
            tk::genEsuelTet( m_inpoel, tk::genEsup(m_inpoel,4) ),
            m_inpoel, m_coord ),
          "Refined mesh partition leaky" );
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

  // Resize node coordinates, global ids, and solution vector
  auto npoin = ref.size();
  auto nprop = g_inputdeck.get<tag::component>().nprop();
  x.resize( npoin );
  y.resize( npoin );
  z.resize( npoin );
  m_gid.resize( npoin, std::numeric_limits< std::size_t >::max() );
  if (!m_initial) m_u.resize( npoin, nprop );

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
        // generate solution for newly added node
        if (!m_initial)   // only during t>0 refinement
          for (std::size_t c=0; c<nprop; ++c)
            m_u(r,c,0) = (m_u(p[0],c,0) + m_u(p[1],c,0))/2.0;
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
}

Refiner::BndFaces
Refiner::boundary()
// *****************************************************************************
//  Generate boundary data structure used to update refined boundary faces and
//  nodes
//! \return Map associating a pair of side set id and adjacent tet id (value) to
//!   a partition-boundary triangle face (key) given by three global node IDs.
//! \details The output of this function is used to regenerate physical boundary
//!   face and node data structures after refinement, see updateBndFaces() and
//!   updateBndNodes().
// *****************************************************************************
{
  using Face = tk::UnsMesh::Face;
  using Tet = tk::UnsMesh::Tet;

  // Generate the inverse of AMR's active inpoel
  std::unordered_map< Tet,
                      std::size_t,
                      tk::UnsMesh::Hash<4>,
                      tk::UnsMesh::Eq<4> > invtets;
  //for (const auto& t : m_refiner.tet_store.tets) invtets[ t.second ] = t.first;
  const auto& refinpoel = m_refiner.tet_store.get_active_inpoel();
  for (std::size_t e=0; e<refinpoel.size()/4; ++e) {
    invtets[ {{ refinpoel[e*4+0], refinpoel[e*4+1],
                refinpoel[e*4+2], refinpoel[e*4+3] }} ] = e;
  }

  // Generate data structure that associates the pair of side set id and
  // adjacent tet id to a boundary triangle face for all boundary faces. After
  // this loop we will have all tets adjacent to boundary faces, where the
  // boundary includes all physical boundaries (where the user may assign
  // boundary conditions via side sets from the mesh file) as well as
  // chare-boundaries due to domain decomposition.
  BndFaces bnd;
  auto oldesuel = tk::genEsuelTet( m_inpoel, tk::genEsup(m_inpoel,4) );
  for (std::size_t e=0; e<oldesuel.size()/4; ++e) {    // for tets on this chare
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f) {
      if (oldesuel[mark+f] == -1) {  // if a face does not have an adjacent tet
        Face b{{ m_ginpoel[ mark+tk::lpofa[f][0] ],
                 m_ginpoel[ mark+tk::lpofa[f][1] ],
                 m_ginpoel[ mark+tk::lpofa[f][2] ] }};
        Tet t{{ m_inpoel[mark+0], m_inpoel[mark+1],
                m_inpoel[mark+2], m_inpoel[mark+3] }};
        // associate tet id adjacent to boundary face to boundary face, at this
        // point we fill in -1 for the side set id, since we don't know it yet
        auto i = invtets.find( t );
        Assert( i != end(invtets), "Tet not found in AMR object: " +
                std::to_string(t[0]) + ',' + std::to_string(t[1]) + ',' +
                std::to_string(t[2]) + ',' + std::to_string(t[3]) );
        bnd[ b ] = { -1, i->second };
      }
    }
  }

  // Assign side set ids to faces on the physical boundary only. Note that
  // m_belem may be empty and that is okay, in that case bnd will contain data
  // on both chare as well as physical boundaries and the side set id in the map
  // value will stay at -1.
  for (const auto& ss : m_belem)  // for all phsyical boundaries (sidesets)
    for (auto f : ss.second) {    // for all faces on this physical boundary
      Face t{{ m_triinpoel[f*3+0], m_triinpoel[f*3+1], m_triinpoel[f*3+2] }};
      auto cf = bnd.find( t );
      if (cf != end(bnd)) cf->second.first = ss.first;
    }

  // Remove chare-boundary faces (to which the above loop did not assign set
  // id), but only if the above loop was run, i.e., if m_belem was not empty.
  // This results in removing the chare-boundary faces keeping only the physical
  // boundary faces (and their associated side set id and tet id).
  if (!m_belem.empty()) {
    auto bndcopy = bnd;
    for (const auto& f : bndcopy)
      if (f.second.first == -1)
        bnd.erase( f.first );
  }

  return bnd;
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
  // generate boundary face data structure used to regenerate boundary face and
  // node data after mesh refinement
  auto bnd = boundary();

  // regerate boundary faces and nodes after mesh refinement
  updateBndFaces( old, ref, bnd );
  updateBndNodes( old, ref, bnd );
}

void
Refiner::updateBndFaces( const std::unordered_set< std::size_t >& old,
                         const std::unordered_set< std::size_t >& ref,
                         const BndFaces& bnd )
// *****************************************************************************
// Regenerate boundary faces after mesh refinement step
//! \param[in] old Unique nodes of the old (unrefined) mesh using local ids
//! \param[in] ref Unique nodes of the refined mesh using local ids
//! \param[in] bnd Map associating a pair of side set id and adjacent tet id
//!   (value) to a partition-boundary triangle face (key) given by three global
//!   node IDs.
// *****************************************************************************
{
  IGNORE(ref);  // to avoid compiler warning when asserts are optimized away

  using Face = tk::UnsMesh::Face;

  // storage for boundary faces associated to side-set IDs of the refined mesh
  decltype(m_belem) belem;              // will become m_belem
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
    triinpoel.insert( end(triinpoel), begin(f), end(f) );
  };

  // Lambda to find the nodes of the parent face of a child face. Argument
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
    // Ensure the number of parent nodes are correct
    Assert( s.size() == 3 || s.size() == 4,
            "Parent node list length must be 3 or 4" );
    // Ensure all parent nodes are part of the old (unrefined) mesh
    Assert( std::all_of( begin(s), end(s), [ &old ]( std::size_t j ){
                           return old.find(j) != end(old); } ),
            "Old face nodes not in old mesh" );
    // Return unique set of nodes of the parent face of child face
    return s;
  };

  // Regenerate boundary faces after refinement step
  for (const auto& f : bnd)     // for all boundary faces
    if (f.second.first != -1) { // for all physical boundary faces
      // construct face of old mesh boundary face using local ids
      Face oldface{{ tk::cref_find( m_lid, f.first[0] ),
                     tk::cref_find( m_lid, f.first[1] ),
                     tk::cref_find( m_lid, f.first[2] ) }};
      // will associate to side set id of old (unrefined) mesh boundary face
      auto& sideface = belem[ f.second.first ];
      // query number of children of boundary tet adjacent to boundary face
      auto nc = m_refiner.tet_store.data( f.second.second ).num_children;
      if (nc == 0) {    // if boundary tet is not refined, add its boundary face
        addBndFace( sideface, f.first );
      } else {          // if boundary tet is refined
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
          std::array<Face,4> face{{{{A,C,B}}, {{A,B,D}}, {{A,D,C}}, {{B,C,D}}}};
          for (const auto& rf : face) {   // for all faces of child tet
            // find nodes of the parent face of child face
            auto parfac = parentFace( rf );
            // if the child shares a face with its parent and all 3 nodes of the
            // parent face of the child's face are the same, the child face is
            // on the same boundary as the parent face
            if ( parfac.size() == 3 &&
                 parfac.find(oldface[0]) != end(parfac) &&
                 parfac.find(oldface[1]) != end(parfac) &&
                 parfac.find(oldface[2]) != end(parfac) )
            {
              addBndFace(sideface, {{m_gid[rf[2]],m_gid[rf[1]],m_gid[rf[0]]}});
            }
          }
        }
      }
    }

  // Update boundary face data structures
  m_belem = std::move(belem);
  m_triinpoel = std::move(triinpoel);
}

void
Refiner::updateBndNodes( const std::unordered_set< std::size_t >& old,
                         const std::unordered_set< std::size_t >& ref,
                         const Refiner::BndFaces& bnd )
// *****************************************************************************
// Update boundary nodes after mesh refinement
//! \param[in] old Unique nodes of the old (unrefined) mesh using local ids
//! \param[in] ref Unique nodes of the refined mesh using local ids
//! \param[in] bnd Map associating a pair of side set id and adjacent tet id
//!   (value) to a partition-boundary triangle face (key) given by three global
//!   node IDs.
// *****************************************************************************
{
  IGNORE(ref);  // to avoid compiler warning when asserts are optimized away

  using Edge = tk::UnsMesh::Edge;

  // storage for boundary nodes associated to side-set IDs of the refined mesh
  decltype(m_bnode) bnode;              // will become m_node

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
  // node was not found.
  auto bndNode = [ &bnodeset ]( std::size_t p ){
    std::vector< int > ss;
    for (const auto& s : bnodeset)  // for all phsyical boundaries (sidesets)
      if (s.second.find(p) != end(s.second))
        ss.push_back( s.first );
    return ss;
  };


  // Regenerate boundary node lists after refinement step
  for (const auto& f : bnd) {
    // query number of children of boundary tet adjacent to boundary face
    auto nc = m_refiner.tet_store.data( f.second.second ).num_children;
    if (nc == 0) {  // if boundary tet is not refined, add its boundary node
      for (auto n : f.first) {
        auto ss = bndNode( n );
        for (auto s : ss) bnode[ s ].push_back( n );
      }
    } else {        // if boundary tet is refined
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
