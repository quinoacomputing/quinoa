// *****************************************************************************
/*!
  \file      src/Inciter/Refiner.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Mesh refiner for interfacing the mesh refinement library
  \see       Refiner.h for more info.
*/
// *****************************************************************************

#include <algorithm>

#include "Refiner.h"
#include "Reorder.h"
#include "AMR/mesh_adapter.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "CGPDE.h"
#include "DGPDE.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< CGPDE > g_cgpde;
extern std::vector< DGPDE > g_dgpde;

} // inciter::

using inciter::Refiner;

Refiner::Refiner( const CProxy_Transporter& transporter,
                  const tk::RefinerCallback& cbr,
                  const std::vector< std::size_t >& ginpoel,
                  const tk::UnsMesh::CoordMap& coordmap ) :
  m_host( transporter ),
  m_cbr( cbr ),
  m_el( tk::global2local( ginpoel ) ),     // fills m_inpoel, m_gid, m_lid
//   m_nedge( 0 ),
//   m_nref( 0 ),
//   m_extra( 1 ),
  m_initref( g_inputdeck.get< tag::amr, tag::init >() )
// *****************************************************************************
//  Constructor
// *****************************************************************************
{
  // Reverse initial mesh refinement type list (will pop from back)
  std::reverse( begin(m_initref), end(m_initref) );

//   if ( !g_inputdeck.get< tag::amr, tag::init >().empty() ||
//        !g_inputdeck.get< tag::amr, tag::edge >().empty() )
//     partref();          // if initial mesh refinement configured, partition
//   else
//     finishref();        // if not, continue
}

// void
// Partitioner::partref()
// // *****************************************************************************
// //  Partition the mesh to NumPes partitions (before an initial refinement step)
// //! \details This function calls the mesh partitioner to partition (or
// //!   re-partition) the (current) mesh as a first step for an initial mesh
// //!   refinement step. The number of partitions always equals the numher of PEs.
// // *****************************************************************************
// {
//   // Generate element IDs for Zoltan
//   std::vector< long > gelemid( m_inpoel.size()/4 );
//   std::iota( begin(gelemid), end(gelemid), 0 );
// 
//   // Partition the mesh using Zoltan to number of PEs parts
//   const auto alg = g_inputdeck.get< tag::selected, tag::partitioner >();
//   const auto pel = tk::zoltan::geomPartMesh( alg,
//                                              centroids( m_inpoel, m_coord ),
//                                              gelemid,
//                                              CkNumPes() );
// 
//   Assert( pel.size() == gelemid.size(), "Size of ownership array (PE of "
//           "elements) after mesh partitioning does not equal the number of mesh "
//           "graph elements" );
// 
//   // Prepare for a step of initial mesh refinement
//   m_extra = 1;
//   m_bndEdges.clear();
//   m_pe.clear();
//   m_edgenodePe.clear();
//   m_edgenode.clear();
//   m_coordmap.clear();
//   auto g = m_ginpoel;
//   m_ginpoel.clear();
// 
//   // Categorize mesh elements (given by their gobal node IDs) by target PE and
//   // distribute to their PEs based on mesh partitioning.
//   distributePe( categorize( pel, g ) );
// }
// 
// void
// Partitioner::distributePe(
//   std::unordered_map< int, std::vector< std::size_t > >&& elems )
// // *****************************************************************************
// // Distribute mesh to their PEs during initial mesh refinement
// //! \param[in] elems Mesh cells (with global node IDs) categorized by target PEs
// // *****************************************************************************
// {
//   // Store own mesh connectivity
//   auto i = elems.find( CkMyPe() );
//   if (i != end(elems)) {
//     // Store our mesh chunk. The receive side also writes to this, so concat.
//     m_ginpoel.insert( end(m_ginpoel), begin(i->second), end(i->second) );
//     // store coordinates associated to global nodes of our mesh chunk
//     auto cm = coordmap( i->second );
//     // the receive side also writes to this, so concatenate
//     m_coordmap.insert( begin(cm), end(cm) );
//     // remove our mesh chunk from list (the rest will be exported)
//     elems.erase( i );
//   }
// 
//   m_nedge = 0;
// 
//   // Export connectivities to other PEs
//   if (elems.empty())
//     contribute( m_cb.get< tag::refdistributed >() );
//   else {
//     m_ndist = 0;
//     m_npeDist = elems.size();
//     for (const auto& p : elems)
//       thisProxy[ p.first ].addPeMesh( CkMyPe(), p.second, coordmap(p.second) );
//   }
// }

// void
// Partitioner::recvPeMesh()
// // *****************************************************************************
// //  Acknowledge received mesh chunk and its nodes
// // *****************************************************************************
// {
//   if (++m_ndist == m_npeDist) contribute( m_cb.get< tag::refdistributed >() );
// }
// 
// void
// Partitioner::bndEdges()
// // *****************************************************************************
// // Generate boundary edges and send them to all PEs
// //! \details This step happens when the mesh chunk on this PE has been
// //!   distributed after partitioning during an initial mesh refinement step. At
// //!   this point we have a contiguous chunk of the mesh on this PE as
// //!   determined by the partitioner. The next step is to extract the edges on
// //!   the boundary only. The boundary edges (shared by multiple PEs) will be
// //!   agreed on a refinement that yields a conforming mesh across PE boundaries.
// // *****************************************************************************
// {
//   Assert( !m_ginpoel.empty(), "No elements are assigned to PE" );
// 
//   // Compute local data from global mesh connectivity (m_inpoel, m_gid, m_lid)
//   m_el = tk::global2local( m_ginpoel );
// 
//   // Generate boundary edges of our mesh chunk
//   tk::UnsMesh::EdgeSet bnded;
//   auto esup = tk::genEsup( m_inpoel, 4 );         // elements surrounding points
//   auto esuel = tk::genEsuelTet( m_inpoel, esup ); // elems surrounding elements
//   for (std::size_t e=0; e<esuel.size()/4; ++e) {
//     auto mark = e*4;
//     for (std::size_t f=0; f<4; ++f) {
//       if (esuel[mark+f] == -1) {
//         auto A = m_gid[ m_inpoel[ mark+tk::lpofa[f][0] ] ];
//         auto B = m_gid[ m_inpoel[ mark+tk::lpofa[f][1] ] ];
//         auto C = m_gid[ m_inpoel[ mark+tk::lpofa[f][2] ] ];
//         bnded.insert( {{{A,B}}} );
//         bnded.insert( {{{B,C}}} );
//         bnded.insert( {{{C,A}}} );
//         Assert( m_lid.find( A ) != end(m_lid), "Local node ID not found" );
//         Assert( m_lid.find( B ) != end(m_lid), "Local node ID not found" );
//         Assert( m_lid.find( C ) != end(m_lid), "Local node ID not found" );
//       }
//     }
//   }
// 
//   // Export boundary edges to all PEs
//   thisProxy.addBndEdges( CkMyPe(), bnded );
// }
// 
// void
// Partitioner::addBndEdges( int frompe, const tk::UnsMesh::EdgeSet& ed )
// // *****************************************************************************
// //! Receive boundary edges from all PEs (including this one)
// //! \param[in] frompe PE call coming from
// //! \param[in] ed Edges on frompe's boundary (with global node IDs)
// // *****************************************************************************
// {
//   // Store incoming boundary edges
//   m_bndEdges[ frompe ].insert( begin(ed), end(ed) );
// 
//   if (++m_nedge == static_cast<std::size_t>(CkNumPes())) {
//     // Compute unique set of PEs that share at least a single edge with this PE
//     const auto& ownedges = tk::cref_find( m_bndEdges, CkMyPe() );
//     for (const auto& p : m_bndEdges)    // for all PEs
//       if (p.first != CkMyPe())          // for all PEs other than this PE
//         for (const auto& e : p.second)  // for all boundary edges
//           if (ownedges.find(e) != end(ownedges))
//             m_pe.insert( p.first );     // if edge is shared, store its PE
// 
//     refine();
//   }
// }
// 
// void
// Partitioner::refine()
// // *****************************************************************************
// //  Do a single step of initial mesh refinement based on user-input
// //! \details This is a single step in a potentially multiple-entry list of
// //!   initial adaptive mesh refinement steps. Distribution of the PE-boundary
// //!   edges has preceded this step, so that boundary edges (shared by multiple
// //!   PEs) can agree on a refinement that yields a conforming mesh across PE
// //!   boundaries.
// // *****************************************************************************
// {
//   // Convert node coordinates associated to global node IDs to a flat vector
//   auto npoin = m_coordmap.size();
//   Assert( m_gid.size() == npoin, "Size mismatch" );
//   m_coord[0].resize( npoin );
//   m_coord[1].resize( npoin );
//   m_coord[2].resize( npoin );
//   for (const auto& c : m_coordmap) {
//     auto i = tk::cref_find( m_lid, c.first );
//     Assert( i < npoin, "Indexing out of coordinate map" );
//     m_coord[0][i] = c.second[0];
//     m_coord[1][i] = c.second[1];
//     m_coord[2][i] = c.second[2];
//   }
// 
//   for (const auto& e : tk::cref_find(m_bndEdges,CkMyPe())) {
//     IGNORE(e);
//     Assert( m_lid.find( e[0] ) != end( m_lid ) &&
//             m_lid.find( e[1] ) != end( m_lid ),
//             "Boundary edge not found before refinement" );
//   }
// 
//   // Refine mesh based on next initial refinement type
//   if (!m_initref.empty()) {
//     auto r = m_initref.back();    // consume (reversed) list from back
//     if (r == ctr::AMRInitialType::UNIFORM)
//       uniformRefine();
//     else if (r == ctr::AMRInitialType::INITIAL_CONDITIONS)
//       errorRefine();
//     else Throw( "Initial AMR type not implemented" );
//   }
// 
//   // Additionally refine mesh based on user explicitly tagging edges
//   userRefine();
// 
//   for (const auto& e : tk::cref_find(m_bndEdges,CkMyPe())) {
//     IGNORE(e);
//     Assert( m_lid.find( e[0] ) != end( m_lid ) &&
//             m_lid.find( e[1] ) != end( m_lid ),
//             "Boundary edge not found after refinement" );
//   }
// 
//   // Ensure valid mesh after refinement
//   Assert( tk::positiveJacobians( m_inpoel, m_coord ),
//           "Refined mesh cell Jacobian non-positive" );
// 
//   // Export added nodes on our mesh chunk boundary to other PEs
//   if (m_pe.empty())
//     contribute( sizeof(std::size_t), &m_extra, CkReduction::max_int,
//                 m_cb.get< tag::matched >() );
//   else {
//     m_nref = 0;
//     for (auto p : m_pe) {       // for all PEs we share at least an edge with
//       // For all boundary edges of PE p, find out if we have added a new
//       // node to it, and if so, export parents->(newid,coords) to p.
//       tk::UnsMesh::EdgeNodeCoord exp;
//       for (const auto& e : tk::cref_find(m_bndEdges,p)) {
//         auto i = m_edgenode.find(e);
//         if (i != end(m_edgenode)) exp[ e ] = i->second;
//       }
//       thisProxy[ p ].addRefBndEdges( CkMyPe(), exp );
//     }
//   }
// }
// 
// void
// Partitioner::addRefBndEdges( int frompe, const tk::UnsMesh::EdgeNodeCoord& ed )
// // *****************************************************************************
// //! Receive newly added mesh node IDs on our PE boundary
// //! \param[in] frompe PE call coming from
// //! \param[in] ed Newly added node IDs associated to parent nodes on PE boundary
// //! \details Receive newly added global node IDs and coordinates associated to
// //!   global parent IDs of edges on our mesh chunk boundary.
// // *****************************************************************************
// {
//   // Save/augment buffer of edge-node (IDs and coords) categorized by sender PE
//   m_edgenodePe[ frompe ].insert( begin(ed), end(ed) );
//   // Acknowledge receipt of PE-boundary edges to sender
//   thisProxy[ frompe ].recvRefBndEdges();
// }
// 
// void
// Partitioner::recvRefBndEdges()
// // *****************************************************************************
// //  Acknowledge received newly added node IDs to edges shared among multiple PEs
// // *****************************************************************************
// {
//   // When we have heard from all PEs we share at least a single edge with,
//   // contribute the number of extra edges that this mesh refinement step has
//   // found that were not refined by this PE but were refined by other PEs this
//   // PE shares the edge with. A global maximum will then be computed on the
//   // number of extra edges appearing in Transporter::matched() and that is then
//   // used to decide if a new correction step is needed. If this is called for
//   // the first time in a given initial mesh refinement step, i.e., not after a
//   // correction step, m_extra = 1 on all PEs, so a correction step is assumed
//   // to be required.
//   if (++m_nref == m_pe.size()) {
//     contribute( sizeof(std::size_t), &m_extra, CkReduction::max_int,
//                 m_cb.get< tag::matched >() );
//   }
// }
// 
// void
// Partitioner::correctref()
// // *****************************************************************************
// //  Correct refinement to arrive at a conforming mesh across PE boundaries
// //! \details This function is called repeatedly until there is not a a single
// //!    edge that needs correction for the whole distributed problem to arrive at
// //!    a conforming mesh across PE boundaries during this initial mesh
// //!    refinement step.
// // *****************************************************************************
// {
//   // Storage for edges that still need a new node to yield a conforming mesh
//   tk::UnsMesh::EdgeSet extra;
// 
//   // Ensure that the same global node ID has been assigned by all PEs and that
//   // the new nodes have the same coordinates generated by potentially multiple
//   // PEs sharing the refined edge. This is done by searching for all edges that
//   // we share with and refined by other PEs: (1) If the incoming edge is found
//   // among our refined ones, we ensure the newly assigned global IDs equal
//   // (independently assigned by multiple PEs) and also that the new node
//   // coordinates equal to machine precision. (2) If the incoming edge is not
//   // found among our refined ones, we need to correct the mesh to make it
//   // conforming since the edge has been refined by the remote PE. We collect
//   // these extra edges, and run a correction refinement, whose result then
//   // needs to be communicated again as the new refinement step may introduce
//   // new edges that other PEs did not refine but are shared.
//   for (const auto& p : m_edgenodePe)        // for all PEs we share edges with
//     for (const auto& e : p.second) {        // for all refined edges on p.first
//       auto i = m_edgenode.find( e.first );  // find refined edge given parents
//       if (i != end(m_edgenode)) {           // found same added node on edge
//         // locally assigned added node ID and coordinates: i->second
//         // remotely assigned added node ID and coordinates: e.second
//         // Ensure global IDs of newly added nodes are the same
//         Assert( std::get< 0 >( i->second ) == std::get< 0 >( e.second ),
//                 "Remotely and locally assigned global ids mismatch" );
//         // Ensure coordinates are the same
//         Assert( std::abs( std::get<1>(i->second) - std::get<1>(e.second) ) <
//                   std::numeric_limits<tk::real>::epsilon() &&
//                 std::abs( std::get<2>(i->second) - std::get<2>(e.second) ) <
//                   std::numeric_limits<tk::real>::epsilon() &&
//                 std::abs( std::get<3>(i->second) - std::get<3>(e.second) ) <
//                   std::numeric_limits<tk::real>::epsilon(),
//                 "Remote and local added node coordinates mismatch" );
//       } else {  // remote PE added node on edge but we did not
//         // Make sure we know about this boundary-PE edge (that we did not refine)
//         Assert( m_bndEdges.find( CkMyPe() )->second.find( e.first ) !=
//                 m_bndEdges.find( CkMyPe() )->second.end(),
//                 "Local node IDs of boundary edge not found" );
//         // Save edge (given by parent global node IDs) to which the remote PE
//         // has added a new node but we did not. Will need to correct the mesh so
//         // it conforms across PEs.
//         extra.insert( {{ { tk::cref_find( m_lid, e.first[0] ),
//                            tk::cref_find( m_lid, e.first[1] ) } }} );
//       }
//     }
// 
//   // Store number of extra edges on this PE which this PE did not add but was
//   // refined by another PE, so now we need to tag and refine them and propagate
//   // reconnection of neighbor cells to arrive at conforming mesh across PE
//   // boundaries.
//   m_extra = extra.size();
// 
//   // Refine mesh triggered by nodes added on PE-boundary edges by other PEs
//   // PEs
//   correctRefine( extra );
// 
//   // Since refining edges that we originally did not but other PEs did may
//   // result in refining new edges that may be shared along PE boundaries, we
//   // now need to communicate these edges and potentially repeat the correction
//   // step. This must happen until all PEs that share edges can agree that there
//   // are no more edges to correct. Only then this refinement step can be
//   // considered complete.
//   if (m_pe.empty())
//     contribute( sizeof(std::size_t), &m_extra, CkReduction::max_int,
//                 m_cb.get< tag::matched >() );
//   else {
//     m_nref = 0;
//     for (auto p : m_pe) {       // for all PEs we share at least an edge with
//       // For all boundary edges of PE p, find out if we have added a new
//       // node to it, and if so, export parents->(newid,coords) to p.
//       tk::UnsMesh::EdgeNodeCoord exp;
//       for (const auto& e : tk::cref_find(m_bndEdges,p)) {
//         auto i = m_edgenode.find(e);
//         if (i != end(m_edgenode)) exp[ e ] = i->second;
//       }
//       thisProxy[ p ].addRefBndEdges( CkMyPe(), exp );
//     }
//   }
// }
// 
// void
// Partitioner::nextref()
// // *****************************************************************************
// // Decide wether to continue with another step of initial mesh refinement
// //! \details At this point the mesh has been refined and all PEs have received
// //!   a map associating the global IDs and the coordinates of a node added to
// //!   an edge during initial mesh refinement associated to all other PEs the
// //!   edges are shared with. Now the mesh is corrected so that it conforms
// //!   across PE-boundaries by tagging those edges for refinement that have been
// //!   refined by at least a PE. This concludes this initial mesh refinement
// //!   step, and we continue if there are more steps configured by the user.
// // *****************************************************************************
// {
//   // Remove initial mesh refinement step from list
//   if (!m_initref.empty()) m_initref.pop_back();
// 
//   if (!m_initref.empty())       // Continue to next initial refinement step
//     partref();
//   else {                        // Finish list of initial mesh refinement steps
//     // Output final mesh after initial mesh refinement
//     tk::UnsMesh refmesh( m_inpoel, m_coord );
//     tk::ExodusIIMeshWriter mw( "initref.final." + std::to_string(CkMyPe()),
//                                tk::ExoWriter::CREATE );
//     mw.writeMesh( refmesh );
//     // Finish initial mesh refinement
//     finishref();
//   }
// }
// 
// void
// Partitioner::finishref()
// // *****************************************************************************
// // Finish initial mesh refinement
// //! \details This function is called as after initial mesh refinement has
// //!   finished. If initial mesh reifnement was not configured by the user, this
// //!   is the point where we continue after the constructor, by computing the
// //!   total number of elements across the whole problem.
// // *****************************************************************************
// {
//   // Compute final number of cells across whole problem
//   auto nelem = m_ginpoel.size()/4;
//   contribute( sizeof(uint64_t), &nelem, CkReduction::sum_int,
//               m_cb.get< tag::refined >() );
// }
// 
// void
// Partitioner::uniformRefine()
// // *****************************************************************************
// // Do uniform mesh refinement
// // *****************************************************************************
// {
//   // Do uniform refinement
//   m_refiner.uniform_refinement();
// 
//   // Update mesh coordinates and connectivity
//   updateMesh();
// }
// 
// void
// Partitioner::errorRefine()
// // *****************************************************************************
// // Do error-based mesh refinement
// // *****************************************************************************
// {
//   // Find number of nodes in old mesh
//   auto npoin = tk::npoin( m_inpoel );
//   // Generate edges surrounding points in old mesh
//   auto esup = tk::genEsup( m_inpoel, 4 );
//   auto psup = tk::genPsup( m_inpoel, 4, esup );
// 
//   // Evaluate initial conditions at mesh nodes
//   auto u = nodeinit( npoin, esup );
// 
//   // Get the indices (in the system of systems) of refinement variables and the
//   // error indicator configured
//   const auto& refidx = g_inputdeck.get< tag::amr, tag::id >();
//   auto errtype = g_inputdeck.get< tag::amr, tag::error >();
// 
//   // Compute errors in ICs and define refinement criteria for edges
//   std::vector< edge_t > edge;
//   std::vector< real_t > crit;
//   AMR::Error error;
//   for (std::size_t p=0; p<npoin; ++p)   // for all mesh nodes on this PE
//     for (auto q : tk::Around(psup,p)) { // for all nodes surrounding p
//        tk::real cmax = 0.0;
//        edge_t e(p,q);
//        for (auto i : refidx) {          // for all refinement variables
//          auto c = error.scalar( u, e, i, m_coord, m_inpoel, esup, errtype );
//          if (c > cmax) cmax = c;        // find max error at edge
//        }
//        if (cmax > 0.0) {         // if nonzero error, will pass edge to refiner
//          edge.push_back( e );
//          crit.push_back( cmax );
//        }
//      }
// 
//   Assert( edge.size() == crit.size(), "Size mismatch" );
// 
//   // Do error-based refinement
//   m_refiner.error_refinement( edge, crit );
// 
//   // Update mesh coordinates and connectivity
//   updateMesh();
// }
// 
// void
// Partitioner::userRefine()
// // *****************************************************************************
// // Do mesh refinement based on user explicitly tagging edges
// // *****************************************************************************
// {
//   // Find number of nodes in old mesh
//   auto npoin = tk::npoin( m_inpoel );
//   // Generate edges surrounding points in old mesh
//   auto esup = tk::genEsup( m_inpoel, 4 );
//   auto psup = tk::genPsup( m_inpoel, 4, esup );
// 
//   // Get user-defined node-pairs (edges) to tag for refinement
//   const auto& edgenodelist = g_inputdeck.get< tag::amr, tag::edge >();
//   tk::UnsMesh::EdgeSet edgeset;
//   for (std::size_t i=0; i<edgenodelist.size()/2; ++i)
//     edgeset.insert( {{ {edgenodelist[i*2+0], edgenodelist[i*2+1]} }} );
// 
//   // Compute errors in ICs and define refinement criteria for edges
//   std::vector< edge_t > edge;
//   std::vector< real_t > crit;
//   for (std::size_t p=0; p<npoin; ++p)        // for all mesh nodes on this PE
//     for (auto q : tk::Around(psup,p)) {      // for all nodes surrounding p
//       tk::UnsMesh::Edge e{{p,q}};
//       if (edgeset.find(e) != end(edgeset)) { // tag edge if on user's list
//         edge.push_back( edge_t(e[0],e[1]) );
//         crit.push_back( 1.0 );
//       }
//     }
// 
//   Assert( edge.size() == crit.size(), "Size mismatch" );
// 
//   // Do error-based refinement
//   m_refiner.error_refinement( edge, crit );
// 
//   // Update mesh coordinates and connectivity
//   updateMesh();
// }
// 
// tk::Fields
// Partitioner::nodeinit( std::size_t npoin,
//                        const std::pair< std::vector< std::size_t >,
//                           std::vector< std::size_t > >& esup )
// // *****************************************************************************
// // Evaluate initial conditions (IC) at mesh nodes
// //! \param[in] npoin Number points in mesh (on this PE)
// //! \param[in] esup Elements surrounding points as linked lists, see tk::genEsup
// //! \return Initial conditions (evaluated at t0) at nodes
// // *****************************************************************************
// {
//   auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();
//   auto nprop = g_inputdeck.get< tag::component >().nprop();
// 
//   // Will store nodal ICs
//   tk::Fields u( m_coord[0].size(), nprop );
// 
//   // Evaluate ICs differently depending on nodal or cell-centered discretization
//   const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
//   if (scheme == ctr::SchemeType::MatCG || scheme == ctr::SchemeType::DiagCG) {
// 
//     // Node-centered: evaluate ICs for all scalar components integrated
//     for (const auto& eq : g_cgpde) eq.initialize( m_coord, u, t0 );
// 
//   } else if (scheme == ctr::SchemeType::DG) {
// 
//     // Initialize cell-centered unknowns
//     tk::Fields ue( m_inpoel.size()/4, nprop );
//     auto geoElem = tk::genGeoElemTet( m_inpoel, m_coord );
//     for (const auto& eq : g_dgpde) eq.initialize( geoElem, ue, t0 );
// 
//     // Transfer initial conditions from cells to nodes
//     for (std::size_t p=0; p<npoin; ++p) {       // for all mesh nodes on this PE
//       std::vector< tk::real > up( nprop, 0.0 );
//       tk::real vol = 0.0;
//       for (auto e : tk::Around(esup,p)) {       // for all cells around node p
//         // compute nodal volume: every element contributes their volume / 4
//         vol += geoElem(e,0,0) / 4.0;
//         // sum cell value to node weighed by cell volume / 4
//         for (std::size_t c=0; c<nprop; ++c)
//           up[c] += ue[e][c] * geoElem(e,0,0) / 4.0;
//       }
//       // store nodal value
//       for (std::size_t c=0; c<nprop; ++c) u(p,c,0) = up[c] / vol;
//     }
// 
//   } else Throw( "Nodal initialization not handled for discretization scheme" );
// 
//   Assert( u.nunk() == m_coord[0].size(), "Size mismatch" );
//   Assert( u.nprop() == nprop, "Size mismatch" );
// 
//   return u;
// }
// 
// void
// Partitioner::correctRefine( const tk::UnsMesh::EdgeSet& extra )
// // *****************************************************************************
// // Do mesh refinement correcting PE-boundary edges
// //! \param[in] extra Unique set of edges that need a new node on PE boundaries
// // *****************************************************************************
// {
//   if (!extra.empty()) {
//     // Generate list of edges that need to be corrected
//     std::vector< edge_t > edge;
//     for (const auto& e : extra) edge.push_back( edge_t(e[0],e[1]) );
//     std::vector< real_t > crit( edge.size(), 1.0 );
//   
//     // Do refinement including edges that need to be corrected
//     m_refiner.error_refinement( edge, crit );
//   
//     // Update mesh coordinates and connectivity
//     updateMesh();
//   }
// }
// 
// void
// Partitioner::updateMesh()
// // *****************************************************************************
// // Update mesh after refinement
// // *****************************************************************************
// {
//   // Get refined mesh connectivity
//   const auto& refinpoel = m_refiner.tet_store.get_active_inpoel();
//   Assert( refinpoel.size()%4 == 0, "Inconsistent refined mesh connectivity" );
// 
//   // Generate unique node lists of old and refined mesh using local ids
//   std::unordered_set< std::size_t > old( m_inpoel.cbegin(), m_inpoel.cend() );
//   std::unordered_set< std::size_t > ref( refinpoel.cbegin(), refinpoel.cend() );
// 
//   updateVolumeMesh( old, ref );
// 
//   updateBoundaryMesh( old, ref );
// 
//   // Update mesh connectivity with local node IDs
//   m_inpoel = refinpoel;
// 
//   // Update mesh connectivity with new global node ids
//   m_ginpoel = m_inpoel;
//   Assert( tk::uniquecopy(m_ginpoel).size() == m_coord[0].size(),
//           "Size mismatch" );
//   for (auto& i : m_ginpoel) i = m_gid[i];
// }
// 
// void
// Partitioner::updateVolumeMesh( const std::unordered_set< std::size_t >& old,
//                                const std::unordered_set< std::size_t >& ref )
// // *****************************************************************************
// //  Update volume mesh after mesh refinement
// //! \param[in] old Unique nodes of the old (unrefined) mesh using local ids
// //! \param[in] ref Unique nodes of the refined mesh using local ids
// // *****************************************************************************
// {
//   auto& x = m_coord[0];
//   auto& y = m_coord[1];
//   auto& z = m_coord[2];
// 
//   // Resize node coordinate storage to accommodate refined mesh nodes
//   auto npoin = ref.size();
//   x.resize( npoin );
//   y.resize( npoin );
//   z.resize( npoin );
//   m_gid.resize( npoin, std::numeric_limits< std::size_t >::max() );
// 
//   // Generate coordinates and ids to newly added nodes after refinement step
//   for (auto r : ref) {               // for all unique nodes of the refined mesh
//     if (old.find(r) == end(old)) {   // if node is newly added (in this step)
//       // get (local) parent ids of newly added node
//       auto p = m_refiner.node_connectivity.get( r );
//       Assert( old.find(p[0]) != end(old) && old.find(p[1]) != end(old),
//               "Parent(s) not in old mesh" );
//       Assert( r >= old.size(), "Attempting to overwrite node with added one" );
//       // generate coordinates for newly added node
//       x[r] = (x[p[0]] + x[p[1]])/2.0;
//       y[r] = (y[p[0]] + y[p[1]])/2.0;
//       z[r] = (z[p[0]] + z[p[1]])/2.0;
//       decltype(p) gp{{ m_gid[p[0]], m_gid[p[1]] }}; // global parent ids
//       // generate new global ID for newly added node
//       auto g = tk::UnsMesh::EdgeHash()( gp );
//       // ensure newly generated node id has not yet been used
//       Assert( g >= old.size(), "Hashed id overwriting old id" );
//       Assert( m_coordmap.find(g) == end(m_coordmap),
//               "Hash collision: ID already exist" );
//       // assign new global ids to local->global and to global->local maps
//       m_gid[r] = g;
//       Assert( m_lid.find(g) == end(m_lid),
//               "Overwriting entry global->local node ID map" );
//       m_lid[g] = r;
//       // assign new coordinates to new global node id
//       Assert( m_coordmap.find(g) == end(m_coordmap),
//               "Overwriting entry coordmap" );
//       m_coordmap.insert( {g, {{x[r], y[r], z[r]}}} );
//       // assign new coordinates and new global node id to global parent id pair
//       m_edgenode[ gp ] = std::make_tuple( g, x[r], y[r], z[r] );
//     }
//   }
// 
//   Assert( m_gid.size() == m_lid.size(), "Size mismatch" );
// 
//   Assert( std::none_of( begin(m_gid), end(m_gid), [](std::size_t i){
//             return i == std::numeric_limits< std::size_t >::max(); } ),
//           "Not all local->global node IDs have been assigned" );
// }
// 
// void
// Partitioner::updateBoundaryMesh( const std::unordered_set< std::size_t >& old,
//                                  const std::unordered_set< std::size_t >& ref )
// // *****************************************************************************
// // Update boundary data structures after mesh refinement
// //! \param[in] old Unique nodes of the old (unrefined) mesh using local ids
// //! \param[in] ref Unique nodes of the refined mesh using local ids
// // *****************************************************************************
// {
//   IGNORE(ref);  // to avoid compiler warning when asserts are optimized away
// 
//   using Face = tk::UnsMesh::Face;
// 
//   // Will associate a pair of side set id and adjacent tet id to a boundary
//   // triangle face
//   using BndFaces = std::unordered_map< Face,
//                                        std::pair< int, std::size_t >,
//                                        tk::UnsMesh::FaceHasher,
//                                        tk::UnsMesh::FaceEq >;
// 
//   // Generate data structure that associates the pair of side set id and
//   // adjacent tet id to a boundary triangle face for all boundary faces. After
//   // this loop we will have all tets adjacent to boundary faces, where the
//   // "boundary" includes all physical boundaries (where the user may assign
//   // boundary conditions via side sets from the mesh file) as well as
//   // PE-boundaries due to domain decomposition.
//   BndFaces bnd;
//   auto oldesuel = tk::genEsuelTet( m_inpoel, tk::genEsup(m_inpoel,4) );
//   for (std::size_t e=0; e<oldesuel.size()/4; ++e) {     // for tets on this PE
//     auto mark = e*4;
//     for (std::size_t f=0; f<4; ++f) {
//       if (oldesuel[mark+f] == -1) {  // if a face does not have an adjacent  tet
//         Face t{{ m_gid[ m_inpoel[ mark+tk::lpofa[f][0] ] ],
//                  m_gid[ m_inpoel[ mark+tk::lpofa[f][1] ] ],
//                  m_gid[ m_inpoel[ mark+tk::lpofa[f][2] ] ] }};
//         // associate tet id adjacent to boundary face to boundary face, at this
//         // point we fill in -1 for the side set id, since we don't know it yet
//         bnd[t] = {-1,e};
//       }
//     }
//   }
// 
//   // Assign side set ids to faces on the physical boundary only
//   for (const auto& ss : m_bface)  // for all phsyical boundaries (sidesets)
//     for (auto f : ss.second) {    // for all faces on this physical boundary
//       Face t{{ m_triinpoel[f*3+0], m_triinpoel[f*3+1], m_triinpoel[f*3+2] }};
//       auto cf = bnd.find( t );
//       if (cf != end(bnd)) cf->second.first = ss.first;
//     }
// 
//   // Remove PE-boundary faces (to which the above loop did not assign a set id)
//   auto bndcopy = bnd;
//   for (const auto& f : bndcopy) if (f.second.first == -1) bnd.erase( f.first );
//   tk::destroy( bndcopy );
// 
//   // Now in bnd we have a pair of side set ids of all tet ids that are adjacent
//   // to all physical boundary faces that are associated to a side set in the
//   // mesh file and only along faces on the physical boundary and not along faces
//   // on the PE-boundary.
// 
//   // storage for boundary faces associated to side-set IDs of the refined mesh
//   decltype(m_bface) bface;              // will become m_bface
//   // storage for boundary faces-node connectivity of the refined mesh
//   decltype(m_triinpoel) triinpoel;      // will become m_triinpoel
//   // face id counter
//   std::size_t facecnt = 0;
// 
//   // Lambda to associate a boundary face and connectivity to a side set.
//   // Parameter 's' is the list of faces (ids) to add new face to. Parameter 'f'
//   // is the triangle face connecttivity to add.
//   auto addBndFace = [ &facecnt, &triinpoel ]
//                     ( std::vector< std::size_t >& s, const Face& f )
//   {
//     s.push_back( facecnt++ );
//     triinpoel.push_back( f[0] );
//     triinpoel.push_back( f[1] );
//     triinpoel.push_back( f[2] );
//   };
// 
//   // Lambda to find the nodes of the parent face of a child face. Parameter
//   // 'face' denotes the node ids of the child face whose parent face we are
//   // looking for. This search may find 3 or 4 parent nodes, depending on whether
//   // the child shares or does not share a face with the parent tet,
//   // respectively.
//   auto parentFace = [ &old, this ]( const Face& face ){
//     std::unordered_set< std::size_t > s;// will store nodes of parent face
//     for (auto n : face) {               // for all 3 nodes of the face
//       if (old.find(n) != end(old))      // if child node found in old mesh,
//         s.insert( n );                  // that node is also in the parent face
//       else {                            // if child node is a newly added one
//         // find its parent nodes and store both uniquely
//         auto p = this->m_refiner.node_connectivity.get( n );
//         s.insert( begin(p), end(p) );
//       }
//     }
//     // Ensure all parent nodes are part of the old (unrefined) mesh
//     Assert( std::all_of( begin(s), end(s), [ &old ]( std::size_t j ){
//                            return old.find(j) != end(old); } ),
//             "Old mesh nodes not in old mesh" );
//     // Return unique set of nodes of the parent face of child face
//     return s;
//   };
// 
//   // Generate boundary face data structures after refinement step
//   for (const auto& f : bnd) {
//     // construct face of old mesh boundary face using local ids
//     Face oldface{{ tk::cref_find( m_lid, f.first[0] ),
//                    tk::cref_find( m_lid, f.first[1] ),
//                    tk::cref_find( m_lid, f.first[2] ) }};
//     // will associate to side set id of old (unrefined) mesh boundary face
//     auto& side = bface[ f.second.first ];
//     // query number of children of boundary tet adjacent to boundary face
//     auto nc = m_refiner.tet_store.data( f.second.second ).num_children;
//     if (nc == 0) {      // if boundary tet is not refined, add its boundary face
//       addBndFace( side, f.first );
//     } else {            // if boundary tet is refined
//       const auto& tets = m_refiner.tet_store.tets;
//       for (decltype(nc) i=0; i<nc; ++i ) {      // for all child tets
//         // get child tet id
//         auto childtet = m_refiner.tet_store.get_child_id( f.second.second, i );
//         auto ct = tets.find( childtet );
//         Assert( ct != end(tets), "Child tet not found" );
//         // ensure all nodes of child tet are in refined mesh
//         Assert( ref.find(ct->second[0]) != end(ref) &&
//                 ref.find(ct->second[1]) != end(ref) &&
//                 ref.find(ct->second[2]) != end(ref) &&
//                 ref.find(ct->second[3]) != end(ref),
//                 "Boundary child tet node id not found in refined mesh" );
//         // get nodes of child tet
//         auto A = ct->second[0];
//         auto B = ct->second[1];
//         auto C = ct->second[2];
//         auto D = ct->second[3];
//         // form all 4 faces of child tet
//         std::array<Face, 4> face{{{{A,C,B}}, {{A,B,D}}, {{A,D,C}}, {{B,C,D}}}};
//         for (const auto& rf : face) {   // for all faces of child tet
//           // find nodes of the parent face of child face
//           auto parfac = parentFace( rf );
//           // if the child shares a face with its parent and all 3 nodes of the
//           // parent face of the child's face are the same, the child face is on
//           // the same boundary as the parent face
//           if ( parfac.size() == 3 &&
//                parfac.find(oldface[0]) != end(parfac) &&
//                parfac.find(oldface[1]) != end(parfac) &&
//                parfac.find(oldface[2]) != end(parfac) )
//           {
//             addBndFace( side, {{m_gid[rf[0]], m_gid[rf[1]], m_gid[rf[2]]}} );
//           }
//         }
//       }
//     }
//   }
// 
//   // Update boundary face data structures
//   m_bface = std::move(bface);
//   m_triinpoel = std::move(triinpoel);
// }

#include "NoWarning/refiner.def.h"
