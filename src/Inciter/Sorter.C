// *****************************************************************************
/*!
  \file      src/Inciter/Sorter.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Mesh sorter for global distributed mesh reordering
  \see       Sorter.h for more info.
*/
// *****************************************************************************

#include <vector>
#include <algorithm>

#include "Sorter.h"
#include "Reorder.h"
#include "DerivedData.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::Sorter;

Sorter::Sorter( const CProxy_Transporter& transporter,
                const tk::CProxy_Solver& solver,
                const tk::SorterCallback& cbs,
                const Scheme& scheme,
                const std::vector< std::size_t >& ginpoel,
                const tk::UnsMesh::CoordMap& coordmap,
                const std::map< int, std::vector< std::size_t > >& bface,
                const std::vector< std::size_t >& triinpoel,
                const std::map< int, std::vector< std::size_t > >& bnode,
                int nchare ) :
  m_host( transporter ),
  m_solver( solver ),
  m_cbs( cbs ),
  m_scheme( scheme ),
  m_ginpoel( ginpoel ),
  m_coordmap( coordmap ),
  m_bface( bface ),
  m_triinpoel( triinpoel ),
  m_bnode( bnode ),
  m_nchare( nchare ),
  m_nodeset( begin(ginpoel), end(ginpoel) ),
  m_nquery( 0 ),
  m_nmask( 0 ),
//   m_lower( 0 ),
//   m_upper( 0 ),
  m_ncomm(),
  m_ncommunication()
// *****************************************************************************
//  Constructor: prepare owned mesh node IDs for reordering
//! \param[in] bface Face lists mapped to side set ids
//! \param[in] triinpoel Interconnectivity of points and boundary-face
// *****************************************************************************
{
  // Find chare-boundary nodes
  std::vector< std::size_t > chbnode;
  auto el = tk::global2local( ginpoel );      // generate local mesh data
  const auto& inpoel = std::get< 0 >( el );   // local connectivity
  const auto& gid = std::get< 1 >( el );      // local->global node ids
  auto esup = tk::genEsup( inpoel, 4 );       // elements surrounding points
  auto esuel = tk::genEsuelTet( inpoel, esup ); // elems surrounding elements
  for (std::size_t e=0; e<esuel.size()/4; ++e) {
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f) {
      if (esuel[mark+f] == -1) {
        chbnode.push_back( gid[ inpoel[ mark+tk::lpofa[f][0] ] ] );
        chbnode.push_back( gid[ inpoel[ mark+tk::lpofa[f][1] ] ] );
        chbnode.push_back( gid[ inpoel[ mark+tk::lpofa[f][2] ] ] );
      }
    }
  }
  // Make boundary nodes unique
  tk::unique( chbnode );

  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.chbnd();

  // Send chare boundary nodes to all fellow chares in a broadcast fashion
  thisProxy.query( thisIndex, chbnode );
}

void
Sorter::query( int fromch, const std::vector< std::size_t >& bnodes )
// *****************************************************************************
//  Query mesh nodes to identify if they are shared
//! \param[in] fromch Querying chare
//! \param[in] nodes List of global mesh node IDs to query
//! \details Note that every chare calls this function in a broadcast fashion,
//!   including our own. However, to compute the correct result, this would
//!   only be necessary for chare whose ID is higher than ours. However, the
//!   broadcast (calling everyone) is more efficient. This also results in a
//!   simpler logic, because every chare goes through this single call path.
//!   The returned mask is a vector of node IDs found on this chare.
// *****************************************************************************
{
  std::vector< std::size_t > found;
  for (auto j : bnodes) {
    const auto it = m_nodeset.find( j );
    if (it != end(m_nodeset)) found.push_back( j );
  }

  if ( g_inputdeck.get< tag::cmd, tag::feedback >() && ++m_nquery == m_nchare )
    m_host.chquery();

  // Return nodes found to sender
  thisProxy[ fromch ].mask( thisIndex, found );
}

void
Sorter::mask( int fromch, const std::vector< std::size_t >& found )
// *****************************************************************************
//  Receive mesh node IDs found on sender chare
//! \param[in] c The chare that has found the mesh nodes passed in
//! \param[in] found Mesh nodes shared with chare c
//! \details Note that every chare will call this function, since query() was
//!   called in a broadcast fashion and query() answers to every chare once.
//!   This is more efficient than calling only the chares from which we would
//!   have to receive results from. Thus the incoming results are only
//!   interesting from chares with lower IDs than ours.
// *****************************************************************************
{
  // Store the global mesh node IDs associated to chare IDs bordering our mesh
  // chunk. This loop computes m_msum, a symmetric chare-node communication map,
  // that associates a unique set of global node IDs to chare IDs we share these
  // nodes with.
  if (fromch != thisIndex) m_msum[ fromch ].insert( begin(found), end(found) );

  if (++m_nmask == m_nchare) {
    create();   // shortcut
   }

//   // Associate global mesh node IDs to lower PEs we will need to receive from
//   // during node reordering. The choice of associated container is std::map,
//   // which is ordered (vs. unordered, hash-map). This is required by the
//   // following operation that makes the mesh node IDs unique in the
//   // communication map. (We are called in an unordered fashion, so we need to
//   // collect from all PEs and then we need to make the node IDs unique, keeping
//   // only the lowest PEs a node ID is associated with.) Note that m_ncomm is an
//   // asymmetric PE-node communication map, associating a set of unique global
//   // node IDs to PE IDs from which we (this PE) will need to receive newly
//   // assigned node IDs during global mesh node reordering. This map is
//   // asymmetric, beacuse of the agreement of the reordering that if a mesh node
//   // is shared by multiple PEs, the PE with the lowest ID gets to assign a new
//   // ID to it and all others must receive it instead of assigning it.
// 
//   if (p < CkMyPe()) {
//     auto& id = m_ncomm[ p ];
//     for (const auto& h : cn) id.insert( h.first );
//   }
// 
//   if (++m_nmask == static_cast<std::size_t>(CkNumPes())) {
//     // Make sure we have received all we need
//     Assert( m_ncomm.size() == static_cast<std::size_t>(CkMyPe()),
//             "Communication map size on PE " +
//             std::to_string(CkMyPe()) + " must equal " +
//             std::to_string(CkMyPe()) );
//     // Fill new hash-map, keeping only unique node IDs obtained from the
//     // lowest possible PEs
//     for (auto c=m_ncomm.cbegin(); c!=m_ncomm.cend(); ++c) {
//       auto& n = m_ncommunication[ c->first ];
//       for (auto j : c->second)
//         if (std::none_of( m_ncomm.cbegin(), c,
//              [ j ]( const typename decltype(m_ncomm)::value_type& s )
//              { return s.second.find(j) != end(s.second); } )) {
//           n.insert(j);
//         }
//       if (n.empty()) m_ncommunication.erase( c->first );
//     }
//     // Free storage of temporary communication map used to receive global
//     // mesh node IDs as it is no longer needed once the final communication
//     // map is generated.
//     tk::destroy( m_ncomm );
//     // Count up total number of nodes we will need receive during reordering
//     std::size_t nrecv = 0;
//     for (const auto& u : m_ncommunication) nrecv += u.second.size();
// 
//     if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pemask();
// 
//     // Compute number of mesh node IDs we will assign IDs to
//     auto nuniq = m_nodeset.size() - nrecv;
// 
//     // Start computing PE offsets for node reordering
//     thisProxy.offset( CkMyPe(), nuniq );
//   }
}


// void
// Sorter::offset( int p, std::size_t u )
// // *****************************************************************************
// //  Receive number of uniquely assigned global mesh node IDs from lower PEs
// //! \param[in] p PE ID
// //! \param[in] u Number of mesh node IDs PE p will assign IDs to
// //! \details This function computes the offset each PE will need to start
// //!   assigning its new node IDs from (for those nodes that are not assigned
// //!   new IDs by any PEs with lower indices). The offset for a PE is the
// //!   offset for the previous PE plus the number of node IDs the previous PE
// //!   (uniquely) assigns new IDs for minus the number of node IDs the
// //!   previous PE receives from others (lower PEs). This is computed here in
// //!   a parallel/distributed fashion by each PE sending its number of node
// //!   IDs (that it uniquely assigns) to all PEs. Note that each PE would
// //!   only need to send this information to higher PEs, but instead this
// //!   function is called in a broadcast fashion, because that is more
// //!   efficient than individual calls to only the higher PEs. Therefore when
// //!   computing the offsets, we only count the lower PEs. When this is done,
// //!   we have the precise communication map as well as the start offset on
// //!   all PEs and so we can start the distributed global mesh node ID
// //!   reordering.
// // *****************************************************************************
// {
//   if (p < CkMyPe()) m_start += u;
//   if (++m_noffset == static_cast<std::size_t>(CkNumPes())) reorder();
// }
// 
// 
// void
// Sorter::request( int p, const std::unordered_set< std::size_t >& nd )
// // *****************************************************************************
// //  Request new global node IDs for old node IDs
// //! \param[in] p PE request coming from and to which we send new IDs to
// //! \param[in] nd Set of old node IDs whose new IDs are requested
// // *****************************************************************************
// {
//   // Queue up requesting PE and node IDs
//   m_reqNodes.push_back( { p, nd } );
//   // Trigger SDAG wait signaling that node IDs have been requested from us
//   nodes_requested_complete();
// }
// 
// void
// Sorter::neworder( const std::unordered_map< std::size_t,
//                         std::tuple< std::size_t, tk::UnsMesh::Coord > >& nodes )
// // *****************************************************************************
// //  Receive new (reordered) global node IDs
// //! \param[in] nd Map associating new to old node IDs
// // *****************************************************************************
// {
//   // Signal to the runtime system that we have participated in reordering
//   participated_complete();
// 
//   // Store new node IDs associated to old ones, and node coordinates associated
//   // to new node IDs. Since multiple chares can contribute to a single node, we
//   // store such shared node coordinates for all chares that contribute.
//   for (const auto& p : nodes) {
//     auto id = std::get< 0 >( p.second );
//     auto coord = std::get< 1 >( p.second );
//     m_linnodes[ p.first ] = id;
//     for (auto c : tk::cref_find(m_nodech,p.first))
//       m_chcoordmap[ c ].emplace( id, coord );
//   }
// 
//   // If all our nodes have new IDs assigned, reorder complete on this PE
//   if (m_linnodes.size() == m_nodeset.size()) reordered();
// }
// 
// void
// Sorter::lower( std::size_t low )
// // *****************************************************************************
// //  Receive lower bound of node IDs our PE operates on after reordering
// //! \param[in] low Lower bound of node IDs assigned to us
// // *****************************************************************************
// {
//   m_lower = low;
//   lower_complete();
// }
// 
// void
// Sorter::stdCost( tk::real av )
// // *****************************************************************************
// //  Compute the variance of the communication cost
// //! \param[in] av Average of the communication cost
// //! \details Computing the standard deviation is done via computing and
// //!   summing up the variances on each PE and asynchronously reducing the
// //!   sum to our host.
// // *****************************************************************************
// {
//   tk::real var = (m_cost-av)*(m_cost-av);
//   contribute( sizeof(tk::real), &var, CkReduction::sum_double,
//               m_cbp.get< tag::stdcost >() );
// }
// 
// tk::real
// Sorter::cost( std::size_t l, std::size_t u )
// // *****************************************************************************
// //  Compute communication cost on our PE
// //! \param[in] l Lower global node ID this PE works on
// //! \param[in] u Upper global node ID this PE works on
// //! \return Communication cost for our PE
// //! \details The cost is a real number between 0 and 1, defined as the
// //!   number of mesh points we do not own, i.e., need to send to some other
// //!   PE, divided by the total number of points we contribute to. The lower
// //!   the better.
// // *****************************************************************************
// {
//   std::size_t ownpts = 0, compts = 0;
//   for (auto p : m_nodeset) if (p >= l && p < u) ++ownpts; else ++compts;
// 
//   // Free storage of unique global node IDs chares on our PE will contribute to
//   // as it is no longer needed after computing the communication cost.
//   tk::destroy( m_nodeset );
// 
//   return static_cast<tk::real>(compts) / static_cast<tk::real>(ownpts + compts);
// }
// 
// void
// Sorter::reorder()
// // *****************************************************************************
// //  Reorder global mesh node IDs
// // *****************************************************************************
// {
//   // Activate SDAG waits for having requests arrive from other PEs for some
//   // of our node IDs; and for computing/receiving lower and upper bounds of
//   // global node IDs our PE operates on after reordering
//   thisProxy[ CkMyPe() ].wait4prep();
//   thisProxy[ CkMyPe() ].wait4bounds();
// 
//   // In serial signal to the runtime system that we have participated in
//   // reordering. This is required here because this is only triggered if
//   // communication is required during mesh node reordering. See also
//   // particioner.ci.
//   if (CkNumPes() == 1) participated_complete();
// 
//   // Send out request for new global node IDs for nodes we do not reorder
//   for (const auto& c : m_ncommunication)
//     thisProxy[ c.first ].request( CkMyPe(), c.second );
// 
//   // Lambda to decide if node ID is being assigned a new ID by us
//   auto ownnode = [ this ]( std::size_t p ) {
//     using Set = typename std::remove_reference<
//                   decltype(m_ncommunication) >::type::value_type;
//     return !std::any_of( m_ncommunication.cbegin(), m_ncommunication.cend(),
//                          [&](const Set& s)
//                          { return s.second.find(p) != s.second.cend(); } );
//   };
// 
//   // Reorder our chunk of the mesh node IDs by looping through all of our node
//   // IDs. We test if this PE is to assign a new ID to a node ID, and if so, we
//   // assign a new ID, i.e., reorder, by constructing a map associating new to
//   // old IDs. We also count up the reordered nodes, which also serves as the new
//   // node id. Also, we store the coordinates associated to the new node ID for
//   // each chare on this PE. Since multiple chares can contribute to a single
//   // node, we store such shared node coordinates for all chares that contribute.
//   for (auto p : m_nodeset)
//     if (ownnode(p)) {
//       m_linnodes[ p ] = m_start;
//       auto coord = tk::cref_find( m_coordmap, p );
//       for (auto c : tk::cref_find(m_nodech,p))
//         m_chcoordmap[ c ].emplace( m_start, coord );
//       ++m_start;
//     }
// 
//   // Trigger SDAG wait indicating that reordering own node IDs are complete
//   reorderowned_complete();
// 
//   // If all our nodes have new IDs assigned, reordering complete on this PE
//   if (m_linnodes.size() == m_nodeset.size()) reordered();
// }
// 
// void
// Sorter::prepare()
// // *****************************************************************************
// //  Associate new node IDs to old ones and return them to the requestor(s)
// // *****************************************************************************
// {
//   // Signal to the runtime system that we have participated in reordering
//   participated_complete();
// 
//   // Find and return new node IDs to sender
//   for (const auto& r : m_reqNodes) {
//     std::unordered_map< std::size_t,
//       std::tuple< std::size_t, tk::UnsMesh::Coord > > n;
//     for (auto p : r.second)
//       n.emplace( p, std::make_tuple( tk::cref_find(m_linnodes,p),
//                                      tk::cref_find(m_coordmap,p) ) );
//     thisProxy[ r.first ].neworder( n );
//     tk::destroy( n );
//   }
// 
//   tk::destroy( m_reqNodes ); // Clear queue of requests just fulfilled
// 
//   // Re-enable SDAG wait for preparing new node requests
//   thisProxy[ CkMyPe() ].wait4prep();
// 
//   // Re-enable trigger signaling that reordering of owned node IDs are
//   // complete right away
//   reorderowned_complete();
// }
// 
// void
// Sorter::reordered()
// // *****************************************************************************
// //  Compute final result of reordering
// //! \details At this point the node coordinates on all PEs have been updated to
// //!   be consistent with the new ordering. We continue by updating other data,
// //!   such as mesh connectivities of chares.
// // *****************************************************************************
// {
//   tk::destroy( m_bnodechares );
//   tk::destroy( m_ncommunication );
// 
//   // Construct maps associating old node IDs (as in file) to new node IDs
//   // (as in producing contiguous-row-id linear system contributions)
//   // associated to chare IDs (outer key).
//   for (const auto& c : m_chinpoel) {
//     auto& nodes = m_chfilenodes[ c.first ];
//     for (auto p : c.second)
//       nodes[ tk::cref_find(m_linnodes,p) ] = p;
//   }
// 
//   // Update chare-categorized elem connectivities with the reordered node IDs
//   for (auto& c : m_chinpoel)
//     for (auto& p : c.second)
//        p = tk::cref_find( m_linnodes, p );
// 
//   // Update chare-categorized chare-mesh-nodes comm map with the reordered IDs
//   for (auto& c : m_msum)
//     for (auto& s : c.second) {
//       decltype(s.second) n;
//       for (auto p : s.second) {
//         n.insert( tk::cref_find( m_linnodes, p ) );
//       }
//       s.second = std::move( n );
//     }
// 
//   // Update unique global node IDs chares on our PE will contribute to with
//   // the reordered node IDs
//   m_nodeset.clear();
//   for (const auto& c : m_chinpoel)
//     for (auto i : c.second)
//       m_nodeset.insert( i );
// 
//   if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pereordered();
// 
//   // Compute lower and upper bounds of reordered node IDs our PE operates on
//   bounds();
// }
// 
// void
// Sorter::bounds()
// // *****************************************************************************
// // Compute lower and upper bounds of reordered node IDs our PE operates on
// // \details This function computes the global row IDs at which the linear
// //   system will have a PE boundary. We simply find the largest node ID
// //   assigned on each PE by the reordering and use that as the upper global
// //   row index. Note that while this rarely results in equal number of rows
// //   assigned to PEs, potentially resulting in some load-imbalance, it
// //   yields a pretty good division reducing communication costs during the
// //   assembly of the linear system, which is more important than a slight
// //   (FLOP) load imbalance. Since the upper index for PE 1 is the same as
// //   the lower index for PE 2, etc., we find the upper indices and then the
// //   lower indices for all PEs are communicated.
// // *****************************************************************************
// {
//   m_upper = 0;
// 
//   using P1 = std::pair< const std::size_t, std::size_t >;
//   for (const auto& c : m_chfilenodes) {
//     auto x = std::max_element( begin(c.second), end(c.second),
//              [](const P1& a, const P1& b){ return a.first < b.first; } );
//     if (x->first > m_upper) m_upper = x->first;
//   }
// 
//   // The bounds are the dividers (global mesh point indices) at which the
//   // linear system assembly is divided among PEs. However, Hypre and thus
//   // Solver expect exclusive upper indices, so we increase the last one by
//   // one here. Note that the cost calculation, Sorter::cost() also
//   // expects exclusive upper indices.
//   if (CkMyPe() == CkNumPes()-1) ++m_upper;
// 
//   // Tell the runtime system that the upper bound has been computed
//   upper_complete();
// 
//   // Set lower index for PE 0 as 0
//   if (CkMyPe() == 0) lower(0);
// 
//   // All PEs except the last one send their upper indices as the lower index for
//   // PE+1
//   if (CkMyPe() < CkNumPes()-1) thisProxy[ CkMyPe()+1 ].lower( m_upper );
// }

void
Sorter::create()
// *****************************************************************************
// Create chare array elements on this PE and assign the global mesh element IDs
// they will operate on
// *****************************************************************************
{
//   if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.pebounds();
// 
//   // Initiate asynchronous reduction across all Sorter objects computing
//   // the average communication cost of merging the linear system
//   m_cost = cost( m_lower, m_upper );
//   contribute( sizeof(tk::real), &m_cost, CkReduction::sum_double,
//               m_cbp.get< tag::avecost >() );

  // Create worker chare array elements
  createDiscWorkers();

//   // Broadcast our bounds of global node IDs to all matrix solvers
//   const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
//   if (scheme == ctr::SchemeType::MatCG || scheme == ctr::SchemeType::DiagCG)
//     m_solver.bounds( CkMyPe(), m_lower, m_upper );
//   else // if no MatCG, no matrix solver, continue
//    contribute( m_cbp.get< tag::coord >() );
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
//   auto dist = distribution( m_nchare );
// 
//   for (int c=0; c<dist[1]; ++c) {
//     // Compute chare ID
//     auto cid = CkMyPe() * dist[0] + c;
//     // Guard those searches that operate on empty containers in serial
//     typename decltype(m_msum)::mapped_type msum;
//     if (!m_msum.empty()) msum = tk::cref_find( m_msum, cid );
//     // Create worker array element using Charm++ dynamic chare array element
//     // insertion: 1st arg: chare id, last arg: PE chare is created on, middle
//     // args: Discretization ctor args. See also Charm++ manual, Sec. "Dynamic
//     // Insertion".
//     m_scheme.discInsert( cid, m_host, tk::cref_find(m_chinpoel,cid),
//       tk::cref_find(m_chcoordmap,cid), msum, m_nchare, CkMyPe() );
//   }

  // Create worker array element using Charm++ dynamic chare array element
  // insertion: 1st arg: chare id, last arg: PE chare is created on, middle
  // args: Discretization ctor args. See also Charm++ manual, Sec. "Dynamic
  // Insertion".
  m_scheme.discInsert( thisIndex, m_host, m_ginpoel, m_coordmap, m_msum,
                       m_nchare, CkMyPe() );

  contribute( m_cbs.get< tag::discinserted >() );

//   // Free storage of global mesh node IDs associated to chare IDs bordering
//   // the mesh chunk held by and associated to chare IDs we own as it is no
//   // longer needed after creating the workers.
//   tk::destroy( m_msum );
}

void
Sorter::createWorkers()
// *****************************************************************************
//  Create worker chare array element
// *****************************************************************************
{
//   // Generate map associating new(value) to file(key) node ids for this chare
//   decltype(m_linnodes) newnodes;
//   const auto& chfilenodes = tk::cref_find( m_chfilenodes, cid );
//   for (auto i : chinpoel) newnodes[ tk::cref_find(chfilenodes,i) ] = i;

//   // Generate set of alll mesh faces
//   tk::UnsMesh::FaceSet faceset;
//   for (std::size_t e=0; e<m_ginpoel.size()/4; ++e) { // for all tets
//     auto mark = e*4;
//     for (std::size_t f=0; f<4; ++f) // for all tet faces
//       faceset.insert( {{{ m_ginpoel[ mark + tk::lpofa[f][0] ],
//                           m_ginpoel[ mark + tk::lpofa[f][1] ],
//                           m_ginpoel[ mark + tk::lpofa[f][2] ] }}} );
//   }

  // Generate boundary face and node ids (after mesh node reordering)
//   std::vector< std::size_t > chtriinpoel;
//   std::unordered_map< int, std::vector< std::size_t > > chbface;
//   std::size_t cnt = 0;

  // Extract our portion of the boundary node lists
  decltype(m_bnode) bnode;
  for (const auto& s : m_bnode) {
    auto& n = bnode[ s.first ];
    for (auto p : s.second) {
      if (m_nodeset.find(p) != end(m_nodeset))
        n.push_back( p );
    }
  }

//   for (const auto& ss : m_bface)  // for all phsyical boundaries (sidesets)
//     for (auto f : ss.second) {    // for all faces on this physical boundary
//       // attempt to find face nodes on this chare
//       tk::UnsMesh::Face t{{ m_triinpoel[f*3+0],
//                             m_triinpoel[f*3+1],
//                             m_triinpoel[f*3+2] }};
//       auto f1 = m_nodeset.find( t[0] );
//       auto f2 = m_nodeset.find( t[1] );
//       auto f3 = m_nodeset.find( t[2] );
//       // if face node on this chare, assign to side set
//       auto& n = bnode[ ss.first ];
//       if (f1 != end(m_nodeset)) n.push_back( t[0] );
//       if (f2 != end(m_nodeset)) n.push_back( t[1] );
//       if (f3 != end(m_nodeset)) n.push_back( t[2] );
// 
// //       // if all 3 nodes of the physical boundary face are on this chare
// //       if (f1 != end(newnodes) && f2 != end(newnodes) && f3 != end(newnodes)) {
// //         // Create face with new node ids (after mesh node reordering)
// //         std::array< std::size_t, 3 > t{{f1->second, f2->second, f3->second}};
// //         // if this boundary face is on this chare
// //         if (faceset.find(t) != end(faceset)) {
// //           // store face connectivity with new (global) node ids of this chare
// //           chtriinpoel.insert( end(chtriinpoel), begin(t), end(t) );
// //           // generate/store physical boundary face id associated to sideset id
// //           chbface[ ss.first ].push_back( cnt++ );
// //         }
// //       }
//     }

//   // Make boundary node IDs unique for each physical boundary (side set)
//   for (auto& s : bnode) tk::unique( s.second );
//
  // Face data class
  FaceData fd( m_ginpoel, m_bface, bnode, m_triinpoel );

  // Make sure (bound) base is already created and accessible
  Assert( m_scheme.get()[thisIndex].ckLocal() != nullptr,
          "About to pass nullptr" );

  // Create worker array element using Charm++ dynamic chare array element
  // insertion: 1st arg: chare id, last arg: PE chare is created on, middle
  // args: Discretization's child ctor args. See also Charm++ manual, Sec.
  // "Dynamic Insertion".
  m_scheme.insert( thisIndex, m_scheme.get(), m_solver, fd, CkMyPe() );

  contribute( m_cbs.get< tag::workinserted >() );
}

#include "NoWarning/sorter.def.h"
