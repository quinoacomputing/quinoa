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
#include "Solver.h"
#include "HashMapReducer.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

static CkReduction::reducerType BndNodeMerger;

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
  m_noffset( 0 ),
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
//! \param[in] transporter Transporter (host) Charm++ proxy
//! \param[in] solver Linear system solver Charm++ proxy
//! \param[in] cbs Charm++ callbacks for Sorter
//! \param[in] scheme Discretization scheme
//! \param[in] ginpoel Mesh connectivity (this chare) using global node IDs
//! \param[in] coordmap Mesh node coordinates (this chare) for global node IDs
//! \param[in] bface Face lists mapped to side set ids
//! \param[in] triinpoel Interconnectivity of points and boundary-faces
//! \param[in] bnode Node ids mapped to side set ids
//! \param[in] nchare Total number of Charm++ Refiner chares
// *****************************************************************************
{
  // Ensure boundary face ids will not index out of face connectivity
  Assert( std::all_of( begin(m_bface), end(m_bface),
            [&](const decltype(m_bface)::value_type& s)
            { return std::all_of( begin(s.second), end(s.second),
                                  [&](decltype(s.second)::value_type f)
                                  { return f*3+2 < m_triinpoel.size(); } ); } ),
          "Boundary face data structures inconsistent" );

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

  if (g_inputdeck.get< tag::cmd, tag::feedback >()) m_host.chbnd();

  // Aggregate boundary nodes across all Sorter chares
  std::unordered_map< int, std::vector< std::size_t > >
    bnd{{ thisIndex, std::move(chbnode) }};
  auto stream = tk::serialize( bnd );
  contribute( stream.first, stream.second.get(), BndNodeMerger,
    CkCallback(CkIndex_Sorter::comm(nullptr),thisProxy) );
}

void
Sorter::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types
//! \details Since this is a [nodeinit] routine, the runtime system executes the
//!   routine exactly once on every logical node early on in the Charm++ init
//!   sequence. Must be static as it is called without an object. See also:
//!   Section "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  BndNodeMerger = CkReduction::addReducer(
                    tk::mergeHashMap< int, std::vector< std::size_t > > );
}

void
Sorter::comm( CkReductionMsg* msg )
// *****************************************************************************
//  Receive aggregated chare boundary nodes associated to chares
//! \param[in] msg Charm++ reduction message containing aggregate boundary nodes
//! \details This is a reduction target receiving the aggregated chare boundary
//!    nodes associated to chare IDs across the whole problem. Here we setup
//!    the chare-node communication map.
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::feedback >()) m_host.chcomm();

   // Unpack final result of aggregating boundary nodes associated to chares
   std::unordered_map< int, std::vector< std::size_t > > bnd;
   PUP::fromMem creator( msg->getData() );
   creator | bnd;
   delete msg;

  // Store the global mesh node IDs associated to chare IDs bordering our mesh
  // chunk. This loop computes m_msum, a symmetric chare-node communication map,
  // that associates a unique set of global node IDs to chare IDs we share these
  // nodes with.
   for (const auto& c : bnd)
     if (c.first != thisIndex)
       for (auto i : c.second) {
         const auto it = m_nodeset.find( i );
         if (it != end(m_nodeset)) m_msum[ c.first ].insert( i );
       }

  if (g_inputdeck.get< tag::discr, tag::reorder >())
    mask();   // continue with mesh node reordering if requested (or required)
  else
    create(); // skip mesh node reordering
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
      for (auto j : c->second)
        if (std::none_of( m_msum.cbegin(), c,
             [ j ]( const decltype(m_msum)::value_type& s )
             { return s.second.find(j) != end(s.second); } )) {
          n.insert(j);
        }
      if (n.empty()) m_reordcomm.erase( c->first );
    }

  // Count up total number of nodes this chare will need to receive
  std::size_t nrecv = 0;
  for (const auto& u : m_reordcomm) nrecv += u.second.size();

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
  thisProxy[ thisIndex ].wait4bounds();

  // Send out request for new global node IDs for nodes we do not reorder
  for (const auto& c : m_reordcomm)
    thisProxy[ c.first ].request( thisIndex, c.second );

  // Lambda to decide if node is assigned a new ID by this chare. If node is not
  // found in the asymmetric communication map, it is owned, i.e., this chare
  // assigns its new id.
  auto ownnode = [ this ]( std::size_t p ) {
    return std::all_of( m_reordcomm.cbegin(), m_reordcomm.cend(),
                        [&](const decltype(m_reordcomm)::value_type& s)
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
  for (const auto& r : m_reqnodes) {
    std::unordered_map< std::size_t,
      std::tuple< std::size_t, tk::UnsMesh::Coord > > n;
    for (auto p : r.second) {
      auto newid = tk::cref_find( m_newnodes, p );
      n.emplace( p,
        std::make_tuple( newid, tk::cref_find(m_newcoordmap,newid) ) );
    }
    thisProxy[ r.first ].neworder( n );
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
  for (const auto& p : nodes) {
    auto id = std::get< 0 >( p.second );
    m_newnodes[ p.first ] = id;
    m_newcoordmap.emplace( id, std::get<1>(p.second) );
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
  for (auto& p : m_ginpoel) p = tk::cref_find( m_newnodes, p );

  // Update node coordinate map with the reordered IDs
  m_coordmap = m_newcoordmap;

  // Update symmetric chare-node communication map with the reordered IDs
  for (auto& c : m_msum) {
    decltype(c.second) n;
    for (auto p : c.second) n.insert( tk::cref_find(m_newnodes,p) );
    c.second = std::move( n );
  }

  // Update unique global node IDs of this chare with the reordered node IDs
  m_nodeset.clear();
  m_nodeset.insert( begin(m_ginpoel), end(m_ginpoel) );

  // Progress report to host
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) m_host.chreordered();

  // Compute lower and upper bounds of reordered node IDs on this chare
  bounds();
}

void
Sorter::bounds()
// *****************************************************************************
// Compute lower and upper bounds of reordered node IDs for this chare
//! \details This function computes the bounds that this chare will contribute
//!   to in a linear system solve. We find the largest node ID assigned on each
//!   chare by the reordering and use that as the upper global row index for
//!   this chare. Note that while this rarely results in equal number of rows
//!   assigned to chares, potentially resulting in some load-imbalance, it
//!   yields a pretty good division reducing communication costs during the
//!   assembly of the linear system, which is more important than a slight
//!   (FLOP) load imbalance. Since the upper index for chare 1 is the same as
//!   the lower index for chare 2, etc., we find the upper indices and then the
//!   lower indices for all chares are communicated.
// *****************************************************************************
{
  m_upper = *std::max_element( begin(m_nodeset), end(m_nodeset) );

  // The bounds are the dividers (global mesh point indices) at which the linear
  // system assembly is divided among PEs. However, Solver expect exclusive
  // upper indices, so we increase the last one by one here.
  if (thisIndex == m_nchare-1) ++m_upper;

  // Tell the runtime system that the upper bound has been computed
  upper_complete();

  // Set lower index for chare 0 as 0
  if (thisIndex == 0) lower( 0 );

  // All chares except the last one send their upper bound as the lower index for
  // the chare with thisIndex+1
  if (thisIndex < m_nchare-1) thisProxy[ thisIndex+1 ].lower( m_upper );
}

void
Sorter::lower( std::size_t low )
// *****************************************************************************
//  Receive lower bound of node IDs for this chare
//! \param[in] low Lower bound of node IDs assigned to this chare
// *****************************************************************************
{
  m_lower = low;
  lower_complete();
}

void
Sorter::create()
// *****************************************************************************
// Create chare array elements on this PE and assign the global mesh element IDs
// they will operate on
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::feedback >()) m_host.chbounds();

  if ( g_inputdeck.get< tag::discr, tag::scheme >() == ctr::SchemeType::MatCG)
    // broadcast this chare's bounds of global node IDs to matrix solvers
    m_solver.ckLocalBranch()->chbounds( m_lower, m_upper );
  else // if no MatCG, no matrix solver, continue
    createDiscWorkers();
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
  // Create worker array element using Charm++ dynamic chare array element
  // insertion: 1st arg: chare id, last arg: PE chare is created on, middle
  // args: Discretization ctor args. See also Charm++ manual, Sec. "Dynamic
  // Insertion".
  m_scheme.discInsert( thisIndex, m_host, m_ginpoel, m_coordmap, m_msum,
                       m_nchare, CkMyPe() );

  contribute( m_cbs.get< tag::discinserted >() );
}

void
Sorter::createWorkers()
// *****************************************************************************
//  Create worker chare array element
// *****************************************************************************
{
  // If there was no reordering, assign a one-to-one node-map
  if (m_newnodes.empty()) for (auto n : m_nodeset) m_newnodes[ n ] = n;

  // Extract this chare's portion of the boundary node lists
  decltype(m_bnode) chbnode;
  for (const auto& s : m_bnode) {
    auto& n = chbnode[ s.first ];
    for (auto p : s.second) {
      auto q = m_newnodes.find(p);
      if (q != end(m_newnodes)) n.push_back( q->second );
    }
    if (n.empty()) chbnode.erase( s.first );
  }
  // Make boundary node IDs unique for each physical boundary (side set)
  for (auto& s : chbnode) tk::unique( s.second );


  // Generate set of all mesh faces
  tk::UnsMesh::FaceSet faceset;
  for (std::size_t e=0; e<m_ginpoel.size()/4; ++e) { // for all tets
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f) // for all tet faces
      faceset.insert( {{{ m_ginpoel[ mark + tk::lpofa[f][0] ],
                          m_ginpoel[ mark + tk::lpofa[f][1] ],
                          m_ginpoel[ mark + tk::lpofa[f][2] ] }}} );
  }

  // Extract this chare's portion of the boundary faces and their connectivity
  decltype(m_bface) chbface;
  decltype(m_triinpoel) chtriinpoel;
  std::size_t cnt = 0;

  // Generate boundary 
  for (const auto& ss : m_bface)  // for all phsyical boundaries (sidesets)
    for (auto f : ss.second) {    // for all faces on this physical boundary
      // attempt to find face nodes on this chare
      auto f1 = m_newnodes.find( m_triinpoel[f*3+0] );
      auto f2 = m_newnodes.find( m_triinpoel[f*3+1] );
      auto f3 = m_newnodes.find( m_triinpoel[f*3+2] );
      // if all 3 nodes of the physical boundary face are on this chare
      if (f1!=end(m_newnodes) && f2!=end(m_newnodes) && f3!=end(m_newnodes)) {
        // Create face with new node ids (after mesh node reordering)
        std::array< std::size_t, 3 > n{{f1->second, f2->second, f3->second}};
        // if this boundary face is on this chare
        if (faceset.find(n) != end(faceset)) {
          // store face connectivity with new (global) node ids of this chare
          chtriinpoel.insert( end(chtriinpoel), begin(n), end(n) );
          // generate/store physical boundary face id associated to sideset id
          chbface[ ss.first ].push_back( cnt++ );
        }
      }
    }

  // Create face data
  FaceData fd( m_ginpoel, chbface, chbnode, chtriinpoel );

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
