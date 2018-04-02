// *****************************************************************************
/*!
  \file      src/Inciter/DG.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     DG advances a system of PDEs with the discontinuous Galerkin scheme
  \details   DG advances a system of partial differential equations (PDEs) using
    discontinuous Galerkin (DG) finite element (FE) spatial discretization (on
    tetrahedron elements) combined with Runge-Kutta (RK) time stepping.
  \see The documentation in DG.h.
*/
// *****************************************************************************

#include <algorithm>
#include <numeric>

#include "DG.h"
#include "Discretization.h"
#include "DGPDE.h"
#include "Solver.h"
#include "DiagReducer.h"
#include "DerivedData.h"
#include "ElemDiagnostics.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "ExodusIIMeshWriter.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< DGPDE > g_dgpde;

} // inciter::

using inciter::DG;

DG::DG( const CProxy_Discretization& disc,
        const tk::CProxy_Solver&,
        const FaceData& fd ) :
  m_ncomfac( 0 ),
  m_nadj( 0 ),
  m_nsol( 0 ),
  m_itf( 0 ),
  m_disc( disc ),
  m_fd( fd ),
  m_u( m_disc[thisIndex].ckLocal()->Inpoel().size()/4,
       g_inputdeck.get< tag::component >().nprop() ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_vol( 0.0 ),
  m_geoFace( tk::genGeoFaceTri( fd.Ntfac(), fd.Inpofa(),
                                m_disc[thisIndex].ckLocal()->Coord()) ),
  m_geoElem( tk::genGeoElemTet( m_disc[thisIndex].ckLocal()->Inpoel(),
                                m_disc[thisIndex].ckLocal()->Coord() ) ),
  m_lhs( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_nfac( fd.Inpofa().size()/3 ),
  m_nunk( m_u.nunk() ),
  m_msumset( msumset() ),
  m_esuelTet( tk::genEsuelTet( m_disc[thisIndex].ckLocal()->Inpoel(),
                tk::genEsup( m_disc[thisIndex].ckLocal()->Inpoel(), 4 ) ) ),
  m_bndFace(),
  m_ghostData(),
  m_ghostReq( 0 ),
  m_exptNbface( 0 ),
  m_ghost(),
  m_exptGhost(),
  m_recvGhost(),
  m_diag()
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] Face data structures
// *****************************************************************************
{
  // Perform leak test on mesh partition
  Assert( !leakyPartition(), "Mesh partition leaky" );

  // Activate SDAG waits for face adjacency map (ghost data) calculation
  wait4ghost();

  auto d = Disc();

  const auto& gid = d->Gid();
  const auto& inpoel = d->Inpoel();
  const auto& inpofa = fd.Inpofa();

  // Invert inpofa to enable searching for faces based on (global) node triplets
  Assert( inpofa.size() % 3 == 0, "Inpofa must contain triplets" );
  for (std::size_t f=0; f<inpofa.size()/3; ++f)
    m_ipface.insert( {{{ gid[ inpofa[f*3+0] ],
                         gid[ inpofa[f*3+1] ],
                         gid[ inpofa[f*3+2] ] }}} );

  // At this point ipface has node-id-triplets (faces) on the internal
  // chare-domain and on the physical boundary but not on chare boundaries,
  // hence the name internal + physical boundary faces.

  // Build a set of faces (each face given by 3 global node IDs) associated to
  // chares we potentially share boundary faces with.
  tk::UnsMesh::FaceSet potbndface;
  for (std::size_t e=0; e<m_esuelTet.size()/4; ++e) {   // for all our tets
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f)     // for all tet faces
      if (m_esuelTet[mark+f] == -1) {   // if face has no outside-neighbor tet
        // if does not exist in among the internal and physical boundary faces,
        // store as a potential chare-boundary face
        tk::UnsMesh::Face t{{ gid[ inpoel[ mark + tk::lpofa[f][0] ] ],
                              gid[ inpoel[ mark + tk::lpofa[f][1] ] ],
                              gid[ inpoel[ mark + tk::lpofa[f][2] ] ] }};
        if (m_ipface.find(t) == end(m_ipface)) {
          Assert( ++m_exptNbface, "Sum up expected number of boundary faces" );
          potbndface.insert( t );
        }
      }
  }

  // In the following we assume that the size of the (potential) boundary-face
  // adjacency map above does not necessarily equal to that of the node
  // adjacency map (m_msumset). This is because while a node can be shared at a
  // single corner or along an edge, that does not necessarily share a face as
  // well (in other words, shared nodes or edges can exist that are not part of
  // a shared face). So the chares we communicate with across faces are not
  // necessarily the same as the chares we would communicate nodes with.
  //
  // Since the sizes of the node and face adjacency maps are not the same, while
  // sending the faces on chare boundaries would be okay, however, the receiver
  // would not necessarily know how many chares it must receive from. To solve
  // this problem we send to chares which we share at least a single node with,
  // i.e., rely on the node-adjacency map. Note that to all chares we share at
  // least a single node with we send all our potential chare-boundary faces.
  // This is the same list of faces to all chares we send.
  //
  // Another underlying assumption here is, of course, that the size of the face
  // adjacency map is always smaller than or equal to that of the node adjacency
  // map, which is always true. Since the receive side already knows how many
  // fellow chares it must receive shared node ids from, we use that to detect
  // completion of the number of receives in comfac(). This simplifies the
  // communication pattern and code.

  // Send sets of faces adjacent to chare boundaries to fellow workers (if any)
  if (m_msumset.empty())        // in serial, skip setting up ghosts altogether
    adj();
  else
    for (const auto& c : m_msumset) {   // for all chares we share nodes with
      thisProxy[ c.first ].comfac( thisIndex, potbndface );
    }
}

bool
DG::leakyPartition()
// *****************************************************************************
// Perform leak-test on mesh partition
//! \details This function computes a surface integral over the boundary of the
//!   incoming mesh partition. A non-zero vector result indicates a leak, e.g.,
//!   a hole in the partition, which indicates an error upstream of this code,
//!   either in the mesh geometry, mesh partitioning, or in the data structures
//!   that represent faces.
//! \return True if our chare partition leaks.
// *****************************************************************************
{
  auto d = Disc();
  const auto& inpoel = d->Inpoel();
  const auto& coord = d->Coord();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Storage for surface integral over our mesh partition
  std::array< tk::real, 3 > s{{ 0.0, 0.0, 0.0}};

  for (std::size_t e=0; e<m_esuelTet.size()/4; ++e) {   // for all our tets
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f)     // for all tet faces
      if (m_esuelTet[mark+f] == -1) {   // if face has no outside-neighbor tet
        // 3 local node IDs of face
        auto A = inpoel[ mark + tk::lpofa[f][0] ];
        auto B = inpoel[ mark + tk::lpofa[f][1] ];
        auto C = inpoel[ mark + tk::lpofa[f][2] ];
        // Compute geometry data for face
        auto geoface = tk::geoFaceTri( {{x[A], x[B], x[C]}},
                                       {{y[A], y[B], y[C]}},
                                       {{z[A], z[B], z[C]}} );
        // Sum up face area * face unit-normal
        s[0] += geoface(0,0,0) * geoface(0,1,0);
        s[1] += geoface(0,0,0) * geoface(0,2,0);
        s[2] += geoface(0,0,0) * geoface(0,3,0);
      }
  }

  auto eps = std::numeric_limits< tk::real >::epsilon() * 100;
  return std::abs(s[0]) > eps || std::abs(s[1]) > eps || std::abs(s[2]) > eps;
}

bool
DG::leakyAdjacency()
// *****************************************************************************
// Perform leak-test on chare boundary faces
//! \details This function computes a surface integral over the boundary of the
//!   faces after the face adjacency communication map is completed.. A non-zero
//!   vector result indicates a leak, e.g., a hole in the partition (covered by
//!   the faces of the face adjacency communication map), which indicates an
//!   error upstream in the code that sets up the face communication data
//!   structures.
//! \note Compared to leakyPartition() this function performs the leak-test on
//!   the face geometry data structure enlarged by ghost faces on this partition
//!   by computing a discrete surface integral considering the physical and
//!   chare boundary faces, which should be equal to zero for a closed domain.
//! \return True if our chare face adjacency leaks.
// *****************************************************************************
{
  // Storage for surface integral over our chunk of the adjacency
  std::array< tk::real, 3 > s{{ 0.0, 0.0, 0.0}};

  // physical boundary faces
  for (std::size_t f=0; f<m_fd.Nbfac(); ++f) {
    s[0] += m_geoFace(f,0,0) * m_geoFace(f,1,0);
    s[1] += m_geoFace(f,0,0) * m_geoFace(f,2,0);
    s[2] += m_geoFace(f,0,0) * m_geoFace(f,3,0);
  }

  // chare-boundary faces
  for (std::size_t f=m_fd.Ntfac(); f<m_fd.Esuf().size()/2; ++f) {
    s[0] += m_geoFace(f,0,0) * m_geoFace(f,1,0);
    s[1] += m_geoFace(f,0,0) * m_geoFace(f,2,0);
    s[2] += m_geoFace(f,0,0) * m_geoFace(f,3,0);
  }

  auto eps = std::numeric_limits< tk::real >::epsilon() * 100;
  return std::abs(s[0]) > eps || std::abs(s[1]) > eps || std::abs(s[2]) > eps;
}

std::unordered_map< int, std::unordered_set< std::size_t > >
DG::msumset() const
// *****************************************************************************
// Convert chare-node adjacency map to hold sets instead of vectors
//! \return Chare-node adjacency map that holds sets instead of vectors
// *****************************************************************************
{
  auto d = Disc();

  std::unordered_map< int, std::unordered_set< std::size_t > > m;
  for (const auto& n : d->Msum())
    m[ n.first ].insert( n.second.cbegin(), n.second.cend() );

  Assert( m.find( thisIndex ) == m.cend(),
          "Chare-node adjacency map should not contain data for own chare ID" );

  return m;
}

void
DG::comfac( int fromch, const tk::UnsMesh::FaceSet& infaces )
// *****************************************************************************
//  Receive unique set of faces we potentially share with/from another chare
//! \param[in] fromch Sender chare id
//! \param[in] infaces Unique set of faces we potentially share with fromch
// *****************************************************************************
{
  // Attempt to find sender chare among chares we potentially share faces with.
  // Note that it is feasible that a sender chare called us but we do not have a
  // set of faces associated to that chare. This can happen if we only share a
  // single node or an edge but note a face with that chare.
  auto& bndface = m_bndFace[ fromch ];  // will associate to sender chare
  // Try to find incoming faces on our chare boundary with other chares. If
  // found, generate and assign new local ID to face, associated to sender
  // chare.
  auto d = Disc();
  const auto& gid = d->Gid();
  const auto& inpoel = d->Inpoel();
  for (std::size_t e=0; e<m_esuelTet.size()/4; ++e) {  // for all our tets
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f) {  // for all cell faces
      if (m_esuelTet[mark+f] == -1) {  // if face has no outside-neighbor tet
        tk::UnsMesh::Face t{{ gid[ inpoel[ mark + tk::lpofa[f][0] ] ],
                              gid[ inpoel[ mark + tk::lpofa[f][1] ] ],
                              gid[ inpoel[ mark + tk::lpofa[f][2] ] ] }};
        // if found among the incoming faces and if not one of our internal nor
        // physical boundary faces
        if ( infaces.find(t) != end(infaces) &&
             m_ipface.find(t) == end(m_ipface) )
          bndface[t][0] = m_nfac++;    // assign new local face ID
      }
    }
  }
  // If at this point we have not found any face among our faces we potentially
  // share with fromch, there is no need to keep an empty set of faces
  // associated to fromch as we only share nodes or edges with it, but not
  // faces.
  if (bndface.empty()) m_bndFace.erase( fromch );

  // if we have heard from all fellow chares that we share at least a single
  // node, edge, or face with
  if (++m_ncomfac == m_msumset.size()) {

    tk::destroy(m_ipface);

    // Error checking on the number of expected vs received/found chare-boundary
    // faces
    Assert( m_exptNbface ==
             std::accumulate( m_bndFace.cbegin(), m_bndFace.cend(),
               std::size_t(0),
               []( std::size_t acc, const decltype(m_bndFace)::value_type& b )
                 { return acc + b.second.size(); } ),
            "Expected and received number of boundary faces mismatch" );

    // Basic error checking on chare-boundary-face map
    Assert( m_bndFace.find( thisIndex ) == m_bndFace.cend(),
            "Face-communication map should not contain data for own chare ID" );

    // Store (local) tet ID adjacent to our chare boundary from the inside
    for (std::size_t e=0; e<m_esuelTet.size()/4; ++e) {  // for all our tets
      auto mark = e*4;
      for (std::size_t f=0; f<4; ++f) {  // for all cell faces
        if (m_esuelTet[mark+f] == -1) {  // if face has no outside-neighbor tet
          tk::UnsMesh::Face t{{ gid[ inpoel[ mark + tk::lpofa[f][0] ] ],
                                gid[ inpoel[ mark + tk::lpofa[f][1] ] ],
                                gid[ inpoel[ mark + tk::lpofa[f][2] ] ] }};
          auto c = findchare( t );
          if (c > -1) {
            auto& lbndface = tk::ref_find( m_bndFace, c );
            auto& face = tk::ref_find( lbndface, t );
            face[1] = e;  // store (local) inner tet ID adjacent to face
          }
        }
      }
    }

    // At this point m_bndFace is complete on this PE. This means that starting
    // from the sets of faces we potentially share with fellow chares we now
    // only have those faces we actually share faces with (through which we need
    // to communicate later). Also, m_bndFace not only has the unique faces
    // associated to fellow chares, but also a newly assigned local face ID as
    // well as the local id of the inner tet adjacent to the face. Continue by
    // starting setting up ghost data
    setupGhost();
    // Besides setting up our own ghost data, we also issue requests (for ghost
    // data) to those chares which we share faces with. Note that similar to
    // comfac() we are calling reqGhost() by going through msumset instead,
    // which may send requests to those chare we do not share faces with. This
    // is so that we can test for completing by querying the size of the already
    // complete msumset in reqGhost. Requests in sendGhost will only be
    // fullfilled based on m_ghostData.
    for (const auto& c : m_msumset)     // for all chares we share nodes with
      thisProxy[ c.first ].reqGhost();
  }
}

void
DG::setupGhost()
// *****************************************************************************
// Setup own ghost data on this chare
// *****************************************************************************
{
  auto d = Disc();
  const auto& gid = d->Gid();
  const auto& inpoel = d->Inpoel();

  // Enlarge elements surrounding faces data structure for ghosts
  m_fd.Esuf().resize( 2*m_nfac, -1 );
  // Enlarge face geometry data structure for ghosts
  m_geoFace.resize( m_nfac, 0.0 );

  // Collect tet ids, their face connectivity (given by 3 global node IDs, each
  // triplet for potentially multiple faces on the chare boundary), and their
  // elem geometry data (see GhostData) associated to fellow chares adjacent to
  // chare boundaries. Once received by fellow chares, these tets will become
  // known as ghost elements and their data as ghost data.
  for (std::size_t e=0; e<m_esuelTet.size()/4; ++e) {  // for all our tets
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f) {  // for all cell faces
      if (m_esuelTet[mark+f] == -1) {  // if face has no outside-neighbor tet
        tk::UnsMesh::Face t{{ gid[ inpoel[ mark + tk::lpofa[f][0] ] ],
                              gid[ inpoel[ mark + tk::lpofa[f][1] ] ],
                              gid[ inpoel[ mark + tk::lpofa[f][2] ] ] }};
        auto c = findchare( t );
        // It is possible that we do not find the chare for this face. We are
        // looping through all of our tets and interrogating all faces that do
        // not have neighboring tets but we only care about chare-boundary faces
        // here as only those need ghost data. (m_esuelTet may also contain
        // physical boundary faces.)
        if (c > -1) {
          // Will store ghost data associated to neighbor chare
          auto& ghost = m_ghostData[ c ];
          // Store tet id adjacent to chare boundary as key for ghost data
          auto& tuple = ghost[ e ];
          // If tetid e has not yet been encountered, store geometry (only once)
          auto& nodes = std::get< 0 >( tuple );
          if (nodes.empty()) std::get< 1 >( tuple ) = m_geoElem[ e ];
          // (Always) store face node IDs on chare boundary, even if tetid e has
          // already been stored. Thus we store potentially multiple faces along
          // the same chare-boundary. This happens, e.g., when the boundary
          // between chares is zig-zaggy enough to have 2 or even 3 faces of the
          // same tet.
          nodes.push_back( t[0] );
          nodes.push_back( t[1] );
          nodes.push_back( t[2] );
          Assert( nodes.size() <= 4*3, "Overflow of faces/tet to send" );
        }
      }
    }
  }

  // Basic error checking on local ghost data
  Assert( m_ghostData.find( thisIndex ) == m_ghostData.cend(),
          "Chare-node adjacency map should not contain data for own chare ID" );

  // More in-depth error checking on local ghost data
  for (const auto& c : m_ghostData)
    for (const auto& t : c.second) {
      IGNORE(t);
      Assert( !std::get< 0 >( t.second ).empty(),
              "Emtpy face vector in ghost data" );
      Assert( std::get< 0 >( t.second ).size() % 3 == 0,
              "Face node IDs must be triplets" );
      Assert( std::get< 0 >( t.second ).size() <= 4*3,    // <= 4*3 (4*numfaces)
              "Max number of faces for a single ghost tet is 4" );
      Assert( !std::get< 1 >( t.second ).empty(),
              "No elem geometry data for ghost" );
      Assert( std::get< 1 >( t.second ).size() == m_geoElem.nprop(),
              "Elem geometry data for ghost must be for single tet" );
    }

  ownghost_complete();
}

void
DG::reqGhost()
// *****************************************************************************
// Receive requests for ghost data
// *****************************************************************************
{
  // If every chare we communicate with has requested ghost data from us, we may
  // fulfill the requests, but only if we have already setup our ghost data.
  if (++m_ghostReq == m_msumset.size()) reqghost_complete();
}

void
DG::sendGhost()
// *****************************************************************************
// Send all of our ghost data to fellow chares
// *****************************************************************************
{
  for (const auto& c : m_ghostData)
    thisProxy[ c.first ].comGhost( thisIndex, c.second );
}

int
DG::findchare( const tk::UnsMesh::Face& t )
// *****************************************************************************
// Find any chare for face (given by 3 global node IDs)
//! \param[in] t Face given by three global node IDs
//! \return Chare ID if found among any of the chares we communicate along
//!   faces with, -1 if the face cannot be found.
// *****************************************************************************
{
  for (const auto& cf : m_bndFace)
    if (cf.second.find(t) != end(cf.second))
      return cf.first;
  return -1;
}

void
DG::comGhost( int fromch, const GhostData& ghost )
// *****************************************************************************
// Receive ghost data on chare boundaries from fellow chare
//! \param[in] fromch Caller chare ID
//! \param[in] ghost Ghost data, see Inciter/FaceData.h for the type
// *****************************************************************************
{
  // nodelist with fromch, currently only used for an assert
  const auto& nl = tk::cref_find( m_msumset, fromch );
  IGNORE(nl);

  auto& ghostelem = m_ghost[ fromch ];  // will associate to sender chare

  // Store ghost data coming from chare
  for (const auto& g : ghost) {  // loop over incoming ghost data
    auto e = g.first;  // remote/ghost tet id outside of chare boundary
    const auto& nodes = std::get< 0 >( g.second );  // node IDs of face(s)
    const auto& geo = std::get< 1 >( g.second );    // ghost elem geometry data
    Assert( nodes.size() % 3 == 0, "Face node IDs must be triplets" );
    Assert( nodes.size() <= 4*3, "Overflow of faces/tet received" );
    Assert( geo.size() % 4 == 0, "Ghost geometry size mismatch" );
    Assert( geo.size() == m_geoElem.nprop(), "Ghost geometry number mismatch" );
    for (std::size_t n=0; n<nodes.size()/3; ++n) {  // face(s) of ghost e
      // node IDs of face on chare boundary
      tk::UnsMesh::Face t{{ nodes[n*3+0], nodes[n*3+1], nodes[n*3+2] }};
      // must find t in nodelist of chare-boundary adjacent to fromch
      Assert( nl.find(t[0]) != end(nl) &&
              nl.find(t[1]) != end(nl) &&
              nl.find(t[2]) != end(nl),
           "Ghost face not found in chare-node adjacency map on receiving end" );
      // must find face in boundary-face adjacency map for fromch
      Assert( tk::cref_find(m_bndFace,fromch).find( t ) !=
              tk::cref_find(m_bndFace,fromch).cend(), "Ghost face not "
              "found in boundary-face adjacency map on receiving end" );
      // find local face & tet ids for t
      auto id = tk::cref_find( tk::cref_find(m_bndFace,fromch), t );
      // compute face geometry for chare-boundary face
      addGeoFace( t, id );
      // if ghost tet id not yet encountered on boundary with fromch
      auto i = ghostelem.find( e );
      if (i != end(ghostelem))
        addEsuf( id, i->second );  // fill in elements surrounding face
      else {
        addEsuf( id, m_nunk );     // fill in elements surrounding face
        ghostelem[e] = m_nunk;     // assign new local tet id to remote ghost id
        m_geoElem.push_back( geo );// store ghost elem geometry
        ++m_nunk;                  // increase number of unknowns on this chare
      }
    }
  }

  // Signal the runtime system that all workers have received their adjacency
  if (++m_nadj == m_ghostData.size()) adj();
}

void
DG::addEsuf( const std::array< std::size_t, 2 >& id, std::size_t ghostid )
// *****************************************************************************
// Fill elements surrounding a face along chare boundary
//! \param[in] id Local face and (inner) tet id adjacent to it
//! \param[in] ghostid Local ID for ghost tet
//! \details This function extends and fills in the elements surrounding faces
//!   data structure (esuf) so that the left and right element id is filled
//!   in correctly on chare boundaries to contain the correct inner tet id and
//!   the local tet id for the outer (ghost) tet, both adjacent to the given
//1   chare-face boundary. Prior to this function, this data structure does not
//!   have yet face-element connectivity adjacent to chare-boundary faces, only
//!   for physical boundaries and internal faces that are not on the chare
//!   boundary (this latter purely as a result of mesh partitioning). The remote
//!   element id of the ghost is stored in a location that is local to our own
//!   esuf. The face numbering is such that esuf stores the element-face
//!   connectivity first for the physical-boundary faces, followed by that of
//!   the internal faces, followed by the chare-boundary faces. As a result,
//!   esuf can be used by physics algorithms in exactly the same way as would be
//!   used in serial. In serial, of course, this data structure is not extended
//!   at the end by the chare-boundaries.
// *****************************************************************************
{
  auto& esuf = m_fd.Esuf();
  Assert( 2*id[0]+1 < esuf.size(), "Indexing out of esuf" );

  // put in inner tet id
  esuf[ 2*id[0]+0 ] = static_cast< int >( id[1] );
  // put in local id for outer/ghost tet
  esuf[ 2*id[0]+1 ] = static_cast< int >( ghostid );
}

void
DG::addGeoFace( const tk::UnsMesh::Face& t,
                const std::array< std::size_t, 2 >& id )
// *****************************************************************************
// Fill face-geometry data along chare boundary
//! \param[in] t Face (given by 3 global node IDs) on the chare boundary
//! \param[in] id Local face and (inner) tet id adjacent to face t
//! \details This function fills in the face geometry data along a chare
//!    boundary.
// *****************************************************************************
{
  auto d = Disc();
  auto coord = d->Coord();
  auto lid = d->Lid();

  // get global node IDs reversing order to get outward-pointing normal
  auto A = tk::cref_find( lid, t[2] );
  auto B = tk::cref_find( lid, t[1] );
  auto C = tk::cref_find( lid, t[0] );
  auto geochf = tk::geoFaceTri( {{coord[0][A], coord[0][B], coord[0][C]}},
                                {{coord[1][A], coord[1][B], coord[1][C]}},
                                {{coord[2][A], coord[2][B], coord[2][C]}} );

  for (std::size_t i=0; i<7; ++i)
    m_geoFace(id[0],i,0) = geochf(0,i,0);
}

void
DG::adj()
// *****************************************************************************
// Continue after face adjacency communication map completed on this chare
//! \details At this point the face/ghost communication map has been established
//!    on this chare.
// *****************************************************************************
{
  tk::destroy(m_bndFace);

  // Ensure that all elements surrounding faces (are correct) including those at
  // chare boundaries
  for (std::size_t f=0; f<m_nfac; ++f) {
    Assert( m_fd.Esuf()[2*f] > -1,
            "Left element in esuf cannot be physical ghost" );
    if (f >= m_fd.Nbfac())
      Assert( m_fd.Esuf()[2*f+1] > -1,
           "Right element in esuf for internal/chare faces cannot be a ghost" );
  }

  // Error checking on ghost data
  for(const auto& n : m_ghostData)
    for(const auto& i : n.second) {
      Assert( i.first < m_esuelTet.size()/4, "Sender contains ghost tet id " );
      IGNORE(i);
    }

  // Perform leak test on face geometry data structure enlarged by ghosts
  Assert( !leakyAdjacency(), "Face adjacency leaky" );

  // Resize solution vectors, lhs, and rhs by the number of ghost tets
  m_u.resize( m_nunk );
  m_un.resize( m_nunk );
  m_lhs.resize( m_nunk );
  m_rhs.resize( m_nunk );

  // Ensure that we also have all the geometry data (including those of ghosts)
  Assert( m_geoElem.nunk() == m_u.nunk(), "GeoElem unknowns size mismatch" );

  // Basic error checking on ghost tet ID map
  Assert( m_ghost.find( thisIndex ) == m_ghost.cend(),
          "Ghost id map should not contain data for own chare ID" );

  // Store expected ghost tet IDs
  for (const auto& c : m_ghost)
    for (const auto& g : c.second) {
      IGNORE(g);
      Assert( m_exptGhost.insert( g.second ).second,
              "Failed to store local tetid as exptected ghost id" );
    }

  // Signal the runtime system that all workers have received their adjacency
  auto d = Disc();
  d->contribute(CkCallback(CkReductionTarget(Transporter,comfinal), d->Tr()));
}

void
DG::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types
//! \details Since this is a [nodeinit] routine, the runtime system executes the
//!   routine exactly once on every logical node early on in the Charm++ init
//!   sequence. Must be static as it is called without an object. See also:
//!   Section "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  ElemDiagnostics::registerReducers();
}

void
DG::setup( tk::real v )
// *****************************************************************************
// Setup rows, query boundary conditions, output mesh, etc.
//! \param[in] v Total mesh volume
// *****************************************************************************
{
  tk::destroy(m_msumset);

  auto d = Disc();

  // Store total mesh volume
  m_vol = v;
  // Output chare mesh to file
  d->writeMesh();
  // Output fields metadata to output file
  d->writeElemMeta();

  // Basic error checking on element geometry data size
  Assert( m_geoElem.nunk() == m_lhs.nunk(), "Size mismatch in DG::setup()" );

  // Compute left-hand side of discrete PDEs
  lhs();

  // Set initial conditions for all PDEs
  for (const auto& eq : g_dgpde) eq.initialize( m_geoElem, m_u,  d->T() );
  m_un = m_u;

  // Output initial conditions to file (regardless of whether it was requested)
  if ( !g_inputdeck.get< tag::cmd, tag::benchmark >() ) writeFields( d->T() );

  // Start timer measuring time stepping wall clock time
  d->Timer().zero();

  // Start time stepping
  dt();
}

void
DG::dt()
// *****************************************************************************
// Compute time step size
// *****************************************************************************
{
  auto mindt = std::numeric_limits< tk::real >::max();

  auto const_dt = g_inputdeck.get< tag::discr, tag::dt >();
  auto def_const_dt = g_inputdeck_defaults.get< tag::discr, tag::dt >();
  auto eps = std::numeric_limits< tk::real >::epsilon();

  // use constant dt if configured
  if (std::abs(const_dt - def_const_dt) > eps) {

    mindt = const_dt;

  } else {      // compute dt based on CFL

    // find the minimum dt across all PDEs integrated
    // ...
    mindt = 0.1;        // stub for now to overwrite numeric_limits::max

    // Scale smallest dt with CFL coefficient
    mindt *= g_inputdeck.get< tag::discr, tag::cfl >();

  }

  // Enable SDAG wait for building the solution vector
  wait4sol();

  // Contribute to minimum dt across all chares then advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(DG,advance), thisProxy) );
}

void
DG::advance( tk::real newdt )
// *****************************************************************************
// Advance equations to next time step
//! \param[in] newdt Size of this new time step
// *****************************************************************************
{
  auto d = Disc();

  // Set new time step size
  d->setdt( newdt );

  // communicate solution ghost data (if any)
  if (m_ghostData.empty())
    comsol_complete();
  else
    for(const auto& n : m_ghostData) {
      std::vector< std::size_t > tetid;
      std::vector< std::vector< tk::real > > u;
      for(const auto& i : n.second) {
        Assert( i.first < m_esuelTet.size()/4, "Sending solution ghost data" );
        tetid.push_back( i.first );
        u.push_back( m_u[i.first] );
      }
      thisProxy[ n.first ].comsol( thisIndex, tetid, u );
    }

  ownsol_complete();
}

void
DG::comsol( int fromch,
            const std::vector< std::size_t >& tetid,
            const std::vector< std::vector< tk::real > >& u )
// *****************************************************************************
//  Receive chare-boundary solution ghost data from neighboring chares
//! \param[in] tetid Ghost tet ids we receive solution data for
//! \param[in] u Solution ghost data
//! \details This function receives contributions to m_u from fellow chares.
// *****************************************************************************
{
  Assert( u.size() == tetid.size(), "Size mismatch in DG::comsol()" );

  // Find local-to-ghost tet id map for sender chare
  const auto& n = tk::cref_find( m_ghost, fromch );

  for (std::size_t i=0; i<tetid.size(); ++i) {
    auto j = tk::cref_find( n, tetid[i] );
    Assert( j >= m_esuelTet.size()/4, "Receiving solution non-ghost data" );
    Assert( j < m_u.nunk(), "Indexing out of bounds in DG::comsol()" );
    Assert( m_recvGhost.insert( j ).second,
            "Failed to store local tetid of received ghost tetid" );
    for (std::size_t c=0; c<m_u.nprop(); ++c)
      m_u(j,c,0) = u[i][c];
  }

  // if we have received all solution ghost contributions from those chares we
  // communicate along chare-boundary faces with, solve the system
  if (++m_nsol == m_ghostData.size()) {
    Assert( m_exptGhost == m_recvGhost,
            "Expected/received ghost tet id mismatch" );
    m_recvGhost.clear();
    m_nsol = 0;
    comsol_complete();
  }
}

void
DG::writeFields( tk::real time )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] time Physical time
// *****************************************************************************
{
  auto d = Disc();

  // Save time stamp at which the last field write happened
  d->LastFieldWriteTime() = time;

  // Increase field output iteration count
  ++m_itf;

  // Collect element field output
  std::vector< std::vector< tk::real > > elemfields;
  auto u = m_u;   // make a copy as eq::output() may overwrite its arg
  for (const auto& eq : g_dgpde) {
    auto output = eq.fieldOutput( time, m_vol, m_geoElem, u );
    // cut off ghost elements
    for (auto& o : output) o.resize( m_esuelTet.size()/4 );
    elemfields.insert( end(elemfields), begin(output), end(output) );
  }

  // Create ExodusII writer
  tk::ExodusIIMeshWriter ew( d->OutFilename(), tk::ExoWriter::OPEN );
  // Write time stamp
  ew.writeTimeStamp( m_itf, time );
  // Write node fields to file
  d->writeElemSolution( ew, m_itf, elemfields );
}

void
DG::out()
// *****************************************************************************
// Output mesh field data
// *****************************************************************************
{
  auto d = Disc();

  // Optionally output field and particle data
  if ( !((d->It()+1) % g_inputdeck.get< tag::interval, tag::field >()) &&
       !g_inputdeck.get< tag::cmd, tag::benchmark >() )
  {
    writeFields( d->T() + d->Dt() );
  }

  // Output final field data to file (regardless of whether it was requested)
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  if ( (std::fabs(d->T() + d->Dt() - term) < eps || (d->It()+1) >= nstep) &&
       (!g_inputdeck.get< tag::cmd, tag::benchmark >()) )
  {
    writeFields( d->T()+d->Dt() );
  }
}

void
DG::lhs()
// *****************************************************************************
// Compute left-hand side of discrete transport equations
// *****************************************************************************
{
  // Compute left-hand side matrix for all equations
  for (const auto& eq : g_dgpde)
    eq.lhs( m_geoElem, m_lhs );
}

void
DG::solve()
// *****************************************************************************
// Compute right-hand side of discrete transport equations
// *****************************************************************************
{
  auto d = Disc();

  for (const auto& eq : g_dgpde)
    eq.rhs( d->T(), m_geoFace, m_fd, m_u, m_rhs );

  // Explicit time-stepping using forward Euler to discretize time-derivative
  //m_u = m_un + d->Dt() * m_rhs/m_lhs;
  //m_un = m_u;
  m_u += d->Dt() * m_rhs/m_lhs;

  // Output field data to file
  out();
  // Compute diagnostics, e.g., residuals
  auto diag = m_diag.compute( *d, m_u.nunk()-m_esuelTet.size()/4, m_geoElem, m_u );
  // Increase number of iterations and physical time
  d->next();
  // Output one-liner status report
  d->status();

  // Evaluate whether to continue with next step
  if (!diag) eval();
}

void
DG::eval()
// *****************************************************************************
// Evaluate whether to continue with next step
// *****************************************************************************
{
  auto d = Disc();

  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();

  // If neither max iterations nor max time reached, continue, otherwise finish
  if (std::fabs(d->T()-term) > eps && d->It() < nstep)
    dt();
  else
    contribute( CkCallback( CkReductionTarget(Transporter,finish), d->Tr() ) );
}

#include "NoWarning/dg.def.h"
