// *****************************************************************************
/*!
  \file      src/Inciter/DG.C
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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
#include "DiagReducer.h"
#include "DerivedData.h"
#include "ElemDiagnostics.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Refiner.h"
#include "Limiter.h"
#include "Reorder.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< DGPDE > g_dgpde;

//! Runge-Kutta coefficients
static const std::array< std::array< tk::real, 3 >, 2 >
  rkcoef{{ {{ 0.0, 3.0/4.0, 1.0/3.0 }}, {{ 1.0, 1.0/4.0, 2.0/3.0 }} }};

} // inciter::

using inciter::DG;

DG::DG( const CProxy_Discretization& disc,
        const std::map< int, std::vector< std::size_t > >& bface,
        const std::map< int, std::vector< std::size_t > >& /* bnode */,
        const std::vector< std::size_t >& triinpoel ) :
  m_disc( disc ),
  m_ncomfac( 0 ),
  m_nadj( 0 ),
  m_nsol( 0 ),
  m_ninitsol( 0 ),
  m_nlim( 0 ),
  m_fd( Disc()->Inpoel(), bface, tk::remap(triinpoel,Disc()->Lid()) ),
  m_u( Disc()->Inpoel().size()/4,
       g_inputdeck.get< tag::discr, tag::ndof >()*
       g_inputdeck.get< tag::component >().nprop() ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_geoFace( tk::genGeoFaceTri( m_fd.Nipfac(), m_fd.Inpofa(), Disc()->Coord()) ),
  m_geoElem( tk::genGeoElemTet( Disc()->Inpoel(), Disc()->Coord() ) ),
  m_lhs( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_limFunc( limFunc(Disc()->Inpoel().size()/4) ),
  m_nfac( m_fd.Inpofa().size()/3 ),
  m_nunk( m_u.nunk() ),
  m_ncoord( Disc()->Coord()[0].size() ),
  m_msumset( Disc()->msumset() ),
  m_bndFace(),
  m_ghostData(),
  m_ghostReq( 0 ),
  m_exptNbface( 0 ),
  m_ghost(),
  m_exptGhost(),
  m_recvGhost(),
  m_diag(),
  m_stage( 0 ),
  m_initial( 1 ),
  m_refined( 0 )
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
// *****************************************************************************
{
  usesAtSync = true;    // enable migration at AtSync

  // Size communication buffers and setup ghost data
  resizeComm();
}

void
DG::resizeComm()
// *****************************************************************************
//  Start sizing communication buffers and setting up ghost data
// *****************************************************************************
{
  auto d = Disc();

  const auto& gid = d->Gid();
  const auto& inpoel = d->Inpoel();
  const auto& inpofa = m_fd.Inpofa();
  const auto& esuel = m_fd.Esuel();

  // Perform leak test on mesh partition
  Assert( !tk::leakyPartition( esuel, inpoel, d->Coord() ),
          "Mesh partition leaky" );

  // Activate SDAG waits for face adjacency map (ghost data) calculation
  thisProxy[ thisIndex ].wait4ghost();

  // Enable SDAG wait for initially building the solution vector and limiting
  if (m_initial) {
    thisProxy[ thisIndex ].wait4initlim();
    thisProxy[ thisIndex ].wait4sol();
    thisProxy[ thisIndex ].wait4lim();
  }

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
  for (std::size_t e=0; e<esuel.size()/4; ++e) {   // for all our tets
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f)     // for all tet faces
      if (esuel[mark+f] == -1) {        // if face has no outside-neighbor tet
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

  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) d->Tr().chbndface();

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

tk::Fields
DG::limFunc( std::size_t nelem ) const
// *****************************************************************************
// Size and create limiter function data container
//! \param[in] nelem Number elements in mesh
// *****************************************************************************
{
  auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
  auto nprop = g_inputdeck.get< tag::component >().nprop();

  return tk::Fields( (ndof == 1 ? 0 : nelem), (ndof-1)*nprop );
}

bool
DG::leakyAdjacency()
// *****************************************************************************
// Perform leak-test on chare boundary faces
//! \details This function computes a surface integral over the boundary of the
//!   faces after the face adjacency communication map is completed. A non-zero
//!   vector result indicates a leak, e.g., a hole in the partition (covered by
//!   the faces of the face adjacency communication map), which indicates an
//!   error upstream in the code that sets up the face communication data
//!   structures.
//! \note Compared to tk::leakyPartition() this function performs the leak-test
//!   on the face geometry data structure enlarged by ghost faces on this
//!   partition by computing a discrete surface integral considering the
//!   physical and chare boundary faces, which should be equal to zero for a
//!   closed domain.
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
  for (std::size_t f=m_fd.Nipfac(); f<m_fd.Esuf().size()/2; ++f) {
    s[0] += m_geoFace(f,0,0) * m_geoFace(f,1,0);
    s[1] += m_geoFace(f,0,0) * m_geoFace(f,2,0);
    s[2] += m_geoFace(f,0,0) * m_geoFace(f,3,0);
  }

  auto eps = std::numeric_limits< tk::real >::epsilon() * 100;
  return std::abs(s[0]) > eps || std::abs(s[1]) > eps || std::abs(s[2]) > eps;
}

bool
DG::faceMatch()
// *****************************************************************************
// Check if esuf of chare-boundary faces matches
//! \details This function checks each chare-boundary esuf entry for the left
//!   and right elements. Then, it tries to match all vertices of these
//!   elements. Exactly three of these vertices must match if the esuf entry
//!   has been updated correctly at chare-boundaries.
//! \return True if chare-boundary faces match.
// *****************************************************************************
{
  const auto& esuf = m_fd.Esuf();
  const auto& inpoel = Disc()->Inpoel();
  const auto& coord = Disc()->Coord();
  bool match(true);

  auto eps = std::numeric_limits< tk::real >::epsilon() * 100;

  for (auto f=m_fd.Nipfac(); f<esuf.size()/2; ++f)
  {
    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    std::size_t count = 0;

    for (std::size_t i=0; i<4; ++i)
    {
      auto ip = inpoel[4*el+i];
      for (std::size_t j=0; j<4; ++j)
      {
        auto jp = inpoel[4*er+j];
        auto xdiff = std::abs( coord[0][ip] - coord[0][jp] );
        auto ydiff = std::abs( coord[1][ip] - coord[1][jp] );
        auto zdiff = std::abs( coord[2][ip] - coord[2][jp] );

        if ( xdiff<=eps && ydiff<=eps && zdiff<=eps ) ++count;
      }
    }

    match = (match && count == 3);
  }

  return match;
}

void
DG::comfac( int fromch, const tk::UnsMesh::FaceSet& infaces )
// *****************************************************************************
//  Receive unique set of faces we potentially share with/from another chare
//! \param[in] fromch Sender chare id
//! \param[in] infaces Unique set of faces we potentially share with fromch
// *****************************************************************************
{
  const auto& esuel = m_fd.Esuel();

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
  for (std::size_t e=0; e<esuel.size()/4; ++e) {  // for all our tets
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f) {  // for all cell faces
      if (esuel[mark+f] == -1) {  // if face has no outside-neighbor tet
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
    m_ncomfac = 0;
    if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) d->Tr().chcomfac();
    tk::destroy(m_ipface);

    // Ensure correct number of expected vs received/found chare-boundary faces
    Assert( m_exptNbface == tk::sumvalsize(m_bndFace), 
            "Expected and received number of boundary faces mismatch" );
    m_exptNbface = 0;

    // Basic error checking on chare-boundary-face map
    Assert( m_bndFace.find( thisIndex ) == m_bndFace.cend(),
            "Face-communication map should not contain data for own chare ID" );

    // Store (local) tet ID adjacent to our chare boundary from the inside
    for (std::size_t e=0; e<esuel.size()/4; ++e) {  // for all our tets
      auto mark = e*4;
      for (std::size_t f=0; f<4; ++f) {  // for all cell faces
        if (esuel[mark+f] == -1) {  // if face has no outside-neighbor tet
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
  const auto& coord = d->Coord();

  // Enlarge elements surrounding faces data structure for ghosts
  m_fd.Esuf().resize( 2*m_nfac, -2 );
  m_fd.Inpofa().resize( 3*m_nfac, 0 );
  // Enlarge face geometry data structure for ghosts
  m_geoFace.resize( m_nfac, 0.0 );

  const auto& esuel = m_fd.Esuel();

  // Collect tet ids, their face connectivity (given by 3 global node IDs, each
  // triplet for potentially multiple faces on the chare boundary), and their
  // elem geometry data (see GhostData) associated to fellow chares adjacent to
  // chare boundaries. Once received by fellow chares, these tets will become
  // known as ghost elements and their data as ghost data.
  for (std::size_t e=0; e<esuel.size()/4; ++e) {  // for all our tets
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f) {  // for all cell faces
      if (esuel[mark+f] == -1) {  // if face has no outside-neighbor tet
        tk::UnsMesh::Face t{{ gid[ inpoel[ mark + tk::lpofa[f][0] ] ],
                              gid[ inpoel[ mark + tk::lpofa[f][1] ] ],
                              gid[ inpoel[ mark + tk::lpofa[f][2] ] ] }};
        auto c = findchare( t );
        // It is possible that we do not find the chare for this face. We are
        // looping through all of our tets and interrogating all faces that do
        // not have neighboring tets but we only care about chare-boundary faces
        // here as only those need ghost data. (esuel may also contain
        // physical boundary faces.)
        if (c > -1) {
          // Will store ghost data associated to neighbor chare
          auto& ghost = m_ghostData[ c ];
          // Store tet id adjacent to chare boundary as key for ghost data
          auto& tuple = ghost[ e ];
          // If tetid e has not yet been encountered, store geometry (only once)
          auto& nodes = std::get< 0 >( tuple );
          if (nodes.empty()) {
            std::get< 1 >( tuple ) = m_geoElem[ e ];

            auto& ncoord = std::get< 2 >( tuple );
            ncoord[0] = coord[0][ inpoel[ mark+f ] ];
            ncoord[1] = coord[1][ inpoel[ mark+f ] ];
            ncoord[2] = coord[2][ inpoel[ mark+f ] ];

            std::get< 3 >( tuple ) = f;

            std::get< 4 >( tuple ) = {{ gid[ inpoel[ mark ] ],
                                        gid[ inpoel[ mark+1 ] ],
                                        gid[ inpoel[ mark+2 ] ],
                                        gid[ inpoel[ mark+3 ] ] }};
          }
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
    for ([[maybe_unused]] const auto& t : c.second) {
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
      Assert( !std::get< 2 >( t.second ).empty(),
              "No nodal coordinate data for ghost" );
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
  if (++m_ghostReq == m_msumset.size()) {
    m_ghostReq = 0;
    reqghost_complete();
  }
}

void
DG::sendGhost()
// *****************************************************************************
// Send all of our ghost data to fellow chares
// *****************************************************************************
{
  for (const auto& c : m_ghostData)
    thisProxy[ c.first ].comGhost( thisIndex, c.second );

  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) Disc()->Tr().chghost();
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
    // cppcheck-suppress useStlAlgorithm
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
  auto d = Disc();
  const auto& lid = d->Lid();
  auto& inpofa = m_fd.Inpofa();
  auto& inpoel = d->Inpoel();
  auto& coord = d->Coord();
  auto ncoord = coord[0].size();

  // nodelist with fromch, currently only used for an assert
  [[maybe_unused]] const auto& nl = tk::cref_find( m_msumset, fromch );

  auto& ghostelem = m_ghost[ fromch ];  // will associate to sender chare

  // Store ghost data coming from chare
  for (const auto& g : ghost) {  // loop over incoming ghost data
    auto e = g.first;  // remote/ghost tet id outside of chare boundary
    const auto& nodes = std::get< 0 >( g.second );  // node IDs of face(s)
    const auto& geo = std::get< 1 >( g.second );    // ghost elem geometry data
    const auto& coordg = std::get< 2 >( g.second );  // coordinate of ghost node
    const auto& inpoelg = std::get< 4 >( g.second ); // inpoel of ghost tet

    Assert( nodes.size() % 3 == 0, "Face node IDs must be triplets" );
    Assert( nodes.size() <= 4*3, "Overflow of faces/tet received" );
    Assert( geo.size() % 4 == 0, "Ghost geometry size mismatch" );
    Assert( geo.size() == m_geoElem.nprop(), "Ghost geometry number mismatch" );
    Assert( coordg.size() == 3, "Incorrect ghost node coordinate size" );
    Assert( inpoelg.size() == 4, "Incorrect ghost inpoel size" );

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
      // add node-triplet to node-face connectivity
      inpofa[3*id[0]+0] = tk::cref_find( lid, t[2] );
      inpofa[3*id[0]+1] = tk::cref_find( lid, t[1] );
      inpofa[3*id[0]+2] = tk::cref_find( lid, t[0] );

      // if ghost tet id not yet encountered on boundary with fromch
      auto i = ghostelem.find( e );
      if (i != end(ghostelem)) {
        addEsuf( id, i->second );  // fill in elements surrounding face
        addEsuel( id, i->second, t ); // fill in elements surrounding element
      } else {
        addEsuf( id, m_nunk );     // fill in elements surrounding face
        addEsuel( id, m_nunk, t ); // fill in elements surrounding element
        ghostelem[e] = m_nunk;     // assign new local tet id to remote ghost id
        m_geoElem.push_back( geo );// store ghost elem geometry
        ++m_nunk;                  // increase number of unknowns on this chare
        std::size_t counter = 0;
        for (std::size_t gp=0; gp<4; ++gp) {
          auto it = lid.find( inpoelg[gp] );
          std::size_t lp;
          if (it != end(lid))
            lp = it->second;
          else {
            Assert( nodes.size() == 3, "Expected node not found in lid" );
            Assert( gp == std::get< 3 >( g.second ),
                    "Ghost node not matching correct entry in ghost inpoel" );
            lp = ncoord;
            ++counter;
          }
          inpoel.push_back( lp );       // store ghost element connectivity
        }
        // only a single or no ghost node should be found
        Assert( counter <= 1, "Incorrect number of ghost nodes detected. "
                "Detected "+ std::to_string(counter) +" ghost nodes" );
        if (counter == 1) {
          coord[0].push_back( coordg[0] ); // store ghost node coordinate
          coord[1].push_back( coordg[1] );
          coord[2].push_back( coordg[2] );
          Assert( inpoel[ 4*(m_nunk-1)+std::get< 3 >( g.second ) ] == ncoord,
                  "Mismatch in extended inpoel for ghost element" );
          ++ncoord;                // increase number of nodes on this chare
        }
      }

      // additional tests to ensure that entries in inpoel and t/inpofa match
      Assert( nodetripletMatch(id, t) == 3, "Mismatch/Overmatch in inpoel and "
              "inpofa at chare-boundary face" );
    }
  }

  // Signal the runtime system that all workers have received their adjacency
  if (++m_nadj == m_ghostData.size()) adj();
}

std::size_t
DG::nodetripletMatch( const std::array< std::size_t, 2 >& id,
                      const tk::UnsMesh::Face& t )
// *****************************************************************************
// Check if entries in inpoel, inpofa and node-triplet are consistent
//! \param[in] id Local face and (inner) tet id adjacent to it
//! \param[in] t node-triplet associated with the chare boundary face
//! \return number of nodes in inpoel that matched with t and inpofa
// *****************************************************************************
{
  const auto& lid = Disc()->Lid();
  const auto& inpoel = Disc()->Inpoel();
  const auto& esuf = m_fd.Esuf();
  const auto& inpofa = m_fd.Inpofa();

  std::size_t counter = 0;
  for (std::size_t k=0; k<4; ++k)
  {
    auto el = esuf[ 2*id[0] ];
    auto ip = inpoel[ 4*static_cast< std::size_t >( el )+k ];
    Assert( el == static_cast< int >( id[1] ), "Mismatch in id and esuf" );
    for (std::size_t j=0; j<3; ++j)
    {
      auto jp = tk::cref_find( lid, t[j] );
      auto fp = inpofa[ 3*id[0]+(2-j) ];
      if (ip == jp && ip == fp) ++counter;
    }
  }

  return counter;
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
  Assert( esuf[ 2*id[0] ] == -2 && esuf[ 2*id[0]+1 ] == -2, "Updating esuf at "
          "wrong location instead of chare-boundary" );
  esuf[ 2*id[0]+0 ] = static_cast< int >( id[1] );
  // put in local id for outer/ghost tet
  esuf[ 2*id[0]+1 ] = static_cast< int >( ghostid );
}

void
DG::addEsuel( const std::array< std::size_t, 2 >& id,
              std::size_t ghostid,
              const tk::UnsMesh::Face& t )
// *****************************************************************************
// Fill elements surrounding a element along chare boundary
//! \param[in] id Local face and (inner) tet id adjacent to it
//! \param[in] ghostid Local ID for ghost tet
//! \param[in] t node-triplet associated with the chare boundary face
//! \details This function updates the elements surrounding element (esuel) data
//    structure for the (inner) tets adjacent to the chare-boundaries. It fills
//    esuel of this inner tet with the local tet-id that has been assigned to
//    the outer ghost tet in DG::comGhost in place of the -1 before.
// *****************************************************************************
{
  auto d = Disc();
  const auto& inpoel = d->Inpoel();
  [[maybe_unused]] const auto& esuf = m_fd.Esuf();
  const auto& lid = d->Lid();

  std::array< tk::UnsMesh::Face, 4 > face;
  for (std::size_t f = 0; f<4; ++f)
    for (std::size_t i = 0; i<3; ++i)
      face[f][i] = inpoel[ id[1]*4 + tk::lpofa[f][i] ];

  tk::UnsMesh::Face tl{{ tk::cref_find( lid, t[0] ),
                         tk::cref_find( lid, t[1] ),
                         tk::cref_find( lid, t[2] ) }};

  auto& esuel = m_fd.Esuel();
  std::size_t i(0), nmatch(0);
  for (const auto& f : face) {
    if (tk::UnsMesh::Eq< 3 >()( tl, f )) {
      Assert( esuel[ id[1]*4 + i ] == -1, "Incorrect boundary element found in "
             "esuel");
      esuel[ id[1]*4 + i ] = static_cast<int>(ghostid);
      ++nmatch;
      Assert( esuel[ id[1]*4 + i ] == esuf[ 2*id[0]+1 ], "Incorrect boundary "
             "element entered in esuel" );
      Assert( static_cast<int>(id[1]) == esuf[ 2*id[0]+0 ], "Boundary "
             "element entered in incorrect esuel location" );
    }
    ++i;
  }

  // ensure that exactly one face matched
  Assert( nmatch == 1, "Incorrect number of node-triplets (faces) matched for "
         "updating esuel; matching faces = "+ std::to_string(nmatch) );
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
  const auto& coord = d->Coord();
  const auto& lid = d->Lid();

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
  m_nadj = 0;

  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) Disc()->Tr().chadj();

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

  // Ensure that all elements surrounding elements are correct including those
  // at chare boundaries
  const auto& esuel = m_fd.Esuel();
  std::size_t nbound = 0;
  for (std::size_t e=0; e<esuel.size()/4; ++e) {
    for (std::size_t f=0; f<4; ++f)
      if (esuel[4*e+f] == -1) ++nbound;
  }
  Assert( nbound == m_fd.Nbfac(), "Incorrect number of ghost-element -1's in "
         "updated esuel" );

  // Error checking on ghost data
  for(const auto& n : m_ghostData)
    for([[maybe_unused]] const auto& i : n.second)
      Assert( i.first < m_fd.Esuel().size()/4, "Sender contains ghost tet id " );

  // Perform leak test on face geometry data structure enlarged by ghosts
  Assert( !leakyAdjacency(), "Face adjacency leaky" );
  Assert( faceMatch(), "Chare-boundary element-face connectivity (esuf) does "
         "not match" );

  // Resize solution vectors, lhs, rhs and limiter function by the number of
  // ghost tets
  m_u.resize( m_nunk );
  m_un.resize( m_nunk );
  m_lhs.resize( m_nunk );
  m_rhs.resize( m_nunk );
  m_limFunc.resize( m_nunk );

  // Ensure that we also have all the geometry and connectivity data 
  // (including those of ghosts)
  Assert( m_geoElem.nunk() == m_u.nunk(), "GeoElem unknowns size mismatch" );
  Assert( Disc()->Inpoel().size()/4 == m_u.nunk(), "Inpoel size mismatch" );

  // Basic error checking on ghost tet ID map
  Assert( m_ghost.find( thisIndex ) == m_ghost.cend(),
          "Ghost id map should not contain data for own chare ID" );

  // Store expected ghost tet IDs
  for (const auto& c : m_ghost)
    for ([[maybe_unused]] const auto& g : c.second)
      Assert( m_exptGhost.insert( g.second ).second,
              "Failed to store local tetid as exptected ghost id" );

  // Signal the runtime system that all workers have received their adjacency
  contribute( sizeof(int), &m_initial, CkReduction::sum_int,
    CkCallback(CkReductionTarget(Transporter,comfinal), Disc()->Tr()) );
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
DG::ResumeFromSync()
// *****************************************************************************
//  Return from migration
//! \details This is called when load balancing (LB) completes. The presence of
//!   this function does not affect whether or not we block on LB.
// *****************************************************************************
{
  if (Disc()->It() == 0) Throw( "it = 0 in ResumeFromSync()" );

  if (!g_inputdeck.get< tag::cmd, tag::nonblocking >()) next();
}

void
DG::setup( tk::real )
// *****************************************************************************
// Set initial conditions, generate lhs, output mesh
// *****************************************************************************
{
  tk::destroy(m_msumset);

  auto d = Disc();

  // Basic error checking on sizes of element geometry data and connectivity
  Assert( m_geoElem.nunk() == m_lhs.nunk(), "Size mismatch in DG::setup()" );
  Assert( d->Inpoel().size()/4 == m_lhs.nunk(),
          "Size mismatch in DG::setup()" );

  // Compute left-hand side of discrete PDEs
  lhs();

  // Set initial conditions for all PDEs
  for (const auto& eq : g_dgpde) 
    eq.initialize( m_lhs, d->Inpoel(), d->Coord(), m_u, d->T(),
                   m_fd.Esuel().size()/4 );
  m_un = m_u;

  for (std::size_t e=0; e<m_fd.Esuel().size()/4; ++e)
    for (std::size_t c=0; c<m_limFunc.nprop(); ++c)
      m_limFunc(e,c,0)=1.0;

  // Communicate for initial solution limiting
  contribute( CkCallback(CkReductionTarget(Transporter,sendinit), d->Tr()) );
}

void
DG::limitIC()
// *****************************************************************************
//  Limit initial solution and prepare for time stepping
//! \details This function applies limiter to initial solution and then proceeds
//!   to communicate this limited solution and begin time stepping
// *****************************************************************************
{
  // Limit initial solution
  const auto limiter = g_inputdeck.get< tag::discr, tag::limiter >();
  if (limiter == ctr::LimiterType::WENOP1)
  {
    WENO_P1( m_fd.Esuel(), 0, m_u, m_limFunc );

    const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
    const auto ncomp= m_u.nprop()/ndof;
    for (inciter::ncomp_t c=0; c<ncomp; ++c)
    {
      auto mark = c*ndof;
      auto lmark = c*(ndof-1);
      for (std::size_t e=0; e<m_u.nunk(); ++e)
      {
        // limit P1 dofs
        m_u( e, mark+1, 0 ) = m_limFunc( e, lmark  , 0 ) * m_u( e, mark+1, 0 );
        m_u( e, mark+2, 0 ) = m_limFunc( e, lmark+1, 0 ) * m_u( e, mark+2, 0 );
        m_u( e, mark+3, 0 ) = m_limFunc( e, lmark+2, 0 ) * m_u( e, mark+3, 0 );
      }
    }
  }

  // Output initial conditions to file (regardless of whether it was requested)
  writeFields( CkCallback(CkIndex_DG::start(), thisProxy[thisIndex]) );
}

void
DG::start()
// *****************************************************************************
//  Start time stepping
// *****************************************************************************
{
  auto d = Disc();

  // Start timer measuring time stepping wall clock time
  d->Timer().zero();

  tk::real fdt = 0.0;
  // Start time stepping
  contribute( sizeof(tk::real), &fdt, CkReduction::nop,
              CkCallback(CkReductionTarget(Transporter,advance), d->Tr()) );
}

void
DG::sendinit()
// *****************************************************************************
// Send own chare-boundary data to neighboring chares
// *****************************************************************************
{
  // communicate solution ghost data (if any)
  if (m_ghostData.empty())
    cominit_complete();
  else
    for(const auto& n : m_ghostData) {
      std::vector< std::size_t > tetid;
      std::vector< std::vector< tk::real > > u;
      for(const auto& i : n.second) {
        Assert( i.first < m_fd.Esuel().size()/4, "Sending solution ghost data" );
        tetid.push_back( i.first );
        u.push_back( m_u[i.first] );
      }
      thisProxy[ n.first ].cominit( thisIndex, tetid, u );
    }

  owninit_complete();
}

void
DG::cominit( int fromch,
             const std::vector< std::size_t >& tetid,
             const std::vector< std::vector< tk::real > >& u )
// *****************************************************************************
//  Receive chare-boundary solution ghost data from neighboring chares
//! \param[in] fromch Sender chare id
//! \param[in] tetid Ghost tet ids we receive solution data for
//! \param[in] u Solution ghost data
//! \details This function receives contributions to m_u from fellow chares.
// *****************************************************************************
{
  Assert( u.size() == tetid.size(), "Size mismatch in DG::cominit()" );

  // Find local-to-ghost tet id map for sender chare
  const auto& n = tk::cref_find( m_ghost, fromch );

  for (std::size_t i=0; i<tetid.size(); ++i) {
    auto j = tk::cref_find( n, tetid[i] );
    Assert( j >= m_fd.Esuel().size()/4, "Receiving solution non-ghost data" );
    Assert( j < m_u.nunk(), "Indexing out of bounds in DG::cominit()" );
    Assert( m_recvGhost.insert( j ).second,
            "Failed to store local tetid of received ghost tetid" );
    for (std::size_t c=0; c<m_u.nprop(); ++c)
      m_u(j,c,0) = u[i][c];
  }

  // if we have received all solution ghost contributions from those chares we
  // communicate along chare-boundary faces with, solve the system
  if (++m_ninitsol == m_ghostData.size()) {
    Assert( m_exptGhost == m_recvGhost,
            "Expected/received ghost tet id mismatch" );
    m_exptGhost.clear();
    m_recvGhost.clear();
    m_ninitsol = 0;
    cominit_complete();
  }
}

void
DG::dt()
// *****************************************************************************
// Compute time step size
// *****************************************************************************
{
  auto mindt = std::numeric_limits< tk::real >::max();

  auto d = Disc();

  if (m_stage == 0)
  {
    auto const_dt = g_inputdeck.get< tag::discr, tag::dt >();
    auto def_const_dt = g_inputdeck_defaults.get< tag::discr, tag::dt >();
    auto eps = std::numeric_limits< tk::real >::epsilon();

    // use constant dt if configured
    if (std::abs(const_dt - def_const_dt) > eps) {

      mindt = const_dt;

    } else {      // compute dt based on CFL

      // find the minimum dt across all PDEs integrated
      for (const auto& eq : g_dgpde) {
        auto eqdt = eq.dt( d->Coord(), d->Inpoel(), m_fd, m_geoFace, m_geoElem,
                           m_limFunc, m_u );
        if (eqdt < mindt) mindt = eqdt;
      }

      // Scale smallest dt with CFL coefficient
      mindt *= g_inputdeck.get< tag::discr, tag::cfl >();

    }
  }
  else
  {
    mindt = d->Dt();
  }

  // Contribute to minimum dt across all chares then advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(DG,solve), thisProxy) );
}

void
DG::advance( tk::real )
// *****************************************************************************
// Advance equations to next time step
// *****************************************************************************
{
  // communicate solution ghost data (if any)
  if (m_ghostData.empty())
    comsol_complete();
  else
    for(const auto& n : m_ghostData) {
      std::vector< std::size_t > tetid;
      std::vector< std::vector< tk::real > > u;
      for(const auto& i : n.second) {
        Assert( i.first < m_fd.Esuel().size()/4, "Sending solution ghost data" );
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
//! \param[in] fromch Sender chare id
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
    Assert( j >= m_fd.Esuel().size()/4, "Receiving solution non-ghost data" );
    Assert( j < m_u.nunk(), "Indexing out of bounds in DG::comsol()" );
    for (std::size_t c=0; c<m_u.nprop(); ++c)
      m_u(j,c,0) = u[i][c];
  }

  // if we have received all solution ghost contributions from those chares we
  // communicate along chare-boundary faces with, solve the system
  if (++m_nsol == m_ghostData.size()) {
    m_nsol = 0;
    comsol_complete();
  }
}

void
DG::writeFields( CkCallback c )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  auto d = Disc();

  const auto& esuel = m_fd.Esuel();

  // Copy mesh form Discretization object and chop off ghosts for dump
  auto inpoel = d->Inpoel();
  inpoel.resize( esuel.size() );
  auto coord = d->Coord();
  for (std::size_t i=0; i<3; ++i) coord[i].resize( m_ncoord );

  // Query fields names from all PDEs integrated
  std::vector< std::string > names;
  for (const auto& eq : g_dgpde) {
    auto n = eq.fieldNames();
    names.insert( end(names), begin(n), end(n) );
  }

  // Collect element field solution
  std::vector< std::vector< tk::real > > fields;
  auto u = m_u;
  for (const auto& eq : g_dgpde) {
    auto o =
      eq.fieldOutput( d->T(), m_geoElem, u );
    // cut off ghost elements
    for (auto& field : o) field.resize( esuel.size()/4 );
    fields.insert( end(fields), begin(o), end(o) );
  }

  // // Collect node field solution
  // std::vector< std::vector< tk::real > > nodefields;
  // for (const auto& eq : g_dgpde) {
  //   auto fields =
  //     eq.avgElemToNode( d->Inpoel(), d->Coord(), m_geoElem, m_limFunc, m_u );
  //   nodefields.insert( end(nodefields), begin(fields), end(fields) );
  // }

  // Output chare mesh and fields metadata to file
  d->write( inpoel, coord, m_fd.Bface(), {}, m_fd.Triinpoel(), names,
            fields, tk::Centering::ELEM, c );
}

void
DG::lhs()
// *****************************************************************************
// Compute left-hand side of discrete transport equations
// *****************************************************************************
{
  for (const auto& eq : g_dgpde) eq.lhs( m_geoElem, m_lhs );

  if (!m_initial) stage();
}

void
DG::lim()
// *****************************************************************************
// Compute limiter function
// *****************************************************************************
{
  if (g_inputdeck.get< tag::discr, tag::ndof >() > 1) {
  
    Assert( m_u.nunk() == m_limFunc.nunk(), "Number of unknowns in solution "
            "vector and limiter at recent time step incorrect" );
  
    Assert( m_u.nprop() == m_limFunc.nprop() + 
            (m_u.nprop()/g_inputdeck.get< tag::discr, tag::ndof >()),
            "Number of components in solution vector and limiter at recent "
            "time step incorrect" );

    for (std::size_t e=0; e<m_fd.Esuel().size()/4; ++e)
      for (std::size_t c=0; c<m_limFunc.nprop(); ++c)
        m_limFunc(e,c,0)=1.0;
  
    const auto limiter = g_inputdeck.get< tag::discr, tag::limiter >();
    if (limiter == ctr::LimiterType::WENOP1)
      WENO_P1( m_fd.Esuel(), 0, m_u, m_limFunc );
  
    // communicate solution ghost data (if any)
    if (m_ghostData.empty())
      comlim_complete();
    else
      for(const auto& n : m_ghostData) {
        std::vector< std::size_t > tetid;
        std::vector< std::vector< tk::real > > limfunc;
        for(const auto& i : n.second) {
          Assert( i.first < m_fd.Esuel().size()/4, "Sending limiter ghost data" );
          tetid.push_back( i.first );
          limfunc.push_back( m_limFunc[i.first] );
        }
        thisProxy[ n.first ].comlim( thisIndex, tetid, limfunc );
      }

  } else {

    comlim_complete();

  }

  ownlim_complete();
}

void
DG::comlim( int fromch,
            const std::vector< std::size_t >& tetid,
            const std::vector< std::vector< tk::real > >& lfn )
// *****************************************************************************
//  Receive chare-boundary limiter ghost data from neighboring chares
//! \param[in] fromch Sender chare id
//! \param[in] tetid Ghost tet ids we receive solution data for
//! \param[in] lfn Limiter function ghost data
//! \details This function receives contributions to m_limFunc from fellow
//    chares.
// *****************************************************************************
{
  Assert( lfn.size() == tetid.size(), "Size mismatch in DG::comlim()" );

  // Find local-to-ghost tet id map for sender chare
  const auto& n = tk::cref_find( m_ghost, fromch );

  for (std::size_t i=0; i<tetid.size(); ++i) {
    auto j = tk::cref_find( n, tetid[i] );
    Assert( j >= m_fd.Esuel().size()/4, "Receiving solution non-ghost data" );
    Assert( j < m_limFunc.nunk(), "Indexing out of bounds in DG::comlim()" );
    Assert( lfn[i].size() == m_limFunc.nprop(), "size mismatch in received "
           "limfunc in DG::comlim()" );
    for (std::size_t c=0; c<m_limFunc.nprop(); ++c)
      m_limFunc(j,c,0) = lfn[i][c];
  }

  // if we have received all solution ghost contributions from those chares we
  // communicate along chare-boundary faces with, solve the system
  if (++m_nlim == m_ghostData.size()) {
    m_nlim = 0;
    comlim_complete();
  }
}

void
DG::solve( tk::real newdt )
// *****************************************************************************
// Compute right-hand side of discrete transport equations
//! \param[in] newdt Size of this new time step
// *****************************************************************************
{
  // Enable SDAG wait for building the solution vector during the next stage
  thisProxy[ thisIndex ].wait4sol();
  thisProxy[ thisIndex ].wait4lim();

  auto d = Disc();

  // Set new time step size
  d->setdt( newdt );

  for (const auto& eq : g_dgpde)
    eq.rhs( d->T(), m_geoFace, m_geoElem, m_fd, d->Inpoel(), d->Coord(), m_u,
            m_limFunc, m_rhs );

  // Explicit time-stepping using RK3 to discretize time-derivative
  m_u =  rkcoef[0][m_stage] * m_un
       + rkcoef[1][m_stage] * ( m_u + d->Dt() * m_rhs/m_lhs );

  if (m_stage < 2) {

    // continue with next tims step stage
    stage();

  } else {

    thisProxy[ thisIndex ].wait4recompghost();

    // Compute diagnostics, e.g., residuals
    auto diag_computed =
      m_diag.compute( *d, m_u.nunk()-m_fd.Esuel().size()/4, m_geoElem, m_u );
    // Increase number of iterations and physical time
    d->next();
    // Update Un
    m_un = m_u;
    // Signal that diagnostics have been computed (or in this case, skipped)
    if (!diag_computed) diag();
    // Optionally refine mesh
    refine();

  }
}

void
DG::resized()
// *****************************************************************************
// Resizing data sutrctures after mesh refinement has been completed
// *****************************************************************************
{
  resize_complete();
}

void
DG::diag()
// *****************************************************************************
// Signal the runtime system that diagnostics have been computed
// *****************************************************************************
{
  diag_complete();
}

void
DG::refine()
// *****************************************************************************
// Optionally refine/derefine mesh
// *****************************************************************************
{
  auto d = Disc();

  auto dtref = g_inputdeck.get< tag::amr, tag::dtref >();
  auto dtfreq = g_inputdeck.get< tag::amr, tag::dtfreq >();

  // if t>0 refinement enabled and we hit the frequency
  if (dtref && !(d->It() % dtfreq)) {   // refine

    d->Ref()->dtref( m_fd.Bface(), {}, tk::remap(m_fd.Triinpoel(),d->Gid()) );
    m_refined = 1;

  } else {      // do not refine

    ref_complete();
    resize_complete();
    m_refined = 0;

  }
}

void
DG::resizeAfterRefined(
  const std::vector< std::size_t >& /*ginpoel*/,
  const tk::UnsMesh::Chunk& chunk,
  const tk::UnsMesh::Coords& coord,
  const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& /*addedNodes*/,
  const std::unordered_map< std::size_t, std::size_t >& addedTets,
  const std::unordered_map< int, std::vector< std::size_t > >& msum,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& /* bnode */,
  const std::vector< std::size_t >& triinpoel )
// *****************************************************************************
//  Receive new mesh from refiner
//! \param[in] ginpoel Mesh connectivity with global node ids
//! \param[in] chunk New mesh chunk (connectivity and global<->local id maps)
//! \param[in] coord New mesh node coordinates
//! \param[in] addedNodes Newly added mesh nodes and their parents (local ids)
//! \param[in] addedTets Newly added mesh cells and their parents (local ids)
//! \param[in] msum New node communication map
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] bnode Boundary-node lists mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
// *****************************************************************************
{
  auto d = Disc();

  // Set flag that indicates that we are during time stepping
  m_initial = 0;

  // Zero field output iteration count between two mesh refinement steps
  d->Itf() = 0;

  // Increase number of iterations with mesh refinement
  ++d->Itr();

  // Save old number of elements
  [[maybe_unused]] auto old_nelem = d->Inpoel().size()/4;

  // Resize mesh data structures
  d->resize( chunk, coord, msum );

  // Update state
  auto nelem = d->Inpoel().size()/4;
  auto nprop = m_u.nprop();
  m_u.resize( nelem, nprop );
  m_un.resize( nelem, nprop );
  m_lhs.resize( nelem, nprop );
  m_rhs.resize( nelem, nprop );


  m_fd = FaceData( d->Inpoel(), bface, tk::remap(triinpoel,d->Lid()) );

  m_geoFace =
    tk::Fields( tk::genGeoFaceTri( m_fd.Nipfac(), m_fd.Inpofa(), coord ) );
  m_geoElem = tk::Fields( tk::genGeoElemTet( d->Inpoel(), coord ) );

  m_limFunc = limFunc( nelem );

  m_nfac = m_fd.Inpofa().size()/3;
  m_nunk = nelem;
  m_ncoord = coord[0].size();
  m_msumset = d->msumset();
  m_bndFace.clear();
  m_ghostData.clear();
  m_ghost.clear();

  // Update solution on new mesh, P0 (cell center value) only for now
  m_un = m_u;
  for (const auto& e : addedTets) {
    Assert( e.first < nelem, "Indexing out of new solution vector" );
    Assert( e.second < old_nelem, "Indexing out of old solution vector" );
    for (std::size_t c=0; c<nprop; ++c)
      m_u(e.first,c,0) = m_un(e.second,c,0);
  }
  m_un = m_u;

  ref_complete();

  contribute( CkCallback(CkReductionTarget(Transporter,workresized), d->Tr()) );
}

void
DG::recompGhostRefined()
// *****************************************************************************
// Start recomputing ghost data after a mesh refinement step
// *****************************************************************************
{
  if (m_refined) resizeComm(); else stage();
}

void
DG::next()
// *****************************************************************************
// Continue to next time step stage
// *****************************************************************************
{
  auto d = Disc();

  tk::real fdt = 0.0;
  d->contribute( sizeof(tk::real), &fdt, CkReduction::nop,
                 CkCallback(CkReductionTarget(Transporter,advance), d->Tr()) );
}

void
DG::out()
// *****************************************************************************
// Output mesh field data
// *****************************************************************************
{
  auto d = Disc();

  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto fieldfreq = g_inputdeck.get< tag::interval, tag::field >();

  // output field data if field iteration count is reached or in the last time
  // step, otherwise continue to next time step
  if ( !((d->It()) % fieldfreq) ||
       (std::fabs(d->T()-term) < eps || d->It() >= nstep) )
    writeFields( CkCallback(CkIndex_DG::step(), thisProxy[thisIndex]) );
  else
    step();
}

void
DG::stage()
// *****************************************************************************
// Evaluate whether to continue with next time step stage
// *****************************************************************************
{
  // Increment Runge-Kutta stage counter
  ++m_stage;

  // if not all Runge-Kutta stages complete, continue to next time stage,
  // otherwise output field data to file(s)
  if (m_stage < 3) next(); else out();
}

void
DG::step()
// *****************************************************************************
// Evaluate wether to continue with next time step
// *****************************************************************************
{
  auto d = Disc();

  // Output one-liner status report to screen
  d->status();
  // Reset Runge-Kutta stage counter
  m_stage = 0;

  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto lbfreq = g_inputdeck.get< tag::cmd, tag::lbfreq >();
  const auto nonblocking = g_inputdeck.get< tag::cmd, tag::nonblocking >();

  // If neither max iterations nor max time reached, continue, otherwise finish
  if (std::fabs(d->T()-term) > eps && d->It() < nstep) {

    if ( (d->It()) % lbfreq == 0 ) {
      AtSync();
      if (nonblocking) next();
    }
    else {
      next();
    }

  } else {
    contribute(CkCallback( CkReductionTarget(Transporter,finish), d->Tr() ));
  }
}

#include "NoWarning/dg.def.h"
