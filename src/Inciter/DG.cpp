// *****************************************************************************
/*!
  \file      src/Inciter/DG.cpp
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
#include <sstream>

#include "DG.hpp"
#include "Discretization.hpp"
#include "DGPDE.hpp"
#include "DiagReducer.hpp"
#include "DerivedData.hpp"
#include "ElemDiagnostics.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Refiner.hpp"
#include "Limiter.hpp"
#include "PrefIndicator.hpp"
#include "Reorder.hpp"
#include "Vector.hpp"

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
       g_inputdeck.get< tag::discr, tag::rdof >()*
       g_inputdeck.get< tag::component >().nprop() ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_p( m_u.nunk(),
       g_inputdeck.get< tag::discr, tag::rdof >()*
         std::accumulate( begin(g_dgpde), end(g_dgpde), 0u,
           [](std::size_t s, const DGPDE& eq){ return s + eq.nprim(); } ) ),
  m_geoFace( tk::genGeoFaceTri( m_fd.Nipfac(), m_fd.Inpofa(), Disc()->Coord()) ),
  m_geoElem( tk::genGeoElemTet( Disc()->Inpoel(), Disc()->Coord() ) ),
  m_lhs( m_u.nunk(),
         g_inputdeck.get< tag::discr, tag::ndof >()*
         g_inputdeck.get< tag::component >().nprop() ),
  m_rhs( m_u.nunk(), m_lhs.nprop() ),
  m_nfac( m_fd.Inpofa().size()/3 ),
  m_nunk( m_u.nunk() ),
  m_ncoord( Disc()->Coord()[0].size() ),
  m_msumset( Disc()->msumset() ),
  m_bndFace(),
  m_ghostData(),
  m_ghostReq( 0 ),
  m_ghost(),
  m_exptGhost(),
  m_recvGhost(),
  m_diag(),
  m_stage( 0 ),
  m_ndof(),
  m_bid(),
  m_uc(),
  m_pc(),
  m_ndofc(),
  m_initial( 1 ),
  m_expChBndFace()
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
// *****************************************************************************
{
  usesAtSync = true;    // enable migration at AtSync

  // Enable SDAG wait for setting up chare boundary faces
  thisProxy[ thisIndex ].wait4fac();

  // Ensure that mesh partition is not leaky
  Assert( !tk::leakyPartition(m_fd.Esuel(), Disc()->Inpoel(), Disc()->Coord()),
          "Input mesh to DG leaky" );

  // Ensure mesh physical boundary for the entire problem not leaky,
  // effectively checking if the user has specified boundary conditions on all
  // physical boundary faces
  bndIntegral();
}

void
DG::bndIntegral()
// *****************************************************************************
//  Compute partial boundary surface integral and sum across all chares
//! \details This function computes a partial surface integral over the boundary
//!   of the faces of this mesh partition then sends its contribution to perform
//!   the integral acorss the total problem boundary. After the global sum a
//!   non-zero vector result indicates a leak, e.g., a hole in the boundary
//!   which indicates an error in the boundary face data structures used to
//!   compute the partial surface integrals.
// *****************************************************************************
{
  // Storage for surface integral over our mesh chunk physical boundary
  std::vector< tk::real > s{{ 0.0, 0.0, 0.0 }};

  // Integrate over all physical boundary faces
  for (std::size_t f=0; f<m_fd.Nbfac(); ++f) {
    s[0] += m_geoFace(f,0,0) * m_geoFace(f,1,0);
    s[1] += m_geoFace(f,0,0) * m_geoFace(f,2,0);
    s[2] += m_geoFace(f,0,0) * m_geoFace(f,3,0);
  }

  s.push_back( 1.0 );  // positive: call-back to resizeComm() after reduction

  // Send contribution to host summing partial surface integrals
  contribute( s, CkReduction::sum_double,
    CkCallback(CkReductionTarget(Transporter,bndint), Disc()->Tr()) );
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
        // if does not exist among the internal and physical boundary faces,
        // store as a potential chare-boundary face
        tk::UnsMesh::Face t{{ gid[ inpoel[ mark + tk::lpofa[f][0] ] ],
                              gid[ inpoel[ mark + tk::lpofa[f][1] ] ],
                              gid[ inpoel[ mark + tk::lpofa[f][2] ] ] }};
        if (m_ipface.find(t) == end(m_ipface)) {
          Assert( m_expChBndFace.insert(t).second,
                  "Store expected chare-boundary face" );
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

  ownfac_complete();
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
  std::array< tk::real, 3 > s{{ 0.0, 0.0, 0.0 }};

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
  // Buffer up incoming data
  m_infaces[ fromch ] = infaces;

  // if we have heard from all fellow chares that we share at least a single
  // node, edge, or face with
  if (++m_ncomfac == m_msumset.size()) {
    m_ncomfac = 0;
    comfac_complete();
  }
}

void
DG::bndFaces()
// *****************************************************************************
// Compute chare-boundary faces
//! \details This is called when both send and receives are completed on a
//!  chare and thus we are ready to compute chare-boundary faces and ghost data.
// *****************************************************************************
{
  auto d = Disc();
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) d->Tr().chcomfac();
  const auto& esuel = m_fd.Esuel();
  const auto& gid = d->Gid();
  const auto& inpoel = d->Inpoel();

  for (const auto& in : m_infaces) {
    // Find sender chare among chares we potentially share faces with. Note that
    // it is feasible that a sender chare called us but we do not have a set of
    // faces associated to that chare. This can happen if we only share a single
    // node or an edge but not a face with that chare.
    auto& bndface = m_bndFace[ in.first ];  // will associate to sender chare
    // Try to find incoming faces on our chare boundary with other chares. If
    // found, generate and assign new local face ID, associated to sender chare.
    for (std::size_t e=0; e<esuel.size()/4; ++e) {  // for all our tets
      auto mark = e*4;
      for (std::size_t f=0; f<4; ++f) {  // for all cell faces
        if (esuel[mark+f] == -1) {  // if face has no outside-neighbor tet
          tk::UnsMesh::Face t{{ gid[ inpoel[ mark + tk::lpofa[f][0] ] ],
                                gid[ inpoel[ mark + tk::lpofa[f][1] ] ],
                                gid[ inpoel[ mark + tk::lpofa[f][2] ] ] }};
          // if found among the incoming faces and if not one of our internal
          // nor physical boundary faces
          if ( in.second.find(t) != end(in.second) &&
               m_ipface.find(t) == end(m_ipface) ) {
            bndface[t][0] = m_nfac++;    // assign new local face ID
          }
        }
      }
    }
    // If at this point if we have not found any face among our faces we
    // potentially share with fromch, there is no need to keep an empty set of
    // faces associated to fromch as we only share nodes or edges with it, but
    // not faces.
    if (bndface.empty()) m_bndFace.erase( in.first );
  }

  tk::destroy(m_ipface);
  tk::destroy(m_infaces);

  // Ensure all expected faces have been received
  Assert( receivedChBndFaces(),
          "Expected and received chare boundary faces mismatch" );

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

bool
DG::receivedChBndFaces()
// *****************************************************************************
// Verify that all chare-boundary faces have been received
//! \return True if all chare-boundary faces have been received
// *****************************************************************************
{
  auto d = Disc();
  tk::UnsMesh::FaceSet recvBndFace;

  // Collect chare-boundary faces that have been received and expected
  for (const auto& c : m_bndFace)
    for (const auto& f : c.second)
      if (m_expChBndFace.find(f.first) != end(m_expChBndFace))
        recvBndFace.insert(f.first);

   // Collect info on expected but not received faces
   std::stringstream msg;
   for (const auto& f : m_expChBndFace)
     if (recvBndFace.find(f) == end(recvBndFace)) {
       const auto& coord = d->Coord();
       const auto& x = coord[0];
       const auto& y = coord[1];
       const auto& z = coord[2];
       auto A = tk::cref_find( d->Lid(), f[0] );
       auto B = tk::cref_find( d->Lid(), f[1] );
       auto C = tk::cref_find( d->Lid(), f[2] );
       msg << '{' << A << ',' << B << ',' << C << "}:("
           << x[A] << ',' << y[A] << ',' << z[A] << ' '
           << x[B] << ',' << y[B] << ',' << z[B] << ' '
           << x[C] << ',' << y[C] << ',' << z[C] << ") ";
     }

  tk::destroy( m_expChBndFace );

  // Error out with info on missing faces
  auto s = msg.str();
  if (!s.empty()) {
    Throw( "DG chare " + std::to_string(thisIndex) +
           " missing face(s) {local node ids} (node coords): " + s );
  } else {
    return true;
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
        // physical boundary faces)
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

  // Resize solution vectors, lhs and rhs by the number of ghost tets
  m_u.resize( m_nunk );
  m_un.resize( m_nunk );
  m_p.resize( m_nunk );
  m_lhs.resize( m_nunk );
  m_rhs.resize( m_nunk );

  // Create a mapping between local ghost tet ids and zero-based boundary ids
  std::vector< std::size_t > c( tk::sumvalsize( m_ghost ) );
  std::size_t j = 0;
  for (const auto& n : m_ghost) {
    for(const auto& i : n.second) {
      c[j++] = i.second;
    }
  }
  m_bid = tk::assignLid( c );

  // Size communication buffer that receives number of degrees of freedom
  for (auto& n : m_ndofc) n.resize( m_bid.size() );
  for (auto& u : m_uc) u.resize( m_bid.size() );
  for (auto& p : m_pc) p.resize( m_bid.size() );

  // Initialize number of degrees of freedom in mesh elements
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();
  if( pref )
  {
    const auto ndofmax = g_inputdeck.get< tag::pref, tag::ndofmax >();
    m_ndof.resize( m_nunk, ndofmax );
  }
  else
  {
    const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
    m_ndof.resize( m_nunk, ndof );
  }

  // Ensure that we also have all the geometry and connectivity data 
  // (including those of ghosts)
  Assert( m_geoElem.nunk() == m_u.nunk(), "GeoElem unknowns size mismatch" );
  Assert( Disc()->Inpoel().size()/4 == m_u.nunk(), "Inpoel size mismatch" );

  // Basic error checking on ghost tet ID map
  Assert( m_ghost.find( thisIndex ) == m_ghost.cend(),
          "Ghost id map should not contain data for own chare ID" );

  // Store expected ghost tet IDs
  for (const auto& n : m_ghost)
    for ([[maybe_unused]] const auto& g : n.second)
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
//! \details Since this is a [initnode] routine, the runtime system executes the
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
DG::setup()
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
  {
    eq.initialize( m_lhs, d->Inpoel(), d->Coord(), m_u, d->T(),
                   m_fd.Esuel().size()/4 );
    eq.updatePrimitives( m_u, m_p, m_fd.Esuel().size()/4 );
  }

  m_un = m_u;

  // Start timer measuring time stepping wall clock time
  d->Timer().zero();

  // Output initial conditions to file (regardless of whether it was requested)
  writeFields( CkCallback(CkIndex_DG::next(), thisProxy[thisIndex]) );
}

void
DG::next()
// *****************************************************************************
// Advance equations to next time step
// *****************************************************************************
{
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

  auto d = Disc();

  if (pref && m_stage == 0 && d->T() > 0)
    eval_ndof( m_nunk, Disc()->Coord(), Disc()->Inpoel(), m_fd, m_u,
               g_inputdeck.get< tag::pref, tag::indicator >(),
               g_inputdeck.get< tag::discr, tag::ndof >(),
               g_inputdeck.get< tag::pref, tag::ndofmax >(),
               g_inputdeck.get< tag::pref, tag::tolref >(),
               m_ndof );

  // communicate solution ghost data (if any)
  if (m_ghostData.empty())
    comsol_complete();
  else
    for(const auto& n : m_ghostData) {
      std::vector< std::size_t > tetid( n.second.size() );
      std::vector< std::vector< tk::real > > u( n.second.size() ),
                                             prim( n.second.size() );
      std::vector< std::size_t > ndof;
      std::size_t j = 0;
      for(const auto& i : n.second) {
        Assert( i.first < m_fd.Esuel().size()/4, "Sending solution ghost data" );
        tetid[j] = i.first;
        u[j] = m_u[i.first];
        prim[j] = m_p[i.first];
        if (pref && m_stage == 0) ndof.push_back( m_ndof[i.first] );
        ++j;
      }
      thisProxy[ n.first ].comsol( thisIndex, m_stage, tetid, u, prim, ndof );
    }

  ownsol_complete();
}

void
DG::comsol( int fromch,
            std::size_t fromstage,
            const std::vector< std::size_t >& tetid,
            const std::vector< std::vector< tk::real > >& u,
            const std::vector< std::vector< tk::real > >& prim,
            const std::vector< std::size_t >& ndof )
// *****************************************************************************
//  Receive chare-boundary solution ghost data from neighboring chares
//! \param[in] fromch Sender chare id
//! \param[in] fromstage Sender chare time step stage
//! \param[in] tetid Ghost tet ids we receive solution data for
//! \param[in] u Solution ghost data
//! \param[in] prim Primitive variables in ghost cells
//! \param[in] ndof Number of degrees of freedom for chare-boundary elements
//! \details This function receives contributions to the unlimited solution
//!   from fellow chares.
// *****************************************************************************
{
  Assert( u.size() == tetid.size(), "Size mismatch in DG::comsol()" );
  Assert( prim.size() == tetid.size(), "Size mismatch in DG::comsol()" );

  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

  if (pref && fromstage == 0)
    Assert( ndof.size() == tetid.size(), "Size mismatch in DG::comsol()" );

  // Find local-to-ghost tet id map for sender chare
  const auto& n = tk::cref_find( m_ghost, fromch );

  for (std::size_t i=0; i<tetid.size(); ++i) {
    auto j = tk::cref_find( n, tetid[i] );
    Assert( j >= m_fd.Esuel().size()/4, "Receiving solution non-ghost data" );
    auto b = tk::cref_find( m_bid, j );
    Assert( b < m_uc[0].size(), "Indexing out of bounds" );
    m_uc[0][b] = u[i];
    m_pc[0][b] = prim[i];
    if (pref && fromstage == 0) {
      Assert( b < m_ndofc[0].size(), "Indexing out of bounds" );
      m_ndofc[0][b] = ndof[i];
    }
  }

  // if we have received all solution ghost contributions from those chares we
  // communicate along chare-boundary faces with, solve the system
  if (++m_nsol == m_ghostData.size()) {
    m_nsol = 0;
    comsol_complete();
  }
}

void
DG::writeFields( CkCallback c ) const
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
  std::vector< std::string > elemfieldnames;
  for (const auto& eq : g_dgpde) {
    auto n = eq.fieldNames();
    elemfieldnames.insert( end(elemfieldnames), begin(n), end(n) );
  }

  // Collect element field solution
  std::vector< std::vector< tk::real > > elemfields;
  auto u = m_u;
  for (const auto& eq : g_dgpde) {
    auto o = eq.fieldOutput( d->T(), m_geoElem, u );

    // cut off ghost elements
    for (auto& f : o) f.resize( esuel.size()/4 );
    elemfields.insert( end(elemfields), begin(o), end(o) );
  }

  // Add adaptive indicator array to element-centered field output
  std::vector<tk::real> ndof( begin(m_ndof), end(m_ndof) );
  ndof.resize( esuel.size()/4 );  // cut off ghosts
  elemfields.push_back( ndof );

  // // Collect node field solution
  // std::vector< std::vector< tk::real > > nodefields;
  // for (const auto& eq : g_dgpde) {
  //   auto fields =
  //     eq.avgElemToNode( d->Inpoel(), d->Coord(), m_geoElem, m_u );
  //   nodefields.insert( end(nodefields), begin(fields), end(fields) );
  // }

  // Output chare mesh and fields metadata to file
  d->write( inpoel, coord, m_fd.Bface(), {}, m_fd.Triinpoel(), elemfieldnames,
            {}, elemfields, {}, c );
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
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();

  // Combine own and communicated contributions of unlimited solution and
  // degrees of freedom in cells (if p-adaptive)
  for (const auto& b : m_bid) {
    Assert( m_uc[0][b.second].size() == m_u.nprop(), "ncomp size mismatch" );
    Assert( m_pc[0][b.second].size() == m_p.nprop(), "ncomp size mismatch" );
    for (std::size_t c=0; c<m_u.nprop(); ++c) {
      m_u(b.first,c,0) = m_uc[0][b.second][c];
    }
    for (std::size_t c=0; c<m_p.nprop(); ++c) {
      m_p(b.first,c,0) = m_pc[0][b.second][c];
    }
    if (pref && m_stage == 0) {
      m_ndof[ b.first ] = m_ndofc[0][ b.second ];
    }
  }

  if (pref && m_stage==0) propagate_ndof();

  if (rdof > 1) {

    auto d = Disc();

    // Reconstruct second-order solution and primitive quantities
    // if P0P1
    if (rdof == 4 && g_inputdeck.get< tag::discr, tag::ndof >() == 1)
      for (const auto& eq : g_dgpde)
        eq.reconstruct( d->T(), m_geoFace, m_geoElem, m_fd, d->Inpoel(),
                        d->Coord(), m_u, m_p );

    for (const auto& eq : g_dgpde)
      eq.limit( d->T(), m_geoFace, m_geoElem, m_fd, d->Inpoel(), d->Coord(),
                m_ndof, m_u, m_p );
  }

  // Send limited solution to neighboring chares
  if (m_ghostData.empty())
    comlim_complete();
  else
    for(const auto& n : m_ghostData) {
      std::vector< std::size_t > tetid( n.second.size() );
      std::vector< std::vector< tk::real > > u( n.second.size() ),
                                             prim( n.second.size() );
      std::vector< std::size_t > ndof;
      std::size_t j = 0;
      for(const auto& i : n.second) {
        Assert( i.first < m_fd.Esuel().size()/4, "Sending limiter ghost data" );
        tetid[j] = i.first;
        u[j] = m_u[i.first];
        prim[j] = m_p[i.first];
        if (pref && m_stage == 0) ndof.push_back( m_ndof[i.first] );
        ++j;
      }
      thisProxy[ n.first ].comlim( thisIndex, tetid, u, prim, ndof );
    }

  ownlim_complete();
}

void
DG::propagate_ndof()
// *****************************************************************************
//  p-refine all elements that are adjacent to p-refined elements
//! \details This function p-refines all the neighbors of an element that has
//!   been p-refined as a result of an error indicator.
// *****************************************************************************
{
  const auto& esuf = m_fd.Esuf();

  // Copy number of degrees of freedom for each cell
  auto ndof = m_ndof;

  // p-refine all neighboring elements of elements that have been p-refined as a
  // result of error indicators
  for( auto f=m_fd.Nbfac(); f<esuf.size()/2; ++f )
  {
    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    if (m_ndof[el] > m_ndof[er])
      ndof[er] = m_ndof[el];

    if (m_ndof[el] < m_ndof[er])
      ndof[el] = m_ndof[er];
  }

  // Update number of degrees of freedom for each cell
  m_ndof = ndof;
}

void
DG::comlim( int fromch,
            const std::vector< std::size_t >& tetid,
            const std::vector< std::vector< tk::real > >& u,
            const std::vector< std::vector< tk::real > >& prim,
            const std::vector< std::size_t >& ndof )
// *****************************************************************************
//  Receive chare-boundary limiter ghost data from neighboring chares
//! \param[in] fromch Sender chare id
//! \param[in] tetid Ghost tet ids we receive solution data for
//! \param[in] u Limited high-order solution
//! \param[in] prim Limited high-order primitive quantities
//! \param[in] ndof Number of degrees of freedom for chare-boundary elements
//! \details This function receives contributions to the limited solution from
//!   fellow chares.
// *****************************************************************************
{
  Assert( u.size() == tetid.size(), "Size mismatch in DG::comlim()" );
  Assert( prim.size() == tetid.size(), "Size mismatch in DG::comlim()" );

  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

  if (pref && m_stage == 0)
    Assert( ndof.size() == tetid.size(), "Size mismatch in DG::comlim()" );

  // Find local-to-ghost tet id map for sender chare
  const auto& n = tk::cref_find( m_ghost, fromch );

  for (std::size_t i=0; i<tetid.size(); ++i) {
    auto j = tk::cref_find( n, tetid[i] );
    Assert( j >= m_fd.Esuel().size()/4, "Receiving solution non-ghost data" );
    auto b = tk::cref_find( m_bid, j );
    Assert( b < m_uc[1].size(), "Indexing out of bounds" );
    Assert( b < m_pc[1].size(), "Indexing out of bounds" );
    m_uc[1][b] = u[i];
    m_pc[1][b] = prim[i];
    if (pref && m_stage == 0) {
      Assert( b < m_ndofc[1].size(), "Indexing out of bounds" );
      m_ndofc[1][b] = ndof[i];
    }
  }

  // if we have received all solution ghost contributions from those chares we
  // communicate along chare-boundary faces with, solve the system
  if (++m_nlim == m_ghostData.size()) {
    m_nlim = 0;
    comlim_complete();
  }
}

void
DG::dt()
// *****************************************************************************
// Compute time step size
// *****************************************************************************
{
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

  auto d = Disc();

  // Combine own and communicated contributions of limited solution and degrees
  // of freedom in cells (if p-adaptive)
  for (const auto& b : m_bid) {
    Assert( m_uc[1][b.second].size() == m_u.nprop(), "ncomp size mismatch" );
    Assert( m_pc[1][b.second].size() == m_p.nprop(), "ncomp size mismatch" );
    for (std::size_t c=0; c<m_u.nprop(); ++c) {
      m_u(b.first,c,0) = m_uc[1][b.second][c];
    }
    for (std::size_t c=0; c<m_p.nprop(); ++c) {
      m_p(b.first,c,0) = m_pc[1][b.second][c];
    }
    if (pref && m_stage == 0) {
      m_ndof[ b.first ] = m_ndofc[1][ b.second ];
    }
  }

  auto mindt = std::numeric_limits< tk::real >::max();

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
        auto eqdt =
          eq.dt( d->Coord(), d->Inpoel(), m_fd, m_geoFace, m_geoElem, m_ndof, m_u );
        if (eqdt < mindt) mindt = eqdt;
      }

      auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
      tk::real dgp = 0.0;

      if (ndof == 4)
      {
        dgp = 1.0;
      }
      else if (ndof == 10)
      {
        dgp = 2.0;
      }

      // Scale smallest dt with CFL coefficient and the CFL is scaled by (2*p+1)
      // where p is the order of the DG polynomial by linear stability theory.
      mindt *= g_inputdeck.get< tag::discr, tag::cfl >() / (2.0*dgp + 1.0);

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
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
  const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
  const auto neq = m_u.nprop()/rdof;

  // Set new time step size
  if (m_stage == 0) d->setdt( newdt );

  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();
  if (pref && m_stage == 0)
  {
    // When the element are coarsened, high order terms should be zero
    for(std::size_t e = 0; e < m_nunk; e++)
    {
      const auto ncomp= m_u.nprop()/rdof;
      if(m_ndof[e] == 1)
      {
        for (std::size_t c=0; c<ncomp; ++c)
        {
          auto mark = c*rdof;
          m_u(e, mark+1, 0) = 0.0;
          m_u(e, mark+2, 0) = 0.0;
          m_u(e, mark+3, 0) = 0.0;
        }
      } else if(m_ndof[e] == 4)
      {
        for (std::size_t c=0; c<ncomp; ++c)
        {
          auto mark = c*ndof;
          m_u(e, mark+4, 0) = 0.0;
          m_u(e, mark+5, 0) = 0.0;
          m_u(e, mark+6, 0) = 0.0;
          m_u(e, mark+7, 0) = 0.0;
          m_u(e, mark+8, 0) = 0.0;
          m_u(e, mark+9, 0) = 0.0;
        }
      }
    }
  }

  // Update Un
  if (m_stage == 0) m_un = m_u;

  for (const auto& eq : g_dgpde)
    eq.rhs( d->T(), m_geoFace, m_geoElem, m_fd, d->Inpoel(), d->Coord(), m_u,
            m_p, m_ndof, m_rhs );

  // Explicit time-stepping using RK3 to discretize time-derivative
  for(std::size_t e=0; e<m_nunk; ++e)
    for(std::size_t c=0; c<neq; ++c)
      for (std::size_t k=0; k<ndof; ++k)
      {
        auto rmark = c*rdof+k;
        auto mark = c*ndof+k;
        m_u(e, rmark, 0) =  rkcoef[0][m_stage] * m_un(e, rmark, 0)
          + rkcoef[1][m_stage] * ( m_u(e, rmark, 0)
            + d->Dt() * m_rhs(e, mark, 0)/m_lhs(e, mark, 0) );
      }

  // Update primitives based on the evolved solution
  for (const auto& eq : g_dgpde)
    eq.updatePrimitives( m_u, m_p, m_fd.Esuel().size()/4 );

  if (m_stage < 2) {

    // continue with next tims step stage
    stage();

  } else {

    // Compute diagnostics, e.g., residuals
    auto diag_computed = m_diag.compute( *d, m_u.nunk()-m_fd.Esuel().size()/4,
                                         m_geoElem, m_ndof, m_u );

    // Increase number of iterations and physical time
    d->next();

    // Continue to mesh refinement (if configured)
    if (!diag_computed) refine();

  }
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

  // if t>0 refinement enabled and we hit the dtref frequency
  if (dtref && !(d->It() % dtfreq)) {   // refine

    d->startvol();
    d->Ref()->dtref( m_fd.Bface(), {}, tk::remap(m_fd.Triinpoel(),d->Gid()) );
    d->refined() = 1;

  } else {      // do not refine

    d->refined() = 0;
    stage();

  }
}

void
DG::resizePostAMR(
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
//! \param[in] chunk New mesh chunk (connectivity and global<->local id maps)
//! \param[in] coord New mesh node coordinates
//! \param[in] addedTets Newly added mesh cells and their parents (local ids)
//! \param[in] msum New node communication map
//! \param[in] bface Boundary-faces mapped to side set ids
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
  d->resizePostAMR( chunk, coord, msum );

  // Update state
  auto nelem = d->Inpoel().size()/4;
  auto nprop = m_p.nprop();
  m_p.resize( nelem, nprop );
  nprop = m_u.nprop();
  m_u.resize( nelem, nprop );
  m_un.resize( nelem, nprop );
  m_lhs.resize( nelem, nprop );
  m_rhs.resize( nelem, nprop );

  m_fd = FaceData( d->Inpoel(), bface, tk::remap(triinpoel,d->Lid()) );

  m_geoFace =
    tk::Fields( tk::genGeoFaceTri( m_fd.Nipfac(), m_fd.Inpofa(), coord ) );
  m_geoElem = tk::Fields( tk::genGeoElemTet( d->Inpoel(), coord ) );

  m_nfac = m_fd.Inpofa().size()/3;
  m_nunk = nelem;
  m_ncoord = coord[0].size();
  m_msumset = d->msumset();
  m_bndFace.clear();
  m_ghostData.clear();
  m_ghost.clear();

  // Update solution on new mesh, P0 (cell center value) only for now
  m_un = m_u;
  auto pn = m_p;
  for (const auto& e : addedTets) {
    Assert( e.first < nelem, "Indexing out of new solution vector" );
    Assert( e.second < old_nelem, "Indexing out of old solution vector" );
    for (std::size_t c=0; c<nprop; ++c)
      m_u(e.first,c,0) = m_un(e.second,c,0);
    for (std::size_t c=0; c<m_p.nprop(); ++c)
      m_p(e.first,c,0) = pn(e.second,c,0);
  }
  m_un = m_u;

  // Enable SDAG wait for setting up chare boundary faces
  thisProxy[ thisIndex ].wait4fac();

  // Resize communication buffers
  resizeComm();
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
DG::evalLB()
// *****************************************************************************
// Evaluate whether to do load balancing
// *****************************************************************************
{
  auto d = Disc();

  const auto lbfreq = g_inputdeck.get< tag::cmd, tag::lbfreq >();
  const auto nonblocking = g_inputdeck.get< tag::cmd, tag::nonblocking >();

  // Load balancing if user frequency is reached or after the second time-step
  if ( (d->It()) % lbfreq == 0 || d->It() == 2 ) {

    AtSync();
    if (nonblocking) next();

  } else {

    next();

  }
}

void
DG::evalRestart()
// *****************************************************************************
// Evaluate whether to save checkpoint/restart
// *****************************************************************************
{
  auto d = Disc();

  const auto rsfreq = g_inputdeck.get< tag::cmd, tag::rsfreq >();

  if ( (d->It()) % rsfreq == 0 ) {

    std::vector< tk::real > t{{ static_cast<tk::real>(d->It()), d->T() }};
    d->contribute( t, CkReduction::nop,
      CkCallback(CkReductionTarget(Transporter,checkpoint), d->Tr()) );

  } else {

    evalLB();

  }
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

  // If neither max iterations nor max time reached, continue, otherwise finish
  if (std::fabs(d->T()-term) > eps && d->It() < nstep) {

    evalRestart();
 
  } else {

    std::vector< tk::real > t{{ static_cast<tk::real>(d->It()), d->T() }};
    d->contribute( t, CkReduction::nop,
      CkCallback(CkReductionTarget(Transporter,finish), d->Tr()) );

  }
}

#include "NoWarning/dg.def.h"
