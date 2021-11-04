// *****************************************************************************
/*!
  \file      src/Inciter/DG.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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
#include "Reorder.hpp"
#include "Vector.hpp"
#include "Around.hpp"
#include "Integrate/Basis.hpp"
#include "FieldOutput.hpp"
#include "ChareStateCollector.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< DGPDE > g_dgpde;

//! Runge-Kutta coefficients
static const std::array< std::array< tk::real, 3 >, 2 >
  rkcoef{{ {{ 0.0, 3.0/4.0, 1.0/3.0 }}, {{ 1.0, 1.0/4.0, 2.0/3.0 }} }};

} // inciter::

extern tk::CProxy_ChareStateCollector stateProxy;

using inciter::DG;

DG::DG( const CProxy_Discretization& disc,
        const std::map< int, std::vector< std::size_t > >& bface,
        const std::map< int, std::vector< std::size_t > >& /* bnode */,
        const std::vector< std::size_t >& triinpoel ) :
  m_disc( disc ),
  m_ndof_NodalExtrm( 3 ), // for the first order derivatives in 3 directions
  m_ncomfac( 0 ),
  m_nadj( 0 ),
  m_ncomEsup( 0 ),
  m_nsol( 0 ),
  m_ninitsol( 0 ),
  m_nlim( 0 ),
  m_nnod( 0 ),
  m_nreco( 0 ),
  m_nnodalExtrema( 0 ),
  m_inpoel( Disc()->Inpoel() ),
  m_coord( Disc()->Coord() ),
  m_fd( m_inpoel, bface, tk::remap(triinpoel,Disc()->Lid()) ),
  m_u( m_inpoel.size()/4,
       g_inputdeck.get< tag::discr, tag::rdof >()*
       g_inputdeck.get< tag::component >().nprop() ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_p( m_u.nunk(),
       g_inputdeck.get< tag::discr, tag::rdof >()*
         std::accumulate( begin(g_dgpde), end(g_dgpde), 0u,
           [](std::size_t s, const DGPDE& eq){ return s + eq.nprim(); } ) ),
  m_geoFace( tk::genGeoFaceTri( m_fd.Nipfac(), m_fd.Inpofa(), m_coord) ),
  m_geoElem( tk::genGeoElemTet( m_inpoel, m_coord ) ),
  m_lhs( m_u.nunk(),
         g_inputdeck.get< tag::discr, tag::ndof >()*
         g_inputdeck.get< tag::component >().nprop() ),
  m_rhs( m_u.nunk(), m_lhs.nprop() ),
  m_uNodalExtrm(),
  m_pNodalExtrm(),
  m_uNodalExtrmc(),
  m_pNodalExtrmc(),
  m_nfac( m_fd.Inpofa().size()/3 ),
  m_nunk( m_u.nunk() ),
  m_npoin( m_coord[0].size() ),
  m_ipface(),
  m_bndFace(),
  m_ghostData(),
  m_sendGhost(),
  m_ghostReq( 0 ),
  m_ghost(),
  m_exptGhost(),
  m_recvGhost(),
  m_diag(),
  m_stage( 0 ),
  m_ndof(),
  m_numEqDof(),
  m_bid(),
  m_uc(),
  m_pc(),
  m_ndofc(),
  m_initial( 1 ),
  m_expChBndFace(),
  m_infaces(),
  m_esup(),
  m_esupc(),
  m_elemfields(),
  m_nodefields(),
  m_nodefieldsc(),
  m_outmesh(),
  m_boxelems(),
  m_shockmarker(m_u.nunk())
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::chare >() ||
      g_inputdeck.get< tag::cmd, tag::quiescence >())
    stateProxy.ckLocalBranch()->insert( "DG", thisIndex, CkMyPe(), Disc()->It(),
                                        "DG" );

  // assign number of dofs for each equation in all pde systems
  for (const auto& eq : g_dgpde) {
    eq.numEquationDofs(m_numEqDof);
  }

  // Allocate storage for the vector of nodal extrema
  m_uNodalExtrm.resize( Disc()->Bid().size(), std::vector<tk::real>( 2*
    m_ndof_NodalExtrm*g_inputdeck.get< tag::component >().nprop() ) );
  m_pNodalExtrm.resize( Disc()->Bid().size(), std::vector<tk::real>( 2*
    m_ndof_NodalExtrm*m_p.nprop()/g_inputdeck.get< tag::discr, tag::rdof >()));

  // Initialization for the buffer vector of nodal extrema
  resizeNodalExtremac();

  usesAtSync = true;    // enable migration at AtSync

  // Enable SDAG wait for setting up chare boundary faces
  thisProxy[ thisIndex ].wait4fac();

  // Ensure that mesh partition is not leaky
  Assert( !tk::leakyPartition(m_fd.Esuel(), m_inpoel, m_coord),
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
  s.push_back( static_cast< tk::real >( Disc()->MeshId() ) );

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
  const auto& inpofa = m_fd.Inpofa();
  const auto& esuel = m_fd.Esuel();

  // Perform leak test on mesh partition
  Assert( !tk::leakyPartition( esuel, m_inpoel, m_coord ),
          "Mesh partition leaky" );

  // Activate SDAG waits for face adjacency map (ghost data) calculation
  thisProxy[ thisIndex ].wait4ghost();
  thisProxy[ thisIndex ].wait4esup();

  // Enable SDAG wait for initially building the solution vector and limiting
  if (m_initial) {
    thisProxy[ thisIndex ].wait4sol();
    thisProxy[ thisIndex ].wait4reco();
    thisProxy[ thisIndex ].wait4nodalExtrema();
    thisProxy[ thisIndex ].wait4lim();
    thisProxy[ thisIndex ].wait4nod();
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
        tk::UnsMesh::Face t{{ gid[ m_inpoel[ mark + tk::lpofa[f][0] ] ],
                              gid[ m_inpoel[ mark + tk::lpofa[f][1] ] ],
                              gid[ m_inpoel[ mark + tk::lpofa[f][2] ] ] }};
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
  // adjacency map. This is because while a node can be shared at a single
  // corner or along an edge, that does not necessarily share a face as well
  // (in other words, shared nodes or edges can exist that are not part of a
  // shared face). So the chares we communicate with across faces are not
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
  if (d->NodeCommMap().empty())  // in serial, skip setting up ghosts altogether
    faceAdj();
  else
    // for all chares we share nodes with
    for (const auto& c : d->NodeCommMap()) {
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
  bool match(true);

  auto eps = std::numeric_limits< tk::real >::epsilon() * 100;

  for (auto f=m_fd.Nipfac(); f<esuf.size()/2; ++f)
  {
    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    std::size_t count = 0;

    for (std::size_t i=0; i<4; ++i)
    {
      auto ip = m_inpoel[4*el+i];
      for (std::size_t j=0; j<4; ++j)
      {
        auto jp = m_inpoel[4*er+j];
        auto xdiff = std::abs( m_coord[0][ip] - m_coord[0][jp] );
        auto ydiff = std::abs( m_coord[1][ip] - m_coord[1][jp] );
        auto zdiff = std::abs( m_coord[2][ip] - m_coord[2][jp] );

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
  if (g_inputdeck.get< tag::cmd, tag::chare >() ||
      g_inputdeck.get< tag::cmd, tag::quiescence >())
    stateProxy.ckLocalBranch()->insert( "DG", thisIndex, CkMyPe(), Disc()->It(),
                                        "comfac" );

  // Buffer up incoming data
  m_infaces[ fromch ] = infaces;

  // if we have heard from all fellow chares that we share at least a single
  // node, edge, or face with
  if (++m_ncomfac == Disc()->NodeCommMap().size()) {
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
          tk::UnsMesh::Face t{{ gid[ m_inpoel[ mark + tk::lpofa[f][0] ] ],
                                gid[ m_inpoel[ mark + tk::lpofa[f][1] ] ],
                                gid[ m_inpoel[ mark + tk::lpofa[f][2] ] ] }};
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
        tk::UnsMesh::Face t{{ gid[ m_inpoel[ mark + tk::lpofa[f][0] ] ],
                              gid[ m_inpoel[ mark + tk::lpofa[f][1] ] ],
                              gid[ m_inpoel[ mark + tk::lpofa[f][2] ] ] }};
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
  // comfac() we are calling reqGhost() by going through the node communication
  // map instead, which may send requests to those chare we do not share faces
  // with. This is so that we can test for completing by querying the size of
  // the already complete node commincation map in reqGhost. Requests in
  // sendGhost will only be fullfilled based on m_ghostData.
  for (const auto& c : d->NodeCommMap())  // for all chares we share nodes with
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
       const auto& x = m_coord[0];
       const auto& y = m_coord[1];
       const auto& z = m_coord[2];
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
        tk::UnsMesh::Face t{{ gid[ m_inpoel[ mark + tk::lpofa[f][0] ] ],
                              gid[ m_inpoel[ mark + tk::lpofa[f][1] ] ],
                              gid[ m_inpoel[ mark + tk::lpofa[f][2] ] ] }};
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
            ncoord[0] = m_coord[0][ m_inpoel[ mark+f ] ];
            ncoord[1] = m_coord[1][ m_inpoel[ mark+f ] ];
            ncoord[2] = m_coord[2][ m_inpoel[ mark+f ] ];

            std::get< 3 >( tuple ) = f;

            std::get< 4 >( tuple ) = {{ gid[ m_inpoel[ mark ] ],
                                        gid[ m_inpoel[ mark+1 ] ],
                                        gid[ m_inpoel[ mark+2 ] ],
                                        gid[ m_inpoel[ mark+3 ] ] }};
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
  if (g_inputdeck.get< tag::cmd, tag::chare >() ||
      g_inputdeck.get< tag::cmd, tag::quiescence >())
    stateProxy.ckLocalBranch()->insert( "DG", thisIndex, CkMyPe(), Disc()->It(),
                                        "reqGhost" );

  // If every chare we communicate with has requested ghost data from us, we may
  // fulfill the requests, but only if we have already setup our ghost data.
  if (++m_ghostReq == Disc()->NodeCommMap().size()) {
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
  if (g_inputdeck.get< tag::cmd, tag::chare >() ||
      g_inputdeck.get< tag::cmd, tag::quiescence >())
    stateProxy.ckLocalBranch()->insert( "DG", thisIndex, CkMyPe(), Disc()->It(),
                                        "sendGhost" );

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
  if (g_inputdeck.get< tag::cmd, tag::chare >() ||
      g_inputdeck.get< tag::cmd, tag::quiescence >())
    stateProxy.ckLocalBranch()->insert( "DG", thisIndex, CkMyPe(), Disc()->It(),
                                        "comGhost" );

  auto d = Disc();
  const auto& lid = d->Lid();
  auto& inpofa = m_fd.Inpofa();
  auto ncoord = m_coord[0].size();

  // nodelist with fromch, currently only used for an assert
  [[maybe_unused]] const auto& nl = tk::cref_find( d->NodeCommMap(), fromch );

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
    Assert( geo.size() % 5 == 0, "Ghost geometry size mismatch" );
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
          m_inpoel.push_back( lp );       // store ghost element connectivity
        }
        // only a single or no ghost node should be found
        Assert( counter <= 1, "Incorrect number of ghost nodes detected. "
                "Detected "+ std::to_string(counter) +" ghost nodes" );
        if (counter == 1) {
          m_coord[0].push_back( coordg[0] ); // store ghost node coordinate
          m_coord[1].push_back( coordg[1] );
          m_coord[2].push_back( coordg[2] );
          Assert( m_inpoel[ 4*(m_nunk-1)+std::get< 3 >( g.second ) ] == ncoord,
                  "Mismatch in extended inpoel for ghost element" );
          ++ncoord;                // increase number of nodes on this chare
        }
      }

      // additional tests to ensure that entries in inpoel and t/inpofa match
      Assert( nodetripletMatch(id, t) == 3, "Mismatch/Overmatch in inpoel and "
              "inpofa at chare-boundary face" );
    }
  }

  // Signal the runtime system that all workers have received their
  // face-adjacency
  if (++m_nadj == m_ghostData.size()) faceAdj();
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
  const auto& esuf = m_fd.Esuf();
  const auto& inpofa = m_fd.Inpofa();

  std::size_t counter = 0;
  for (std::size_t k=0; k<4; ++k)
  {
    auto el = esuf[ 2*id[0] ];
    auto ip = m_inpoel[ 4*static_cast< std::size_t >( el )+k ];
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
  [[maybe_unused]] const auto& esuf = m_fd.Esuf();
  const auto& lid = d->Lid();

  std::array< tk::UnsMesh::Face, 4 > face;
  for (std::size_t f = 0; f<4; ++f)
    for (std::size_t i = 0; i<3; ++i)
      face[f][i] = m_inpoel[ id[1]*4 + tk::lpofa[f][i] ];

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
  const auto& lid = d->Lid();

  // get global node IDs reversing order to get outward-pointing normal
  auto A = tk::cref_find( lid, t[2] );
  auto B = tk::cref_find( lid, t[1] );
  auto C = tk::cref_find( lid, t[0] );
  auto geochf = tk::geoFaceTri( {{m_coord[0][A], m_coord[0][B], m_coord[0][C]}},
                                {{m_coord[1][A], m_coord[1][B], m_coord[1][C]}},
                                {{m_coord[2][A], m_coord[2][B], m_coord[2][C]}} );

  for (std::size_t i=0; i<7; ++i)
    m_geoFace(id[0],i,0) = geochf(0,i,0);
}

void
DG::faceAdj()
// *****************************************************************************
// Continue after face adjacency communication map completed on this chare
//! \details At this point the face communication map has been established
//!    on this chare. Proceed to set up the nodal-comm map.
// *****************************************************************************
{
  m_nadj = 0;

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

  // Create new map of elements along chare boundary which are ghosts for
  // neighboring chare, associated with that chare ID
  for (const auto& [cid, cgd] : m_ghostData)
  {
    auto& sg = m_sendGhost[cid];
    for (const auto& e : cgd)
    {
      Assert(sg.find(e.first) == sg.end(), "Repeating element found in "
        "ghost data");
      sg.insert(e.first);
    }
    Assert(sg.size() == cgd.size(), "Incorrect size for sendGhost");
  }
  Assert(m_sendGhost.size() == m_ghostData.size(), "Incorrect number of "
    "chares in sendGhost");

  // Error checking on ghost data
  for(const auto& n : m_sendGhost)
    for([[maybe_unused]] const auto& i : n.second)
      Assert( i < m_fd.Esuel().size()/4, "Sender contains ghost tet id. " );

  // Generate and store Esup data-structure in a map
  auto esup = tk::genEsup(m_inpoel, 4);
  for (std::size_t p=0; p<Disc()->Gid().size(); ++p)
  {
    for (auto e : tk::Around(esup, p))
    {
      // since inpoel has been augmented with the face-ghost cell previously,
      // esup also contains cells which are not on this mesh-chunk, hence the
      // following test
      if (e < m_fd.Esuel().size()/4) m_esup[p].push_back(e);
    }
  }

  // Error checking on Esup map
  for(const auto& p : m_esup)
    for([[maybe_unused]] const auto& e : p.second)
      Assert( e < m_fd.Esuel().size()/4, "Esup contains tet id greater than "
      + std::to_string(m_fd.Esuel().size()/4-1) +" : "+ std::to_string(e) );

  auto meshid = Disc()->MeshId();
  contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
    CkCallback(CkReductionTarget(Transporter,startEsup), Disc()->Tr()) );
}

void
DG::nodeNeighSetup()
// *****************************************************************************
// Setup node-neighborhood (esup)
//! \details At this point the face-ghost communication map has been established
//!    on this chare. This function begins generating the node-ghost comm map.
// *****************************************************************************
{
  if (Disc()->NodeCommMap().empty())
  // in serial, skip setting up node-neighborhood
  { comesup_complete(); }
  else
  {
    const auto& nodeCommMap = Disc()->NodeCommMap();

    // send out node-neighborhood map
    for (const auto& [cid, nlist] : nodeCommMap)
    {
      std::unordered_map< std::size_t, std::vector< std::size_t > > bndEsup;
      std::unordered_map< std::size_t, std::vector< tk::real > > nodeBndCells;
      for (const auto& p : nlist)
      {
        auto pl = tk::cref_find(Disc()->Lid(), p);
        // fill in the esup for the chare-boundary
        const auto& pesup = tk::cref_find(m_esup, pl);
        bndEsup[p] = pesup;

        // fill a map with the element ids from esup as keys and geoElem as
        // values, and another map containing these elements associated with
        // the chare id with which they are node-neighbors.
        for (const auto& e : pesup)
        {
          nodeBndCells[e] = m_geoElem[e];

          // add these esup-elements into map of elements along chare boundary
          Assert( e < m_fd.Esuel().size()/4, "Sender contains ghost tet id." );
          m_sendGhost[cid].insert(e);
        }
      }

      thisProxy[cid].comEsup(thisIndex, bndEsup, nodeBndCells);
    }
  }

  ownesup_complete();
}

void
DG::comEsup( int fromch,
  const std::unordered_map< std::size_t, std::vector< std::size_t > >& bndEsup,
  const std::unordered_map< std::size_t, std::vector< tk::real > >&
    nodeBndCells )
// *****************************************************************************
//! \brief Receive elements-surrounding-points data-structure for points on
//    common boundary between receiving and sending neighbor chare, and the
//    element geometries for these new elements
//! \param[in] fromch Sender chare id
//! \param[in] bndEsup Elements-surrounding-points data-structure from fromch
//! \param[in] nodeBndCells Map containing element geometries associated with
//!   remote element IDs in the esup
// *****************************************************************************
{
  auto& chghost = m_ghost[fromch];

  // Extend remote-local element id map and element geometry array
  for (const auto& e : nodeBndCells)
  {
    // need to check following, because 'e' could have been added previously in
    // remote-local element id map as a part of face-communication, i.e. as a
    // face-ghost element
    if (chghost.find(e.first) == chghost.end())
    {
      chghost[e.first] = m_nunk;
      m_geoElem.push_back(e.second);
      ++m_nunk;
    }
  }

  // Store incoming data in comm-map buffer for Esup
  for (const auto& [node, elist] : bndEsup)
  {
    auto pl = tk::cref_find(Disc()->Lid(), node);
    auto& pesup = m_esupc[pl];
    for (auto e : elist)
    {
      auto el = tk::cref_find(chghost, e);
      pesup.push_back(el);
    }
  }

  // if we have heard from all fellow chares that we share at least a single
  // node, edge, or face with
  if (++m_ncomEsup == Disc()->NodeCommMap().size()) {
    m_ncomEsup = 0;
    comesup_complete();
  }
}

void
DG::adj()
// *****************************************************************************
// Finish up with adjacency maps, and do a global-sync to begin problem setup
//! \details At this point, the nodal- and face-adjacency has been set up. This
//    function does some error checking on the nodal-adjacency and prepares
//    for problem setup.
// *****************************************************************************
{
  // combine own and communicated contributions to elements surrounding points
  for (auto& [p, elist] : m_esupc)
  {
    auto& pesup = tk::ref_find(m_esup, p);
    for ([[maybe_unused]] auto e : elist)
    {
      Assert( e >= m_fd.Esuel().size()/4, "Non-ghost element received from "
        "esup buffer." );
    }
    tk::concat< std::size_t >(std::move(elist), pesup);
  }

  tk::destroy(m_ghostData);
  tk::destroy(m_esupc);

  if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) Disc()->Tr().chadj();

  // Error checking on ghost data
  for(const auto& n : m_sendGhost)
    for([[maybe_unused]] const auto& i : n.second)
      Assert( i < m_fd.Esuel().size()/4, "Sender contains ghost tet id. ");

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

  // Basic error checking on ghost tet ID map
  Assert( m_ghost.find( thisIndex ) == m_ghost.cend(),
          "Ghost id map should not contain data for own chare ID" );

  // Store expected ghost tet IDs
  for (const auto& n : m_ghost)
    for ([[maybe_unused]] const auto& g : n.second)
      Assert( m_exptGhost.insert( g.second ).second,
              "Failed to store local tetid as exptected ghost id" );

  // Signal the runtime system that all workers have received their adjacency
  std::vector< std::size_t > meshdata{ m_initial, Disc()->MeshId() };
  contribute( meshdata, CkReduction::sum_ulong,
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
  if (g_inputdeck.get< tag::cmd, tag::chare >() ||
      g_inputdeck.get< tag::cmd, tag::quiescence >())
    stateProxy.ckLocalBranch()->insert( "DG", thisIndex, CkMyPe(), Disc()->It(),
                                        "setup" );

  auto d = Disc();

  // Basic error checking on sizes of element geometry data and connectivity
  Assert( m_geoElem.nunk() == m_lhs.nunk(), "Size mismatch in DG::setup()" );

  // Compute left-hand side of discrete PDEs
  lhs();

  // Determine elements inside user-defined IC box
  for (auto& eq : g_dgpde)
    eq.IcBoxElems( m_geoElem, m_fd.Esuel().size()/4, m_boxelems );

  // Compute volume of user-defined box IC
  d->boxvol( {} );      // punt for now

  // Query time history field output labels from all PDEs integrated
  const auto& hist_points = g_inputdeck.get< tag::history, tag::point >();
  if (!hist_points.empty()) {
    std::vector< std::string > histnames;
    for (const auto& eq : g_dgpde) {
      auto n = eq.histNames();
      histnames.insert( end(histnames), begin(n), end(n) );
    }
    d->histheader( std::move(histnames) );
  }
}

void
DG::box( tk::real v )
// *****************************************************************************
// Receive total box IC volume and set conditions in box
//! \param[in] v Total volume within user-specified box
// *****************************************************************************
{
  auto d = Disc();

  // Store user-defined box IC volume
  d->Boxvol() = v;

  // Set initial conditions for all PDEs
  for (const auto& eq : g_dgpde)
  {
    eq.initialize( m_lhs, m_inpoel, m_coord, m_boxelems, m_u, d->T(),
      m_fd.Esuel().size()/4 );
    eq.updatePrimitives( m_u, m_lhs, m_geoElem, m_p, m_fd.Esuel().size()/4 );
  }

  m_un = m_u;

  // Output initial conditions to file (regardless of whether it was requested)
  startFieldOutput( CkCallback(CkIndex_DG::start(), thisProxy[thisIndex]) );
}

void
DG::start()
// *****************************************************************************
//  Start time stepping
// *****************************************************************************
{
  // Free memory storing output mesh
  m_outmesh.destroy();

  // Start timer measuring time stepping wall clock time
  Disc()->Timer().zero();
  // Zero grind-timer
  Disc()->grindZero();
  // Start time stepping by computing the size of the next time step)
  next();
}

void
DG::startFieldOutput( CkCallback c )
// *****************************************************************************
// Start preparing fields for output to file
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  // No field output in benchmark mode or if field output frequency not hit
  if (g_inputdeck.get< tag::cmd, tag::benchmark >() || !fieldOutput()) {

    c.send();

  } else {

    // Optionally refine mesh for field output
    auto d = Disc();

    if (refinedOutput()) {

      const auto& tr = tk::remap( m_fd.Triinpoel(), d->Gid() );
      d->Ref()->outref( m_fd.Bface(), {}, tr, c );

    } else {

      // cut off ghosts from mesh connectivity and coordinates
      const auto& tr = tk::remap( m_fd.Triinpoel(), d->Gid() );
      extractFieldOutput( {}, d->Chunk(), d->Coord(), {}, {},
                          d->NodeCommMap(), m_fd.Bface(), {}, tr, c );

    }

  }
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
    for (const auto& eq : g_dgpde)
      eq.eval_ndof( m_nunk, m_coord, m_inpoel, m_fd, m_u,
                    g_inputdeck.get< tag::pref, tag::indicator >(),
                    g_inputdeck.get< tag::discr, tag::ndof >(),
                    g_inputdeck.get< tag::pref, tag::ndofmax >(),
                    g_inputdeck.get< tag::pref, tag::tolref >(),
                    m_ndof );

  // communicate solution ghost data (if any)
  if (m_sendGhost.empty())
    comsol_complete();
  else
    for(const auto& [cid, ghostdata] : m_sendGhost) {
      std::vector< std::size_t > tetid( ghostdata.size() );
      std::vector< std::vector< tk::real > > u( ghostdata.size() ),
                                             prim( ghostdata.size() );
      std::vector< std::size_t > ndof;
      std::size_t j = 0;
      for(const auto& i : ghostdata) {
        Assert( i < m_fd.Esuel().size()/4, "Sending solution ghost data" );
        tetid[j] = i;
        u[j] = m_u[i];
        prim[j] = m_p[i];
        if (pref && m_stage == 0) ndof.push_back( m_ndof[i] );
        ++j;
      }
      thisProxy[ cid ].comsol( thisIndex, m_stage, tetid, u, prim, ndof );
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

  // if we have received all solution ghost contributions from neighboring
  // chares (chares we communicate along chare-boundary faces with), and
  // contributed our solution to these neighbors, proceed to reconstructions
  if (++m_nsol == m_sendGhost.size()) {
    m_nsol = 0;
    comsol_complete();
  }
}

void
DG::extractFieldOutput(
  const std::vector< std::size_t >& /*ginpoel*/,
  const tk::UnsMesh::Chunk& chunk,
  const tk::UnsMesh::Coords& coord,
  const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& /*addedNodes*/,
  const std::unordered_map< std::size_t, std::size_t >& addedTets,
  const tk::NodeCommMap& nodeCommMap,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& /* bnode */,
  const std::vector< std::size_t >& triinpoel,
  CkCallback c )
// *****************************************************************************
// Extract field output going to file
//! \param[in] chunk Field-output mesh chunk (connectivity and global<->local
//!    id maps)
//! \param[in] coord Field-output mesh node coordinates
//! \param[in] addedTets Field-output mesh cells and their parents (local ids)
//! \param[in] nodeCommMap Field-output mesh node communication map
//! \param[in] bface Field-output meshndary-faces mapped to side set ids
//! \param[in] triinpoel Field-output mesh boundary-face connectivity
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  m_outmesh.chunk = chunk;
  m_outmesh.coord = coord;
  m_outmesh.triinpoel = triinpoel;
  m_outmesh.bface = bface;
  m_outmesh.nodeCommMap = nodeCommMap;

  const auto& inpoel = std::get< 0 >( chunk );
  auto nelem = inpoel.size() / 4;

  // Evaluate element solution on incoming mesh
  auto [ue,pe,un,pn] = evalSolution( inpoel, coord, addedTets );

  // Collect field output from numerical solution requested by user
  m_elemfields = numericFieldOutput( ue, tk::Centering::ELEM, pe );
  m_nodefields = numericFieldOutput( un, tk::Centering::NODE, pn );

  // Collect field output from analytical solutions (if exist)
  auto geoElem = tk::genGeoElemTet( inpoel, coord );
  auto t = Disc()->T();
  for (const auto& eq : g_dgpde) {
    analyticFieldOutput( eq, tk::Centering::ELEM, geoElem.extract(1,0),
      geoElem.extract(2,0), geoElem.extract(3,0), t, m_elemfields );
    analyticFieldOutput( eq, tk::Centering::NODE, coord[0], coord[1], coord[2],
      t, m_nodefields );
  }

  // Add adaptive indicator array to element-centered field output
  if (g_inputdeck.get< tag::pref, tag::pref >()) {
    std::vector< tk::real > ndof( begin(m_ndof), end(m_ndof) );
    ndof.resize( nelem );
    for (const auto& [child,parent] : addedTets)
      ndof[child] = static_cast< tk::real >( m_ndof[parent] );
    m_elemfields.push_back( ndof );
  }

  // Add shock detection marker array to element-centered field output
  std::vector< tk::real > shockmarker( begin(m_shockmarker), end(m_shockmarker) );
  // Here m_shockmarker has a size of m_u.nunk() which is the number of the
  // elements within this partition (nelem) plus the ghost partition cells. In
  // terms of output purpose, we only need the solution data within this
  // partition. Therefore, resizing it to nelem removes the extra partition
  // boundary allocations in the shockmarker vector. Since the code assumes that
  // the boundary elements are on the top, the resize operation keeps the lower
  // portion.
  shockmarker.resize( nelem );
  for (const auto& [child,parent] : addedTets)
    shockmarker[child] = static_cast< tk::real >(m_shockmarker[parent]);
  m_elemfields.push_back( shockmarker );

  // Send node fields contributions to neighbor chares
  if (nodeCommMap.empty())
    comnodeout_complete();
  else {
    const auto& lid = std::get< 2 >( chunk );
    auto esup = tk::genEsup( inpoel, 4 );
    for(const auto& [ch,nodes] : nodeCommMap) {
      // Pack node field data in chare boundary nodes
      std::vector< std::vector< tk::real > >
        l( m_nodefields.size(), std::vector< tk::real >( nodes.size() ) );
      for (std::size_t f=0; f<m_nodefields.size(); ++f) {
        std::size_t j = 0;
        for (auto g : nodes)
          l[f][j++] = m_nodefields[f][ tk::cref_find(lid,g) ];
      }
      // Pack (partial) number of elements surrounding chare boundary nodes
      std::vector< std::size_t > nesup( nodes.size() );
      std::size_t j = 0;
      for (auto g : nodes) {
        auto i = tk::cref_find( lid, g );
        nesup[j++] = esup.second[i+1] - esup.second[i];
      }
      thisProxy[ch].comnodeout(
        std::vector<std::size_t>(begin(nodes),end(nodes)), nesup, l );
    }
  }

  ownnod_complete( c );
}

std::tuple< tk::Fields, tk::Fields, tk::Fields, tk::Fields >
DG::evalSolution(
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const std::unordered_map< std::size_t, std::size_t >& addedTets )
// *****************************************************************************
// Evaluate solution on incomping (a potentially refined) mesh
//! \param[in] inpoel Incoming (potentially refined field-output) mesh
//!   connectivity
//! \param[in] coord Incoming (potentially refined Field-output) mesh node
//!   coordinates
//! \param[in] addedTets Field-output mesh cells and their parents (local ids)
//! \details This function evaluates the solution on the incoming mesh. The
//!   incoming mesh can be refined but can also be just the mesh the numerical
//!   solution is computed on.
//! \note If the incoming mesh is refined (for field putput) compared to the
//!   mesh the numerical solution is computed on, the solution is evaluated in
//!   cells as wells as in nodes. If the solution is not refined, the solution
//!   is evaluated in nodes.
//! \return Solution in cells, primitive variables in cells, solution in nodes,
//!   primitive variables in nodes of incoming mesh.
// *****************************************************************************
{
  using tk::dot;
  using tk::real;

  const auto nelem = inpoel.size()/4;
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
  const auto uncomp = m_u.nprop() / rdof;
  const auto pncomp = m_p.nprop() / rdof;
  auto ue = m_u;
  auto pe = m_p;

  // If mesh is not refined for field output, cut off ghosts from element
  // solution. (No need to output ghosts and writer would error.) If mesh is
  // refined for field output, resize element solution fields to refined mesh.
  ue.resize( nelem );
  pe.resize( nelem );

  auto npoin = coord[0].size();
  tk::Fields un( npoin, m_u.nprop()/rdof );
  tk::Fields pn( npoin, m_p.nprop()/rdof );
  un.fill(0.0);
  pn.fill(0.0);

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // If mesh is not refined for output, evaluate solution in nodes
  if (addedTets.empty()) {

    for (std::size_t e=0; e<nelem; ++e) {
      auto e4 = e*4;
      // Extract element node coordinates
      std::array< std::array< real, 3>, 4 > ce{{
        {{ x[inpoel[e4  ]], y[inpoel[e4  ]], z[inpoel[e4  ]] }},
        {{ x[inpoel[e4+1]], y[inpoel[e4+1]], z[inpoel[e4+1]] }},
        {{ x[inpoel[e4+2]], y[inpoel[e4+2]], z[inpoel[e4+2]] }},
        {{ x[inpoel[e4+3]], y[inpoel[e4+3]], z[inpoel[e4+3]] }} }};
      // Compute inverse Jacobian
      auto J = tk::inverseJacobian( ce[0], ce[1], ce[2], ce[3] );
      // Evaluate solution in child nodes
      for (std::size_t j=0; j<4; ++j) {
        std::array< real, 3 >
           h{{ce[j][0]-ce[0][0], ce[j][1]-ce[0][1], ce[j][2]-ce[0][2] }};
        auto Bn = tk::eval_basis( m_ndof[e],
                                  dot(J[0],h), dot(J[1],h), dot(J[2],h) );
        auto u = eval_state( uncomp, 0, rdof, m_ndof[e], e, m_u, Bn, {0, uncomp-1} );
        auto p = eval_state( pncomp, 0, rdof, m_ndof[e], e, m_p, Bn, {0, pncomp-1} );
        // Assign child node solution
        for (std::size_t i=0; i<uncomp; ++i) un(inpoel[e4+j],i,0) += u[i];
        for (std::size_t i=0; i<pncomp; ++i) pn(inpoel[e4+j],i,0) += p[i];
      }
    }

  // If mesh is refed for output, evaluate solution in elements and nodes of
  // refined mesh
  } else {

    const auto& pinpoel = Disc()->Inpoel();  // unrefined (parent) mesh

    for ([[maybe_unused]] const auto& [child,parent] : addedTets) {
      Assert( child < nelem, "Indexing out of new solution vector" );
      Assert( parent < pinpoel.size()/4,
              "Indexing out of old solution vector" );
    }

    for (const auto& [child,parent] : addedTets) {
      // Extract parent element's node coordinates
      auto p4 = 4*parent;
      std::array< std::array< real, 3>, 4 > cp{{
        {{ x[pinpoel[p4  ]], y[pinpoel[p4  ]], z[pinpoel[p4  ]] }},
        {{ x[pinpoel[p4+1]], y[pinpoel[p4+1]], z[pinpoel[p4+1]] }},
        {{ x[pinpoel[p4+2]], y[pinpoel[p4+2]], z[pinpoel[p4+2]] }},
        {{ x[pinpoel[p4+3]], y[pinpoel[p4+3]], z[pinpoel[p4+3]] }} }};
      // Evaluate inverse Jacobian of the parent
      auto Jp = tk::inverseJacobian( cp[0], cp[1], cp[2], cp[3] );
      // Compute child cell centroid
      auto c4 = 4*child;
      auto cx = (x[inpoel[c4  ]] + x[inpoel[c4+1]] +
                 x[inpoel[c4+2]] + x[inpoel[c4+3]]) / 4.0;
      auto cy = (y[inpoel[c4  ]] + y[inpoel[c4+1]] +
                 y[inpoel[c4+2]] + y[inpoel[c4+3]]) / 4.0;
      auto cz = (z[inpoel[c4  ]] + z[inpoel[c4+1]] +
                 z[inpoel[c4+2]] + z[inpoel[c4+3]]) / 4.0;
      // Compute solution in child centroid
      std::array< real, 3 > h{{cx-cp[0][0], cy-cp[0][1], cz-cp[0][2] }};
      auto B = tk::eval_basis( m_ndof[parent],
                               dot(Jp[0],h), dot(Jp[1],h), dot(Jp[2],h) );
      auto u = eval_state( uncomp, 0, rdof, m_ndof[parent], parent, m_u, B, {0, uncomp-1} );
      auto p = eval_state( pncomp, 0, rdof, m_ndof[parent], parent, m_p, B, {0, pncomp-1} );
      // Assign cell center solution from parent to child
      for (std::size_t i=0; i<uncomp; ++i) ue(child,i*rdof,0) = u[i];
      for (std::size_t i=0; i<pncomp; ++i) pe(child,i*rdof,0) = p[i];
      // Extract child element's node coordinates
      std::array< std::array< real, 3>, 4 > cc{{
        {{ x[inpoel[c4  ]], y[inpoel[c4  ]], z[inpoel[c4  ]] }},
        {{ x[inpoel[c4+1]], y[inpoel[c4+1]], z[inpoel[c4+1]] }},
        {{ x[inpoel[c4+2]], y[inpoel[c4+2]], z[inpoel[c4+2]] }},
        {{ x[inpoel[c4+3]], y[inpoel[c4+3]], z[inpoel[c4+3]] }} }};
      // Evaluate solution in child nodes
      for (std::size_t j=0; j<4; ++j) {
        std::array< real, 3 >
           hn{{cc[j][0]-cp[0][0], cc[j][1]-cp[0][1], cc[j][2]-cp[0][2] }};
        auto Bn = tk::eval_basis( m_ndof[parent],
                                  dot(Jp[0],hn), dot(Jp[1],hn), dot(Jp[2],hn) );
        auto cnu = eval_state(uncomp, 0, rdof, m_ndof[parent], parent, m_u, Bn, {0, uncomp-1});
        auto cnp = eval_state(pncomp, 0, rdof, m_ndof[parent], parent, m_p, Bn, {0, pncomp-1});
        // Assign child node solution
        for (std::size_t i=0; i<uncomp; ++i) un(inpoel[c4+j],i,0) += cnu[i];
        for (std::size_t i=0; i<pncomp; ++i) pn(inpoel[c4+j],i,0) += cnp[i];
      }
    }
  }

  return { ue, pe, un, pn };
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
DG::reco()
// *****************************************************************************
// Compute reconstructions
// *****************************************************************************
{
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();

  // Combine own and communicated contributions of unreconstructed solution and
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
    for (const auto& eq : g_dgpde)
      eq.reconstruct( d->T(), m_geoFace, m_geoElem, m_fd, m_esup, m_inpoel,
                      m_coord, m_u, m_p );
  }

  // Send reconstructed solution to neighboring chares
  if (m_sendGhost.empty())
    comreco_complete();
  else
    for(const auto& [cid, ghostdata] : m_sendGhost) {
      std::vector< std::size_t > tetid( ghostdata.size() );
      std::vector< std::vector< tk::real > > u( ghostdata.size() ),
                                             prim( ghostdata.size() );
      std::vector< std::size_t > ndof;
      std::size_t j = 0;
      for(const auto& i : ghostdata) {
        Assert( i < m_fd.Esuel().size()/4, "Sending reconstructed ghost "
          "data" );
        tetid[j] = i;
        u[j] = m_u[i];
        prim[j] = m_p[i];
        if (pref && m_stage == 0) ndof.push_back( m_ndof[i] );
        ++j;
      }
      thisProxy[ cid ].comreco( thisIndex, tetid, u, prim, ndof );
    }

  ownreco_complete();
}

void
DG::comreco( int fromch,
             const std::vector< std::size_t >& tetid,
             const std::vector< std::vector< tk::real > >& u,
             const std::vector< std::vector< tk::real > >& prim,
             const std::vector< std::size_t >& ndof )
// *****************************************************************************
//  Receive chare-boundary reconstructed ghost data from neighboring chares
//! \param[in] fromch Sender chare id
//! \param[in] tetid Ghost tet ids we receive solution data for
//! \param[in] u Reconstructed high-order solution
//! \param[in] prim Limited high-order primitive quantities
//! \param[in] ndof Number of degrees of freedom for chare-boundary elements
//! \details This function receives contributions to the reconstructed solution
//!   from fellow chares.
// *****************************************************************************
{
  Assert( u.size() == tetid.size(), "Size mismatch in DG::comreco()" );
  Assert( prim.size() == tetid.size(), "Size mismatch in DG::comreco()" );

  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

  if (pref && m_stage == 0)
    Assert( ndof.size() == tetid.size(), "Size mismatch in DG::comreco()" );

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

  // if we have received all solution ghost contributions from neighboring
  // chares (chares we communicate along chare-boundary faces with), and
  // contributed our solution to these neighbors, proceed to limiting
  if (++m_nreco == m_sendGhost.size()) {
    m_nreco = 0;
    comreco_complete();
  }
}

void
DG::nodalExtrema()
// *****************************************************************************
// Compute nodal extrema at chare-boundary nodes. Extrema at internal nodes
// are calculated in limiter function.
// *****************************************************************************
{
  auto d = Disc();
  auto gid = d->Gid();
  auto bid = d->Bid();
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();
  const auto ncomp = m_u.nprop() / rdof;
  const auto nprim = m_p.nprop() / rdof;

  // Combine own and communicated contributions of unlimited solution, and
  // if a p-adaptive algorithm is used, degrees of freedom in cells
  for (const auto& [boundary, localtet] : m_bid) {
    Assert( m_uc[1][localtet].size() == m_u.nprop(), "ncomp size mismatch" );
    Assert( m_pc[1][localtet].size() == m_p.nprop(), "ncomp size mismatch" );
    for (std::size_t c=0; c<m_u.nprop(); ++c) {
      m_u(boundary,c,0) = m_uc[1][localtet][c];
    }
    for (std::size_t c=0; c<m_p.nprop(); ++c) {
      m_p(boundary,c,0) = m_pc[1][localtet][c];
    }
    if (pref && m_stage == 0) {
      m_ndof[ boundary ] = m_ndofc[1][ localtet ];
    }
  }

  // Initialize nodal extrema vector
  auto large = std::numeric_limits< tk::real >::max();
  for(std::size_t i = 0; i<bid.size(); i++)
  {
    for (std::size_t c=0; c<ncomp; ++c)
    {
      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
      {
        auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
        auto min_mark = max_mark + 1;
        m_uNodalExtrm[i][max_mark] = -large;
        m_uNodalExtrm[i][min_mark] =  large;
      }
    }
    for (std::size_t c=0; c<nprim; ++c)
    {
      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
      {
        auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
        auto min_mark = max_mark + 1;
        m_pNodalExtrm[i][max_mark] = -large;
        m_pNodalExtrm[i][min_mark] =  large;
      }
    }
  }

  // Evaluate the max/min value for the chare-boundary nodes
  if(rdof > 4) {
      evalNodalExtrm(ncomp, nprim, m_ndof_NodalExtrm, d->bndel(), m_inpoel,
        m_coord, gid, bid, m_u, m_p, m_uNodalExtrm, m_pNodalExtrm);
  }

  // Communicate extrema at nodes to other chares on chare-boundary
  if (d->NodeCommMap().empty())        // in serial we are done
    comnodalExtrema_complete();
  else  // send nodal extrema to chare-boundary nodes to fellow chares
  {
    for (const auto& [c,n] : d->NodeCommMap()) {
      std::vector< std::vector< tk::real > > g1( n.size() ), g2( n.size() );
      std::size_t j = 0;
      for (auto i : n)
      {
        auto p = tk::cref_find(d->Bid(),i);
        g1[ j   ] = m_uNodalExtrm[ p ];
        g2[ j++ ] = m_pNodalExtrm[ p ];
      }
      thisProxy[c].comnodalExtrema( std::vector<std::size_t>(begin(n),end(n)),
        g1, g2 );
    }
  }
  ownnodalExtrema_complete();
}

void
DG::comnodalExtrema( const std::vector< std::size_t >& gid,
                     const std::vector< std::vector< tk::real > >& G1,
                     const std::vector< std::vector< tk::real > >& G2 )
// *****************************************************************************
//  Receive contributions to nodal extrema on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive grad contributions
//! \param[in] G1 Partial contributions of extrema for conservative variables to
//!   chare-boundary nodes
//! \param[in] G2 Partial contributions of extrema for primitive variables to
//!   chare-boundary nodes
//! \details This function receives contributions to m_uNodalExtrm/m_pNodalExtrm
//!   , which stores nodal extrems at mesh chare-boundary nodes. While
//!   m_uNodalExtrm/m_pNodalExtrm stores own contributions, m_uNodalExtrmc
//!   /m_pNodalExtrmc collects the neighbor chare contributions during
//!   communication.
// *****************************************************************************
{
  Assert( G1.size() == gid.size(), "Size mismatch" );
  Assert( G2.size() == gid.size(), "Size mismatch" );

  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
  const auto ncomp = m_u.nprop() / rdof;
  const auto nprim = m_p.nprop() / rdof;

  for (std::size_t i=0; i<gid.size(); ++i)
  {
    auto& u = m_uNodalExtrmc[gid[i]];
    auto& p = m_pNodalExtrmc[gid[i]];
    for (std::size_t c=0; c<ncomp; ++c)
    {
      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
      {
        auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
        auto min_mark = max_mark + 1;
        u[max_mark] = std::max( G1[i][max_mark], u[max_mark] );
        u[min_mark] = std::min( G1[i][min_mark], u[min_mark] );
      }
    }
    for (std::size_t c=0; c<nprim; ++c)
    {
      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
      {
        auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
        auto min_mark = max_mark + 1;
        p[max_mark] = std::max( G2[i][max_mark], p[max_mark] );
        p[min_mark] = std::min( G2[i][min_mark], p[min_mark] );
      }
    }
  }

  if (++m_nnodalExtrema == Disc()->NodeCommMap().size())
  {
    m_nnodalExtrema = 0;
    comnodalExtrema_complete();
  }
}

void DG::resizeNodalExtremac()
// *****************************************************************************
//  Resize the buffer vector of nodal extrema
// *****************************************************************************
{
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
  const auto ncomp = m_u.nprop() / rdof;
  const auto nprim = m_p.nprop() / rdof;

  auto large = std::numeric_limits< tk::real >::max();
  for (const auto& [c,n] : Disc()->NodeCommMap())
  {
    for (auto i : n) {
      auto& u = m_uNodalExtrmc[i];
      auto& p = m_pNodalExtrmc[i];
      u.resize( 2*m_ndof_NodalExtrm*ncomp, large );
      p.resize( 2*m_ndof_NodalExtrm*nprim, large );

      // Initialize the minimum nodal extrema
      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
      {
        for(std::size_t k = 0; k < ncomp; k++)
          u[2*k*m_ndof_NodalExtrm+2*idof] = -large;
        for(std::size_t k = 0; k < nprim; k++)
          p[2*k*m_ndof_NodalExtrm+2*idof] = -large;
      }
    }
  }
}

void DG::evalNodalExtrm( const std::size_t ncomp,
                         const std::size_t nprim,
                         const std::size_t ndof_NodalExtrm,
                         const std::vector< std::size_t >& bndel,
                         const std::vector< std::size_t >& inpoel,
                         const tk::UnsMesh::Coords& coord,
                         const std::vector< std::size_t >& gid,
                         const std::unordered_map< std::size_t, std::size_t >&
                           bid,
                         const tk::Fields& U,
                         const tk::Fields& P,
                         std::vector< std::vector<tk::real> >& uNodalExtrm,
                         std::vector< std::vector<tk::real> >& pNodalExtrm )
// *****************************************************************************
//  Compute the nodal extrema for chare-boundary nodes
//! \param[in] ncomp Number of conservative variables
//! \param[in] nprim Number of primitive variables
//! \param[in] ndof_NodalExtrm Degree of freedom for nodal extrema
//! \param[in] bndel List of elements contributing to chare-boundary nodes
//! \param[in] inpoel Element-node connectivity for element e
//! \param[in] coord Array of nodal coordinates
//! \param[in] gid Local->global node id map
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!   global node ids (key)
//! \param[in] U Vector of conservative variables
//! \param[in] P Vector of primitive variables
//! \param[in,out] uNodalExtrm Chare-boundary nodal extrema for conservative
//!   variables
//! \param[in,out] pNodalExtrm Chare-boundary nodal extrema for primitive
//!   variables
// *****************************************************************************
{
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();

  for (auto e : bndel)
  {
    // access node IDs
    const std::vector<std::size_t> N
      { inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };

    // Loop over nodes of element e
    for(std::size_t ip=0; ip<4; ++ip)
    {
      auto i = bid.find( gid[N[ip]] );
      if (i != end(bid))      // If ip is the chare boundary point
      {
        // If DG(P2) is applied, find the nodal extrema of the gradients of
        // conservative/primitive variables in the physical domain

        // Vector used to store the first order derivatives for both
        // conservative and primitive variables
        std::vector< std::array< tk::real, 3 > > gradc(ncomp, {0.0, 0.0, 0.0});
        std::vector< std::array< tk::real, 3 > > gradp(ncomp, {0.0, 0.0, 0.0});

        const auto& cx = coord[0];
        const auto& cy = coord[1];
        const auto& cz = coord[2];

        std::array< std::array< tk::real, 3>, 4 > coordel {{
          {{ cx[ N[0] ], cy[ N[0] ], cz[ N[0] ] }},
          {{ cx[ N[1] ], cy[ N[1] ], cz[ N[1] ] }},
          {{ cx[ N[2] ], cy[ N[2] ], cz[ N[2] ] }},
          {{ cx[ N[3] ], cy[ N[3] ], cz[ N[3] ] }}
        }};

        auto jacInv = tk::inverseJacobian( coordel[0], coordel[1],
          coordel[2], coordel[3] );

        // Compute the derivatives of basis functions
        auto dBdx = tk::eval_dBdx_p1( rdof, jacInv );

        std::array< std::vector< tk::real >, 3 > center;
        center[0].resize(1, 0.25);
        center[1].resize(1, 0.25);
        center[2].resize(1, 0.25);
        tk::eval_dBdx_p2(0, center, jacInv, dBdx);

        // Evaluate the first order derivative in physical domain
        for(std::size_t icomp = 0; icomp < ncomp; icomp++)
        {
          auto mark = icomp * rdof;
          for(std::size_t idir = 0; idir < 3; idir++)
          {
            gradc[icomp][idir] = 0;
            for(std::size_t idof = 1; idof < rdof; idof++)
              gradc[icomp][idir] += U(e, mark+idof, 0) * dBdx[idir][idof];
          }
        }
        for(std::size_t icomp = 0; icomp < nprim; icomp++)
        {
          auto mark = icomp * rdof;
          for(std::size_t idir = 0; idir < 3; idir++)
          {
            gradp[icomp][idir] = 0;
            for(std::size_t idof = 1; idof < rdof; idof++)
              gradp[icomp][idir] += P(e, mark+idof, 0) * dBdx[idir][idof];
          }
        }

        // Store the extrema for the gradients
        for (std::size_t c=0; c<ncomp; ++c)
        {
          for (std::size_t idof = 0; idof < ndof_NodalExtrm; idof++)
          {
            auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
            auto min_mark = max_mark + 1;
            auto& ex = uNodalExtrm[i->second];
            ex[max_mark] = std::max(ex[max_mark], gradc[c][idof-1]);
            ex[min_mark] = std::min(ex[min_mark], gradc[c][idof-1]);
          }
        }
        for (std::size_t c=0; c<nprim; ++c)
        {
          for (std::size_t idof = 0; idof < ndof_NodalExtrm; idof++)
          {
            auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
            auto min_mark = max_mark + 1;
            auto& ex = pNodalExtrm[i->second];
            ex[max_mark] = std::max(ex[max_mark], gradp[c][idof-1]);
            ex[min_mark] = std::min(ex[min_mark], gradp[c][idof-1]);
          }
        }
      }
    }
  }
}

void
DG::lim()
// *****************************************************************************
// Compute limiter function
// *****************************************************************************
{
  auto d = Disc();
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
  const auto ncomp = m_u.nprop() / rdof;
  const auto nprim = m_p.nprop() / rdof;

  // Combine own and communicated contributions to nodal extrema
  for (const auto& [gid,g] : m_uNodalExtrmc) {
    auto bid = tk::cref_find( d->Bid(), gid );
    for (ncomp_t c=0; c<ncomp; ++c)
    {
      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
      {
        auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
        auto min_mark = max_mark + 1;
        m_uNodalExtrm[bid][max_mark] =
          std::max(g[max_mark], m_uNodalExtrm[bid][max_mark]);
        m_uNodalExtrm[bid][min_mark] =
          std::min(g[min_mark], m_uNodalExtrm[bid][min_mark]);
      }
    }
  }
  for (const auto& [gid,g] : m_pNodalExtrmc) {
    auto bid = tk::cref_find( d->Bid(), gid );
    for (ncomp_t c=0; c<nprim; ++c)
    {
      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
      {
        auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
        auto min_mark = max_mark + 1;
        m_pNodalExtrm[bid][max_mark] =
          std::max(g[max_mark], m_pNodalExtrm[bid][max_mark]);
        m_pNodalExtrm[bid][min_mark] =
          std::min(g[min_mark], m_pNodalExtrm[bid][min_mark]);
      }
    }
  }

  // clear gradients receive buffer
  tk::destroy(m_uNodalExtrmc);
  tk::destroy(m_pNodalExtrmc);

  if (rdof > 1)
    for (const auto& eq : g_dgpde)
      eq.limit( d->T(), m_geoFace, m_geoElem, m_fd, m_esup, m_inpoel, m_coord,
                m_ndof, d->Gid(), d->Bid(), m_uNodalExtrm, m_pNodalExtrm, m_u,
                m_p, m_shockmarker );

  // Send limited solution to neighboring chares
  if (m_sendGhost.empty())
    comlim_complete();
  else
    for(const auto& [cid, ghostdata] : m_sendGhost) {
      std::vector< std::size_t > tetid( ghostdata.size() );
      std::vector< std::vector< tk::real > > u( ghostdata.size() ),
                                             prim( ghostdata.size() );
      std::vector< std::size_t > ndof;
      std::size_t j = 0;
      for(const auto& i : ghostdata) {
        Assert( i < m_fd.Esuel().size()/4, "Sending limiter ghost data" );
        tetid[j] = i;
        u[j] = m_u[i];
        prim[j] = m_p[i];
        if (pref && m_stage == 0) ndof.push_back( m_ndof[i] );
        ++j;
      }
      thisProxy[ cid ].comlim( thisIndex, tetid, u, prim, ndof );
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
    Assert( b < m_uc[2].size(), "Indexing out of bounds" );
    Assert( b < m_pc[2].size(), "Indexing out of bounds" );
    m_uc[2][b] = u[i];
    m_pc[2][b] = prim[i];
    if (pref && m_stage == 0) {
      Assert( b < m_ndofc[2].size(), "Indexing out of bounds" );
      m_ndofc[2][b] = ndof[i];
    }
  }

  // if we have received all solution ghost contributions from neighboring
  // chares (chares we communicate along chare-boundary faces with), and
  // contributed our solution to these neighbors, proceed to limiting
  if (++m_nlim == m_sendGhost.size()) {
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
    Assert( m_uc[2][b.second].size() == m_u.nprop(), "ncomp size mismatch" );
    Assert( m_pc[2][b.second].size() == m_p.nprop(), "ncomp size mismatch" );
    for (std::size_t c=0; c<m_u.nprop(); ++c) {
      m_u(b.first,c,0) = m_uc[2][b.second][c];
    }
    for (std::size_t c=0; c<m_p.nprop(); ++c) {
      m_p(b.first,c,0) = m_pc[2][b.second][c];
    }
    if (pref && m_stage == 0) {
      m_ndof[ b.first ] = m_ndofc[2][ b.second ];
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
          eq.dt( m_coord, m_inpoel, m_fd, m_geoFace, m_geoElem, m_ndof,
            m_u, m_p, m_fd.Esuel().size()/4 );
        if (eqdt < mindt) mindt = eqdt;
      }

      mindt *= g_inputdeck.get< tag::discr, tag::cfl >();
    }
  }
  else
  {
    mindt = d->Dt();
  }

  // Resize the buffer vector of nodal extrema
  resizeNodalExtremac();

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
  thisProxy[ thisIndex ].wait4reco();
  thisProxy[ thisIndex ].wait4nodalExtrema();
  thisProxy[ thisIndex ].wait4lim();
  thisProxy[ thisIndex ].wait4nod();

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
    eq.rhs( d->T(), m_geoFace, m_geoElem, m_fd, m_inpoel, m_boxelems, m_coord,
            m_u, m_p, m_ndof, m_rhs );

  // Explicit time-stepping using RK3 to discretize time-derivative
  for(std::size_t e=0; e<m_nunk; ++e)
    for(std::size_t c=0; c<neq; ++c)
    {
      for (std::size_t k=0; k<m_numEqDof[c]; ++k)
      {
        auto rmark = c*rdof+k;
        auto mark = c*ndof+k;
        m_u(e, rmark, 0) =  rkcoef[0][m_stage] * m_un(e, rmark, 0)
          + rkcoef[1][m_stage] * ( m_u(e, rmark, 0)
            + d->Dt() * m_rhs(e, mark, 0)/m_lhs(e, mark, 0) );
        if(fabs(m_u(e, rmark, 0)) < 1e-16)
          m_u(e, rmark, 0) = 0;
      }
      // zero out unused/reconstructed dofs of equations using reduced dofs
      // (see DGMultiMat::numEquationDofs())
      if (m_numEqDof[c] < rdof) {
        for (std::size_t k=m_numEqDof[c]; k<rdof; ++k)
        {
          auto rmark = c*rdof+k;
          m_u(e, rmark, 0) = 0.0;
        }
      }
    }

  // Update primitives based on the evolved solution
  for (const auto& eq : g_dgpde)
  {
    eq.updateInterfaceCells( m_u, m_fd.Esuel().size()/4, m_ndof );
    eq.updatePrimitives( m_u, m_lhs, m_geoElem, m_p, m_fd.Esuel().size()/4 );
    eq.cleanTraceMaterial( m_geoElem, m_u, m_p, m_fd.Esuel().size()/4 );
  }

  if (m_stage < 2) {

    // continue with next time step stage
    stage();

  } else {

    // Compute diagnostics, e.g., residuals
    auto diag_computed = m_diag.compute( *d, m_u.nunk()-m_fd.Esuel().size()/4,
                                         m_geoElem, m_ndof, m_u );

    // Increase number of iterations and physical time
    d->next();

    // Continue to mesh refinement (if configured)
    if (!diag_computed) refine( std::vector< tk::real >( m_u.nprop(), 0.0 ) );

  }
}

void
DG::refine( [[maybe_unused]] const std::vector< tk::real >& l2res )
// *****************************************************************************
// Optionally refine/derefine mesh
//! \param[in] l2res L2-norms of the residual for each scalar component
//!   computed across the whole problem
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
  const std::set< std::size_t >& /*removedNodes*/,
  const tk::NodeCommMap& nodeCommMap,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& /* bnode */,
  const std::vector< std::size_t >& triinpoel )
// *****************************************************************************
//  Receive new mesh from Refiner
//! \param[in] chunk New mesh chunk (connectivity and global<->local id maps)
//! \param[in] coord New mesh node coordinates
//! \param[in] addedTets Newly added mesh cells and their parents (local ids)
//! \param[in] nodeCommMap New node communication map
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
  [[maybe_unused]] auto old_nelem = m_inpoel.size()/4;

  // Resize mesh data structures
  d->resizePostAMR( chunk, coord, nodeCommMap );

  // Update state
  m_inpoel = d->Inpoel();
  m_coord = d->Coord();
  auto nelem = m_inpoel.size()/4;
  m_p.resize( nelem );
  m_u.resize( nelem );
  m_un.resize( nelem );
  m_lhs.resize( nelem );
  m_rhs.resize( nelem );
  m_uNodalExtrm.resize( Disc()->Bid().size(), std::vector<tk::real>( 2*
    m_ndof_NodalExtrm*g_inputdeck.get< tag::component >().nprop() ) );
  m_pNodalExtrm.resize( Disc()->Bid().size(), std::vector<tk::real>( 2*
    m_ndof_NodalExtrm*m_p.nprop()/g_inputdeck.get< tag::discr, tag::rdof >()));

  // Resize the buffer vector of nodal extrema
  resizeNodalExtremac();

  m_fd = FaceData( m_inpoel, bface, tk::remap(triinpoel,d->Lid()) );

  m_geoFace =
    tk::Fields( tk::genGeoFaceTri( m_fd.Nipfac(), m_fd.Inpofa(), coord ) );
  m_geoElem = tk::Fields( tk::genGeoElemTet( m_inpoel, coord ) );

  m_nfac = m_fd.Inpofa().size()/3;
  m_nunk = nelem;
  m_npoin = coord[0].size();
  m_bndFace.clear();
  m_exptGhost.clear();
  m_sendGhost.clear();
  m_ghost.clear();
  m_esup.clear();

  // Update solution on new mesh, P0 (cell center value) only for now
  m_un = m_u;
  auto pn = m_p;
  auto unprop = m_u.nprop();
  auto pnprop = m_p.nprop();
  for (const auto& [child,parent] : addedTets) {
    Assert( child < nelem, "Indexing out of new solution vector" );
    Assert( parent < old_nelem, "Indexing out of old solution vector" );
    for (std::size_t i=0; i<unprop; ++i) m_u(child,i,0) = m_un(parent,i,0);
    for (std::size_t i=0; i<pnprop; ++i) m_p(child,i,0) = pn(parent,i,0);
  }
  m_un = m_u;

  // Enable SDAG wait for setting up chare boundary faces
  thisProxy[ thisIndex ].wait4fac();

  // Resize communication buffers
  resizeComm();
}

bool
DG::fieldOutput() const
// *****************************************************************************
// Decide wether to output field data
//! \return True if field data is output in this step
// *****************************************************************************
{
  auto d = Disc();

  // Output field data
  return d->fielditer() or d->fieldtime() or d->fieldrange() or d->finished();
}

bool
DG::refinedOutput() const
// *****************************************************************************
// Decide if we write field output using a refined mesh
//! \return True if field output will use a refined mesh
// *****************************************************************************
{
  return g_inputdeck.get< tag::cmd, tag::io, tag::refined >() &&
         g_inputdeck.get< tag::discr, tag::scheme >() != ctr::SchemeType::DG;
}

void
DG::writeFields( CkCallback c )
// *****************************************************************************
// Output mesh field data
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  auto d = Disc();

  // Output time history
  if (d->histiter() or d->histtime() or d->histrange()) {
    std::vector< std::vector< tk::real > > hist;
    for (const auto& eq : g_dgpde) {
      auto h = eq.histOutput( d->Hist(), m_inpoel, m_coord, m_u, m_p );
      hist.insert( end(hist), begin(h), end(h) );
    }
    d->history( std::move(hist) );
  }

  const auto& inpoel = std::get< 0 >( m_outmesh.chunk );
  auto esup = tk::genEsup( inpoel, 4 );

  // Combine own and communicated contributions and finish averaging of node
  // field output in chare boundary nodes
  const auto& lid = std::get< 2 >( m_outmesh.chunk );
  for (const auto& [g,f] : m_nodefieldsc) {
    Assert( m_nodefields.size() == f.first.size(), "Size mismatch" );
    auto p = tk::cref_find( lid, g );
    for (std::size_t i=0; i<f.first.size(); ++i) {
      m_nodefields[i][p] += f.first[i];
      m_nodefields[i][p] /= static_cast< tk::real >(
                             esup.second[p+1] - esup.second[p] + f.second );
    }
  }
  tk::destroy( m_nodefieldsc );

  // Lambda to decide if a node (global id) is on a chare boundary of the field
  // output mesh. p - global node id, return true if node is on the chare
  // boundary.
  auto chbnd = [ this ]( std::size_t p ) {
    return
      std::any_of( m_outmesh.nodeCommMap.cbegin(), m_outmesh.nodeCommMap.cend(),
        [&](const auto& s) { return s.second.find(p) != s.second.cend(); } );
  };

  // Finish computing node field output averages in internal nodes
  auto npoin = m_outmesh.coord[0].size();
  auto& gid = std::get< 1 >( m_outmesh.chunk );
  for (std::size_t p=0; p<npoin; ++p) {
    if (!chbnd(gid[p])) {
      auto n = static_cast< tk::real >( esup.second[p+1] - esup.second[p] );
      for (auto& f : m_nodefields) f[p] /= n;
    }
  }

  // Query fields names requested by user
  auto elemfieldnames = numericFieldNames( tk::Centering::ELEM );
  auto nodefieldnames = numericFieldNames( tk::Centering::NODE );

  // Collect field output names for analytical solutions
  for (const auto& eq : g_dgpde) {
    analyticFieldNames( eq, tk::Centering::ELEM, elemfieldnames );
    analyticFieldNames( eq, tk::Centering::NODE, nodefieldnames );
  }

  if (g_inputdeck.get< tag::pref, tag::pref >())
    elemfieldnames.push_back( "NDOF" );

  elemfieldnames.push_back( "shock_marker" );

  Assert( elemfieldnames.size() == m_elemfields.size(), "Size mismatch" );
  Assert( nodefieldnames.size() == m_nodefields.size(), "Size mismatch" );

  // Output chare mesh and fields metadata to file
  const auto& triinpoel = m_outmesh.triinpoel;
  d->write( inpoel, m_outmesh.coord, m_outmesh.bface, {},
            tk::remap( triinpoel, lid ), elemfieldnames, nodefieldnames,
            {}, m_elemfields, m_nodefields, {}, c );
}

void
DG::comnodeout( const std::vector< std::size_t >& gid,
                const std::vector< std::size_t >& nesup,
                const std::vector< std::vector< tk::real > >& L )
// *****************************************************************************
//  Receive chare-boundary nodal solution (for field output) contributions from
//  neighboring chares
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] nesup Number of elements surrounding points
//! \param[in] L Partial contributions of node fields to chare-boundary nodes
// *****************************************************************************
{
  Assert( gid.size() == nesup.size(), "Size mismatch" );
  for (std::size_t f=0; f<L.size(); ++f)
    Assert( gid.size() == L[f].size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i) {
    auto& nf = m_nodefieldsc[ gid[i] ];
    nf.first.resize( L.size() );
    for (std::size_t f=0; f<L.size(); ++f) nf.first[f] += L[f][i];
    nf.second += nesup[i];
  }

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_nnod == Disc()->NodeCommMap().size()) {
    m_nnod = 0;
    comnodeout_complete();
  }
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
  // otherwise prepare for nodal field output
  if (m_stage < 3)
    next();
  else
    startFieldOutput( CkCallback(CkIndex_DG::step(), thisProxy[thisIndex]) );
}

void
DG::evalLB( int nrestart )
// *****************************************************************************
// Evaluate whether to do load balancing
//! \param[in] nrestart Number of times restarted
// *****************************************************************************
{
  auto d = Disc();

  // Detect if just returned from a checkpoint and if so, zero timers
  d->restarted( nrestart );

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
  const auto benchmark = g_inputdeck.get< tag::cmd, tag::benchmark >();

  if ( !benchmark && (d->It()) % rsfreq == 0 ) {

    std::vector< std::size_t > meshdata{ /* finished = */ 0, d->MeshId() };
    contribute( meshdata, CkReduction::nop,
      CkCallback(CkReductionTarget(Transporter,checkpoint), d->Tr()) );

  } else {

    evalLB( /* nrestart = */ -1 );

  }
}

void
DG::step()
// *****************************************************************************
// Evaluate wether to continue with next time step
// *****************************************************************************
{
  // Free memory storing output mesh
  m_outmesh.destroy();

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

    auto meshid = d->MeshId();
    d->contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
                   CkCallback(CkReductionTarget(Transporter,finish), d->Tr()) );

  }
}

#include "NoWarning/dg.def.h"
