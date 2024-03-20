// *****************************************************************************
/*!
  \file      src/Inciter/Ghosts.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Definitions file for generating ghost data structures
  \details   Definitions file for asynchronous distributed
             ghost data structures using Charm++.
*/
// *****************************************************************************

#include "Ghosts.hpp"
#include "DerivedData.hpp"
#include "Reorder.hpp"
#include "Around.hpp"
#include "ChareStateCollector.hpp"

extern tk::CProxy_ChareStateCollector stateProxy;

using inciter::Ghosts;

Ghosts::Ghosts( const CProxy_Discretization& disc,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel,
  std::size_t nunk,
  CkCallback cbDone ) :
  m_disc( disc ),
  m_nunk( nunk ),
  m_inpoel( Disc()->Inpoel() ),
  m_coord( Disc()->Coord() ),
  m_fd( m_inpoel, bface, tk::remap(triinpoel,Disc()->Lid()) ),
  m_geoFace( tk::genGeoFaceTri( m_fd.Nipfac(), m_fd.Inpofa(), m_coord) ),
  m_geoElem( tk::genGeoElemTet( m_inpoel, m_coord ) ),
  m_nfac( m_fd.Inpofa().size()/3 ),
  m_bndFace(),
  m_sendGhost(),
  m_ghost(),
  m_exptGhost(),
  m_bid(),
  m_esup(),
  m_initial( 1 ),
  m_ncomfac( 0 ),
  m_nadj( 0 ),
  m_ncomEsup( 0 ),
  m_ipface(),
  m_ghostData(),
  m_ghostReq( 0 ),
  m_expChBndFace(),
  m_infaces(),
  m_esupc(),
  m_cbAfterDone( cbDone )
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
//! \param[in] nunk Number of unknowns
//! \param[in] cbDone Function to continue with when Ghosts have been computed
// *****************************************************************************
{
  if (g_inputdeck.get< newtag::cmd, tag::chare >() ||
      g_inputdeck.get< newtag::cmd, tag::quiescence >())
    stateProxy.ckLocalBranch()->insert( "Ghosts", thisIndex, CkMyPe(), Disc()->It(),
                                        "Ghosts" );
}

void
Ghosts::startCommSetup()
// *****************************************************************************
//  Start setup of communication maps for cell-centered schemes
// *****************************************************************************
{
  // Ensure that mesh partition is not leaky
  Assert( !tk::leakyPartition(m_fd.Esuel(), m_inpoel, m_coord),
    "Input mesh to Ghosts leaky" );

  // Ensure mesh physical boundary for the entire problem not leaky,
  // effectively checking if the user has specified boundary conditions on all
  // physical boundary faces
  bndIntegral();
}

void
Ghosts::bndIntegral()
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
    s[0] += m_geoFace(f,0) * m_geoFace(f,1);
    s[1] += m_geoFace(f,0) * m_geoFace(f,2);
    s[2] += m_geoFace(f,0) * m_geoFace(f,3);
  }

  s.push_back( 1.0 );  // positive: call-back to resizeComm() after reduction
  s.push_back( static_cast< tk::real >( Disc()->MeshId() ) );

  // Send contribution to host summing partial surface integrals
  contribute( s, CkReduction::sum_double,
    CkCallback(CkReductionTarget(Transporter,bndint), Disc()->Tr()) );
}

void
Ghosts::resizeComm()
// *****************************************************************************
//  Start sizing communication buffers and setting up ghost data
// *****************************************************************************
{
  // Enable SDAG wait for setting up chare boundary faces
  thisProxy[ thisIndex ].wait4fac();

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

  if ( g_inputdeck.get< newtag::cmd, tag::feedback >() ) d->Tr().chbndface();

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

void
Ghosts::comfac( int fromch, const tk::UnsMesh::FaceSet& infaces )
// *****************************************************************************
//  Receive unique set of faces we potentially share with/from another chare
//! \param[in] fromch Sender chare id
//! \param[in] infaces Unique set of faces we potentially share with fromch
// *****************************************************************************
{
  if (g_inputdeck.get< newtag::cmd, tag::chare >() ||
      g_inputdeck.get< newtag::cmd, tag::quiescence >())
    stateProxy.ckLocalBranch()->insert( "Ghosts", thisIndex, CkMyPe(), Disc()->It(),
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
Ghosts::bndFaces()
// *****************************************************************************
// Compute chare-boundary faces
//! \details This is called when both send and receives are completed on a
//!  chare and thus we are ready to compute chare-boundary faces and ghost data.
// *****************************************************************************
{
  auto d = Disc();
  if ( g_inputdeck.get< newtag::cmd, tag::feedback >() ) d->Tr().chcomfac();
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
        auto c = findchare(t);
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

void
Ghosts::setupGhost()
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
        auto c = findchare(t);
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
Ghosts::reqGhost()
// *****************************************************************************
// Receive requests for ghost data
// *****************************************************************************
{
  if (g_inputdeck.get< newtag::cmd, tag::chare >() ||
      g_inputdeck.get< newtag::cmd, tag::quiescence >())
    stateProxy.ckLocalBranch()->insert( "Ghosts", thisIndex, CkMyPe(), Disc()->It(),
                                        "reqGhost" );

  // If every chare we communicate with has requested ghost data from us, we may
  // fulfill the requests, but only if we have already setup our ghost data.
  if (++m_ghostReq == Disc()->NodeCommMap().size()) {
    m_ghostReq = 0;
    reqghost_complete();
  }
}

void
Ghosts::sendGhost()
// *****************************************************************************
// Send all of our ghost data to fellow chares
// *****************************************************************************
{
  if (g_inputdeck.get< newtag::cmd, tag::chare >() ||
      g_inputdeck.get< newtag::cmd, tag::quiescence >())
    stateProxy.ckLocalBranch()->insert( "Ghosts", thisIndex, CkMyPe(), Disc()->It(),
                                        "sendGhost" );

  for (const auto& c : m_ghostData)
    thisProxy[ c.first ].comGhost( thisIndex, c.second );

  if ( g_inputdeck.get< newtag::cmd, tag::feedback >() ) Disc()->Tr().chghost();
}

void
Ghosts::comGhost( int fromch, const GhostData& ghost )
// *****************************************************************************
// Receive ghost data on chare boundaries from fellow chare
//! \param[in] fromch Caller chare ID
//! \param[in] ghost Ghost data, see Inciter/FaceData.h for the type
// *****************************************************************************
{
  if (g_inputdeck.get< newtag::cmd, tag::chare >() ||
      g_inputdeck.get< newtag::cmd, tag::quiescence >())
    stateProxy.ckLocalBranch()->insert( "Ghosts", thisIndex, CkMyPe(), Disc()->It(),
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
      addGeoFace(t, id);
      // add node-triplet to node-face connectivity
      inpofa[3*id[0]+0] = tk::cref_find( lid, t[2] );
      inpofa[3*id[0]+1] = tk::cref_find( lid, t[1] );
      inpofa[3*id[0]+2] = tk::cref_find( lid, t[0] );

      // if ghost tet id not yet encountered on boundary with fromch
      auto i = ghostelem.find( e );
      if (i != end(ghostelem)) {
        // fill in elements surrounding face
        addEsuf(id, i->second);
        // fill in elements surrounding element
        addEsuel(id, i->second, t);
      } else {
        // fill in elements surrounding face
        addEsuf(id, m_nunk);
        // fill in elements surrounding element
        addEsuel(id, m_nunk, t);
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
      Assert( nodetripletMatch(id, t) == 3,
        "Mismatch/Overmatch in inpoel and inpofa at chare-boundary face" );
    }
  }

  // Signal the runtime system that all workers have received their
  // face-adjacency
  if (++m_nadj == m_ghostData.size()) faceAdj();
}

void
Ghosts::faceAdj()
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
  Assert( faceMatch(), "Chare-boundary element-face "
    "connectivity (esuf) does not match" );

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
Ghosts::nodeNeighSetup()
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
Ghosts::comEsup( int fromch,
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
Ghosts::adj()
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

  if ( g_inputdeck.get< newtag::cmd, tag::feedback >() ) Disc()->Tr().chadj();

  // Error checking on ghost data
  for(const auto& n : m_sendGhost)
    for([[maybe_unused]] const auto& i : n.second)
      Assert( i < m_fd.Esuel().size()/4, "Sender contains ghost tet id. ");

  // Create a mapping between local ghost tet ids and zero-based boundary ids
  std::vector< std::size_t > c( tk::sumvalsize( m_ghost ) );
  std::size_t j = 0;
  for (const auto& n : m_ghost) {
    for(const auto& i : n.second) {
      c[j++] = i.second;
    }
  }
  m_bid = tk::assignLid( c );

  // Basic error checking on ghost tet ID map
  Assert( m_ghost.find( thisIndex ) == m_ghost.cend(),
          "Ghost id map should not contain data for own chare ID" );

  // Store expected ghost tet IDs
  for (const auto& n : m_ghost)
    for ([[maybe_unused]] const auto& g : n.second)
      Assert( m_exptGhost.insert( g.second ).second,
              "Failed to store local tetid as exptected ghost id" );

  // Callback function from DG/FV after ghost-setup is done
  m_cbAfterDone.send();
}

bool
Ghosts::leakyAdjacency()
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
    s[0] += m_geoFace(f,0) * m_geoFace(f,1);
    s[1] += m_geoFace(f,0) * m_geoFace(f,2);
    s[2] += m_geoFace(f,0) * m_geoFace(f,3);
  }

  // chare-boundary faces
  for (std::size_t f=m_fd.Nipfac(); f<m_fd.Esuf().size()/2; ++f) {
    s[0] += m_geoFace(f,0) * m_geoFace(f,1);
    s[1] += m_geoFace(f,0) * m_geoFace(f,2);
    s[2] += m_geoFace(f,0) * m_geoFace(f,3);
  }

  auto eps = std::numeric_limits< tk::real >::epsilon() * 100;
  return std::abs(s[0]) > eps || std::abs(s[1]) > eps || std::abs(s[2]) > eps;
}

bool
Ghosts::faceMatch()
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

bool
Ghosts::receivedChBndFaces()
// *****************************************************************************
// Verify that all chare-boundary faces have been received
//! \return True if all chare-boundary faces have been received
// *****************************************************************************
{
  const auto& lid = Disc()->Lid();
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
       auto A = tk::cref_find( lid, f[0] );
       auto B = tk::cref_find( lid, f[1] );
       auto C = tk::cref_find( lid, f[2] );
       msg << '{' << A << ',' << B << ',' << C << "}:("
           << x[A] << ',' << y[A] << ',' << z[A] << ' '
           << x[B] << ',' << y[B] << ',' << z[B] << ' '
           << x[C] << ',' << y[C] << ',' << z[C] << ") ";
     }

  tk::destroy( m_expChBndFace );

  // Error out with info on missing faces
  auto s = msg.str();
  if (!s.empty()) {
    Throw( "Ghosts chare " + std::to_string(thisIndex) +
           " missing face(s) {local node ids} (node coords): " + s );
  } else {
    return true;
  }
}

int
Ghosts::findchare( const tk::UnsMesh::Face& t )
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

std::size_t
Ghosts::nodetripletMatch(
  const std::array< std::size_t, 2 >& id,
  const tk::UnsMesh::Face& t )
// *****************************************************************************
// Check if entries in inpoel, inpofa and node-triplet are consistent
//! \param[in] id Local face and (inner) tet id adjacent to it
//! \param[in] t node-triplet associated with the chare boundary face
//! \return number of nodes in inpoel that matched with t and inpofa
// *****************************************************************************
{
  const auto& esuf = m_fd.Esuf();
  const auto& inpofa = m_fd.Inpofa();
  const auto& lid = Disc()->Lid();

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
Ghosts::addEsuf(
  const std::array< std::size_t, 2 >& id,
  std::size_t ghostid )
// *****************************************************************************
// Fill elements surrounding a face along chare boundary
//! \param[in] id Local face and (inner) tet id adjacent to it
//! \param[in] ghostid Local ID for ghost tet
//! \details This function extends and fills in the elements surrounding faces
//!   data structure (esuf) so that the left and right element id is filled
//!   in correctly on chare boundaries to contain the correct inner tet id and
//!   the local tet id for the outer (ghost) tet, both adjacent to the given
//!   chare-face boundary. Prior to this function, this data structure does not
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
Ghosts::addEsuel(
  const std::array< std::size_t, 2 >& id,
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
//    the outer ghost tet in Ghosts::comGhost in place of the -1 before.
// *****************************************************************************
{
  const auto& lid = Disc()->Lid();
  [[maybe_unused]] const auto& esuf = m_fd.Esuf();
  auto& esuel = m_fd.Esuel();

  std::array< tk::UnsMesh::Face, 4 > face;
  for (std::size_t f = 0; f<4; ++f)
    for (std::size_t i = 0; i<3; ++i)
      face[f][i] = m_inpoel[ id[1]*4 + tk::lpofa[f][i] ];

  tk::UnsMesh::Face tl{{ tk::cref_find( lid, t[0] ),
                         tk::cref_find( lid, t[1] ),
                         tk::cref_find( lid, t[2] ) }};

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
Ghosts::addGeoFace(
  const tk::UnsMesh::Face& t,
  const std::array< std::size_t, 2 >& id )
// *****************************************************************************
// Fill face-geometry data along chare boundary
//! \param[in] t Face (given by 3 global node IDs) on the chare boundary
//! \param[in] id Local face and (inner) tet id adjacent to face t
//! \details This function fills in the face geometry data along a chare
//!    boundary.
// *****************************************************************************
{
  const auto& lid = Disc()->Lid();

  // get global node IDs reversing order to get outward-pointing normal
  auto A = tk::cref_find( lid, t[2] );
  auto B = tk::cref_find( lid, t[1] );
  auto C = tk::cref_find( lid, t[0] );
  auto geochf = tk::geoFaceTri( {{m_coord[0][A], m_coord[0][B], m_coord[0][C]}},
                                {{m_coord[1][A], m_coord[1][B], m_coord[1][C]}},
                                {{m_coord[2][A], m_coord[2][B], m_coord[2][C]}} );

  for (std::size_t i=0; i<7; ++i)
    m_geoFace(id[0],i) = geochf(0,i);
}

#include "NoWarning/ghosts.def.h"
