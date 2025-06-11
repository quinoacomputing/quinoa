// *****************************************************************************
/*!
  \file      src/Inciter/Discretization.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \details   Data and functionality common to all discretization schemes
  \see       Discretization.h and Discretization.C for more info.
*/
// *****************************************************************************

#include "Tags.hpp"
#include "Reorder.hpp"
#include "Vector.hpp"
#include "DerivedData.hpp"
#include "Discretization.hpp"
#include "MeshWriter.hpp"
#include "DiagWriter.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Inciter/Options/Scheme.hpp"
#include "Print.hpp"
#include "Around.hpp"
#include "QuinoaBuildConfig.hpp"
#include "ConjugateGradients.hpp"
#include "ALE.hpp"

#include "M2MTransfer.hpp"

namespace inciter {

static CkReduction::reducerType PDFMerger;
extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;

} // inciter::

using inciter::Discretization;

Discretization::Discretization(
  std::size_t meshid,
  const std::vector< CProxy_Discretization >& disc,
  const CProxy_ALE& aleproxy,
  const tk::CProxy_ConjugateGradients& conjugategradientsproxy,
  const CProxy_Transporter& transporter,
  const tk::CProxy_MeshWriter& meshwriter,
  const tk::UnsMesh::CoordMap& coordmap,
  const tk::UnsMesh::Chunk& el,
  const tk::CommMaps& msum,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel,
  const std::unordered_map< std::size_t, std::set< std::size_t > >& elemblockid,
  int nc ) :
  m_meshid( meshid ),
  m_transfer( g_inputdeck.get< tag::transfer >() ),
  m_disc( disc ),
  m_nchare( nc ),
  m_it( 0 ),
  m_itr( 0 ),
  m_itf( 0 ),
  m_initial( 1 ),
  m_t( g_inputdeck.get< tag::t0 >() ),
  m_lastDumpTime( -std::numeric_limits< tk::real >::max() ),
  m_physFieldFloor( 0.0 ),
  m_physHistFloor( 0.0 ),
  m_rangeFieldFloor( 0.0 ),
  m_rangeHistFloor( 0.0 ),
  m_dt( g_inputdeck.get< tag::dt >() ),
  m_dtn( m_dt ),
  m_nvol( 0 ),
  m_nxfer( 0 ),
  m_ale( aleproxy ),
  m_transporter( transporter ),
  m_meshwriter( meshwriter ),
  m_el( el ),     // fills m_inpoel, m_gid, m_lid
  m_coord( setCoord( coordmap ) ),
  m_coordn( m_coord ),
  m_nodeCommMap(),
  m_edgeCommMap(),
  m_meshvol( 0.0 ),
  m_v( m_gid.size(), 0.0 ),
  m_vol( m_gid.size(), 0.0 ),
  m_volc(),
  m_voln( m_vol ),
  m_vol0( m_inpoel.size()/4, 0.0 ),
  m_bid(),
  m_timer(),
  m_refined( 0 ),
  m_prevstatus( std::chrono::high_resolution_clock::now() ),
  m_nrestart( 0 ),
  m_histdata(),
  m_nsrc( 0 ),
  m_ndst( 0 ),
  m_meshvel( 0, 3 ),
  m_meshvel_converged( true ),
  m_bface( bface ),
  m_triinpoel( triinpoel ),
  m_elemblockid( elemblockid )
// *****************************************************************************
//  Constructor
//! \param[in] meshid Mesh ID
//! \param[in] disc All Discretization proxies (one per mesh)
//! \param[in] aleproxy Distributed ALE proxy
//! \param[in] conjugategradientsproxy Distributed Conjugrate Gradients linear
//!   solver proxy
//! \param[in] transporter Host (Transporter) proxy
//! \param[in] meshwriter Mesh writer proxy
//! \param[in] coordmap Coordinates of mesh nodes and their global IDs
//! \param[in] el Elements of the mesh chunk we operate on
//! \param[in] msum Communication maps associated to chare IDs bordering the
//!   mesh chunk we operate on
//! \param[in] bface Face lists mapped to side set ids
//! \param[in] triinpoel Interconnectivity of points and boundary-faces
//! \param[in] elemblockid Local tet ids associated with mesh block ids
//! \param[in] nc Total number of Discretization chares
// *****************************************************************************
{
  Assert( !m_inpoel.empty(), "No elements assigned to Discretization chare" );
  Assert( tk::positiveJacobians( m_inpoel, m_coord ),
          "Jacobian in input mesh to Discretization non-positive" );
  #if not defined(__INTEL_COMPILER) || defined(NDEBUG)
  // The above ifdef skips running the conformity test with the intel compiler
  // in debug mode only. This is necessary because in tk::conforming(), filling
  // up the map can fail with some meshes (only in parallel), e.g., tube.exo,
  // used by some regression tests, due to the intel compiler generating some
  // garbage incorrect code - only in debug, only in parallel, only with that
  // mesh.
  Assert( tk::conforming( m_inpoel, m_coord ),
          "Input mesh to Discretization not conforming" );
  #endif

  // Store communication maps
  for (const auto& [ c, maps ] : msum) {
    m_nodeCommMap[c] = maps.get< tag::node >();
    m_edgeCommMap[c] = maps.get< tag::edge >();
  }

  // Get ready for computing/communicating nodal volumes
  startvol();

  // Get chare-boundary node-id map
  m_bid = genBid();

  // Find host elements of user-specified points where time histories are
  // saved, and save the shape functions evaluated at the point locations
  const auto& pt = g_inputdeck.get< tag::history_output, tag::point >();
  for (std::size_t p=0; p<pt.size(); ++p) {
    std::array< tk::real, 4 > N;
    const auto& l = pt[p].get< tag::coord >();
    const auto& id = pt[p].get< tag::id >();
    for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
      if (tk::intet( m_coord, m_inpoel, l, e, N )) {
        m_histdata.push_back( HistData{{ id, e, {l[0],l[1],l[2]}, N }} );
        break;
      }
    }
  }

  // Insert ConjugrateGradients solver chare array element if needed
  if (g_inputdeck.get< tag::ale, tag::ale >()) {
    m_ale[ thisIndex ].insert( conjugategradientsproxy,
                               m_coord, m_inpoel,
                               m_gid, m_lid, m_nodeCommMap );
  } else {
    m_meshvel.resize( m_gid.size() );
  }

  // Register mesh with mesh-transfer lib
  if (m_disc.size() == 1 || m_transfer.empty()) {
    // skip transfer if single mesh or if not involved in coupling
    transferInit();
  } else {
    if (thisIndex == 0) {
      exam2m::addMesh( thisProxy, m_nchare,
        CkCallback( CkIndex_Discretization::transferInit(), thisProxy ) );
      //std::cout << "Disc: " << m_meshid << " m2m::addMesh()\n";
    }
  }
}

std::unordered_map< std::size_t, std::size_t >
Discretization::genBid()
// *****************************************************************************
// Generate the Bid data-structure based on the node communication-map
// *****************************************************************************
{
  // Count the number of mesh nodes at which we receive data from other chares
  // and compute map associating boundary-chare node ID to global node ID
  std::vector< std::size_t > c( tk::sumvalsize( m_nodeCommMap ) );
  std::size_t j = 0;
  for (const auto& [ch,n] : m_nodeCommMap) for (auto i : n) c[j++] = i;
  tk::unique( c );
  return tk::assignLid( c );
}

void
Discretization::transferInit()
// *****************************************************************************
// Our mesh has been registered with the mesh-to-mesh transfer library (if
// coupled to other solver)
// *****************************************************************************
{
  // Compute number of mesh points owned
  std::size_t npoin = m_gid.size();
  for (auto g : m_gid) if (tk::slave(m_nodeCommMap,g,thisIndex)) --npoin;

  // Tell the RTS that the Discretization chares have been created and compute
  // the total number of mesh points across the distributed mesh
  std::vector< std::size_t > meshdata{ m_meshid, npoin };
  contribute( meshdata, CkReduction::sum_ulong,
    CkCallback( CkReductionTarget(Transporter,disccreated), m_transporter ) );
}

void
Discretization::meshvelStart(
  const tk::UnsMesh::Coords vel,
  const std::vector< tk::real >& soundspeed,
  const std::unordered_map< int,
    std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& bnorm,
  tk::real adt,
  CkCallback done ) const
// *****************************************************************************
// Start computing new mesh velocity for ALE mesh motion
//! \param[in] vel Fluid velocity at mesh nodes
//! \param[in] soundspeed Speed of sound at mesh nodes
//! \param[in] bnorm Face normals in boundary points associated to side sets
//! \param[in] adt alpha*dt of the RK time step
//! \param[in] done Function to continue with when mesh velocity has been
//!   computed
// *****************************************************************************
{
  if (g_inputdeck.get< tag::ale, tag::ale >())
    m_ale[ thisIndex ].ckLocal()->start( vel, soundspeed, done,
      m_coord, m_coordn, m_vol0, m_vol, bnorm, m_initial, m_it, m_t, adt );
  else
    done.send();
}

const tk::Fields&
Discretization::meshvel() const
// *****************************************************************************
//! Query the mesh velocity
//! \return Mesh velocity
// *****************************************************************************
{
  if (g_inputdeck.get< tag::ale, tag::ale >())
    return m_ale[ thisIndex ].ckLocal()->meshvel();
  else
    return m_meshvel;
}

void
Discretization::meshvelBnd(
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& bnode,
  const std::vector< std::size_t >& triinpoel) const
// *****************************************************************************
// Query ALE mesh velocity boundary condition node lists and node lists at
// which ALE moves boundaries
// *****************************************************************************
{
  if (g_inputdeck.get< tag::ale, tag::ale >())
    m_ale[ thisIndex ].ckLocal()->meshvelBnd( bface, bnode, triinpoel );
}

void
Discretization::meshvelConv()
// *****************************************************************************
//! Assess and record mesh velocity linear solver convergence
// *****************************************************************************
{
  auto smoother = g_inputdeck.get< tag::ale, tag::smoother >();

  if (g_inputdeck.get< tag::ale, tag::ale >() &&
      (smoother == ctr::MeshVelocitySmootherType::LAPLACE or
       smoother == ctr::MeshVelocitySmootherType::HELMHOLTZ))
  {
    m_meshvel_converged &= m_ale[ thisIndex ].ckLocal()->converged();
  }
}

void
Discretization::comfinal()
// *****************************************************************************
// Finish setting up communication maps and solution transfer callbacks
// *****************************************************************************
{
  // Generate own subset of solver/mesh transfer list
  for (const auto& t : m_transfer) {
    if (t.src == m_meshid || t.dst == m_meshid) {
      m_mytransfer.push_back( t );
    }
  }

  // Signal the runtime system that the workers have been created
  std::vector< std::size_t > meshdata{ /* initial = */ 1, m_meshid };
  contribute( meshdata, CkReduction::sum_ulong,
    CkCallback(CkReductionTarget(Transporter,comfinal), m_transporter) );
}

void
Discretization::transfer(
  tk::Fields& u,
  std::size_t dirn,
  CkCallback cb )
// *****************************************************************************
//  Start solution transfer (if coupled)
//! \param[in,out] u Solution to transfer from/to
//! \param[in] dirn Direction of solution transfer. 0: from background to
//!   overset, 1: from overset to background
//! \param[in] cb Callback to call when back and forth transfers complete.
//! \details This function initiates the solution transfer (direction dependent
//!   on 'dirn') between meshes. It invokes a reduction to Transporter when the
//!   transfer in one direction is complete (dirn == 0), or calls back the
//!   'cb' function in Scheme when transfers both directions are complete.
//!   The function relies on 'dirn' to make this decision.
// *****************************************************************************
{
  if (m_mytransfer.empty()) {   // skip transfer if not involved in coupling

    cb.send();

  } else {

    m_transfer_complete = cb;

    // determine source and destination mesh depending on direction of transfer
    std::size_t fromMesh(0), toMesh(0);
    CkCallback cb_xfer;
    if (dirn == 0) {
      fromMesh = m_mytransfer[m_nsrc].src;
      toMesh = m_mytransfer[m_ndst].dst;
      cb_xfer = CkCallback( CkIndex_Discretization::to_complete(), thisProxy[thisIndex] );
    }
    else {
      fromMesh = m_mytransfer[m_nsrc].dst;
      toMesh = m_mytransfer[m_ndst].src;
      cb_xfer = CkCallback( CkIndex_Discretization::from_complete(), thisProxy[thisIndex] );
    }

    // Pass source and destination meshes to mesh transfer lib (if coupled)
    Assert( m_nsrc < m_mytransfer.size(), "Indexing out of mytransfer[src]" );
    if (fromMesh == m_meshid) {
      exam2m::setSourceTets( thisProxy, thisIndex, &m_inpoel, &m_coord, u );
      ++m_nsrc;
    } else {
      m_nsrc = 0;
    }
    Assert( m_ndst < m_mytransfer.size(), "Indexing out of mytransfer[dst]" );
    if (toMesh == m_meshid) {
      exam2m::setDestPoints( thisProxy, thisIndex, &m_coord, u,
        cb_xfer );
      ++m_ndst;
    } else {
      m_ndst = 0;
    }

  }

  m_nsrc = 0;
  m_ndst = 0;
}

void Discretization::to_complete()
// *****************************************************************************
//! Solution transfer from background to overset mesh completed (from ExaM2M)
//! \brief This is called by ExaM2M on the destination mesh when the
//!   transfer completes. Since this is called only on the destination, we find
//!   and notify the corresponding source of the completion.
// *****************************************************************************
{
  // Lookup the source disc and notify it of completion
  for (auto& t : m_transfer) {
    if (m_meshid == t.dst) {
      m_disc[ t.src ][ thisIndex ].transfer_complete();
    }
  }

  thisProxy[ thisIndex ].transfer_complete();
}

void Discretization::from_complete()
// *****************************************************************************
//! Solution transfer from overset to background mesh completed (from ExaM2M)
//! \brief This is called by ExaM2M on the destination mesh when the
//!   transfer completes. Since this is called only on the destination, we find
//!   and notify the corresponding source of the completion.
// *****************************************************************************
{
  // Lookup the source disc and notify it of completion
  for (auto& t : m_transfer) {
    if (m_meshid == t.src) {
      m_disc[ t.dst ][ thisIndex ].transfer_complete_from_dest();
    }
  }

  m_transfer_complete.send();
}

void Discretization::transfer_complete_from_dest()
// *****************************************************************************
//! Solution transfer completed (from dest Discretization)
//! \details Called on the source only by the destination when a back and forth
//!   transfer step completes.
// *****************************************************************************
{
  m_transfer_complete.send();
}

void Discretization::transfer_complete()
// *****************************************************************************
//! Solution transfer completed (one-way)
//! \note Single exit point after solution transfer between meshes
// *****************************************************************************
{
  contribute( sizeof(nullptr), nullptr, CkReduction::nop,
    CkCallback(CkReductionTarget(Transporter,solutionTransferred),
    m_transporter) );
}

std::vector< std::size_t >
Discretization::bndel() const
// *****************************************************************************
// Find elements along our mesh chunk boundary
//! \return List of local element ids that have at least a single node
//!   contributing to a chare boundary
// *****************************************************************************
{
  // Lambda to find out if a mesh node is shared with another chare
  auto shared = [this]( std::size_t i ){
    for (const auto& [c,n] : m_nodeCommMap)
      if (n.find(i) != end(n)) return true;
    return false;
  };

  // Find elements along our mesh chunk boundary
  std::vector< std::size_t > e;
  for (std::size_t n=0; n<m_inpoel.size(); ++n)
    if (shared( m_gid[ m_inpoel[n] ] )) e.push_back( n/4 );
  tk::unique( e );

  return e;
}

void
Discretization::resizePostAMR(
  const tk::UnsMesh::Chunk& chunk,
  const tk::UnsMesh::Coords& coord,
  const std::unordered_map< std::size_t, std::size_t >& /*amrNodeMap*/,
  const tk::NodeCommMap& nodeCommMap,
  const std::set< std::size_t >& /*removedNodes*/,
  const std::unordered_map< std::size_t, std::set< std::size_t > >& elemblockid )
// *****************************************************************************
//  Resize mesh data structures after mesh refinement
//! \param[in] chunk New mesh chunk (connectivity and global<->local id maps)
//! \param[in] coord New mesh node coordinates
//! \param[in] amrNodeMap Node id map after amr (local ids)
//! \param[in] nodeCommMap New node communication map
//! \param[in] removedNodes Newly removed mesh node local ids
//! \param[in] elemblockid New local tet ids associated with mesh block ids
// *****************************************************************************
{
  m_el = chunk;         // updates m_inpoel, m_gid, m_lid
  m_nodeCommMap.clear();
  m_nodeCommMap = nodeCommMap;        // update node communication map
  m_elemblockid.clear();
  m_elemblockid = elemblockid;

  // Update mesh volume container size
  m_vol.resize( m_gid.size(), 0.0 );
  if (!m_voln.empty()) m_voln.resize( m_gid.size(), 0.0 );

  // Regenerate bid data
  tk::destroy(m_bid);
  m_bid = genBid();

  // update mesh node coordinates
  m_coord = coord;

  // we are no longer during setup
  m_initial = 0;
}

void
Discretization::startvol()
// *****************************************************************************
//  Get ready for (re-)computing/communicating nodal volumes
// *****************************************************************************
{
  m_nvol = 0;
  thisProxy[ thisIndex ].wait4vol();

  // Zero out mesh volume container
  std::fill( begin(m_vol), end(m_vol), 0.0 );

  // Clear receive buffer that will be used for collecting nodal volumes
  m_volc.clear();
}

void
Discretization::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types
//!  \details Since this is a [initnode] routine, see the .ci file, the
//!   Charm++ runtime system executes the routine exactly once on every
//!   logical node early on in the Charm++ init sequence. Must be static as
//!   it is called without an object. See also: Section "Initializations at
//!   Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  PDFMerger = CkReduction::addReducer( tk::mergeUniPDFs );
}

tk::UnsMesh::Coords
Discretization::setCoord( const tk::UnsMesh::CoordMap& coordmap )
// *****************************************************************************
// Set mesh coordinates based on coordinates map
// *****************************************************************************
{
  Assert( coordmap.size() == m_gid.size(), "Size mismatch" );
  Assert( coordmap.size() == m_lid.size(), "Size mismatch" );

  tk::UnsMesh::Coords coord;
  coord[0].resize( coordmap.size() );
  coord[1].resize( coordmap.size() );
  coord[2].resize( coordmap.size() );

  for (const auto& [ gid, coords ] : coordmap) {
    auto i = tk::cref_find( m_lid, gid );
    coord[0][i] = coords[0];
    coord[1][i] = coords[1];
    coord[2][i] = coords[2];
  }

  return coord;
}

void
Discretization::remap(
  const std::unordered_map< std::size_t, std::size_t >& map )
// *****************************************************************************
//  Remap mesh data based on new local ids
//! \param[in] map Mapping of old->new local ids
// *****************************************************************************
{
  // Remap connectivity containing local IDs
  for (auto& l : m_inpoel) l = tk::cref_find(map,l);

  // Remap global->local id map
  for (auto& [g,l] : m_lid) l = tk::cref_find(map,l);

  // Remap global->local id map
  auto maxid = std::numeric_limits< std::size_t >::max();
  std::vector< std::size_t > newgid( m_gid.size(), maxid );
  for (const auto& [o,n] : map) newgid[n] = m_gid[o];
  m_gid = std::move( newgid );

  Assert( std::all_of( m_gid.cbegin(), m_gid.cend(),
            [=](std::size_t i){ return i < maxid; } ),
          "Not all gid have been remapped" );

  // Remap nodal volumes (with contributions along chare-boundaries)
  std::vector< tk::real > newvol( m_vol.size(), 0.0 );
  for (const auto& [o,n] : map) newvol[n] = m_vol[o];
  m_vol = std::move( newvol );

  // Remap nodal volumes (without contributions along chare-boundaries)
  std::vector< tk::real > newv( m_v.size(), 0.0 );
  for (const auto& [o,n] : map) newv[n] = m_v[o];
  m_v = std::move( newv );

  // Remap locations of node coordinates
  tk::UnsMesh::Coords newcoord;
  auto npoin = m_coord[0].size();
  newcoord[0].resize( npoin );
  newcoord[1].resize( npoin );
  newcoord[2].resize( npoin );
  for (const auto& [o,n] : map) {
    newcoord[0][n] = m_coord[0][o];
    newcoord[1][n] = m_coord[1][o];
    newcoord[2][n] = m_coord[2][o];
  }
  m_coord = std::move( newcoord );
}

void
Discretization::setRefiner( const CProxy_Refiner& ref )
// *****************************************************************************
//  Set Refiner Charm++ proxy
//! \param[in] ref Incoming refiner proxy to store
// *****************************************************************************
{
  m_refiner = ref;
}

void
Discretization::vol()
// *****************************************************************************
// Sum mesh volumes to nodes, start communicating them on chare-boundaries
// *****************************************************************************
{
  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  // Compute nodal volumes on our chunk of the mesh
  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ m_inpoel[e*4+0], m_inpoel[e*4+1],
                                           m_inpoel[e*4+2], m_inpoel[e*4+3] }};
    // compute element Jacobi determinant * 5/120 = element volume / 4
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da ) * 5.0 / 120.0;
    ErrChk( J > 0, "Element Jacobian non-positive: PE:" +
                   std::to_string(CkMyPe()) + ", node IDs: " +
                   std::to_string(m_gid[N[0]]) + ',' +
                   std::to_string(m_gid[N[1]]) + ',' +
                   std::to_string(m_gid[N[2]]) + ',' +
                   std::to_string(m_gid[N[3]]) + ", coords: (" +
                   std::to_string(x[N[0]]) + ", " +
                   std::to_string(y[N[0]]) + ", " +
                   std::to_string(z[N[0]]) + "), (" +
                   std::to_string(x[N[1]]) + ", " +
                   std::to_string(y[N[1]]) + ", " +
                   std::to_string(z[N[1]]) + "), (" +
                   std::to_string(x[N[2]]) + ", " +
                   std::to_string(y[N[2]]) + ", " +
                   std::to_string(z[N[2]]) + "), (" +
                   std::to_string(x[N[3]]) + ", " +
                   std::to_string(y[N[3]]) + ", " +
                   std::to_string(z[N[3]]) + ')' );
    // scatter add V/4 to nodes
    for (std::size_t j=0; j<4; ++j) m_vol[N[j]] += J;

    // save element volumes at t=t0
    if (m_it == 0) m_vol0[e] = J * 4.0;
  }

  // Store nodal volumes without contributions from other chares on
  // chare-boundaries
  m_v = m_vol;

  // Send our nodal volume contributions to neighbor chares
  if (m_nodeCommMap.empty())
   totalvol();
  else
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< tk::real > v( n.size() );
      std::size_t j = 0;
      for (auto i : n) v[ j++ ] = m_vol[ tk::cref_find(m_lid,i) ];
      thisProxy[c].comvol( std::vector<std::size_t>(begin(n), end(n)), v );
    }

  ownvol_complete();
}

void
Discretization::comvol( const std::vector< std::size_t >& gid,
                        const std::vector< tk::real >& nodevol )
// *****************************************************************************
//  Receive nodal volumes on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive volume contributions
//! \param[in] nodevol Partial sums of nodal volume contributions to
//!    chare-boundary nodes
//! \details This function receives contributions to m_vol, which stores the
//!   nodal volumes. While m_vol stores own contributions, m_volc collects the
//!   neighbor chare contributions during communication. This way work on m_vol
//!   and m_volc is overlapped. The contributions are applied in totalvol().
// *****************************************************************************
{
  Assert( nodevol.size() == gid.size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i)
    m_volc[ gid[i] ] += nodevol[i];

  if (++m_nvol == m_nodeCommMap.size()) {
    m_nvol = 0;
    comvol_complete();
  }
}

void
Discretization::totalvol()
// *****************************************************************************
// Sum mesh volumes and contribute own mesh volume to total volume
// *****************************************************************************
{
  // Add received contributions to nodal volumes
  for (const auto& [gid, vol] : m_volc)
    m_vol[ tk::cref_find(m_lid,gid) ] += vol;

  // Clear receive buffer
  tk::destroy(m_volc);

  // Sum mesh volume to host
  std::vector< tk::real > tvol{ 0.0,
                                static_cast<tk::real>(m_initial),
                                static_cast<tk::real>(m_meshid) };
  for (auto v : m_v) tvol[0] += v;
  contribute( tvol, CkReduction::sum_double,
    CkCallback(CkReductionTarget(Transporter,totalvol), m_transporter) );
}

void
Discretization::stat( tk::real mesh_volume )
// *****************************************************************************
// Compute mesh cell statistics
//! \param[in] mesh_volume Total mesh volume
// *****************************************************************************
{
  // Store total mesh volume
  m_meshvol = mesh_volume;

  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  auto MIN = -std::numeric_limits< tk::real >::max();
  auto MAX = std::numeric_limits< tk::real >::max();
  std::vector< tk::real > min{ MAX, MAX, MAX };
  std::vector< tk::real > max{ MIN, MIN, MIN };
  std::vector< tk::real > sum{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  tk::UniPDF edgePDF( 1e-4 );
  tk::UniPDF volPDF( 1e-4 );
  tk::UniPDF ntetPDF( 1e-4 );

  // Compute points surrounding points
  auto psup = tk::genPsup( m_inpoel, 4, tk::genEsup(m_inpoel,4) );
  Assert( psup.second.size()-1 == m_gid.size(),
          "Number of mesh points and number of global IDs unequal" );

  // Compute edge length statistics
  // Note that while the min and max edge lengths are independent of the number
  // of CPUs (by the time they are aggregated across all chares), the sum of
  // the edge lengths and the edge length PDF are not. This is because the
  // edges on the chare-boundary are counted multiple times and we
  // conscientiously do not make an effort to precisely compute this, because
  // that would require communication and more complex logic. Since these
  // statistics are intended as simple average diagnostics, we ignore these
  // small differences. For reproducible average edge lengths and edge length
  // PDFs, run the mesh in serial.
  for (std::size_t p=0; p<m_gid.size(); ++p)
    for (auto i : tk::Around(psup,p)) {
       const auto dx = x[ i ] - x[ p ];
       const auto dy = y[ i ] - y[ p ];
       const auto dz = z[ i ] - z[ p ];
       const auto length = std::sqrt( dx*dx + dy*dy + dz*dz );
       if (length < min[0]) min[0] = length;
       if (length > max[0]) max[0] = length;
       sum[0] += 1.0;
       sum[1] += length;
       edgePDF.add( length );
    }

  // Compute mesh cell volume statistics
  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ m_inpoel[e*4+0], m_inpoel[e*4+1],
                                           m_inpoel[e*4+2], m_inpoel[e*4+3] }};
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto L = std::cbrt( tk::triple( ba, ca, da ) / 6.0 );
    if (L < min[1]) min[1] = L;
    if (L > max[1]) max[1] = L;
    sum[2] += 1.0;
    sum[3] += L;
    volPDF.add( L );
  }

  // Contribute stats of number of tetrahedra (ntets)
  sum[4] = 1.0;
  min[2] = max[2] = sum[5] = static_cast< tk::real >( m_inpoel.size() / 4 );
  ntetPDF.add( min[2] );

  min.push_back( static_cast<tk::real>(m_meshid) );
  max.push_back( static_cast<tk::real>(m_meshid) );
  sum.push_back( static_cast<tk::real>(m_meshid) );

  // Contribute to mesh statistics across all Discretization chares
  contribute( min, CkReduction::min_double,
    CkCallback(CkReductionTarget(Transporter,minstat), m_transporter) );
  contribute( max, CkReduction::max_double,
    CkCallback(CkReductionTarget(Transporter,maxstat), m_transporter) );
  contribute( sum, CkReduction::sum_double,
    CkCallback(CkReductionTarget(Transporter,sumstat), m_transporter) );

  // Serialize PDFs to raw stream
  auto stream = tk::serialize( m_meshid, { edgePDF, volPDF, ntetPDF } );
  // Create Charm++ callback function for reduction of PDFs with
  // Transporter::pdfstat() as the final target where the results will appear.
  CkCallback cb( CkIndex_Transporter::pdfstat(nullptr), m_transporter );
  // Contribute serialized PDF of partial sums to host via Charm++ reduction
  contribute( stream.first, stream.second.get(), PDFMerger, cb );
}

void
Discretization::boxvol(
  const std::vector< std::unordered_set< std::size_t > >& nodes,
  const std::unordered_map< std::size_t, std::set< std::size_t > >& nodeblk,
  std::size_t nuserblk )
// *****************************************************************************
// Compute total box IC volume
//! \param[in] nodes Node list contributing to box IC volume (for each IC box)
//! \param[in] nodeblk Node list associated to mesh blocks contributing to block
//!   volumes (for each IC box)
//! \param[in] nuserblk Number of user IC mesh blocks
// *****************************************************************************
{
  // Compute partial box IC volume (just add up all boxes)
  tk::real boxvol = 0.0;
  for (const auto& b : nodes) for (auto i : b) boxvol += m_v[i];

  // Compute partial IC mesh block volume
  std::vector< tk::real > blockvols;
  if (nuserblk > 0) {
    blockvols.resize(nuserblk,0.0);
    for (const auto& [blid, ndset] : nodeblk) {
      // The following if-test makes sure we access volumes only of mesh blocks
      // with user-specified ICs
      if (blid < nuserblk) {
        for (const auto& n : ndset) blockvols[blid] += m_v[n];
      }
    }
  }

  // Sum up box IC volume across all chares
  auto meshdata = blockvols;
  meshdata.push_back(boxvol);
  meshdata.push_back(static_cast<tk::real>(m_meshid));
  contribute( meshdata, CkReduction::sum_double,
    CkCallback(CkReductionTarget(Transporter,boxvol), m_transporter) );
}

void
Discretization::write(
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& bnode,
  const std::vector< std::size_t >& triinpoel,
  const std::vector< std::string>& elemfieldnames,
  const std::vector< std::string>& nodefieldnames,
  const std::vector< std::string>& elemsurfnames,
  const std::vector< std::string>& nodesurfnames,
  const std::vector< std::vector< tk::real > >& elemfields,
  const std::vector< std::vector< tk::real > >& nodefields,
  const std::vector< std::vector< tk::real > >& elemsurfs,
  const std::vector< std::vector< tk::real > >& nodesurfs,
  CkCallback c )
// *****************************************************************************
//  Output mesh and fields data (solution dump) to file(s)
//! \param[in] inpoel Mesh connectivity for the mesh chunk to be written
//! \param[in] coord Node coordinates of the mesh chunk to be written
//! \param[in] bface Map of boundary-face lists mapped to corresponding side set
//!   ids for this mesh chunk
//! \param[in] bnode Map of boundary-node lists mapped to corresponding side set
//!   ids for this mesh chunk
//! \param[in] triinpoel Interconnectivity of points and boundary-face in this
//!   mesh chunk
//! \param[in] elemfieldnames Names of element fields to be output to file
//! \param[in] nodefieldnames Names of node fields to be output to file
//! \param[in] elemsurfnames Names of elemental surface fields to be output to
//!   file
//! \param[in] nodesurfnames Names of node surface fields to be output to file
//! \param[in] elemfields Field data in mesh elements to output to file
//! \param[in] nodefields Field data in mesh nodes to output to file
//! \param[in] elemsurfs Surface field data in mesh elements to output to file
//! \param[in] nodesurfs Surface field data in mesh nodes to output to file
//! \param[in] c Function to continue with after the write
//! \details Since m_meshwriter is a Charm++ chare group, it never migrates and
//!   an instance is guaranteed on every PE. We index the first PE on every
//!   logical compute node. In Charm++'s non-SMP mode, a node is the same as a
//!   PE, so the index is the same as CkMyPe(). In SMP mode the index is the
//!   first PE on every logical node. In non-SMP mode this yields one or more
//!   output files per PE with zero or non-zero virtualization, respectively. If
//!   there are multiple chares on a PE, the writes are serialized per PE, since
//!   only a single entry method call can be executed at any given time. In SMP
//!   mode, still the same number of files are output (one per chare), but the
//!   output is serialized through the first PE of each compute node. In SMP
//!   mode, channeling multiple files via a single PE on each node is required
//!   by NetCDF and HDF5, as well as ExodusII, since none of these libraries are
//!   thread-safe.
// *****************************************************************************
{
  // If the previous iteration refined (or moved) the mesh or this is called
  // before the first time step, we also output the mesh.
  bool meshoutput = m_itf == 0 ? true : false;

  auto eps = std::numeric_limits< tk::real >::epsilon();
  bool fieldoutput = false;

  // Output field data only if there is no dump at this physical time yet
  if (std::abs(m_lastDumpTime - m_t) > eps ) {
    m_lastDumpTime = m_t;
    ++m_itf;
    fieldoutput = true;
  }

  // set of sidesets where fieldoutput is required
  std::set< int > outsets;
  const auto& osv = g_inputdeck.get< tag::field_output, tag::sideset >();
  outsets.insert(osv.begin(), osv.end());

  m_meshwriter[ CkNodeFirst( CkMyNode() ) ].
    write( m_meshid, meshoutput, fieldoutput, m_itr, m_itf, m_t, thisIndex,
           g_inputdeck.get< tag::cmd, tag::io, tag::output >(),
           inpoel, coord, bface, bnode, triinpoel, elemfieldnames,
           nodefieldnames, elemsurfnames, nodesurfnames, elemfields, nodefields,
           elemsurfs, nodesurfs, outsets, c );
}

void
Discretization::setdt( tk::real newdt )
// *****************************************************************************
// Set time step size
//! \param[in] newdt Size of the new time step
// *****************************************************************************
{
  m_dtn = m_dt;
  m_dt = newdt;

  // Truncate the size of last time step
  const auto term = g_inputdeck.get< tag::term >();
  if (m_t+m_dt > term) m_dt = term - m_t;
}

void
Discretization::next()
// *****************************************************************************
// Prepare for next step
// *****************************************************************************
{
  // Update floor of physics time divided by output interval times
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto ft = g_inputdeck.get< tag::field_output, tag::time_interval >();
  if (ft > eps) m_physFieldFloor = std::floor( m_t / ft );
  const auto ht = g_inputdeck.get< tag::history_output, tag::time_interval >();
  if (ht > eps) m_physHistFloor = std::floor( m_t / ht );

  // Update floors of physics time divided by output interval times for ranges
  const auto& rf = g_inputdeck.get< tag::field_output, tag::time_range >();
  if (!rf.empty()) {
    if (m_t > rf[0] and m_t < rf[1])
      m_rangeFieldFloor = std::floor( m_t / rf[2] );
  }
  const auto& rh = g_inputdeck.get< tag::history_output, tag::time_range >();
  if (!rh.empty()) {
    if (m_t > rh[0] and m_t < rh[1])
      m_rangeHistFloor = std::floor( m_t / rh[2] );
  }

  ++m_it;
  m_t += m_dt;
}

void
Discretization::grindZero()
// *****************************************************************************
//  Zero grind-time
// *****************************************************************************
{
  m_prevstatus = std::chrono::high_resolution_clock::now();

  if (thisIndex == 0 && m_meshid == 0) {
    const auto verbose = g_inputdeck.get< tag::cmd, tag::verbose >();
    const auto& def =
      g_inputdeck_defaults.get< tag::cmd, tag::io, tag::screen >();
    tk::Print print( g_inputdeck.get< tag::cmd >().logname( def, m_nrestart ),
                     verbose ? std::cout : std::clog,
                     std::ios_base::app );
    print.diag( "Starting time stepping ..." );
  }
}

bool
Discretization::restarted( int nrestart )
// *****************************************************************************
//  Detect if just returned from a checkpoint and if so, zero timers
//! \param[in] nrestart Number of times restarted
//! \return True if restart detected
// *****************************************************************************
{
  // Detect if just restarted from checkpoint:
  //   nrestart == -1 if there was no checkpoint this step
  //   d->Nrestart() == nrestart if there was a checkpoint this step
  //   if both false, just restarted from a checkpoint
  bool restarted = nrestart != -1 and m_nrestart != nrestart;

   // If just restarted from checkpoint
  if (restarted) {
    // Update number of restarts
    m_nrestart = nrestart;
    // Start timer measuring time stepping wall clock time
    m_timer.zero();
    // Zero grind-timer
    grindZero();
  }

  return restarted;
}

std::string
Discretization::histfilename( const std::string& id,
                              std::streamsize precision )
// *****************************************************************************
//  Construct history output filename
//! \param[in] id History point id
//! \param[in] precision Floating point precision to use for output
//! \return History file name
// *****************************************************************************
{
  auto of = g_inputdeck.get< tag::cmd, tag::io, tag::output >();
  std::stringstream ss;

  auto mid =
    m_disc.size() > 1 ? std::string( '.' + std::to_string(m_meshid) ) : "";
  ss << std::setprecision(static_cast<int>(precision)) << of << mid << ".hist." << id;

  return ss.str();
}

void
Discretization::histheader( std::vector< std::string >&& names )
// *****************************************************************************
//  Output headers for time history files (one for each point)
//! \param[in] names History output variable names
// *****************************************************************************
{
  for (const auto& h : m_histdata) {
    auto prec = g_inputdeck.get< tag::history_output, tag::precision >();
    tk::DiagWriter hw( histfilename( h.get< tag::id >(), prec ),
                       g_inputdeck.get< tag::history_output, tag::format >(),
                       prec );
    hw.header( names );
  }
}

void
Discretization::history( std::vector< std::vector< tk::real > >&& data )
// *****************************************************************************
//  Output time history for a time step
//! \param[in] data Time history data for all variables and equations integrated
// *****************************************************************************
{
  Assert( data.size() == m_histdata.size(), "Size mismatch" );

  std::size_t i = 0;
  for (const auto& h : m_histdata) {
    auto prec = g_inputdeck.get< tag::history_output, tag::precision >();
    tk::DiagWriter hw( histfilename( h.get< tag::id >(), prec ),
                       g_inputdeck.get< tag::history_output, tag::format >(),
                       prec,
                       std::ios_base::app );
    hw.diag( m_it, m_t, m_dt, data[i] );
    ++i;
  }
}

bool
Discretization::fielditer() const
// *****************************************************************************
//  Decide if field output iteration count interval is hit
//! \return True if field output iteration count interval is hit
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::benchmark >()) return false;

  return m_it % g_inputdeck.get< tag::field_output, tag::interval >() == 0;
}

bool
Discretization::fieldtime() const
// *****************************************************************************
//  Decide if field output physics time interval is hit
//! \return True if field output physics time interval is hit
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::benchmark >()) return false;

  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto ft = g_inputdeck.get< tag::field_output, tag::time_interval >();

  if (ft < eps) return false;

  return std::floor(m_t/ft) - m_physFieldFloor > eps;
}

bool
Discretization::fieldrange() const
// *****************************************************************************
//  Decide if physics time falls into a field output time range
//! \return True if physics time falls into a field output time range
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::benchmark >()) return false;

  const auto eps = std::numeric_limits< tk::real >::epsilon();

  bool output = false;

  const auto& rf = g_inputdeck.get< tag::field_output, tag::time_range >();
  if (!rf.empty()) {
    if (m_t > rf[0] and m_t < rf[1])
      output |= std::floor(m_t/rf[2]) - m_rangeFieldFloor > eps;
  }

  return output;
}

bool
Discretization::histiter() const
// *****************************************************************************
//  Decide if history output iteration count interval is hit
//! \return True if history output iteration count interval is hit
// *****************************************************************************
{
  const auto hist = g_inputdeck.get< tag::history_output, tag::interval >();
  const auto& hist_points = g_inputdeck.get< tag::history_output, tag::point >();

  return m_it % hist == 0 and not hist_points.empty();
}

bool
Discretization::histtime() const
// *****************************************************************************
//  Decide if history output physics time interval is hit
//! \return True if history output physics time interval is hit
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::benchmark >()) return false;

  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto ht = g_inputdeck.get< tag::history_output, tag::time_interval >();

  if (ht < eps) return false;

  return std::floor(m_t/ht) - m_physHistFloor > eps;
}

bool
Discretization::histrange() const
// *****************************************************************************
//  Decide if physics time falls into a history output time range
//! \return True if physics time falls into a history output time range
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::benchmark >()) return false;

  const auto eps = std::numeric_limits< tk::real >::epsilon();

  bool output = false;

  const auto& rh = g_inputdeck.get< tag::history_output, tag::time_range >();
  if (!rh.empty()) {
    if (m_t > rh[0] and m_t < rh[1])
      output |= std::floor(m_t/rh[2]) - m_rangeHistFloor > eps;
  }

  return output;
}

bool
Discretization::finished() const
// *****************************************************************************
//  Decide if this is the last time step
//! \return True if this is the last time step
// *****************************************************************************
{
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto nstep = g_inputdeck.get< tag::nstep >();
  const auto term = g_inputdeck.get< tag::term >();

  return std::abs(m_t-term) < eps or m_it >= nstep;
}

void
Discretization::status()
// *****************************************************************************
// Output one-liner status report
// *****************************************************************************
{
  // Query after how many time steps user wants TTY dump
  const auto tty = g_inputdeck.get< tag::ttyi >();

  // estimate grind time (taken between this and the previous time step)
  using std::chrono::duration_cast;
  using ms = std::chrono::milliseconds;
  using clock = std::chrono::high_resolution_clock;
  auto grind_time = duration_cast< ms >(clock::now() - m_prevstatus).count();
  m_prevstatus = clock::now();

  if (thisIndex==0 and m_meshid == 0 and not (m_it%tty)) {

    const auto term = g_inputdeck.get< tag::term >();
    const auto t0 = g_inputdeck.get< tag::t0 >();
    const auto nstep = g_inputdeck.get< tag::nstep >();
    const auto diag = g_inputdeck.get< tag::diagnostics, tag::interval >();
    const auto lbfreq = g_inputdeck.get< tag::cmd, tag::lbfreq >();
    const auto rsfreq = g_inputdeck.get< tag::cmd, tag::rsfreq >();
    const auto verbose = g_inputdeck.get< tag::cmd, tag::verbose >();
    const auto benchmark = g_inputdeck.get< tag::cmd, tag::benchmark >();
    const auto steady = g_inputdeck.get< tag::steady_state >();

    // estimate time elapsed and time for accomplishment
    tk::Timer::Watch ete, eta;
    if (not steady) m_timer.eta( term-t0, m_t-t0, nstep, m_it, ete, eta );

    const auto& def =
      g_inputdeck_defaults.get< tag::cmd, tag::io, tag::screen >();
    tk::Print print( g_inputdeck.get< tag::cmd >().logname( def, m_nrestart ),
                     verbose ? std::cout : std::clog,
                     std::ios_base::app );

    // Output one-liner
    print << std::setfill(' ') << std::setw(8) << m_it << "  "
          << std::scientific << std::setprecision(6)
          << std::setw(12) << m_t << "  "
          << m_dt << "  "
          << std::setfill('0')
          << std::setw(3) << ete.hrs.count() << ":"
          << std::setw(2) << ete.min.count() << ":"
          << std::setw(2) << ete.sec.count() << "  "
          << std::setw(3) << eta.hrs.count() << ":"
          << std::setw(2) << eta.min.count() << ":"
          << std::setw(2) << eta.sec.count() << "  "
          << std::scientific << std::setprecision(6) << std::setfill(' ')
          << std::setw(9) << grind_time << "  ";

    // Augment one-liner status with output indicators
    if (fielditer() or fieldtime() or fieldrange()) print << 'f';
    if (not (m_it % diag)) print << 'd';
    if (histiter() or histtime() or histrange()) print << 't';
    if (m_refined) print << 'h';
    if (not (m_it % lbfreq) && not finished()) print << 'l';
    if (not benchmark && (not (m_it % rsfreq) || finished())) print << 'r';

    if (not m_meshvel_converged) print << 'a';
    m_meshvel_converged = true; // get ready for next time step

    print << std::endl;
  }
}

#include "NoWarning/discretization.def.h"
