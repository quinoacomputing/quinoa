// *****************************************************************************
/*!
  \file      src/Inciter/FV.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     FV advances a system of PDEs with the finite volume scheme
  \details   FV advances a system of partial differential equations (PDEs) using
    the finite volume (FV) (on tetrahedron elements).
  \see The documentation in FV.h.
*/
// *****************************************************************************

#include <algorithm>
#include <numeric>
#include <sstream>

#include "FV.hpp"
#include "Discretization.hpp"
#include "FVPDE.hpp"
#include "DiagReducer.hpp"
#include "DerivedData.hpp"
#include "ElemDiagnostics.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Refiner.hpp"
#include "Limiter.hpp"
#include "PrefIndicator.hpp"
#include "Reorder.hpp"
#include "Vector.hpp"
#include "Around.hpp"
#include "Integrate/Basis.hpp"
#include "FieldOutput.hpp"
#include "ChareStateCollector.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< FVPDE > g_fvpde;

//! Runge-Kutta coefficients
static const std::array< std::array< tk::real, 3 >, 2 >
  rkcoef{{ {{ 0.0, 3.0/4.0, 1.0/3.0 }}, {{ 1.0, 1.0/4.0, 2.0/3.0 }} }};

} // inciter::

extern tk::CProxy_ChareStateCollector stateProxy;

using inciter::FV;

FV::FV( const CProxy_Discretization& disc,
        const CProxy_Ghosts& ghostsproxy,
        const std::map< int, std::vector< std::size_t > >& bface,
        const std::map< int, std::vector< std::size_t > >& /* bnode */,
        const std::vector< std::size_t >& triinpoel ) :
  m_disc( disc ),
  m_ghosts( ghostsproxy ),
  m_nsol( 0 ),
  m_ninitsol( 0 ),
  m_nlim( 0 ),
  m_nnod( 0 ),
  m_u( Disc()->Inpoel().size()/4,
       g_inputdeck.get< tag::discr, tag::rdof >()*
       g_inputdeck.get< tag::component >().nprop() ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_p( m_u.nunk(),
       g_inputdeck.get< tag::discr, tag::rdof >()*
         std::accumulate( begin(g_fvpde), end(g_fvpde), 0u,
           [](std::size_t s, const FVPDE& eq){ return s + eq.nprim(); } ) ),
  m_lhs( m_u.nunk(),
         g_inputdeck.get< tag::component >().nprop() ),
  m_rhs( m_u.nunk(), m_lhs.nprop() ),
  m_nunk( m_u.nunk() ),
  m_npoin( Disc()->Coord()[0].size() ),
  m_diag(),
  m_stage( 0 ),
  m_uc(),
  m_pc(),
  m_initial( 1 ),
  m_elemfields(),
  m_nodefields(),
  m_nodefieldsc(),
  m_outmesh(),
  m_boxelems()
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::chare >() ||
      g_inputdeck.get< tag::cmd, tag::quiescence >())
    stateProxy.ckLocalBranch()->insert( "FV", thisIndex, CkMyPe(), Disc()->It(),
                                        "FV" );

  usesAtSync = true;    // enable migration at AtSync

  // Enable SDAG wait for initially building the solution vector and limiting
  if (m_initial) {
    thisProxy[ thisIndex ].wait4sol();
    thisProxy[ thisIndex ].wait4lim();
    thisProxy[ thisIndex ].wait4nod();
  }

  // Ensure that mesh partition is not leaky
  Assert( !tk::leakyPartition(myGhosts()->m_fd.Esuel(), myGhosts()->m_inpoel,
          myGhosts()->m_coord), "Input mesh to FV leaky" );

  m_ghosts[thisIndex].insert(m_disc, bface, triinpoel, m_nunk);

  // global-sync to call doneInserting on m_ghosts
  auto meshid = Disc()->MeshId();
  contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
    CkCallback(CkReductionTarget(Transporter,doneInsertingGhosts),
    Disc()->Tr()) );
}

void
FV::registerReducers()
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
FV::ResumeFromSync()
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
FV::setup()
// *****************************************************************************
// Set initial conditions, generate lhs, output mesh
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::chare >() ||
      g_inputdeck.get< tag::cmd, tag::quiescence >())
    stateProxy.ckLocalBranch()->insert( "FV", thisIndex, CkMyPe(), Disc()->It(),
                                        "setup" );

  // Resize solution vectors, lhs and rhs by the number of ghost tets
  m_nunk = myGhosts()->m_nunk;
  m_u.resize( m_nunk );
  m_un.resize( m_nunk );
  m_p.resize( m_nunk );
  m_lhs.resize( m_nunk );
  m_rhs.resize( m_nunk );

  // Size communication buffer for solution
  for (auto& u : m_uc) u.resize( myGhosts()->m_bid.size() );
  for (auto& p : m_pc) p.resize( myGhosts()->m_bid.size() );

  // Ensure that we also have all the geometry and connectivity data
  // (including those of ghosts)
  Assert( myGhosts()->m_geoElem.nunk() == m_u.nunk(),
    "GeoElem unknowns size mismatch" );

  auto d = Disc();

  // Basic error checking on sizes of element geometry data and connectivity
  Assert( myGhosts()->m_geoElem.nunk() == m_lhs.nunk(),
    "Size mismatch in FV::setup()" );

  // Compute left-hand side of discrete PDEs
  lhs();

  // Determine elements inside user-defined IC box
  for (auto& eq : g_fvpde)
    eq.IcBoxElems( myGhosts()->m_geoElem, myGhosts()->m_fd.Esuel().size()/4,
      m_boxelems );

  // Compute volume of user-defined box IC
  d->boxvol( {} );      // punt for now

  // Query time history field output labels from all PDEs integrated
  const auto& hist_points = g_inputdeck.get< tag::history, tag::point >();
  if (!hist_points.empty()) {
    std::vector< std::string > histnames;
    for (const auto& eq : g_fvpde) {
      auto n = eq.histNames();
      histnames.insert( end(histnames), begin(n), end(n) );
    }
    d->histheader( std::move(histnames) );
  }
}

void
FV::box( tk::real v )
// *****************************************************************************
// Receive total box IC volume and set conditions in box
//! \param[in] v Total volume within user-specified box
// *****************************************************************************
{
  auto d = Disc();

  // Store user-defined box IC volume
  d->Boxvol() = v;

  // Set initial conditions for all PDEs
  for (const auto& eq : g_fvpde)
  {
    eq.initialize( m_lhs, myGhosts()->m_inpoel, myGhosts()->m_coord, m_boxelems,
      m_u, d->T(), myGhosts()->m_fd.Esuel().size()/4 );
    eq.updatePrimitives( m_u, m_p, myGhosts()->m_fd.Esuel().size()/4 );
  }

  m_un = m_u;

  // Output initial conditions to file (regardless of whether it was requested)
  startFieldOutput( CkCallback(CkIndex_FV::start(), thisProxy[thisIndex]) );
}

void
FV::start()
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
FV::startFieldOutput( CkCallback c )
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

      const auto& tr = tk::remap( myGhosts()->m_fd.Triinpoel(), d->Gid() );
      d->Ref()->outref( myGhosts()->m_fd.Bface(), {}, tr, c );

    } else {

      // cut off ghosts from mesh connectivity and coordinates
      const auto& tr = tk::remap( myGhosts()->m_fd.Triinpoel(), d->Gid() );
      extractFieldOutput( {}, d->Chunk(), d->Coord(), {}, {},
                          d->NodeCommMap(), myGhosts()->m_fd.Bface(), {}, tr, c );

    }

  }
}

void
FV::next()
// *****************************************************************************
// Advance equations to next time step
// *****************************************************************************
{
  // communicate solution ghost data (if any)
  if (myGhosts()->m_sendGhost.empty())
    comsol_complete();
  else
    for(const auto& [cid, ghostdata] : myGhosts()->m_sendGhost) {
      std::vector< std::size_t > tetid( ghostdata.size() );
      std::vector< std::vector< tk::real > > u( ghostdata.size() ),
                                             prim( ghostdata.size() );
      std::size_t j = 0;
      for(const auto& i : ghostdata) {
        Assert( i < myGhosts()->m_fd.Esuel().size()/4,
          "Sending solution ghost data" );
        tetid[j] = i;
        u[j] = m_u[i];
        prim[j] = m_p[i];
        ++j;
      }
      thisProxy[ cid ].comsol( thisIndex, tetid, u, prim );
    }

  ownsol_complete();
}

void
FV::comsol( int fromch,
            const std::vector< std::size_t >& tetid,
            const std::vector< std::vector< tk::real > >& u,
            const std::vector< std::vector< tk::real > >& prim )
// *****************************************************************************
//  Receive chare-boundary solution ghost data from neighboring chares
//! \param[in] fromch Sender chare id
//! \param[in] tetid Ghost tet ids we receive solution data for
//! \param[in] u Solution ghost data
//! \param[in] prim Primitive variables in ghost cells
//! \details This function receives contributions to the unlimited solution
//!   from fellow chares.
// *****************************************************************************
{
  Assert( u.size() == tetid.size(), "Size mismatch in FV::comsol()" );
  Assert( prim.size() == tetid.size(), "Size mismatch in FV::comsol()" );

  // Find local-to-ghost tet id map for sender chare
  const auto& n = tk::cref_find( myGhosts()->m_ghost, fromch );

  for (std::size_t i=0; i<tetid.size(); ++i) {
    auto j = tk::cref_find( n, tetid[i] );
    Assert( j >= myGhosts()->m_fd.Esuel().size()/4,
      "Receiving solution non-ghost data" );
    auto b = tk::cref_find( myGhosts()->m_bid, j );
    Assert( b < m_uc[0].size(), "Indexing out of bounds" );
    m_uc[0][b] = u[i];
    m_pc[0][b] = prim[i];
  }

  // if we have received all solution ghost contributions from neighboring
  // chares (chares we communicate along chare-boundary faces with), and
  // contributed our solution to these neighbors, proceed to reconstructions
  if (++m_nsol == myGhosts()->m_sendGhost.size()) {
    m_nsol = 0;
    comsol_complete();
  }
}

void
FV::extractFieldOutput(
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

  // Evaluate element solution on incoming mesh
  auto [ue,pe,un,pn] = evalSolution( inpoel, coord, addedTets );

  // Collect field output from numerical solution requested by user
  m_elemfields = numericFieldOutput( ue, tk::Centering::ELEM, pe );
  m_nodefields = numericFieldOutput( un, tk::Centering::NODE, pn );

  // Collect field output from analytical solutions (if exist)
  auto geoElem = tk::genGeoElemTet( inpoel, coord );
  auto t = Disc()->T();
  for (const auto& eq : g_fvpde) {
    analyticFieldOutput( eq, tk::Centering::ELEM, geoElem.extract(1,0),
      geoElem.extract(2,0), geoElem.extract(3,0), t, m_elemfields );
    analyticFieldOutput( eq, tk::Centering::NODE, coord[0], coord[1], coord[2],
      t, m_nodefields );
  }

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
FV::evalSolution(
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
        auto Bn = tk::eval_basis( 1, dot(J[0],h), dot(J[1],h), dot(J[2],h) );
        auto u = eval_state( uncomp, 0, rdof, 1, e, m_u, Bn, {0, uncomp-1} );
        auto p = eval_state( pncomp, 0, rdof, 1, e, m_p, Bn, {0, pncomp-1} );
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
      auto B = tk::eval_basis( 1, dot(Jp[0],h), dot(Jp[1],h), dot(Jp[2],h) );
      auto u = eval_state( uncomp, 0, rdof, 1, parent, m_u, B, {0, uncomp-1} );
      auto p = eval_state( pncomp, 0, rdof, 1, parent, m_p, B, {0, pncomp-1} );
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
        auto Bn = tk::eval_basis( 1, dot(Jp[0],hn), dot(Jp[1],hn), dot(Jp[2],hn) );
        auto cnu = eval_state(uncomp, 0, rdof, 1, parent, m_u, Bn,
          {0, uncomp-1});
        auto cnp = eval_state(pncomp, 0, rdof, 1, parent, m_p, Bn,
          {0, pncomp-1});
        // Assign child node solution
        for (std::size_t i=0; i<uncomp; ++i) un(inpoel[c4+j],i,0) += cnu[i];
        for (std::size_t i=0; i<pncomp; ++i) pn(inpoel[c4+j],i,0) += cnp[i];
      }
    }
  }

  return { ue, pe, un, pn };
}

void
FV::lhs()
// *****************************************************************************
// Compute left-hand side of discrete transport equations
// *****************************************************************************
{
  for (const auto& eq : g_fvpde) eq.lhs( myGhosts()->m_geoElem, m_lhs );

  if (!m_initial) stage();
}

void
FV::reco()
// *****************************************************************************
// Compute reconstructions
// *****************************************************************************
{
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();

  // Combine own and communicated contributions of unreconstructed solution and
  // degrees of freedom in cells (if p-adaptive)
  for (const auto& b : myGhosts()->m_bid) {
    Assert( m_uc[0][b.second].size() == m_u.nprop(), "ncomp size mismatch" );
    Assert( m_pc[0][b.second].size() == m_p.nprop(), "ncomp size mismatch" );
    for (std::size_t c=0; c<m_u.nprop(); ++c) {
      m_u(b.first,c,0) = m_uc[0][b.second][c];
    }
    for (std::size_t c=0; c<m_p.nprop(); ++c) {
      m_p(b.first,c,0) = m_pc[0][b.second][c];
    }
  }

  if (rdof > 1) {
    // Reconstruct second-order solution and primitive quantities
    for (const auto& eq : g_fvpde)
      eq.reconstruct( myGhosts()->m_geoElem, myGhosts()->m_fd,
        myGhosts()->m_esup, myGhosts()->m_inpoel, myGhosts()->m_coord, m_u, m_p );
  }

  // start limiting
  lim();
}

void
FV::lim()
// *****************************************************************************
// Compute limiter function
// *****************************************************************************
{
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();

  if (rdof > 1)
    for (const auto& eq : g_fvpde)
      eq.limit( myGhosts()->m_geoElem, myGhosts()->m_fd, myGhosts()->m_esup,
        myGhosts()->m_inpoel, myGhosts()->m_coord, m_u, m_p );

  // Send limited solution to neighboring chares
  if (myGhosts()->m_sendGhost.empty())
    comlim_complete();
  else
    for(const auto& [cid, ghostdata] : myGhosts()->m_sendGhost) {
      std::vector< std::size_t > tetid( ghostdata.size() );
      std::vector< std::vector< tk::real > > u( ghostdata.size() ),
                                             prim( ghostdata.size() );
      std::size_t j = 0;
      for(const auto& i : ghostdata) {
        Assert( i < myGhosts()->m_fd.Esuel().size()/4,
          "Sending limiter ghost data" );
        tetid[j] = i;
        u[j] = m_u[i];
        prim[j] = m_p[i];
        ++j;
      }
      thisProxy[ cid ].comlim( thisIndex, tetid, u, prim );
    }

  ownlim_complete();
}

void
FV::comlim( int fromch,
            const std::vector< std::size_t >& tetid,
            const std::vector< std::vector< tk::real > >& u,
            const std::vector< std::vector< tk::real > >& prim )
// *****************************************************************************
//  Receive chare-boundary limiter ghost data from neighboring chares
//! \param[in] fromch Sender chare id
//! \param[in] tetid Ghost tet ids we receive solution data for
//! \param[in] u Limited high-order solution
//! \param[in] prim Limited high-order primitive quantities
//! \details This function receives contributions to the limited solution from
//!   fellow chares.
// *****************************************************************************
{
  Assert( u.size() == tetid.size(), "Size mismatch in FV::comlim()" );
  Assert( prim.size() == tetid.size(), "Size mismatch in FV::comlim()" );

  // Find local-to-ghost tet id map for sender chare
  const auto& n = tk::cref_find( myGhosts()->m_ghost, fromch );

  for (std::size_t i=0; i<tetid.size(); ++i) {
    auto j = tk::cref_find( n, tetid[i] );
    Assert( j >= myGhosts()->m_fd.Esuel().size()/4,
      "Receiving solution non-ghost data" );
    auto b = tk::cref_find( myGhosts()->m_bid, j );
    Assert( b < m_uc[2].size(), "Indexing out of bounds" );
    Assert( b < m_pc[2].size(), "Indexing out of bounds" );
    m_uc[2][b] = u[i];
    m_pc[2][b] = prim[i];
  }

  // if we have received all solution ghost contributions from neighboring
  // chares (chares we communicate along chare-boundary faces with), and
  // contributed our solution to these neighbors, proceed to limiting
  if (++m_nlim == myGhosts()->m_sendGhost.size()) {
    m_nlim = 0;
    comlim_complete();
  }
}

void
FV::dt()
// *****************************************************************************
// Compute time step size
// *****************************************************************************
{
  auto d = Disc();

  // Combine own and communicated contributions of limited solution and degrees
  // of freedom in cells (if p-adaptive)
  for (const auto& b : myGhosts()->m_bid) {
    Assert( m_uc[2][b.second].size() == m_u.nprop(), "ncomp size mismatch" );
    Assert( m_pc[2][b.second].size() == m_p.nprop(), "ncomp size mismatch" );
    for (std::size_t c=0; c<m_u.nprop(); ++c) {
      m_u(b.first,c,0) = m_uc[2][b.second][c];
    }
    for (std::size_t c=0; c<m_p.nprop(); ++c) {
      m_p(b.first,c,0) = m_pc[2][b.second][c];
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
      for (const auto& eq : g_fvpde) {
        auto eqdt =
          eq.dt( myGhosts()->m_fd, myGhosts()->m_geoFace, myGhosts()->m_geoElem,
            m_u, m_p, myGhosts()->m_fd.Esuel().size()/4 );
        if (eqdt < mindt) mindt = eqdt;
      }

      mindt *= g_inputdeck.get< tag::discr, tag::cfl >();
    }
  }
  else
  {
    mindt = d->Dt();
  }

  // Contribute to minimum dt across all chares then advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(FV,solve), thisProxy) );
}

void
FV::solve( tk::real newdt )
// *****************************************************************************
// Compute right-hand side of discrete transport equations
//! \param[in] newdt Size of this new time step
// *****************************************************************************
{
  // Enable SDAG wait for building the solution vector during the next stage
  thisProxy[ thisIndex ].wait4sol();
  thisProxy[ thisIndex ].wait4lim();
  thisProxy[ thisIndex ].wait4nod();

  auto d = Disc();
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
  const auto neq = m_u.nprop()/rdof;

  // Set new time step size
  if (m_stage == 0) d->setdt( newdt );

  // Update Un
  if (m_stage == 0) m_un = m_u;

  for (const auto& eq : g_fvpde)
    eq.rhs( d->T(), myGhosts()->m_geoFace, myGhosts()->m_geoElem,
      myGhosts()->m_fd, myGhosts()->m_inpoel, myGhosts()->m_coord, m_u, m_p,
      m_rhs );

  // Explicit time-stepping using RK3 to discretize time-derivative
  for (std::size_t e=0; e<m_nunk; ++e)
    for (std::size_t c=0; c<neq; ++c)
    {
      auto rmark = c*rdof;
      m_u(e, rmark, 0) =  rkcoef[0][m_stage] * m_un(e, rmark, 0)
        + rkcoef[1][m_stage] * ( m_u(e, rmark, 0)
          + d->Dt() * m_rhs(e, c, 0)/m_lhs(e, c, 0) );
      // zero out reconstructed dofs of equations using reduced dofs
      if (rdof > 1) {
        for (std::size_t k=1; k<rdof; ++k)
        {
          rmark = c*rdof+k;
          m_u(e, rmark, 0) = 0.0;
        }
      }
    }

  // Update primitives based on the evolved solution
  for (const auto& eq : g_fvpde)
  {
    eq.updatePrimitives( m_u, m_p, myGhosts()->m_fd.Esuel().size()/4 );
    eq.cleanTraceMaterial( myGhosts()->m_geoElem, m_u, m_p,
      myGhosts()->m_fd.Esuel().size()/4 );
  }

  if (m_stage < 2) {

    // continue with next time step stage
    stage();

  } else {

    //// Compute diagnostics, e.g., residuals
    //auto diag_computed = m_diag.compute( *d, m_u.nunk()-myGhosts()->m_fd.Esuel().size()/4,
    //                                     myGhosts()->m_geoElem, m_ndof, m_u );

    // Increase number of iterations and physical time
    d->next();

    // Continue to mesh refinement (if configured)
    /*if (!diag_computed)*/ refine( std::vector< tk::real >( m_u.nprop(), 0.0 ) );

  }
}

void
FV::refine( [[maybe_unused]] const std::vector< tk::real >& l2res )
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
    d->Ref()->dtref( myGhosts()->m_fd.Bface(), {},
      tk::remap(myGhosts()->m_fd.Triinpoel(),d->Gid()) );
    d->refined() = 1;

  } else {      // do not refine

    d->refined() = 0;
    stage();

  }
}

void
FV::resizePostAMR(
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
  [[maybe_unused]] auto old_nelem = myGhosts()->m_inpoel.size()/4;

  // Resize mesh data structures
  d->resizePostAMR( chunk, coord, nodeCommMap );

  // Update state
  myGhosts()->m_inpoel = d->Inpoel();
  myGhosts()->m_coord = d->Coord();
  auto nelem = myGhosts()->m_inpoel.size()/4;
  m_p.resize( nelem );
  m_u.resize( nelem );
  m_un.resize( nelem );
  m_lhs.resize( nelem );
  m_rhs.resize( nelem );

  myGhosts()->m_fd = FaceData( myGhosts()->m_inpoel, bface,
    tk::remap(triinpoel,d->Lid()) );

  myGhosts()->m_geoFace =
    tk::Fields( tk::genGeoFaceTri( myGhosts()->m_fd.Nipfac(),
    myGhosts()->m_fd.Inpofa(), coord ) );
  myGhosts()->m_geoElem = tk::Fields( tk::genGeoElemTet( myGhosts()->m_inpoel,
    coord ) );

  myGhosts()->m_nfac = myGhosts()->m_fd.Inpofa().size()/3;
  m_nunk = nelem;
  m_npoin = coord[0].size();
  myGhosts()->m_bndFace.clear();
  myGhosts()->m_exptGhost.clear();
  myGhosts()->m_sendGhost.clear();
  myGhosts()->m_ghost.clear();
  myGhosts()->m_esup.clear();

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

  // Enable SDAG wait for initially building the solution vector and limiting
  if (m_initial) {
    thisProxy[ thisIndex ].wait4sol();
    thisProxy[ thisIndex ].wait4lim();
    thisProxy[ thisIndex ].wait4nod();
  }

  // Resize communication buffers
  m_ghosts[thisIndex].resizeComm();
}

bool
FV::fieldOutput() const
// *****************************************************************************
// Decide wether to output field data
//! \return True if field data is output in this step
// *****************************************************************************
{
  auto d = Disc();

  // output field data if field iteration count is reached or if the field
  // physics time output frequency is hit or in the last time step
  return d->fielditer() or d->fieldtime() or d->finished();
}

bool
FV::refinedOutput() const
// *****************************************************************************
// Decide if we write field output using a refined mesh
//! \return True if field output will use a refined mesh
// *****************************************************************************
{
  return g_inputdeck.get< tag::cmd, tag::io, tag::refined >() &&
         g_inputdeck.get< tag::discr, tag::scheme >() != ctr::SchemeType::FV;
}

void
FV::writeFields( CkCallback c )
// *****************************************************************************
// Output mesh field data
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  auto d = Disc();

  // Output time history if we hit its output frequency
  if (d->histiter() or d->histtime()) {
    std::vector< std::vector< tk::real > > hist;
    for (const auto& eq : g_fvpde) {
      auto h = eq.histOutput( d->Hist(), myGhosts()->m_inpoel,
        myGhosts()->m_coord, m_u, m_p );
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
  for (const auto& eq : g_fvpde) {
    analyticFieldNames( eq, tk::Centering::ELEM, elemfieldnames );
    analyticFieldNames( eq, tk::Centering::NODE, nodefieldnames );
  }

  Assert( elemfieldnames.size() == m_elemfields.size(), "Size mismatch" );
  Assert( nodefieldnames.size() == m_nodefields.size(), "Size mismatch" );

  // Output chare mesh and fields metadata to file
  const auto& triinpoel = m_outmesh.triinpoel;
  d->write( inpoel, m_outmesh.coord, m_outmesh.bface, {},
            tk::remap( triinpoel, lid ), elemfieldnames, nodefieldnames,
            {}, m_elemfields, m_nodefields, {}, c );
}

void
FV::comnodeout( const std::vector< std::size_t >& gid,
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
FV::stage()
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
    startFieldOutput( CkCallback(CkIndex_FV::step(), thisProxy[thisIndex]) );
}

void
FV::evalLB( int nrestart )
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
FV::evalRestart()
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
FV::step()
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

#include "NoWarning/fv.def.h"
