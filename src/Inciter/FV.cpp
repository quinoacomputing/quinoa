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
  m_p( m_u.nunk(), g_inputdeck.get< tag::discr, tag::rdof >()*
    g_fvpde[Disc()->MeshId()].nprim() ),
  m_lhs( m_u.nunk(),
         g_inputdeck.get< tag::component >().nprop() ),
  m_rhs( m_u.nunk(), m_lhs.nprop() ),
  m_npoin( Disc()->Coord()[0].size() ),
  m_diag(),
  m_stage( 0 ),
  m_uc(),
  m_pc(),
  m_initial( 1 ),
  m_uElemfields(m_u.nunk(), g_inputdeck.get< tag::component >().nprop()),
  m_pElemfields(m_u.nunk(),
    m_p.nprop()/g_inputdeck.get< tag::discr, tag::rdof >()),
  m_uNodefields(m_npoin, g_inputdeck.get< tag::component >().nprop()),
  m_pNodefields(m_npoin,
    m_p.nprop()/g_inputdeck.get< tag::discr, tag::rdof >()),
  m_uNodefieldsc(),
  m_pNodefieldsc(),
  m_boxelems(),
  m_srcFlag(m_u.nunk(), 1),
  m_nrk(0)
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
// *****************************************************************************
{
  //! Runge-Kutta coefficients
  m_nrk = 2;
  m_rkcoef[0].resize(m_nrk);
  m_rkcoef[1].resize(m_nrk);
  if (m_nrk == 2) {
    m_rkcoef = {{ {{ 0.0, 1.0/2.0 }}, {{ 1.0, 1.0/2.0 }} }};
  }
  else {
    m_rkcoef = {{ {{ 0.0, 3.0/4.0, 1.0/3.0 }}, {{ 1.0, 1.0/4.0, 2.0/3.0 }} }};
  }

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

  m_ghosts[thisIndex].insert(m_disc, bface, triinpoel, m_u.nunk(),
    CkCallback(CkIndex_FV::resizeSolVectors(), thisProxy[thisIndex]));

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
FV::resizeSolVectors()
// *****************************************************************************
// Resize solution vectors after extension due to Ghosts
// *****************************************************************************
{
  // Resize solution vectors, lhs and rhs by the number of ghost tets
  m_u.resize( myGhosts()->m_nunk );
  m_un.resize( myGhosts()->m_nunk );
  m_srcFlag.resize( myGhosts()->m_nunk );
  m_p.resize( myGhosts()->m_nunk );
  m_lhs.resize( myGhosts()->m_nunk );
  m_rhs.resize( myGhosts()->m_nunk );

  // Size communication buffer for solution
  for (auto& u : m_uc) u.resize( myGhosts()->m_bid.size() );
  for (auto& p : m_pc) p.resize( myGhosts()->m_bid.size() );

  // Ensure that we also have all the geometry and connectivity data
  // (including those of ghosts)
  Assert( myGhosts()->m_geoElem.nunk() == m_u.nunk(),
    "GeoElem unknowns size mismatch" );

  // Signal the runtime system that all workers have received their adjacency
  std::vector< std::size_t > meshdata{ myGhosts()->m_initial, Disc()->MeshId() };
  contribute( meshdata, CkReduction::sum_ulong,
    CkCallback(CkReductionTarget(Transporter,comfinal), Disc()->Tr()) );
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

  auto d = Disc();

  // Basic error checking on sizes of element geometry data and connectivity
  Assert( myGhosts()->m_geoElem.nunk() == m_lhs.nunk(),
    "Size mismatch in FV::setup()" );

  // Compute left-hand side of discrete PDEs
  lhs();

  // Determine elements inside user-defined IC box
  g_fvpde[d->MeshId()].IcBoxElems( myGhosts()->m_geoElem,
    myGhosts()->m_fd.Esuel().size()/4, m_boxelems );

  // Compute volume of user-defined box IC
  d->boxvol( {}, {}, 0 );      // punt for now

  // Query time history field output labels from all PDEs integrated
  const auto& hist_points = g_inputdeck.get< tag::history, tag::point >();
  if (!hist_points.empty()) {
    std::vector< std::string > histnames;
    auto n = g_fvpde[d->MeshId()].histNames();
    histnames.insert( end(histnames), begin(n), end(n) );
    d->histheader( std::move(histnames) );
  }
}

void
FV::box( tk::real v, const std::vector< tk::real >& )
// *****************************************************************************
// Receive total box IC volume and set conditions in box
//! \param[in] v Total volume within user-specified box
// *****************************************************************************
{
  auto d = Disc();

  // Store user-defined box IC volume
  d->Boxvol() = v;

  // Set initial conditions for all PDEs
  g_fvpde[d->MeshId()].initialize( m_lhs, myGhosts()->m_inpoel,
    myGhosts()->m_coord, m_boxelems, d->ElemBlockId(), m_u, d->T(),
    myGhosts()->m_fd.Esuel().size()/4 );
  g_fvpde[d->MeshId()].updatePrimitives( m_u, m_p,
    myGhosts()->m_fd.Esuel().size()/4 );

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
      extractFieldOutput( {}, d->Chunk(), d->Coord(), {}, {}, d->NodeCommMap(),
        {}, {}, {}, c );

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
  const std::map< int, std::vector< std::size_t > >& /* bface */,
  const std::map< int, std::vector< std::size_t > >& /* bnode */,
  const std::vector< std::size_t >& /* triinpoel */,
  CkCallback c )
// *****************************************************************************
// Extract field output going to file
//! \param[in] chunk Field-output mesh chunk (connectivity and global<->local
//!    id maps)
//! \param[in] coord Field-output mesh node coordinates
//! \param[in] addedTets Field-output mesh cells and their parents (local ids)
//! \param[in] nodeCommMap Field-output mesh node communication map
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  const auto& inpoel = std::get< 0 >( chunk );

  // Evaluate element solution on incoming mesh
  evalSolution( *Disc(), inpoel, coord, addedTets, std::vector< std::size_t>{},
    m_u, m_p, m_uElemfields, m_pElemfields, m_uNodefields, m_pNodefields );

  // Send node fields contributions to neighbor chares
  if (nodeCommMap.empty())
    comnodeout_complete();
  else {
    const auto& lid = std::get< 2 >( chunk );
    auto esup = tk::genEsup( inpoel, 4 );
    for(const auto& [ch,nodes] : nodeCommMap) {
      // Pack node field data in chare boundary nodes
      std::vector< std::vector< tk::real > >
        lu( m_uNodefields.nprop(), std::vector< tk::real >( nodes.size() ) );
      std::vector< std::vector< tk::real > >
        lp( m_pNodefields.nprop(), std::vector< tk::real >( nodes.size() ) );
      for (std::size_t f=0; f<m_uNodefields.nprop(); ++f) {
        std::size_t j = 0;
        for (auto g : nodes)
          lu[f][j++] = m_uNodefields(tk::cref_find(lid,g),f);
      }
      for (std::size_t f=0; f<m_pNodefields.nprop(); ++f) {
        std::size_t j = 0;
        for (auto g : nodes)
          lp[f][j++] = m_pNodefields(tk::cref_find(lid,g),f);
      }
      // Pack (partial) number of elements surrounding chare boundary nodes
      std::vector< std::size_t > nesup( nodes.size() );
      std::size_t j = 0;
      for (auto g : nodes) {
        auto i = tk::cref_find( lid, g );
        nesup[j++] = esup.second[i+1] - esup.second[i];
      }
      thisProxy[ch].comnodeout(
        std::vector<std::size_t>(begin(nodes),end(nodes)), nesup, lu, lp );
    }
  }

  ownnod_complete( c );
}

void
FV::lhs()
// *****************************************************************************
// Compute left-hand side of discrete transport equations
// *****************************************************************************
{
  g_fvpde[Disc()->MeshId()].lhs( myGhosts()->m_geoElem, m_lhs );

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
      m_u(b.first,c) = m_uc[0][b.second][c];
    }
    for (std::size_t c=0; c<m_p.nprop(); ++c) {
      m_p(b.first,c) = m_pc[0][b.second][c];
    }
  }

  if (rdof > 1) {
    // Reconstruct second-order solution and primitive quantities
    g_fvpde[Disc()->MeshId()].reconstruct( myGhosts()->m_geoElem, myGhosts()->m_fd,
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

  if (rdof > 1) {
    g_fvpde[Disc()->MeshId()].limit( myGhosts()->m_geoFace, myGhosts()->m_fd,
      myGhosts()->m_esup,
      myGhosts()->m_inpoel, myGhosts()->m_coord, m_u, m_p );

    if (g_inputdeck.get< tag::discr, tag::limsol_projection >())
      g_fvpde[Disc()->MeshId()].CPL(m_p, myGhosts()->m_geoElem,
        myGhosts()->m_inpoel, myGhosts()->m_coord, m_u,
        myGhosts()->m_fd.Esuel().size()/4);
  }

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
    Assert( b < m_uc[1].size(), "Indexing out of bounds" );
    Assert( b < m_pc[1].size(), "Indexing out of bounds" );
    m_uc[1][b] = u[i];
    m_pc[1][b] = prim[i];
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
    Assert( m_uc[1][b.second].size() == m_u.nprop(), "ncomp size mismatch" );
    Assert( m_pc[1][b.second].size() == m_p.nprop(), "ncomp size mismatch" );
    for (std::size_t c=0; c<m_u.nprop(); ++c) {
      m_u(b.first,c) = m_uc[1][b.second][c];
    }
    for (std::size_t c=0; c<m_p.nprop(); ++c) {
      m_p(b.first,c) = m_pc[1][b.second][c];
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
      auto eqdt =
        g_fvpde[d->MeshId()].dt( myGhosts()->m_fd, myGhosts()->m_geoFace,
          myGhosts()->m_geoElem, m_u, m_p, myGhosts()->m_fd.Esuel().size()/4,
          m_srcFlag );
      if (eqdt < mindt) mindt = eqdt;

      tk::real coeff(1.0);
      if (d->It() < 100) coeff = 0.01 * static_cast< tk::real >(d->It());

      mindt *= coeff * g_inputdeck.get< tag::discr, tag::cfl >();
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

  // physical time at time-stage for computing exact source terms
  tk::real physT(d->T());
  // 2-stage RK
  if (m_nrk == 2) {
    if (m_stage == 1) {
      physT += d->Dt();
    }
  }
  // 3-stage RK
  else {
    if (m_stage == 1) {
      physT += d->Dt();
    }
    else if (m_stage == 2) {
      physT += 0.5*d->Dt();
    }
  }

  // initialize energy source as not added (modified in eq.rhs appropriately)
  for (auto& fl : m_srcFlag) fl = 0;
  g_fvpde[d->MeshId()].rhs( physT, myGhosts()->m_geoFace, myGhosts()->m_geoElem,
    myGhosts()->m_fd, myGhosts()->m_inpoel, myGhosts()->m_coord,
    d->ElemBlockId(), m_u, m_p, m_rhs, m_srcFlag );

  // Explicit time-stepping using RK3 to discretize time-derivative
  for (std::size_t e=0; e<myGhosts()->m_nunk; ++e)
    for (std::size_t c=0; c<neq; ++c)
    {
      auto rmark = c*rdof;
      m_u(e, rmark) =  m_rkcoef[0][m_stage] * m_un(e, rmark)
        + m_rkcoef[1][m_stage] * ( m_u(e, rmark)
          + d->Dt() * m_rhs(e, c)/m_lhs(e, c) );
      // zero out reconstructed dofs of equations using reduced dofs
      if (rdof > 1) {
        for (std::size_t k=1; k<rdof; ++k)
        {
          rmark = c*rdof+k;
          m_u(e, rmark) = 0.0;
        }
      }
    }

  // Update primitives based on the evolved solution
  g_fvpde[d->MeshId()].updatePrimitives( m_u, m_p,
    myGhosts()->m_fd.Esuel().size()/4 );
  if (!g_inputdeck.get< tag::discr, tag::accuracy_test >()) {
    g_fvpde[d->MeshId()].cleanTraceMaterial( myGhosts()->m_geoElem, m_u, m_p,
      myGhosts()->m_fd.Esuel().size()/4 );
  }

  if (m_stage < m_nrk-1) {

    // continue with next time step stage
    stage();

  } else {

    // Increase number of iterations and physical time
    d->next();

    // Compute diagnostics, e.g., residuals
    auto diag_computed = m_diag.compute( *d,
      m_u.nunk()-myGhosts()->m_fd.Esuel().size()/4, myGhosts()->m_geoElem,
      std::vector< std::size_t>{}, m_u );

    // Continue to mesh refinement (if configured)
    if (!diag_computed) refine( std::vector< tk::real >( m_u.nprop(), 0.0 ) );

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
  const std::set< std::size_t >& removedNodes,
  const std::unordered_map< std::size_t, std::size_t >& amrNodeMap,
  const tk::NodeCommMap& nodeCommMap,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& /* bnode */,
  const std::vector< std::size_t >& triinpoel,
  const std::unordered_map< std::size_t, std::set< std::size_t > >& elemblockid )
// *****************************************************************************
//  Receive new mesh from Refiner
//! \param[in] chunk New mesh chunk (connectivity and global<->local id maps)
//! \param[in] coord New mesh node coordinates
//! \param[in] addedTets Newly added mesh cells and their parents (local ids)
//! \param[in] removedNodes Newly removed mesh node local ids
//! \param[in] amrNodeMap Node id map after amr (local ids)
//! \param[in] nodeCommMap New node communication map
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
//! \param[in] elemblockid Local tet ids associated with mesh block ids
// *****************************************************************************
{
  auto d = Disc();

  // Set flag that indicates that we are during time stepping
  m_initial = 0;
  myGhosts()->m_initial = 0;

  // Zero field output iteration count between two mesh refinement steps
  d->Itf() = 0;

  // Increase number of iterations with mesh refinement
  ++d->Itr();

  // Save old number of elements
  [[maybe_unused]] auto old_nelem = myGhosts()->m_inpoel.size()/4;

  // Resize mesh data structures
  d->resizePostAMR( chunk, coord, amrNodeMap, nodeCommMap, removedNodes,
    elemblockid );

  // Update state
  myGhosts()->m_inpoel = d->Inpoel();
  myGhosts()->m_coord = d->Coord();
  auto nelem = myGhosts()->m_inpoel.size()/4;
  m_p.resize( nelem );
  m_u.resize( nelem );
  m_srcFlag.resize( nelem );
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
  myGhosts()->m_nunk = nelem;
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
    for (std::size_t i=0; i<unprop; ++i) m_u(child,i) = m_un(parent,i);
    for (std::size_t i=0; i<pnprop; ++i) m_p(child,i) = pn(parent,i);
  }
  m_un = m_u;

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

  // Output field data
  return d->fielditer() or d->fieldtime() or d->fieldrange() or d->finished();
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

  const auto& inpoel = std::get< 0 >( d->Chunk() );
  auto esup = tk::genEsup( inpoel, 4 );

  // Combine own and communicated contributions and finish averaging of node
  // field output in chare boundary nodes
  const auto& lid = std::get< 2 >( d->Chunk() );
  for (const auto& [g,f] : m_uNodefieldsc) {
    Assert( m_uNodefields.nprop() == f.first.size(), "Size mismatch" );
    auto p = tk::cref_find( lid, g );
    for (std::size_t i=0; i<f.first.size(); ++i) {
      m_uNodefields(p,i) += f.first[i];
      m_uNodefields(p,i) /= static_cast< tk::real >(
                              esup.second[p+1] - esup.second[p] + f.second );
    }
  }
  tk::destroy( m_uNodefieldsc );
  for (const auto& [g,f] : m_pNodefieldsc) {
    Assert( m_pNodefields.nprop() == f.first.size(), "Size mismatch" );
    auto p = tk::cref_find( lid, g );
    for (std::size_t i=0; i<f.first.size(); ++i) {
      m_pNodefields(p,i) += f.first[i];
      m_pNodefields(p,i) /= static_cast< tk::real >(
                              esup.second[p+1] - esup.second[p] + f.second );
    }
  }
  tk::destroy( m_pNodefieldsc );

  // Lambda to decide if a node (global id) is on a chare boundary of the field
  // output mesh. p - global node id, return true if node is on the chare
  // boundary.
  auto chbnd = [ this ]( std::size_t p ) {
    return
      std::any_of( Disc()->NodeCommMap().cbegin(), Disc()->NodeCommMap().cend(),
        [&](const auto& s) { return s.second.find(p) != s.second.cend(); } );
  };

  // Finish computing node field output averages in internal nodes
  auto npoin = d->Coord()[0].size();
  auto& gid = std::get< 1 >( d->Chunk() );
  for (std::size_t p=0; p<npoin; ++p) {
    if (!chbnd(gid[p])) {
      auto n = static_cast< tk::real >( esup.second[p+1] - esup.second[p] );
      for (std::size_t i=0; i<m_uNodefields.nprop(); ++i)
        m_uNodefields(p,i) /= n;
      for (std::size_t i=0; i<m_pNodefields.nprop(); ++i)
        m_pNodefields(p,i) /= n;
    }
  }

  // Collect field output from numerical solution requested by user
  auto elemfields = numericFieldOutput( m_uElemfields, tk::Centering::ELEM,
    m_pElemfields );
  auto nodefields = numericFieldOutput( m_uNodefields, tk::Centering::NODE,
    m_pNodefields );

  // Collect field output from analytical solutions (if exist)
  const auto& coord = d->Coord();
  auto geoElem = tk::genGeoElemTet( inpoel, coord );
  auto t = Disc()->T();
  analyticFieldOutput( g_fvpde[d->MeshId()], tk::Centering::ELEM,
    geoElem.extract_comp(1), geoElem.extract_comp(2), geoElem.extract_comp(3),
    t, elemfields );
  analyticFieldOutput( g_fvpde[d->MeshId()], tk::Centering::NODE, coord[0],
    coord[1], coord[2], t, nodefields );

  // Query fields names requested by user
  auto elemfieldnames = numericFieldNames( tk::Centering::ELEM );
  auto nodefieldnames = numericFieldNames( tk::Centering::NODE );

  // Collect field output names for analytical solutions
  analyticFieldNames( g_fvpde[d->MeshId()], tk::Centering::ELEM, elemfieldnames );
  analyticFieldNames( g_fvpde[d->MeshId()], tk::Centering::NODE, nodefieldnames );

  Assert( elemfieldnames.size() == elemfields.size(), "Size mismatch" );
  Assert( nodefieldnames.size() == nodefields.size(), "Size mismatch" );

  // Collect surface output names
  auto surfnames = g_fvpde[d->MeshId()].surfNames();

  // Collect surface field solution
  const auto& fd = myGhosts()->m_fd;
  auto elemsurfs = g_fvpde[d->MeshId()].surfOutput(fd, m_u, m_p);

  // Output chare mesh and fields metadata to file
  const auto& triinpoel = tk::remap( fd.Triinpoel(), d->Gid() );
  d->write( inpoel, d->Coord(), fd.Bface(), {},
            tk::remap( triinpoel, lid ), elemfieldnames, nodefieldnames,
            surfnames, {}, elemfields, nodefields, elemsurfs, {}, c );
}

void
FV::comnodeout( const std::vector< std::size_t >& gid,
                const std::vector< std::size_t >& nesup,
                const std::vector< std::vector< tk::real > >& Lu,
                const std::vector< std::vector< tk::real > >& Lp )
// *****************************************************************************
//  Receive chare-boundary nodal solution (for field output) contributions from
//  neighboring chares
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] nesup Number of elements surrounding points
//! \param[in] Lu Partial contributions of solution nodal fields to
//!   chare-boundary nodes
//! \param[in] Lp Partial contributions of primitive quantity nodal fields to
//!   chare-boundary nodes
// *****************************************************************************
{
  Assert( gid.size() == nesup.size(), "Size mismatch" );
  Assert(Lu.size() == m_uNodefields.nprop(), "Fields size mismatch");
  Assert(Lp.size() == m_pNodefields.nprop(), "Fields size mismatch");
  for (std::size_t f=0; f<Lu.size(); ++f)
    Assert( gid.size() == Lu[f].size(), "Size mismatch" );
  for (std::size_t f=0; f<Lp.size(); ++f)
    Assert( gid.size() == Lp[f].size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i) {
    auto& nfu = m_uNodefieldsc[ gid[i] ];
    nfu.first.resize( Lu.size() );
    for (std::size_t f=0; f<Lu.size(); ++f) nfu.first[f] += Lu[f][i];
    nfu.second += nesup[i];
    auto& nfp = m_pNodefieldsc[ gid[i] ];
    nfp.first.resize( Lp.size() );
    for (std::size_t f=0; f<Lp.size(); ++f) nfp.first[f] += Lp[f][i];
    nfp.second += nesup[i];
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
  if (m_stage < m_nrk)
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
  auto d = Disc();

  // Output time history
  if (d->histiter() or d->histtime() or d->histrange()) {
    std::vector< std::vector< tk::real > > hist;
    auto h = g_fvpde[d->MeshId()].histOutput( d->Hist(), myGhosts()->m_inpoel,
      myGhosts()->m_coord, m_u, m_p );
    hist.insert( end(hist), begin(h), end(h) );
    d->history( std::move(hist) );
  }

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
