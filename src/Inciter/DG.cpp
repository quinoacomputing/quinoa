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
        const CProxy_Ghosts& ghostsproxy,
        const std::map< int, std::vector< std::size_t > >& bface,
        const std::map< int, std::vector< std::size_t > >& /* bnode */,
        const std::vector< std::size_t >& triinpoel ) :
  m_disc( disc ),
  m_ghosts( ghostsproxy ),
  m_ndof_NodalExtrm( 3 ), // for the first order derivatives in 3 directions
  m_nsol( 0 ),
  m_ninitsol( 0 ),
  m_nlim( 0 ),
  m_nnod( 0 ),
  m_nreco( 0 ),
  m_nnodalExtrema( 0 ),
  m_u( Disc()->Inpoel().size()/4,
       g_inputdeck.get< tag::discr, tag::rdof >()*
       g_inputdeck.get< tag::component >().nprop() ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_p( m_u.nunk(),
       g_inputdeck.get< tag::discr, tag::rdof >()*
         std::accumulate( begin(g_dgpde), end(g_dgpde), 0u,
           [](std::size_t s, const DGPDE& eq){ return s + eq.nprim(); } ) ),
  m_lhs( m_u.nunk(),
         g_inputdeck.get< tag::discr, tag::ndof >()*
         g_inputdeck.get< tag::component >().nprop() ),
  m_rhs( m_u.nunk(), m_lhs.nprop() ),
  m_uNodalExtrm(),
  m_pNodalExtrm(),
  m_uNodalExtrmc(),
  m_pNodalExtrmc(),
  m_npoin( Disc()->Coord()[0].size() ),
  m_diag(),
  m_stage( 0 ),
  m_ndof(),
  m_numEqDof(),
  m_uc(),
  m_pc(),
  m_ndofc(),
  m_initial( 1 ),
  m_uElemfields(m_u.nunk(), g_inputdeck.get< tag::component >().nprop()),
  m_pElemfields(m_u.nunk(),
    m_p.nprop()/g_inputdeck.get< tag::discr, tag::rdof >()),
  m_uNodefields(m_npoin, g_inputdeck.get< tag::component >().nprop()),
  m_pNodefields(m_npoin,
    m_p.nprop()/g_inputdeck.get< tag::discr, tag::rdof >()),
  m_uNodefieldsc(),
  m_pNodefieldsc(),
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

  // Enable SDAG wait for initially building the solution vector and limiting
  if (m_initial) {
    thisProxy[ thisIndex ].wait4sol();
    thisProxy[ thisIndex ].wait4lim();
    thisProxy[ thisIndex ].wait4nod();
    thisProxy[ thisIndex ].wait4reco();
    thisProxy[ thisIndex ].wait4nodalExtrema();
  }

  m_ghosts[thisIndex].insert(m_disc, bface, triinpoel, m_u.nunk(),
    CkCallback(CkIndex_DG::resizeSolVectors(), thisProxy[thisIndex]));

  // global-sync to call doneInserting on m_ghosts
  auto meshid = Disc()->MeshId();
  contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
    CkCallback(CkReductionTarget(Transporter,doneInsertingGhosts),
    Disc()->Tr()) );
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
DG::resizeSolVectors()
// *****************************************************************************
// Resize solution vectors after extension due to Ghosts and continue with setup
// *****************************************************************************
{
  // Resize solution vectors, lhs and rhs by the number of ghost tets
  m_u.resize( myGhosts()->m_nunk );
  m_un.resize( myGhosts()->m_nunk );
  m_p.resize( myGhosts()->m_nunk );
  m_lhs.resize( myGhosts()->m_nunk );
  m_rhs.resize( myGhosts()->m_nunk );

  // Size communication buffer for solution and number of degrees of freedom
  for (auto& n : m_ndofc) n.resize( myGhosts()->m_bid.size() );
  for (auto& u : m_uc) u.resize( myGhosts()->m_bid.size() );
  for (auto& p : m_pc) p.resize( myGhosts()->m_bid.size() );

  // Initialize number of degrees of freedom in mesh elements
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();
  if( pref )
  {
    const auto ndofmax = g_inputdeck.get< tag::pref, tag::ndofmax >();
    m_ndof.resize( myGhosts()->m_nunk, ndofmax );
  }
  else
  {
    const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
    m_ndof.resize( myGhosts()->m_nunk, ndof );
  }

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
  Assert( myGhosts()->m_geoElem.nunk() == m_lhs.nunk(),
    "Size mismatch in DG::setup()" );

  // Compute left-hand side of discrete PDEs
  lhs();

  // Determine elements inside user-defined IC box
  for (auto& eq : g_dgpde)
    eq.IcBoxElems( myGhosts()->m_geoElem, myGhosts()->m_fd.Esuel().size()/4,
      m_boxelems );

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
    eq.initialize( m_lhs, myGhosts()->m_inpoel, myGhosts()->m_coord, m_boxelems,
      m_u, d->T(), myGhosts()->m_fd.Esuel().size()/4 );
    eq.updatePrimitives( m_u, m_lhs, myGhosts()->m_geoElem, m_p,
      myGhosts()->m_fd.Esuel().size()/4 );
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
DG::next()
// *****************************************************************************
// Advance equations to next time step
// *****************************************************************************
{
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

  auto d = Disc();

  if (pref && m_stage == 0 && d->T() > 0)
    for (const auto& eq : g_dgpde)
      eq.eval_ndof( myGhosts()->m_nunk, myGhosts()->m_coord, myGhosts()->m_inpoel,
                    myGhosts()->m_fd, m_u, m_p,
                    g_inputdeck.get< tag::pref, tag::indicator >(),
                    g_inputdeck.get< tag::discr, tag::ndof >(),
                    g_inputdeck.get< tag::pref, tag::ndofmax >(),
                    g_inputdeck.get< tag::pref, tag::tolref >(),
                    m_ndof );

  // communicate solution ghost data (if any)
  if (myGhosts()->m_sendGhost.empty())
    comsol_complete();
  else
    for(const auto& [cid, ghostdata] : myGhosts()->m_sendGhost) {
      std::vector< std::size_t > tetid( ghostdata.size() );
      std::vector< std::vector< tk::real > > u( ghostdata.size() ),
                                             prim( ghostdata.size() );
      std::vector< std::size_t > ndof;
      std::size_t j = 0;
      for(const auto& i : ghostdata) {
        Assert( i < myGhosts()->m_fd.Esuel().size()/4,
          "Sending solution ghost data" );
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
  const auto& n = tk::cref_find( myGhosts()->m_ghost, fromch );

  for (std::size_t i=0; i<tetid.size(); ++i) {
    auto j = tk::cref_find( n, tetid[i] );
    Assert( j >= myGhosts()->m_fd.Esuel().size()/4,
      "Receiving solution non-ghost data" );
    auto b = tk::cref_find( myGhosts()->m_bid, j );
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
  if (++m_nsol == myGhosts()->m_sendGhost.size()) {
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

  // Evaluate element solution on incoming mesh
  evalSolution( inpoel, coord, addedTets );

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
          lu[f][j++] = m_uNodefields(tk::cref_find(lid,g),f,0);
      }
      for (std::size_t f=0; f<m_pNodefields.nprop(); ++f) {
        std::size_t j = 0;
        for (auto g : nodes)
          lp[f][j++] = m_pNodefields(tk::cref_find(lid,g),f,0);
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

  ownnod_complete( c, addedTets );
}

void
DG::evalSolution(
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const std::unordered_map< std::size_t, std::size_t >& addedTets )
// *****************************************************************************
// Evaluate solution on incoming (a potentially refined) mesh
//! \param[in] inpoel Incoming (potentially refined field-output) mesh
//!   connectivity
//! \param[in] coord Incoming (potentially refined Field-output) mesh node
//!   coordinates
//! \param[in] addedTets Field-output mesh cells and their parents (local ids)
//! \details This function evaluates the solution on the incoming mesh, and
//!   stores it in uElemfields, pElemfields, uNodefields, and pNodefields
//!   appropriately. The incoming mesh can be refined but can also be just the
//!   mesh the numerical solution is computed on.
//! \note If the incoming mesh is refined (for field output) compared to the
//!   mesh the numerical solution is computed on, the solution is evaluated in
//!   cells as wells as in nodes. If the solution is not refined, the solution
//!   is evaluated in nodes.
// *****************************************************************************
{
  using tk::dot;
  using tk::real;

  const auto nelem = inpoel.size()/4;
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
  const auto uncomp = m_u.nprop() / rdof;
  const auto pncomp = m_p.nprop() / rdof;

  // If mesh is not refined for field output, cut off ghosts from element
  // solution. (No need to output ghosts and writer would error.) If mesh is
  // refined for field output, resize element solution fields to refined mesh.
  m_uElemfields.resize( nelem );
  m_pElemfields.resize( nelem );

  auto npoin = coord[0].size();
  m_uNodefields.resize( npoin );
  m_pNodefields.resize( npoin );
  m_uNodefields.fill(0.0);
  m_pNodefields.fill(0.0);

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Assign values to element-fields
  for (std::size_t e=0; e<m_u.nunk(); ++e) {
    if (e < nelem) {
      for (std::size_t i=0; i<uncomp; ++i) {
        m_uElemfields(e,i,0) = m_u(e,rdof*i,0);
      }
      for (std::size_t i=0; i<pncomp; ++i) {
        m_pElemfields(e,i,0) = m_p(e,rdof*i,0);
      }
    }
  }

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
        for (std::size_t i=0; i<uncomp; ++i) m_uNodefields(inpoel[e4+j],i,0) += u[i];
        for (std::size_t i=0; i<pncomp; ++i) m_pNodefields(inpoel[e4+j],i,0) += p[i];
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
      for (std::size_t i=0; i<uncomp; ++i) m_uElemfields(child,i,0) = u[i];
      for (std::size_t i=0; i<pncomp; ++i) m_pElemfields(child,i,0) = p[i];
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
        for (std::size_t i=0; i<uncomp; ++i) m_uNodefields(inpoel[c4+j],i,0) += cnu[i];
        for (std::size_t i=0; i<pncomp; ++i) m_pNodefields(inpoel[c4+j],i,0) += cnp[i];
      }
    }
  }
}

void
DG::lhs()
// *****************************************************************************
// Compute left-hand side of discrete transport equations
// *****************************************************************************
{
  for (const auto& eq : g_dgpde) eq.lhs( myGhosts()->m_geoElem, m_lhs );

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
  for (const auto& b : myGhosts()->m_bid) {
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
      eq.reconstruct( d->T(), myGhosts()->m_geoFace, myGhosts()->m_geoElem,
        myGhosts()->m_fd, myGhosts()->m_esup, myGhosts()->m_inpoel,
        myGhosts()->m_coord, m_u, m_p );
  }

  // Send reconstructed solution to neighboring chares
  if (myGhosts()->m_sendGhost.empty())
    comreco_complete();
  else
    for(const auto& [cid, ghostdata] : myGhosts()->m_sendGhost) {
      std::vector< std::size_t > tetid( ghostdata.size() );
      std::vector< std::vector< tk::real > > u( ghostdata.size() ),
                                             prim( ghostdata.size() );
      std::vector< std::size_t > ndof;
      std::size_t j = 0;
      for(const auto& i : ghostdata) {
        Assert( i < myGhosts()->m_fd.Esuel().size()/4, "Sending reconstructed ghost "
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
    if (pref && m_stage == 0) {
      Assert( b < m_ndofc[1].size(), "Indexing out of bounds" );
      m_ndofc[1][b] = ndof[i];
    }
  }

  // if we have received all solution ghost contributions from neighboring
  // chares (chares we communicate along chare-boundary faces with), and
  // contributed our solution to these neighbors, proceed to limiting
  if (++m_nreco == myGhosts()->m_sendGhost.size()) {
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
  for (const auto& [boundary, localtet] : myGhosts()->m_bid) {
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
      evalNodalExtrm(ncomp, nprim, m_ndof_NodalExtrm, d->bndel(),
        myGhosts()->m_inpoel, myGhosts()->m_coord, gid, bid, m_u, m_p,
        m_uNodalExtrm, m_pNodalExtrm);
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
    for (const auto& eq : g_dgpde) {
      eq.limit( d->T(), myGhosts()->m_geoFace, myGhosts()->m_geoElem,
                myGhosts()->m_fd, myGhosts()->m_esup, myGhosts()->m_inpoel,
                myGhosts()->m_coord, m_ndof, d->Gid(), d->Bid(), m_uNodalExtrm,
                m_pNodalExtrm, m_u, m_p, m_shockmarker );

      eq.Correct_Conserv(m_p, myGhosts()->m_geoElem, m_u,
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
      std::vector< std::size_t > ndof;
      std::size_t j = 0;
      for(const auto& i : ghostdata) {
        Assert( i < myGhosts()->m_fd.Esuel().size()/4,
          "Sending limiter ghost data" );
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
  const auto& esuf = myGhosts()->m_fd.Esuf();

  // Copy number of degrees of freedom for each cell
  auto ndof = m_ndof;

  // p-refine all neighboring elements of elements that have been p-refined as a
  // result of error indicators
  for( auto f=myGhosts()->m_fd.Nbfac(); f<esuf.size()/2; ++f )
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
    if (pref && m_stage == 0) {
      Assert( b < m_ndofc[2].size(), "Indexing out of bounds" );
      m_ndofc[2][b] = ndof[i];
    }
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
DG::dt()
// *****************************************************************************
// Compute time step size
// *****************************************************************************
{
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();
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
          eq.dt( myGhosts()->m_coord, myGhosts()->m_inpoel, myGhosts()->m_fd,
            myGhosts()->m_geoFace, myGhosts()->m_geoElem, m_ndof, m_u, m_p,
            myGhosts()->m_fd.Esuel().size()/4 );
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
    for(std::size_t e = 0; e < myGhosts()->m_nunk; e++)
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
    eq.rhs( d->T(), myGhosts()->m_geoFace, myGhosts()->m_geoElem,
      myGhosts()->m_fd, myGhosts()->m_inpoel, m_boxelems, myGhosts()->m_coord,
      m_u, m_p, m_ndof, m_rhs );

  // Explicit time-stepping using RK3 to discretize time-derivative
  for(std::size_t e=0; e<myGhosts()->m_nunk; ++e)
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
    eq.updateInterfaceCells( m_u, myGhosts()->m_fd.Esuel().size()/4, m_ndof );
    eq.updatePrimitives( m_u, m_lhs, myGhosts()->m_geoElem, m_p,
      myGhosts()->m_fd.Esuel().size()/4 );
    if (!g_inputdeck.get< tag::discr, tag::accuracy_test >()) {
      eq.cleanTraceMaterial( myGhosts()->m_geoElem, m_u, m_p,
        myGhosts()->m_fd.Esuel().size()/4 );
    }
  }

  if (m_stage < 2) {

    // continue with next time step stage
    stage();

  } else {

    // Increase number of iterations and physical time
    d->next();

    // Compute diagnostics, e.g., residuals
    auto diag_computed = m_diag.compute( *d,
      m_u.nunk()-myGhosts()->m_fd.Esuel().size()/4, myGhosts()->m_geoElem,
      m_ndof, m_u );

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
    d->Ref()->dtref( myGhosts()->m_fd.Bface(), {},
      tk::remap(myGhosts()->m_fd.Triinpoel(),d->Gid()) );
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
  const std::set< std::size_t >& removedNodes,
  const std::unordered_map< std::size_t, std::size_t >& amrNodeMap,
  const tk::NodeCommMap& nodeCommMap,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& /* bnode */,
  const std::vector< std::size_t >& triinpoel )
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
  d->resizePostAMR( chunk, coord, amrNodeMap, nodeCommMap, removedNodes );

  // Update state
  myGhosts()->m_inpoel = d->Inpoel();
  myGhosts()->m_coord = d->Coord();
  auto nelem = myGhosts()->m_inpoel.size()/4;
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
    for (std::size_t i=0; i<unprop; ++i) m_u(child,i,0) = m_un(parent,i,0);
    for (std::size_t i=0; i<pnprop; ++i) m_p(child,i,0) = pn(parent,i,0);
  }
  m_un = m_u;

  // Resize communication buffers
  m_ghosts[thisIndex].resizeComm();
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
DG::writeFields(
  CkCallback c,
  const std::unordered_map< std::size_t, std::size_t >& addedTets )
// *****************************************************************************
// Output mesh field data
//! \param[in] c Function to continue with after the write
//! \param[in] addedTets Newly added mesh cells and their parents (local ids)
// *****************************************************************************
{
  auto d = Disc();

  // Output time history
  if (d->histiter() or d->histtime() or d->histrange()) {
    std::vector< std::vector< tk::real > > hist;
    for (const auto& eq : g_dgpde) {
      auto h = eq.histOutput( d->Hist(), myGhosts()->m_inpoel,
        myGhosts()->m_coord, m_u, m_p );
      hist.insert( end(hist), begin(h), end(h) );
    }
    d->history( std::move(hist) );
  }

  const auto& inpoel = std::get< 0 >( m_outmesh.chunk );
  auto esup = tk::genEsup( inpoel, 4 );
  auto nelem = inpoel.size() / 4;

  // Combine own and communicated contributions and finish averaging of node
  // field output in chare boundary nodes
  const auto& lid = std::get< 2 >( m_outmesh.chunk );
  for (const auto& [g,f] : m_uNodefieldsc) {
    Assert( m_uNodefields.nprop() == f.first.size(), "Size mismatch" );
    auto p = tk::cref_find( lid, g );
    for (std::size_t i=0; i<f.first.size(); ++i) {
      m_uNodefields(p,i,0) += f.first[i];
      m_uNodefields(p,i,0) /= static_cast< tk::real >(
                              esup.second[p+1] - esup.second[p] + f.second );
    }
  }
  tk::destroy( m_uNodefieldsc );
  for (const auto& [g,f] : m_pNodefieldsc) {
    Assert( m_pNodefields.nprop() == f.first.size(), "Size mismatch" );
    auto p = tk::cref_find( lid, g );
    for (std::size_t i=0; i<f.first.size(); ++i) {
      m_pNodefields(p,i,0) += f.first[i];
      m_pNodefields(p,i,0) /= static_cast< tk::real >(
                              esup.second[p+1] - esup.second[p] + f.second );
    }
  }
  tk::destroy( m_pNodefieldsc );

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
      for (std::size_t i=0; i<m_uNodefields.nprop(); ++i)
        m_uNodefields(p,i,0) /= n;
      for (std::size_t i=0; i<m_pNodefields.nprop(); ++i)
        m_pNodefields(p,i,0) /= n;
    }
  }

  // Collect field output from numerical solution requested by user
  auto elemfields = numericFieldOutput( m_uElemfields, tk::Centering::ELEM,
    m_pElemfields );
  auto nodefields = numericFieldOutput( m_uNodefields, tk::Centering::NODE,
    m_pNodefields );

  // Collect field output from analytical solutions (if exist)
  const auto& coord = m_outmesh.coord;
  auto geoElem = tk::genGeoElemTet( inpoel, coord );
  auto t = Disc()->T();
  for (const auto& eq : g_dgpde) {
    analyticFieldOutput( eq, tk::Centering::ELEM, geoElem.extract(1,0),
      geoElem.extract(2,0), geoElem.extract(3,0), t, elemfields );
    analyticFieldOutput( eq, tk::Centering::NODE, coord[0], coord[1], coord[2],
      t, nodefields );
  }

  // Add adaptive indicator array to element-centered field output
  if (g_inputdeck.get< tag::pref, tag::pref >()) {
    std::vector< tk::real > ndof( begin(m_ndof), end(m_ndof) );
    ndof.resize( nelem );
    for (const auto& [child,parent] : addedTets)
      ndof[child] = static_cast< tk::real >( m_ndof[parent] );
    elemfields.push_back( ndof );
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
  elemfields.push_back( shockmarker );

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

  Assert( elemfieldnames.size() == elemfields.size(), "Size mismatch" );
  Assert( nodefieldnames.size() == nodefields.size(), "Size mismatch" );

  // Output chare mesh and fields metadata to file
  const auto& triinpoel = m_outmesh.triinpoel;
  d->write( inpoel, m_outmesh.coord, m_outmesh.bface, {},
            tk::remap( triinpoel, lid ), elemfieldnames, nodefieldnames,
            {}, elemfields, nodefields, {}, c );
}

void
DG::comnodeout( const std::vector< std::size_t >& gid,
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

  if (not benchmark and not (d->It() % rsfreq)) {

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
