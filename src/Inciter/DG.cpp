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
#include "CGPDE.hpp"
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
#include "PDE/MultiMat/MultiMatIndexing.hpp"
#include "Integrate/Quadrature.hpp"

#include <fstream>

// ignore old-style-casts required for lapack/blas calls
#if defined(__clang__)
  #pragma clang diagnostic ignored "-Wold-style-cast"
#endif

// Lapacke forward declarations
extern "C" {

using lapack_int = long;

#define LAPACK_ROW_MAJOR 101

extern lapack_int LAPACKE_dgesv( int, lapack_int, lapack_int, double*,
  lapack_int, lapack_int*, double*, lapack_int );

}

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< DGPDE > g_dgpde;

//! Runge-Kutta coefficients
static const std::array< std::array< tk::real, 3 >, 2 >
  rkcoef{{ {{ 0.0, 3.0/4.0, 1.0/3.0 }}, {{ 1.0, 1.0/4.0, 2.0/3.0 }} }};

//! Implicit-Explicit Runge-Kutta Coefficients
[[maybe_unused]] static const tk::real rk_gamma = (2.0-std::sqrt(2.0))/2.0;
[[maybe_unused]] static const tk::real rk_delta = -2.0*std::sqrt(2.0)/3.0;
static const tk::real c2 =
  (27.0 + std::pow(2187.0-1458.0*std::sqrt(2.0),1.0/3.0)
   + 9.0*std::pow(3.0+2.0*std::sqrt(2.0),1.0/3.0))/54.0;
static const tk::real c3 = c2/(6.0*std::pow(c2,2.0)-3.0*c2+1.0);
static const tk::real b2 = (3.0*c2-1.0)/(6.0*std::pow(c2,2.0));
static const tk::real b3 =
  (6.0*std::pow(c2,2.0)-3.0*c2+1.0)/(6.0*std::pow(c2,2.0));
static const tk::real a22_impl = c2;
static const tk::real a21_expl = c2;
static const tk::real a32_expl = c3;
static const tk::real a33_impl =
  (1.0/6.0-b2*std::pow(c2,2.0)-b3*c2*c3)/(b3*(c3-c2));
static const tk::real a32_impl = a33_impl-c3;
static const std::array< std::array< tk::real, 3 >, 2 >
  expl_rkcoef{{ {{ 0.0, 0.0, b2 }},
                {{ a21_expl, a32_expl, b3 }} }};
static const std::array< std::array< tk::real, 3 >, 2>
  impl_rkcoef{{ {{ 0.0, a32_impl, b2 }},
                {{ a22_impl, a33_impl, b3}} }};
static const std::array< tk::real, 10 > mass_dubiner( tk::massMatrixDubiner() );

} // inciter::

extern tk::CProxy_ChareStateCollector stateProxy;

using inciter::DG;

DG::DG( const CProxy_Discretization& disc,
        const CProxy_Ghosts& ghostsproxy,
        const std::map< int, std::vector< std::size_t > >& bface,
        const std::map< int, std::vector< std::size_t > >& bnode,
        const std::vector< std::size_t >& triinpoel ) :
  m_disc( disc ),
  m_ghosts( ghostsproxy ),
  m_nsol( 0 ),
  m_ninitsol( 0 ),
  m_nlim( 0 ),
  m_nnod( 0 ),
  m_nrefine( 0 ),
  m_nsmooth( 0 ),
  m_nreco( 0 ),
  m_nale( 0 ),
  m_nbnorm( 0 ),
  m_nstiffeq( g_dgpde[Disc()->MeshId()].nstiffeq() ),
  m_nnonstiffeq( g_dgpde[Disc()->MeshId()].nnonstiffeq() ),
  m_u( Disc()->Inpoel().size()/4,
       g_inputdeck.get< tag::rdof >()*
       g_inputdeck.get< tag::ncomp >() ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_p( m_u.nunk(), g_inputdeck.get< tag::rdof >()*
    g_dgpde[Disc()->MeshId()].nprim() ),
  m_geoElemk( m_u.nunk(), 5 ),
  m_geoElemn( m_u.nunk(), 5 ),
  m_rhs( m_u.nunk(),
         g_inputdeck.get< tag::ndof >()*
         g_inputdeck.get< tag::ncomp >() ),
  m_rhsprev( m_u.nunk(), m_rhs.nprop() ),
  m_stiffrhs( m_u.nunk(), g_inputdeck.get< tag::ndof >()*
              g_dgpde[Disc()->MeshId()].nstiffeq() ),
  m_stiffrhsprev( m_u.nunk(), g_inputdeck.get< tag::ndof >()*
                  g_dgpde[Disc()->MeshId()].nstiffeq() ),
  m_stiffEqIdx( g_dgpde[Disc()->MeshId()].nstiffeq() ),
  m_nonStiffEqIdx( g_dgpde[Disc()->MeshId()].nnonstiffeq() ),
  m_mtInv(
    tk::invMassMatTaylorRefEl(g_inputdeck.get< tag::rdof >()) ),
  m_npoin( Disc()->Coord()[0].size() ),
  m_diag(),
  m_nstage( 3 ),
  m_stage( 0 ),
  m_ndof(),
  m_interface(),
  m_numEqDof(),
  m_uc(),
  m_pc(),
  m_ndofc(),
  m_interfacec(),
  m_initial( 1 ),
  m_uElemfields( m_u.nunk(),
                 g_inputdeck.get< tag::ncomp >() ),
  m_pElemfields( m_u.nunk(),
                 m_p.nprop() / g_inputdeck.get< tag::rdof >() ),
  m_uNodefields( m_npoin,
                 g_inputdeck.get< tag::ncomp >() ),
  m_pNodefields( m_npoin,
                 m_p.nprop() / g_inputdeck.get< tag::rdof >() ),
  m_uNodefieldsc(),
  m_pNodefieldsc(),
  m_outmesh(),
  m_boxelems(),
  m_shockmarker(m_u.nunk(), 1),
  m_nodevel( {{ std::vector<tk::real>(Disc()->Lid().size(), 0.0),
                std::vector<tk::real>(Disc()->Lid().size(), 0.0),
                std::vector<tk::real>(Disc()->Lid().size(), 0.0) }} ),
  m_bnode( bnode ),
  m_bface( bface ),
  m_triinpoel( tk::remap( triinpoel, Disc()->Lid() ) ),
  m_bnorm(),
  m_bnormc(),
  m_dte(m_u.nunk(), 0.0),
  m_finished(0)
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
  g_dgpde[Disc()->MeshId()].numEquationDofs(m_numEqDof);

  usesAtSync = true;    // enable migration at AtSync

  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

  // Enable SDAG wait for initially building the solution vector and limiting
  if (m_initial) {
    thisProxy[ thisIndex ].wait4sol();
    if (pref) thisProxy[ thisIndex ].wait4refine();
    thisProxy[ thisIndex ].wait4smooth();
    thisProxy[ thisIndex ].wait4lim();
    thisProxy[ thisIndex ].wait4ale();
    thisProxy[ thisIndex ].wait4nod();
    thisProxy[ thisIndex ].wait4reco();
  }

  m_ghosts[thisIndex].insert(m_disc, bface, triinpoel, m_u.nunk(),
    CkCallback(CkIndex_DG::resizeSolVectors(), thisProxy[thisIndex]));

  // Query ALE mesh velocity boundary condition node lists
  Disc()->meshvelBnd( m_bface, m_bnode, m_triinpoel );

  // insert array-element into the implicit solver chare array
  if (g_inputdeck.get< tag::implicit_timestepping >()) {
    // Single-stage BDF1 for implicit solver
    m_nstage = 1;

    const auto& inpoel = myGhosts()->m_inpoel;
    // TODO: linear solver:
    //  modify CSR to handle element-based structures (or create new one)
    tk::CSR A(m_rhs.nprop(), tk::genPsup(inpoel,4,tk::genEsup(inpoel,4)));
    std::vector< tk::real > x(m_u.nunk()*m_rhs.nprop(), 0.0),
      b(m_u.nunk()*m_rhs.nprop(), 0.0);

    Disc()->ImplicitSolver()[ thisIndex ].insert(std::move(A), std::move(x),
      std::move(b), Disc()->Gid(), Disc()->Lid(), Disc()->NodeCommMap());
  }

  // global-sync to call doneInserting on m_ghosts
  auto meshid = Disc()->MeshId();
  contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
    CkCallback(CkReductionTarget(Transporter,doneInsertingGhosts),
    Disc()->Tr()) );

  // Array elements must not use the chare_objs table
  chareIdx = -1;
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
  m_rhs.resize( myGhosts()->m_nunk );
  m_rhsprev.resize( myGhosts()->m_nunk );
  m_stiffrhs.resize( myGhosts()->m_nunk );
  m_stiffrhsprev.resize( myGhosts()->m_nunk );
  for (std::size_t i=0; i<3; ++i)
    m_nodevel[i].resize( Disc()->Coord()[0].size() );
  m_dte.resize( myGhosts()->m_nunk );

  // Size communication buffer for solution and number of degrees of freedom
  for (auto& n : m_ndofc) n.resize( myGhosts()->m_bid.size() );
  for (auto& u : m_uc) u.resize( myGhosts()->m_bid.size() );
  for (auto& p : m_pc) p.resize( myGhosts()->m_bid.size() );
  for (auto& i : m_interfacec) i.resize( myGhosts()->m_bid.size() );

  // Initialize number of degrees of freedom in mesh elements
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();
  if( pref )
  {
    const auto ndofmax = g_inputdeck.get< tag::pref, tag::ndofmax >();
    m_ndof.resize( myGhosts()->m_nunk, ndofmax );
  }
  else
  {
    const auto ndof = g_inputdeck.get< tag::ndof >();
    m_ndof.resize( myGhosts()->m_nunk, ndof );
  }
  m_interface.resize( myGhosts()->m_nunk, 0 );

  // Ensure that we also have all the geometry and connectivity data
  // (including those of ghosts)
  Assert( myGhosts()->m_geoElem.nunk() == m_u.nunk(),
    "GeoElem unknowns size mismatch" );

  // Signal the runtime system that all workers have received their adjacency
  std::vector< std::size_t > meshdata{ m_initial, Disc()->MeshId() };
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

  // Determine elements inside user-defined IC box
  g_dgpde[d->MeshId()].IcBoxElems( myGhosts()->m_geoElem,
    myGhosts()->m_fd.Esuel().size()/4, m_boxelems );

  // Compute volume of user-defined box IC
  d->boxvol( {}, {}, 0 );      // punt for now

  // Query time history field output labels from all PDEs integrated
  const auto& hist_points = g_inputdeck.get< tag::history_output, tag::point >();
  if (!hist_points.empty()) {
    std::vector< std::string > histnames;
    auto n = g_dgpde[d->MeshId()].histNames();
    histnames.insert( end(histnames), begin(n), end(n) );
    d->histheader( std::move(histnames) );
  }

  // If working with IMEX-RK, Store stiff equations into m_stiffEqIdx
  if (g_inputdeck.get< tag::imex_runge_kutta >())
  {
    g_dgpde[Disc()->MeshId()].setStiffEqIdx(m_stiffEqIdx);
    g_dgpde[Disc()->MeshId()].setNonStiffEqIdx(m_nonStiffEqIdx);
  }
}

void
DG::box( tk::real v, const std::vector< tk::real >& )
// *****************************************************************************
// Receive total box IC volume and set conditions in box
//! \param[in] v Total volume within user-specified box
// *****************************************************************************
{
  auto d = Disc();

  // Store user-defined box IC volume
  d->Boxvol() = v;

  // Store previous time step and stage element volumes for GCL
  m_geoElemk = myGhosts()->m_geoElem;
  m_geoElemn = myGhosts()->m_geoElem;

  // Set initial conditions for all PDEs
  g_dgpde[d->MeshId()].initialize( myGhosts()->m_geoElem, myGhosts()->m_inpoel,
    myGhosts()->m_coord, m_boxelems, d->ElemBlockId(), m_u, d->T(),
    myGhosts()->m_fd.Esuel().size()/4 );
  g_dgpde[d->MeshId()].updatePrimitives( m_u, myGhosts()->m_geoElem, m_p,
    myGhosts()->m_fd.Esuel().size()/4, m_ndof );

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
    g_dgpde[d->MeshId()].eval_ndof( myGhosts()->m_nunk, myGhosts()->m_coord,
                  myGhosts()->m_inpoel,
                  myGhosts()->m_fd, m_u, m_p,
                  g_inputdeck.get< tag::pref, tag::indicator >(),
                  g_inputdeck.get< tag::ndof >(),
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
      std::vector< std::size_t > interface( ghostdata.size() );
      std::vector< std::size_t > ndof( ghostdata.size() );
      std::size_t j = 0;
      for(const auto& i : ghostdata) {
        Assert( i < myGhosts()->m_fd.Esuel().size()/4,
          "Sending solution ghost data" );
        tetid[j] = i;
        u[j] = m_u[i];
        prim[j] = m_p[i];
        if (pref && m_stage == 0) {
          ndof[j] = m_ndof[i];
          interface[j] = m_interface[i];
        }
        ++j;
      }
      thisProxy[ cid ].comsol( thisIndex, m_stage, tetid, u, prim, interface, ndof );
    }

  ownsol_complete();
}

void
DG::comsol( int fromch,
            std::size_t fromstage,
            const std::vector< std::size_t >& tetid,
            const std::vector< std::vector< tk::real > >& u,
            const std::vector< std::vector< tk::real > >& prim,
            const std::vector< std::size_t >& interface,
            const std::vector< std::size_t >& ndof )
// *****************************************************************************
//  Receive chare-boundary solution ghost data from neighboring chares
//! \param[in] fromch Sender chare id
//! \param[in] fromstage Sender chare time step stage
//! \param[in] tetid Ghost tet ids we receive solution data for
//! \param[in] u Solution ghost data
//! \param[in] prim Primitive variables in ghost cells
//! \param[in] interface Interface marker in ghost cells
//! \param[in] ndof Number of degrees of freedom for chare-boundary elements
//! \details This function receives contributions to the unlimited solution
//!   from fellow chares.
// *****************************************************************************
{
  Assert( u.size() == tetid.size(), "Size mismatch in DG::comsol()" );
  Assert( prim.size() == tetid.size(), "Size mismatch in DG::comsol()" );

  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

  if (pref && fromstage == 0) {
    Assert( ndof.size() == tetid.size(), "Size mismatch in DG::comsol()" );
    Assert( interface.size() == tetid.size(), "Size mismatch in DG::comsol()" );
  }

  // Find local-to-ghost tet id map for sender chare
  const auto& n = tk::cref_find( myGhosts()->m_ghost, fromch );

  for (std::size_t i=0; i<tetid.size(); ++i) {
    auto j = tk::cref_find( n, tetid[i] );
    Assert( j >= myGhosts()->m_fd.Esuel().size()/4,
      "Receiving solution non-ghost data" );
    auto b = tk::cref_find( myGhosts()->m_bid, j );
    Assert( b < m_uc[0].size(), "Indexing out of bounds" );
    Assert( b < m_pc[0].size(), "Indexing out of bounds" );
    m_uc[0][b] = u[i];
    m_pc[0][b] = prim[i];
    if (pref && fromstage == 0) {
      Assert( b < m_ndofc[0].size(), "Indexing out of bounds" );
      m_ndofc[0][b] = ndof[i];
      Assert( b < m_interfacec[0].size(), "Indexing out of bounds" );
      m_interfacec[0][b] = interface[i];
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
  evalSolution( *Disc(), inpoel, coord, addedTets, m_ndof, m_u, m_p,
    m_uElemfields, m_pElemfields, m_uNodefields, m_pNodefields );

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

  ownnod_complete( c, addedTets );
}

void DG::p_refine()
// *****************************************************************************
// Add the protective layer for ndof refinement
// *****************************************************************************
{
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

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
    if (pref && m_stage == 0) {
      m_ndof[ b.first ] = m_ndofc[0][ b.second ];
      m_interface[ b.first ] = m_interfacec[0][ b.second ];
    }
  }

  if (pref && m_stage==0) refine_ndof();

  if (!pref) {
    // if p-refinement is not configured, proceed directly to reconstructions
    reco();
  }
  else {
    // if p-refinement is configured, do refine-smoothing before reconstruction
    if (myGhosts()->m_sendGhost.empty())
      comrefine_complete();
    else
      for(const auto& [cid, ghostdata] : myGhosts()->m_sendGhost) {
        std::vector< std::size_t > tetid( ghostdata.size() );
        std::vector< std::vector< tk::real > > u( ghostdata.size() ),
                                               prim( ghostdata.size() );
        std::vector< std::size_t > ndof( ghostdata.size() );
        std::size_t j = 0;
        for(const auto& i : ghostdata) {
          Assert( i < myGhosts()->m_fd.Esuel().size()/4, "Sending refined ndof  "
            "data" );
          tetid[j] = i;
          if (pref && m_stage == 0) ndof[j] = m_ndof[i];
          ++j;
        }
        thisProxy[ cid ].comrefine( thisIndex, tetid, ndof );
      }

    ownrefine_complete();
  }
}

void
DG::comrefine( int fromch,
               const std::vector< std::size_t >& tetid,
               const std::vector< std::size_t >& ndof )
// *****************************************************************************
//  Receive chare-boundary ghost data from neighboring chares
//! \param[in] fromch Sender chare id
//! \param[in] tetid Ghost tet ids we receive solution data for
//! \param[in] ndof Number of degrees of freedom for chare-boundary elements
//! \details This function receives contributions to the refined ndof data
//!   from fellow chares.
// *****************************************************************************
{
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

  if (pref && m_stage == 0)
    Assert( ndof.size() == tetid.size(), "Size mismatch in DG::comrefine()" );

  // Find local-to-ghost tet id map for sender chare
  const auto& n = tk::cref_find( myGhosts()->m_ghost, fromch );

  for (std::size_t i=0; i<tetid.size(); ++i) {
    auto j = tk::cref_find( n, tetid[i] );
    Assert( j >= myGhosts()->m_fd.Esuel().size()/4,
      "Receiving solution non-ghost data" );
    auto b = tk::cref_find( myGhosts()->m_bid, j );
    if (pref && m_stage == 0) {
      Assert( b < m_ndofc[1].size(), "Indexing out of bounds" );
      m_ndofc[1][b] = ndof[i];
    }
  }

  // if we have received all solution ghost contributions from neighboring
  // chares (chares we communicate along chare-boundary faces with), and
  // contributed our solution to these neighbors, proceed to limiting
  if (++m_nrefine == myGhosts()->m_sendGhost.size()) {
    m_nrefine = 0;
    comrefine_complete();
  }
}

void
DG::smooth()
// *****************************************************************************
// Smooth the refined ndof distribution
// *****************************************************************************
{
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

  for (const auto& b : myGhosts()->m_bid) {
    if (pref && m_stage == 0)
      m_ndof[ b.first ] = m_ndofc[1][ b.second ];
  }

  if (pref && m_stage==0) smooth_ndof();

  if (myGhosts()->m_sendGhost.empty())
    comsmooth_complete();
  else
    for(const auto& [cid, ghostdata] : myGhosts()->m_sendGhost) {
      std::vector< std::size_t > tetid( ghostdata.size() );
      std::vector< std::size_t > ndof;
      std::size_t j = 0;
      for(const auto& i : ghostdata) {
        Assert( i < myGhosts()->m_fd.Esuel().size()/4, "Sending ndof data" );
        tetid[j] = i;
        if (pref && m_stage == 0) ndof.push_back( m_ndof[i] );
        ++j;
      }
      thisProxy[ cid ].comsmooth( thisIndex, tetid, ndof );
    }

  ownsmooth_complete();
}

void
DG::comsmooth( int fromch,
               const std::vector< std::size_t >& tetid,
               const std::vector< std::size_t >& ndof )
// *****************************************************************************
//  Receive chare-boundary ghost data from neighboring chares
//! \param[in] fromch Sender chare id
//! \param[in] tetid Ghost tet ids we receive solution data for
//! \param[in] ndof Number of degrees of freedom for chare-boundary elements
//! \details This function receives contributions to the smoothed ndof data
//!   from fellow chares.
// *****************************************************************************
{
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

  if (pref && m_stage == 0)
    Assert( ndof.size() == tetid.size(), "Size mismatch in DG::comsmooth()" );

  const auto& n = tk::cref_find( myGhosts()->m_ghost, fromch );

  for (std::size_t i=0; i<tetid.size(); ++i) {
    auto j = tk::cref_find( n, tetid[i] );
    Assert( j >= myGhosts()->m_fd.Esuel().size()/4, "Receiving ndof data" );
    auto b = tk::cref_find( myGhosts()->m_bid, j );
    if (pref && m_stage == 0) {
      Assert( b < m_ndofc[2].size(), "Indexing out of bounds" );
      m_ndofc[2][b] = ndof[i];
    }
  }

  if (++m_nsmooth == myGhosts()->m_sendGhost.size()) {
    m_nsmooth = 0;
    comsmooth_complete();
  }
}

void
DG::reco()
// *****************************************************************************
// Compute reconstructions
// *****************************************************************************
{
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();
  const auto rdof = g_inputdeck.get< tag::rdof >();

  // Combine own and communicated contributions of unreconstructed solution and
  // degrees of freedom in cells (if p-adaptive)
  if (pref) {
    for (const auto& b : myGhosts()->m_bid) {
      if (m_stage == 0) {
        m_ndof[ b.first ] = m_ndofc[2][ b.second ];
      }
    }
  }

  auto d = Disc();
  if (pref && m_stage==0) {
    g_dgpde[d->MeshId()].resetAdapSol( myGhosts()->m_fd, m_u, m_p, m_ndof );
  }

  if (rdof > 1)
    // Reconstruct second-order solution and primitive quantities
    g_dgpde[d->MeshId()].reconstruct( d->T(), myGhosts()->m_geoFace,
      myGhosts()->m_geoElem,
      myGhosts()->m_fd, myGhosts()->m_esup, myGhosts()->m_inpoel,
      myGhosts()->m_coord, m_u, m_p, pref, m_ndof );

  // Send reconstructed solution to neighboring chares
  if (myGhosts()->m_sendGhost.empty())
    comreco_complete();
  else
    for(const auto& [cid, ghostdata] : myGhosts()->m_sendGhost) {
      std::vector< std::size_t > tetid( ghostdata.size() );
      std::vector< std::vector< tk::real > > u( ghostdata.size() ),
                                             prim( ghostdata.size() );
      std::size_t j = 0;
      for(const auto& i : ghostdata) {
        Assert( i < myGhosts()->m_fd.Esuel().size()/4, "Sending reconstructed ghost "
          "data" );
        tetid[j] = i;
        u[j] = m_u[i];
        prim[j] = m_p[i];
        ++j;
      }
      thisProxy[ cid ].comreco( thisIndex, tetid, u, prim );
    }

  ownreco_complete();
}

void
DG::comreco( int fromch,
             const std::vector< std::size_t >& tetid,
             const std::vector< std::vector< tk::real > >& u,
             const std::vector< std::vector< tk::real > >& prim )
// *****************************************************************************
//  Receive chare-boundary reconstructed ghost data from neighboring chares
//! \param[in] fromch Sender chare id
//! \param[in] tetid Ghost tet ids we receive solution data for
//! \param[in] u Reconstructed high-order solution
//! \param[in] prim Limited high-order primitive quantities
//! \details This function receives contributions to the reconstructed solution
//!   from fellow chares.
// *****************************************************************************
{
  Assert( u.size() == tetid.size(), "Size mismatch in DG::comreco()" );
  Assert( prim.size() == tetid.size(), "Size mismatch in DG::comreco()" );

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
  if (++m_nreco == myGhosts()->m_sendGhost.size()) {
    m_nreco = 0;
    comreco_complete();
  }
}

void
DG::lim()
// *****************************************************************************
// Compute limiter function
// *****************************************************************************
{
  auto d = Disc();
  auto gid = d->Gid();
  auto bid = d->Bid();
  const auto rdof = g_inputdeck.get< tag::rdof >();
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

  // Combine own and communicated contributions of unlimited solution, and
  // if a p-adaptive algorithm is used, degrees of freedom in cells
  for (const auto& [boundary, localtet] : myGhosts()->m_bid) {
    Assert( m_uc[1][localtet].size() == m_u.nprop(), "ncomp size mismatch" );
    Assert( m_pc[1][localtet].size() == m_p.nprop(), "ncomp size mismatch" );
    for (std::size_t c=0; c<m_u.nprop(); ++c) {
      m_u(boundary,c) = m_uc[1][localtet][c];
    }
    for (std::size_t c=0; c<m_p.nprop(); ++c) {
      m_p(boundary,c) = m_pc[1][localtet][c];
    }
  }

  if (rdof > 1) {
    g_dgpde[d->MeshId()].limit( d->T(), pref, myGhosts()->m_geoFace,
              myGhosts()->m_geoElem, myGhosts()->m_fd, myGhosts()->m_esup,
              myGhosts()->m_inpoel, myGhosts()->m_coord, m_ndof, d->Gid(),
              d->Bid(), m_mtInv, m_u, m_p, m_shockmarker );

    if (g_inputdeck.get< tag::limsol_projection >())
      g_dgpde[d->MeshId()].CPL(m_p, myGhosts()->m_geoElem,
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
      std::vector< std::size_t > ndof;
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
DG::refine_ndof()
// *****************************************************************************
//  p-refine all elements that are adjacent to p-refined elements
//! \details This function p-refines all the neighbors of an element that has
//!   been p-refined as a result of an error indicator.
// *****************************************************************************
{
  auto d = Disc();
  const auto& coord = d->Coord();
  const auto& inpoel = d->Inpoel();
  const auto npoin = coord[0].size();
  const auto nelem = myGhosts()->m_fd.Esuel().size()/4;
  std::vector<std::size_t> node_ndof(npoin, 1);

  // Mark the max ndof for each node and store in node_ndof
  for(std::size_t ip=0; ip<npoin; ip++)
  {
    const auto& pesup = tk::cref_find(myGhosts()->m_esup, ip);
    for(auto er : pesup)
      node_ndof[ip] = std::max(m_ndof[er], node_ndof[ip]);
  }

  for(std::size_t e = 0; e < nelem; e++)
  {
    // Find if any node of this element has p1/p2 ndofs
    std::size_t counter_p2(0);
    std::size_t counter_p1(0);
    for(std::size_t inode = 0; inode < 4; inode++)
    {
      auto node = inpoel[4*e+inode];
      if(node_ndof[node] == 10)
        counter_p2++;
      else if (node_ndof[node] == 4)
        counter_p1++;
    }

    // If there is at least one node with p1/p2 ndofs, all of the elements
    // around this node are refined to p1/p2.
    if(counter_p2 > 0 && m_ndof[e] < 10)
    {
      if(m_ndof[e] == 4)
        m_ndof[e] = 10;
      if(m_ndof[e] == 1)
        m_ndof[e] = 4;
    }
    else if(counter_p1 > 0 && m_ndof[e] < 4)
      m_ndof[e] = 4;
  }
}

void DG::smooth_ndof()
// *****************************************************************************
//  Smooth the refined ndof distribution to avoid zigzag refinement
// *****************************************************************************
{
  auto d = Disc();
  const auto& inpoel = d->Inpoel();
  const auto& coord = d->Coord();
  const auto npoin = coord[0].size();
  const auto nelem = myGhosts()->m_fd.Esuel().size()/4;
  std::vector<std::size_t> node_ndof(npoin, 1);

  // Mark the max ndof for each node and store in node_ndof
  for(std::size_t ip=0; ip<npoin; ip++)
  {
    const auto& pesup = tk::cref_find(myGhosts()->m_esup, ip);
    for(auto er : pesup)
      node_ndof[ip] = std::max(m_ndof[er], node_ndof[ip]);
  }

  for(std::size_t e = 0; e < nelem; e++)
  {
    // Find if any node of this element has p1/p2 ndofs
    std::size_t counter_p2(0);
    std::size_t counter_p1(0);
    for(std::size_t inode = 0; inode < 4; inode++)
    {
      auto node = inpoel[4*e+inode];
      if(node_ndof[node] == 10)
        counter_p2++;
      else if (node_ndof[node] == 4)
        counter_p1++;
    }

    // If all the nodes in the element are p1/p2, this element is refined to
    // p1/p2.
    if(counter_p2 == 4 && m_ndof[e] == 4)
      m_ndof[e] = 10;
    else if(counter_p1 == 4 && m_ndof[e] == 1)
      m_ndof[e] = 4;
  }
}

void
DG::comlim( int fromch,
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
  Assert( u.size() == tetid.size(), "Size mismatch in DG::comlim()" );
  Assert( prim.size() == tetid.size(), "Size mismatch in DG::comlim()" );

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
DG::updateChareBoundaryGeoFace()
// *****************************************************************************
// Recompute chare-boundary face geometry after ALE ghost updates arrive
// *****************************************************************************
{
  auto d = Disc();
  auto g = myGhosts();
  const auto& esuf = g->m_fd.Esuf();

  for (std::size_t f = g->m_fd.Nipfac(); f < esuf.size()/2; ++f) {
    std::size_t el = static_cast< std::size_t >( esuf[2*f] );
    tk::UnsMesh::Face t{{
      d->Gid()[ g->m_fd.Inpofa()[3*f+2] ],
      d->Gid()[ g->m_fd.Inpofa()[3*f+1] ],
      d->Gid()[ g->m_fd.Inpofa()[3*f+0] ]
    }};
    std::array< std::size_t, 2 > id{{ f, el }};
    g->addGeoFace( t, id );
  }
}

void
DG::bnorm()
// *****************************************************************************
//  Compute boundary point normals for mesh velocity symmetry BCs
// *****************************************************************************
{
  auto d = Disc();

  std::unordered_map< int, std::unordered_set< std::size_t > > bcnodes;
  for (const auto& s : g_inputdeck.get< tag::ale, tag::symmetry >()) {
    auto k = m_bface.find(static_cast<int>(s));
    if (k != end(m_bface)) {
      auto& n = bcnodes[ k->first ];
      for (auto f : k->second) {
        n.insert( m_triinpoel[f*3+0] );
        n.insert( m_triinpoel[f*3+1] );
        n.insert( m_triinpoel[f*3+2] );
      }
    }
  }

  m_bnorm = cg::bnorm( m_bface, m_triinpoel, d->Coord(), d->Gid(), bcnodes );

  // Send our nodal normal contributions to neighbor chares
  if (d->NodeCommMap().empty())
    comnorm_complete();
  else
    for (const auto& [neighborchare, sharednodes] : d->NodeCommMap()) {
      decltype(m_bnorm) exp;
      for (auto i : sharednodes) {
        for (const auto& [s,norms] : m_bnorm) {
          auto j = norms.find(i);
          if (j != end(norms)) exp[s][i] = j->second;
        }
      }
      thisProxy[ neighborchare ].comnorm( exp );
    }

  ownnorm_complete();
}

void
DG::comnorm( const std::unordered_map< int,
  std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& innorm )
// *****************************************************************************
// Receive boundary point normals on chare-boundaries
//! \param[in] innorm Incoming partial sums of boundary point normal
//!   contributions to normals (first 3 components), inverse distance squared
//!   (4th component), associated to side set ids
// *****************************************************************************
{
  // Buffer up incoming boundary-point normal vector contributions
  for (const auto& [s,norms] : innorm) {
    auto& bnorms = m_bnormc[s];
    for (const auto& [p,n] : norms) {
      auto& bnorm = bnorms[p];
      bnorm[0] += n[0];
      bnorm[1] += n[1];
      bnorm[2] += n[2];
      bnorm[3] += n[3];
    }
  }

  if (++m_nbnorm == Disc()->NodeCommMap().size()) {
    m_nbnorm = 0;
    comnorm_complete();
  }
}

void
DG::normfinal()
// *****************************************************************************
//  Finish computing boundary point normals
// *****************************************************************************
{
  const auto& lid = Disc()->Lid();

  // Combine own and communicated contributions to boundary point normals
  for (const auto& [s,norms] : m_bnormc) {
    auto& bnorms = m_bnorm[s];
    for (const auto& [p,n] : norms) {
      auto& norm = bnorms[p];
      norm[0] += n[0];
      norm[1] += n[1];
      norm[2] += n[2];
      norm[3] += n[3];
    }
  }
  tk::destroy( m_bnormc );

  // Divide summed point normals by the sum of inverse distance squared
  for (auto& [s,norms] : m_bnorm)
    for (auto& [p,n] : norms) {
      n[0] /= n[3];
      n[1] /= n[3];
      n[2] /= n[3];
      Assert( (n[0]*n[0] + n[1]*n[1] + n[2]*n[2] - 1.0) <
              1.0e+3*std::numeric_limits< tk::real >::epsilon(),
              "Non-unit normal" );
    }

  // Replace global->local ids associated to boundary point normals
  decltype(m_bnorm) bnorm;
  for (auto& [s,norms] : m_bnorm) {
    auto& bnorms = bnorm[s];
    for (auto&& [g,n] : norms)
      bnorms[ tk::cref_find(lid,g) ] = std::move(n);
  }
  m_bnorm = std::move(bnorm);

  meshvelstart();
}

void
DG::ALEComm()
// *****************************************************************************
// Perform ALE mesh update and communicate updated ghost mesh data
// *****************************************************************************
{
  auto d = Disc();
  const auto is_ale = g_inputdeck.get< tag::ale, tag::ale >();

  if (!is_ale) {
    ownale_complete();
    comale_complete();
    return;
  }

  if (m_stage == 0) {
    d->UpdateCoordn();
    m_geoElemn = myGhosts()->m_geoElem;
  }

  // Advance owned mesh coordinates and mirror them into the extended array.
  const auto& meshvel = d->meshvel();
  auto& coord = myGhosts()->m_coord;
  auto& disc_coord = d->Coord();
  const auto dc_size = disc_coord[0].size();
  for (auto j : g_inputdeck.get< tag::ale, tag::mesh_motion_directions >()) {
    for (std::size_t i=0; i<dc_size; ++i) {
      auto x = rkcoef[0][m_stage] * d->Coordn()[j][i]
        + rkcoef[1][m_stage] * ( disc_coord[j][i] + d->Dt() * meshvel(i,j) );
      disc_coord[j][i] = x;
      coord[j][i] = x;
    }
  }

  // Store element volumes at previous stage for GCL consistent RK
  m_geoElemk = myGhosts()->m_geoElem;

  // Recompute internal + physical boundary face geometry from the updated mesh.
  auto gf_temp = tk::genGeoFaceTri( myGhosts()->m_fd.Nipfac(),
    myGhosts()->m_fd.Inpofa(), myGhosts()->m_coord );
  for (std::size_t f=0; f<myGhosts()->m_fd.Nipfac(); ++f)
    for (std::size_t i=0; i<gf_temp.nprop(); ++i)
      myGhosts()->m_geoFace(f,i) = gf_temp(f,i);

  // Recompute element geometries for owned elements only.
  auto ge_temp = tk::genGeoElemTet( d->Inpoel(), disc_coord );
  for (std::size_t e=0; e<ge_temp.nunk(); ++e)
    for (std::size_t i=0; i<ge_temp.nprop(); ++i)
      myGhosts()->m_geoElem(e,i) = ge_temp(e,i);

  if (myGhosts()->m_sendGhost.empty()) {
    updateChareBoundaryGeoFace();
    comale_complete();
  } else {
    for (const auto& [cid, ghostdata] : myGhosts()->m_sendGhost) {
      std::vector< std::size_t > tetid( ghostdata.size() );
      std::vector< std::vector< tk::real > > geoElem( ghostdata.size() );
      std::vector< std::array< tk::real, 3 > > coordg(
        ghostdata.size(), std::array< tk::real, 3 >{{ 0.0, 0.0, 0.0 }} );
      const auto sendnode = myGhosts()->m_sendChBndNode.find( cid );

      std::size_t j = 0;
      for (const auto& e : ghostdata) {
        Assert( e < myGhosts()->m_fd.Esuel().size()/4,
          "Sending ALE ghost data" );
        tetid[j] = e;
        geoElem[j] = myGhosts()->m_geoElem[e];

        if (sendnode != end(myGhosts()->m_sendChBndNode)) {
          const auto n = sendnode->second.find( e );
          if (n != end(sendnode->second)) {
            coordg[j] = {{ coord[0][n->second],
                           coord[1][n->second],
                           coord[2][n->second] }};
          }
        }
        ++j;
      }

      thisProxy[ cid ].comale( thisIndex, tetid, geoElem, coordg );
    }
  }

  ownale_complete();
}

void
DG::comale( int fromch,
            const std::vector< std::size_t >& tetid,
            const std::vector< std::vector< tk::real > >& geoElem,
            const std::vector< std::array< tk::real, 3 > >& coord )
// *****************************************************************************
//  Receive updated ALE ghost mesh data from neighboring chares
//! \param[in] fromch Sender chare id
//! \param[in] tetid Ghost tet ids we receive ALE mesh data for
//! \param[in] geoElem Updated ghost-element geometry
//! \param[in] coord Updated off-face coordinates for face ghosts
// *****************************************************************************
{
  Assert( geoElem.size() == tetid.size(), "Size mismatch in DG::comale()" );
  Assert( coord.size() == tetid.size(), "Size mismatch in DG::comale()" );

  const auto& ghost = tk::cref_find( myGhosts()->m_ghost, fromch );
  const auto recvnode = myGhosts()->m_recvChBndNode.find( fromch );

  for (std::size_t i=0; i<tetid.size(); ++i) {
    auto j = tk::cref_find( ghost, tetid[i] );
    Assert( j >= myGhosts()->m_fd.Esuel().size()/4,
      "Receiving ALE non-ghost data" );
    Assert( geoElem[i].size() == myGhosts()->m_geoElem.nprop(),
      "Geometry size mismatch in DG::comale()" );

    for (std::size_t c=0; c<geoElem[i].size(); ++c)
      myGhosts()->m_geoElem(j,c) = geoElem[i][c];

    if (recvnode != end(myGhosts()->m_recvChBndNode)) {
      const auto n = recvnode->second.find( tetid[i] );
      if (n != end(recvnode->second)) {
        auto p = n->second;
        Assert( p < myGhosts()->m_coord[0].size(),
          "Indexing out of extended ALE ghost coordinates" );
        myGhosts()->m_coord[0][p] = coord[i][0];
        myGhosts()->m_coord[1][p] = coord[i][1];
        myGhosts()->m_coord[2][p] = coord[i][2];
      }
    }
  }

  if (++m_nale == myGhosts()->m_sendGhost.size()) {
    m_nale = 0;
    updateChareBoundaryGeoFace();
    comale_complete();
  }
}

void
DG::dt()
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
      m_u(b.first,c) = m_uc[2][b.second][c];
    }
    for (std::size_t c=0; c<m_p.nprop(); ++c) {
      m_p(b.first,c) = m_pc[2][b.second][c];
    }
  }

  auto mindt = std::numeric_limits< tk::real >::max();

  if (m_stage == 0)
  {
    auto const_dt = g_inputdeck.get< tag::dt >();
    auto eps = std::numeric_limits< tk::real >::epsilon();

    // use constant dt if configured
    if (std::abs(const_dt) > eps) {

      mindt = const_dt;

    } else {      // compute dt based on CFL

      // find the minimum dt across all PDEs integrated
      auto eqdt =
        g_dgpde[d->MeshId()].dt( myGhosts()->m_coord, myGhosts()->m_inpoel,
          myGhosts()->m_fd,
          myGhosts()->m_geoFace, myGhosts()->m_geoElem, m_ndof, m_u, m_p,
          myGhosts()->m_fd.Esuel().size()/4, m_dte );
      if (eqdt < mindt) mindt = eqdt;

      // time-step suppression for unsteady problems
      tk::real coeff(1.0);
      auto ramp_steps = g_inputdeck.get< tag::cfl_ramping_steps >();
      if (g_inputdeck.get< tag::cfl_ramping >() && d->It() < ramp_steps)
        coeff = 1.0/static_cast< tk::real >(ramp_steps)
          * static_cast< tk::real >(d->It()+1);

      if (g_inputdeck.get< tag::steady_state >()) {
        for (auto& edt : m_dte) edt *= coeff * g_inputdeck.get< tag::cfl >();
      }

      mindt *= coeff * g_inputdeck.get< tag::cfl >();

      // time-step restriction based on max volume change
      auto mindtv = std::numeric_limits< tk::real >::max();
      auto dvcfl = g_inputdeck.get< tag::ale, tag::dvcfl >();
      if (d->Dtn() > 1e-12 && dvcfl > 0.0) {
        for (std::size_t e=0; e<myGhosts()->m_nunk; ++e) {
          auto dt_v = dvcfl *
            d->Dtn() * std::min(m_geoElemn(e,0), myGhosts()->m_geoElem(e,0))
            / (std::abs(m_geoElemn(e,0)-myGhosts()->m_geoElem(e,0)) + 1.0e-12);
          mindtv = std::min(mindtv, dt_v);
        }
      }

      mindt = std::min(mindt, mindtv);
    }
  }
  else
  {
    mindt = d->Dt();
  }

  // Set up the reduction target for finding minimum dt across chares.
  // 1. If implicit solver is used, first invoke the solver object via
  // appropriate entry methods, and then proceed to solve.
  // 2. If explicit, directly proceed to solve.
  CkCallback minDtDone;
  if (!g_inputdeck.get< tag::implicit_timestepping >())
    minDtDone = CkCallback(CkReductionTarget(DG,solve), thisProxy);
  else
    minDtDone = CkCallback(CkReductionTarget(DG,initializeLinearSystem),
      thisProxy);

  // Contribute to minimum dt across all chares then advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double, minDtDone );
}

void
DG::initializeLinearSystem( tk::real newdt )
// *****************************************************************************
// Initialize the linear solver via the interface BiCG::init()
//! \param[in] newdt Size of this new time step
// *****************************************************************************
{
  auto d = Disc();

  // Set new time step size
  if (m_stage == 0) d->setdt( newdt );

  // Initialize linear solver, and route to solveLinearSystem()
  // TODO: linear solver:
  // 1. jacobian computation (call to e.g. g_dgpde[d->MeshId()].computeJacobian)
  // 2. the following call is just a stand-in/example- verify correctness
  d->ImplicitSolver()[ thisIndex ].init( m_u.flat(), {}, {}, 1,
    CkCallback(CkIndex_DG::solveLinearSystem(), thisProxy[thisIndex]) );
}

void
DG::solveLinearSystem()
// *****************************************************************************
// Solve the linear system via the interface BiCG::solve()
// *****************************************************************************
{
  auto d = Disc();

  // Get new time step size to pass along to solve()
  auto dt = d->Dt();

  // Solve linear system, and route to solve()
  // TODO: linear solver:
  //  the following call is just a stand-in/example- verify correctness
  d->ImplicitSolver()[ thisIndex ].solve(
     g_inputdeck.get< tag::ale, tag::maxit >(),
     g_inputdeck.get< tag::residual >(),
     CkCallback(CkIndex_DG::solve(dt), thisProxy[thisIndex]) );
}

void
DG::solve( tk::real newdt )
// *****************************************************************************
// Compute right-hand side of discrete transport equations
//! \param[in] newdt Size of this new time step
// *****************************************************************************
{
  const auto pref = g_inputdeck.get< tag::pref, tag::pref >();

  // Enable SDAG wait for building the solution vector during the next stage
  thisProxy[ thisIndex ].wait4sol();
  if (pref) thisProxy[ thisIndex ].wait4refine();
  thisProxy[ thisIndex ].wait4smooth();
  thisProxy[ thisIndex ].wait4reco();
  thisProxy[ thisIndex ].wait4lim();
  thisProxy[ thisIndex ].wait4ale();
  thisProxy[ thisIndex ].wait4nod();

  auto d = Disc();
  const auto rdof = g_inputdeck.get< tag::rdof >();
  const auto ndof = g_inputdeck.get< tag::ndof >();
  const auto neq = m_u.nprop()/rdof;

  // Set new time step size. If implicit solver, time step has already been
  // set in initializeImplicitSystem()
  if (m_stage == 0 && !g_inputdeck.get< tag::implicit_timestepping >())
    d->setdt( newdt );

  // Update Un
  if (m_stage == 0) m_un = m_u;

  // Explicit or IMEX
  const auto imex_runge_kutta = g_inputdeck.get< tag::imex_runge_kutta >();
  const auto implicit_ts = g_inputdeck.get< tag::implicit_timestepping >();

  // physical time at time-stage for computing exact source terms
  tk::real physT(d->T());
  if (m_stage == 1) {
    physT += d->Dt();
  }
  else if (m_stage == 2) {
    physT += 0.5*d->Dt();
  }

  if (imex_runge_kutta) {
    if (m_stage == 0)
    {
      // Initialize m_stiffrhs to zero
      m_stiffrhs.fill(0.0);
      m_stiffrhsprev.fill(0.0);
    }
  }

  if (!imex_runge_kutta || m_stage < m_nstage-1) {
    if (imex_runge_kutta && m_stage < m_nstage-1) m_rhsprev = m_rhs;
    g_dgpde[d->MeshId()].rhs( physT, pref, myGhosts()->m_geoFace,
      myGhosts()->m_geoElem, myGhosts()->m_fd, myGhosts()->m_inpoel, m_boxelems,
      myGhosts()->m_coord, m_u, m_p, d->meshvel(), m_ndof, d->Dt(), m_rhs );
  }

  if (imex_runge_kutta) {
    // Implicit-Explicit time-stepping using RK3 to discretize time-derivative
    DG::imex_integrate();
  }
  else if (implicit_ts) {
    // Implicit time-stepping using BDF1 to discretize time-derivative
    DG::BDF1_integrate();
  }
  else {
    // Explicit time-stepping using RK3 to discretize time-derivative
    const auto steady = g_inputdeck.get< tag::steady_state >();
    for(std::size_t e=0; e<myGhosts()->m_nunk; ++e) {

      // Stage-wise volumes for GCL consistent RK
      auto vole = myGhosts()->m_geoElem(e,0);
      auto vole_n = m_geoElemn(e,0);
      auto vole_k = m_geoElemk(e,0);

      auto dte = d -> Dt();
      if (steady) dte = m_dte[e];
      for(std::size_t c=0; c<neq; ++c)
      {
        for (std::size_t k=0; k<m_numEqDof[c]; ++k)
        {
          if(k < m_ndof[e]) {
            auto rmark = c*rdof+k;
            auto mark = c*ndof+k;

            auto mm_i = vole * mass_dubiner[k];
            auto mm_n = vole_n * mass_dubiner[k];
            auto mm_k = vole_k * mass_dubiner[k];
            m_u(e, rmark) = (
              rkcoef[0][m_stage] * mm_n * m_un(e, rmark)
              + rkcoef[1][m_stage] * ( mm_k * m_u(e, rmark)
              + dte * m_rhs(e, mark) )
              ) / mm_i;
            if(fabs(m_u(e, rmark)) < 1e-16)
              m_u(e, rmark) = 0;
          }
        }
      }
    }
  }

  for(std::size_t e=0; e<myGhosts()->m_nunk; ++e)
    for(std::size_t c=0; c<neq; ++c)
    {
      // zero out unused/reconstructed dofs of equations using reduced dofs
      // (see DGMultiMat::numEquationDofs())
      if (m_numEqDof[c] < rdof) {
        for (std::size_t k=m_numEqDof[c]; k<rdof; ++k)
        {
          auto rmark = c*rdof+k;
          m_u(e, rmark) = 0.0;
        }
      }
    }

  // Evolve damage
  g_dgpde[d->MeshId()].evolveDamage( d->Dt(), myGhosts()->m_geoElem, m_u,
      m_p, myGhosts()->m_fd.Esuel().size()/4 );

  // Update primitives based on the evolved solution
  g_dgpde[d->MeshId()].updateInterfaceCells( m_u,
    myGhosts()->m_fd.Esuel().size()/4, m_ndof, m_interface );
  g_dgpde[d->MeshId()].updatePrimitives( m_u, myGhosts()->m_geoElem, m_p,
    myGhosts()->m_fd.Esuel().size()/4, m_ndof );
  if (!g_inputdeck.get< tag::accuracy_test >()) {
    g_dgpde[d->MeshId()].cleanTraceMaterial( physT, myGhosts()->m_geoElem, m_u,
      m_p, myGhosts()->m_fd.Esuel().size()/4 );
  }

  if (m_stage < m_nstage-1) {

    // continue with next time step stage
    stage();

  } else {

    // Increase number of iterations and physical time
    d->next();

    // Compute diagnostics, e.g., residuals
    auto diag_computed = m_diag.compute( *d,
      m_u.nunk()-myGhosts()->m_fd.Esuel().size()/4, myGhosts()->m_geoElem,
      m_ndof, m_u, m_un );

    // Continue to mesh refinement (if configured)
    if (!diag_computed) refine( std::vector< tk::real >( m_u.nprop(), 1.0 ) );

  }
}

void
DG::refine( const std::vector< tk::real >& l2res )
// *****************************************************************************
// Optionally refine/derefine mesh
//! \param[in] l2res L2-norms of the residual for each scalar component
//!   computed across the whole problem
// *****************************************************************************
{
  auto d = Disc();

  // Assess convergence for steady state
  const auto steady = g_inputdeck.get< tag::steady_state >();
  const auto residual = g_inputdeck.get< tag::residual >();
  const auto rc = g_inputdeck.get< tag::rescomp >() - 1;

  bool converged(false);
  if (steady) converged = l2res[rc] < residual;

  // this is the last time step if max time of max number of time steps
  // reached or the residual has reached its convergence criterion
  if (d->finished() or converged) m_finished = 1;

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
  const std::map< int, std::vector< std::size_t > >& bnode,
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
  m_un.resize( nelem );
  m_rhs.resize( nelem );
  m_rhsprev.resize( nelem );
  m_stiffrhs.resize( nelem );
  m_stiffrhsprev.resize( nelem );
  for (std::size_t i=0; i<3; ++i) m_nodevel[i].resize( coord[0].size() );

  myGhosts()->m_fd = FaceData( myGhosts()->m_inpoel, bface,
    tk::remap(triinpoel,d->Lid()) );

  m_bnode = bnode;
  m_bface = bface;
  m_triinpoel = tk::remap( triinpoel, d->Lid() );
  d->meshvelBnd( m_bface, m_bnode, m_triinpoel );

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
DG::fieldOutput() const
// *****************************************************************************
// Decide wether to output field data
//! \return True if field data is output in this step
// *****************************************************************************
{
  auto d = Disc();

  // Output field data
  return d->fielditer() or d->fieldtime() or d->fieldrange() or m_finished;
}

bool
DG::refinedOutput() const
// *****************************************************************************
// Decide if we write field output using a refined mesh
//! \return True if field output will use a refined mesh
// *****************************************************************************
{
  return g_inputdeck.get< tag::field_output, tag::refined >() &&
         g_inputdeck.get< tag::scheme >() != ctr::SchemeType::DGP0;
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
        m_uNodefields(p,i) /= n;
      for (std::size_t i=0; i<m_pNodefields.nprop(); ++i)
        m_pNodefields(p,i) /= n;
    }
  }

  // Collect field output from numerical solution requested by user
  auto elemfields = numericFieldOutput( m_uElemfields, tk::Centering::ELEM,
    g_dgpde[Disc()->MeshId()].OutVarFn(), m_pElemfields );
  auto nodefields = numericFieldOutput( m_uNodefields, tk::Centering::NODE,
    g_dgpde[Disc()->MeshId()].OutVarFn(), m_pNodefields );

  // Collect field output from analytical solutions (if exist)
  const auto& coord = m_outmesh.coord;
  auto geoElem = tk::genGeoElemTet( inpoel, coord );
  auto t = Disc()->T();
  analyticFieldOutput( g_dgpde[d->MeshId()], tk::Centering::ELEM,
    geoElem.extract_comp(1), geoElem.extract_comp(2), geoElem.extract_comp(3),
    t, elemfields );
  analyticFieldOutput( g_dgpde[d->MeshId()], tk::Centering::NODE, coord[0],
    coord[1], coord[2], t, nodefields );

  // Add adaptive indicator array to element-centered field output
  if (g_inputdeck.get< tag::pref, tag::pref >()) {
    std::vector< tk::real > ndof( begin(m_ndof), end(m_ndof) );
    ndof.resize( nelem );
    for(std::size_t k = 0; k < nelem; k++) {
      // Mark the cell with THINC reconstruction as 0 for output
      if(m_interface[k] == 1) ndof[k] = 0;
    }
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

  // Compute plastic deformation averaged for all materials
  std::vector< tk::real > plasticDeformation(nelem);
  g_dgpde[d->MeshId()].computePlasticDeformation(nelem, m_u, m_p, plasticDeformation);
  for (const auto& [child,parent] : addedTets)
    plasticDeformation[child] = 0.0;
  if (plasticDeformation.size() > 0) elemfields.push_back( plasticDeformation );

  // Query fields names requested by user
  auto elemfieldnames = numericFieldNames( tk::Centering::ELEM );
  auto nodefieldnames = numericFieldNames( tk::Centering::NODE );

  // Collect field output names for analytical solutions
  analyticFieldNames( g_dgpde[d->MeshId()], tk::Centering::ELEM, elemfieldnames );
  analyticFieldNames( g_dgpde[d->MeshId()], tk::Centering::NODE, nodefieldnames );

  if (g_inputdeck.get< tag::pref, tag::pref >()) {
    elemfieldnames.push_back( "NDOF" );
  }

  elemfieldnames.push_back( "shock_marker" );

  if (plasticDeformation.size() > 0)
    elemfieldnames.push_back( "plastic_deformation" );

  //! Lambda to put in a field for output if not empty
  auto add_node_field = [&]( const auto& name, const auto& field ){
    if (not field.empty()) {
      nodefieldnames.push_back( name );
      nodefields.push_back( field );
    }
  };

  // Output mesh velocity if ALE is enabled
  if (g_inputdeck.get< tag::ale, tag::ale >()) {
    add_node_field( "x-mesh-velocity", d->meshvel().extract_comp(0) );
    add_node_field( "y-mesh-velocity", d->meshvel().extract_comp(1) );
    add_node_field( "z-mesh-velocity", d->meshvel().extract_comp(2) );
  }

  Assert( elemfieldnames.size() == elemfields.size(), "Size mismatch" );
  Assert( nodefieldnames.size() == nodefields.size(), "Size mismatch" );

  // Collect surface output names
  auto surfnames = g_dgpde[d->MeshId()].surfNames();

  // Collect surface field solution
  const auto& fd = myGhosts()->m_fd;
  auto elemsurfs = g_dgpde[d->MeshId()].surfOutput(
    fd, myGhosts()->m_geoFace, myGhosts()->m_inpoel, myGhosts()->m_coord,
    m_u, m_p );

  // Output chare mesh and fields metadata to file
  const auto& triinpoel = m_outmesh.triinpoel;
  d->write( inpoel, m_outmesh.coord, m_outmesh.bface, {},
            tk::remap( triinpoel, lid ), elemfieldnames, nodefieldnames,
            surfnames, {}, elemfields, nodefields, elemsurfs, {}, c );
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
  if (m_stage < m_nstage)
    next();
  else {
    // Ensure new field output file if ALE is enabled
    if (g_inputdeck.get< tag::ale, tag::ale >()) {
      Disc()->Itf() = 0;  // Zero field output iteration count if mesh moved
      ++Disc()->Itr();    // Increase number of iterations with a change in the mesh
    }
    startFieldOutput( CkCallback(CkIndex_DG::step(), thisProxy[thisIndex]) );
  }
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
  if (d->restarted( nrestart )) m_finished = 0;

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
  auto d = Disc();

  // Output time history
  if (d->histiter() or d->histtime() or d->histrange()) {
    std::vector< std::vector< tk::real > > hist;
    auto h = g_dgpde[d->MeshId()].histOutput( d->Hist(), myGhosts()->m_inpoel,
      myGhosts()->m_coord, m_u, m_p );
    hist.insert( end(hist), begin(h), end(h) );
    d->history( std::move(hist) );
  }

  // Free memory storing output mesh
  m_outmesh.destroy();

  // Output one-liner status report to screen
  d->status();
  // Reset Runge-Kutta stage counter
  m_stage = 0;

  // If neither max iterations nor max time reached, continue, otherwise finish
  if (not m_finished) {

    evalRestart();
 
  } else {

    auto meshid = d->MeshId();
    d->contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
                   CkCallback(CkReductionTarget(Transporter,finish), d->Tr()) );

  }
}

void
DG::computeBNorm()
// *****************************************************************************
// Start computing the boundary normals for ALE
// *****************************************************************************
{
  if (g_inputdeck.get< tag::ale, tag::ale >() &&
      !g_inputdeck.get< tag::ale, tag::symmetry >().empty())
  {
    thisProxy[ thisIndex ].wait4norm();
    bnorm();
  } else {
    meshvelstart();
  }
}

void
DG::meshvelstart()
// *****************************************************************************
// Start computing the mesh mesh velocity after boundary normals are ready
// *****************************************************************************
{
  auto d = Disc();

  // Compute fluid velocity at nodes
  if (g_inputdeck.get< tag::ale, tag::ale >()) {
    const auto smoother = g_inputdeck.get< tag::ale, tag::smoother >();
    const auto meshveltype =
      g_inputdeck.get< tag::ale, tag::mesh_velocity >();

    if (smoother == ctr::MeshVelocitySmootherType::HELMHOLTZ)
      Throw( "DG-ALE does not yet support the Helmholtz mesh velocity "
             "smoother" );

    if (smoother == ctr::MeshVelocitySmootherType::LAPLACE) {
      if (meshveltype != ctr::MeshVelocityType::FLUID)
        Throw( "DG-ALE Laplace mesh velocity smoothing is currently "
               "supported only with mesh_velocity = \"fluid\"" );

      const auto& meshforce = g_inputdeck.get< tag::ale, tag::meshforce >();
      const auto eps = std::numeric_limits< tk::real >::epsilon();
      if (!std::all_of( begin(meshforce), end(meshforce),
            [eps](tk::real c){ return std::abs(c) <= eps; } ))
        Throw( "DG-ALE Laplace mesh velocity smoothing does not yet support "
               "nonzero meshforce coefficients" );
    }

    g_dgpde[d->MeshId()].nodeVelocity( myGhosts()->m_geoElem,
      myGhosts()->m_esup, myGhosts()->m_inpoel, myGhosts()->m_coord, m_u, m_p,
      m_nodevel );
  }

  // Start computing the mesh velocity for ALE
  const auto adt = rkcoef[1][m_stage] * d->Dt();
  d->meshvelStart( m_nodevel, {}, m_bnorm, adt,
    CkCallback(CkIndex_DG::meshveldone(), thisProxy[thisIndex]) );
}

void
DG::meshveldone()
// *****************************************************************************
// Done with computing the mesh velocity for ALE
// *****************************************************************************
{
  // Assess and record mesh velocity linear solver convergence
  Disc()->meshvelConv();

  m_initial = 0;

  p_refine();
}

void
DG::imex_integrate()
// *****************************************************************************
//  Perform the Implicit-Explicit Runge-Kutta stage update
//
//!  \details
//!    Performs the Implicit-Explicit Runge-Kutta step. Scheme taken from
//!    Cavaglieri, D., & Bewley, T. (2015). Low-storage implicit/explicit
//!    Runge–Kutta schemes for the simulation of stiff high-dimensional ODE
//!    systems. Journal of Computational Physics, 286, 172-193.
//!
//!    Scheme given by equations (25a,b):
//!
//!    u[0] = u[n] + dt * (expl_rkcoef[1,0]*R_ex(u[n])+impl_rkcoef[1,1]*R_im(u[0]))
//!
//!    u[1] = u[n] + dt * (expl_rkcoef[2,1]*R_ex(u[0])+impl_rkcoef[2,1]*R_im(u[0])
//!                                                   +impl_rkcoef[2,2]*R_im(u[1]))
//!
//!    u[n+1] = u[n] + dt * (expl_rkcoef[3,1]*R_ex(u[0])+impl_rkcoef[3,1]*R_im(u[0])
//!                          expl_rkcoef[3,2]*R_ex(u[1])+impl_rkcoef[3,2]*R_im(u[1]))
//!
//!    In order to solve the first two equations we need to solve a series of
//!    systems of non-linear equations:
//!
//!    F1(u[0]) = B1 + R1(u[0]) = 0, and
//!    F2(u[1]) = B2 + R2(u[1]) = 0,
//!
//!    where
//!
//!    B1 = u[n] + dt * expl_rkcoef[1,0]*R_ex(u[n]),
//!    R1 = dt * impl_rkcoef[1,1]*R_im(u[0]) - u([0]),
//!    B2 = u[n] + dt * (expl_rkcoef[2,1]*R_ex(u[0])+impl_rkcoef[2,1]*R_im(u[0])),
//!    R2 = dt * impl_rkcoef[2,2]*R_im(u[1]) - u([1]).
//!
//!    In order to solve the non-linear system F(U) = 0, we employ:
//!    First, Broyden's method.
//!    If that fails, Newton's method (with FD approximation for jacobian).
//!
//!    Broyden's method:
//!    ----------------
//!
//!    Taken from https://en.wikipedia.org/wiki/Broyden%27s_method.
//!    The method consists in obtaining an approximation for the inverse of the
//!    Jacobian H = J^(-1) and advancing in a quasi-newton step:
//!
//!    U[k+1] = U[k] - H[k]*F(U[k]),
//!
//!    until F(U) is close enough to zero.
//!
//!    The approximation H[k] is improved at every iteration following
//!
//!    H[k] = H[k-1] + (DU[k]-H[k-1]*DF[k])/(DU[k]^T*H[k-1]*DF[k]) * DU[k]^T*H[k-1],
//!
//!    where DU[k] = U[k] - U[k-1] and DF[k] = F(U[k]) - F(U[k-1)).
//!
//!    Newton's method:
//!    ----------------
//!
//!    Taken from https://en.wikipedia.org/wiki/Newton%27s_method
//!    The method consists in inverting the jacobian
//!    Jacobian J and advancing in a Newton step:
//!
//!    U[k+1] = U[k] - J^(-1)[k]*F(U[k]),
//!
//!    until F(U) is close enough to zero.
//!
//!
//!    This function performs the following main algorithmic steps:
//!    - If stage == 0 or stage == 1:
//!      - Take Initial value:
//!        U[0] = U[n] + dt * expl_rkcoef[1,0]*R_ex(U[n]) (for stage 0)
//!        U[1] = U[n] + dt * (expl_rkcoef[2,1]*R_ex(U[0])
//!                           +impl_rkcoef[2,1]*R_im(U[0])) (for stage 1)
//!      - Loop over the Elements (e++)
//!        - Broyden steps:
//!        - Initialize Jacobian inverse approximation using FD
//!        - Compute implicit right-hand-side (F_im) with current U
//!        - Iterate for the solution (iter++)
//!          - Perform line search prior to solution update
//!          - Compute new solution U[k+1] = U[k] - H[k]*F(U[k])
//!          - Compute implicit right-hand-side (F_im) with current U
//!          - Compute DU and DF
//!          - Update inverse Jacobian approximation by:
//!            - Compute V1 = H[k-1]*DF[k] and V2 = DU[k]^T*H[k-1]
//!            - Compute d = DU[k]^T*V1 and V3 = DU[k]-V1
//!            - Compute V4 = V3/d
//!            - Update H[k] = H[k-1] + V4*V2
//!          - Save old U and F
//!          - Compute absolute and relative errors
//!          - Break iterations if error < tol or iter == max_iter
//!        - Newton steps:
//!          - Initialize Jacobian using FD approximation.
//!          - Compute implicit right-hand-side (F_im) with current U
//!          - Iterate for the solution (iter++)
//!          - Perform line search prior to solution update
//!          - Compute new solution U[k+1] = U[k] - J^(-1)[k]*F(U[k])
//!          - Compute implicit right-hand-side (F_im) with current U
//!          - Compute DU and DF
//!          - Save old U and F
//!          - Compute absolute and relative errors
//!          - Break iterations if error < tol or iter == max_iter
//!       - Update explicit equations using only the explicit terms.
//!    - Else if stage == 2:
//!       - Update explicit equations using only the explicit terms.
//!       - Update implicit equations using:
//!       u[n+1] = u[n]+dt*(expl_rkcoef[3,1]*R_ex(u[0])+impl_rkcoef[3,1]*R_im(u[0])
//!                         expl_rkcoef[3,2]*R_ex(u[1])+impl_rkcoef[3,2]*R_im(u[1]))
// *****************************************************************************
{
  auto d = Disc();
  const auto rdof = g_inputdeck.get< tag::rdof >();
  const auto ndof = g_inputdeck.get< tag::ndof >();
  if (m_stage < m_nstage-1) {
    // Save previous stiff_rhs
    m_stiffrhsprev = m_stiffrhs;

    // Compute the imex update
    const auto nelem = myGhosts()->m_fd.Esuel().size()/4;
    const auto neq = m_u.nprop()/rdof;

    for (std::size_t e=0; e<nelem; ++e) {
      auto vole = myGhosts()->m_geoElem(e,0);
      auto vole_n = m_geoElemn(e,0);
      // Integrate explicitly on all equations
      for (std::size_t c=0; c<neq; ++c)
      {
        for (std::size_t k=0; k<m_numEqDof[c]; ++k)
        {
          auto rmark = c*rdof+k;
          auto mark = c*ndof+k;

          auto mm_i = vole * mass_dubiner[k];
          auto mm_n = vole_n * mass_dubiner[k];

          m_u(e, rmark) = ( mm_n * m_un(e, rmark) + d->Dt() * (
            expl_rkcoef[0][m_stage] * m_rhsprev(e, mark)
            + expl_rkcoef[1][m_stage] * m_rhs(e, mark) ) ) / mm_i;
          if(fabs(m_u(e, rmark)) < 1e-16)
            m_u(e, rmark) = 0;
        }
      }
      // Integrate previous implicit step, which is now explicit
      for (std::size_t c=0; c<m_nstiffeq; ++c)
      {
        for (std::size_t k=0; k<m_numEqDof[c]; ++k)
        {
          auto rmark = m_stiffEqIdx[c]*rdof+k;

          auto mm_i = vole * mass_dubiner[k];

          m_u(e, rmark) += d->Dt() *
            ( impl_rkcoef[0][m_stage]
            * m_stiffrhsprev(e,c*ndof+k) / mm_i );
          if(fabs(m_u(e, rmark)) < 1e-16)
            m_u(e, rmark) = 0;
        }
      }
    }

    // Solve for implicit-explicit equations
    for (std::size_t e=0; e<nelem; ++e)
    {

      // Non-linear solver solves for x.
      // Copy the relevant variables from the state array into x.
      std::vector< tk::real > x(m_nstiffeq*ndof, 0.0);
      for (size_t ieq=0; ieq<m_nstiffeq; ++ieq)
        for (size_t idof=0; idof<m_numEqDof[ieq]; ++idof)
        {
          auto stiffrmark = m_stiffEqIdx[ieq]*rdof+idof;
          x[ieq*ndof+idof] = m_u(e, stiffrmark);
        }

      // Save all the values of m_u at stiffEqIdx as x_star,
      // They will serve to balance the energy exchange
      // from the implicit step
      auto x_star = x;

      // Solve nonlinear system, first try broyden
      bool solver_failed = false;
      x = DG::nonlinear_broyden(e, x, solver_failed);

      // If solver_failed, do newton
      if (solver_failed) {
        solver_failed = false;
        x = DG::nonlinear_newton(e, x, solver_failed);
      }

      // If newton failed, crash
      if (solver_failed)
        Throw("At element " + std::to_string(e) +
              " nonlinear solvers was not able to converge");

      // Balance energy
      g_dgpde[d->MeshId()].balance_plastic_energy(e, x_star, x, m_un);

      // Update the state u with the converged vector x.
      for (size_t ieq=0; ieq<m_nstiffeq; ++ieq)
        for (size_t idof=0; idof<m_numEqDof[ieq]; ++idof)
        {
          auto stiffrmark = m_stiffEqIdx[ieq]*rdof+idof;
          m_u(e, stiffrmark) = x[ieq*ndof+idof];
        }

    }

  }
  else {
    // For last stage just use all previously computed stages
    const auto nelem = myGhosts()->m_fd.Esuel().size()/4;
    const auto neq = m_u.nprop()/rdof;
    for (std::size_t e=0; e<nelem; ++e)
    {
      auto vole = myGhosts()->m_geoElem(e,0);
      auto vole_n = m_geoElemn(e,0);
      // First integrate explicitly on all equations
      for (std::size_t c=0; c<neq; ++c)
      {
        for (std::size_t k=0; k<m_numEqDof[c]; ++k)
        {
          auto rmark = c*rdof+k;
          auto mark = c*ndof+k;

          auto mm_i = vole * mass_dubiner[k];
          auto mm_n = vole_n * mass_dubiner[k];

          m_u(e, rmark) =  ( mm_n * m_un(e, rmark) + d->Dt() * (
            expl_rkcoef[0][m_stage] * m_rhsprev(e, mark)
            + expl_rkcoef[1][m_stage] * m_rhs(e, mark)) ) / mm_i;
          if(fabs(m_u(e, rmark)) < 1e-16)
            m_u(e, rmark) = 0;
        }
      }
      // Then, integrate the implicit part
      for (std::size_t ieq=0; ieq<m_nstiffeq; ++ieq)
        for (std::size_t idof=0; idof<m_numEqDof[ieq]; ++idof)
        {
          auto rmark = m_stiffEqIdx[ieq]*rdof+idof;

          auto mm_i = vole * mass_dubiner[idof];

          m_u(e, rmark) +=
            d->Dt() * ( impl_rkcoef[0][m_stage]
                      * m_stiffrhsprev(e,ieq*ndof+idof) / mm_i
                      + impl_rkcoef[1][m_stage]
                      * m_stiffrhs(e,ieq*ndof+idof) / mm_i );
          if(fabs(m_u(e, rmark)) < 1e-16)
            m_u(e, rmark) = 0;
        }
    }
  }
}

void
DG::BDF1_integrate()
// *****************************************************************************
//  Perform the BDF1 update
//! \details This function updates the solution using the BDF1 (backward Euler)
//!   time discretization.
// *****************************************************************************
{
  //TODO: implicit solver:
  // update solution m_u
}

std::vector< tk::real > DG::nonlinear_func(std::size_t e,
                                           std::vector< tk::real > x)
// *****************************************************************************
// Evaluate the stiff RHS and stiff equations f = b - A(x)
//! \param[in] e Element number
//! \param[in,out] x Array of unknowns to solve for
//! \details
//!   Defines the F(x) function that the non-linear solvers
//!   look to minimize. Deals with properly calling the stiff
//!   RHS functions.
// *****************************************************************************
{
  auto d = Disc();
  const auto rdof = g_inputdeck.get< tag::rdof >();
  const auto ndof = g_inputdeck.get< tag::ndof >();
  std::size_t n = x.size();

  // m_u <- x
  for (size_t ieq=0; ieq<m_nstiffeq; ++ieq)
    for (size_t idof=0; idof<m_numEqDof[ieq]; ++idof)
    {
      auto stiffrmark = m_stiffEqIdx[ieq]*rdof+idof;
      m_u(e, stiffrmark) = x[ieq*ndof+idof];
    }

  auto vole = myGhosts()->m_geoElem(e,0);
  auto vole_n = m_geoElemn(e,0);

  // Compute explicit terms (Should be computed once)
  std::vector< tk::real > expl_terms(n, 0.0);
  for (size_t ieq=0; ieq<m_nstiffeq; ++ieq)
    for (size_t idof=0; idof<m_numEqDof[ieq]; ++idof)
    {
      auto stiffmark = m_stiffEqIdx[ieq]*ndof+idof;
      auto stiffrmark = m_stiffEqIdx[ieq]*rdof+idof;
      auto mm_i = vole * mass_dubiner[idof];
      auto mm_n = vole_n * mass_dubiner[idof];
      expl_terms[ieq*ndof+idof] = ( mm_n * m_un(e, stiffrmark)
        + d->Dt() * ( expl_rkcoef[0][m_stage]
        * m_rhsprev(e,stiffmark)
        + expl_rkcoef[1][m_stage]
        * m_rhs(e,stiffmark)
        + impl_rkcoef[0][m_stage]
        * m_stiffrhsprev(e,ieq*ndof+idof)) ) / mm_i;
    }

  // Compute stiff_rhs
  g_dgpde[d->MeshId()].stiff_rhs( e, myGhosts()->m_geoElem,
    m_u, m_ndof, m_stiffrhs );

  // Store f
  std::vector< tk::real > f(n, 0.0);
  for (std::size_t ieq=0; ieq<m_nstiffeq; ++ieq)
    for (std::size_t idof=0; idof<m_numEqDof[ieq]; ++idof)
    {
      auto stiffrmark = m_stiffEqIdx[ieq]*rdof+idof;
      auto mm_i = vole * mass_dubiner[idof];
      f[ieq*ndof+idof] = expl_terms[ieq*ndof+idof]
        + d->Dt() * impl_rkcoef[1][m_stage]
        * m_stiffrhs(e,ieq*ndof+idof) / mm_i
        - m_u(e, stiffrmark);
    }

  return f;
}

std::vector< tk::real > DG::nonlinear_broyden(std::size_t e,
                                              std::vector< tk::real > x,
                                              bool solver_failed )
// *****************************************************************************
// Performs Broyden's method to solve a non-linear system on
// element e.
//! \param[in] e Element number
//! \param[in,out] x Array of unknowns to solve for
//! \param[out] solver_failed Returns true if solver did not converge
//! \details
//!    Taken from https://en.wikipedia.org/wiki/Broyden%27s_method.
//!    The method consists in obtaining an approximation for the inverse of the
//!    Jacobian H = J^(-1) and advancing in a quasi-newton step:
//!
//!    U[k+1] = U[k] - H[k]*F(U[k]),
//!
//!    until F(U) is close enough to zero.
// *****************************************************************************
{
  // Broyden's method
  // Control parameters
  std::size_t max_iter = g_inputdeck.get< tag::imex_maxiter >();
  tk::real rel_tol = g_inputdeck.get< tag::imex_reltol >();
  tk::real abs_tol = g_inputdeck.get< tag::imex_abstol >();
  tk::real rel_err = rel_tol+1;
  tk::real abs_err = abs_tol+1;
  std::size_t n = x.size();

  // Compute f with initial guess
  std::vector< tk::real > f = DG::nonlinear_func(e, x);

  // Initialize x_old and f_old
  std::vector< tk::real > x_old(n, 0.0), f_old(n, 0.0);
  for (std::size_t i=0; i<n; ++i)
  {
    x_old[i] = x[i];
    f_old[i] = f[i];
  }

  // Initialize delta_x and delta_f
  std::vector< tk::real > delta_x(n, 0.0), delta_f(n, 0.0);

  // Store the norm of f initially, for relative error measure
  tk::real err0 = 0.0;
  for (std::size_t i=0; i<n; ++i)
    err0 += f[i]*f[i];
  err0 = std::sqrt(err0);

  // Iterate for the solution if err0 > 0
  solver_failed = false;
  tk::real alpha_jacob = 1.0;
  if (err0 > abs_tol) {

    // Evaluate finite difference based jacobian
    std::vector< double > jacob(n*n);
    tk::real dx = 0.0;
    for (std::size_t i=0; i<n; ++i)
      for (std::size_t j=0; j<n; ++j)
      {
        // Set dx in the order 1% of the unknown
        dx = std::max(std::abs(0.1*x[j]),1.0e-06);
        // Derivative of f[i] with respect to x[j]
        auto x_perturb = x;
        x_perturb[j] += dx;
        auto f_perturb = DG::nonlinear_func(e, x_perturb);
        jacob[i*n+j] = (f_perturb[i]-f[i])/dx;
      }

    // Initialize Jacobian to be the inverse of this jacobian
    lapack_int ln = static_cast< lapack_int >(n);
    std::vector< lapack_int > ipiv(n);

    #ifndef NDEBUG
    lapack_int ierr =
    #endif
      LAPACKE_dgetrf(LAPACK_ROW_MAJOR, ln, ln, jacob.data(), ln, ipiv.data());
    Assert(ierr==0, "Lapack error in LU factorization of FD Jacobian");

    #ifndef NDEBUG
    lapack_int jerr =
    #endif
      LAPACKE_dgetri(LAPACK_ROW_MAJOR, ln, jacob.data(), ln, ipiv.data());
    Assert(jerr==0, "Lapack error in inverting FD Jacobian");

    std::vector< std::vector< tk::real > >
      approx_jacob(n, std::vector< tk::real >(n, 0.0));
    for (std::size_t i=0; i<n; ++i)
      for (std::size_t j=0; j<n; ++j)
        approx_jacob[i][j] = jacob[i*n+j];

    for (size_t iter=0; iter<max_iter; ++iter)
    {

      // Scale inverse of jacobian if things are not going well
      for (std::size_t i=0; i<n; ++i)
        for (std::size_t j=0; j<n; ++j)
          approx_jacob[i][j] *= alpha_jacob;

      // Compute new solution
      std::vector < tk::real > delta(n, 0.0);
      for (std::size_t i=0; i<n; ++i)
      {
        delta[i] = 0.0;
        for (std::size_t j=0; j<n; ++j)
          delta[i] -= approx_jacob[i][j] * f[j];
      }

      // Update x using line search
      bool ls_failed = false;
      tk::real alpha_ls = 1.0E+00;
      std::size_t nline = 25;
      auto xtest = x;
      for (std::size_t iline = 0; iline<nline; ++iline)
      {
        // Evaluate xtest
        for (std::size_t i=0; i<n; ++i)
          xtest[i] = x[i] + alpha_ls*delta[i];

        // Compute new f(x)
        f = DG::nonlinear_func(e, xtest);

        tk::real err = 0.0;
        for (std::size_t i=0; i<n; ++i)
          err += f[i]*f[i];
        abs_err = std::sqrt(err);

        // If 1. The error went up
        // or 2. The function f flipped in sign
        // Reduce the factor alpha_ls
        bool flipped_sign = false;
        for (std::size_t i=0; i<n; ++i)
          if (f_old[i]*f[i] < 0.0) {
            flipped_sign = true;
            break;
          }

        if (!flipped_sign)
        {
          break;
        }
        else
        {
          alpha_ls *= 0.5;
        }
        if (iline == nline-1) {
          // Try again by reducing the jacobian,
          // but only a few times, otherwise give up
          alpha_jacob *= 0.5;
          if (alpha_jacob < 1.0E-03)
            solver_failed = true;
          else
            ls_failed = true;
        }
      }

      if (solver_failed) {
        break;
      }

      if (!ls_failed) {
        // Save x
        for (std::size_t i=0; i<n; ++i)
          x[i] = xtest[i];

        // Compute delta_x and delta_f
        for (std::size_t i=0; i<n; ++i)
        {
          delta_x[i] = x[i] - x_old[i];
          delta_f[i] = f[i] - f_old[i];
        }

        // Update inverse Jacobian approximation

        // 1. Compute approx_jacob*delta_f and delta_x*jacob_approx
        tk::real sum1, sum2;
        std::vector< tk::real > auxvec1(n, 0.0), auxvec2(n, 0.0);
        for (std::size_t i=0; i<n; ++i)
        {
          sum1 = 0.0;
          sum2 = 0.0;
          for (std::size_t j=0; j<n; ++j)
          {
            sum1 += approx_jacob[i][j]*delta_f[j];
            sum2 += delta_x[j]*approx_jacob[j][i];
          }
          auxvec1[i] = sum1;
          auxvec2[i] = sum2;
        }

        // 2. Compute delta_x*approx_jacob*delta_f
        // and delta_x-approx_jacob*delta_f
        tk::real denom = 0.0;
        for (std::size_t i=0; i<n; ++i)
        {
          denom += delta_x[i]*auxvec1[i];
          auxvec1[i] = delta_x[i]-auxvec1[i];
        }

        // 3. Divide delta_u+approx_jacob*delta_f
        // by delta_x*(approx_jacob*delta_f)
        if (std::abs(denom) < 1.0e-18)
        {
          if (denom < 0.0)
          {
            for (std::size_t i=0; i<n; ++i)
              auxvec1[i] /= -1.0e-18;
          }
          else
          {
            for (std::size_t i=0; i<n; ++i)
              auxvec1[i] /= 1.0e-18;
          }
        }
        else
        {
          for (std::size_t i=0; i<n; ++i)
            auxvec1[i] /= denom;
        }

        // 4. Perform outter product between the two arrays and
        // add that quantity to the new jacobian approximation
        for (std::size_t i=0; i<n; ++i)
          for (std::size_t j=0; j<n; ++j)
            approx_jacob[i][j] += auxvec1[i] * auxvec2[j];

        // Compute a measure of error, use norm of f
        tk::real err = 0.0;
        for (std::size_t i=0; i<n; ++i)
          err += f[i]*f[i];
        abs_err = std::sqrt(err);
        rel_err = abs_err/err0;

        // Save solution and f
        for (std::size_t i=0; i<n; ++i)
        {
          x_old[i] = x[i];
          f_old[i] = f[i];
        }

        // check if error condition is met and loop back
        if (rel_err < rel_tol || abs_err < abs_tol)
          break;

        // If we did not converge, print a message and keep going
        if (iter == max_iter-1)
        {
          solver_failed = true;
        }
      }
    }
  }

  return x;
}

std::vector< tk::real > DG::nonlinear_newton(std::size_t e,
                                             std::vector< tk::real > x,
                                             bool solver_failed )
// *****************************************************************************
// Performs Newton's method to solve a non-linear system on
// element e.
//! \param[in] e Element number
//! \param[in,out] x Array of unknowns to solve for
//! \param[out] solver_failed Returns true if solver did not converge
//! \details
//!    Taken from https://en.wikipedia.org/wiki/Newton%27s_method
//!    The method consists in inverting the jacobian
//!    Jacobian J and advancing in a Newton step:
//!
//!    U[k+1] = U[k] - J^(-1)[k]*F(U[k]),
//!
//!    until F(U) is close enough to zero.
// *****************************************************************************
{
  // Newton's method
  // Control parameters
  std::size_t max_iter = g_inputdeck.get< tag::imex_maxiter >();
  tk::real rel_tol = g_inputdeck.get< tag::imex_reltol >();
  tk::real abs_tol = g_inputdeck.get< tag::imex_abstol >();
  tk::real rel_err = rel_tol+1;
  tk::real abs_err = abs_tol+1;
  std::size_t n = x.size();

  // Define jacobian
  std::vector< double > jacob(n*n);

  // Compute f with initial guess
  std::vector< tk::real > f = DG::nonlinear_func(e, x);

  // Store the norm of f initially, for relative error measure
  tk::real err0 = 0.0;
  for (std::size_t i=0; i<n; ++i)
    err0 += f[i]*f[i];
  err0 = std::sqrt(err0);
  auto abs_err_old = err0;

  // Iterate for the solution if err0 > 0
  solver_failed = false;
  tk::real alpha_jacob = 1.0;
  if (err0 > abs_tol)
    for (std::size_t iter=0; iter<max_iter; ++iter)
    {

      // Evaluate jacobian
      tk::real dx = 0.0;
      for (std::size_t i=0; i<n; ++i)
        for (std::size_t j=0; j<n; ++j)
        {
          // Set dx in the order 1% of the unknown
          dx = alpha_jacob*std::max(std::abs(0.1*x[j]),1.0e-06);
          // Derivative of f[i] with respect to x[j]
          auto x_perturb = x;
          x_perturb[j] += dx;
          auto f_perturb = DG::nonlinear_func(e, x_perturb);
          jacob[i*n+j] = (f_perturb[i]-f[i])/dx;
        }

      // Compute new solution by solving linear system J*dx = -f
      lapack_int ln = static_cast< lapack_int >(n);
      std::vector< double > delta(n);
      for (std::size_t i=0; i<n; ++i)
        delta[i] = -f[i];
      lapack_int info;
      std::vector< lapack_int > ipiv(n);
      info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, ln, 1, jacob.data(), ln, ipiv.data(), delta.data(), 1);

      if (info != 0) {
        printf("Failed with info: %ld\n", info);
      }

      // Save f as fold
      std::vector< tk::real > fold(n);
      for (std::size_t i=0; i<n; ++i)
        fold[i] = f[i];

      // Update x using line search
      bool ls_failed = false;
      tk::real alpha_ls = 1.0E+00;
      std::size_t nline = 25;
      auto xtest = x;
      for (std::size_t iline = 0; iline<nline; ++iline)
      {
        // Evaluate xtest
        for (std::size_t i=0; i<n; ++i)
          xtest[i] = x[i] + alpha_ls*delta[i];

        // Compute new f(x)
        f = DG::nonlinear_func(e, xtest);

        tk::real err = 0.0;
        for (std::size_t i=0; i<n; ++i)
          err += f[i]*f[i];
        abs_err = std::sqrt(err);

        // If 1. The error went up
        // or 2. The function f flipped in sign
        // Reduce the factor alpha_ls
        bool flipped_sign = false;
        for (std::size_t i=0; i<n; ++i)
          if (fold[i]*f[i] < 0.0) {
            flipped_sign = true;
            break;
          }

        if (abs_err < abs_err_old && !flipped_sign)
        {
          break;
        }
        else
        {
          alpha_ls *= 0.5;
        }
        if (iline == nline-1) {
          //printf("Line search failed to decrease f\n");
          // Try again by reducing the jacobian,
          // but only a few times, otherwise give up
          alpha_jacob *= 0.5;
          if (alpha_jacob < 1.0E-03)
            solver_failed = true;
          else
            ls_failed = true;
        }
      }

      if (solver_failed) {
        f = DG::nonlinear_func(e, x);
        printf("\nIMEX-RK: Non-linear solver did not converge in %lu iterations\n", iter+1);
        printf("Element #%lu\n", e);
        printf("Relative error: %e\n", rel_err);
        printf("Absolute error: %e\n\n", abs_err);
        break;
      }

      if (!ls_failed) {
        // Save x
        for (std::size_t i=0; i<n; ++i)
          x[i] = xtest[i];

        // Compute a measure of error, use norm of f
        tk::real err = 0.0;
        for (std::size_t i=0; i<n; ++i)
          err += f[i]*f[i];
        abs_err = std::sqrt(err);
        rel_err = abs_err/err0;

        // check if error condition is met and loop back
        if (rel_err < rel_tol || abs_err < abs_tol) {
          break;
        }

        // If we did not converge, print a message and keep going
        if (iter == max_iter-1)
        {
          printf("\nIMEX-RK: Non-linear solver did not converge in %lu iterations\n", max_iter);
          printf("Element #%lu\n", e);
          printf("Relative error: %e\n", rel_err);
          printf("Absolute error: %e\n\n", abs_err);
        }
      }
    }

  return x;

}

//------------------------------------------------------------------------------
// Unused Nodal Extrema code
//------------------------------------------------------------------------------
// The following code computes the 'nodal extrema' of solutions that can be used
// for high-order limiting purposes. However, these are currently unused, due to
// more effective limiting methods in use. Hence the code is commented. The code
// itself is quite complex and it is worthwhile to keep it.
//------------------------------------------------------------------------------
// 1) in DG.hpp:
//    void pup( PUP::er &p ) override {
//      p | m_ndof_NodalExtrm;
//      p | m_nnodalExtrema;
//      p | m_uNodalExtrm;
//      p | m_pNodalExtrm;
//      p | m_uNodalExtrmc;
//      p | m_pNodalExtrmc;
//    }
//    //! \brief Degree of freedom for nodal extrema vector. When DGP1 is applied,
//    //!   there is one degree of freedom for cell average variable. When DGP2 is
//    //!   applied, the degree of freedom is 4 which refers to cell average and
//    //!   gradients in three directions
//    std::size_t m_ndof_NodalExtrm;
//    //! \brief Counter signaling that we have received all our nodal extrema from
//    //!   ghost chare partitions
//    std::size_t m_nnodalExtrema;
//    //! Vector of nodal extrema for conservative variables
//    std::vector< std::vector<tk::real> > m_uNodalExtrm;
//    //! Vector of nodal extrema for primitive variables
//    std::vector< std::vector<tk::real> > m_pNodalExtrm;
//    //! Buffer for vector of nodal extrema for conservative variables
//    std::unordered_map< std::size_t, std::vector< tk::real > > m_uNodalExtrmc;
//    //! Buffer for vector of nodal extrema for primitive variables
//    std::unordered_map< std::size_t, std::vector< tk::real > > m_pNodalExtrmc;
//
// 2) in dg.ci:
//      entry void comnodalExtrema( const std::vector< std::size_t >& gid,
//                                  const std::vector< std::vector< tk::real > >& G1,
//                                  const std::vector< std::vector< tk::real > >& G2 );
//
//      >> Call nodalExtrema() when wait4reco() is done.
//
//      entry void wait4nodalExtrema() {
//        when ownnodalExtrema_complete(), comnodalExtrema_complete()
//        serial "nodalExtrema" { lim(); } }
//
//      entry void ownnodalExtrema_complete();
//      entry void comnodalExtrema_complete();
//
// 3) in DG.cpp:
//  m_ndof_NodalExtrm( 3 ), // for the first order derivatives in 3 directions
//
//  // Allocate storage for the vector of nodal extrema in DG::ctor
//  m_uNodalExtrm.resize( Disc()->Bid().size(),
//    std::vector<tk::real>( 2 * m_ndof_NodalExtrm *
//    g_inputdeck.get< tag::ncomp >() ) );
//  m_pNodalExtrm.resize( Disc()->Bid().size(),
//    std::vector<tk::real>( 2 * m_ndof_NodalExtrm *
//    m_p.nprop() / g_inputdeck.get< tag::rdof >() ) );
//
//  // Initialization for the buffer vector of nodal extrema in DG::ctor
//  resizeNodalExtremac();
//
//  // In appropriate locations
//  thisProxy[ thisIndex ].wait4nodalExtrema();
//
//void
//DG::nodalExtrema()
//// *****************************************************************************
//// Compute nodal extrema at chare-boundary nodes. Extrema at internal nodes
//// are calculated in limiter function.
//// *****************************************************************************
//{
//  // Initialize nodal extrema vector
//  auto large = std::numeric_limits< tk::real >::max();
//  for(std::size_t i = 0; i<bid.size(); i++)
//  {
//    for (std::size_t c=0; c<ncomp; ++c)
//    {
//      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
//      {
//        auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
//        auto min_mark = max_mark + 1;
//        m_uNodalExtrm[i][max_mark] = -large;
//        m_uNodalExtrm[i][min_mark] =  large;
//      }
//    }
//    for (std::size_t c=0; c<nprim; ++c)
//    {
//      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
//      {
//        auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
//        auto min_mark = max_mark + 1;
//        m_pNodalExtrm[i][max_mark] = -large;
//        m_pNodalExtrm[i][min_mark] =  large;
//      }
//    }
//  }
//
//  // Evaluate the max/min value for the chare-boundary nodes
//  if(rdof > 4) {
//      evalNodalExtrmRefEl(ncomp, nprim, m_ndof_NodalExtrm, d->bndel(),
//        myGhosts()->m_inpoel, gid, bid, m_u, m_p, m_uNodalExtrm, m_pNodalExtrm);
//  }
//
//  // Communicate extrema at nodes to other chares on chare-boundary
//  if (d->NodeCommMap().empty())        // in serial we are done
//    comnodalExtrema_complete();
//  else  // send nodal extrema to chare-boundary nodes to fellow chares
//  {
//    for (const auto& [c,n] : d->NodeCommMap()) {
//      std::vector< std::vector< tk::real > > g1( n.size() ), g2( n.size() );
//      std::size_t j = 0;
//      for (auto i : n)
//      {
//        auto p = tk::cref_find(d->Bid(),i);
//        g1[ j   ] = m_uNodalExtrm[ p ];
//        g2[ j++ ] = m_pNodalExtrm[ p ];
//      }
//      thisProxy[c].comnodalExtrema( std::vector<std::size_t>(begin(n),end(n)),
//        g1, g2 );
//    }
//  }
//  ownnodalExtrema_complete();
//}
//
//void
//DG::comnodalExtrema( const std::vector< std::size_t >& gid,
//                     const std::vector< std::vector< tk::real > >& G1,
//                     const std::vector< std::vector< tk::real > >& G2 )
//// *****************************************************************************
////  Receive contributions to nodal extrema on chare-boundaries
////! \param[in] gid Global mesh node IDs at which we receive grad contributions
////! \param[in] G1 Partial contributions of extrema for conservative variables to
////!   chare-boundary nodes
////! \param[in] G2 Partial contributions of extrema for primitive variables to
////!   chare-boundary nodes
////! \details This function receives contributions to m_uNodalExtrm/m_pNodalExtrm
////!   , which stores nodal extrems at mesh chare-boundary nodes. While
////!   m_uNodalExtrm/m_pNodalExtrm stores own contributions, m_uNodalExtrmc
////!   /m_pNodalExtrmc collects the neighbor chare contributions during
////!   communication.
//// *****************************************************************************
//{
//  Assert( G1.size() == gid.size(), "Size mismatch" );
//  Assert( G2.size() == gid.size(), "Size mismatch" );
//
//  const auto rdof = g_inputdeck.get< tag::rdof >();
//  const auto ncomp = m_u.nprop() / rdof;
//  const auto nprim = m_p.nprop() / rdof;
//
//  for (std::size_t i=0; i<gid.size(); ++i)
//  {
//    auto& u = m_uNodalExtrmc[gid[i]];
//    auto& p = m_pNodalExtrmc[gid[i]];
//    for (std::size_t c=0; c<ncomp; ++c)
//    {
//      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
//      {
//        auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
//        auto min_mark = max_mark + 1;
//        u[max_mark] = std::max( G1[i][max_mark], u[max_mark] );
//        u[min_mark] = std::min( G1[i][min_mark], u[min_mark] );
//      }
//    }
//    for (std::size_t c=0; c<nprim; ++c)
//    {
//      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
//      {
//        auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
//        auto min_mark = max_mark + 1;
//        p[max_mark] = std::max( G2[i][max_mark], p[max_mark] );
//        p[min_mark] = std::min( G2[i][min_mark], p[min_mark] );
//      }
//    }
//  }
//
//  if (++m_nnodalExtrema == Disc()->NodeCommMap().size())
//  {
//    m_nnodalExtrema = 0;
//    comnodalExtrema_complete();
//  }
//}
//
//void DG::resizeNodalExtremac()
//// *****************************************************************************
////  Resize the buffer vector of nodal extrema
//// *****************************************************************************
//{
//  const auto rdof = g_inputdeck.get< tag::rdof >();
//  const auto ncomp = m_u.nprop() / rdof;
//  const auto nprim = m_p.nprop() / rdof;
//
//  auto large = std::numeric_limits< tk::real >::max();
//  for (const auto& [c,n] : Disc()->NodeCommMap())
//  {
//    for (auto i : n) {
//      auto& u = m_uNodalExtrmc[i];
//      auto& p = m_pNodalExtrmc[i];
//      u.resize( 2*m_ndof_NodalExtrm*ncomp, large );
//      p.resize( 2*m_ndof_NodalExtrm*nprim, large );
//
//      // Initialize the minimum nodal extrema
//      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
//      {
//        for(std::size_t k = 0; k < ncomp; k++)
//          u[2*k*m_ndof_NodalExtrm+2*idof] = -large;
//        for(std::size_t k = 0; k < nprim; k++)
//          p[2*k*m_ndof_NodalExtrm+2*idof] = -large;
//      }
//    }
//  }
//}
//
//void DG::evalNodalExtrmRefEl(
//  const std::size_t ncomp,
//  const std::size_t nprim,
//  const std::size_t ndof_NodalExtrm,
//  const std::vector< std::size_t >& bndel,
//  const std::vector< std::size_t >& inpoel,
//  const std::vector< std::size_t >& gid,
//  const std::unordered_map< std::size_t, std::size_t >& bid,
//  const tk::Fields& U,
//  const tk::Fields& P,
//  std::vector< std::vector<tk::real> >& uNodalExtrm,
//  std::vector< std::vector<tk::real> >& pNodalExtrm )
//// *****************************************************************************
////  Compute the nodal extrema of ref el derivatives for chare-boundary nodes
////! \param[in] ncomp Number of conservative variables
////! \param[in] nprim Number of primitive variables
////! \param[in] ndof_NodalExtrm Degree of freedom for nodal extrema
////! \param[in] bndel List of elements contributing to chare-boundary nodes
////! \param[in] inpoel Element-node connectivity for element e
////! \param[in] gid Local->global node id map
////! \param[in] bid Local chare-boundary node ids (value) associated to
////!   global node ids (key)
////! \param[in] U Vector of conservative variables
////! \param[in] P Vector of primitive variables
////! \param[in,out] uNodalExtrm Chare-boundary nodal extrema for conservative
////!   variables
////! \param[in,out] pNodalExtrm Chare-boundary nodal extrema for primitive
////!   variables
//// *****************************************************************************
//{
//  const auto rdof = g_inputdeck.get< tag::rdof >();
//
//  for (auto e : bndel)
//  {
//    // access node IDs
//    const std::vector<std::size_t> N
//      { inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };
//
//    // Loop over nodes of element e
//    for(std::size_t ip=0; ip<4; ++ip)
//    {
//      auto i = bid.find( gid[N[ip]] );
//      if (i != end(bid))      // If ip is the chare boundary point
//      {
//        // If DG(P2) is applied, find the nodal extrema of the gradients of
//        // conservative/primitive variables in the reference element
//
//        // Vector used to store the first order derivatives for both
//        // conservative and primitive variables
//        std::vector< std::array< tk::real, 3 > > gradc(ncomp, {0.0, 0.0, 0.0});
//        std::vector< std::array< tk::real, 3 > > gradp(ncomp, {0.0, 0.0, 0.0});
//
//        // Derivatives of the Dubiner basis
//        std::array< tk::real, 3 > center {{0.25, 0.25, 0.25}};
//        auto dBdxi = tk::eval_dBdxi(rdof, center);
//
//        // Evaluate the first order derivative
//        for(std::size_t icomp = 0; icomp < ncomp; icomp++)
//        {
//          auto mark = icomp * rdof;
//          for(std::size_t idir = 0; idir < 3; idir++)
//          {
//            gradc[icomp][idir] = 0;
//            for(std::size_t idof = 1; idof < rdof; idof++)
//              gradc[icomp][idir] += U(e, mark+idof) * dBdxi[idir][idof];
//          }
//        }
//        for(std::size_t icomp = 0; icomp < nprim; icomp++)
//        {
//          auto mark = icomp * rdof;
//          for(std::size_t idir = 0; idir < 3; idir++)
//          {
//            gradp[icomp][idir] = 0;
//            for(std::size_t idof = 1; idof < rdof; idof++)
//              gradp[icomp][idir] += P(e, mark+idof) * dBdxi[idir][idof];
//          }
//        }
//
//        // Store the extrema for the gradients
//        for (std::size_t c=0; c<ncomp; ++c)
//        {
//          for (std::size_t idof = 0; idof < ndof_NodalExtrm; idof++)
//          {
//            auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
//            auto min_mark = max_mark + 1;
//            auto& ex = uNodalExtrm[i->second];
//            ex[max_mark] = std::max(ex[max_mark], gradc[c][idof]);
//            ex[min_mark] = std::min(ex[min_mark], gradc[c][idof]);
//          }
//        }
//        for (std::size_t c=0; c<nprim; ++c)
//        {
//          for (std::size_t idof = 0; idof < ndof_NodalExtrm; idof++)
//          {
//            auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
//            auto min_mark = max_mark + 1;
//            auto& ex = pNodalExtrm[i->second];
//            ex[max_mark] = std::max(ex[max_mark], gradp[c][idof]);
//            ex[min_mark] = std::min(ex[min_mark], gradp[c][idof]);
//          }
//        }
//      }
//    }
//  }
//}
//
// >> Once communication into buffers is complete, and the next function is
// called store communicated extrema values there as follows:
//
//  // Combine own and communicated contributions to nodal extrema
//  for (const auto& [gid,g] : m_uNodalExtrmc) {
//    auto bid = tk::cref_find( d->Bid(), gid );
//    for (ncomp_t c=0; c<ncomp; ++c)
//    {
//      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
//      {
//        auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
//        auto min_mark = max_mark + 1;
//        m_uNodalExtrm[bid][max_mark] =
//          std::max(g[max_mark], m_uNodalExtrm[bid][max_mark]);
//        m_uNodalExtrm[bid][min_mark] =
//          std::min(g[min_mark], m_uNodalExtrm[bid][min_mark]);
//      }
//    }
//  }
//  for (const auto& [gid,g] : m_pNodalExtrmc) {
//    auto bid = tk::cref_find( d->Bid(), gid );
//    for (ncomp_t c=0; c<nprim; ++c)
//    {
//      for(std::size_t idof=0; idof<m_ndof_NodalExtrm; idof++)
//      {
//        auto max_mark = 2*c*m_ndof_NodalExtrm + 2*idof;
//        auto min_mark = max_mark + 1;
//        m_pNodalExtrm[bid][max_mark] =
//          std::max(g[max_mark], m_pNodalExtrm[bid][max_mark]);
//        m_pNodalExtrm[bid][min_mark] =
//          std::min(g[min_mark], m_pNodalExtrm[bid][min_mark]);
//      }
//    }
//  }
//
//  // clear gradients receive buffer
//  tk::destroy(m_uNodalExtrmc);
//  tk::destroy(m_pNodalExtrmc);
//
//  >> Resize these buffers in dt() prior to the contribute-call by
//  calling resizeNodalExtremac();
//
//------------------------------------------------------------------------------

#include "NoWarning/dg.def.h"
