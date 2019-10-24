// *****************************************************************************
/*!
  \file      src/Inciter/DiagCG.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     DiagCG for a PDE system with continuous Galerkin without a matrix
  \details   DiagCG advances a system of partial differential equations (PDEs)
    using continuous Galerkin (CG) finite element (FE) spatial discretization
    (using linear shapefunctions on tetrahedron elements) combined with a time
    stepping scheme that is equivalent to the Lax-Wendroff (LW) scheme within
    the unstructured-mesh FE context and treats discontinuities with
    flux-corrected transport (FCT). Only the diagonal entries of the left-hand
    side matrix are non-zero thus it does not need a matrix-based linear solver.
  \see The documentation in DiagCG.h.
*/
// *****************************************************************************

#include "QuinoaConfig.hpp"
#include "DiagCG.hpp"
#include "Vector.hpp"
#include "Reader.hpp"
#include "ContainerUtil.hpp"
#include "UnsMesh.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "DerivedData.hpp"
#include "CGPDE.hpp"
#include "Discretization.hpp"
#include "DistFCT.hpp"
#include "DiagReducer.hpp"
#include "NodeBC.hpp"
#include "Refiner.hpp"
#include "Reorder.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< CGPDE > g_cgpde;

} // inciter::

using inciter::DiagCG;

DiagCG::DiagCG( const CProxy_Discretization& disc,
                const std::map< int, std::vector< std::size_t > >& bface,
                const std::map< int, std::vector< std::size_t > >& bnode,
                const std::vector< std::size_t >& triinpoel ) :
  m_disc( disc ),
  m_initial( 1 ),
  m_nsol( 0 ),
  m_nlhs( 0 ),
  m_nrhs( 0 ),
  m_nnorm( 0 ),
  m_bnode( bnode ),
  m_u( m_disc[thisIndex].ckLocal()->Gid().size(),
       g_inputdeck.get< tag::component >().nprop() ),
  m_ul( m_u.nunk(), m_u.nprop() ),
  m_du( m_u.nunk(), m_u.nprop() ),
  m_ue( m_disc[thisIndex].ckLocal()->Inpoel().size()/4, m_u.nprop() ),
  m_lhs( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_bcdir(),
  m_lhsc(),
  m_rhsc(),
  m_difc(),
  m_vol(),
  m_bnorm(),
  m_bnormc(),
  m_diag()
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] bnode Boundary-node lists mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
// *****************************************************************************
{
  usesAtSync = true;    // enable migration at AtSync

  // Activate SDAG wait
  thisProxy[ thisIndex ].wait4norm();
  thisProxy[ thisIndex ].wait4lhs();

  // Query nodes at which symmetry BCs are specified
  std::unordered_set< std::size_t > symbcnodes;
  for (const auto& eq : g_cgpde) eq.symbcnodes( bface, triinpoel, symbcnodes );

  // Compute boundary point normals
  bnorm( bface, triinpoel, std::move(symbcnodes) );
}

void
DiagCG::bnorm( const std::map< int, std::vector< std::size_t > >& bface,
               const std::vector< std::size_t >& triinpoel,
               std::unordered_set< std::size_t >&& symbcnodes )
// *****************************************************************************
//  Compute boundary point normals
// *****************************************************************************
{
  auto d = Disc();

  const auto& coord = d->Coord();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // Lambda to compute the inverse distance squared between boundary face
  // centroid and boundary point. Here p is the global node id and g is the
  // geometry of the boundary face, see tk::geoFaceTri().
  auto invdistsq = [&]( const tk::Fields& g, std::size_t p ){
    return 1.0 / ( (g(0,4,0) - x[p])*(g(0,4,0) - x[p]) +
                   (g(0,5,0) - y[p])*(g(0,5,0) - y[p]) +
                   (g(0,6,0) - z[p])*(g(0,6,0) - z[p]) );
  };

  // Compute boundary point normals on all side sets summing inverse distance
  // weighted face normals to points. This is only a partial sum at shared
  // boundary points in parallel.
  const auto& lid = d->Lid();
  const auto& gid = d->Gid();
  for (const auto& [ setid, faceids ] : bface) {
    for (auto f : faceids) {
      tk::UnsMesh::Face
        face{ tk::cref_find( lid, triinpoel[f*3+0] ),
              tk::cref_find( lid, triinpoel[f*3+1] ),
              tk::cref_find( lid, triinpoel[f*3+2] ) };
      std::array< tk::real, 3 > fx{ x[face[0]], x[face[1]], x[face[2]] };
      std::array< tk::real, 3 > fy{ y[face[0]], y[face[1]], y[face[2]] };
      std::array< tk::real, 3 > fz{ z[face[0]], z[face[1]], z[face[2]] };
      auto g = tk::geoFaceTri( fx, fy, fz );
      for (auto p : face) {
        auto i = symbcnodes.find( gid[p] );
        if (i != end(symbcnodes)) {     // only if user set symbc on node
          tk::real r = invdistsq( g, p );
          auto& n = m_bnorm[ gid[p] ];  // associate global node id
          n[0] += r*g(0,1,0);
          n[1] += r*g(0,2,0);
          n[2] += r*g(0,3,0);
          n[3] += r;
        }
      }
    }
  }

  // Send our nodal normal contributions to neighbor chares
  if (d->Msum().empty())
   comnorm_complete();
  else
    for (const auto& [ neighborchare, sharednodes ] : d->Msum()) {
      std::unordered_map< std::size_t, std::array< tk::real, 4 > > exp;
      for (auto i : sharednodes) {      // symmetry BCs may be specified on only
        auto j = m_bnorm.find(i);       // a subset of chare boundary nodes
        if (j != end(m_bnorm)) exp[i] = j->second;
      }
      thisProxy[ neighborchare ].comnorm( exp );
    }

  ownnorm_complete();
}

void
DiagCG::comnorm(
  const std::unordered_map< std::size_t, std::array< tk::real, 4 > >& innorm )
// *****************************************************************************
// Receive boundary point normals on chare-boundaries
//! \param[in] innorm Incoming partial sums of boundary point normal
//!   contributions to normals (first 3 components), inverse distance squared
//!   (4th component)
// *****************************************************************************
{
  // Buffer up inccoming contributions
  for (const auto& [ p, n ] : innorm) {
    auto& bnorm = m_bnormc[ p ];
    bnorm[0] += n[0];
    bnorm[1] += n[1];
    bnorm[2] += n[2];
    bnorm[3] += n[3];
  }

  if (++m_nnorm == Disc()->Msum().size()) {
    m_nnorm = 0;
    comnorm_complete();
  }
}

void
DiagCG::normfinal()
// *****************************************************************************
//  Finish computing boundary point normals
// *****************************************************************************
{
  // Combine own and communicated contributions to boundary point normals
  for (auto& [ p, n ] : m_bnormc) {
    auto j = m_bnorm.find(p);
    if (j != end(m_bnorm)) {
      auto& norm = j->second;
      norm[0] += n[0];
      norm[1] += n[1];
      norm[2] += n[2];
      norm[3] += n[3];
    }
  }

  // Divie summed point normals by the sum of inverse distance squared
  for (auto& [ p, n ] : m_bnorm) {
    n[0] /= n[3];
    n[1] /= n[3];
    n[2] /= n[3];
    Assert( (n[0]*n[0] + n[1]*n[1] + n[2]*n[2] - 1.0) <
            std::numeric_limits< tk::real >::epsilon(), "Non-unit normal" );
  }

  // Replace global->local ids associated to boundary normals of symbc nodes
  decltype(m_bnorm) bnorm;
  for (auto&& [g,n] : m_bnorm) {
    bnorm[ tk::cref_find( Disc()->Lid(), g ) ] = std::move(n);
  }
  m_bnorm = std::move(bnorm);

  // Signal the runtime system that the workers have been created
  contribute( sizeof(int), &m_initial, CkReduction::sum_int,
    CkCallback(CkReductionTarget(Transporter,comfinal), Disc()->Tr()) );
}

void
DiagCG::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types initiated from this chare array
//! \details Since this is a [initnode] routine, the runtime system executes the
//!   routine exactly once on every logical node early on in the Charm++ init
//!   sequence. Must be static as it is called without an object. See also:
//!   Section "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  NodeDiagnostics::registerReducers();
}

void
DiagCG::ResumeFromSync()
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
DiagCG::setup()
// *****************************************************************************
// Set and output initial conditions and mesh to file
// *****************************************************************************
{
  auto d = Disc();

  // Set initial conditions for all PDEs
  for (const auto& eq : g_cgpde) eq.initialize( d->Coord(), m_u, d->T() );

  // Apply symmetry boundary conditions on initial conditions
  for (const auto& eq : g_cgpde) eq.symbc( m_u, m_bnorm );

  // Output initial conditions to file (regardless of whether it was requested)
  writeFields( CkCallback(CkIndex_DiagCG::init(), thisProxy[thisIndex]) );
}

void
DiagCG::init()
// *****************************************************************************
// Initially compute left hand side diagonal matrix
// *****************************************************************************
{
  lhs();
}

void
DiagCG::start()
// *****************************************************************************
//  Start time stepping
// *****************************************************************************
{
  // Start timer measuring time stepping wall clock time
  Disc()->Timer().zero();

  // Start time stepping by computing the size of the next time step)
  next();
}

void
DiagCG::next()
// *****************************************************************************
// Continue to next time step
// *****************************************************************************
{
  dt();
}

void
DiagCG::lhs()
// *****************************************************************************
// Compute the left-hand side of transport equations
// *****************************************************************************
{
  auto d = Disc();

  // Compute lumped mass lhs required for both high and low order solutions
  m_lhs = d->FCT()->lump( *d );

  if (d->Msum().empty())
    comlhs_complete();
  else // send contributions of lhs to chare-boundary nodes to fellow chares
    for (const auto& n : d->Msum()) {
      std::vector< std::vector< tk::real > > l( n.second.size() );
      std::size_t j = 0;
      for (auto i : n.second) l[ j++ ] = m_lhs[ tk::cref_find(d->Lid(),i) ];
      thisProxy[ n.first ].comlhs( n.second, l );
    }

  ownlhs_complete();
}

void
DiagCG::comlhs( const std::vector< std::size_t >& gid,
                const std::vector< std::vector< tk::real > >& L )
// *****************************************************************************
//  Receive contributions to left-hand side diagonal matrix on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive LHS contributions
//! \param[in] L Partial contributions of LHS to chare-boundary nodes
//! \details This function receives contributions to m_lhs, which stores the
//!   diagonal (lumped) mass matrix at mesh nodes. While m_lhs stores
//!   own contributions, m_lhsc collects the neighbor chare contributions during
//!   communication. This way work on m_lhs and m_lhsc is overlapped. The two
//!   are combined in lhsmerge().
// *****************************************************************************
{
  Assert( L.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  for (std::size_t i=0; i<gid.size(); ++i)
    m_lhsc[ gid[i] ] += L[i];

  if (++m_nlhs == Disc()->Msum().size()) {
    m_nlhs = 0;
    comlhs_complete();
  }
}

void
DiagCG::lhsmerge()
// *****************************************************************************
// The own and communication portion of the left-hand side is complete
// *****************************************************************************
{
  // Combine own and communicated contributions to left hand side
  for (const auto& b : m_lhsc) {
    auto lid = tk::cref_find( Disc()->Lid(), b.first );
    for (ncomp_t c=0; c<m_lhs.nprop(); ++c)
      m_lhs(lid,c,0) += b.second[c];
  }

  // Clear receive buffer
  tk::destroy(m_lhsc);

  // Continue after lhs is complete
  if (m_initial) start(); else lhs_complete();
}

void
DiagCG::dt()
// *****************************************************************************
// Compute time step size
// *****************************************************************************
{
  tk::real mindt = std::numeric_limits< tk::real >::max();

  auto const_dt = g_inputdeck.get< tag::discr, tag::dt >();
  auto def_const_dt = g_inputdeck_defaults.get< tag::discr, tag::dt >();
  auto eps = std::numeric_limits< tk::real >::epsilon();

  auto d = Disc();

  // use constant dt if configured
  if (std::abs(const_dt - def_const_dt) > eps) {

    mindt = const_dt;

  } else {      // compute dt based on CFL

    // find the minimum dt across all PDEs integrated
    for (const auto& eq : g_cgpde) {
      auto eqdt = eq.dt( d->Coord(), d->Inpoel(), m_u );
      if (eqdt < mindt) mindt = eqdt;
    }

    // Scale smallest dt with CFL coefficient
    mindt *= g_inputdeck.get< tag::discr, tag::cfl >();

  }

  // Actiavate SDAG waits for time step
  thisProxy[ thisIndex ].wait4rhs();
  thisProxy[ thisIndex ].wait4out();

  // Activate SDAG-waits for FCT
  d->FCT()->next();

  // Contribute to minimum dt across all chares the advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(DiagCG,advance), thisProxy) );
}

void
DiagCG::advance( tk::real newdt )
// *****************************************************************************
// Advance equations to next time step
//! \param[in] newdt Size of this new time step
// *****************************************************************************
{
  auto d = Disc();

  // Set new time step size
  d->setdt( newdt );

  // Compute rhs for next time step
  rhs();
}

void
DiagCG::rhs()
// *****************************************************************************
// Compute right-hand side of transport equations
// *****************************************************************************
{
  auto d = Disc();

  // Compute right-hand side and query Dirichlet BCs for all equations
  for (const auto& eq : g_cgpde)
    eq.rhs( d->T(), d->Dt(), d->Coord(), d->Inpoel(), m_u, m_ue, m_rhs );

  // Compute mass diffusion rhs contribution required for the low order solution
  auto dif = d->FCT()->diff( *d, m_u );

  // Query and match user-specified boundary conditions to side sets
  m_bcdir = match( m_u.nprop(), d->T(), d->Dt(), d->Coord(),
                   d->Lid(), m_bnode );

  if (d->Msum().empty())
    comrhs_complete();
  else // send contributions of rhs to chare-boundary nodes to fellow chares
    for (const auto& n : d->Msum()) {
      std::vector< std::vector< tk::real > > r( n.second.size() );
      std::vector< std::vector< tk::real > > D( n.second.size() );
      std::size_t j = 0;
      for (auto i : n.second) {
        auto lid = tk::cref_find( d->Lid(), i );
        r[ j ] = m_rhs[ lid ];
        D[ j ] = dif[ lid ];
        ++j;
      }
      thisProxy[ n.first ].comrhs( n.second, r, D );
    }

  ownrhs_complete( dif );
}

void
DiagCG::comrhs( const std::vector< std::size_t >& gid,
                const std::vector< std::vector< tk::real > >& R,
                const std::vector< std::vector< tk::real > >& D )
// *****************************************************************************
//  Receive contributions to right-hand side vector on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive RHS contributions
//! \param[in] R Partial contributions of RHS to chare-boundary nodes
//! \param[in] D Partial contributions to chare-boundary nodes
//! \details This function receives contributions to m_rhs, which stores the
//!   right hand side vector at mesh nodes. While m_rhs stores own
//!   contributions, m_rhsc collects the neighbor chare contributions during
//!   communication. This way work on m_rhs and m_rhsc is overlapped. The two
//!   are combined in solve(). This function also receives contributions to
//!   mass diffusion term of the right hand side vector at mesh nodes.
// *****************************************************************************
{
  Assert( R.size() == gid.size(), "Size mismatch" );
  Assert( D.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  for (std::size_t i=0; i<gid.size(); ++i) {
    m_rhsc[ gid[i] ] += R[i];
    m_difc[ gid[i] ] += D[i];
  }

  if (++m_nrhs == Disc()->Msum().size()) {
    m_nrhs = 0;
    comrhs_complete();
  }
}

void
DiagCG::solve( tk::Fields& dif )
// *****************************************************************************
//  Solve low and high order diagonal systems
//! \param[in,out] dif Mass diffusion own contribution
// *****************************************************************************
{
  const auto ncomp = m_rhs.nprop();

  auto d = Disc();

  // Combine own and communicated contributions to rhs
  for (const auto& b : m_rhsc) {
    auto lid = tk::cref_find( d->Lid(), b.first );
    for (ncomp_t c=0; c<ncomp; ++c) m_rhs(lid,c,0) += b.second[c];
  }

  // Combine own and communicated contributions to mass diffusion
  for (const auto& b : m_difc) {
    auto lid = tk::cref_find( d->Lid(), b.first );
    for (ncomp_t c=0; c<ncomp; ++c) dif(lid,c,0) += b.second[c];
  }

  // Clear receive buffers
  tk::destroy(m_rhsc);
  tk::destroy(m_difc);

  // Set Dirichlet BCs for lhs and both low and high order rhs vectors. Note
  // that the low order rhs (more precisely the mass-diffusion term) is set to
  // zero instead of the solution increment at Dirichlet BCs, because for the
  // low order solution the right hand side is the sum of the high order right
  // hand side and mass diffusion so the low order system is L = R + D, where L
  // is the lumped mass matrix, R is the high order RHS, and D is
  // mass diffusion, and R already will have the Dirichlet BC set.
  for (const auto& [b,bc] : m_bcdir) {
    for (ncomp_t c=0; c<ncomp; ++c) {
      if (bc[c].first) {
        m_lhs( b, c, 0 ) = 1.0;
        m_rhs( b, c, 0 ) = bc[c].second;
        dif( b, c, 0 ) = 0.0;
      }
    }
  }

  // Solve low and high order diagonal systems and update low order solution
  auto dul = (m_rhs + dif) / m_lhs;

  // Apply symmetry BCs on low order solution increment

  m_ul = m_u + dul;
  m_du = m_rhs / m_lhs;

  // Apply symmetry BCs on high order solution increment
  for (const auto& eq : g_cgpde) {
    eq.symbc( dul, m_bnorm );
    eq.symbc( m_ul, m_bnorm );
    eq.symbc( m_du, m_bnorm );
  }

  // Continue with FCT
  d->FCT()->aec( *d, m_du, m_u, m_bcdir, m_bnorm );
  d->FCT()->alw( m_u, m_ul, std::move(dul), thisProxy );
}

void
DiagCG::writeFields( CkCallback c ) const
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  auto d = Disc();

  // Query and collect field names from PDEs integrated
  std::vector< std::string > nodefieldnames;
  for (const auto& eq : g_cgpde) {
    auto n = eq.fieldNames();
    nodefieldnames.insert( end(nodefieldnames), begin(n), end(n) );
  }

  // Collect node field solution
  auto u = m_u;
  std::vector< std::vector< tk::real > > nodefields;
  for (const auto& eq : g_cgpde) {
    auto o = eq.fieldOutput( d->T(), d->meshvol(), d->Coord(), d->V(), u );
    nodefields.insert( end(nodefields), begin(o), end(o) );
  }

  std::vector< std::string > elemfieldnames;
  std::vector< std::vector< tk::real > > elemfields;

  // Query refinement data
  auto dtref = g_inputdeck.get< tag::amr, tag::dtref >();

  std::tuple< std::vector< std::string >,
              std::vector< std::vector< tk::real > >,
              std::vector< std::string >,
              std::vector< std::vector< tk::real > > > r;
  if (dtref) r = d->Ref()->refinementFields();

  const auto& refinement_elemfieldnames = std::get< 0 >( r );
  const auto& refinement_elemfields = std::get< 1 >( r );
  const auto& refinement_nodefieldnames = std::get< 2 >( r );
  const auto& refinement_nodefields = std::get< 3 >( r );

  nodefieldnames.insert( end(nodefieldnames),
    begin(refinement_nodefieldnames), end(refinement_nodefieldnames) );
  nodefields.insert( end(nodefields),
    begin(refinement_nodefields), end(refinement_nodefields) );

  elemfieldnames.insert( end(elemfieldnames),
    begin(refinement_elemfieldnames), end(refinement_elemfieldnames) );
  elemfields.insert( end(elemfields),
    begin(refinement_elemfields), end(refinement_elemfields) );

  // Collect FCT field data (for debugging)
  auto f = d->FCT()->fields();

  const auto& fct_elemfieldnames = std::get< 0 >( f );
  const auto& fct_elemfields = std::get< 1 >( f );
  const auto& fct_nodefieldnames = std::get< 2 >( f );
  const auto& fct_nodefields = std::get< 3 >( f );

  nodefieldnames.insert( end(nodefieldnames),
    begin(fct_nodefieldnames), end(fct_nodefieldnames) );
  nodefields.insert( end(nodefields),
    begin(fct_nodefields), end(fct_nodefields) );

  elemfieldnames.insert( end(elemfieldnames),
    begin(fct_elemfieldnames), end(fct_elemfieldnames) );
  elemfields.insert( end(elemfields),
    begin(fct_elemfields), end(fct_elemfields) );

  // Send mesh and fields data (solution dump) for output to file
  d->write( d->Inpoel(), d->Coord(), {}, tk::remap(m_bnode,d->Lid()), {},
            elemfieldnames, nodefieldnames, elemfields, nodefields, c );
}

void
DiagCG::update( const tk::Fields& a, [[maybe_unused]] tk::Fields&& dul )
// *****************************************************************************
// Prepare for next step
//! \param[in] a Limited antidiffusive element contributions
//! \param[in] dul Low order solution increment
// *****************************************************************************
{
  auto d = Disc();

  // Verify that the change in the solution at those nodes where Dirichlet
  // boundary conditions are set is exactly the amount the BCs prescribe
  Assert( correctBC(a,dul,m_bcdir), "Dirichlet boundary condition incorrect" );

  // Apply limited antidiffusive element contributions to low order solution
  if (g_inputdeck.get< tag::discr, tag::fct >())
    m_u = m_ul + a;
  else
    m_u = m_u + m_du;

  // Set symmetry BCs
  //for (const auto& eq : g_cgpde) eq.symbc( m_u, m_bnorm );

  // Compute diagnostics, e.g., residuals
  auto diag_computed = m_diag.compute( *d, m_u );
  // Increase number of iterations and physical time
  d->next();
  // Continue to mesh refinement (if configured)
  if (!diag_computed) refine();
}

void
DiagCG::refine()
// *****************************************************************************
// Optionally refine/derefine mesh
// *****************************************************************************
{
  auto d = Disc();

  auto dtref = g_inputdeck.get< tag::amr, tag::dtref >();
  auto dtfreq = g_inputdeck.get< tag::amr, tag::dtfreq >();

  // if t>0 refinement enabled and we hit the dtref frequency
  if (dtref && !(d->It() % dtfreq)) {   // refine

    // Activate SDAG waits for re-computing the left-hand side
    thisProxy[ thisIndex ].wait4lhs();

    d->startvol();
    d->Ref()->dtref( {}, m_bnode, {} );
    d->refined() = 1;

  } else {      // do not refine

    d->refined() = 0;
    lhs_complete();
    resized();

  }
}

void
DiagCG::resizePostAMR(
  const std::vector< std::size_t >& /*ginpoel*/,
  const tk::UnsMesh::Chunk& chunk,
  const tk::UnsMesh::Coords& coord,
  const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& addedNodes,
  const std::unordered_map< std::size_t, std::size_t >& /*addedTets*/,
  const std::unordered_map< int, std::vector< std::size_t > >& msum,
  const std::map< int, std::vector< std::size_t > >& /*bface*/,
  const std::map< int, std::vector< std::size_t > >& bnode,
  const std::vector< std::size_t >& /*triinpoel*/ )
// *****************************************************************************
//  Receive new mesh from Refiner
//! \param[in] ginpoel Mesh connectivity with global node ids
//! \param[in] chunk New mesh chunk (connectivity and global<->local id maps)
//! \param[in] coord New mesh node coordinates
//! \param[in] addedNodes Newly added mesh nodes and their parents (local ids)
//! \param[in] addedTets Newly added mesh cells and their parents (local ids)
//! \param[in] msum New node communication map
//! \param[in] bnode Boundary-node lists mapped to side set ids
// *****************************************************************************
{
  auto d = Disc();

  // Set flag that indicates that we are during time stepping
  m_initial = 0;

  // Zero field output iteration count between two mesh refinement steps
  d->Itf() = 0;

  // Increase number of iterations with mesh refinement
  ++d->Itr();

  // Resize mesh data structures
  d->resizePostAMR( chunk, coord, msum );

  // Resize auxiliary solution vectors
  auto nelem = d->Inpoel().size()/4;
  auto npoin = coord[0].size();
  auto nprop = m_u.nprop();
  m_u.resize( npoin, nprop );
  m_ul.resize( npoin, nprop );
  m_du.resize( npoin, nprop );
  m_ue.resize( nelem, nprop );
  m_lhs.resize( npoin, nprop );
  m_rhs.resize( npoin, nprop );

  // Update solution on new mesh
  for (const auto& n : addedNodes) {
    for (std::size_t c=0; c<nprop; ++c) {
      m_u(n.first,c,0) = (m_u(n.second[0],c,0) + m_u(n.second[1],c,0))/2.0;
    }
  }

  // Update physical-boundary node lists
  m_bnode = bnode;

  // Resize FCT data structures
  d->FCT()->resize( npoin, msum, d->Bid(), d->Lid(), d->Inpoel() );

  contribute( CkCallback(CkReductionTarget(Transporter,resized), d->Tr()) );
}

void
DiagCG::resized()
// *****************************************************************************
// Resizing data sutrctures after mesh refinement has been completed
// *****************************************************************************
{
  resize_complete();
}

void
DiagCG::out()
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
    writeFields( CkCallback(CkIndex_DiagCG::step(), thisProxy[thisIndex]) );
  else
    step();
}

void
DiagCG::evalLB()
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
DiagCG::evalRestart()
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
DiagCG::step()
// *****************************************************************************
// Evaluate whether to continue with next time step
// *****************************************************************************
{
  auto d = Disc();

  // Output one-liner status report to screen
  d->status();

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

#include "NoWarning/diagcg.def.h"
