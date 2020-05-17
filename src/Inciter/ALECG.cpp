// *****************************************************************************
/*!
  \file      src/Inciter/ALECG.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     ALECG for a PDE system with continuous Galerkin + ALE + RK
  \details   ALECG advances a system of partial differential equations (PDEs)
    using a continuous Galerkin (CG) finite element (FE) spatial discretization
    (using linear shapefunctions on tetrahedron elements) combined with a
    Runge-Kutta (RK) time stepping scheme in the arbitrary Eulerian-Lagrangian
    reference frame.
  \see The documentation in ALECG.h.
*/
// *****************************************************************************

#include "QuinoaConfig.hpp"
#include "ALECG.hpp"
#include "Vector.hpp"
#include "Reader.hpp"
#include "ContainerUtil.hpp"
#include "UnsMesh.hpp"
#include "ExodusIIMeshWriter.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "DerivedData.hpp"
#include "CGPDE.hpp"
#include "Discretization.hpp"
#include "DiagReducer.hpp"
#include "NodeBC.hpp"
#include "Refiner.hpp"
#include "Reorder.hpp"
#include "Around.hpp"
#include "CGPDE.hpp"
#include "Integrate/Mass.hpp"

#ifdef HAS_ROOT
  #include "RootMeshWriter.hpp"
#endif

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< CGPDE > g_cgpde;

//! Runge-Kutta coefficients
static const std::array< tk::real, 3 > rkcoef{{ 1.0/3.0, 1.0/2.0, 1.0 }};

} // inciter::

using inciter::ALECG;

ALECG::ALECG( const CProxy_Discretization& disc,
              const std::map< int, std::vector< std::size_t > >& bface,
              const std::map< int, std::vector< std::size_t > >& bnode,
              const std::vector< std::size_t >& triinpoel ) :
  m_disc( disc ),
  m_initial( 1 ),
  m_nsol( 0 ),
  m_nlhs( 0 ),
  m_ngrad( 0 ),
  m_nrhs( 0 ),
  m_nbnorm( 0 ),
  m_ndfnorm( 0 ),
  m_bnode( bnode ),
  m_bface( bface ),
  m_triinpoel( tk::remap( triinpoel, Disc()->Lid() ) ),
  m_bndel( bndel() ),
  m_dfnorm(),
  m_dfnormc(),
  m_dfn(),
  m_esup( tk::genEsup( Disc()->Inpoel(), 4 ) ),
  m_psup( tk::genPsup( Disc()->Inpoel(), 4, m_esup ) ),
  m_u( m_disc[thisIndex].ckLocal()->Gid().size(),
       g_inputdeck.get< tag::component >().nprop() ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_lhs( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_grad( Disc()->Bid().size(), m_u.nprop()*3 ),
  m_bcdir(),
  m_lhsc(),
  m_gradc(),
  m_rhsc(),
  m_diag(),
  m_bnorm(),
  m_bnormc(),
  m_symbctri(),
  m_symbcnodes(),
  m_farfieldbcnodes(),
  m_stage( 0 ),
  m_boxnodes(),
  m_dtp( m_u.nunk(), 0.0 ),
  m_tp( m_u.nunk(), g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_finished( 0 )
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] bface Boundary-faces mapped to side sets used in the input file
//! \param[in] bnode Boundary-node lists mapped to side sets used in input file
//! \param[in] triinpoel Boundary-face connectivity where BCs set (global ids)
// *****************************************************************************
//! [Constructor]
{
  usesAtSync = true;    // enable migration at AtSync

  // Perform optional operator-access-pattern mesh node reordering
  if (g_inputdeck.get< tag::discr, tag::operator_reorder >()) {

    auto d = Disc();

    // Create new local ids based on access pattern of PDE operators
    std::unordered_map< std::size_t, std::size_t > map;
    std::size_t n = 0;

    for (std::size_t p=0; p<m_u.nunk(); ++p) {  // for each point p
      if (map.find(p) == end(map)) map[p] = n++;
      for (auto q : tk::Around(m_psup,p)) {     // for each edge p-q
        if (map.find(q) == end(map)) map[q] = n++;
      }
    }

    Assert( map.size() == d->Gid().size(), "Map size mismatch" );

    // Remap data in bound Discretization object
    d->remap( map );
    // Recompute elements surrounding points
    m_esup = tk::genEsup( d->Inpoel(), 4 );
    // Recompute points surrounding points
    m_psup = tk::genPsup( d->Inpoel(), 4, m_esup );
    // Remap boundary triangle face connectivity
    tk::remap( m_triinpoel, map );
  }

  // Activate SDAG wait for initially computing the left-hand side and normals
  thisProxy[ thisIndex ].wait4lhs();

  // Signal the runtime system that the workers have been created
  contribute( sizeof(int), &m_initial, CkReduction::sum_int,
    CkCallback(CkReductionTarget(Transporter,comfinal), Disc()->Tr()) );
}
//! [Constructor]

void
ALECG::norm()
// *****************************************************************************
// Start (re-)computing boundare point-, and dual-face normals
// *****************************************************************************
{
  auto d = Disc();

  // Query nodes at which symmetry BCs are specified
  d->bcnodes< tag::bcsym >( m_bface, m_triinpoel, m_symbcnodes );
  // Query nodes at which farfield BCs are specified
  d->bcnodes< tag::bcfarfield >( m_bface, m_triinpoel, m_farfieldbcnodes );

  std::unordered_set< std::size_t > bcnodes;
  std::merge( m_symbcnodes.cbegin(), m_symbcnodes.cend(),
              m_farfieldbcnodes.cbegin(), m_farfieldbcnodes.cend(),
              std::inserter( bcnodes, begin(bcnodes) ) );

  // Compute boundary point normals
  bnorm( std::move(bcnodes) );

  // Compute dual-face normals associated to edges
  dfnorm();
}

std::vector< std::size_t >
ALECG::bndel() const
// *****************************************************************************
// Find elements along our mesh chunk boundary
//! \return List of local element ids that have at least a single node
//!   contributing to a chare boundary
// *****************************************************************************
{
  // Lambda to find out if a mesh node is shared with another chare
  auto shared = [this]( std::size_t i ){
    for (const auto& [c,n] : Disc()->NodeCommMap())
      if (n.find(i) != end(n)) return true;
    return false;
  };

  // Find elements along our mesh chunk boundary
  std::vector< std::size_t > e;
  const auto& inpoel = Disc()->Inpoel();
  const auto gid = Disc()->Gid();
  for (std::size_t n=0; n<inpoel.size(); ++n)
    if (shared( gid[ inpoel[n] ] )) e.push_back( n/4 );
  tk::unique( e );

  return e;
}

std::array< tk::real, 3 >
ALECG::edfnorm( const tk::UnsMesh::Edge& edge,
                const std::unordered_map< tk::UnsMesh::Edge,
                        std::vector< std::size_t >,
                        tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> >& esued )
// *****************************************************************************
//  Compute normal of dual-mesh associated to edge
//! \param[in] edge Edge whose dual-face normal to compute given by local ids
//! \param[in] inpoel Mesh element connectivity
//! \param[in] esued Elements surrounding edges
//! \return Dual-face normal for edge
// *****************************************************************************
{
  auto d = Disc();
  const auto& inpoel = d->Inpoel();
  const auto& coord = d->Coord();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  std::array< tk::real, 3 > n{ 0.0, 0.0, 0.0 };

  for (auto e : tk::cref_find(esued,edge)) {
    // access node IDs
    const std::array< std::size_t, 4 >
      N{ inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };
    // compute element Jacobi determinant
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );        // J = 6V
    Assert( J > 0, "Element Jacobian non-positive" );
    // shape function derivatives, nnode*ndim [4][3]
    std::array< std::array< tk::real, 3 >, 4 > grad;
    grad[1] = tk::crossdiv( ca, da, J );
    grad[2] = tk::crossdiv( da, ba, J );
    grad[3] = tk::crossdiv( ba, ca, J );
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];
    // sum normal contributions
    auto J48 = J/48.0;
    for (const auto& [a,b] : tk::lpoed) {
      auto s = tk::orient( {N[a],N[b]}, edge );
      for (std::size_t j=0; j<3; ++j)
        n[j] += J48 * s * (grad[a][j] - grad[b][j]);
    }
  }

  return n;
}

void
ALECG::dfnorm()
// *****************************************************************************
// Compute dual-face normals associated to edges
// *****************************************************************************
{
  auto d = Disc();
  const auto& inpoel = d->Inpoel();
  const auto& gid = d->Gid();

  // compute derived data structures
  auto esued = tk::genEsued( inpoel, 4, tk::genEsup( inpoel, 4 ) );

  // Compute dual-face normals for domain edges
  for (std::size_t p=0; p<gid.size(); ++p)    // for each point p
    for (auto q : tk::Around(m_psup,p))       // for each edge p-q
      if (gid[p] < gid[q])
        m_dfnorm[{gid[p],gid[q]}] = edfnorm( {p,q}, esued );

  // Send our dual-face normal contributions to neighbor chares
  if (d->EdgeCommMap().empty())
    comdfnorm_complete();
  else {
    for (const auto& [c,edges] : d->EdgeCommMap()) {
      decltype(m_dfnorm) exp;
      for (const auto& e : edges) exp[e] = tk::cref_find(m_dfnorm,e);
      thisProxy[c].comdfnorm( exp );
    }
  }

  owndfnorm_complete();
}

void
ALECG::comdfnorm( const std::unordered_map< tk::UnsMesh::Edge,
                    std::array< tk::real, 3 >,
                    tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> >& dfnorm )
// *****************************************************************************
// Receive contributions to dual-face normals on chare-boundaries
//! \param[in] dfnorm Incoming partial sums of dual-face normals associated to
//!   chare-boundary edges
// *****************************************************************************
{
  // Buffer up inccoming contributions to dual-face normals
  for (const auto& [e,n] : dfnorm) {
    auto& dfn = m_dfnormc[e];
    dfn[0] += n[0];
    dfn[1] += n[1];
    dfn[2] += n[2];
  }

  if (++m_ndfnorm == Disc()->EdgeCommMap().size()) {
    m_ndfnorm = 0;
    comdfnorm_complete();
  }
}

void
ALECG::bnorm( std::unordered_set< std::size_t >&& bcnodes )
// *****************************************************************************
//  Compute boundary point normals
//! \param[in] Local node ids at which BCs are set that require normals
// *****************************************************************************
{
  auto d = Disc();

  m_bnorm = cg::bnorm( m_bface, m_triinpoel, d->Coord(), d->Gid(), bcnodes );

  // Send our nodal normal contributions to neighbor chares
  if (d->NodeCommMap().empty())
   comnorm_complete();
  else
    for (const auto& [ neighborchare, sharednodes ] : d->NodeCommMap()) {
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
ALECG::comnorm(
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

  if (++m_nbnorm == Disc()->NodeCommMap().size()) {
    m_nbnorm = 0;
    comnorm_complete();
  }
}

void
ALECG::registerReducers()
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
ALECG::ResumeFromSync()
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
ALECG::setup()
// *****************************************************************************
// Setup rows, query boundary conditions, output mesh, etc.
// *****************************************************************************
{
  auto d = Disc();

  // Set initial conditions for all PDEs
  for (auto& eq : g_cgpde) eq.initialize( d->Coord(), m_u, d->T(), m_boxnodes );

  // Compute volume of user-defined box IC
  d->boxvol( m_boxnodes );

  // Query time history field output labels from all PDEs integrated
  const auto& hist_points = g_inputdeck.get< tag::history, tag::point >();
  if (!hist_points.empty()) {
    std::vector< std::string > histnames;
    for (const auto& eq : g_cgpde) {
      auto n = eq.histNames();
      histnames.insert( end(histnames), begin(n), end(n) );
    }
    d->histheader( std::move(histnames) );
  }
}

void
ALECG::boxvol( tk::real v )
// *****************************************************************************
// Receive total box IC volume
//! \param[in] v Total volume within user-specified box
// *****************************************************************************
{
  auto d = Disc();

  // Store user-defined box IC volume
  d->Boxvol() = v;

  // Update density in user-defined IC box based on box volume
  for (const auto& eq : g_cgpde) eq.box( d->Boxvol(), m_boxnodes, m_u );

  // Compute left-hand side of PDEs
  lhs();
}

//! [start]
void
ALECG::start()
// *****************************************************************************
// Start time stepping
// *****************************************************************************
{
  // Start timer measuring time stepping wall clock time
  Disc()->Timer().zero();
  // Zero grind-timer
  Disc()->grindZero();
  // Continue to next time step
  next();
}
//! [start]

//! [Compute own and send lhs on chare-boundary]
void
ALECG::lhs()
// *****************************************************************************
// Compute the left-hand side of transport equations
//! \details Also (re-)compute all data structures if the mesh changed.
// *****************************************************************************
{
  auto d = Disc();

  // Compute lumped mass lhs
  m_lhs = tk::lump( m_u.nprop(), d->Coord(), d->Inpoel() );

  if (d->NodeCommMap().empty())        // in serial we are done
    comlhs_complete();
  else // send contributions of lhs to chare-boundary nodes to fellow chares
    for (const auto& [c,n] : d->NodeCommMap()) {
      std::vector< std::vector< tk::real > > l( n.size() );
      std::size_t j = 0;
      for (auto i : n) l[ j++ ] = m_lhs[ tk::cref_find(d->Lid(),i) ];
      thisProxy[c].comlhs( std::vector<std::size_t>(begin(n),end(n)), l );
    }

  ownlhs_complete();

  // (Re-)compute boundary point-, and dual-face normals
  norm();
}
//! [Compute own and send lhs on chare-boundary]

//! [Receive lhs on chare-boundary]
void
ALECG::comlhs( const std::vector< std::size_t >& gid,
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

  auto d = Disc();

  for (std::size_t i=0; i<gid.size(); ++i) {
    m_lhsc[ gid[i] ] += L[i];
  }

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_nlhs == d->NodeCommMap().size()) {
    m_nlhs = 0;
    comlhs_complete();
  }
}
//! [Receive lhs on chare-boundary]

//! [Merge lhs and continue]
void
ALECG::lhsmerge()
// *****************************************************************************
// The own and communication portion of the left-hand side is complete
// *****************************************************************************
{
  // Combine own and communicated contributions to left hand side
  auto d = Disc();

  // Combine own and communicated contributions to LHS and ICs
  for (const auto& b : m_lhsc) {
    auto lid = tk::cref_find( d->Lid(), b.first );
    for (ncomp_t c=0; c<m_lhs.nprop(); ++c)
      m_lhs(lid,c,0) += b.second[c];
  }

  // Clear receive buffer
  tk::destroy(m_lhsc);

  // Combine own and communicated contributions of normals
  normfinal();

  // Continue after lhs is complete
  if (m_initial) {
    // Output initial conditions to file
    writeFields( CkCallback(CkIndex_ALECG::start(), thisProxy[thisIndex]) );
  } else {
    lhs_complete();
  }
}
//! [Merge lhs and continue]

void
ALECG::normfinal()
// *****************************************************************************
//  Finish computing dual-face and boundary point normals
// *****************************************************************************
{
  auto d = Disc();
  const auto& lid = d->Lid();

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
  tk::destroy( m_bnormc );

  // Divie summed point normals by the sum of inverse distance squared
  for (auto& [p,n] : m_bnorm) {
    n[0] /= n[3];
    n[1] /= n[3];
    n[2] /= n[3];
    Assert( (n[0]*n[0] + n[1]*n[1] + n[2]*n[2] - 1.0) <
            std::numeric_limits< tk::real >::epsilon(), "Non-unit normal" );
  }

  // Replace global->local ids associated to boundary point normals
  decltype(m_bnorm) bnorm;
  for (auto&& [g,n] : m_bnorm) bnorm[ tk::cref_find(lid,g) ] = std::move(n);
  m_bnorm = std::move(bnorm);

  // Apply symmetry boundary conditions on initial conditions
  for (const auto& eq : g_cgpde)
    eq.symbc( m_u, m_bnorm, m_symbcnodes );
  // Apply farfield boundary conditions on initial conditions
  for (const auto& eq : g_cgpde)
    eq.farfieldbc( m_u, m_bnorm, m_farfieldbcnodes );

  // Flatten boundary normal data structure
  m_symbctri.resize( m_triinpoel.size()/3, 0 );
  for (std::size_t e=0; e<m_triinpoel.size()/3; ++e)
    if (m_symbcnodes.find(m_triinpoel[e*3+0]) != end(m_symbcnodes))
      m_symbctri[e] = 1;

  // Count contributions to chare-boundary edges
  std::unordered_map< tk::UnsMesh::Edge, std::size_t,
    tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> > edge_node_count;
  for (const auto& [c,edges] : d->EdgeCommMap())
    for (const auto& e : edges)
      ++edge_node_count[e];

  // Combine and weigh communicated contributions to dual-face normals
  for (auto& [e,n] : m_dfnormc) {
    const auto& dfn = tk::cref_find( m_dfnorm, e );
    n[0] += dfn[0];
    n[1] += dfn[1];
    n[2] += dfn[2];
    auto count = static_cast< tk::real >( tk::cref_find( edge_node_count, e ) );
    auto factor = 1.0/(count + 1.0);
    for (auto & x : n) x *= factor;
  }
  // Convert dual-face normals to streamable (and vectorizable) data structure
  m_dfn.resize( m_psup.first.size() * 6 );      // 2 vectors per access
  const auto& gid = d->Gid();
  for (std::size_t p=0,k=0; p<m_u.nunk(); ++p)
    for (auto q : tk::Around(m_psup,p)) {
      std::array< std::size_t, 2 > e{ gid[p], gid[q] };
      auto n = tk::cref_find( m_dfnorm, e );
      // figure out if this is an edge on the parallel boundary
      auto nit = m_dfnormc.find( e );
      auto m = ( nit != m_dfnormc.end() ) ? nit->second : n;
      // orient normals
      if (gid[p] > gid[q]) { tk::flip(n); tk::flip(m); }
      m_dfn[k+0] = n[0];
      m_dfn[k+1] = n[1];
      m_dfn[k+2] = n[2];
      m_dfn[k+3] = m[0];
      m_dfn[k+4] = m[1];
      m_dfn[k+5] = m[2];
      k += 6;
    }

  tk::destroy( m_dfnorm );
  tk::destroy( m_dfnormc );
}

void
ALECG::next()
// *****************************************************************************
// Continue to next time step
// *****************************************************************************
{
  dt();
}

void
ALECG::dt()
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

    //! [Find the minimum dt across all PDEs integrated]
    if (g_inputdeck.get< tag::discr, tag::steady_state >()) {

      // compute new dt for each mesh point
      for (const auto& eq : g_cgpde)
        eq.dt( d->It(), d->Vol(), m_u, m_dtp );

      // find the smallest dt of all nodes on this chare
      mindt = *std::min_element( begin(m_dtp), end(m_dtp) );

    } else {    // compute new dt for this chare

      // find the smallest dt of all equations on this chare
      for (const auto& eq : g_cgpde) {
        auto eqdt = eq.dt( d->Coord(), d->Inpoel(), m_u );
        if (eqdt < mindt) mindt = eqdt;
      }

    }
    //! [Find the minimum dt across all PDEs integrated]

  }

  //! [Advance]
  // Actiavate SDAG waits for next time step stage
  thisProxy[ thisIndex ].wait4grad();
  thisProxy[ thisIndex ].wait4rhs();
  thisProxy[ thisIndex ].wait4stage();

  // Contribute to minimum dt across all chares the advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(ALECG,advance), thisProxy) );
  //! [Advance]
}

void
ALECG::advance( tk::real newdt )
// *****************************************************************************
// Advance equations to next time step
//! \param[in] newdt The smallest dt across the whole problem
// *****************************************************************************
{
  auto d = Disc();

  // Set new time step size
  if (m_stage == 0) d->setdt( newdt );

  // Compute gradients for next time step
  grad();
}

void
ALECG::grad()
// *****************************************************************************
// Compute nodal gradients
// *****************************************************************************
{
  auto d = Disc();

  // Compute own portion of gradients for all equations
  for (const auto& eq : g_cgpde)
    eq.grad(d->Coord(), d->Inpoel(), m_bndel, d->Gid(), d->Bid(), m_u, m_grad);

  // Communicate gradients to other chares on chare-boundary
  if (d->NodeCommMap().empty())        // in serial we are done
    comgrad_complete();
  else // send gradients contributions to chare-boundary nodes to fellow chares
    for (const auto& [c,n] : d->NodeCommMap()) {
      std::vector< std::vector< tk::real > > g( n.size() );
      std::size_t j = 0;
      for (auto i : n) g[ j++ ] = m_grad[ tk::cref_find(d->Bid(),i) ];
      thisProxy[c].comgrad( std::vector<std::size_t>(begin(n),end(n)), g );
    }

  owngrad_complete();
}

void
ALECG::comgrad( const std::vector< std::size_t >& gid,
                const std::vector< std::vector< tk::real > >& G )
// *****************************************************************************
//  Receive contributions to nodal gradients on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive grad contributions
//! \param[in] G Partial contributions of gradients to chare-boundary nodes
//! \details This function receives contributions to m_grad, which stores the
//!   nodal gradients at mesh nodes. While m_grad stores own
//!   contributions, m_gradc collects the neighbor chare contributions during
//!   communication. This way work on m_grad and m_gradc is overlapped. The two
//!   are combined in rhs().
// *****************************************************************************
{
  Assert( G.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  for (std::size_t i=0; i<gid.size(); ++i) m_gradc[ gid[i] ] += G[i];

  if (++m_ngrad == Disc()->NodeCommMap().size()) {
    m_ngrad = 0;
    comgrad_complete();
  }
}

void
ALECG::rhs()
// *****************************************************************************
// Compute right-hand side of transport equations
// *****************************************************************************
{
  auto d = Disc();

  // Combine own and communicated contributions to nodal gradients
  for (const auto& [gid,g] : m_gradc) {
    auto bid = tk::cref_find( d->Bid(), gid );
    for (ncomp_t c=0; c<m_grad.nprop(); ++c) m_grad(bid,c,0) += g[c];
  }

  // clear gradients receive buffer
  tk::destroy(m_gradc);

  const auto steady = g_inputdeck.get< tag::discr, tag::steady_state >();

  // Compute own portion of right-hand side for all equations
  auto prev_rkcoef = m_stage == 0 ? 0.0 : rkcoef[m_stage-1];
  if (steady)
    for (std::size_t p=0; p<m_tp.size(); ++p) m_tp[p] += prev_rkcoef * m_dtp[p];
  for (const auto& eq : g_cgpde)
    eq.rhs( d->T() + prev_rkcoef * d->Dt(), d->Coord(), d->Inpoel(),
            m_triinpoel, d->Bid(), d->Lid(), m_dfn, m_psup, m_esup,
            m_symbctri, d->Vol(), m_grad, m_u, m_tp, m_rhs );
  if (steady)
    for (std::size_t p=0; p<m_tp.size(); ++p) m_tp[p] -= prev_rkcoef * m_dtp[p];

  // Query and match user-specified boundary conditions to side sets
  if (steady) for (auto& deltat : m_dtp) deltat *= rkcoef[m_stage];
  m_bcdir = match( m_u.nprop(), d->T(), rkcoef[m_stage] * d->Dt(),
                   m_tp, m_dtp, d->Coord(), d->Lid(), m_bnode );
  if (steady) for (auto& deltat : m_dtp) deltat /= rkcoef[m_stage];

  // Communicate rhs to other chares on chare-boundary
  if (d->NodeCommMap().empty())        // in serial we are done
    comrhs_complete();
  else // send contributions of rhs to chare-boundary nodes to fellow chares
    for (const auto& [c,n] : d->NodeCommMap()) {
      std::vector< std::vector< tk::real > > r( n.size() );
      std::size_t j = 0;
      for (auto i : n) r[ j++ ] = m_rhs[ tk::cref_find(d->Lid(),i) ];
      thisProxy[c].comrhs( std::vector<std::size_t>(begin(n),end(n)), r );
    }

  ownrhs_complete();
}

void
ALECG::comrhs( const std::vector< std::size_t >& gid,
               const std::vector< std::vector< tk::real > >& R )
// *****************************************************************************
//  Receive contributions to right-hand side vector on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive RHS contributions
//! \param[in] R Partial contributions of RHS to chare-boundary nodes
//! \details This function receives contributions to m_rhs, which stores the
//!   right hand side vector at mesh nodes. While m_rhs stores own
//!   contributions, m_rhsc collects the neighbor chare contributions during
//!   communication. This way work on m_rhs and m_rhsc is overlapped. The two
//!   are combined in solve().
// *****************************************************************************
{
  Assert( R.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  for (std::size_t i=0; i<gid.size(); ++i) m_rhsc[ gid[i] ] += R[i];

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_nrhs == Disc()->NodeCommMap().size()) {
    m_nrhs = 0;
    comrhs_complete();
  }
}

void
ALECG::solve()
// *****************************************************************************
//  Solve low and high order diagonal systems
// *****************************************************************************
{
  const auto ncomp = m_rhs.nprop();

  auto d = Disc();

  // Combine own and communicated contributions to rhs
  for (const auto& b : m_rhsc) {
    auto lid = tk::cref_find( d->Lid(), b.first );
    for (ncomp_t c=0; c<ncomp; ++c) m_rhs(lid,c,0) += b.second[c];
  }

  // clear receive buffer
  tk::destroy(m_rhsc);

  const auto steady = g_inputdeck.get< tag::discr, tag::steady_state >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();

  // Set Dirichlet BCs for lhs and rhs
  for (const auto& [b,bc] : m_bcdir)
    for (ncomp_t c=0; c<ncomp; ++c)
      if (bc[c].first) {
        m_lhs(b,c,0) = 1.0;
        auto deltat = steady ? m_dtp[b] : d->Dt();
        m_rhs(b,c,0) = bc[c].second / deltat / rkcoef[m_stage];
      }

  // Update Un
  if (m_stage == 0) m_un = m_u;

  // Solve the sytem
  if (steady) {

    for (std::size_t i=0; i<m_u.nunk(); ++i)
      for (ncomp_t c=0; c<m_u.nprop(); ++c)
        m_u(i,c,0) = m_un(i,c,0)
          + rkcoef[m_stage] * m_dtp[i] * m_rhs(i,c,0) / m_lhs(i,c,0);

  } else {

    m_u = m_un + rkcoef[m_stage] * d->Dt() * m_rhs / m_lhs;

  }

  // Apply farfield boundary conditions on new solution
  for (const auto& eq : g_cgpde)
    eq.farfieldbc( m_u, m_bnorm, m_farfieldbcnodes );

  //! [Continue after solve]
  if (m_stage < 2) {

    // Activate SDAG wait for next time step stage
    thisProxy[ thisIndex ].wait4grad();
    thisProxy[ thisIndex ].wait4rhs();

    // continue with next time step stage
    stage();

  } else {

    // Compute diagnostics, e.g., residuals
    auto diag_computed = m_diag.compute( *d, m_u, m_un );
    // Increase number of iterations and physical time
    d->next();

    // Advance physical time for local time stepping
    if (steady) {
      const auto term = g_inputdeck.get< tag::discr, tag::term >();
      for (std::size_t i=0; i<m_u.nunk(); ++i) {
        if (std::abs(m_dtp[i]) > eps)
          m_tp[i] += m_dtp[i];
        else
          m_tp[i] = term;
      }
    }

    // Continue to mesh refinement (if configured)
    if (!diag_computed) refine( 1.0 );
  }
  //! [Continue after solve]
}

void
ALECG::refine( tk::real l2res )
// *****************************************************************************
// Optionally refine/derefine mesh
//! \param[in] l2res L2-norm of the residual across the whole problem
// *****************************************************************************
{
  //! [Refine]
  auto d = Disc();

  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto steady = g_inputdeck.get< tag::discr, tag::steady_state >();
  const auto residual = g_inputdeck.get< tag::discr, tag::residual >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();

  if (steady) {

    // this is the last time step if max time of max number of time steps
    // reached or the residual has reached its convergence criterion
    if (std::abs(d->T()-term) < eps || d->It() >= nstep || l2res < residual)
      m_finished = 1;

  } else {

    // this is the last time step if max time or max iterations reached
    if (std::abs(d->T()-term) < eps || d->It() >= nstep) m_finished = 1;

  }

  auto dtref = g_inputdeck.get< tag::amr, tag::dtref >();
  auto dtfreq = g_inputdeck.get< tag::amr, tag::dtfreq >();

  // if t>0 refinement enabled and we hit the frequency
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
  //! [Refine]
}

//! [Resize]
void
ALECG::resizePostAMR(
  const std::vector< std::size_t >& /*ginpoel*/,
  const tk::UnsMesh::Chunk& chunk,
  const tk::UnsMesh::Coords& coord,
  const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& addedNodes,
  const std::unordered_map< std::size_t, std::size_t >& /*addedTets*/,
  const tk::NodeCommMap& nodeCommMap,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& bnode,
  const std::vector< std::size_t >& triinpoel )
// *****************************************************************************
//  Receive new mesh from Refiner
//! \param[in] ginpoel Mesh connectivity with global node ids
//! \param[in] chunk New mesh chunk (connectivity and global<->local id maps)
//! \param[in] coord New mesh node coordinates
//! \param[in] addedNodes Newly added mesh nodes and their parents (local ids)
//! \param[in] addedTets Newly added mesh cells and their parents (local ids)
//! \param[in] nodeCommMap New node communication map
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

  // Resize mesh data structures
  d->resizePostAMR( chunk, coord, nodeCommMap );

  // Resize auxiliary solution vectors
  auto npoin = coord[0].size();
  auto nprop = m_u.nprop();
  m_u.resize( npoin, nprop );
  m_un.resize( npoin, nprop );
  m_lhs.resize( npoin, nprop );
  m_rhs.resize( npoin, nprop );
  m_grad.resize( d->Bid().size(), nprop*3 );

  // Update solution on new mesh
  for (const auto& n : addedNodes)
    for (std::size_t c=0; c<nprop; ++c)
      m_u(n.first,c,0) = (m_u(n.second[0],c,0) + m_u(n.second[1],c,0))/2.0;

  // Update physical-boundary node-, face-, and element lists
  m_bnode = bnode;
  m_bface = bface;
  m_triinpoel = tk::remap( triinpoel, d->Lid() );

  contribute( CkCallback(CkReductionTarget(Transporter,resized), d->Tr()) );
}
//! [Resize]

void
ALECG::resized()
// *****************************************************************************
// Resizing data sutrctures after mesh refinement has been completed
// *****************************************************************************
{
  resize_complete();
}

void
ALECG::stage()
// *****************************************************************************
// Evaluate whether to continue with next time step stage
// *****************************************************************************
{
  // Increment Runge-Kutta stage counter
  ++m_stage;

  // if not all Runge-Kutta stages complete, continue to next time stage,
  // otherwise output field data to file(s)
  if (m_stage < 3) grad(); else out();
}

void
ALECG::writeFields( CkCallback c ) const
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::benchmark >()) {

    c.send();

  } else {

    auto d = Disc();

    // Query and collect block and surface field names from PDEs integrated
    std::vector< std::string > nodefieldnames;
    std::vector< std::string > nodesurfnames;
    for (const auto& eq : g_cgpde) {
      auto n = eq.fieldNames();
      nodefieldnames.insert( end(nodefieldnames), begin(n), end(n) );
      auto s = eq.surfNames();
      nodesurfnames.insert( end(nodesurfnames), begin(s), end(s) );
    }

    // Collect node block and surface field solution
    auto u = m_u;
    std::vector< std::vector< tk::real > > nodefields;
    std::vector< std::vector< tk::real > > nodesurfs;
    for (const auto& eq : g_cgpde) {
      auto o = eq.fieldOutput( d->T(), d->meshvol(), d->Coord()[0].size(),
                               d->Coord(), d->V(), u );
      nodefields.insert( end(nodefields), begin(o), end(o) );
      auto s = eq.surfOutput( tk::bfacenodes(m_bface,m_triinpoel), u );
      nodesurfs.insert( end(nodesurfs), begin(s), end(s) );
    }

    // Create volume node field assigning 1 where symmetry BC is set
    // std::unordered_set< std::size_t > symbcnodes;
    // d->bcnodes< tag::bcsym >( m_bface, m_triinpoel, bcnodes );
    // nodefieldnames.push_back( "bc_type" );
    // nodefields.push_back( std::vector<tk::real>(d->Coord()[0].size(),0.0) );
    // for (auto i : symbcnodes) nodefields.back()[i] = 1.0;

    // nodefieldnames.push_back( "tp" );
    // nodefieldnames.push_back( "dtp" );
    // nodefields.push_back( m_tp );
    // nodefields.push_back( m_dtp );

    Assert( nodefieldnames.size() == nodefields.size(), "Size mismatch" );

    // Send mesh and fields data (solution dump) for output to file
    d->write( d->Inpoel(), d->Coord(), m_bface, tk::remap(m_bnode,d->Lid()),
              m_triinpoel, {}, nodefieldnames, nodesurfnames, {}, nodefields,
              nodesurfs, c );

  }
}

void
ALECG::out()
// *****************************************************************************
// Output mesh field data
// *****************************************************************************
{
  auto d = Disc();

  // Output time history if we hit its output frequency
  const auto histfreq = g_inputdeck.get< tag::interval, tag::history >();
  if ( !((d->It()) % histfreq) ) {
    std::vector< std::vector< tk::real > > hist;
    for (const auto& eq : g_cgpde) {
      auto h = eq.histOutput( d->Hist(), d->Inpoel(), m_u );
      hist.insert( end(hist), begin(h), end(h) );
    }
    d->history( std::move(hist) );
  }

  const auto fieldfreq = g_inputdeck.get< tag::interval, tag::field >();

  // output field data if field iteration count is reached or in the last time
  // step
  if ( !((d->It()) % fieldfreq) || m_finished )
    writeFields( CkCallback(CkIndex_ALECG::step(), thisProxy[thisIndex]) );
  else
    step();
}

void
ALECG::evalLB( int nrestart )
// *****************************************************************************
// Evaluate whether to do load balancing
//! \param[in] nrestart Number of times restarted
// *****************************************************************************
{
  auto d = Disc();

  // Detect if just returned from a checkpoint and if so, zero timers and
  // finished flag
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
ALECG::evalRestart()
// *****************************************************************************
// Evaluate whether to save checkpoint/restart
// *****************************************************************************
{
  auto d = Disc();

  const auto rsfreq = g_inputdeck.get< tag::cmd, tag::rsfreq >();
  const auto benchmark = g_inputdeck.get< tag::cmd, tag::benchmark >();

  if ( !benchmark && (d->It()) % rsfreq == 0 ) {

    int finished = 0;
    d->contribute( sizeof(int), &finished, CkReduction::nop,
      CkCallback(CkReductionTarget(Transporter,checkpoint), d->Tr()) );

  } else {

    evalLB( /* nrestart = */ -1 );

  }
}

void
ALECG::step()
// *****************************************************************************
// Evaluate whether to continue with next time step
// *****************************************************************************
{
  auto d = Disc();

  // Output one-liner status report to screen
  d->status();
  // Reset Runge-Kutta stage counter
  m_stage = 0;

  if (not m_finished) {

    evalRestart();

  } else {

    d->contribute( CkCallback(CkReductionTarget(Transporter,finish), d->Tr()) );

  }
}

#include "NoWarning/alecg.def.h"
