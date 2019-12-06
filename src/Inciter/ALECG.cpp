// *****************************************************************************
/*!
  \file      src/Inciter/ALECG.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
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
#include "HashMapReducer.hpp"
#include "Integrate/Mass.hpp"

#ifdef HAS_ROOT
  #include "RootMeshWriter.hpp"
#endif

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< CGPDE > g_cgpde;

//! Runge-Kutta coefficients
static const std::array< std::array< tk::real, 3 >, 2 >
  rkcoef{{ {{ 0.0, 3.0/4.0, 1.0/3.0 }}, {{ 1.0, 1.0/4.0, 2.0/3.0 }} }};

static CkReduction::reducerType BndEdgeMerger;

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
  m_nnorm( 0 ),
  m_ndfnorm( 0 ),
  m_bnode( bnode ),
  m_triinpoel( tk::remap(triinpoel,Disc()->Lid()) ),
  m_esued(tk::genEsued(Disc()->Inpoel(), 4, tk::genEsup(Disc()->Inpoel(),4))),
  m_psup(tk::genPsup(Disc()->Inpoel(), 4, tk::genEsup(Disc()->Inpoel(),4))),
  m_dfnorm(),
  m_u( m_disc[thisIndex].ckLocal()->Gid().size(),
       g_inputdeck.get< tag::component >().nprop() ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_lhs( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_grad( m_u.nunk(), m_u.nprop()*3 ),
  m_bndEdges(),
  m_bcdir(),
  m_lhsc(),
  m_gradc(),
  m_rhsc(),
  m_diag(),
  m_bnorm(),
  m_bnormc(),
  m_dfnormc(),
  m_stage( 0 )
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] bface Boundary-faces mapped to side set ids where BCs set
//! \param[in] bnode Boundary-node lists mapped to side set ids where BCs set
//! \param[in] triinpoel Boundary-face connectivity where BCs set (global ids)
// *****************************************************************************
//! [Constructor]
{
  usesAtSync = true;    // enable migration at AtSync

  // Activate SDAG wait for initially computing the left-hand side and normals
  thisProxy[ thisIndex ].wait4norm();
  thisProxy[ thisIndex ].wait4lhs();

  // Query nodes at which symmetry BCs are specified
  std::unordered_set< std::size_t > symbcnodes;
  for (const auto& eq : g_cgpde) eq.symbcnodes( bface, triinpoel, symbcnodes );

  // Compute boundary point normals
  bnorm( bface, triinpoel, std::move(symbcnodes) );

  auto d = Disc();
  const auto& gid = d->Gid();
  const auto& inpoel = d->Inpoel();

  // Generate boundary edges of our mesh chunk
  tk::UnsMesh::EdgeSet bnded;
  auto esup = tk::genEsup( inpoel, 4 );         // elements surrounding points
  auto esuel = tk::genEsuelTet( inpoel, esup ); // elems surrounding elements
  for (std::size_t e=0; e<esuel.size()/4; ++e) {
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f) {
      if (esuel[mark+f] == -1) {
        auto A = gid[ inpoel[ mark+tk::lpofa[f][0] ] ];
        auto B = gid[ inpoel[ mark+tk::lpofa[f][1] ] ];
        auto C = gid[ inpoel[ mark+tk::lpofa[f][2] ] ];
        bnded.insert( {A,B} );
        bnded.insert( {B,C} );
        bnded.insert( {C,A} );
      }
    }
  }

  // Aggregate boundary edges across all chares
  decltype(m_bndEdges) bnd{{ thisIndex, std::move(bnded) }};
  auto stream = tk::serialize( bnd );
  contribute( stream.first, stream.second.get(), BndEdgeMerger,
    CkCallback(CkIndex_ALECG::addBndEdges(nullptr),thisProxy) );
}
//! [Constructor]

void
ALECG::addBndEdges( CkReductionMsg* msg )
// *****************************************************************************
//! Receive boundary edges from all refiner chares (including this one)
//! \param[in] msg Charm++ message containing the aggregated map of bnd edges
// *****************************************************************************
{
  decltype(m_bndEdges) allBndEdges;
  PUP::fromMem creator( msg->getData() );
  creator | allBndEdges;
  delete msg;

  // Compute unique set of chares that share at least a single edge with us
  const auto& ownedges = tk::cref_find( allBndEdges, thisIndex );
  for (const auto& [c,edges] : allBndEdges)
    if (c != thisIndex)
      for (const auto& e : edges)
        if (ownedges.find(e) != end(ownedges))
          m_bndEdges[c].insert(e);

  // Compute dual-face normals associated to edges
  dfnorm();
}

void
ALECG::dfnorm()
// *****************************************************************************
// Compute dual-face normals associated to edges
// *****************************************************************************
{
  auto d = Disc();
  const auto& inpoel = d->Inpoel();
  const auto& coord = d->Coord();
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];
  const auto& gid = d->Gid();

  // compute own portion of dual-face normals
  for (std::size_t p=0; p<m_u.nunk(); ++p) {  // for each point p
    for (auto q : tk::Around(m_psup,p)) {     // for each edge p-q
      if (gid[p] < gid[q]) {
        // compute normal of dual-mesh associated to edge p-q
        auto& n = m_dfnorm[ {gid[p],gid[q]} ];
        n = { 0.0, 0.0, 0.0 };
        for (auto e : tk::cref_find(m_esued,{p,q})) {
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
          // sum normal contributions to nodes
          auto J48 = J/48.0;
          for (const auto& [a,b] : tk::lpoed) {
            auto s = tk::orient( {N[a],N[b]}, {p,q} );
            for (std::size_t j=0; j<3; ++j)
              n[j] += J48 * s * (grad[a][j] - grad[b][j]);
          }
        }
        Assert( std::abs(tk::dot(n,n)) >
                  std::numeric_limits< tk::real >::epsilon(),
                "Dual-face normal zero length" );
      }
    }
  }

  // Send our dual-face normal contributions to neighbor chares
  if (m_bndEdges.empty())
    comdfnorm_complete();
  else
    for (const auto& [c,edges] : m_bndEdges) {
      decltype(m_dfnorm) exp;
      for (const auto& e : edges) exp[e] = tk::cref_find(m_dfnorm,e);
      thisProxy[c].comdfnorm( exp );
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

  if (++m_ndfnorm == m_bndEdges.size()) {
    m_ndfnorm = 0;
    comdfnorm_complete();
  }
}

void
ALECG::bnorm( const std::map< int, std::vector< std::size_t > >& bface,
              const std::vector< std::size_t >& triinpoel,
              std::unordered_set< std::size_t >&& symbcnodes )
// *****************************************************************************
//  Compute boundary point normals
//! \param[in] bface Boundary faces side-set information
//! \param[in] triinpoel Boundary triangle face connecitivity
//! \param[in] Node ids at which symmetry BCs are set
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

  if (++m_nnorm == Disc()->NodeCommMap().size()) {
    m_nnorm = 0;
    comnorm_complete();
  }
}

void
ALECG::normfinal()
// *****************************************************************************
//  Finish computing dual-face and boundary point normals
// *****************************************************************************
{
  const auto& lid = Disc()->Lid();

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
  for (auto& [ p, n ] : m_bnorm) {
    n[0] /= n[3];
    n[1] /= n[3];
    n[2] /= n[3];
    Assert( (n[0]*n[0] + n[1]*n[1] + n[2]*n[2] - 1.0) <
            std::numeric_limits< tk::real >::epsilon(), "Non-unit normal" );
  }

  // Replace global->local ids associated to boundary normals of symbc nodes
  decltype(m_bnorm) bnorm;
  for (auto&& [g,n] : m_bnorm) bnorm[ tk::cref_find(lid,g) ] = std::move(n);
  m_bnorm = std::move(bnorm);

  // Combine own and communicated contributions to dual-face normals
  for (auto& [e,n] : m_dfnormc) {
    auto& dfn = tk::ref_find(m_dfnorm,e);
    dfn[0] += n[0];
    dfn[1] += n[1];
    dfn[2] += n[2];
  }
  tk::destroy( m_dfnormc );

  // Normalize dual-face normals
  for (auto& [e,n] : m_dfnorm) tk::unit(n);

  // Signal the runtime system that the workers have been created
  contribute( sizeof(int), &m_initial, CkReduction::sum_int,
    CkCallback(CkReductionTarget(Transporter,comfinal), Disc()->Tr()) );
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

  BndEdgeMerger = CkReduction::addReducer(
                    tk::mergeHashMap< decltype(m_bndEdges)::key_type,
                                      decltype(m_bndEdges)::mapped_type > );
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
  for (const auto& eq : g_cgpde) eq.initialize( d->Coord(), m_u, d->T() );

  // Apply symmetry boundary conditions on initial conditions
  //for (const auto& eq : g_cgpde) eq.symbc( m_u, m_bnorm );

  // Output initial conditions to file (regardless of whether it was requested)
  writeFields( CkCallback(CkIndex_ALECG::init(), thisProxy[thisIndex]) );
}

//! [init and lhs]
void
ALECG::init()
// *****************************************************************************
// Initially compute left hand side diagonal matrix
// *****************************************************************************
{
  // Compute left-hand side of PDEs
  lhs();
}
//! [init and lhs]


void
ALECG::start()
// *****************************************************************************
//  Start time stepping
// *****************************************************************************
{
  // Start timer measuring time stepping wall clock time
  Disc()->Timer().zero();

  // Start time stepping by computing the size of the next time step)
  next();
}

//! [Compute own and send lhs on chare-boundary]
void
ALECG::lhs()
// *****************************************************************************
// Compute the left-hand side of transport equations
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

  // Continue after lhs is complete
  if (m_initial) start(); else lhs_complete();
}
//! [Merge lhs and continue]

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
    for (const auto& eq : g_cgpde) {
      auto eqdt = eq.dt( d->Coord(), d->Inpoel(), m_u );
      if (eqdt < mindt) mindt = eqdt;
    }

    // Scale smallest dt with CFL coefficient
    mindt *= g_inputdeck.get< tag::discr, tag::cfl >();
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
//! \param[in] newdt Size of this new time step
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
    eq.grad( d->Coord(), d->Inpoel(), m_u, m_grad );

  // Communicate gradients to other chares on chare-boundary
  if (d->NodeCommMap().empty())        // in serial we are done
    comgrad_complete();
  else // send gradients contributions to chare-boundary nodes to fellow chares
    for (const auto& [c,n] : d->NodeCommMap()) {
      std::vector< std::vector< tk::real > > g( n.size() );
      std::size_t j = 0;
      for (auto i : n) g[ j++ ] = m_grad[ tk::cref_find(d->Lid(),i) ];
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
  for (const auto& b : m_gradc) {
    auto lid = tk::cref_find( d->Lid(), b.first );
    for (ncomp_t c=0; c<m_grad.nprop(); ++c) m_grad(lid,c,0) += b.second[c];
  }

  // divide weak result in gradients by nodal volume
  for (std::size_t p=0; p<m_grad.nunk(); ++p)
    for (std::size_t c=0; c<m_grad.nprop(); ++c)
       m_grad(p,c,0) /= d->Vol()[p];

  // Compute own portion of right-hand side for all equations
  for (const auto& eq : g_cgpde)
    eq.rhs( d->T(), d->Coord(), d->Inpoel(), m_esued, m_psup, m_triinpoel,
            d->Gid(), m_dfnorm, m_grad, m_u, m_rhs );

  // Query and match user-specified boundary conditions to side sets
  m_bcdir = match( m_u.nprop(), d->T(), rkcoef[1][m_stage] * d->Dt(),
                   d->Coord(),  d->Lid(), m_bnode );

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

  for (std::size_t i=0; i<gid.size(); ++i) {
    m_rhsc[ gid[i] ] += R[i];
  }

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

  // Set Dirichlet BCs for lhs and rhs
  for (const auto& [b,bc] : m_bcdir) {
    for (ncomp_t c=0; c<ncomp; ++c) {
      if (bc[c].first) {
        m_lhs( b, c, 0 ) = 1.0;
        m_rhs( b, c, 0 ) = bc[c].second / (rkcoef[1][m_stage]*d->Dt());
      }
    }
  }

  // Update Un
  if (m_stage == 0) m_un = m_u;

  // Solve sytem
  m_u = rkcoef[0][m_stage] * m_un
    + rkcoef[1][m_stage] * (m_u + d->Dt() * m_rhs / m_lhs);

  // Apply symmetry BCs
  //for (const auto& eq : g_cgpde) eq.symbc( m_u, m_bnorm );

  //! [Continue after solve]
  if (m_stage < 2) {

    // Activate SDAG wait for next time step stage
    thisProxy[ thisIndex ].wait4grad();
    thisProxy[ thisIndex ].wait4rhs();

    // continue with next time step stage
    stage();

  } else {

    // Compute diagnostics, e.g., residuals
    auto diag_computed = m_diag.compute( *d, m_u );
    // Increase number of iterations and physical time
    d->next();
    // Continue to mesh refinement (if configured)
    if (!diag_computed) refine();
  }
  //! [Continue after solve]
}

void
ALECG::refine()
// *****************************************************************************
// Optionally refine/derefine mesh
// *****************************************************************************
{
  //! [Refine]
  auto d = Disc();

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
  const std::map< int, std::vector< std::size_t > >& /*bface*/,
  const std::map< int, std::vector< std::size_t > >& bnode,
  const std::vector< std::size_t >& /* triinpoel */ )
// *****************************************************************************
//  Receive new mesh from Refiner
//! \param[in] ginpoel Mesh connectivity with global node ids
//! \param[in] chunk New mesh chunk (connectivity and global<->local id maps)
//! \param[in] coord New mesh node coordinates
//! \param[in] addedNodes Newly added mesh nodes and their parents (local ids)
//! \param[in] addedTets Newly added mesh cells and their parents (local ids)
//! \param[in] nodeCommMap New node communication map
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
  d->resizePostAMR( chunk, coord, nodeCommMap );

  // Recompute derived data structures
  m_esued = tk::genEsued( d->Inpoel(), 4, tk::genEsup(d->Inpoel(),4) );
  m_psup = tk::genPsup( d->Inpoel(), 4, tk::genEsup( d->Inpoel(),4) );

  // Resize auxiliary solution vectors
  auto npoin = coord[0].size();
  auto nprop = m_u.nprop();
  m_u.resize( npoin, nprop );
  m_un.resize( npoin, nprop );
  m_lhs.resize( npoin, nprop );
  m_rhs.resize( npoin, nprop );

  // Update solution on new mesh
  for (const auto& n : addedNodes)
    for (std::size_t c=0; c<nprop; ++c)
      m_u(n.first,c,0) = (m_u(n.second[0],c,0) + m_u(n.second[1],c,0))/2.0;

  // Update physical-boundary node lists
  m_bnode = bnode;

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

  // Send mesh and fields data (solution dump) for output to file
  d->write( d->Inpoel(), d->Coord(), {}, tk::remap(m_bnode,d->Lid()), {}, {},
            nodefieldnames, {}, nodefields, c );
}

void
ALECG::out()
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
  // step
  if ( !((d->It()) % fieldfreq) ||
       (std::fabs(d->T()-term) < eps || d->It() >= nstep) )
    writeFields( CkCallback(CkIndex_ALECG::step(), thisProxy[thisIndex]) );
  else
    step();
}

void
ALECG::evalLB()
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
ALECG::evalRestart()
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

#include "NoWarning/alecg.def.h"
