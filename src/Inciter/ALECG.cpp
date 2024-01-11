// *****************************************************************************
/*!
  \file      src/Inciter/ALECG.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     ALECG for a PDE system with continuous Galerkin + ALE + RK
  \details   ALECG advances a system of partial differential equations (PDEs)
    using a continuous Galerkin (CG) finite element (FE) spatial discretization
    (using linear shapefunctions on tetrahedron elements) combined with a
    Runge-Kutta (RK) time stepping scheme in the arbitrary Eulerian-Lagrangian
    reference frame.
  \see The documentation in ALECG.hpp.
*/
// *****************************************************************************

#include "QuinoaBuildConfig.hpp"
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
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< CGPDE > g_cgpde;

//! Runge-Kutta coefficients
static const std::array< tk::real, 3 > rkcoef{{ 1.0/3.0, 1.0/2.0, 1.0 }};

} // inciter::

using inciter::ALECG;

ALECG::ALECG( const CProxy_Discretization& disc,
              const CProxy_Ghosts&,
              const std::map< int, std::vector< std::size_t > >& bface,
              const std::map< int, std::vector< std::size_t > >& bnode,
              const std::vector< std::size_t >& triinpoel ) :
  m_disc( disc ),
  m_nsol( 0 ),
  m_ngrad( 0 ),
  m_nrhs( 0 ),
  m_nbnorm( 0 ),
  m_ndfnorm( 0 ),
  m_nmblk( 0 ),
  m_bnode( bnode ),
  m_bface( bface ),
  m_triinpoel( tk::remap( triinpoel, Disc()->Lid() ) ),
  m_bndel( Disc()->bndel() ),
  m_dfnorm(),
  m_dfnormc(),
  m_dfn(),
  m_esup( tk::genEsup( Disc()->Inpoel(), 4 ) ),
  m_psup( tk::genPsup( Disc()->Inpoel(), 4, m_esup ) ),
  m_u( Disc()->Gid().size(),
       g_inputdeck.get< tag::component >().nprop( Disc()->MeshId() ) ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_rhsc(),
  m_chBndGrad( Disc()->Bid().size(), m_u.nprop()*3 ),
  m_dirbc(),
  m_chBndGradc(),
  m_diag(),
  m_bnorm(),
  m_bnormc(),
  m_symbcnodes(),
  m_farfieldbcnodes(),
  m_symbctri(),
  m_spongenodes(),
  m_timedepbcnodes(),
  m_timedepbcFn(),
  m_stage( 0 ),
  m_boxnodes(),
  m_edgenode(),
  m_edgeid(),
  m_dtp( m_u.nunk(), 0.0 ),
  m_tp( m_u.nunk(), g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_finished( 0 ),
  m_newmesh( 0 ),
  m_refinedmesh( 0 ),
  m_nusermeshblk( 0 ),
  m_nodeblockid(),
  m_nodeblockidc()
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

  auto d = Disc();

  // Perform optional operator-access-pattern mesh node reordering
  if (g_inputdeck.get< tag::discr, tag::operator_reorder >()) {

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

  // Query/update boundary-conditions-related data structures from user input
  queryBnd();

  // Activate SDAG wait for initially computing normals, and mesh blocks
  thisProxy[ thisIndex ].wait4norm();
  thisProxy[ thisIndex ].wait4meshblk();

  d->comfinal();

}
//! [Constructor]

void
ALECG::queryBnd()
// *****************************************************************************
// Query/update boundary-conditions-related data structures from user input
// *****************************************************************************
{
  auto d = Disc();

  // Query and match user-specified Dirichlet boundary conditions to side sets
  const auto steady = g_inputdeck.get< tag::discr, tag::steady_state >();
  if (steady) for (auto& deltat : m_dtp) deltat *= rkcoef[m_stage];
  m_dirbc = match( d->MeshId(), m_u.nprop(), d->T(), rkcoef[m_stage] * d->Dt(),
                   m_tp, m_dtp, d->Coord(), d->Lid(), m_bnode,
                   /* increment = */ false );
  if (steady) for (auto& deltat : m_dtp) deltat /= rkcoef[m_stage];

  // Prepare unique set of symmetry BC nodes
  auto sym = d->bcnodes< tag::bc, tag::bcsym >( m_bface, m_triinpoel );
  for (const auto& [s,nodes] : sym)
    m_symbcnodes.insert( begin(nodes), end(nodes) );

  // Prepare unique set of farfield BC nodes
  auto far = d->bcnodes< tag::bc, tag::bcfarfield >( m_bface, m_triinpoel );
  for (const auto& [s,nodes] : far)
    m_farfieldbcnodes.insert( begin(nodes), end(nodes) );

  // If farfield BC is set on a node, will not also set symmetry BC
  for (auto fn : m_farfieldbcnodes) m_symbcnodes.erase(fn);

  // Prepare boundary nodes contiguously accessible from a triangle-face loop
  m_symbctri.resize( m_triinpoel.size()/3, 0 );
  for (std::size_t e=0; e<m_triinpoel.size()/3; ++e)
    if (m_symbcnodes.find(m_triinpoel[e*3+0]) != end(m_symbcnodes))
      m_symbctri[e] = 1;

  // Prepare unique set of sponge nodes
  auto sponge = d->bcnodes< tag::sponge, tag::sideset >( m_bface, m_triinpoel );
  for (const auto& [s,nodes] : sponge)
    m_spongenodes.insert( begin(nodes), end(nodes) );

  // Prepare unique set of time dependent BC nodes
  m_timedepbcnodes.clear();
  m_timedepbcFn.clear();
  const auto& timedep =
    g_inputdeck.get< tag::param, tag::compflow, tag::bctimedep >();
  if (!timedep.empty()) {
    m_timedepbcnodes.resize(timedep.size());
    m_timedepbcFn.resize(timedep.size());
    std::size_t ib=0;
    for (const auto& bndry : timedep) {
      std::unordered_set< std::size_t > nodes;
      for (const auto& s : bndry.template get< tag::sideset >()) {
        auto k = m_bnode.find( std::stoi(s) );
        if (k != end(m_bnode)) {
          for (auto g : k->second) {      // global node ids on side set
            nodes.insert( tk::cref_find(d->Lid(),g) );
          }
        }
      }
      m_timedepbcnodes[ib].insert( begin(nodes), end(nodes) );

      // Store user defined discrete function in time. This is done in the same
      // loop as the BC nodes, so that the indices for the two vectors
      // m_timedepbcnodes and m_timedepbcFn are consistent with each other
      auto fn = bndry.template get< tag::fn >();
      for (std::size_t ir=0; ir<fn.size()/6; ++ir) {
        m_timedepbcFn[ib].push_back({{ fn[ir*6+0], fn[ir*6+1], fn[ir*6+2],
          fn[ir*6+3], fn[ir*6+4], fn[ir*6+5] }});
      }
      ++ib;
    }
  }

  Assert(m_timedepbcFn.size() == m_timedepbcnodes.size(), "Incorrect number of "
    "time dependent functions.");

  // Query ALE mesh velocity boundary condition node lists and node lists at
  // which ALE moves boundaries
  d->meshvelBnd( m_bface, m_bnode, m_triinpoel );
}

void
ALECG::norm()
// *****************************************************************************
// Start (re-)computing boundary point-, and dual-face normals
// *****************************************************************************
{
  auto d = Disc();

  // Query nodes at which symmetry BCs are specified
  auto bn = d->bcnodes< tag::bc, tag::bcsym >( m_bface, m_triinpoel );

  // Query nodes at which farfield BCs are specified
  auto far = d->bcnodes< tag::bc, tag::bcfarfield >( m_bface, m_triinpoel );
  // Merge BC data where boundary-point normals are required
  for (const auto& [s,n] : far) bn[s].insert( begin(n), end(n) );

  // Query nodes at which mesh velocity symmetry BCs are specified
  std::unordered_map<int, std::unordered_set< std::size_t >> ms;
  for (const auto& s : g_inputdeck.get< tag::ale, tag::bcsym >()) {
    auto k = m_bface.find( std::stoi(s) );
    if (k != end(m_bface)) {
      auto& n = ms[ k->first ];
      for (auto f : k->second) {
        n.insert( m_triinpoel[f*3+0] );
        n.insert( m_triinpoel[f*3+1] );
        n.insert( m_triinpoel[f*3+2] );
      }
    }
  }
  // Merge BC data where boundary-point normals are required
  for (const auto& [s,n] : ms) bn[s].insert( begin(n), end(n) );

  // Compute boundary point normals
  bnorm( bn );

  // Compute dual-face normals associated to edges
  dfnorm();
}

std::array< tk::real, 3 >
ALECG::edfnorm( const tk::UnsMesh::Edge& edge,
                const std::unordered_map< tk::UnsMesh::Edge,
                        std::vector< std::size_t >,
                        tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> >& esued )
const
// *****************************************************************************
//  Compute normal of dual-mesh associated to edge
//! \param[in] edge Edge whose dual-face normal to compute given by local ids
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
    // The constant 1/48: Eq (12) from Waltz et al. Computers & fluids (92) 2014
    // The result of the integral of shape function N on a tet is V/4.
    // This can be written as J/(6*4). Eq (12) has a 1/2 multiplier.
    // This leads to J/48.
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
ALECG::bnorm( const std::unordered_map< int,
                std::unordered_set< std::size_t > >& bcnodes )
// *****************************************************************************
//  Compute boundary point normals
//! \param[in] bcnodes Local node ids associated to side set ids at which BCs
//!    are set that require normals
//*****************************************************************************
{
  auto d = Disc();

  m_bnorm = cg::bnorm( m_bface, m_triinpoel, d->Coord(), d->Gid(), bcnodes );

  // Send our nodal normal contributions to neighbor chares
  if (d->NodeCommMap().empty())
    comnorm_complete();
  else
    for (const auto& [ neighborchare, sharednodes ] : d->NodeCommMap()) {
      std::unordered_map< int,
        std::unordered_map< std::size_t, std::array< tk::real, 4 > > > exp;
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
ALECG::comnorm( const std::unordered_map< int,
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

//! [setup]
void
ALECG::setup()
// *****************************************************************************
// Start setup for solution
// *****************************************************************************
{
  auto d = Disc();

  // Determine nodes inside user-defined IC box
  g_cgpde[d->MeshId()].IcBoxNodes( d->Coord(), d->Inpoel(),
    d->ElemBlockId(), m_boxnodes, m_nodeblockid, m_nusermeshblk );

  // Communicate mesh block nodes to other chares on chare-boundary
  if (d->NodeCommMap().empty())        // in serial we are done
    comblk_complete();
  else // send mesh block information to chare-boundary nodes to fellow chares
    for (const auto& [c,n] : d->NodeCommMap()) {
      // data structure assigning block ids (set of values) to nodes (index).
      // although nodeblockid is a map with key-blockid and value-nodeid, the
      // sending data structure has to be inverted, because of how communication
      // data is handled.
      std::vector< std::set< std::size_t > > mb( n.size() );
      std::size_t j = 0;
      for (auto i : n) {
        for (const auto& [blid, ndset] : m_nodeblockid) {
          // if node was found in a block, add to send-data
          if (ndset.find(tk::cref_find(d->Lid(),i)) != ndset.end())
            mb[j].insert(blid);
        }
        if (m_nusermeshblk > 0)
          Assert(mb[j].size() > 0, "Sending no block data for node");
        ++j;
      }
      thisProxy[c].comblk( std::vector<std::size_t>(begin(n),end(n)), mb );
    }

  ownblk_complete();
}

void
ALECG::comblk( const std::vector< std::size_t >& gid,
               const std::vector< std::set< std::size_t > >& mb )
// *****************************************************************************
//  Receive mesh block information for nodes on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive RHS contributions
//! \param[in] mb Block ids for each node on chare-boundaries
//! \details This function receives mesh block information for nodes on chare
//!   boundaries. While m_nodeblockid stores block information for own nodes,
//!   m_nodeblockidc collects the neighbor chare information during
//!   communication. This way work on m_nodeblockid and m_nodeblockidc is
//!   overlapped. The two are combined in continueSetup().
// *****************************************************************************
{
  Assert( mb.size() == gid.size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i) {
    for (const auto& blid : mb[i]) {
      m_nodeblockidc[blid].insert(gid[i]);
    }
  }

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_nmblk == Disc()->NodeCommMap().size()) {
    m_nmblk = 0;
    comblk_complete();
  }
}

void
ALECG::continueSetup()
// *****************************************************************************
// Continue setup for solution, after communication for mesh blocks
// *****************************************************************************
{
  auto d = Disc();

  // Combine own and communicated mesh block information
  for (const auto& [blid, ndset] : m_nodeblockidc) {
    for (const auto& i : ndset) {
      auto lid = tk::cref_find(d->Lid(), i);
      m_nodeblockid[blid].insert(lid);
    }
  }

  // clear receive buffer
  tk::destroy(m_nodeblockidc);

  // Compute volume of user-defined box IC
  d->boxvol( m_boxnodes, m_nodeblockid, m_nusermeshblk );

  // Query time history field output labels from all PDEs integrated
  const auto& hist_points = g_inputdeck.get< tag::history, tag::point >();
  if (!hist_points.empty()) {
    std::vector< std::string > histnames;
    auto n = g_cgpde[d->MeshId()].histNames();
    histnames.insert( end(histnames), begin(n), end(n) );
    d->histheader( std::move(histnames) );
  }
}
//! [setup]

void
ALECG::volumetric( tk::Fields& u, const std::vector< tk::real >& v )
// *****************************************************************************
//  Multiply solution with mesh volume
//! \param[in,out] u Solution vector
//! \param[in] v Volume to multiply with
// *****************************************************************************
{
  Assert( v.size() == u.nunk(), "Size mismatch" );

  for (std::size_t i=0; i<u.nunk(); ++i)
    for (ncomp_t c=0; c<u.nprop(); ++c)
      u(i,c) *= v[i];
}

void
ALECG::conserved( tk::Fields& u, const std::vector< tk::real >& v )
// *****************************************************************************
//  Divide solution with mesh volume
//! \param[in,out] u Solution vector
//! \param[in] v Volume to divide with
// *****************************************************************************
{
  Assert( v.size() == u.nunk(), "Size mismatch" );

  for (std::size_t i=0; i<u.nunk(); ++i)
    for (ncomp_t c=0; c<u.nprop(); ++c) {
      u(i,c) /= v[i];
    }
}

void
ALECG::box( tk::real v, const std::vector< tk::real >& blkvols )
// *****************************************************************************
// Receive total box IC volume and set conditions in box
//! \param[in] v Total volume within user-specified box
//! \param[in] blkvols Vector of mesh block discrete volumes with user ICs
// *****************************************************************************
{
  Assert(blkvols.size() == m_nusermeshblk,
    "Incorrect size of block volume vector");
  auto d = Disc();

  // Store user-defined box/block IC volume
  d->Boxvol() = v;
  d->MeshBlkVol() = blkvols;

  // Set initial conditions for all PDEs
  g_cgpde[d->MeshId()].initialize( d->Coord(), m_u, d->T(), d->Boxvol(),
    m_boxnodes, d->MeshBlkVol(), m_nodeblockid );

  // Multiply conserved variables with mesh volume
  volumetric( m_u, Disc()->Vol() );

  // Initialize nodal mesh volumes at previous time step stage
  d->Voln() = d->Vol();

  // Start computing the mesh mesh velocity for ALE
  meshvelstart();
}

void
ALECG::meshvelstart()
// *****************************************************************************
// Start computing the mesh mesh velocity for ALE
// *****************************************************************************
{
  auto d = Disc();

  // Apply boundary conditions on numerical solution
  BC();

  conserved( m_u, d->Vol() );

  // query fluid velocity across all systems integrated
  tk::UnsMesh::Coords vel;
  g_cgpde[d->MeshId()].velocity( m_u, vel );
  // query speed of sound in mesh nodes across all systems integrated
  std::vector< tk::real > soundspeed;
  g_cgpde[d->MeshId()].soundspeed( m_u, soundspeed );

  volumetric( m_u, d->Vol() );

  // Start computing the mesh mesh velocity for ALE
  d->meshvelStart( vel, soundspeed, m_bnorm, rkcoef[m_stage] * d->Dt(),
    CkCallback(CkIndex_ALECG::meshveldone(), thisProxy[thisIndex]) );
}

void
ALECG::meshveldone()
// *****************************************************************************
// Done with computing the mesh velocity for ALE
// *****************************************************************************
{
  // Assess and record mesh velocity linear solver conergence
  Disc()->meshvelConv();

  // Continue
  if (Disc()->Initial()) {

    conserved( m_u, Disc()->Vol() );

    // Initiate IC transfer (if coupled)
    Disc()->transfer( m_u, 0,
      CkCallback(CkIndex_ALECG::transfer_complete(), thisProxy[thisIndex]) );

    lhs();

  } else {

    ale();

  }
}

//! [start]
void
ALECG::start()
// *****************************************************************************
// Start time stepping
// *****************************************************************************
{
  // Set flag that indicates that we are now during time stepping
  Disc()->Initial( 0 );
  // Start timer measuring time stepping wall clock time
  Disc()->Timer().zero();
  // Zero grind-timer
  Disc()->grindZero();
  // Continue to first time step
  next();
}
//! [start]

//! [Compute lhs]
void
ALECG::lhs()
// *****************************************************************************
// Compute the left-hand side of transport equations
//! \details Also (re-)compute all data structures if the mesh changed.
// *****************************************************************************
{
  // No need for LHS in ALECG

  // (Re-)compute boundary point-, and dual-face normals
  norm();
}
//! [Compute lhs]

//! [Merge normals and continue]
void
ALECG::mergelhs()
// *****************************************************************************
// The own and communication portion of the left-hand side is complete
// *****************************************************************************
{
  // Combine own and communicated contributions of normals
  normfinal();

  if (Disc()->Initial()) {
    volumetric( m_u, Disc()->Vol() );
    // Output initial conditions to file
    writeFields( CkCallback(CkIndex_ALECG::start(), thisProxy[thisIndex]) );
  } else {
    norm_complete();
  }
}
//! [Merge normals and continue]

void
ALECG::normfinal()
// *****************************************************************************
//  Finish computing dual-face and boundary point normals
// *****************************************************************************
{
  auto d = Disc();
  const auto& lid = d->Lid();

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

  // Generate list of unique edges
  tk::UnsMesh::EdgeSet uedge;
  for (std::size_t p=0; p<m_u.nunk(); ++p)
    for (auto q : tk::Around(m_psup,p))
      uedge.insert( {p,q} );

  // Flatten edge list
  m_edgenode.resize( uedge.size() * 2 );
  std::size_t f = 0;
  const auto& gid = d->Gid();
  for (auto&& [p,q] : uedge) {
    if (gid[p] > gid[q]) {
      m_edgenode[f+0] = std::move(q);
      m_edgenode[f+1] = std::move(p);
    } else {
      m_edgenode[f+0] = std::move(p);
      m_edgenode[f+1] = std::move(q);
    }
    f += 2;
  }
  tk::destroy(uedge);

  // Convert dual-face normals to streamable (and vectorizable) data structure
  m_dfn.resize( m_edgenode.size() * 3 );      // 2 vectors per access
  std::unordered_map< tk::UnsMesh::Edge, std::size_t,
                      tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> > eid;
  for (std::size_t e=0; e<m_edgenode.size()/2; ++e) {
    auto p = m_edgenode[e*2+0];
    auto q = m_edgenode[e*2+1];
    eid[{p,q}] = e;
    std::array< std::size_t, 2 > g{ gid[p], gid[q] };
    auto n = tk::cref_find( m_dfnorm, g );
    // figure out if this is an edge on the parallel boundary
    auto nit = m_dfnormc.find( g );
    auto m = ( nit != m_dfnormc.end() ) ? nit->second : n;
    m_dfn[e*6+0] = n[0];
    m_dfn[e*6+1] = n[1];
    m_dfn[e*6+2] = n[2];
    m_dfn[e*6+3] = m[0];
    m_dfn[e*6+4] = m[1];
    m_dfn[e*6+5] = m[2];
  }

  tk::destroy( m_dfnorm );
  tk::destroy( m_dfnormc );

  // Flatten edge id data structure
  m_edgeid.resize( m_psup.first.size() );
  for (std::size_t p=0,k=0; p<m_u.nunk(); ++p)
    for (auto q : tk::Around(m_psup,p))
      m_edgeid[k++] = tk::cref_find( eid, {p,q} );
}

void
ALECG::BC()
// *****************************************************************************
// Apply boundary conditions
// \details The following BC enforcement changes the initial condition or
//!   updated solution (dependending on when it is called) to ensure strong
//!   imposition of the BCs. This is a matter of choice. Another alternative is
//!   to only apply BCs when computing fluxes at boundary faces, thereby only
//!   weakly enforcing the BCs. The former is conventionally used in continunous
//!   Galerkin finite element methods (such as ALECG implements), whereas the
//!   latter, in finite volume methods.
// *****************************************************************************
{
  auto d = Disc();
  const auto& coord = d->Coord();

  conserved( m_u, d->Vol() );

  // Apply Dirichlet BCs
  for (const auto& [b,bc] : m_dirbc)
    for (ncomp_t c=0; c<m_u.nprop(); ++c)
      if (bc[c].first) m_u(b,c) = bc[c].second;

  // Apply symmetry BCs
  g_cgpde[d->MeshId()].symbc( m_u, coord, m_bnorm, m_symbcnodes );

  // Apply farfield BCs
  g_cgpde[d->MeshId()].farfieldbc( m_u, coord, m_bnorm, m_farfieldbcnodes );

  // Apply sponge conditions
  g_cgpde[d->MeshId()].sponge( m_u, coord, m_spongenodes );

  // Apply user defined time dependent BCs
  g_cgpde[d->MeshId()].timedepbc( d->T(), m_u, m_timedepbcnodes,
    m_timedepbcFn );

  volumetric( m_u, d->Vol() );
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
    conserved( m_u, Disc()->Vol() );
    if (g_inputdeck.get< tag::discr, tag::steady_state >()) {

      // compute new dt for each mesh point
      g_cgpde[d->MeshId()].dt( d->It(), d->Vol(), m_u, m_dtp );

      // find the smallest dt of all nodes on this chare
      mindt = *std::min_element( begin(m_dtp), end(m_dtp) );

    } else {    // compute new dt for this chare

      // find the smallest dt of all equations on this chare
      auto eqdt = g_cgpde[d->MeshId()].dt( d->Coord(), d->Inpoel(), d->T(),
        d->Dtn(), m_u, d->Vol(), d->Voln() );
      if (eqdt < mindt) mindt = eqdt;

    }
    volumetric( m_u, Disc()->Vol() );
    //! [Find the minimum dt across all PDEs integrated]

  }

  //! [Advance]
  // Actiavate SDAG waits for next time step stage
  thisProxy[ thisIndex ].wait4grad();
  thisProxy[ thisIndex ].wait4rhs();

  // Contribute to minimum dt across all chares and advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(ALECG,advance), thisProxy) );
  //! [Advance]
}

void
ALECG::advance( tk::real newdt, tk::real )
// *****************************************************************************
// Advance equations to next time step
//! \param[in] newdt The smallest dt across the whole problem
// *****************************************************************************
{
  auto d = Disc();

  // Set new time step size
  if (m_stage == 0) d->setdt( newdt );

  // Compute gradients for next time step
  chBndGrad();
}

void
ALECG::chBndGrad()
// *****************************************************************************
// Compute nodal gradients at chare-boundary nodes. Gradients at internal nodes
// are calculated locally as needed and are not stored.
// *****************************************************************************
{
  auto d = Disc();

  // Divide solution with mesh volume
  conserved( m_u, Disc()->Vol() );
  // Compute own portion of gradients for all equations
  g_cgpde[d->MeshId()].chBndGrad( d->Coord(), d->Inpoel(), m_bndel, d->Gid(),
    d->Bid(), m_u, m_chBndGrad );
  // Multiply solution with mesh volume
  volumetric( m_u, Disc()->Vol() );

  // Communicate gradients to other chares on chare-boundary
  if (d->NodeCommMap().empty())        // in serial we are done
    comgrad_complete();
  else // send gradients contributions to chare-boundary nodes to fellow chares
    for (const auto& [c,n] : d->NodeCommMap()) {
      std::vector< std::vector< tk::real > > g( n.size() );
      std::size_t j = 0;
      for (auto i : n) g[ j++ ] = m_chBndGrad[ tk::cref_find(d->Bid(),i) ];
      thisProxy[c].comChBndGrad( std::vector<std::size_t>(begin(n),end(n)), g );
    }

  owngrad_complete();
}

void
ALECG::comChBndGrad( const std::vector< std::size_t >& gid,
                     const std::vector< std::vector< tk::real > >& G )
// *****************************************************************************
//  Receive contributions to nodal gradients on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive grad contributions
//! \param[in] G Partial contributions of gradients to chare-boundary nodes
//! \details This function receives contributions to m_chBndGrad, which stores
//!   nodal gradients at mesh chare-boundary nodes. While m_chBndGrad stores
//!   own contributions, m_chBndGradc collects the neighbor chare
//!   contributions during communication. This way work on m_chBndGrad and
//!   m_chBndGradc is overlapped. The two are combined in rhs().
// *****************************************************************************
{
  Assert( G.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  for (std::size_t i=0; i<gid.size(); ++i) m_chBndGradc[ gid[i] ] += G[i];

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
  for (const auto& [gid,g] : m_chBndGradc) {
    auto bid = tk::cref_find( d->Bid(), gid );
    for (ncomp_t c=0; c<m_chBndGrad.nprop(); ++c)
      m_chBndGrad(bid,c) += g[c];
  }

  // clear gradients receive buffer
  tk::destroy(m_chBndGradc);

  const auto steady = g_inputdeck.get< tag::discr, tag::steady_state >();

  // Compute own portion of right-hand side for all equations
  auto prev_rkcoef = m_stage == 0 ? 0.0 : rkcoef[m_stage-1];
  if (steady)
    for (std::size_t p=0; p<m_tp.size(); ++p) m_tp[p] += prev_rkcoef * m_dtp[p];
  conserved( m_u, Disc()->Vol() );
  g_cgpde[d->MeshId()].rhs( d->T() + prev_rkcoef * d->Dt(), d->Coord(), d->Inpoel(),
          m_triinpoel, d->Gid(), d->Bid(), d->Lid(), m_dfn, m_psup, m_esup,
          m_symbctri, m_spongenodes, d->Vol(), m_edgenode, m_edgeid,
          m_boxnodes, m_chBndGrad, m_u, d->meshvel(), m_tp, d->Boxvol(),
          m_rhs );
  volumetric( m_u, Disc()->Vol() );
  if (steady)
    for (std::size_t p=0; p<m_tp.size(); ++p) m_tp[p] -= prev_rkcoef * m_dtp[p];

  // Query/update boundary-conditions-related data structures from user input
  queryBnd();

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
//  Advance systems of equations
// *****************************************************************************
{
  auto d = Disc();

  // Combine own and communicated contributions to rhs
  for (const auto& b : m_rhsc) {
    auto lid = tk::cref_find( d->Lid(), b.first );
    for (ncomp_t c=0; c<m_rhs.nprop(); ++c) m_rhs(lid,c) += b.second[c];
  }

  // clear receive buffer
  tk::destroy(m_rhsc);

  // Update state at time n
  if (m_stage == 0) {
    m_un = m_u;
    if (g_inputdeck.get< tag::ale, tag::ale >()) d->UpdateCoordn();
  }

  // Solve the sytem
  if (g_inputdeck.get< tag::discr, tag::steady_state >()) {

    // Advance solution, converging to steady state
    for (std::size_t i=0; i<m_u.nunk(); ++i)
      for (ncomp_t c=0; c<m_u.nprop(); ++c)
        m_u(i,c) = m_un(i,c) + rkcoef[m_stage] * m_dtp[i] * m_rhs(i,c);

  } else {

    auto adt = rkcoef[m_stage] * d->Dt();

    // Advance unsteady solution
    m_u = m_un + adt * m_rhs;

    // Advance mesh if ALE is enabled
    if (g_inputdeck.get< tag::ale, tag::ale >()) {
      auto& coord = d->Coord();
      const auto& w = d->meshvel();
      for (auto j : g_inputdeck.get< tag::ale, tag::mesh_motion >())
        for (std::size_t i=0; i<coord[j].size(); ++i)
          coord[j][i] = d->Coordn()[j][i] + adt * w(i,j);
    }

  }

  m_newmesh = 0;  // recompute normals after ALE (if enabled)
  m_refinedmesh = 0;  // mesh has not been refined by ALE
  // Activate SDAG waits
  thisProxy[ thisIndex ].wait4norm();
  thisProxy[ thisIndex ].wait4mesh();

  //! [Continue after solve]
  // Recompute mesh volumes if ALE is enabled
  if (g_inputdeck.get< tag::ale, tag::ale >()) {

    transfer_complete();
    // Save nodal volumes at previous time step stage
    d->Voln() = d->Vol();
    // Prepare for recomputing the nodal volumes
    d->startvol();
    auto meshid = d->MeshId();
    contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
                CkCallback(CkReductionTarget(Transporter,resized), d->Tr()) );

  } else {

    norm_complete();
    resized();

  }
  //! [Continue after solve]
}

void
ALECG::ale()
// *****************************************************************************
//  Continue after computing the new mesh velocity for ALE
// *****************************************************************************
{
  auto d = Disc();

  if (m_stage < 2) {

    // Activate SDAG wait for next time step stage
    thisProxy[ thisIndex ].wait4grad();
    thisProxy[ thisIndex ].wait4rhs();

    // continue to mesh-to-mesh transfer (if coupled)
    transfer();

  } else {

    // Ensure new field output file if mesh moved if ALE is enabled
    if (g_inputdeck.get< tag::ale, tag::ale >()) {
      d->Itf() = 0;  // Zero field output iteration count if mesh moved
      ++d->Itr();    // Increase number of iterations with a change in the mesh
    }

    // Compute diagnostics, e.g., residuals
    conserved( m_u, Disc()->Vol() );
    conserved( m_un, Disc()->Voln() );
    auto diag_computed = m_diag.compute( *d, m_u, m_un, m_bnorm,
                                         m_symbcnodes, m_farfieldbcnodes );
    volumetric( m_u, Disc()->Vol() );
    volumetric( m_un, Disc()->Voln() );
    // Increase number of iterations and physical time
    d->next();
    // Advance physical time for local time stepping
    if (g_inputdeck.get< tag::discr, tag::steady_state >())
      for (std::size_t i=0; i<m_u.nunk(); ++i) m_tp[i] += m_dtp[i];
    // Continue to mesh refinement (if configured)
    if (!diag_computed) refine( std::vector< tk::real >( m_u.nprop(), 1.0 ) );

  }
}

//! [Refine]
void
ALECG::refine( const std::vector< tk::real >& l2res )
// *****************************************************************************
// Optionally refine/derefine mesh
//! \param[in] l2res L2-norms of the residual for each scalar component
//!   computed across the whole problem
// *****************************************************************************
{
  auto d = Disc();

  const auto steady = g_inputdeck.get< tag::discr, tag::steady_state >();
  const auto residual = g_inputdeck.get< tag::discr, tag::residual >();
  const auto rc = g_inputdeck.get< tag::discr, tag::rescomp >() - 1;

  if (steady) {

    // this is the last time step if max time of max number of time steps
    // reached or the residual has reached its convergence criterion
    if (d->finished() or l2res[rc] < residual) m_finished = 1;

  } else {

    // this is the last time step if max time or max iterations reached
    if (d->finished()) m_finished = 1;

  }

  auto dtref = g_inputdeck.get< tag::amr, tag::dtref >();
  auto dtfreq = g_inputdeck.get< tag::amr, tag::dtfreq >();

  // Activate SDAG waits for re-computing the normals
  m_newmesh = 1;  // recompute normals after AMR (if enabled)
  thisProxy[ thisIndex ].wait4norm();
  thisProxy[ thisIndex ].wait4mesh();

  // if t>0 refinement enabled and we hit the frequency
  if (dtref && !(d->It() % dtfreq)) {   // refine

    // Convert to conserved unknowns, since the next step changes volumes
    conserved(m_u, d->Vol());

    m_refinedmesh = 1;
    d->startvol();
    d->Ref()->dtref( m_bface, m_bnode, m_triinpoel );
    d->refined() = 1;

  } else {      // do not refine

    m_refinedmesh = 0;
    d->refined() = 0;
    norm_complete();
    resized();

  }
}
//! [Refine]

//! [Resize]
void
ALECG::resizePostAMR(
  const std::vector< std::size_t >& /*ginpoel*/,
  const tk::UnsMesh::Chunk& chunk,
  const tk::UnsMesh::Coords& coord,
  const std::unordered_map< std::size_t, tk::UnsMesh::Edge >& addedNodes,
  const std::unordered_map< std::size_t, std::size_t >& /*addedTets*/,
  const std::set< std::size_t >& removedNodes,
  const std::unordered_map< std::size_t, std::size_t >& amrNodeMap,
  const tk::NodeCommMap& nodeCommMap,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& bnode,
  const std::vector< std::size_t >& triinpoel,
  const std::unordered_map< std::size_t, std::set< std::size_t > >& elemblockid )
// *****************************************************************************
//  Receive new mesh from Refiner
//! \param[in] ginpoel Mesh connectivity with global node ids
//! \param[in] chunk New mesh chunk (connectivity and global<->local id maps)
//! \param[in] coord New mesh node coordinates
//! \param[in] addedNodes Newly added mesh nodes and their parents (local ids)
//! \param[in] addedTets Newly added mesh cells and their parents (local ids)
//! \param[in] removedNodes Newly removed mesh node local ids
//! \param[in] amrNodeMap Node id map after amr (local ids)
//! \param[in] nodeCommMap New node communication map
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] bnode Boundary-node lists mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
//! \param[in] elemblockid Local tet ids associated with mesh block ids
// *****************************************************************************
{
  auto d = Disc();

  d->Itf() = 0;  // Zero field output iteration count if AMR
  ++d->Itr();    // Increase number of iterations with a change in the mesh

  // Resize mesh data structures after mesh refinement
  d->resizePostAMR( chunk, coord, amrNodeMap, nodeCommMap, removedNodes,
    elemblockid );

  // Remove newly removed nodes from solution vectors
  m_u.rm(removedNodes);
  m_un.rm(removedNodes);
  m_rhs.rm(removedNodes);

  // Resize auxiliary solution vectors
  auto npoin = coord[0].size();
  auto nprop = m_u.nprop();
  m_u.resize( npoin );
  m_un.resize( npoin );
  m_rhs.resize( npoin );
  m_chBndGrad.resize( d->Bid().size() );
  tk::destroy(m_esup);
  tk::destroy(m_psup);
  m_esup = tk::genEsup( d->Inpoel(), 4 );
  m_psup = tk::genPsup( d->Inpoel(), 4, m_esup );

  // Update solution on new mesh
  for (const auto& n : addedNodes)
    for (std::size_t c=0; c<nprop; ++c) {
      Assert(n.first < m_u.nunk(), "Added node index out of bounds post-AMR");
      Assert(n.second[0] < m_u.nunk() && n.second[1] < m_u.nunk(),
        "Indices of parent-edge nodes out of bounds post-AMR");
      m_u(n.first,c) = (m_u(n.second[0],c) + m_u(n.second[1],c))/2.0;
    }

  // Update physical-boundary node-, face-, and element lists
  m_bnode = bnode;
  m_bface = bface;
  m_triinpoel = tk::remap( triinpoel, d->Lid() );

  auto meshid = d->MeshId();
  contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
              CkCallback(CkReductionTarget(Transporter,resized), d->Tr()) );
}
//! [Resize]

void
ALECG::resized()
// *****************************************************************************
// Resizing data sutrctures after mesh refinement has been completed
// *****************************************************************************
{
  auto d = Disc();

  // Revert to volumetric unknowns, if soln was converted in ALECG::refine()
  auto dtref = g_inputdeck.get< tag::amr, tag::dtref >();
  auto dtfreq = g_inputdeck.get< tag::amr, tag::dtfreq >();
  if (dtref && !(d->It() % dtfreq) && m_refinedmesh==1) {
    volumetric(m_u, d->Vol());
    // Update previous volumes after refinement
    d->Voln() = d->Vol();
  }

  resize_complete();
}

void
ALECG::transfer()
// *****************************************************************************
// Transfer solution to other solver and mesh if coupled
// *****************************************************************************
{
  // Initiate solution transfer (if coupled)

//TODO: enable this for during-timestepping solution transfer
//  Disc()->transfer(m_u, CkCallback(CkIndex_ALECG::stage(), thisProxy[thisIndex]));
  thisProxy[thisIndex].stage();
}

//! [stage]
void
ALECG::stage()
// *****************************************************************************
// Evaluate whether to continue with next time step stage
// *****************************************************************************
{
  transfer_complete();

  // Increment Runge-Kutta stage counter
  ++m_stage;

  // if not all Runge-Kutta stages complete, continue to next time stage,
  // otherwise output field data to file(s)
  if (m_stage < 3) chBndGrad(); else out();
}
//! [stage]

void
ALECG::writeFields( CkCallback c )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  if (g_inputdeck.get< tag::cmd, tag::benchmark >()) {

    c.send();

  } else {

    auto d = Disc();
    const auto& coord = d->Coord();

    // if coupled: depvars: src:'a', dst:'b','c',...
    char depvar = 0;
    if (not d->Transfers().empty()) {
      depvar = 'a' + static_cast< char >( d->MeshId() );
    }

    // Query fields names requested by user
    auto nodefieldnames = numericFieldNames( tk::Centering::NODE, depvar );

    // Collect field output from numerical solution requested by user
    conserved( m_u, Disc()->Vol() );
    auto nodefields = numericFieldOutput( m_u, tk::Centering::NODE, m_u, depvar );
    volumetric( m_u, Disc()->Vol() );

    //! Lambda to put in a field for output if not empty
    auto add_node_field = [&]( const auto& name, const auto& field ){
      if (not field.empty()) {
        nodefieldnames.push_back( name );
        nodefields.push_back( field );
      }
    };

    // Output mesh velocity if ALE is enabled
    if (g_inputdeck.get< tag::ale, tag::ale >()) {
      const auto& w = d->meshvel();
      add_node_field( "x-mesh-velocity", w.extract_comp(0) );
      add_node_field( "y-mesh-velocity", w.extract_comp(1) );
      add_node_field( "z-mesh-velocity", w.extract_comp(2) );
      add_node_field( "volume", d->Vol() );
    }

    // Collect field output names for analytical solutions
    analyticFieldNames( g_cgpde[d->MeshId()], tk::Centering::NODE,
      nodefieldnames );

    // Collect field output from analytical solutions (if exist)
    analyticFieldOutput( g_cgpde[d->MeshId()], tk::Centering::NODE, coord[0],
      coord[1], coord[2], d->T(), nodefields );

    // Query and collect nodal block and surface field names from PDEs integrated
    std::vector< std::string > nodesurfnames;
    auto sn = g_cgpde[d->MeshId()].surfNames();
    nodesurfnames.insert( end(nodesurfnames), begin(sn), end(sn) );

    // Collect nodal block and surface field solution
    std::vector< std::vector< tk::real > > nodesurfs;
    conserved( m_u, Disc()->Vol() );
    auto so = g_cgpde[d->MeshId()].surfOutput( tk::bfacenodes(m_bface,
      m_triinpoel), m_u );
    nodesurfs.insert( end(nodesurfs), begin(so), end(so) );

    // Collect elemental block and surface field names from PDEs integrated
    auto elemsurfnames = nodesurfnames;

    // Collect elemental block and surface field solution
    std::vector< std::vector< tk::real > > elemsurfs;
    auto eso = g_cgpde[d->MeshId()].elemSurfOutput( m_bface, m_triinpoel, m_u );
    elemsurfs.insert( end(elemsurfs), begin(eso), end(eso) );

    // Query refinement data
    auto dtref = g_inputdeck.get< tag::amr, tag::dtref >();

    std::tuple< std::vector< std::string >,
                std::vector< std::vector< tk::real > >,
                std::vector< std::string >,
                std::vector< std::vector< tk::real > > > r;
    if (dtref) r = d->Ref()->refinementFields();
    volumetric( m_u, Disc()->Vol() );

    auto& refinement_elemfieldnames = std::get< 0 >( r );
    auto& refinement_elemfields = std::get< 1 >( r );
    auto& refinement_nodefieldnames = std::get< 2 >( r );
    auto& refinement_nodefields = std::get< 3 >( r );

    nodefieldnames.insert( end(nodefieldnames),
      begin(refinement_nodefieldnames), end(refinement_nodefieldnames) );
    nodefields.insert( end(nodefields),
      begin(refinement_nodefields), end(refinement_nodefields) );

    auto elemfieldnames = std::move(refinement_elemfieldnames);
    auto elemfields = std::move(refinement_elemfields);

    Assert( elemfieldnames.size() == elemfields.size(), "Size mismatch" );
    Assert( nodefieldnames.size() == nodefields.size(), "Size mismatch" );

    // Send mesh and fields data (solution dump) for output to file
    d->write( d->Inpoel(), coord, m_bface, tk::remap(m_bnode,d->Lid()),
              m_triinpoel, elemfieldnames, nodefieldnames, elemsurfnames,
              nodesurfnames, elemfields, nodefields, elemsurfs, nodesurfs, c );

  }
}

void
ALECG::out()
// *****************************************************************************
// Output mesh field data
// *****************************************************************************
{
  auto d = Disc();

  // Output time history
  if (d->histiter() or d->histtime() or d->histrange()) {
    std::vector< std::vector< tk::real > > hist;
    conserved( m_u, Disc()->Vol() );
    auto h = g_cgpde[d->MeshId()].histOutput( d->Hist(), d->Inpoel(), m_u );
    hist.insert( end(hist), begin(h), end(h) );
    volumetric( m_u, Disc()->Vol() );
    d->history( std::move(hist) );
  }

  // Output field data
  if (d->fielditer() or d->fieldtime() or d->fieldrange() or m_finished)
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

  if (not benchmark and not (d->It() % rsfreq)) {

    std::vector< std::size_t > meshdata{ /* finished = */ 0, d->MeshId() };
    contribute( meshdata, CkReduction::nop,
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

    auto meshid = d->MeshId();
    d->contribute( sizeof(std::size_t), &meshid, CkReduction::nop,
                   CkCallback(CkReductionTarget(Transporter,finish), d->Tr()) );

  }
}

#include "NoWarning/alecg.def.h"
