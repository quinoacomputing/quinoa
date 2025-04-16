// *****************************************************************************
/*!
  \file      src/Inciter/OversetFE.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     OversetFE for a PDE system with continuous Galerkin FE + RK
  \details   OversetFE advances a system of partial differential equations
    using a continuous Galerkin (CG) finite element (FE) spatial discretization
    (using linear shapefunctions on tetrahedron elements) combined with a
    Runge-Kutta (RK) time stepping scheme and overset grids.
  \see The documentation in OversetFE.hpp.
*/
// *****************************************************************************

#include "QuinoaBuildConfig.hpp"
#include "OversetFE.hpp"
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
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern std::vector< CGPDE > g_cgpde;

//! Runge-Kutta coefficients
static const std::array< tk::real, 3 > rkcoef{{ 1.0/3.0, 1.0/2.0, 1.0 }};

} // inciter::

using inciter::OversetFE;

OversetFE::OversetFE( const CProxy_Discretization& disc,
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
       g_inputdeck.get< tag::ncomp >() ),
  m_uc( m_u.nunk(), m_u.nprop()+1 ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_rhsc(),
  m_chBndGrad( Disc()->Bid().size(), m_u.nprop()*3 ),
  m_dirbc(),
  m_chBndGradc(),
  m_blank( m_u.nunk(), 1.0 ),
  m_diag(),
  m_bnorm(),
  m_bnormc(),
  m_symbcnodes(),
  m_farfieldbcnodes(),
  m_symbctri(),
  m_timedepbcnodes(),
  m_timedepbcFn(),
  m_stage( 0 ),
  m_boxnodes(),
  m_edgenode(),
  m_edgeid(),
  m_dtp( m_u.nunk(), 0.0 ),
  m_tp( m_u.nunk(), g_inputdeck.get< tag::t0 >() ),
  m_finished( 0 ),
  m_movedmesh( 0 ),
  m_nusermeshblk( 0 ),
  m_nodeblockid(),
  m_nodeblockidc(),
  m_ixfer(0),
  m_surfForce({{0, 0, 0}}),
  m_surfTorque({{0, 0, 0}}),
  m_centMass({{0, 0, 0}}),
  m_centMassVel({{0, 0, 0}}),
  m_angVelMesh(0)
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
  if (g_inputdeck.get< tag::operator_reorder >()) {

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

  for (std::size_t i=0; i<3; ++i)
    m_centMass[i] = g_inputdeck.get< tag::mesh >()[d->MeshId()].get<
      tag::center_of_mass >()[i];

  // Query/update boundary-conditions-related data structures from user input
  getBCNodes();

  // Activate SDAG wait for initially computing normals, and mesh blocks
  thisProxy[ thisIndex ].wait4norm();
  thisProxy[ thisIndex ].wait4meshblk();

  if (g_inputdeck.get< tag::steady_state >() &&
    g_inputdeck.get< tag::rigid_body_motion >().get< tag::rigid_body_movt >())
    Throw("Rigid body motion cannot be activated for steady state problem");

  d->comfinal();

}
//! [Constructor]

void
OversetFE::getBCNodes()
// *****************************************************************************
// Query/update boundary-conditions-related data structures from user input
// *****************************************************************************
{
  auto d = Disc();

  // Prepare unique set of symmetry BC nodes
  auto sym = d->bcnodes< tag::symmetry >( m_bface, m_triinpoel );
  for (const auto& [s,nodes] : sym)
    m_symbcnodes.insert( begin(nodes), end(nodes) );

  // Prepare unique set of farfield BC nodes
  auto far = d->bcnodes< tag::farfield >( m_bface, m_triinpoel );
  for (const auto& [s,nodes] : far)
    m_farfieldbcnodes.insert( begin(nodes), end(nodes) );

  // If farfield BC is set on a node, will not also set symmetry BC
  for (auto fn : m_farfieldbcnodes) m_symbcnodes.erase(fn);

  // Prepare boundary nodes contiguously accessible from a triangle-face loop
  m_symbctri.resize( m_triinpoel.size()/3, 0 );
  for (std::size_t e=0; e<m_triinpoel.size()/3; ++e)
    if (m_symbcnodes.find(m_triinpoel[e*3+0]) != end(m_symbcnodes))
      m_symbctri[e] = 1;

  // Prepare unique set of time dependent BC nodes
  m_timedepbcnodes.clear();
  m_timedepbcFn.clear();
  const auto& timedep =
    g_inputdeck.get< tag::bc >()[d->MeshId()].get< tag::timedep >();
  if (!timedep.empty()) {
    m_timedepbcnodes.resize(timedep.size());
    m_timedepbcFn.resize(timedep.size());
    std::size_t ib=0;
    for (const auto& bndry : timedep) {
      std::unordered_set< std::size_t > nodes;
      for (const auto& s : bndry.template get< tag::sideset >()) {
        auto k = m_bnode.find(static_cast<int>(s));
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
}

void
OversetFE::norm()
// *****************************************************************************
// Start (re-)computing boundary point-, and dual-face normals
// *****************************************************************************
{
  auto d = Disc();

  // Query nodes at which symmetry BCs are specified
  auto bn = d->bcnodes< tag::symmetry >( m_bface, m_triinpoel );

  // Query nodes at which farfield BCs are specified
  auto far = d->bcnodes< tag::farfield >( m_bface, m_triinpoel );
  // Merge BC data where boundary-point normals are required
  for (const auto& [s,n] : far) bn[s].insert( begin(n), end(n) );

  // Query nodes at which mesh velocity symmetry BCs are specified
  std::unordered_map<int, std::unordered_set< std::size_t >> ms;
  for (const auto& s : g_inputdeck.get< tag::ale, tag::symmetry >()) {
    auto k = m_bface.find(static_cast<int>(s));
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
OversetFE::edfnorm( const tk::UnsMesh::Edge& edge,
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
OversetFE::dfnorm()
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
OversetFE::comdfnorm( const std::unordered_map< tk::UnsMesh::Edge,
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
OversetFE::bnorm( const std::unordered_map< int,
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
OversetFE::comnorm( const std::unordered_map< int,
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
OversetFE::registerReducers()
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
OversetFE::ResumeFromSync()
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
OversetFE::setup()
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
OversetFE::comblk( const std::vector< std::size_t >& gid,
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
OversetFE::continueSetup()
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
  const auto& hist_points = g_inputdeck.get< tag::history_output, tag::point >();
  if (!hist_points.empty()) {
    std::vector< std::string > histnames;
    auto n = g_cgpde[d->MeshId()].histNames();
    histnames.insert( end(histnames), begin(n), end(n) );
    d->histheader( std::move(histnames) );
  }
}
//! [setup]

void
OversetFE::box( tk::real v, const std::vector< tk::real >& blkvols )
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
  // Initialize overset mesh velocity to zero
  auto& u_mesh = d->MeshVel();
  for (std::size_t p=0; p<u_mesh.nunk(); ++p) {
    for (std::size_t i=0; i<3; ++i) u_mesh(p,i) = 0.0;
  }

  // Initialize nodal mesh volumes at previous time step stage
  d->Voln() = d->Vol();

  // Initiate solution transfer (if coupled)
  transferSol();
}

void
OversetFE::transferSol()
// *****************************************************************************
// Transfer solution to other solver and mesh if coupled
// *****************************************************************************
{
  // Set up transfer-flags for receiving mesh
  if (m_ixfer == 1) {
    applySolTransfer(0);
  }
  setTransferFlags(m_ixfer);
  ++m_ixfer;

  // Initiate IC transfer (if coupled)
  Disc()->transfer( m_uc, m_ixfer-1,
    CkCallback(CkIndex_OversetFE::lhs(), thisProxy[thisIndex]) );
}

//! [Compute lhs]
void
OversetFE::lhs()
// *****************************************************************************
// Compute the left-hand side of transport equations
//! \details Also (re-)compute all data structures if the mesh changed.
// *****************************************************************************
{
  // Do corrections in solution based on incoming transfer
  applySolTransfer(1);
  m_ixfer = 0;

  // No need for LHS in OversetFE

  // If mesh moved: (Re-)compute boundary point- and dual-face normals, and
  //   then proceed to stage()
  // If mesh did not move: shortcut to stage()
  if (m_movedmesh || Disc()->Initial()) norm();
  else stage();
}
//! [Compute lhs]

//! [Merge normals and continue]
void
OversetFE::mergelhs()
// *****************************************************************************
// The own and communication portion of the left-hand side is complete
// *****************************************************************************
{
  // Combine own and communicated contributions of normals
  normfinal();

  // Start with time stepping logic
  if (Disc()->Initial()) {
    // Output initial conditions to file and then start time stepping
    writeFields( CkCallback(CkIndex_OversetFE::start(), thisProxy[thisIndex]) );
  }
  else stage();
}
//! [Merge normals and continue]

//! [start]
void
OversetFE::start()
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

void
OversetFE::applySolTransfer(
  std::size_t dirn )
// *****************************************************************************
// \brief Apply the transferred solution to the solution vector based on
//   transfer flags previously set up
//! \param[in] dirn 0 if called from B to O, 1 if called from O to B
// *****************************************************************************
{
  // Change solution only if:
  //   1. undergoing transfer from B to O, and currently on O
  if (dirn == 0 && Disc()->MeshId() != 0) {

    for (auto i : m_farfieldbcnodes) {
      // overset-BC nodes: use transferred solution and blank nodes.
      // the transfer-flag from m_uc is not used since it has been overwritten
      // by Disc()->transfer() with the flag from B
      for (ncomp_t c=0; c<m_u.nprop(); ++c) { // Loop over number of equations
        m_u(i,c) = m_uc(i,c);
      }
      m_blank[i] = 0.0;
    }

  }
  //   2. undergoing transfer from O to B, and currently on B
  else if (dirn == 1 && Disc()->MeshId() == 0) {

    //TODO: index the flag in a better way
    std::size_t iflag = m_uc.nprop()-1;

    // Zero out solution space for nodes with a specific transfer flag set
    for (std::size_t i=0; i<m_uc.nunk(); ++i) { // Check flag value

      if (std::abs(m_uc(i,iflag) - 1.0) < 1e-4) {
        // overset-BC nodes: use transferred solution and blank nodes
        for (ncomp_t c=0; c<m_u.nprop(); ++c) { // Loop over number of equations
          m_u(i,c) = m_uc(i,c);
        }
        m_blank[i] = 0.0;
      }
      else if (std::abs(m_uc(i,iflag) - 2.0) < 1e-4) {
        // hole: blank nodes
        m_blank[i] = 0.0;
      }
      else {
        // do nothing
        m_blank[i] = 1.0;
      }

    }

  }
}

void
OversetFE::setTransferFlags(
  std::size_t dirn )
// *****************************************************************************
//  Set flags informing solution transfer decisions
//! \param[in] dirn 0 if called from B to O, 1 if called from O to B
// *****************************************************************************
{
  // Copy solution and reset flags
  //TODO: index the flag in a better way
  std::size_t iflag = m_uc.nprop()-1;

  for (std::size_t i=0; i<m_u.nunk(); ++i) {
    for (std::size_t c=0; c<m_u.nprop(); ++c) {
      m_uc(i,c) = m_u(i,c);
    }
    // Reset flags
    m_uc(i,iflag) = 0.0;

    // reset blanking coefficient
    m_blank[i] = 1.0;
  }

  // Transfer flags for O to B are based on block-ids that are hardcoded
  // TODO: remove hardcoding

  // Called from transfer-B-to-O
  if (dirn == 0) {
    if (Disc()->MeshId() != 0) {
      // Overset meshes: assign appropriate values to flag
      for (auto i : m_farfieldbcnodes) m_uc(i,iflag) = 1.0;
    }
  }
  // Called from transfer-O-to-B
  else {
    if (Disc()->MeshId() != 0) {
      // Overset meshes: assign appropriate values to flag
      for (const auto& [blid, ndset] : m_nodeblockid) {
        if (blid == 103) {
          for (auto i : ndset) m_uc(i,iflag) = 1.0;
        }
        else if (blid == 104) {
          for (auto i : ndset) m_uc(i,iflag) = 2.0;
        }
      }
    }
  }
}

void
OversetFE::normfinal()
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
OversetFE::BC()
// *****************************************************************************
// Apply boundary conditions
// \details The following BC enforcement changes the initial condition or
//!   updated solution (dependending on when it is called) to ensure strong
//!   imposition of the BCs. This is a matter of choice. Another alternative is
//!   to only apply BCs when computing fluxes at boundary faces, thereby only
//!   weakly enforcing the BCs. The former is conventionally used in continunous
//!   Galerkin finite element methods (such as OversetFE implements), whereas the
//!   latter, in finite volume methods.
// *****************************************************************************
{
  auto d = Disc();
  const auto& coord = d->Coord();

  const auto& bcmesh = g_inputdeck.get< tag::bc >();

  for (const auto& bci : bcmesh) {
    const auto& bcm = bci.get< tag::mesh >();
    for (const auto& im : bcm) {
      // only if this bc is meant for current mesh
      if (im-1 == d->MeshId()) {

        // Query and match user-specified Dirichlet boundary conditions to side sets
        const auto steady = g_inputdeck.get< tag::steady_state >();
        if (steady) for (auto& deltat : m_dtp) deltat *= rkcoef[m_stage];
        m_dirbc = match( d->MeshId(), m_u.nprop(), d->T(), rkcoef[m_stage] * d->Dt(),
                         m_tp, m_dtp, d->Coord(), d->Lid(), m_bnode,
                       /* increment = */ false );
        if (steady) for (auto& deltat : m_dtp) deltat /= rkcoef[m_stage];

        // Apply Dirichlet BCs
        for (const auto& [b,bc] : m_dirbc)
          for (ncomp_t c=0; c<m_u.nprop(); ++c)
            if (bc[c].first) m_u(b,c) = bc[c].second;

        // Apply symmetry BCs
        g_cgpde[d->MeshId()].symbc( m_u, d->meshvel(), coord, m_bnorm, m_symbcnodes );

        // Apply farfield BCs
        if (bci.get< tag::farfield >().empty() || (d->MeshId() == 0)) {
          g_cgpde[d->MeshId()].farfieldbc( m_u, coord, m_bnorm, m_farfieldbcnodes );
        }

        // Apply user defined time dependent BCs
        g_cgpde[d->MeshId()].timedepbc( d->T(), m_u, m_timedepbcnodes,
          m_timedepbcFn );
      }
    }
  }
}

void
OversetFE::next()
// *****************************************************************************
// Continue to next time step
// *****************************************************************************
{
  dt();
}

void
OversetFE::dt()
// *****************************************************************************
// Compute time step size
// *****************************************************************************
{
  tk::real mindt = std::numeric_limits< tk::real >::max();

  auto const_dt = g_inputdeck.get< tag::dt >();
  auto eps = std::numeric_limits< tk::real >::epsilon();

  auto d = Disc();

  // use constant dt if configured
  if (std::abs(const_dt) > eps) {

    mindt = const_dt;

  } else {      // compute dt based on CFL

    //! [Find the minimum dt across all PDEs integrated]
    if (g_inputdeck.get< tag::steady_state >()) {

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
    //! [Find the minimum dt across all PDEs integrated]

  }

  //! [Advance]
  // Actiavate SDAG waits for next time step stage
  thisProxy[ thisIndex ].wait4grad();
  thisProxy[ thisIndex ].wait4rhs();

  // Compute own portion of force on boundary for overset mesh rigid body motion
  std::vector< tk::real > F(6, 0.0);
  if (g_inputdeck.get< tag::rigid_body_motion >().get< tag::rigid_body_movt >()
    && d->MeshId() > 0) {
    g_cgpde[d->MeshId()].bndPressureInt( d->Coord(), m_triinpoel, m_symbctri,
      m_u, m_centMass, F );
  }

  // Tuple-reduction for min-dt and sum-F
  int tupleSize = 7;
  CkReduction::tupleElement advancingData[] = {
    CkReduction::tupleElement (sizeof(tk::real), &mindt, CkReduction::min_double),
    CkReduction::tupleElement (sizeof(tk::real), &F[0], CkReduction::sum_double),
    CkReduction::tupleElement (sizeof(tk::real), &F[1], CkReduction::sum_double),
    CkReduction::tupleElement (sizeof(tk::real), &F[2], CkReduction::sum_double),
    CkReduction::tupleElement (sizeof(tk::real), &F[3], CkReduction::sum_double),
    CkReduction::tupleElement (sizeof(tk::real), &F[4], CkReduction::sum_double),
    CkReduction::tupleElement (sizeof(tk::real), &F[5], CkReduction::sum_double)
  };
  CkReductionMsg* advMsg =
    CkReductionMsg::buildFromTuple(advancingData, tupleSize);

  // Contribute to minimum dt across all chares, find minimum dt across all
  // meshes, and eventually broadcast to OversetFE::advance()
  CkCallback cb(CkReductionTarget(Transporter,collectDtAndForces), d->Tr());
  advMsg->setCallback(cb);
  contribute(advMsg);
  //! [Advance]
}

void
OversetFE::advance( tk::real newdt, std::array< tk::real, 6 > F )
// *****************************************************************************
// Advance equations to next time step
//! \param[in] newdt The smallest dt across the whole problem
//! \param[in] F Total surface force on this mesh
// *****************************************************************************
{
  auto d = Disc();

  // Set new time step size
  if (m_stage == 0) d->setdt( newdt );

  for (std::size_t i=0; i<3; ++i) {
    m_surfForce[i] = F[i];
    m_surfTorque[i] = F[i+3];
  }

  // Compute gradients for next time step
  chBndGrad();
}

void
OversetFE::chBndGrad()
// *****************************************************************************
// Compute nodal gradients at chare-boundary nodes. Gradients at internal nodes
// are calculated locally as needed and are not stored.
// *****************************************************************************
{
  auto d = Disc();

  // Compute own portion of gradients for all equations
  g_cgpde[d->MeshId()].chBndGrad( d->Coord(), d->Inpoel(), m_bndel, d->Gid(),
    d->Bid(), m_u, m_chBndGrad );

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
OversetFE::comChBndGrad( const std::vector< std::size_t >& gid,
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
OversetFE::rhs()
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

  const auto steady = g_inputdeck.get< tag::steady_state >();

  // Compute own portion of right-hand side for all equations
  auto prev_rkcoef = m_stage == 0 ? 0.0 : rkcoef[m_stage-1];
  if (steady)
    for (std::size_t p=0; p<m_tp.size(); ++p) m_tp[p] += prev_rkcoef * m_dtp[p];
  g_cgpde[d->MeshId()].rhs( d->T() + prev_rkcoef * d->Dt(), d->Coord(), d->Inpoel(),
          m_triinpoel, d->Gid(), d->Bid(), d->Lid(), m_dfn, m_psup, m_esup,
          m_symbctri, d->Vol(), m_edgenode, m_edgeid,
          m_boxnodes, m_chBndGrad, m_u, d->meshvel(), m_tp, d->Boxvol(),
          m_rhs );
  if (steady)
    for (std::size_t p=0; p<m_tp.size(); ++p) m_tp[p] -= prev_rkcoef * m_dtp[p];

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
OversetFE::comrhs( const std::vector< std::size_t >& gid,
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
OversetFE::solve()
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
    d->UpdateCoordn();
  }

  // Explicit time-stepping using RK3
  const auto steady = g_inputdeck.get< tag::steady_state >();
  for (std::size_t i=0; i<m_u.nunk(); ++i) {
    // time-step
    auto dtp = d->Dt();
    if (steady) dtp = m_dtp[i];

    for (ncomp_t c=0; c<m_u.nprop(); ++c)
      m_u(i,c) = m_un(i,c) + m_blank[i] * rkcoef[m_stage] * dtp * m_rhs(i,c)
        / d->Vol()[i];
  }

  // Kinematics for rigid body motion of overset meshes
  if (g_inputdeck.get< tag::rigid_body_motion >().get< tag::rigid_body_movt >()
    && d->MeshId() > 0) {

    // Remove symmetry directions if 3 DOF motion
    if (g_inputdeck.get< tag::rigid_body_motion >().get< tag::rigid_body_dof >()
      == 3) {

      auto sym_dir =
        g_inputdeck.get< tag::rigid_body_motion >().get< tag::symmetry_plane >();

      m_surfForce[sym_dir] = 0.0;
      for (std::size_t i=0; i<3; ++i) {
        if (i != sym_dir) m_surfTorque[i] = 0.0;
      }
    }

    // Mark if mesh moved
    if (std::sqrt(tk::dot(m_surfForce, m_surfForce)) > 1e-12 ||
      std::sqrt(tk::dot(m_surfTorque, m_surfTorque)) > 1e-12)
      m_movedmesh = 1;
    else
      m_movedmesh = 0;

    if (m_movedmesh == 1) {
      auto mass_mesh =
        g_inputdeck.get< tag::mesh >()[d->MeshId()].get< tag::mass >();
      auto mI_mesh = g_inputdeck.get< tag::mesh >()[d->MeshId()].get<
        tag::moment_of_inertia >();
      auto dtp = rkcoef[m_stage] * d->Dt();
      auto sym_dir =
        g_inputdeck.get< tag::rigid_body_motion >().get< tag::symmetry_plane >();

      auto pi = 4.0*std::atan(1.0);

      // mesh acceleration
      std::array< tk::real, 3 > a_mesh;
      for (std::size_t i=0; i<3; ++i) a_mesh[i] = m_surfForce[i] / mass_mesh;
      auto alpha_mesh = m_surfTorque[sym_dir]/mI_mesh; // angular acceleration

      auto& u_mesh = d->MeshVel();

      for (std::size_t p=0; p<u_mesh.nunk(); ++p) {

        // rotation (this is currently only configured for planar motion)
        // ---------------------------------------------------------------------
        std::array< tk::real, 3 > rCM{{
          d->Coord()[0][p] - m_centMass[0],
          d->Coord()[1][p] - m_centMass[1],
          d->Coord()[2][p] - m_centMass[2] }};

        // obtain tangential velocity
        tk::real r_mag(0.0);
        for (std::size_t i=0; i<3; ++i) {
          if (i != sym_dir) r_mag += rCM[i]*rCM[i];
        }
        r_mag = std::sqrt(r_mag);
        auto a_tgt = alpha_mesh*r_mag;

        // get the other two directions
        auto i1 = (sym_dir+1)%3;
        auto i2 = (sym_dir+2)%3;

        // project tangential velocity to these two directions
        auto theta = std::atan2(rCM[i2],rCM[i1]);
        auto a1 = a_tgt*std::cos((pi/2.0)+theta);
        auto a2 = a_tgt*std::sin((pi/2.0)+theta);

        // angle of rotation
        auto dtheta = m_angVelMesh*dtp + 0.5*alpha_mesh*dtp*dtp;

        // rectilinear motion
        // ---------------------------------------------------------------------

        // add contribution of rotation to mesh displacement
        std::array< tk::real, 3 > angles{{ 0, 0, 0 }};
        angles[sym_dir] = dtheta * 180.0/pi;
        tk::rotatePoint(angles, rCM);

        tk::real dsT, dsR;

        for (std::size_t i=0; i<3; ++i) {
          // mesh displacement from translation
          dsT = m_centMassVel[i]*dtp + 0.5*a_mesh[i]*dtp*dtp;
          // mesh displacement from rotation
          dsR = rCM[i] + m_centMass[i] - d->Coord()[i][p];
          // add both contributions
          d->Coord()[i][p] = d->Coordn()[i][p] + dsT + dsR;
          // mesh velocity change from translation
          u_mesh(p,i) += a_mesh[i]*dtp;
        }

        // add contribution of rotation to mesh velocity
        u_mesh(p,i1) += a1*dtp;
        u_mesh(p,i2) += a2*dtp;
      }

      // update angular velocity
      m_angVelMesh += alpha_mesh*dtp;

      // move center of mass
      for (std::size_t i=0; i<3; ++i) {
        m_centMass[i] += m_centMassVel[i]*dtp + 0.5*a_mesh[i]*dtp*dtp;
        m_centMassVel[i] += a_mesh[i]*dtp;  // no rotational component
      }
    }

  }

  // Apply boundary-conditions
  BC();

  // Increment Runge-Kutta stage counter
  ++m_stage;

  // Activate SDAG wait for next time step stage
  thisProxy[ thisIndex ].wait4grad();
  thisProxy[ thisIndex ].wait4rhs();

  // Compute diagnostics, and finish-up time step (if m_stage == 3)
  bool diag_computed(false);
  if (m_stage == 3) {
    // Compute diagnostics, e.g., residuals
    diag_computed = m_diag.compute( *d, m_u, m_un, m_bnorm,
                                    m_symbcnodes, m_farfieldbcnodes );
    // Increase number of iterations and physical time
    d->next();
    // Advance physical time for local time stepping
    if (g_inputdeck.get< tag::steady_state >())
      for (std::size_t i=0; i<m_u.nunk(); ++i) m_tp[i] += m_dtp[i];
  }
  // Continue to finish-up time-step-stage
  // Note: refine is called via a bcast if diag_computed == true
  if (!diag_computed) refine( std::vector< tk::real >( m_u.nprop(), 1.0 ) );
}

//! [Refine]
void
OversetFE::refine( const std::vector< tk::real >& l2res )
// *****************************************************************************
// Finish up end of time-step procedures and continue to moving mesh
//! \param[in] l2res L2-norms of the residual for each scalar component
//!   computed across the whole problem
// *****************************************************************************
{
  auto d = Disc();

  if (m_stage == 3) {
    const auto steady = g_inputdeck.get< tag::steady_state >();
    const auto residual = g_inputdeck.get< tag::residual >();
    const auto rc = g_inputdeck.get< tag::rescomp >() - 1;

    if (m_movedmesh) {
      d->Itf() = 0;  // Zero field output iteration count if mesh moved
      ++d->Itr();    // Increase number of iterations with a change in the mesh
    }

    if (steady) {

      // this is the last time step if max time of max number of time steps
      // reached or the residual has reached its convergence criterion
      if (d->finished() or l2res[rc] < residual) m_finished = 1;

    } else {

      // this is the last time step if max time or max iterations reached
      if (d->finished()) m_finished = 1;

    }
  }

  if (m_movedmesh) {
    // Normals need to be recomputed if overset mesh has been moved
    thisProxy[ thisIndex ].wait4norm();
  }

  // Start solution transfer
  transferSol();
}
//! [Refine]

//! [stage]
void
OversetFE::stage()
// *****************************************************************************
// Evaluate whether to continue with next time step stage
// *****************************************************************************
{
  // if not all Runge-Kutta stages complete, continue to next time stage,
  // otherwise start next time step
  if (m_stage == 3) {
    // output field data and start with next time step
    out();
  }
  else {
    // start with next time-step stage
    chBndGrad();
  }
}
//! [stage]

void
OversetFE::writeFields( CkCallback c )
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

    //// if coupled: depvars: src:'a', dst:'b','c',...
    //char depvar = 0;
    //if (not d->Transfers().empty()) {
    //  depvar = 'a' + static_cast< char >( d->MeshId() );
    //}

    // Query fields names requested by user
    auto nodefieldnames = numericFieldNames( tk::Centering::NODE );

    // Collect field output from numerical solution requested by user
    auto nodefields = numericFieldOutput( m_u, tk::Centering::NODE,
      g_cgpde[Disc()->MeshId()].OutVarFn(), m_u );

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
    auto so = g_cgpde[d->MeshId()].surfOutput( tk::bfacenodes(m_bface,
      m_triinpoel), m_u );
    nodesurfs.insert( end(nodesurfs), begin(so), end(so) );

    // Collect elemental block and surface field names from PDEs integrated
    auto elemsurfnames = nodesurfnames;

    // Collect elemental block and surface field solution
    std::vector< std::vector< tk::real > > elemsurfs;
    auto eso = g_cgpde[d->MeshId()].elemSurfOutput( m_bface, m_triinpoel, m_u );
    elemsurfs.insert( end(elemsurfs), begin(eso), end(eso) );

    Assert( nodefieldnames.size() == nodefields.size(), "Size mismatch" );

    // Send mesh and fields data (solution dump) for output to file
    d->write( d->Inpoel(), coord, m_bface, tk::remap(m_bnode,d->Lid()),
              m_triinpoel, {}, nodefieldnames, elemsurfnames,
              nodesurfnames, {}, nodefields, elemsurfs, nodesurfs, c );

  }
}

void
OversetFE::out()
// *****************************************************************************
// Output mesh field data and continue to next time step
// *****************************************************************************
{
  auto d = Disc();

  // Output time history
  if (d->histiter() or d->histtime() or d->histrange()) {
    std::vector< std::vector< tk::real > > hist;
    auto h = g_cgpde[d->MeshId()].histOutput( d->Hist(), d->Inpoel(), m_u );
    hist.insert( end(hist), begin(h), end(h) );
    d->history( std::move(hist) );
  }

  // Output field data
  if (d->fielditer() or d->fieldtime() or d->fieldrange() or m_finished)
    writeFields(CkCallback( CkIndex_OversetFE::step(), thisProxy[thisIndex]) );
  else
    step();
}

void
OversetFE::evalLB( int nrestart )
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
OversetFE::evalRestart()
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
OversetFE::step()
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

#include "NoWarning/oversetfe.def.h"
