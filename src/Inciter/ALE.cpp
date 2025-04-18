// *****************************************************************************
/*!
  \file      src/Inciter/ALE.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Definitions file for distributed ALE mesh motion
  \details   Definitions file for asynchronous distributed
             arbitrary Lagrangian-Eulerian (ALE) mesh motion using Charm++.
*/
// *****************************************************************************

#include "ALE.hpp"
#include "DerivedData.hpp"
#include "Vector.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::ALE;

ALE::ALE( const tk::CProxy_BiCG& bicgproxy,
          const std::array< std::vector< tk::real >, 3 >& coord,
          const std::vector< std::size_t >& inpoel,
          const std::vector< std::size_t >& gid,
          const std::unordered_map< std::size_t, std::size_t >& lid,
          const tk::NodeCommMap& nodeCommMap ) :
  m_bicg( bicgproxy ),
  m_done(),
  m_nvort( 0 ),
  m_ndiv( 0 ),
  m_npot( 0 ),
  m_nwf( 0 ),
  m_nodeCommMap( nodeCommMap ),
  m_lid( lid ),
  m_coord0( coord ),
  m_inpoel( inpoel ),
  m_vol0(),
  m_vol(),
  m_it( 0 ),
  m_t( 0.0 ),
  m_w( gid.size(), 3 ),
  m_wf( gid.size(), 3 ),
  m_wfc(),
  m_veldiv(),
  m_veldivc(),
  m_gradpot(),
  m_gradpotc(),
  m_vorticity(),
  m_vorticityc(),
  m_meshveldirbcnodes(),
  m_meshvelsymbcnodes(),
  m_move( moveCfg() )
// *****************************************************************************
//  Constructor
//! \param[in] bicgproxy Distributed Conjugrate Gradients linear
//!   solver proxy (bound to ALE proxy)
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] gid Local->global node id map
//! \param[in] lid Global->locbal node id map
//! \param[in] nodeCommMap Node communication map
// *****************************************************************************
{
  auto smoother = g_inputdeck.get< tag::ale, tag::smoother >();

  if (smoother == ctr::MeshVelocitySmootherType::LAPLACE)
    m_bicg[ thisIndex ].        // solve for mesh velocity
      insert( Laplacian( 3, coord ), gid, m_lid, m_nodeCommMap );
  else if (smoother == ctr::MeshVelocitySmootherType::HELMHOLTZ)
    m_bicg[ thisIndex ].        // solve for scalar potential
      insert( Laplacian( 1, coord ), gid, m_lid, m_nodeCommMap );

  // Zero ALE mesh velocity
  m_w.fill( 0.0 );

  // Activate SDAG wait for initially computing prerequisites for ALE
  thisProxy[ thisIndex ].wait4vel();
  thisProxy[ thisIndex ].wait4pot();
  thisProxy[ thisIndex ].wait4for();
}

std::tuple< tk::CSR, std::vector< tk::real >, std::vector< tk::real > >
ALE::Laplacian( std::size_t ncomp,
                const std::array< std::vector< tk::real >, 3 >& coord ) const
// *****************************************************************************
// Generate {A,x,b} for Laplacian mesh velocity smoother
//! \param[in] ncomp Number of scalar components
//! \param[in] coord Mesh node coordinates
//! \return {A,x,b} with a Laplacian, unknown, and rhs initialized with zeros
//! \see Waltz, et al. "A three-dimensional finite element arbitrary
//!   Lagrangian-Eulerian method for shock hydrodynamics on unstructured
//!   grids", Computers& Fluids, 2013, and Bakosi, et al. "Improved ALE mesh
//!   velocities for complex flows, International Journal for Numerical Methods
//!   in Fluids, 2017.
// *****************************************************************************
{
  tk::CSR A( ncomp, tk::genPsup(m_inpoel,4,tk::genEsup(m_inpoel,4)) );

  const auto& X = coord[0];
  const auto& Y = coord[1];
  const auto& Z = coord[2];

  // fill matrix with Laplacian
  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    // access node IDs
    const std::array< std::size_t, 4 >
      N{{ m_inpoel[e*4+0], m_inpoel[e*4+1], m_inpoel[e*4+2], m_inpoel[e*4+3] }};
    // compute element Jacobi determinant
    const std::array< tk::real, 3 >
      ba{{ X[N[1]]-X[N[0]], Y[N[1]]-Y[N[0]], Z[N[1]]-Z[N[0]] }},
      ca{{ X[N[2]]-X[N[0]], Y[N[2]]-Y[N[0]], Z[N[2]]-Z[N[0]] }},
      da{{ X[N[3]]-X[N[0]], Y[N[3]]-Y[N[0]], Z[N[3]]-Z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );        // J = 6V
    Assert( J > 0, "Element Jacobian non-positive" );

    // shape function derivatives, nnode*ndim [4][3]
    std::array< std::array< tk::real, 3 >, 4 > grad;
    grad[1] = tk::crossdiv( ca, da, J );
    grad[2] = tk::crossdiv( da, ba, J );
    grad[3] = tk::crossdiv( ba, ca, J );
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

    for (std::size_t a=0; a<4; ++a)
      for (std::size_t k=0; k<3; ++k)
         for (std::size_t b=0; b<4; ++b)
           for (std::size_t i=0; i<ncomp; ++i)
             A(N[a],N[b],i) -= J/6 * grad[a][k] * grad[b][k];
  }

  auto npoin = coord[0].size();
  std::vector< tk::real > b(npoin*ncomp,0.0), x(npoin*ncomp,0.0);

  return { std::move(A), std::move(x), std::move(b) };
}

decltype(ALE::m_move)
ALE::moveCfg()
// *****************************************************************************
// Initialize user-defined functions for ALE moving sides
//! \details This function fills in only part of the data structure
//!   returned, containing the user-defined functions in discrete form that will
//!   be sampled in time. The node lists will be initialized later.
// *****************************************************************************
{
  decltype(m_move) cfg;

  for (const auto& m : g_inputdeck.get< tag::ale, tag::move >()) {
    const auto& fn = m.get< tag::fn >();
    Assert( fn.size() % 4 == 0, "Incomplete user-defined function" );
    cfg.emplace_back();
    // store user-defined function type
    std::get<0>(cfg.back()) = m.get< tag::fntype >();
    // store user-defined function discrete data
    for (std::size_t i=0; i<fn.size()/4; ++i)
      std::get<1>(cfg.back()).
        push_back( {{ fn[i*4+0], fn[i*4+1], fn[i*4+2], fn[i*4+3] }} );
  }

  return cfg;
}

void
ALE::meshvelBnd(
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::map< int, std::vector< std::size_t > >& bnode,
  const std::vector< std::size_t >& triinpoel )
// *****************************************************************************
//  Query mesh velocity boundary conditions node lists and node list at which
//  ALE moves boundaries
//! \param[in] bface Boundary-faces mapped to side sets used in the input file
//! \param[in] bnode Boundary-node lists mapped to side sets used in input file
//! \param[in] triinpoel Boundary-face connectivity where BCs set
// *****************************************************************************
{
  // Prepare unique set of mesh velocity Dirichlet BC nodes
  tk::destroy( m_meshveldirbcnodes );
  std::unordered_map<int, std::unordered_set< std::size_t >> meshveldirbcnodes;
  for (const auto& s : g_inputdeck.template get< tag::ale, tag::dirichlet >()) {
    auto k = bface.find(static_cast<int>(s));
    if (k != end(bface)) {
      auto& n = meshveldirbcnodes[ k->first ];  // associate set id
      for (auto f : k->second) {                // face ids on side set
        n.insert( triinpoel[f*3+0] );
        n.insert( triinpoel[f*3+1] );
        n.insert( triinpoel[f*3+2] );
      }
    }
  }
  for (const auto& [s,nodes] : meshveldirbcnodes)
    m_meshveldirbcnodes.insert( begin(nodes), end(nodes) );

  // Prepare unique set of mesh velocity symmetry BC nodes. Note that somewhat
  // counter-intuitively, we interrogate the boundary nodes instead of boundary
  // faces here. This is because if we query the boundary faces, then we will
  // get the mathematically correctly defined finite discrete surfaces
  // (triangles) where mesh velocity symmetry BCs are configured by the user.
  // However, in parallel, decomposing the domain and the boundary in various
  // ways can produce situations on the boundary where boundary nodes are part
  // of the given side set for mesh velocity symmetry BCs but not a full
  // triangle face because, not all 3 nodes lie on the boundary. Thus
  // interrogating the boundary nodes will be a superset and will include those
  // nodes that are part of imposing symmetry BCs on nodes of faces that are
  // only partial due to domain decomposition.
  tk::destroy( m_meshvelsymbcnodes );
  std::unordered_map<int, std::unordered_set< std::size_t >> meshvelsymbcnodes;
  for (const auto& s : g_inputdeck.template get< tag::ale, tag::symmetry >()) {
    auto k = bnode.find(static_cast<int>(s));
    if (k != end(bnode)) {
      auto& n = meshvelsymbcnodes[ k->first ];  // associate set id
      for (auto g : k->second) {                // node ids on side set
        n.insert( tk::cref_find(m_lid,g) );     // store local ids
      }
    }
  }
  for (const auto& [s,nodes] : meshvelsymbcnodes)
    m_meshvelsymbcnodes.insert( begin(nodes), end(nodes) );

  // Prepare unique sets of boundary nodes at which ALE moves the boundary
  // based on user-defined functions.
  std::unordered_map< int, std::unordered_set< std::size_t > > movenodes;
  std::size_t i = 0;
  for (const auto& m : g_inputdeck.get< tag::ale, tag::move >()) {
    for (const auto& s : m.get< tag::sideset >()) {
      auto k = bnode.find(static_cast<int>(s));
      if (k != end(bnode)) {
        auto& n = movenodes[ k->first ];        // associate set id
        for (auto g : k->second) {              // node ids on side set
          n.insert( tk::cref_find(m_lid,g) );   // store local ids
        }
      }
    }
    // store all nodes from multiple side sets moved by this usrdef fn
    auto& n = std::get<2>(m_move[i]);
    n.clear();
    for (const auto& [s,nodes] : movenodes) n.insert(begin(nodes), end(nodes));
    // increment move ... end configuration block counter
    ++i;
  }
}

bool
ALE::move( std::size_t i ) const
// *****************************************************************************
// Find Dirichlet BCs on mesh velocity with prescribed movement
//! \param[in] i Local node id to check
//! \return True of node falls on a boundary that is prescribed to move
// *****************************************************************************
{
  for (const auto& m : m_move)
    if (std::get<2>(m).find(i) != end(std::get<2>(m)))
      return true;

  return false;
}

bool
ALE::converged() const
// *****************************************************************************
//  Query the solution of the Conjugrate Gradients linear solver
//! \return Solution to the Conjugate Gradients linear solve
// *****************************************************************************
{
  return m_bicg[ thisIndex ].ckLocal()->converged();
}

void
ALE::start(
  const tk::UnsMesh::Coords vel,
  const std::vector< tk::real >& soundspeed,
  CkCallback done,
  const std::array< std::vector< tk::real >, 3 >& coord,
  const tk::UnsMesh::Coords coordn,
  const std::vector< tk::real >& vol0,
  const std::vector< tk::real >& vol,
  const std::unordered_map< int,
    std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& bnorm,
  std::size_t initial,
  std::size_t it,
  tk::real t,
  tk::real adt )
// *****************************************************************************
// Start computing new mesh velocity for ALE mesh motion
//! \param[in] vel Fluid velocity at mesh nodes
//! \param[in] soundspeed Speed of sound at mesh nodes
//! \param[in] done Function to continue with when mesh velocity has been
//!   computed
//! \param[in] coord Mesh node coordinates
//! \param[in] coordn Mesh node coordinates at the previous time step
//! \param[in] vol0 Nodal mesh volumes at t=t0
//! \param[in] vol Nodal mesh volumes
//! \param[in] bnorm Face normals in boundary points associated to side sets
//! \param[in] initial Nonzero during the first time step stage, zero otherwise
//! \param[in] it Iteration count
//! \param[in] t Physics time
//! \param[in] adt alpha*dt of the RK time step
// *****************************************************************************
{
  m_done = done;

  m_coord = coord;
  m_soundspeed = soundspeed;
  m_vol0 = vol0;
  m_vol = vol;
  m_bnorm = bnorm;
  m_it = it;
  m_t = t;
  m_adt = adt;

  // assign mesh velocity
  auto meshveltype = g_inputdeck.get< tag::ale, tag::mesh_velocity >();
  if (meshveltype == ctr::MeshVelocityType::SINE) {

    // prescribe mesh velocity with a sine function during setup
    if (initial)
      for (std::size_t i=0; i<m_w.nunk(); ++i)
        m_w(i,0) = std::pow( std::sin(coord[0][i]*M_PI), 2.0 );

  } else if (meshveltype == ctr::MeshVelocityType::FLUID) {

    // equate mesh velocity with fluid velocity
    for (auto j : g_inputdeck.get< tag::ale, tag::mesh_motion >())
      for (std::size_t i=0; i<vel[j].size(); ++i)
        m_w(i,j) = vel[j][i];

  } else if (meshveltype == ctr::MeshVelocityType::USER_DEFINED) {

    // assign mesh velocity to sidesets from user-defined functions
    for (const auto& m : m_move)
      if (std::get<0>(m) == tk::ctr::UserTableType::VELOCITY) {
        auto meshvel = tk::sample<3>( t, std::get<1>(m) );
        for (auto i : std::get<2>(m))
          for (auto j : g_inputdeck.get< tag::ale, tag::mesh_motion >())
            m_w(i,j) = meshvel[j];
      } else if (std::get<0>(m) == tk::ctr::UserTableType::POSITION) {
        auto eps = std::numeric_limits< tk::real >::epsilon();
        if (adt > eps) {      // dt == 0 during setup
          auto pos = tk::sample<3>( t+adt, std::get<1>(m) );
          for (auto i : std::get<2>(m))
            for (auto j : g_inputdeck.get< tag::ale, tag::mesh_motion >())
              m_w(i,j) = (m_coord0[j][i] + pos[j] - coordn[j][i]) / adt;
        }
      }

  }

  // start computing the fluid vorticity
  m_vorticity = tk::curl( coord, m_inpoel, vel );
  // communicate vorticity sums to other chares on chare-boundary
  if (m_nodeCommMap.empty()) {
    comvort_complete();
  } else {
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::array< tk::real, 3 > > v( n.size() );
      std::size_t j = 0;
      for (auto i : n) {
        auto lid = tk::cref_find( m_lid, i );
        v[j][0] = m_vorticity[0][lid];
        v[j][1] = m_vorticity[1][lid];
        v[j][2] = m_vorticity[2][lid];
        ++j;
      }
      thisProxy[c].comvort( std::vector<std::size_t>(begin(n),end(n)), v );
    }
  }
  ownvort_complete();

  // start computing the fluid velocity divergence
  m_veldiv = tk::div( coord, m_inpoel, vel );
  // communicate vorticity sums to other chares on chare-boundary
  if (m_nodeCommMap.empty()) {
    comdiv_complete();
  } else {
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< tk::real > v( n.size() );
      std::size_t j = 0;
      for (auto i : n) v[j++] = m_veldiv[ tk::cref_find( m_lid, i ) ];
      thisProxy[c].comdiv( std::vector<std::size_t>(begin(n),end(n)), v );
    }
  }
  owndiv_complete();
}

void
ALE::comvort( const std::vector< std::size_t >& gid,
              const std::vector< std::array< tk::real, 3 > >& v )
// *****************************************************************************
//  Receive contributions to vorticity on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] v Partial contributions to chare-boundary nodes
// *****************************************************************************
{
  Assert( v.size() == gid.size(), "Size mismatch" );
  using tk::operator+=;
  for (std::size_t i=0; i<gid.size(); ++i) m_vorticityc[ gid[i] ] += v[i];

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_nvort == m_nodeCommMap.size()) {
    m_nvort = 0;
    comvort_complete();
  }
}

void
ALE::comdiv( const std::vector< std::size_t >& gid,
             const std::vector< tk::real >& v )
// *****************************************************************************
//  Receive contributions to velocity divergence on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive RHS contributions
//! \param[in] v Partial contributions to chare-boundary nodes
// *****************************************************************************
{
  Assert( v.size() == gid.size(), "Size mismatch" );
  using tk::operator+=;
  for (std::size_t i=0; i<gid.size(); ++i) m_veldivc[ gid[i] ] += v[i];

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_ndiv == m_nodeCommMap.size()) {
    m_ndiv = 0;
    comdiv_complete();
  }
}

void
ALE::mergevel()
// *****************************************************************************
// Finalize computing fluid vorticity and velocity divergence for ALE
// *****************************************************************************
{
  // combine own and communicated contributions to vorticity
  for (const auto& [g,v] : m_vorticityc) {
    auto lid = tk::cref_find( m_lid, g );
    m_vorticity[0][lid] += v[0];
    m_vorticity[1][lid] += v[1];
    m_vorticity[2][lid] += v[2];
  }
  // clear fluid vorticity receive buffer
  tk::destroy(m_vorticityc);

  // finish computing vorticity dividing the weak sum by the nodal volumes
  for (std::size_t j=0; j<3; ++j)
    for (std::size_t p=0; p<m_vorticity[j].size(); ++p)
      m_vorticity[j][p] /= m_vol[p];

  // compute vorticity magnitude
  for (std::size_t p=0; p<m_vorticity[0].size(); ++p)
    m_vorticity[0][p] =
      tk::length( m_vorticity[0][p], m_vorticity[1][p], m_vorticity[2][p] );

  // get rid of the y and z vorticity components, since we just overwrote the
  // x component with the magnitude
  tk::destroy( m_vorticity[1] );
  tk::destroy( m_vorticity[2] );

  // compute max vorticity magnitude
  auto maxv =
    *std::max_element( m_vorticity[0].cbegin(), m_vorticity[0].cend() );

  // combine own and communicated contributions to velocidy divergence
  for (const auto& [g,v] : m_veldivc)
    m_veldiv[ tk::cref_find( m_lid, g ) ] += v;
  // clear velocity divergence receive buffer
  tk::destroy(m_veldivc);

  // finish computing velocity divergence dividing weak sum by the nodal volumes
  for (std::size_t p=0; p<m_veldiv.size(); ++p) m_veldiv[p] /= m_vol[p];

  // Compute max vorticity magnitude across all chares
  contribute( sizeof(tk::real), &maxv, CkReduction::max_double,
              CkCallback(CkReductionTarget(ALE,meshvelbc), thisProxy) );
}

void
ALE::meshvelbc( tk::real maxv )
// *****************************************************************************
// Apply mesh velocity smoother boundary conditions for ALE mesh motion
//! \param[in] maxv The largest vorticity magnitude across the whole problem
// *****************************************************************************
{
  std::size_t ignorebc = false;

  // smooth mesh velocity if needed
  auto smoother = g_inputdeck.get< tag::ale, tag::smoother >();

  if (smoother == ctr::MeshVelocitySmootherType::LAPLACE) {

    // scale mesh velocity with a function of the fluid vorticity
    if (maxv > 1.0e-8) {
      auto mult = g_inputdeck.get< tag::ale, tag::vortmult >();
      for (auto j : g_inputdeck.get< tag::ale, tag::mesh_motion >())
        for (std::size_t p=0; p<m_vorticity[0].size(); ++p)
          m_w(p,j) *= std::max( 0.0, 1.0 - mult*m_vorticity[0][p]/maxv );
    }

    // Set mesh velocity smoother linear solve boundary conditions
    std::unordered_map< std::size_t,
      std::vector< std::pair< bool, tk::real > > > wbc;

    // Dirichlet BCs on mesh velocity with prescribed movement
    for (const auto& m : m_move)
      if (std::get<0>(m) == tk::ctr::UserTableType::VELOCITY) {
        auto meshvel = tk::sample<3>( m_t, std::get<1>(m) );
        for (auto i : std::get<2>(m))
          for (auto j : g_inputdeck.get< tag::ale, tag::mesh_motion >())
            m_w(i,j) = meshvel[j];
      } else if (std::get<0>(m) == tk::ctr::UserTableType::POSITION) {
        auto eps = std::numeric_limits< tk::real >::epsilon();
        if (m_adt > eps) {
          ignorebc = m_it > 0;
          for (auto i : std::get<2>(m))
            if (m_meshveldirbcnodes.find(i) != end(m_meshveldirbcnodes))
              wbc[i] = {{ {false,0}, {false,0}, {false,0} }};
        }
      }

    // Dirichlet BCs where user specified mesh velocity BCs
    for (auto i : m_meshveldirbcnodes)
      if (not move(i)) wbc[i] = {{ {true,0}, {true,0}, {true,0} }};

    // initialize mesh velocity smoother linear solver
    m_bicg[ thisIndex ].ckLocal()->
      init( m_w.flat(), {}, wbc, ignorebc,
            CkCallback(CkIndex_ALE::applied(nullptr), thisProxy[thisIndex]) );

  } else if (smoother == ctr::MeshVelocitySmootherType::HELMHOLTZ) {

    // Set scalar potential linear solve boundary conditions
    std::unordered_map< std::size_t,
      std::vector< std::pair< bool, tk::real > > > pbc;

    // Dirichlet BCs where user specified mesh velocity BCs
    for (auto i : m_meshveldirbcnodes) pbc[i] = {{ {true,0} }};

    // prepare velocity divergence as weak sum required for Helmholtz solve
    decltype(m_veldiv) wveldiv = m_veldiv;
    for (std::size_t p=0; p<wveldiv.size(); ++p) wveldiv[p] *= m_vol[p];

    // initialize Helmholtz decomposition linear solver
    m_bicg[ thisIndex ].ckLocal()->
      init( {}, wveldiv, pbc, ignorebc,
            CkCallback(CkIndex_ALE::applied(nullptr), thisProxy[thisIndex]) );

  } else {

    // continue
    applied();

  }
}

void
ALE::applied( [[maybe_unused]] CkDataMsg* msg )
// *****************************************************************************
// Solve mesh velocity linear solve for ALE mesh motion
// *****************************************************************************
{
  //if (msg != nullptr) {
  //  auto *norm = static_cast< tk::real * >( msg->getData() );
  //  std::cout << "applied: " << *norm << '\n';
  //}

  auto smoother = g_inputdeck.get< tag::ale, tag::smoother >();

  if (smoother == ctr::MeshVelocitySmootherType::LAPLACE) {

    m_bicg[ thisIndex ].ckLocal()->solve(
       g_inputdeck.get< tag::ale, tag::maxit >(),
       g_inputdeck.get< tag::ale, tag::tolerance >(),
       CkCallback(CkIndex_ALE::solved(nullptr), thisProxy[thisIndex]) );

  } else if (smoother == ctr::MeshVelocitySmootherType::HELMHOLTZ) {

    m_bicg[ thisIndex ].ckLocal()->solve(
       g_inputdeck.get< tag::ale, tag::maxit >(),
       g_inputdeck.get< tag::ale, tag::tolerance >(),
       CkCallback(CkIndex_ALE::helmholtz(nullptr), thisProxy[thisIndex]) );

  } else {

    solved();

  }
}

void
ALE::helmholtz( [[maybe_unused]] CkDataMsg* msg )
// *****************************************************************************
//  Compute the gradient of the scalar potential for ALE
// *****************************************************************************
{
  //if (msg != nullptr) {
  //  auto *norm = static_cast< tk::real * >( msg->getData() );
  //  std::cout << "solved: " << *norm << '\n';
  //}

  // compute gradient of scalar potential for ALE (own portion)
  m_gradpot = tk::grad( m_coord, m_inpoel,
                 m_bicg[ thisIndex ].ckLocal()->solution() );

  // communicate scalar potential sums to other chares on chare-boundary
  if (m_nodeCommMap.empty())
    compot_complete();
  else
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::array< tk::real, 3 > > v( n.size() );
      std::size_t j = 0;
      for (auto i : n) {
        auto lid = tk::cref_find( m_lid, i );
        v[j][0] = m_gradpot[0][lid];
        v[j][1] = m_gradpot[1][lid];
        v[j][2] = m_gradpot[2][lid];
        ++j;
      }
      thisProxy[c].compot( std::vector<std::size_t>(begin(n),end(n)), v );
    }
  ownpot_complete();
}

void
ALE::compot( const std::vector< std::size_t >& gid,
               const std::vector< std::array< tk::real, 3 > >& v )
// *****************************************************************************
//  Receive contributions to scalar potential gradient on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] v Partial contributions to chare-boundary nodes
// *****************************************************************************
{
  Assert( v.size() == gid.size(), "Size mismatch" );
  using tk::operator+=;
  for (std::size_t i=0; i<gid.size(); ++i) m_gradpotc[ gid[i] ] += v[i];

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_npot == m_nodeCommMap.size()) {
    m_npot = 0;
    compot_complete();
  }
}

void
ALE::gradpot()
// *****************************************************************************
// Finalize computing gradient of the scalar potential for ALE
// *****************************************************************************
{
  // combine own and communicated contributions to scalar potential gradient
  for (const auto& [g,v] : m_gradpotc) {
    auto lid = tk::cref_find( m_lid, g );
    m_gradpot[0][lid] += v[0];
    m_gradpot[1][lid] += v[1];
    m_gradpot[2][lid] += v[2];
  }
  // clear receive buffer
  tk::destroy(m_gradpotc);

  // finish computing the gradient dividing weak sum by the nodal volumes
  for (auto j : g_inputdeck.get< tag::ale, tag::mesh_motion >())
    for (std::size_t p=0; p<m_gradpot[j].size(); ++p)
      m_gradpot[j][p] /= m_vol[p];

  solved();
}

void
ALE::solved( [[maybe_unused]] CkDataMsg* msg )
// *****************************************************************************
//  Mesh smoother linear solver converged
// *****************************************************************************
{
  //if (msg != nullptr) {
  //  auto *normres = static_cast< tk::real * >( msg->getData() );
  //  std::cout << "solved: " << *normres << '\n';
  //}

  auto smoother = g_inputdeck.get< tag::ale, tag::smoother >();

  if (smoother == ctr::MeshVelocitySmootherType::LAPLACE) {

    // Read out linear solution
    auto w = m_bicg[ thisIndex ].ckLocal()->solution();

    // Assign mesh velocity from linear solution skipping dimensions that are
    // not allowed to move
    for (auto j : g_inputdeck.get< tag::ale, tag::mesh_motion >())
      for (std::size_t i=0; i<m_w.nunk(); ++i)
        m_w(i,j) = w[i*m_w.nprop()+j];

  } else if (smoother == ctr::MeshVelocitySmootherType::HELMHOLTZ) {

    auto a1 = g_inputdeck.get< tag::ale, tag::vortmult >();
    for (auto j : g_inputdeck.get< tag::ale, tag::mesh_motion >())
      for (std::size_t p=0; p<m_w.nunk(); ++p)
        m_w(p,j) += a1 * (m_gradpot[j][p] - m_w(p,j));

  }

  // continue to applying a mesh force to the mesh velocity
  startforce();
}

void
ALE::startforce()
// *****************************************************************************
//  Compute mesh force for the ALE mesh velocity
//! \details Compute mesh forces. See Sec.4 in Bakosi, Waltz, Morgan, Improved
//!   ALE mesh velocities for complex flows, International Journal for Numerical
//!   Methods in Fluids, 2017.
// *****************************************************************************
{
  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  // compute pseudo-pressure gradient for mesh force
  const auto& f = g_inputdeck.get< tag::ale, tag::meshforce >();
  m_wf.fill( 0.0 );
  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    // access node IDs
    const std::array< std::size_t, 4 >
      N{{m_inpoel[e*4+0], m_inpoel[e*4+1], m_inpoel[e*4+2], m_inpoel[e*4+3]}};
    // compute element Jacobi determinant
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );        // J = 6V
    auto J24 = J/24.0;
    Assert( J > 0, "Element Jacobian non-positive" );

    // shape function derivatives, nnode*ndim [4][3]
    std::array< std::array< tk::real, 3 >, 4 > grad;
    grad[1] = tk::crossdiv( ca, da, J );
    grad[2] = tk::crossdiv( da, ba, J );
    grad[3] = tk::crossdiv( ba, ca, J );
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

    // max sound speed across nodes of element
    tk::real c = 1.0e-3;        // floor on sound speed
    for (std::size_t a=0; a<4; ++a) c = std::max( c, m_soundspeed[N[a]] );

    // mesh force in nodes
    auto V = J/6.0;
    auto L = std::cbrt(V);  // element length scale, L=V^(1/3)
    auto dv = m_vol0[e] / V;
    std::array< tk::real, 4 > q{0,0,0,0};
    for (std::size_t a=0; a<4; ++a) {
      auto du = L * m_veldiv[N[a]];
      q[a] = - du*(f[0]*c + f[1]*std::abs(du) + f[2]*du*du)
             + f[3]*c*c*std::abs(dv-1.0);
    }

    // scatter add pseudo-pressure gradient
    for (std::size_t a=0; a<4; ++a)
      for (std::size_t b=0; b<4; ++b)
        for (std::size_t i=0; i<3; ++i)
          m_wf(N[a],i) += J24 * grad[b][i] * q[b];
  }

  // communicate mesh force sums to other chares on chare-boundary
  if (m_nodeCommMap.empty()) {
    comfor_complete();
  } else {
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::array< tk::real, 3 > > w( n.size() );
      std::size_t j = 0;
      for (auto i : n) {
        auto lid = tk::cref_find( m_lid, i );
        w[j][0] = m_wf(lid,0);
        w[j][1] = m_wf(lid,1);
        w[j][2] = m_wf(lid,2);
        ++j;
      }
      thisProxy[c].comfor(std::vector<std::size_t>(begin(n),end(n)), w);
    }
  }
  ownfor_complete();
}

void
ALE::comfor( const std::vector< std::size_t >& gid,
             const std::vector< std::array< tk::real, 3 > >& w )
// *****************************************************************************
//  Receive contributions to ALE mesh force on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] w Partial contributions to chare-boundary nodes
// *****************************************************************************
{
  Assert( w.size() == gid.size(), "Size mismatch" );
  using tk::operator+=;
  for (std::size_t i=0; i<gid.size(); ++i) m_wfc[ gid[i] ] += w[i];

  // When we have heard from all chares we communicate with, this chare is done
  if (++m_nwf == m_nodeCommMap.size()) {
    m_nwf = 0;
    comfor_complete();
  }
}

void
ALE::meshforce()
// *****************************************************************************
// Apply mesh force
// *****************************************************************************
{
  // combine own and communicated contributions to mesh force
  for (const auto& [g,w] : m_wfc) {
    auto lid = tk::cref_find( m_lid, g );
    m_wf(lid,0) += w[0];
    m_wf(lid,1) += w[1];
    m_wf(lid,2) += w[2];
  }
  // clear receive buffer
  tk::destroy(m_wfc);

  // finish computing the mesh force by dviding weak sum by the nodal volumes
  for (std::size_t j=0; j<3; ++j)
    for (std::size_t p=0; p<m_wf.nunk(); ++p)
      m_wf(p,j) /= m_vol[p];

  // advance mesh velocity in time due to pseudo-pressure gradient mesh force
  for (auto j : g_inputdeck.get< tag::ale, tag::mesh_motion >())
    for (std::size_t i=0; i<m_w.nunk(); ++i)
       // This is likely incorrect. It should be m_w = m_w0 + ...
       m_w(i,j) += m_adt * m_wf(i,j);

  // Enforce mesh velocity Dirichlet BCs where user specfied but did not
  // prescribe a move
  for (auto i : m_meshveldirbcnodes)
    if (not move(i)) m_w(i,0) = m_w(i,1) = m_w(i,2) = 0.0;

  // On meshvel symmetry BCs remove normal component of mesh velocity
  const auto& sbc = g_inputdeck.get< tag::ale, tag::symmetry >();
  for (auto p : m_meshvelsymbcnodes) {
    for (const auto& s : sbc) {
      auto j = m_bnorm.find(static_cast<int>(s));
      if (j != end(m_bnorm)) {
        auto i = j->second.find(p);
        if (i != end(j->second)) {
          std::array< tk::real, 3 >
            n{ i->second[0], i->second[1], i->second[2] },
            v{ m_w(p,0), m_w(p,1), m_w(p,2) };
          auto v_dot_n = tk::dot( v, n );
          // symbc: remove normal component of mesh velocity
          m_w(p,0) -= v_dot_n * n[0];
          m_w(p,1) -= v_dot_n * n[1];
          m_w(p,2) -= v_dot_n * n[2];
        }
      }
    }
  }

  // Activate SDAG wait for re-computing prerequisites for ALE
  thisProxy[ thisIndex ].wait4vel();
  thisProxy[ thisIndex ].wait4pot();
  thisProxy[ thisIndex ].wait4for();

  m_done.send();
}

#include "NoWarning/ale.def.h"
