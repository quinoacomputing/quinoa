// *****************************************************************************
/*!
  \file      src/PDE/CGPDE.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions common to ALECG
  \details   Functions common to ALECG.
*/
// *****************************************************************************

#include <array>
#include <vector>
#include <unordered_map>

#include "Vector.hpp"
#include "DerivedData.hpp"
#include "Exception.hpp"
#include "Around.hpp"
#include "Fields.hpp"
#include "CGPDE.hpp"
#include "FunctionPrototypes.hpp"

namespace inciter {
namespace cg {

void
chbgrad( ncomp_t ncomp,
         ncomp_t offset,
         const std::array< std::vector< tk::real >, 3 >& coord,
         const std::vector< std::size_t >& inpoel,
         const std::vector< std::size_t >& bndel,
         const std::vector< std::size_t >& gid,
         const std::unordered_map< std::size_t, std::size_t >& bid,
         const std::tuple< std::vector< tk::real >,
                           std::vector< tk::real > >& stag,
         const tk::Fields& U,
         tk::ElemGradFn egrad,
         tk::Fields& G )
// *****************************************************************************
//  Compute nodal gradients of primitive variables for ALECG on chare boundary
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] bndel List of elements contributing to chare-boundary nodes
//! \param[in] gid Local->global node id map
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!    global node ids (key)
//! \param[in] stag Stagnation point BC configuration
//! \param[in] U Solution vector at recent time step
//! \param[in] egrad Function to compute element contribution to nodal gradients
//! \param[in,out] G Nodal gradients of primitive variables
// *****************************************************************************
{
  Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
          "vector at recent time step incorrect" );

  // compute gradients of primitive variables in points
  G.fill( 0.0 );

  for (auto e : bndel) {  // for elements contributing to chare boundary
    const auto [N,g,u,J] = egrad( ncomp, offset, e, coord, inpoel, stag, U );
    auto J24 = J/24.0;
    for (std::size_t a=0; a<4; ++a) {
      auto i = bid.find( gid[N[a]] );
      if (i != end(bid))    // only contribute to chare-boundary node
        for (std::size_t b=0; b<4; ++b)
          for (std::size_t j=0; j<3; ++j)
            for (std::size_t c=0; c<ncomp; ++c)
              G(i->second,c*3+j,offset) += J24 * g[b][j] * u[c][b];
    }
  }
}

tk::Fields
nodegrad( ncomp_t offset,
          const std::array< std::vector< tk::real >, 3 >& coord,
          const std::vector< std::size_t >& inpoel,
          const std::unordered_map< std::size_t, std::size_t >& lid,
          const std::unordered_map< std::size_t, std::size_t >& bid,
          const std::vector< tk::real >& vol,
          const std::tuple< std::vector< tk::real >,
                            std::vector< tk::real > >& stag,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& esup,
          const tk::Fields& U,
          const tk::Fields& G,
          tk::ElemGradFn egrad )
// *****************************************************************************
//  Compute/assemble nodal gradients of primitive variables for ALECG in all
//  points
//! \param[in] offset Offset this PDE system operates from
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] lid Global->local node ids
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!    global node ids (key)
//! \param[in] vol Nodal volumes
//! \param[in] stag Stagnation point BC configuration
//! \param[in] U Solution vector at recent time step
//! \param[in] G Nodal gradients of primitive variables in chare-boundary nodes
//! \param[in] egrad Function to compute element contribution to nodal gradients
//! \return Gradients of primitive variables in all mesh points
// *****************************************************************************
{
  // allocate storage for nodal gradients of primitive variables
  tk::Fields Grad( U.nunk(), 15 );
  Grad.fill( 0.0 );

  // access node cooordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  // compute gradients of primitive variables in points
  auto npoin = U.nunk();
  #pragma omp simd
  for (std::size_t p=0; p<npoin; ++p)
    for (auto e : tk::Around(esup,p)) {
      // access node IDs
      std::size_t N[4] =
        { inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };
      // compute element Jacobi determinant, J = 6V
      tk::real bax = x[N[1]]-x[N[0]];
      tk::real bay = y[N[1]]-y[N[0]];
      tk::real baz = z[N[1]]-z[N[0]];
      tk::real cax = x[N[2]]-x[N[0]];
      tk::real cay = y[N[2]]-y[N[0]];
      tk::real caz = z[N[2]]-z[N[0]];
      tk::real dax = x[N[3]]-x[N[0]];
      tk::real day = y[N[3]]-y[N[0]];
      tk::real daz = z[N[3]]-z[N[0]];
      auto J = tk::triple( bax, bay, baz, cax, cay, caz, dax, day, daz );
      auto J24 = J/24.0;
      // shape function derivatives, nnode*ndim [4][3]
      tk::real g[4][3];
      tk::crossdiv( cax, cay, caz, dax, day, daz, J,
                    g[1][0], g[1][1], g[1][2] );
      tk::crossdiv( dax, day, daz, bax, bay, baz, J,
                    g[2][0], g[2][1], g[2][2] );
      tk::crossdiv( bax, bay, baz, cax, cay, caz, J,
                    g[3][0], g[3][1], g[3][2] );
      for (std::size_t i=0; i<3; ++i)
        g[0][i] = -g[1][i] - g[2][i] - g[3][i];
      // compute primitive variables at nodes
      tk::real u[5];
      for (std::size_t b=0; b<4; ++b) {
        u[0] = U(N[b],0,offset);
        u[1] = U(N[b],1,offset)/u[0];
        u[2] = U(N[b],2,offset)/u[0];
        u[3] = U(N[b],3,offset)/u[0];
        u[4] = U(N[b],4,offset)/u[0] - 0.5*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);
        for (std::size_t c=0; c<5; ++c)
          for (std::size_t i=0; i<3; ++i)
            Grad(p,c*3+i,0) += J24 * g[b][i] * u[c];
      }
    }

  // put in nodal gradients of chare-boundary points
  for (const auto& [g,b] : bid) {
    auto i = tk::cref_find( lid, g );
    for (ncomp_t c=0; c<Grad.nprop(); ++c)
      Grad(i,c,0) = G(b,c,0);
  }

  // divide weak result in gradients by nodal volume
  for (std::size_t p=0; p<npoin; ++p)
    for (std::size_t c=0; c<15; ++c)
      Grad(p,c,0) /= vol[p];

  return Grad;
}

std::array< tk::real, 3 >
edfnorm( const tk::UnsMesh::Edge& edge,
         const std::array< std::vector< tk::real >, 3 >&  coord,
         const std::vector< std::size_t >& inpoel,
         const std::unordered_map< tk::UnsMesh::Edge,
                  std::vector< std::size_t >,
                  tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> >& esued )
// *****************************************************************************
//  Compute normal of dual-mesh associated to edge
//! \param[in] edge Edge whose dual-face normal to compute given by local ids
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] esued Elements surrounding edges
//! \return Dual-face normal for edge
// *****************************************************************************
{
  std::array< tk::real, 3 > n{ 0.0, 0.0, 0.0 };

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

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

std::vector< tk::real >
solinc( tk::ncomp_t system, tk::ncomp_t ncomp, tk::real x, tk::real y,
        tk::real z, tk::real t, tk::real dt, tk::SolutionFn solution )
// *****************************************************************************
// Evaluate the increment from t to t+dt of the analytical solution at (x,y,z)
// for all components
//! \param[in] system Equation system index, i.e., which equation system we
//!   operate on among the systems of PDEs configured by the user
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution increment starting from
//! \param[in] dt Time increment at which evaluate the solution increment to
//! \param[in] solution Function used to evaluate the solution
//! \return Increment in values of all components evaluated at (x,y,z,t+dt)
// *****************************************************************************
{
  int inbox = 0;
  auto st1 = solution( system, ncomp, x, y, z, t, inbox );
  auto st2 = solution( system, ncomp, x, y, z, t+dt, inbox );

  std::transform( begin(st1), end(st1), begin(st2), begin(st2),
                  []( tk::real s, tk::real& d ){ return d -= s; } );

  return st2;
}

} // cg::
} // inciter::
