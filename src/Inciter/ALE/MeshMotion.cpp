// *****************************************************************************
/*!
  \file      src/Inciter/ALE/MeshMotion.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh motion functionality for ALE
  \details   Functionality for arbitrary Lagrangian-Eulerian (ALE) mesh motion.
*/
// *****************************************************************************

#include <limits>
#include <iostream>     // NOT NEEDED

#include "MeshMotion.hpp"
#include "Mesh/Around.hpp"
#include "Vector.hpp"

namespace inciter {

void
sine( const tk::UnsMesh::Coords& coord, tk::Fields& w )
// *****************************************************************************
//  Prescribe mesh motion as a sine wave for testing purposes
//! \param[in] coord Mesh node coordinates
//! \param[in,out] w Mesh velocity assigned
// *****************************************************************************
{
  for (std::size_t i=0; i<w.nunk(); ++i)
    w(i,0,0) = std::pow( std::sin(coord[0][i]*M_PI), 2.0 );
}

void
fluid( const tk::UnsMesh::Coords& v, tk::Fields& w )
// *****************************************************************************
//! Prescribe mesh velocity as the fluid velocity
//! \param[in] v Fluid velocity
//! \param[in,out] w Mesh velocity assigned
// *****************************************************************************
{
  for (std::size_t j=0; j<3; ++j)
    for (std::size_t i=0; i<w.nunk(); ++i)
      w(i,j,0) = v[j][i];
}

void
meshvel( ctr::MeshVelocityType m,
         const tk::UnsMesh::Coords& coord,
         const tk::UnsMesh::Coords& v,
         tk::Fields& w )
// *****************************************************************************
//  Assign mesh velocity based on user config
//! \param[in] m Mesh velocity type
//! \param[in] coord Mesh node coordinates
//! \param[in] vel Fluid velocity
//! \param[in,out] w Mesh velocity assigned
// *****************************************************************************
{
  Assert( w.nprop() == 3, "The mesh velocity must have 3 scalar components" );

  if (m == ctr::MeshVelocityType::SINE)
    sine( coord, w );
  else if (m == ctr::MeshVelocityType::FLUID)
    fluid( v, w );
  else
    Throw( "Mesh velocity not implemented" );
}

void
vortscale( const std::array< std::vector< tk::real >, 3 >& coord,
           const std::vector< std::size_t >& inpoel,
           const std::pair< std::vector< std::size_t >,
                            std::vector< std::size_t > >& esup,
           const std::vector< tk::real >& vol,
           const tk::UnsMesh::Coords& vel,
           tk::real c1,
           tk::real c2,
           tk::Fields& w )
// *****************************************************************************
//  Scale the mesh velocity with a function of the fluid vorticity for ALE
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] esup Elements surrounding points
//! \param[in] vol Nodal volumes
//! \param[in] vel Fluid velocity in mesh points
//! \param[in] w Mesh velocity to scale
// *****************************************************************************
{
   // access node cooordinates
   const auto& x = coord[0];
   const auto& y = coord[1];
   const auto& z = coord[2];
   // access velocity components
   const auto& vx = vel[0];
   const auto& vy = vel[1];
   const auto& vz = vel[2];

   auto npoin = x.size();
   std::vector< tk::real > v( npoin, 0.0 );
   auto maxv = -std::numeric_limits< tk::real >::max();

   #pragma omp simd
   for (std::size_t p=0; p<npoin; ++p) {
     tk::real vort[3] = { 0.0, 0.0, 0.0 };
     for (auto e : tk::Around(esup,p)) {
       // access node IDs
       std::size_t N[4] =
         { inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };
       // compute element Jacobi determinant, J = 6V
       auto bax = x[N[1]]-x[N[0]];
       auto bay = y[N[1]]-y[N[0]];
       auto baz = z[N[1]]-z[N[0]];
       auto cax = x[N[2]]-x[N[0]];
       auto cay = y[N[2]]-y[N[0]];
       auto caz = z[N[2]]-z[N[0]];
       auto dax = x[N[3]]-x[N[0]];
       auto day = y[N[3]]-y[N[0]];
       auto daz = z[N[3]]-z[N[0]];
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
       for (std::size_t b=0; b<4; ++b) {
         tk::real curl[3];
         tk::cross( g[b][0], g[b][1], g[b][2],
                    vx[N[b]], vy[N[b]], vz[N[b]],
                    curl[0], curl[1], curl[2] );
 
// if (p==0)
//   std::cout << "coords: " << x[N[b]] << ',' << y[N[b]] << ',' << z[N[b]]
//             << ": e:" << e
//             << " g: " << g[b][0] << ',' << g[b][1] << ',' << g[b][2]
//             << " v: " << vx[N[b]] << ',' << vy[N[b]] << ',' << vz[N[b]] << ' '
//             << " c: " << curl[0] << ',' << curl[1] << ',' << curl[2]
//             << '\n';

          vort[0] += J24 * curl[0];
          vort[1] += J24 * curl[1];
          vort[2] += J24 * curl[2];

          //for (std::size_t i=0; i<3; ++i)
          //  for (std::size_t j=0; j<3; ++j)
          //    for (std::size_t k=0; k<3; ++k)
          //      vort[i] += J24 * tk::perm(i,j,k) * g[b][j] * vel[k][N[b]];
       }
     }
     // divide weak result by nodal volume
     vort[0] /= vol[p];
     vort[1] /= vol[p];
     vort[2] /= vol[p];
std::cout << p << ": " << vort[0] << ',' << vort[1] << ',' << vort[2] << '\n';
     v[p] = tk::length( vort[0], vort[1], vort[2] );
     maxv = std::max( maxv, v[p] );
   }
 
//   for (std::size_t p=0; p<npoin; ++p) {
//     std::cout << "coords: " << x[p] << ',' << y[p] << ',' << z[p] << ": "
//               << v[p] << '\n';
//   }
//std::cout << "maxv: " << maxv << '\n';

   if (maxv > 1.0e-8) {
     for (std::size_t j=0; j<3; ++j) {
       #pragma omp simd
       for (std::size_t i=0; i<npoin; ++i)
         w(i,j,0) *= c1 * std::max( 0.0, 1.0 - c2*v[i]/maxv );
     }
   }
}

} // inciter::
