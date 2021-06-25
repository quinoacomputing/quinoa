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
#include "Mesh/DerivedData.hpp"
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
  else if (m == ctr::MeshVelocityType::FLUID ||
           m == ctr::MeshVelocityType::LAGRANGE)
    w = v;
  else
    Throw( "Mesh velocity not implemented" );
}

void
vortscale( const std::array< std::vector< tk::real >, 3 >& coord,
           const std::vector< std::size_t >& inpoel,
           const std::vector< tk::real >& vol,
           const tk::UnsMesh::Coords& vel,
           tk::real c1,
           tk::real c2,
           tk::Fields& w )
// *****************************************************************************
//  Scale the mesh velocity with a function of the fluid vorticity for ALE
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] vol Nodal volumes
//! \param[in] vel Fluid velocity in mesh points
//! \param[in] w Mesh velocity to scale
//! \see J. Waltz, N.R. Morgan, T.R. Canfield, M.R.J. Charest, L.D. Risinger,
//!   J.G. Wohlbier, A three-dimensional finite element arbitrary
//!   Lagrangian–Eulerian method for shock hydrodynamics on unstructured grids,
//!   Computers & Fluids, 92: 172-187, 2014.
// *****************************************************************************
{
   // compute vorticity
   auto vort = tk::curl( coord, inpoel, vol, vel );

   auto npoin = coord[0].size();

   // compute max vorticity
   std::vector< tk::real > v( npoin, 0.0 );
   auto maxv = -std::numeric_limits< tk::real >::max();
   for (std::size_t p=0; p<npoin; ++p) {
     v[p] = tk::length( vort[0][p], vort[1][p], vort[2][p] );
     maxv = std::max( maxv, v[p] );
   }
//std::cout << "maxv: " << maxv << '\n';
 
   // scale mesh velocity with a function of the fluid vorticity
   if (maxv > 1.0e-8) {
     for (std::size_t j=0; j<3; ++j) {
       #pragma omp simd
       for (std::size_t p=0; p<npoin; ++p) {
// std::cout << v[p] << '/' << maxv << '>';
         w(p,j,0) *= c1 * std::max( 0.0, 1.0 - c2*v[p]/maxv );
// std::cout << w(p,j,0) << ' ';
       }
     }
   }
}

} // inciter::
