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
//! \param[in] v Fluid velocity
//! \param[in,out] w Mesh velocity assigned
// *****************************************************************************
{
  Assert( w.nprop() == 3, "The mesh velocity must have 3 scalar components" );

  if (m == ctr::MeshVelocityType::SINE)
    sine( coord, w );
  else if (m == ctr::MeshVelocityType::FLUID ||
           m == ctr::MeshVelocityType::LAGRANGE ||
           m == ctr::MeshVelocityType::HELMHOLTZ)
    w = v;
  else
    Throw( "Mesh velocity not implemented" );
}

void
vortscale( const std::vector< tk::real >& vort,
           tk::real vmult,
           tk::real maxv,
           tk::Fields& w )
// *****************************************************************************
//  Scale the mesh velocity with a function of the fluid vorticity for ALE
//! \param[in] vort Vorticity magnitude in mesh points
//! \param[in] vmult Vorticity multiplier
//! \param[in] maxv Largest vorticity magnitude
//! \param[in] w Mesh velocity to scale
//! \see J. Waltz, N.R. Morgan, T.R. Canfield, M.R.J. Charest, L.D. Risinger,
//!   J.G. Wohlbier, A three-dimensional finite element arbitrary
//!   Lagrangian–Eulerian method for shock hydrodynamics on unstructured grids,
//!   Computers & Fluids, 92: 172-187, 2014.
// *****************************************************************************
{
  Assert( w.nunk() == vort.size(), "Size mismatch" );

  // scale mesh velocity with a function of the fluid vorticity
  if (maxv > 1.0e-8)
    for (std::size_t j=0; j<3; ++j)
      for (std::size_t p=0; p<vort.size(); ++p)
        w(p,j,0) *= std::max( 0.0, 1.0 - vmult*vort[p]/maxv );
}

} // inciter::
