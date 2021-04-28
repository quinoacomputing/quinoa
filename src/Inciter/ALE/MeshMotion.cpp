// *****************************************************************************
/*!
  \file      src/Inciter/ALE/MeshMotion.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh motion options for ALE
  \details   Mesh motion options for ALE.
*/
// *****************************************************************************

#include "MeshMotion.hpp"

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

} // inciter::
