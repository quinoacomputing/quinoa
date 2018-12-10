// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Quadrature.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Quadrature coordinates and weights for numerical integration
  \details   This file contains functions that provide Gauss quadrature
    coordinates and weights for numerical integration.
*/
// *****************************************************************************
#ifndef Quadrature_h
#define Quadrature_h

#include <array>

#include "Types.h"

namespace tk {

//! Five Gaussian quadrature points locations and weights for a tetrahedron
void
GaussQuadratureTet( std::array< std::array< real, 5 >, 3 >& coordgp,
                    std::array< real, 5 >& wgp );

//! Four Gaussian quadrature points locations and weights for a tetrahedron
void
GaussQuadratureTet( std::array< std::array< real, 4 >, 3 >& coordgp,
                    std::array< real, 4 >& wgp );

//! Three Gaussian quadrature points locations and weights for a triangle
void
GaussQuadratureTri( std::array< std::array< real, 3 >, 2 >& coordgp,
                    std::array< real, 3 >& wgp );

} // tk::

#endif // Quadrature_h