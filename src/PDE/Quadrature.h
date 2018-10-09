// *****************************************************************************
/*!
  \file      src/PDE/Quadrature.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Quadrature coordinates and weights for numerical integration
  \details   This file contains functions that provide Gauss quadrature
    coordinates and weights for numerical integration.
*/
// *****************************************************************************
#ifndef Quadrature_h
#define Quadrature_h

#include <array>

#include "Types.h"

namespace inciter {

//! Gaussian quadrature points locations and weights for a tetrahedron
void
GaussQuadratureTet( std::array< std::array< tk::real, 5 >, 3 >& coordgp,
                    std::array< tk::real, 5 >& wgp );


//! Gaussian quadrature points locations and weights for a triangle
void
GaussQuadratureTri( std::array< std::array< tk::real, 3 >, 2 >& coordgp,
                    std::array< tk::real, 3 >& wgp );

} // inciter::

#endif // Quadrature_h
