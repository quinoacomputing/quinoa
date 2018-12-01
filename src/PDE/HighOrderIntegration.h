// *****************************************************************************
/*!
  \file      src/PDE/HighOrderIntegration.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for high order numerical integration
  \details   This file contains functions that provide Gauss quadrature
    coordinates and weights for and functionality for high order numerical
    integration.
*/
// *****************************************************************************
#ifndef HighOrderIntegration_h
#define HighOrderIntegration_h

#include <array>

#include "Types.h"

namespace inciter {

//! Five Gaussian quadrature points locations and weights for a tetrahedron
void
GaussQuadratureTet( std::array< std::array< tk::real, 5 >, 3 >& coordgp,
                    std::array< tk::real, 5 >& wgp );

//! Four Gaussian quadrature points locations and weights for a tetrahedron
void
GaussQuadratureTet( std::array< std::array< tk::real, 4 >, 3 >& coordgp,
                    std::array< tk::real, 4 >& wgp );

//! Three Gaussian quadrature points locations and weights for a triangle
void
GaussQuadratureTri( std::array< std::array< tk::real, 3 >, 2 >& coordgp,
                    std::array< tk::real, 3 >& wgp );

} // inciter::

#endif // HighOrderIntegration_h
