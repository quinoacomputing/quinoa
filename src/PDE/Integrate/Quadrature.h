// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Quadrature.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Quadrature coordinates and weights for numerical integration
  \details   This file contains functions that provide Gauss quadrature
    coordinates and weights for numerical integration.
*/
// *****************************************************************************
#ifndef Quadrature_h
#define Quadrature_h

#include <array>
#include <vector>

#include "Types.h"
#include "Exception.h"

namespace tk {

//! Initialization of number of gauss points for face integration
//! \param[in] ndof Number of degree of freedom
constexpr std::size_t NGfa( const std::size_t ndof ) {
  return ndof == 1 ? 1 :
         ndof == 4 ? 3 :
         ndof == 10 ? 6 :
         throw std::logic_error("ndof must be one of 1,4,10");
}

//! Fourteen Gaussian quadrature points locations and weights for a tetrahedron
void
GaussQuadratureTet( std::array< std::array< real, 14 >, 3 >& coordgp,
                    std::array< real, 14 >& wgp );

//! Eleven Gaussian quadrature points locations and weights for a tetrahedron
void
GaussQuadratureTet( std::array< std::array< real, 11 >, 3 >& coordgp,
                    std::array< real, 11 >& wgp );

//! Five Gaussian quadrature points locations and weights for a tetrahedron
void
GaussQuadratureTet( std::array< std::array< real, 5 >, 3 >& coordgp,
                    std::array< real, 5 >& wgp );

//! Four Gaussian quadrature points locations and weights for a tetrahedron
void
GaussQuadratureTet( std::array< std::array< real, 4 >, 3 >& coordgp,
                    std::array< real, 4 >& wgp );

//! Initialize Gaussian quadrature points locations and weights for a triangle
void
GaussQuadratureTri( std::size_t NG,
                    std::array< std::vector< real >, 2 >& coordgp,
                    std::vector< real >& wgp );
} // tk::

#endif // Quadrature_h
