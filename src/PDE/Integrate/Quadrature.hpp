// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Quadrature.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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
#include <stdexcept>

#include "Types.hpp"
#include "Exception.hpp"

namespace tk {
//! Initialization of number of Gauss points for volume integration
//! \param[in] ndof Number of degrees of freedom
KOKKOS_INLINE_FUNCTION
constexpr std::size_t NGvol( const std::size_t ndof ) {
  return ndof == 1 ? 1 :
         ndof == 4 ? 5 :
         ndof == 10 ? 11 :
         throw std::logic_error("ndof must be one of 1,4,10");
}


//! Initialization of number of Gauss points for face integration
//! \param[in] ndof Number of degrees of freedom
constexpr std::size_t NGfa( const std::size_t ndof ) {
  return ndof == 1 ? 1 :
         ndof == 4 ? 3 :
         ndof == 10 ? 6 :
         throw std::logic_error("ndof must be one of 1,4,10");
}

//! \brief Initialization of number of Gauss points for volume integration
//!   in error estimation
//! \param[in] ndof Number of degrees of freedom
constexpr std::size_t NGdiag( const std::size_t ndof ) {
  return ndof == 1 ? 1 :
         ndof == 4 ? 4 :
         ndof == 10 ? 14 :
         throw std::logic_error("ndof must be one of 1,4,10");
}

//! \brief Initialization of number of Gauss points for volume integration
//!   in DG initialization
//! \param[in] ndof Number of degrees of freedom
constexpr std::size_t NGinit( const std::size_t ndof ) {
  return ndof == 1 ? 1 :
         ndof == 4 ? 14 :
         ndof == 10 ? 14 :
         throw std::logic_error("ndof must be one of 1,4,10");
}

//! Initialize Gaussian quadrature points locations and weights for a tetrahedron
void
GaussQuadratureTet( std::size_t NG,
                    std::array< std::vector< real >, 3 >& coordgp,
                    std::vector< real >& wgp );

//! Kokkos version of GaussQuadratureTet
KOKKOS_INLINE_FUNCTION 
void tk::GaussQuadratureTet( const std::size_t NG,
              Kokkos::View<real**, memory_space> coordgp,
              Kokkos::View<real**, memory_space> wgp);

//! Initialize Gaussian quadrature points locations and weights for a triangle
void
GaussQuadratureTri( std::size_t NG,
                    std::array< std::vector< real >, 2 >& coordgp,
                    std::vector< real >& wgp );
} // tk::

#endif // Quadrature_h
