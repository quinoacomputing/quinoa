// *****************************************************************************
/*!
  \file      src/Inciter/AMR/Error.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Class for computing error estimates for mesh refinement
  \details   Class for computing error estimates for mesh refinement.
*/
// *****************************************************************************

#include <cmath>
#include <limits>

#include "Exception.hpp"
#include "Error.hpp"
#include "Vector.hpp"
#include "Gradients.hpp"
#include "Inciter/Options/AMRError.hpp"

using AMR::Error;

tk::real
Error::scalar( const tk::Fields& u,
               const edge_t& edge,
               ncomp_t c,
               const std::array< std::vector< tk::real >, 3 >& coord,
               const std::vector< std::size_t >& inpoel,
               const std::pair< std::vector< std::size_t >,
                                std::vector< std::size_t > >& esup,
               inciter::ctr::AMRErrorType err ) const
// *****************************************************************************
//  Estimate error for scalar quantity
//! \param[in] u Solution vector
//! \param[in] edge Edge defined by its two end-point IDs
//! \param[in] c Scalar component to compute error of
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] esup Linked lists storing elements surrounding points, see
//!    tk::genEsup()
//! \param[in] err AMR Error indicator type
//! \return Error indicator: a real number between [0...1] inclusive
// *****************************************************************************
{
  if (err == inciter::ctr::AMRErrorType::JUMP)
    return error_jump( u, edge, c );
  else if (err == inciter::ctr::AMRErrorType::HESSIAN)
    return error_hessian( u, edge, c, coord, inpoel, esup );
  else
    Throw( "No such AMR error indicator type" );
}

tk::real
Error::error_jump( const tk::Fields& u,
                   const edge_t& edge,
                   ncomp_t c ) const
// *****************************************************************************
//  Estimate error for scalar quantity on edge based on jump in solution
//! \param[in] u Solution vector
//! \param[in] edge Edge defined by its two end-point IDs
//! \param[in] c Scalar component to compute error of
//! \return Error indicator: a real number between [0...1] inclusive
// *****************************************************************************
{
  auto a = edge.first();
  auto b = edge.second();

  // Normalization factor with a noise filter
  auto norm = std::abs(u(a,c,0)) + std::abs(u(b,c,0)) + 1e-3;

  // Normalized error
  return std::abs( u(a,c,0) - u(b,c,0) ) / norm;
}

tk::real
Error::error_hessian( const tk::Fields& u,
                      const edge_t& edge,
                      ncomp_t c,
                      const std::array< std::vector< tk::real >, 3 >& coord,
                      const std::vector< std::size_t >& inpoel,
                      const std::pair< std::vector< std::size_t >,
                                       std::vector< std::size_t > >& esup )
const
// *****************************************************************************
//  Estimate error for scalar quantity on edge based on Hessian of solution
//! \param[in] u Solution vector
//! \param[in] edge Edge defined by its two end-point IDs
//! \param[in] c Scalar component to compute error of
//! \param[in] coord Mesh node coordinates
//! \param[in] inpoel Mesh element connectivity
//! \param[in] esup Linked lists storing elements surrounding points, see
//!    tk::genEsup()
//! \return Error indicator: a real number between [0...1] inclusive
// *****************************************************************************
{
  const tk::real small = std::numeric_limits< tk::real >::epsilon();

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];
  auto a = edge.first();
  auto b = edge.second();

  // Compute edge vector
  std::array< tk::real, 3 > h {{ x[a]-x[b], y[a]-y[b], z[a]-z[b] }};

  // Compute gradients at edge-end points
  auto ga = nodegrad( a, coord, inpoel, esup, u, c );
  auto gb = nodegrad( b, coord, inpoel, esup, u, c );

  // Compute dot products of gradients and edge vectors
  auto dua = tk::dot( ga, h );
  auto dub = tk::dot( gb, h );

  // If the normalization factor is zero, return zero error
  auto norm = std::abs(dua) + std::abs(dub);
  if (norm < small) return 0.0;

  return std::abs(dub-dua) / norm;
}
