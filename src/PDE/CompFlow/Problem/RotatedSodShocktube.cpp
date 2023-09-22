/*!
  \file      src/PDE/CompFlow/Problem/RotatedSodShocktube.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a Problem policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problem.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************

#include "Vector.hpp"
#include "RotatedSodShocktube.hpp"
#include "EoS/EOS.hpp"

using inciter::CompFlowProblemRotatedSodShocktube;

tk::InitializeFn::result_type
CompFlowProblemRotatedSodShocktube::initialize(
  ncomp_t ncomp,
  const std::vector< EOS >& mat_blk,
  tk::real x,
  tk::real y,
  tk::real z,
  tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time at which to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  // Assume the domain is rotated by 45 degrees about the X, Y, and then Z
  // axis compared to the original tube with largest dimension in X
  tk::real a = -45.0*M_PI/180.0;
  auto c = tk::rotateX( tk::rotateY( tk::rotateZ( {{ x, y, z }}, a ), a ), a );
  return CompFlowProblemSodShocktube::
           initialize( ncomp, mat_blk, c[0], c[1], c[2], t );
}
