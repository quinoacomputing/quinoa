// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/RichtmyerMeshkov.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a Problem policy class for the multi-material
    compressible flow equations, defined in PDE/MultiMat/MultiMat.hpp. See
    PDE/MultiMat/Problem.hpp for general requirements on Problem policy classes
    for MultiMat.
*/
// *****************************************************************************
#include <iostream>
#include "MixingLayer.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "PDE/MultiSpecies/Mixture/Mixture.hpp"
namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiSpeciesProblemMixingLayer;

tk::InitializeFn::result_type
MultiSpeciesProblemMixingLayer::initialize( ncomp_t ncomp,
  const std::vector< EOS >& mat_blk,
  tk::real x, tk::real y, tk::real,
  tk::real )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
//! \note The function signature must follow tk::InitializeFn
//! \details This function only initializes the mixing layer
//!   problem, but does not actually give the analytical solution at time
//!   greater than 0.
// *****************************************************************************
{
  // see also Control/Inciter/InputDeck/Grammar.hpp
  Assert( ncomp == 5, "Number of scalar components must be 5" );

  auto nspec = g_inputdeck.get< eq, tag::nspec >();
  //auto alphamin = g_inputdeck.get< eq, tag::min_volumefrac >();

  std::vector< tk::real > s( ncomp, 0.0 );
  tk::real p, T, u, v, w;

  // impose analytical solution instead of initial condition

  //logistic curve fit of analytical solution
  auto U1 = 10;
  auto mu = 1.716e-5;
  auto rho0 = 1.225;
  auto nu = mu/rho0;
  auto f0 = 0.5;
  auto f1 = 1.0;
  auto c = 1.4498905297373783;
  auto eta0 = -0.10009499242392805;
  auto x0 = 200;
  auto eta = y * std::sqrt(U1 / (2*x0*nu));
  auto f = 0.0;
  if (eta < -5) {
    f = f0;
  }
  else if (eta > 5) {
    f = f1;
  }
  else {
    f = f0 + (f1 - f0) / (1 + std::exp(-c * (eta - eta0)));
  }

  auto denom = (1 + std::exp(-c * std::abs(eta - eta0)));
  p = 101325;
  T = 273.15;
  u = f*U1;
  v = 0;
  w = 0;
  auto rhob = 0.0;
  // density
  for (std::size_t k = 0; k<nspec; ++k) {
    auto rho = mat_blk[k].compute< EOS::density >(p, T);
    rhob += rho;
    s[multispecies::densityIdx(nspec,k)] = rho; 
    s[multispecies::energyIdx(nspec,k)] =
      mat_blk[k].compute< EOS::totalenergy >( rho, u, v, w, p, 1.0);
 }

  // momentum
  s[multispecies::momentumIdx(nspec, 0)] = rhob*u;
  s[multispecies::momentumIdx(nspec, 1)] = rhob*v;
  s[multispecies::momentumIdx(nspec, 2)] = rhob*w;
 
  return s;
}
