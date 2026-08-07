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
#include <algorithm>

#include "MixingLayer.hpp"
#include "MixingLayerSpline.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "PDE/MultiSpecies/Mixture/Mixture.hpp"
namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiSpeciesProblemMixingLayer;

namespace {

tk::real mixingLayerVelocityProfile( tk::real eta )
// ****************************************************************************
//! Evaluate u/U1 using the tabulated similarity-solution cubic spline
// ****************************************************************************
{
  const auto& spline = inciter::mixing_layer::velocity_spline;

  const auto evaluate = []( const auto& segment, tk::real location ) {
    const auto deta = location - segment.eta_left;
    return ((segment.c3*deta + segment.c2)*deta + segment.c1)*deta +
      segment.c0;
  };

  if (eta <= spline.front().eta_left) {
    return spline.front().c0;
  }

  if (eta >= spline.back().eta_right) {
    return evaluate( spline.back(), spline.back().eta_right );
  }

  const auto segment = std::upper_bound(
    begin(spline), end(spline), eta,
    []( tk::real location, const auto& interval ) {
      return location < interval.eta_right;
    } );

  return evaluate( *segment, eta );
}

} // anonymous namespace

tk::InitializeFn::result_type
MultiSpeciesProblemMixingLayer::initialize( ncomp_t ncomp,
  const std::vector< EOS >& mat_blk,
  tk::real x, tk::real y, tk::real,
  tk::real)
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
//! \note The function signature must follow tk::InitializeFn
//! \details This function initializes the mixing layer problem using the same
//!  smooth velocity profile imposed at the boundaries.
// *****************************************************************************
{
  // see also Control/Inciter/InputDeck/Grammar.hpp
  Assert( ncomp == 5, "Number of scalar components must be 5" );

  auto nspec = g_inputdeck.get< eq, tag::nspec >();

  std::vector< tk::real > s( ncomp+1, 0.0 );
  tk::real u, v, w;

  // impose analytical solution instead of initial condition

  // Tabulated cubic-spline interpolation of the similarity solution
  const auto p = 100000.0;
  const auto T = 400.0;
  const auto U1 = 10.0;
  const auto mu = 1.47e-5;
  const auto rho0 =
    mat_blk[0].compute< EOS::density >( p, T );
  const auto nu = mu/rho0;
  const auto x0 = 200.0; // choose far downstream for thick mixing layer
  const auto eta = y * std::sqrt(U1 / (2*x0*nu));
  const auto f = mixingLayerVelocityProfile( eta );
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
    s[ncomp + multispecies::temperatureIdx(nspec,k)] = T;
 }

  // momentum
  s[multispecies::momentumIdx(nspec, 0)] = rhob*u;
  s[multispecies::momentumIdx(nspec, 1)] = rhob*v;
  s[multispecies::momentumIdx(nspec, 2)] = rhob*w;

 
  return s;
}
