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

// T = 400 K is temp of Mixing Layer used in all following methods.
constexpr tk::real T = 400.0;

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
  Assert( nspec == 1, "Mixing layer currently supports one species" );

  std::vector< tk::real > s( ncomp, 0.0 );
  tk::real u, v, w;

  // impose analytical solution instead of initial condition

  // Tabulated cubic-spline interpolation of the similarity solution
  const auto p = 100000.0;
  const auto U1 = 10.0;
  const auto mu = 1.47e-5;

  const std::vector< tk::real > Ys{ 1.0 };
  inciter::Mixture mix( nspec, Ys, p, T, mat_blk );
  tk::real rhob = mix.get_mix_density();

  const auto nu = mu/rhob;
  const auto x0 = 200.0; // choose far downstream for thick mixing layer
  const auto eta = y * std::sqrt(U1 / (2*x0*nu));
  const auto f = mixingLayerVelocityProfile( eta );
  u = f*U1;
  v = 0;
  w = 0;

  // Partial species densities
  for (std::size_t k = 0; k<nspec; ++k) {
    s[multispecies::densityIdx(nspec,k)] = Ys[k] * rhob;
  }

  // Mixture total energy
  s[multispecies::energyIdx(nspec,0)] =
    mix.totalenergy( rhob, u, v, w, T, mat_blk );

  // momentum
  s[multispecies::momentumIdx(nspec, 0)] = rhob*u;
  s[multispecies::momentumIdx(nspec, 1)] = rhob*v;
  s[multispecies::momentumIdx(nspec, 2)] = rhob*w;

 
  return s;
}


std::vector< tk::real >
MultiSpeciesProblemMixingLayer::analyticSolution(
  ncomp_t ncomp, const std::vector< EOS >& mat_blk,
  tk::real x, tk::real y, tk::real z, tk::real t)
// *****************************************************************************
//! Evaluate the mixing layer boundary solution at (x,y,z,t)
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] mat_blk Equations of state for all species
//! \param[in] y Vertical coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
// *****************************************************************************
{ 
  const auto nspec = g_inputdeck.get< eq, tag::nspec >();

  auto s = initialize(ncomp, mat_blk, x, y, z, t);
  s.resize(ncomp + 1);
  s[ncomp + multispecies::temperatureIdx(nspec,0)] = T;

  return s; 
}
