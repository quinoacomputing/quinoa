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
#include <cmath>

#include "MMS.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "PDE/MultiSpecies/Mixture/Mixture.hpp"
namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiSpeciesProblemMMS;

tk::InitializeFn::result_type
MultiSpeciesProblemMMS::initialize( ncomp_t ncomp,
  const std::vector< EOS >& mat_blk,
  tk::real x, tk::real y, tk::real z,
  tk::real t)
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
//! \note The function signature must follow tk::InitializeFn
//! \details This function initializes the method of manufactured solutions
// *****************************************************************************
{
  // see also Control/Inciter/InputDeck/Grammar.hpp
  Assert( ncomp == 5, "Number of scalar components must be 5" );

  auto nspec = g_inputdeck.get< eq, tag::nspec >();
  Assert( nspec == 1, "Method of manufactured solutions currently supports one species" );

  std::vector< tk::real > s( ncomp, 0.0 );
 
  // define the manufactured solution for conserved quantities
  const auto a = 1.0 + 0.1 * std::sin(2*M_PI * (x + y + z) - 2.0 * t);
  const auto b = 1.0 + 0.1 * std::cos(2*M_PI * (x + y + z) - 2.0 * t);
  s[multispecies::densityIdx(nspec, 0)] = a;
  s[multispecies::momentumIdx(nspec, 0)] = a * a; 
  s[multispecies::momentumIdx(nspec, 1)] = a * b;
  s[multispecies::momentumIdx(nspec, 2)] = a * b;
  s[multispecies::energyIdx(nspec, 0)] = a * (b + 0.5 * (a * a + 2*(b * b)));

  return s;
}


std::vector< tk::real >
MultiSpeciesProblemMMS::analyticSolution(
  ncomp_t ncomp, const std::vector< EOS >& mat_blk,
  tk::real x, tk::real y, tk::real z, tk::real t)
// *****************************************************************************
//! Evaluate the method of manufactured solutions at (x,y,z,t)
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] mat_blk Equations of state for all species
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
// *****************************************************************************
{ 
  const auto nspec = g_inputdeck.get< eq, tag::nspec >();

  const auto temperature = 1.0 + 0.1 * std::cos(2*M_PI * (x + y + z) - 2.0 * t);

  auto s = initialize(ncomp, mat_blk, x, y, z, t);
  s.resize(ncomp + 1);
  s[ncomp + multispecies::temperatureIdx(nspec,0)] = temperature;

  return s; 
}

tk::SrcFn::result_type
MultiSpeciesProblemMMS::src(
  ncomp_t,
  const std::vector< EOS >&,
  tk::real x,
  tk::real y,
  tk::real z,
  tk::real t,
  std::vector< tk::real >& sv )
// *****************************************************************************
//! Evaluate the manufactured source for the 3D Danis–Yan traveling wave.
//! The source is added to the right-hand side of the conserved equations.
// *****************************************************************************
{
  // The volume integrator passes nmat as the first argument, not ncomp.
  const auto nspec =
    g_inputdeck.get< tag::multispecies, tag::nspec >();

  Assert( nspec == 1,
          "Danis2022 manufactured solution requires one species" );
  Assert( sv.size() == 5,
          "Danis2022 manufactured solution requires five components" );

  constexpr tk::real gamma = 1.4;
  constexpr tk::real mu = 1.0e-3;
  constexpr tk::real Pr = 0.71;

  const tk::real k = 2.0 * std::acos(-1.0);
  const tk::real theta = k * (x + y + z) - 2.0 * t;

  const auto st = std::sin(theta);
  const auto ct = std::cos(theta);

  // Profiles and derivatives with respect to theta:
  // rho = u = r; v = w = e = b.
  const auto r   = 1.0 + 0.1 * st;
  const auto rp  =       0.1 * ct;
  const auto rpp =      -0.1 * st;

  const auto b   = 1.0 + 0.1 * ct;
  const auto bp  =      -0.1 * st;
  const auto bpp =      -0.1 * ct;

  // Sum of velocity components and its phase derivatives.
  const auto D   = r   + 2.0 * b;
  const auto Dp  = rp  + 2.0 * bp;
  const auto Dpp = rpp + 2.0 * bpp;

  // Specific total energy H and total energy density E.
  const auto H  = b + 0.5 * (r*r + 2.0*b*b);
  const auto Hp = bp + r*rp + 2.0*b*bp;

  const auto E  = r * H;
  const auto Ep = rp * H + r * Hp;

  // Pressure and its phase derivative.
  const auto p  = (gamma - 1.0) * r * b;
  const auto pp = (gamma - 1.0) * (rp*b + r*bp);

  // Momentum densities and their phase derivatives.
  const auto mx  = r * r;
  const auto mxp = 2.0 * r * rp;

  const auto my  = r * b;
  const auto myp = rp*b + r*bp;

  // Viscous momentum-flux divergences.
  const auto visc_x =
    mu * k*k * (3.0*rpp + Dpp/3.0);

  const auto visc_yz =
    mu * k*k * (3.0*bpp + Dpp/3.0);

  // Viscous work and heat-conduction contribution to energy.
  const auto J =
    3.0 * (rp*rp + r*rpp + 2.0*bp*bp + 2.0*b*bpp)
    + (Dp*Dp + D*Dpp) / 3.0
    + 3.0 * gamma / Pr * bpp;

  const auto visc_E = mu * k*k * J;

  // S = Q_t + div(F_c) - div(F_v).
  sv[multispecies::densityIdx(nspec, 0)] =
    -2.0*rp + k * (rp*D + r*Dp);

  sv[multispecies::momentumIdx(nspec, 0)] =
    -2.0*mxp + k * (mxp*D + mx*Dp + pp) - visc_x;

  const auto source_yz =
    -2.0*myp + k * (myp*D + my*Dp + pp) - visc_yz;

  sv[multispecies::momentumIdx(nspec, 1)] = source_yz;
  sv[multispecies::momentumIdx(nspec, 2)] = source_yz;

  sv[multispecies::energyIdx(nspec, 0)] =
    -2.0*Ep
    + k * ((Ep + pp)*D + (E + p)*Dp)
    - visc_E;
}