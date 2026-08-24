// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Problem/CouetteFlow.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for Couette flow
*/
// *****************************************************************************
#include <algorithm>

#include "CouetteFlow.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::MultiSpeciesProblemCouetteFlow;

namespace {

tk::InitializeFn::result_type
couetteState( tk::ncomp_t ncomp, const std::vector< inciter::EOS >& mat_blk,
              tk::real u )
// *****************************************************************************
//! Construct a conservative multi-species state with streamwise velocity u
// *****************************************************************************
{
  Assert( ncomp == 5, "Number of scalar components must be 5" );

  const auto nspec = inciter::g_inputdeck.get<
    tag::multispecies, tag::nspec >();

  constexpr tk::real p = 100000.0;
  constexpr tk::real T = 400.0;
  constexpr tk::real v = 0.0;
  constexpr tk::real w = 0.0;

  std::vector< tk::real > s( ncomp+1, 0.0 );
  tk::real rhob = 0.0;

  for (std::size_t k=0; k<nspec; ++k) {
    const auto rho =
      mat_blk[k].compute< inciter::EOS::density >( p, T );
    rhob += rho;
    s[inciter::multispecies::densityIdx(nspec,k)] = rho;
    s[inciter::multispecies::energyIdx(nspec,k)] =
      mat_blk[k].compute< inciter::EOS::totalenergy >(
        rho, u, v, w, p, 1.0 );
    s[ncomp + inciter::multispecies::temperatureIdx(nspec,k)] = T;
  }

  s[inciter::multispecies::momentumIdx(nspec,0)] = rhob*u;
  s[inciter::multispecies::momentumIdx(nspec,1)] = rhob*v;
  s[inciter::multispecies::momentumIdx(nspec,2)] = rhob*w;

  return s;
}

} // anonymous namespace

tk::InitializeFn::result_type
MultiSpeciesProblemCouetteFlow::initialize(
  ncomp_t ncomp, const std::vector< EOS >& mat_blk,
  tk::real, tk::real, tk::real, tk::real )
// *****************************************************************************
//! Initialize the domain with quiescent fluid
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] mat_blk Equations of state for all species
//! \return Values of all components evaluated at (x,y,z,t)
//! \note The function signature must follow tk::InitializeFn.
//! \details The initially stationary fluid evolves toward Couette flow under
//!   the boundary conditions returned by analyticSolution().
// *****************************************************************************
{
  return couetteState( ncomp, mat_blk, 0.0 );
}

std::vector< tk::real >
MultiSpeciesProblemCouetteFlow::analyticSolution(
  ncomp_t ncomp, const std::vector< EOS >& mat_blk,
  tk::real, tk::real y, tk::real, tk::real )
// *****************************************************************************
//! Evaluate the Couette-flow boundary solution at (x,y,z,t)
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] mat_blk Equations of state for all species
//! \param[in] y Vertical coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
//! \details The domain dimensions match the mixing-layer mesh: x=[0,0.5],
//!   y=[-0.1,0.1], and z=[-0.01,0]. The streamwise velocity varies linearly
//!   from zero at the bottom wall to 10 m/s at the top wall. This solution is
//!   also imposed on every surface selected as Dirichlet in the control file.
// *****************************************************************************
{
  constexpr tk::real ymin = -0.1;
  constexpr tk::real ymax = 0.1;
  constexpr tk::real umax = 10.0;

  const auto vertical_fraction =
    std::clamp( (y-ymin)/(ymax-ymin), 0.0, 1.0 );
  const auto u = umax * vertical_fraction;

  return couetteState( ncomp, mat_blk, u );
}
