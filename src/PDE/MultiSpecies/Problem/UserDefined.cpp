// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Problem/UserDefined.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a Problem policy class for the compressible flow
    equations, defined in PDE/MultiSpecies/MultiSpecies.h. See PDE/MultiSpecies/Problem.h
    for general requirements on Problem policy classes for MultiSpecies.
*/
// *****************************************************************************

#include <limits>

#include "UserDefined.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiSpeciesProblemUserDefined;

tk::InitializeFn::result_type
MultiSpeciesProblemUserDefined::initialize( ncomp_t ncomp,
                                        const std::vector< EOS >& mat_blk,
                                        tk::real,
                                        tk::real,
                                        tk::real,
                                        tk::real )
// *****************************************************************************
//! Set initial conditions
//! \param[in] ncomp Number of scalar components in this PDE system
//! \return Values of all components
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  tk::InitializeFn::result_type s( ncomp, 0.0 );

  auto nspec = g_inputdeck.get< eq, tag::nspec >();

  // Set background ICs
  const auto& ic = g_inputdeck.get< tag::ic >();
  const auto& bgmassfrac = ic.get< tag::mass_fractions >();
  const auto& bgvelic = ic.get< tag::velocity >();
  const auto& bgpreic = ic.get< tag::pressure >();
  const auto& bgtempic = ic.get< tag::temperature >();

  if (bgtempic < 1e-12) Throw( "No background temperature IC" );
  if (bgpreic < 1e-12) Throw( "No background pressure IC" );

  auto alphamin = 1.0e-12;

  // initialize background species states
  auto Ys = bgmassfrac;
  tk::real total_al(0.0);
  for (std::size_t k=0; k<nspec; ++k) {
    Ys[k] = std::max(Ys[k], alphamin);
    total_al += Ys[k];
  }
  for (std::size_t k=0; k<nspec; ++k) Ys[k] /= total_al;

  tk::real u = bgvelic[0];
  tk::real v = bgvelic[1];
  tk::real w = bgvelic[2];

  // Initialize mixture
  Mixture mix(nspec);
  mix.set_massfrac(Ys, bgpreic, bgtempic, mat_blk);

  auto rb = mix.get_mix_density();
  for (std::size_t k=0; k<nspec; ++k) {
    // partial density
    s[multispecies::densityIdx(nspec,k)] = Ys[k] * rb;
  }

  // total specific energy
  s[multispecies::energyIdx(nspec,0)] = mix.totalenergy(rb,
    u, v, w, bgpreic, mat_blk);

  // bulk momentum
  s[multispecies::momentumIdx(nspec,0)] = rb * u;
  s[multispecies::momentumIdx(nspec,1)] = rb * v;
  s[multispecies::momentumIdx(nspec,2)] = rb * w;

  if (bgpreic< 1e-12 || bgtempic< 1e-12)
    Throw("User must specify background pressure and temperature in IC.");

  return s;
}
