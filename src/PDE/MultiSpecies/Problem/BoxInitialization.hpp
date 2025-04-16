// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Problem/BoxInitialization.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     User-defined box initialization
  \details   This file defines functions for initializing solutions for
    compressible multi-species equations inside the user-defined box.
*/
// *****************************************************************************
#ifndef MultiSpeciesBoxInitialization_h
#define MultiSpeciesBoxInitialization_h

#include "Fields.hpp"
#include "EoS/EOS.hpp"
#include "ContainerUtil.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"
#include "MultiSpecies/Mixture/Mixture.hpp"

namespace inciter {

using ncomp_t = tk::ncomp_t;

template< class B >
void initializeBox( const std::vector< EOS >& mat_blk,
                    tk::real /*V_ex*/,
                    tk::real /*t*/,
                    const B& b,
                    std::vector< tk::real >& s )
// *****************************************************************************
// Set the solution in the user-defined IC box/mesh block
//! \tparam B IC-block type to operate, ctr::box, or ctr::meshblock
// //! \param[in] V_ex Exact box volume
// //! \param[in] t Physical time
//! \param[in] b IC box configuration to use
//! \param[in,out] s Solution vector that is set to box ICs
//! \details This function sets the fluid density and total specific energy
//!   within a box initial condition, configured by the user. If the user
//!   is specified a box where mass is specified, we also assume here that
//!   internal energy content (energy per unit volume) is also
//!   specified. Specific internal energy (energy per unit mass) is then
//!   computed here (and added to the kinetic energy) from the internal
//!   energy per unit volume by multiplying it with the total box volume
//!   and dividing it by the total mass of the material in the box.
//!   Example (SI) units of the quantities involved:
//!    * internal energy content (energy per unit volume): J/m^3
//!    * specific energy (internal energy per unit mass): J/kg
// *****************************************************************************
{

  auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

  const auto& initiate = b.template get< tag::initiate >();

  // get species id in box (offset by 1, since input deck uses 1-based ids)
  const auto& boxmassfrac = b.template get< tag::mass_fractions >();
  const auto& boxvel = b.template get< tag::velocity >();
  auto boxpre = b.template get< tag::pressure >();
  auto boxene = b.template get< tag::energy >();
  auto boxtemp = b.template get< tag::temperature >();
  auto boxmas = b.template get< tag::mass >();

  auto alphamin = 1.0e-12;

  // [I] Compute the states inside the IC box/block based on the type of user
  // input.

  // species volume fractions
  auto alphas = boxmassfrac;
  tk::real total_al(0.0);
  for (std::size_t k=0; k<nspec; ++k) {
    alphas[k] = std::max(alphas[k], alphamin);
    total_al += alphas[k];
  }
  for (std::size_t k=0; k<nspec; ++k) alphas[k] /= total_al;

  // material states (density, pressure, velocity)
  tk::real u = 0.0, v = 0.0, w = 0.0, spi(0.0), pr(0.0), rbulk(0.0);
  std::vector< tk::real > rhok(nspec, 0.0);

  // 1. User-specified mass, specific energy (J/m^3) and volume of box
  if (boxmas > 0.0 || initiate != ctr::InitiateType::IMPULSE) {
    Throw( "IC-box initialization type not supported in multispecies" );
  }
  // 2. User-specified temperature, pressure and velocity in box
  else {
    Mixture mix(nspec);
    mix.set_massfrac(alphas, boxpre, boxtemp, mat_blk); // Initialize
    rbulk = mix.get_mix_density();
    for (std::size_t k=0; k<nspec; ++k) {
      rhok[k] = alphas[k]*rbulk;
    }
    if (boxvel.size() == 3) {
      u = boxvel[0];
      v = boxvel[1];
      w = boxvel[2];
    }
    if (boxpre > 0.0) {
      pr = boxpre;
    }
    if (boxene > 0.0) {
      Throw("IC-box with specified energy not set up for multispecies");
    }

    spi = mix.totalenergy(rbulk, u, v, w, pr, mat_blk) / rbulk;
  }

  // [II] Finally initialize the solution vector
  // partial density
  for (std::size_t k=0; k<nspec; ++k) {
    s[multispecies::densityIdx(nspec,k)] = alphas[k] * rhok[k];
  }
  // total specific energy
  s[multispecies::energyIdx(nspec,0)] = rbulk * spi;
  // bulk momentum
  s[multispecies::momentumIdx(nspec,0)] = rbulk * u;
  s[multispecies::momentumIdx(nspec,1)] = rbulk * v;
  s[multispecies::momentumIdx(nspec,2)] = rbulk * w;
}

} //inciter::

#endif // MultiSpeciesBoxInitialization_h
