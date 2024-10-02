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

namespace inciter {

using ncomp_t = tk::ncomp_t;

template< class B >
void initializeBox( const std::vector< EOS >& mat_blk,
                    tk::real V_ex,
                    tk::real t,
                    const B& b,
                    tk::real bgpreic,
                    tk::real bgtempic,
                    std::vector< tk::real >& s )
// *****************************************************************************
// Set the solution in the user-defined IC box/mesh block
//! \tparam B IC-block type to operate, ctr::box, or ctr::meshblock
//! \param[in] V_ex Exact box volume
//! \param[in] t Physical time
//! \param[in] b IC box configuration to use
//! \param[in] bgpreic Background pressure user input
//! \param[in] bgtempic Background temperature user input
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
  std::size_t boxspecid = b.template get< tag::materialid >() - 1;
  const auto& boxvel = b.template get< tag::velocity >();
  auto boxpre = b.template get< tag::pressure >();
  auto boxene = b.template get< tag::energy >();
  auto boxtemp = b.template get< tag::temperature >();
  auto boxmas = b.template get< tag::mass >();
  auto boxenc = b.template get< tag::energy_content >();

  auto alphamin = 1.0e-12;

  // [I] Compute the states inside the IC box/block based on the type of user
  // input.

  // species volume fractions
  std::vector< tk::real > alphas(nspec, alphamin);

  for (std::size_t k=0; k<nspec; ++k) {
    if (k == boxspecid) {
      alphas[k] = 1.0 - (static_cast< tk::real >(nspec-1))*alphamin;
    }
  }
  // material states (density, pressure, velocity)
  tk::real u = 0.0, v = 0.0, w = 0.0, spi(0.0), pr(0.0), tmp(0.0), rbulk(0.0);
  std::vector< tk::real > rhok(nspec, 0.0);

  // 1. User-specified mass, specific energy (J/m^3) and volume of box
  if (boxmas > 0.0) {
    if (boxenc <= 1e-12) Throw( "Box energy content must be nonzero" );
    // determine density and energy of species in the box
    rhok[boxspecid] = boxmas / V_ex;
    spi = boxenc / rhok[boxspecid];

    // Determine pressure and temperature

    // For initiate type 'impulse', pressure and temperature are determined from
    // energy content that needs to be dumped into the box at IC.
    if (initiate == ctr::InitiateType::IMPULSE) {
      pr = mat_blk[0].compute< EOS::pressure >(
        rhok[boxspecid], u, v, w, rhok[boxspecid]*spi );
      tmp = mat_blk[0].compute< EOS::temperature >(
        rhok[boxspecid], u, v, w, rhok[boxspecid]*spi );
    }
    else Throw( "IC box initiate type not implemented for multispecies" );

    // find density of trace species quantities in the box based on pressure
    for (std::size_t k=0; k<nspec; ++k) {
      if (k != boxspecid) {
        rhok[k] = mat_blk[0].compute< EOS::density >(pr, tmp);
      }
    }
  }
  // 2. User-specified temperature, pressure and velocity in box
  else {
    for (std::size_t k=0; k<nspec; ++k) {
      rhok[k] = mat_blk[0].compute< EOS::density >(boxpre, boxtemp);
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
      Throw("IC-box with specified energy not set up for multispec");
    }
  }
  // bulk density
  for (std::size_t k=0; k<nspec; ++k) rbulk += alphas[k]*rhok[k];

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
