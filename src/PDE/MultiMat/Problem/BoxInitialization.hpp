// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/BoxInitialization.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     User-defined box initialization
  \details   This file defines functions for initializing solutions for
    compressible multi-material equations inside the user-defined box.
*/
// *****************************************************************************
#ifndef MultiMatBoxInitialization_h
#define MultiMatBoxInitialization_h

#include "Fields.hpp"
#include "EoS/EOS.hpp"
#include "ContainerUtil.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

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
  auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
  auto alphamin = g_inputdeck.get< tag::multimat, tag::min_volumefrac >();

  const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

  const auto& initiate = b.template get< tag::initiate >();

  // get material id in box (offset by 1, since input deck uses 1-based ids)
  std::size_t boxmatid = b.template get< tag::materialid >() - 1;
  const auto& boxvel = b.template get< tag::velocity >();
  auto boxpre = b.template get< tag::pressure >();
  auto boxene = b.template get< tag::energy >();
  auto boxtemp = b.template get< tag::temperature >();
  auto boxmas = b.template get< tag::mass >();
  auto boxenc = b.template get< tag::energy_content >();

  // [I] Compute the states inside the IC box/block based on the type of user
  // input.

  // material volume fractions
  for (std::size_t k=0; k<nmat; ++k) {
    if (k == boxmatid) {
      s[volfracIdx(nmat,k)] = 1.0 - (static_cast< tk::real >(nmat-1))*alphamin;
    }
    else {
      s[volfracIdx(nmat,k)] = alphamin;
    }
  }
  // material states (density, pressure, velocity)
  tk::real u = 0.0, v = 0.0, w = 0.0, spi(0.0), pr(0.0), tmp(0.0), rbulk(0.0);
  std::vector< tk::real > rhok(nmat, 0.0);
  // 1. User-specified mass, specific energy (J/m^3) and volume of box
  if (boxmas > 0.0) {
    if (boxenc <= 1e-12) Throw( "Box energy content must be nonzero" );
    // determine density and energy of material in the box
    rhok[boxmatid] = boxmas / V_ex;
    spi = boxenc / rhok[boxmatid];

    // Determine pressure and temperature
    auto boxmat_vf = s[volfracIdx(nmat,boxmatid)];

    // Since initiate type 'linear' assigns the background IC values to all
    // nodes within a box at initialization (followed by adding a time-dependent
    // energy source term representing a propagating wave-front), the pressure
    // in the box needs to be set to background pressure.
    if (initiate == ctr::InitiateType::LINEAR && t < 1e-12) {
      if (boxmas <= 1e-12 || boxenc <= 1e-12 || bgpreic <= 1e-12 ||
        bgtempic <= 1e-12)
        Throw("Box mass, energy content, background pressure and background "
          "temperature must be specified for IC with linear propagating source");

      pr = bgpreic;
      auto te = mat_blk[boxmatid].compute< EOS::totalenergy >(
        boxmat_vf*rhok[boxmatid], u, v, w, boxmat_vf*pr, boxmat_vf );
      tmp = mat_blk[boxmatid].compute< EOS::temperature >(
        boxmat_vf*rhok[boxmatid], u, v, w, te, boxmat_vf );

      // if the above density and pressure lead to an invalid thermodynamic
      // state (negative temperature/energy), change temperature to background
      // temperature and use corresponding density.
      if (tmp < 0.0 || te < 0.0) {
        tmp = bgtempic;
        rhok[boxmatid] = mat_blk[boxmatid].compute< EOS::density >(pr, tmp);
        spi = boxenc / rhok[boxmatid];
      }
    }
    // For initiate type 'impulse', pressure and temperature are determined from
    // energy content that needs to be dumped into the box at IC.
    else if (initiate == ctr::InitiateType::IMPULSE) {
      pr = mat_blk[boxmatid].compute< EOS::pressure >(
        boxmat_vf*rhok[boxmatid], u, v, w, boxmat_vf*rhok[boxmatid]*spi,
        boxmat_vf, boxmatid );
      tmp = mat_blk[boxmatid].compute< EOS::temperature >(
        boxmat_vf*rhok[boxmatid], u, v, w, boxmat_vf*rhok[boxmatid]*spi,
        boxmat_vf );
    }
    else Throw( "IC box initiate type not implemented for multimat" );

    // find density of trace material quantities in the box based on pressure
    for (std::size_t k=0; k<nmat; ++k) {
      if (k != boxmatid) {
        rhok[k] = mat_blk[k].compute< EOS::density >(pr, tmp);
      }
    }
  }
  // 2. User-specified temperature, pressure and velocity in box
  else {
    for (std::size_t k=0; k<nmat; ++k) {
      rhok[k] = mat_blk[k].compute< EOS::density >(boxpre, boxtemp);
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
      Throw("IC-box with specified energy not set up for multimat");
    }
  }
  // bulk density
  for (std::size_t k=0; k<nmat; ++k) rbulk += s[volfracIdx(nmat,k)]*rhok[k];

  // [II] Finally initialize the solution vector
  for (std::size_t k=0; k<nmat; ++k) {
    // partial density
    s[densityIdx(nmat,k)] = s[volfracIdx(nmat,k)] * rhok[k];
    // total specific energy
    if (boxmas > 0.0 && k == boxmatid &&
      initiate == ctr::InitiateType::IMPULSE) {
      s[energyIdx(nmat,k)] = s[volfracIdx(nmat,k)] * rhok[k] * spi;
    }
    else {
      // TEMP: Eventually we would need to initialize gk from control file
      std::array< std::array< tk::real, 3 >, 3 > gk;
      if (solidx[k] > 0) {
        for (std::size_t i=0; i<3; ++i) {
          for (std::size_t j=0; j<3; ++j) {
            if (i==j) gk[i][j] = 1.0;
            else gk[i][j] = 0.0;
            s[deformIdx(nmat,solidx[k],i,j)] = gk[i][j];
          }
        }
      }
      else {
        gk = {{}};
      }
      s[energyIdx(nmat,k)] =
        mat_blk[k].compute< EOS::totalenergy >( s[volfracIdx(nmat,k)]*rhok[k],
        u, v, w, s[volfracIdx(nmat,k)]*pr, s[volfracIdx(nmat,k)], gk );
    }
  }
  // bulk momentum
  s[momentumIdx(nmat,0)] = rbulk * u;
  s[momentumIdx(nmat,1)] = rbulk * v;
  s[momentumIdx(nmat,2)] = rbulk * w;
}

} //inciter::

#endif // MultiMatBoxInitialization_h
