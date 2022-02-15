// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/BoxInitialization.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     User-defined box initialization
  \details   This file defines functions for initializing solutions for
    compressible single-material equations inside the user-defined box.
*/
// *****************************************************************************

#include "Control/Inciter/Types.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "BoxInitialization.hpp"
#include "EoS/EoS.hpp"
#include "EoS/MatBlock.hpp"
#include "EoS/StiffenedGas.hpp"
#include "ContainerUtil.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

void initializeBox( std::size_t system,
                    tk::real VRatio,
                    tk::real,
                    const inciter::ctr::box& b,
                    std::vector< tk::real >& s )
// *****************************************************************************
// Set the solution in the user-defined IC box
//! \param[in] system Equation system index
//! \param[in] VRatio Ratio of exact box volume to discrete box volume
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
  auto nmat =
    g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[system];

  std::vector< tk::real > box
    { b.template get< tag::xmin >(), b.template get< tag::xmax >(),
      b.template get< tag::ymin >(), b.template get< tag::ymax >(),
      b.template get< tag::zmin >(), b.template get< tag::zmax >() };

  const auto& initiate = b.template get< tag::initiate >();
  auto inittype = initiate.template get< tag::init >();

  auto boxmatid = b.template get< tag::materialid >();
  const auto& boxvel = b.template get< tag::velocity >();
  auto boxpre = b.template get< tag::pressure >();
  auto boxene = b.template get< tag::energy >();
  auto boxtemp = b.template get< tag::temperature >();
  auto boxmas = b.template get< tag::mass >();
  auto boxenc = b.template get< tag::energy_content >();

  auto alphamin = 1.0e-12;

  // initialize box material volume fraction
  for (std::size_t k=0; k<nmat; ++k) {
    if (k == boxmatid-1) {
      s[volfracIdx(nmat,k)] = 1.0 - (static_cast< tk::real >(nmat-1))*alphamin;
    }
    else {
      s[volfracIdx(nmat,k)] = alphamin;
    }
  }

  // initialize box material states
  tk::real u = 0.0, v = 0.0, w = 0.0;
  std::vector< tk::real > rhok(nmat, 0.0);

  // Initiate type 'impulse' simply assigns the prescribed values to all
  // nodes within a box.
  if (inittype == ctr::InitiateType::IMPULSE) {
    if (boxmas > 0.0) {
      Assert( boxenc > 0.0, "Box energy content must be nonzero" );
      // determine density and energy of material in the box
      auto V_ex = (box[1]-box[0]) * (box[3]-box[2]) * (box[5]-box[4]);
      rhok[boxmatid-1] = boxmas / V_ex;
      auto spi = boxenc * VRatio / rhok[boxmatid-1];

      // based on the density and energy of the material, determine pressure
      // and temperature
      auto boxmat_vf = s[volfracIdx(nmat,boxmatid-1)];
//      auto pr_box = eos_pressure< tag::multimat >(system,
//        boxmat_vf*rhok[boxmatid-1], u, v, w, boxmat_vf*rhok[boxmatid-1]*spi,
//        boxmat_vf, boxmatid-1);
      auto pr_box = m_mat_blk[boxmatid-1]->eos_pressure(system,
        boxmat_vf*rhok[boxmatid-1], u, v, w, boxmat_vf*rhok[boxmatid-1]*spi,
        boxmat_vf, boxmatid-1);
      auto t_box = eos_temperature< tag::multimat >(system,
        boxmat_vf*rhok[boxmatid-1], u, v, w, boxmat_vf*rhok[boxmatid-1]*spi,
        boxmat_vf, boxmatid-1);

      // find density of trace material quantities in the box based on pressure
      for (std::size_t k=0; k<nmat; ++k) {
        if (k != boxmatid-1) {
          rhok[k] = eos_density< tag::multimat >(system, pr_box, t_box, k);
        }
      }

      // initialize box based on above
      auto rb(0.0);
      for (std::size_t k=0; k<nmat; ++k) {
        // partial density
        s[densityIdx(nmat,k)] = s[volfracIdx(nmat,k)] * rhok[k];
        // total specific energy
        if (k == boxmatid-1) {
          s[energyIdx(nmat,k)] = s[volfracIdx(nmat,k)] * rhok[k] * spi;
        }
        else {
          s[energyIdx(nmat,k)] = s[volfracIdx(nmat,k)] *
            eos_totalenergy< tag::multimat >(system, rhok[k], u, v, w, pr_box,
            k);
        }
        // bulk density
        rb += s[densityIdx(nmat,k)];
      }
      // bulk momentum
      s[momentumIdx(nmat,0)] = rb * u;
      s[momentumIdx(nmat,1)] = rb * v;
      s[momentumIdx(nmat,2)] = rb * w;
    } else {
      for (std::size_t k=0; k<nmat; ++k) {
        rhok[k] = eos_density< tag::multimat >(system, boxpre, boxtemp, k);
      }
      if (boxvel.size() == 3) {
        u = boxvel[0];
        v = boxvel[1];
        w = boxvel[2];
      }
      if (boxpre > 0.0) {
        auto rb(0.0);
        for (std::size_t k=0; k<nmat; ++k) {
          // partial density
          s[densityIdx(nmat,k)] = s[volfracIdx(nmat,k)] * rhok[k];
          // total specific energy
          s[energyIdx(nmat,k)] = s[volfracIdx(nmat,k)] *
            eos_totalenergy< tag::multimat >(system, rhok[k], u, v, w, boxpre,
            k);
          // bulk density
          rb += s[densityIdx(nmat,k)];
        }
        // bulk momentum
        s[momentumIdx(nmat,0)] = rb * u;
        s[momentumIdx(nmat,1)] = rb * v;
        s[momentumIdx(nmat,2)] = rb * w;
      }
      if (boxene > 0.0) {
        Throw("IC-box with specified energy not set up for multimat");
      }
    }
  }
  else Throw( "IC box initiate type not implemented" );
}

} //inciter::
