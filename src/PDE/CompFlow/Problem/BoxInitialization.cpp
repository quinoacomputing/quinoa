// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/BoxInitialization.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     User-defined box initialization
  \details   This file defines functions for initializing solutions for
    compressible single-material equations inside the user-defined box.
*/
// *****************************************************************************

#include "BoxInitialization.hpp"
#include "EoS/EoS.hpp"
#include "ContainerUtil.hpp"
#include "Control/Inciter/Types.hpp"

namespace inciter {

void initializeBox( std::size_t system,
  tk::real VRatio,
  tk::real t,
  const inciter::ctr::box& icbox,
  const std::vector< std::vector< tk::real > >& bgpreic,
  const std::vector< std::vector< tk::real > >& cv,
  std::vector< tk::real >& s )
// *****************************************************************************
// Set the solution in the user-defined IC box
//! \param[in] VRatio Ratio of exact box volume to discrete box volume
//! \param[in] t Physical time
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
  // Read all box IC related parameters from inputdeck
  const auto& initiate = icbox.get< tag::initiate >();
  const auto& inittype = initiate.get< tag::init >();

  const auto& boxrho = icbox.get< tag::density >();
  const auto& boxvel = icbox.get< tag::velocity >();
  const auto& boxpre = icbox.get< tag::pressure >();
  const auto& boxene = icbox.get< tag::energy >();
  const auto& boxtem = icbox.get< tag::temperature >();
  const auto& boxmas = icbox.get< tag::mass >();
  const auto& boxenc = icbox.get< tag::energy_content >();
  std::array< tk::real, 6 >
    boxdim{ icbox.get< tag::xmin >(), icbox.get< tag::xmax >(),
            icbox.get< tag::ymin >(), icbox.get< tag::ymax >(),
            icbox.get< tag::zmin >(), icbox.get< tag::zmax >() };

  tk::real rho = 0.0, ru = 0.0, rv = 0.0, rw = 0.0, re = 0.0, spi = 0.0;
  bool boxmassic = false;
  if (boxmas.size() > system && !boxmas[system].empty()) {

    Assert( boxenc.size() > system && !boxenc[system].empty(),
      "Box energy content unspecified in input file" );
    auto V_ex = (boxdim[1]-boxdim[0]) * (boxdim[3]-boxdim[2]) *
      (boxdim[5]-boxdim[4]);
    rho = boxmas[system][0] / V_ex;
    spi = boxenc[system][0] * VRatio / rho;
    boxmassic = true;

  } else {

    if (boxrho.size() > system && !boxrho[system].empty()) {
      rho = boxrho[system][0];
    }
    if (boxvel.size() > system && boxvel[system].size() > 2) {
      ru = rho * boxvel[system][0];
      rv = rho * boxvel[system][1];
      rw = rho * boxvel[system][2];
    }
    if (boxpre.size() > system && !boxpre[system].empty()) {
      re = eos_totalenergy< tag::compflow >
             ( system, rho, ru/rho, rv/rho, rw/rho, boxpre[system][0] );
    }
    if (boxene.size() > system && !boxene[system].empty()) {
      const auto ux = ru/rho, uy = rv/rho, uz = rw/rho;
      const auto ke = 0.5*(ux*ux + uy*uy + uz*uz);
      re = rho * (boxene[system][0] + ke);
    }
    if (boxtem.size() > system && !boxtem[system].empty())
    {
      re = rho * boxtem[system][0] * cv.at(system).at(0);
    }

  }

  // Initiate type 'impulse' simply assigns the prescribed values to all
  // nodes within a box.
  if (inittype[system] == ctr::InitiateType::IMPULSE) {

    // superimpose on existing velocity field
    const auto u = s[1] / s[0],
               v = s[2] / s[0],
               w = s[3] / s[0];
    const auto ke = 0.5*(u*u + v*v + w*w);
    s[0] = rho;
    if (boxmassic) {
      s[1] = rho * u;
      s[2] = rho * v;
      s[3] = rho * w;
      s[4] = rho * (spi + ke);
    } else {
      s[1] = ru;
      s[2] = rv;
      s[3] = rw;
      s[4] = re;
    }
  }

  // Initiate type 'linear' assigns the prescribed values to all
  // nodes within a box. This is followed by adding a time-dependent energy
  // source term representing a planar wave-front propagating along the
  // z-direction with a velocity specified in the IC linear...end block.
  // The wave front is smoothed out to include a couple of mesh nodes.
  // see boxSrc() for further details.
  else if (inittype[system] == ctr::InitiateType::LINEAR && t < 1e-12) {

    // superimpose on existing velocity field
    const auto u = s[1]/s[0],
               v = s[2]/s[0],
               w = s[3]/s[0];
    const auto ke = 0.5*(u*u + v*v + w*w);

    // The linear-propagating source initialization can be done only based
    // on background pressure (not on temperature): The IC box can have a
    // different density than the background, while having the same
    // pressure and temperature as the background. This means, the
    // material in the box has a different specific heat (Cv) than the
    // background material. If such a box has to be initialized based
    // on temperature, the Cv of the box will have to be specified
    // separately. This is not currently supported.
    if (bgpreic.size() > system && !bgpreic[system].empty()) {
      // energy based on box density and background pressure
      spi = eos_totalenergy< tag::compflow >( system, rho, u, v, w,
        bgpreic[system][0] )/rho;
    }
    else Throw("Background pressure must be specified for box-IC with "
               "linear propagating source");

    s[0] = rho;
    s[1] = rho * u;
    s[2] = rho * v;
    s[3] = rho * w;
    s[4] = rho * (spi + ke);

  } else Throw( "IC box initiate type not implemented" );
}

} //inciter::
