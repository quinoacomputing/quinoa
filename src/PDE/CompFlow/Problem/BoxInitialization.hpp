// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/BoxInitialization.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     User-defined box initialization
  \details   This file defines functions for initializing solutions for
    compressible single-material equations inside the user-defined box.
*/
// *****************************************************************************
#ifndef BoxInitialization_h
#define BoxInitialization_h

#include "Fields.hpp"
#include "EoS/EOS.hpp"
#include "ContainerUtil.hpp"
#include "Control/Inciter/Types.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

//! Set the solution in the user-defined IC box

template< class B >
void initializeBox( std::size_t,
                    const std::vector< EOS >& mat_blk,
                    tk::real VRatio,
                    tk::real V_ex,
                    tk::real t,
                    const B& b,
                    tk::real bgpreic,
                    tk::real cv,
                    std::vector< tk::real >& s )
// *****************************************************************************
// Set the solution in the user-defined IC box/block
//! \tparam B IC-block type to operate, ctr::box, or ctr::meshblock
//! \param[in] system Equation system index
//! \param[in] VRatio Ratio of exact box volume to discrete box volume
//! \param[in] V_ex Exact box volume
//! \param[in] t Physical time
//! \param[in] b IC box configuration to use
//! \param[in] bgpreic Background pressure user input
//! \param[in] cv Specific heats ratio user input
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
  const auto& initiate = b.template get< tag::initiate >();
  auto inittype = initiate.template get< tag::init >();

  auto boxrho = b.template get< tag::density >();
  const auto& boxvel = b.template get< tag::velocity >();
  auto boxpre = b.template get< tag::pressure >();
  auto boxene = b.template get< tag::energy >();
  auto boxtem = b.template get< tag::temperature >();
  auto boxmas = b.template get< tag::mass >();
  auto boxenc = b.template get< tag::energy_content >();

  tk::real rho = 0.0, ru = 0.0, rv = 0.0, rw = 0.0, re = 0.0, spi = 0.0;
  bool boxmassic = false;
  if (boxmas > 0.0) {
    Assert( boxenc > 0.0, "Box energy content must be nonzero" );
    rho = boxmas / V_ex;
    spi = boxenc * VRatio / rho;
    boxmassic = true;
  } else {
    if (boxrho > 0.0) rho = boxrho;
    if (boxvel.size() == 3) {
      ru = rho * boxvel[0];
      rv = rho * boxvel[1];
      rw = rho * boxvel[2];
    }
    if (boxpre > 0.0) {
      re = mat_blk[0].compute< EOS::totalenergy >( rho, ru/rho, rv/rho, rw/rho,
        boxpre );
    }
    if (boxene > 0.0) {
      auto ux = ru/rho, uy = rv/rho, uz = rw/rho;
      auto ke = 0.5*(ux*ux + uy*uy + uz*uz);
      re = rho * (boxene + ke);
    }
    if (boxtem > 0.0) re = rho * boxtem * cv;
  }

  // Initiate type 'impulse' simply assigns the prescribed values to all
  // nodes within a box.
  if (inittype == ctr::InitiateType::IMPULSE) {
    // superimpose on existing velocity field
    auto u = s[1]/s[0], v = s[2]/s[0], w = s[3]/s[0];
    auto ke = 0.5*(u*u + v*v + w*w);
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
  // Initiate type 'linear' assigns the prescribed values to all nodes
  // within a box. This is followed by adding a time-dependent energy
  // source term representing a planar wave-front propagating along the
  // z-direction with a velocity specified in the IC linear...end block.
  // The wave front is smoothed out to include a couple of mesh nodes.
  // see boxSrc() for further details.
  else if (inittype == ctr::InitiateType::LINEAR && t < 1e-12) {

    // superimpose on existing velocity field
    const auto u = s[1]/s[0],
               v = s[2]/s[0],
               w = s[3]/s[0];
    const auto ke = 0.5*(u*u + v*v + w*w);

    // The linear-propagating source initialization can be done only
    // based on background pressure (not on temperature): The IC box can
    // have a different density than the background, while having the
    // same pressure and temperature as the background. This means, the
    // material in the box has a different specific heat (Cv) than the
    // background material. If such a box has to be initialized based on
    // temperature, the Cv of the box will have to be specified
    // separately. This is not currently supported.
    if (bgpreic > 0.0) {
      // energy based on box density and background pressure
      spi = mat_blk[0].compute< EOS::totalenergy >(rho, u, v, w, bgpreic) / rho;
    } else Throw( "Background pressure must be specified for box-IC "
                  "with linear propagating source");

    s[0] = rho;
    s[1] = rho * u;
    s[2] = rho * v;
    s[3] = rho * w;
    s[4] = rho * (spi + ke);

  } else Throw( "IC box initiate type not implemented" );
}

} //inciter::

#endif // BoxInitialization_h
