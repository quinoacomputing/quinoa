// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Physics/FVEnergyPill.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Physics policy for the Euler equation governing multi-material flow
    using a finite volume method
  \details   This file defines a Physics policy class for the compressible
    flow equations class fv::MultiMat, defined in PDE/MultiMat/FVMultiMat.h.
    This specific class allows energy pill initialization of a user defined
    box/block for multi-material flow. See PDE/MultiMat/Physics/FV.h for general
    requirements on Physics policy classes for fv::MultiMat.
*/
// *****************************************************************************

#include "FVEnergyPill.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::fv::MultiMatPhysicsEnergyPill;

tk::real
MultiMatPhysicsEnergyPill::dtRestriction( std::size_t system,
  const tk::Fields& geoElem,
  std::size_t nelem,
  const int engSrcAd ) const
// *****************************************************************************
//  Compute the time step size restriction based on this physics
//! \param[in] system Index for equation systems
//! \param[in] geoElem Element geometry array
//! \param[in] nelem Number of elements
//! \param[in] engSrcAd Whether the energy source was added
//! \return Maximum allowable time step based on front propagation speed
// *****************************************************************************
{
  auto mindt = std::numeric_limits< tk::real >::max();
  // determine front propagation speed if relevant energy sources were added
  tk::real v_front(0.0);
  if (engSrcAd == 1) {
    const auto& icmbk = g_inputdeck.get< tag::param, tag::multimat, tag::ic,
      tag::meshblock >();
    if (icmbk.size() > system) {
      for (const auto& b : icmbk[system]) { // for all blocks
        auto inittype = b.template get< tag::initiate, tag::init >();
        if (inittype == ctr::InitiateType::LINEAR) {
          v_front = std::max(v_front,
            b.template get< tag::initiate, tag::velocity >());
        }
      }
    }
  }

  for (std::size_t e=0; e<nelem; ++e)
  {
    // characteristic length (radius of insphere)
    auto dx = std::min(std::cbrt(geoElem(e,0)), geoElem(e,4))
      /std::sqrt(24.0);

    // element dt
    if (std::abs(v_front) > 1e-8) mindt = dx/v_front;
  }

  return mindt;
}
