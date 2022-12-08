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
#include "ContainerUtil.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

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

void MultiMatPhysicsEnergyPill::
physSrc( std::size_t system,
  std::size_t nmat,
  tk::real t,
  const tk::Fields& geoElem,
  const std::unordered_map< std::size_t, std::set< std::size_t > >& elemblkid,
  tk::Fields& R,
  int& engSrcAdded ) const
// *****************************************************************************
//! Compute sources corresponding to a propagating front in user-defined box
//! \param[in] system Index for equation systems
//! \param[in] nmat Number of materials
//! \param[in] t Physical time
//! \param[in] geoElem Element geometry array
//! \param[in] elemblkid Element ids associated with mesh block ids where
//!   user ICs are set
//! \param[in,out] R Right-hand side vector
//! \param[in,out] engSrcAdded Whether the energy source was added
//! \details This function adds the energy source corresponding to a
//!   spherically growing wave-front propagating with a user-specified
//!   velocity, within a user-configured mesh block initial condition.
//!   Example (SI) units of the quantities involved:
//!    * internal energy content (energy per unit volume): J/m^3
//!    * specific energy (internal energy per unit mass): J/kg
// *****************************************************************************
{
  const auto& icmbk = g_inputdeck.get< tag::param, tag::multimat, tag::ic,
    tag::meshblock >();
  if (icmbk.size() > system) {
    for (const auto& mb : icmbk[system]) { // for all blocks
      auto blid = mb.get< tag::blockid >();
      if (elemblkid.find(blid) != elemblkid.end()) { // if elements exist in blk
        const auto& initiate = mb.template get< tag::initiate >();
        auto inittype = initiate.template get< tag::init >();
        if (inittype == ctr::InitiateType::LINEAR) { // if propagating src

          const auto& blkelems = tk::cref_find(elemblkid,blid);
          std::size_t ncomp = R.nprop();

          auto enc = mb.template get< tag::energy_content >();
          Assert( enc > 0.0, "Box energy content must be nonzero" );
          const auto& x0_front = mb.template get< tag::initiate, tag::point >();
          Assert(x0_front.size()==3, "Incorrectly sized front initial location");
          auto blkmatid = mb.template get< tag::materialid >();

          // determine times at which sourcing is initialized and terminated
          auto v_front = mb.template get< tag::initiate, tag::velocity >();
          auto w_front = mb.template get< tag::initiate, tag::front_width >();
          auto tInit = mb.template get< tag::initiate, tag::init_time >();

          if (t >= tInit) {
            // current radius of front
            tk::real r_front = v_front * (t-tInit);
            // arbitrary shape form
            auto amplE = enc * v_front / w_front;

            for (auto e : blkelems) {
              std::array< tk::real, 3 > node{{ geoElem(e,1), geoElem(e,2),
                geoElem(e,3) }};

              auto r_e = std::sqrt(
                (node[0]-x0_front[0])*(node[0]-x0_front[0]) +
                (node[1]-x0_front[1])*(node[1]-x0_front[1]) +
                (node[2]-x0_front[2])*(node[2]-x0_front[2]) );

              // if element centroid lies within spherical shell add sources
              if (r_e >= r_front && r_e <= r_front+w_front) {
                engSrcAdded = 1;
                // Compute the source term variable
                std::vector< tk::real > s(ncomp, 0.0);
                // arbitrary shape form
                s[energyIdx(nmat,blkmatid-1)] = amplE;

                // Add the source term to the rhs
                for (std::size_t c=0; c<ncomp; ++c)
                {
                  R(e, c) += geoElem(e,0) * s[c];
                }
              }
            }
          }

        }
      }
    }
  }
}
