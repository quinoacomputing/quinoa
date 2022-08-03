// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/MiscMultiMatFns.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Misc multi-material system functions
  \details   This file defines functions that required for multi-material
    compressible hydrodynamics.
*/
// *****************************************************************************

#include "MiscMultiMatFns.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "EoS/EoS.hpp"
#include "EoS/StiffenedGas.hpp"
#include "EoS/JWL.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

void initializeMaterialEoS( std::size_t system,
  std::vector< EoS_Base* >& mat_blk )
// *****************************************************************************
//  Initialize the material block with configured EOS
//! \param[in] system Index of system being solved
//! \param[in,out] mat_blk Material block that gets initialized
// *****************************************************************************
{
  // EoS initialization
  auto nmat =
    g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[system];
  const auto& matidxmap = g_inputdeck.get< tag::param, tag::multimat,
    tag::matidxmap >();
  for (std::size_t k=0; k<nmat; ++k) {
    if (matidxmap.get< tag::eosidx >()[k] ==
      static_cast<std::size_t>(inciter::ctr::MaterialType::STIFFENEDGAS)) {
      // query input deck to get gamma, p_c, cv
      auto g = gamma< tag::multimat >(system, k);
      auto ps = pstiff< tag::multimat >(system, k);
      auto c_v = cv< tag::multimat >(system, k);
      mat_blk.push_back(new StiffenedGas(g, ps, c_v));
    }
    else if (matidxmap.get< tag::eosidx >()[k] ==
      static_cast<std::size_t>(inciter::ctr::MaterialType::JWL)) {
      // query input deck to get jwl parameters
      auto g = gamma< tag::multimat >(system, k);
      auto c_v = cv< tag::multimat >(system, k);
      [[maybe_unused]] auto A_jwl =
        getmatprop< tag::multimat, tag::A_jwl >(system, k);
      [[maybe_unused]] auto B_jwl =
        getmatprop< tag::multimat, tag::B_jwl >(system, k);
      [[maybe_unused]] auto C_jwl =
        getmatprop< tag::multimat, tag::C_jwl >(system, k);
      [[maybe_unused]] auto R1_jwl =
        getmatprop< tag::multimat, tag::R1_jwl >(system, k);
      [[maybe_unused]] auto R2_jwl =
        getmatprop< tag::multimat, tag::R2_jwl >(system, k);
      [[maybe_unused]] auto rho0_jwl =
        getmatprop< tag::multimat, tag::rho0_jwl >(system, k);
      // TODO: update eos-ctor call when the JWL class is updated
      mat_blk.push_back(new JWL(g, 0.0, c_v));
    }
  }
}

} //inciter::
