// *****************************************************************************
/*!
  \file      src/PDE/EoS/EOS.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Polymorphic variant-style implementation for equations of state,
    where children implement specific EOS functions.
*/
// *****************************************************************************

#include "EoS/EOS.hpp"
#include "Exception.hpp"
#include "EoS/GetMatProp.hpp"

using inciter::EOS;

//! Constructor
//! \param[in] mattype Material type
//! \param[in] eq Type of PDE being solved
//! \param[in] system Index of system being solved
//! \param[in] k Material index
//! \details Based on the input enum we assign the correct material eos
EOS::EOS( ctr::MaterialType mattype,
  EqType eq,
  std::size_t system,
  std::size_t k )
{
  if (mattype == ctr::MaterialType::STIFFENEDGAS) {
    // query input deck to get gamma, p_c, cv
    tk::real g, ps, c_v;
    if (eq == EqType::compflow) {
      g = getmatprop< tag::compflow, tag::gamma >(system, k);
      ps = getmatprop< tag::compflow, tag::pstiff >(system, k);
      c_v = getmatprop< tag::compflow, tag::cv >(system, k);
    }
    else {
      g = getmatprop< tag::multimat, tag::gamma >(system, k);
      ps = getmatprop< tag::multimat, tag::pstiff >(system, k);
      c_v = getmatprop< tag::multimat, tag::cv >(system, k);
    }
    m_material = StiffenedGas(g, ps, c_v);
  }
  else if (mattype == ctr::MaterialType::JWL) {
    if (eq == EqType::compflow) Throw("JWL not set up for PDE type");
    // query input deck to get jwl parameters
    auto w = getmatprop< tag::multimat, tag::w_gru >(system, k);
    auto c_v = getmatprop< tag::multimat, tag::cv >(system, k);
    auto A_jwl = getmatprop< tag::multimat, tag::A_jwl >(system, k);
    auto B_jwl = getmatprop< tag::multimat, tag::B_jwl >(system, k);
    //[[maybe_unused]] auto C_jwl =
    //  getmatprop< tag::multimat, tag::C_jwl >(system, k);
    auto R1_jwl = getmatprop< tag::multimat, tag::R1_jwl >(system, k);
    auto R2_jwl = getmatprop< tag::multimat, tag::R2_jwl >(system, k);
    auto rho0_jwl = getmatprop< tag::multimat, tag::rho0_jwl >(system, k);
    auto de_jwl = getmatprop< tag::multimat, tag::de_jwl >(system, k);
    auto rhor_jwl = getmatprop< tag::multimat, tag::rhor_jwl >(system, k);
    auto Tr_jwl = getmatprop< tag::multimat, tag::Tr_jwl >(system, k);
    auto Pr_jwl = getmatprop< tag::multimat, tag::Pr_jwl >(system, k);
    m_material = JWL(w, c_v, rho0_jwl, de_jwl, rhor_jwl, Tr_jwl, Pr_jwl, A_jwl,
      B_jwl, R1_jwl, R2_jwl);
  }
  else Throw( "Unknown EOS for material " + std::to_string(k+1) );
}
