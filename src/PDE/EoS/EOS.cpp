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
//! \param[in] k Material index
//! \details Based on the input enum we assign the correct material eos
EOS::EOS( ctr::MaterialType mattype, EqType eq, std::size_t k )
{
  if (mattype == ctr::MaterialType::STIFFENEDGAS) {
    // query input deck to get gamma, p_c, cv
    auto g = getmatprop< newtag::gamma >(k);
    auto ps = getmatprop< newtag::pstiff >(k);
    auto c_v = getmatprop< newtag::cv >(k);
    m_material = StiffenedGas(g, ps, c_v);
  }
  else if (mattype == ctr::MaterialType::JWL) {
    if (eq == EqType::compflow) Throw("JWL not set up for PDE type");
    // query input deck to get jwl parameters
    auto w = getmatprop< newtag::w_gru >(k);
    auto c_v = getmatprop< newtag::cv >(k);
    auto A_jwl = getmatprop< newtag::A_jwl >(k);
    auto B_jwl = getmatprop< newtag::B_jwl >(k);
    //[[maybe_unused]] auto C_jwl =
    //  getmatprop< newtag::multimat, newtag::C_jwl >(k);
    auto R1_jwl = getmatprop< newtag::R1_jwl >(k);
    auto R2_jwl = getmatprop< newtag::R2_jwl >(k);
    auto rho0_jwl = getmatprop< newtag::rho0_jwl >(k);
    auto de_jwl = getmatprop< newtag::de_jwl >(k);
    auto rhor_jwl = getmatprop< newtag::rhor_jwl >(k);
    auto Tr_jwl = getmatprop< newtag::Tr_jwl >(k);
    auto Pr_jwl = getmatprop< newtag::Pr_jwl >(k);
    m_material = JWL(w, c_v, rho0_jwl, de_jwl, rhor_jwl, Tr_jwl, Pr_jwl, A_jwl,
      B_jwl, R1_jwl, R2_jwl);
  }
  else if (mattype == ctr::MaterialType::SMALLSHEARSOLID) {
    if (eq == EqType::compflow)
      Throw("SmallShearSolid not set up for PDE type");
    // query input deck for SmallShearSolid parameters
    auto g = getmatprop< newtag::gamma >(k);
    auto ps = getmatprop< newtag::pstiff >(k);
    auto c_v = getmatprop< newtag::cv >(k);
    auto mu = getmatprop< newtag::mu >(k);
    m_material = SmallShearSolid(g, ps, c_v, mu);
  }
  else Throw( "Unknown EOS for material " + std::to_string(k+1) );
}
