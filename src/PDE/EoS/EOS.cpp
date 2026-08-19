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

#include <new>
#include <utility>

using inciter::EOS;

//! Constructor
//! \param[in] mattype Material type
//! \param[in] eq Type of PDE being solved
//! \param[in] k Material/species index
//! \details Based on the input enum we assign the correct material eos
EOS::EOS( ctr::MaterialType mattype, EqType eq, std::size_t k )
{
  // Multi-material hydrodynamics PDEs
  if (eq == EqType::multimat) {
  if (mattype == ctr::MaterialType::STIFFENEDGAS) {
    // query input deck to get gamma, p_c, cv
    auto g = getmatprop< tag::gamma >(k);
    auto ps = getmatprop< tag::pstiff >(k);
    auto c_v = getmatprop< tag::cv >(k);
    auto mu = getmatprop< tag::mu >(k);
    type = EOSType::StiffenedGas;
    new (&m_material.stiffenedGas)
      StiffenedGas(g, ps, c_v, mu);
    m_active = true;
  }
  else if (mattype == ctr::MaterialType::JWL) {
    // query input deck to get jwl parameters
    auto w = getmatprop< tag::w_gru >(k);
    auto c_v = getmatprop< tag::cv >(k);
    auto A_jwl = getmatprop< tag::A_jwl >(k);
    auto B_jwl = getmatprop< tag::B_jwl >(k);
    //[[maybe_unused]] auto C_jwl =
    //  getmatprop< tag::multimat, tag::C_jwl >(k);
    auto R1_jwl = getmatprop< tag::R1_jwl >(k);
    auto R2_jwl = getmatprop< tag::R2_jwl >(k);
    auto rho0_jwl = getmatprop< tag::rho0_jwl >(k);
    auto de_jwl = getmatprop< tag::de_jwl >(k);
    auto rhor_jwl = getmatprop< tag::rhor_jwl >(k);
    auto Tr_jwl = getmatprop< tag::Tr_jwl >(k);
    auto Pr_jwl = getmatprop< tag::Pr_jwl >(k);
    type = EOSType::JWL;
    new (&m_material.jwl) JWL(w, c_v, rho0_jwl, de_jwl, rhor_jwl, Tr_jwl, Pr_jwl, A_jwl,
      B_jwl, R1_jwl, R2_jwl);
    m_active = true;
  }
  else if (mattype == ctr::MaterialType::SMALLSHEARSOLID) {
    // query input deck for SmallShearSolid parameters
    auto g = getmatprop< tag::gamma >(k);
    auto ps = getmatprop< tag::pstiff >(k);
    auto c_v = getmatprop< tag::cv >(k);
    auto mu = getmatprop< tag::mu >(k);
    type = EOSType::SmallShearSolid;
    new (&m_material.smallShearSolid) SmallShearSolid(g, ps, c_v, mu);
    m_active = true;
  }
  else if (mattype == ctr::MaterialType::LINEARMIEGRUNEISEN) {
    // query input deck for LinearMieGruneisen parameters
    auto gamma0 = getmatprop< tag::w_gru >(k);
    auto alpha = getmatprop< tag::alpha >(k);
    auto c0 = getmatprop< tag::c0 >(k);
    auto s1 = getmatprop< tag::s1 >(k);
    auto c_v = getmatprop< tag::cv >(k);
    auto mu = getmatprop< tag::mu >(k);
    auto rho0_gr = getmatprop< tag::rho0_jwl >(k);
    type = EOSType::LinearMieGruneisen;
    new (&m_material.linearMieGruneisen)
      LinearMieGruneisen(gamma0, rho0_gr, alpha, c0, s1, c_v, mu);
    m_active = true;
  }
  else if (mattype == ctr::MaterialType::WILKINSALUMINUM) {
    // query input deck for Wilkins parameters
    auto g = getmatprop< tag::gamma >(k);
    auto c_v = getmatprop< tag::cv >(k);
    auto mu = getmatprop< tag::mu >(k);
    type = EOSType::WilkinsAluminum;
    new (&m_material.wilkinsAluminum) WilkinsAluminum(g, c_v, mu);
    m_active = true;
  }
  else if (mattype == ctr::MaterialType::GODUNOVROMENSKI) {
    // query input deck for Wilkins parameters
    auto g = getmatprop< tag::gamma >(k);
    auto mu = getmatprop< tag::mu >(k);
    auto rho0_gr = getmatprop< tag::rho0_jwl >(k);
    auto alpha = getmatprop< tag::alpha >(k);
    auto K0 = getmatprop< tag::K0 >(k);
    type = EOSType::GodunovRomenski;
    new (&m_material.godunovRomenski) GodunovRomenski(g, mu, rho0_gr, alpha, K0);
    m_active = true;
  }
  else Throw( "Unknown EOS for material " + std::to_string(k+1) );
  }

  // Multi-species PDEs
  else if (eq == EqType::multispecies) {
    if (mattype == ctr::MaterialType::STIFFENEDGAS) {
      // query input deck to get gamma, p_c, cv
      auto g = getspecprop< tag::gamma >(k);
      auto ps = getspecprop< tag::pstiff >(k);
      auto c_v = getspecprop< tag::cv >(k);
      auto mu = getspecprop< tag::mu >(k);
      type = EOSType::StiffenedGas;
      new (&m_material.stiffenedGas)
        StiffenedGas(g, ps, c_v, mu);
      m_active = true;
    }
    else if (mattype == ctr::MaterialType::THERMALLYPERFECTGAS) {
      // query input deck for ThermallyPerfectGas parameters
      auto R = getspecprop< tag::R >(k);
      // assume only one type of species
      auto cp_coeff =
        g_inputdeck.get< tag::species >()[0].get< tag::cp_coeff >()[k];
      auto t_range =
        g_inputdeck.get< tag::species >()[0].get< tag::t_range >()[k];
      auto dH_ref = getspecprop< tag::dH_ref >(k);
      auto mu = getspecprop< tag::mu >(k);
      auto temp_ref = getspecprop< tag::temp_ref >(k);
      auto mu_ref = getspecprop< tag::mu_ref >(k);
      auto C = getspecprop< tag::C >(k);
      auto Sutherland = g_inputdeck.get< tag::multispecies >().get< tag::Sutherland >();
      type = EOSType::ThermallyPerfectGas;
      new (&m_material.thermallyPerfectGas)
        ThermallyPerfectGas(R, cp_coeff, t_range, dH_ref, mu, temp_ref, mu_ref, C, Sutherland);
      m_active = true;
    }
    else Throw( "Unknown EOS for species " + std::to_string(k+1) );
  }

  // Compressible Euler (single material) PDEs
  else if (eq == EqType::compflow) {
    if (mattype == ctr::MaterialType::STIFFENEDGAS) {
      // query input deck to get gamma, p_c, cv
      auto g = getmatprop< tag::gamma >(k);
      auto ps = getmatprop< tag::pstiff >(k);
      auto c_v = getmatprop< tag::cv >(k);
      type = EOSType::StiffenedGas;
      new (&m_material.stiffenedGas) StiffenedGas(g, ps, c_v);
      m_active = true;
    }
    else Throw( "Unknown EOS for material " + std::to_string(k+1) );
  }
  else
    Throw( "Unknown PDE type encountered in EOS ctor" );
}

/*
// Destroy
EOS::~EOS()
{
  destroy();
}

void EOS::destroy() noexcept
{
  if (!m_active) return;

  switch (type) {
    case EOSType::StiffenedGas:
      m_material.stiffenedGas.~StiffenedGas();
      break;

    case EOSType::JWL:
      m_material.jwl.~JWL();
      break;

    case EOSType::SmallShearSolid:
      m_material.smallShearSolid.~SmallShearSolid();
      break;

    case EOSType::LinearMieGruneisen:
      m_material.linearMieGruneisen.~LinearMieGruneisen();
      break;

    case EOSType::WilkinsAluminum:
      m_material.wilkinsAluminum.~WilkinsAluminum();
      break;

    case EOSType::GodunovRomenski:
      m_material.godunovRomenski.~GodunovRomenski();
      break;

    case EOSType::ThermallyPerfectGas:
      m_material.thermallyPerfectGas.~ThermallyPerfectGas();
      break;
  }

  m_active = false;
}

// Move constructor
void EOS::moveFrom(EOS&& other)
{
  if (!other.m_active) {
    m_active = false;
    return;
  }

  type = other.type;

  switch (type) {
    case EOSType::StiffenedGas:
      new (&m_material.stiffenedGas)
        StiffenedGas(std::move(other.m_material.stiffenedGas));
      break;

    case EOSType::JWL:
      new (&m_material.jwl)
        JWL(std::move(other.m_material.jwl));
      break;

    case EOSType::SmallShearSolid:
      new (&m_material.smallShearSolid)
        SmallShearSolid(std::move(other.m_material.smallShearSolid));
      break;

    case EOSType::LinearMieGruneisen:
      new (&m_material.linearMieGruneisen)
        LinearMieGruneisen(
          std::move(other.m_material.linearMieGruneisen)
        );
      break;

    case EOSType::WilkinsAluminum:
      new (&m_material.wilkinsAluminum)
        WilkinsAluminum(std::move(other.m_material.wilkinsAluminum));
      break;

    case EOSType::GodunovRomenski:
      new (&m_material.godunovRomenski)
        GodunovRomenski(std::move(other.m_material.godunovRomenski));
      break;

    case EOSType::ThermallyPerfectGas:
      new (&m_material.thermallyPerfectGas)
        ThermallyPerfectGas(
          std::move(other.m_material.thermallyPerfectGas)
        );
      break;
  }

  m_active = true;
}

// Copy constructor
void EOS::copyFrom(const EOS& other)
{
  if (!other.m_active) {
    m_active = false;
    return;
  }

  type = other.type;

  switch (type) {
    case EOSType::StiffenedGas:
      new (&m_material.stiffenedGas)
        StiffenedGas(other.m_material.stiffenedGas);
      break;

    case EOSType::JWL:
      new (&m_material.jwl)
        JWL(other.m_material.jwl);
      break;

    case EOSType::SmallShearSolid:
      new (&m_material.smallShearSolid)
        SmallShearSolid(other.m_material.smallShearSolid);
      break;

    case EOSType::LinearMieGruneisen:
      new (&m_material.linearMieGruneisen)
        LinearMieGruneisen(
          other.m_material.linearMieGruneisen
        );
      break;

    case EOSType::WilkinsAluminum:
      new (&m_material.wilkinsAluminum)
        WilkinsAluminum(other.m_material.wilkinsAluminum);
      break;

    case EOSType::GodunovRomenski:
      new (&m_material.godunovRomenski)
        GodunovRomenski(other.m_material.godunovRomenski);
      break;

    case EOSType::ThermallyPerfectGas:
      new (&m_material.thermallyPerfectGas)
        ThermallyPerfectGas(
          other.m_material.thermallyPerfectGas
        );
      break;
  }

  m_active = true;
}

EOS::EOS(const EOS& other)
{
  copyFrom(other);
}

EOS& EOS::operator=(const EOS& other)
{
  if (this != &other) {
    destroy();
    copyFrom(other);
  }

  return *this;
}

EOS::EOS(EOS&& other)
{
  moveFrom(std::move(other));
}

EOS& EOS::operator=(EOS&& other)
{
  if (this != &other) {
    destroy();
    moveFrom(std::move(other));
  }

  return *this;
}
*/
