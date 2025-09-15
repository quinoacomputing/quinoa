// *****************************************************************************
/*!
  \file      src/PDE/EoS/GodunovRomenski.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Godunov-Romenski equation of state for solids
  \details   This file defines functions for the Godunov-Romenski equation of
             state for solids and a hydro EoS for aluminum. These functions were
             taken from Example 1 of Barton, Philip T. "An interface-capturing
             Godunov method for the simulation of compressible solid-fluid
             problems." Journal of Computational Physics 390 (2019): 25-50.
*/
// *****************************************************************************

#include <cmath>
#include <iostream>
#include "Vector.hpp"
#include "EoS/GodunovRomenski.hpp"
#include "EoS/GetMatProp.hpp"

namespace {
  constexpr tk::real A_TRACE = 1e-4;
  constexpr tk::real RHO_FLOOR = 1e-12;     // per-solid density floor [kg/m^3]
  constexpr tk::real DET_FLOOR = 1e-18;     // for SPD guards on g

  inline tk::real alpha_eff(tk::real a) noexcept {
    return (a > A_TRACE) ? a : A_TRACE;     // only for UNWEIGHTING
  }
  inline tk::real safe_pos(tk::real x, tk::real eps) noexcept {
    return std::isfinite(x) && x > eps ? x : eps;
  }
  inline void symmetrize(std::array<std::array<tk::real,3>,3>& g) {
    for (int i=0;i<3;++i) for (int j=i+1;j<3;++j) {
      tk::real s = 0.5*(g[i][j] + g[j][i]); g[i][j]=g[j][i]=s;
    }
  }
  inline std::array<std::array<tk::real,3>,3>
  spherical_from_det(const std::array<std::array<tk::real,3>,3>& g) {
    auto detg = safe_pos(tk::determinant(g), DET_FLOOR);
    tk::real d = std::pow(detg, 1.0/3.0);
    std::array<std::array<tk::real,3>,3> G{};
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) G[i][j] = (i==j)? d : 0.0;
    return G;
  }
}

using inciter::GodunovRomenski;

GodunovRomenski::GodunovRomenski(
  tk::real gamma,
  tk::real mu,
  tk::real rho0,
  tk::real alpha,
  tk::real K0 ) :
  m_gamma(gamma),
  m_mu(mu),
  m_rho0(rho0),
  m_alpha(alpha),
  m_K0(K0)
// *************************************************************************
//  Constructor
//! \param[in] gamma Ratio of specific heats
//! \param[in] mu Constant shear modulus
//! \param[in] rho0 Unstressed density of material
//! \param[in] alpha Alpha parameter
//! \param[in] K0 K0 parameter
// *************************************************************************
{ }

tk::real
GodunovRomenski::density(
  tk::real pr,
  tk::real ) const
// *************************************************************************
//! \brief Calculate density from the material pressure and temperature
//!   using the GodunovRomenski equation of state
//! \param[in] pr Material pressure
//! \return Material density calculated using the cold compression pressure
// *************************************************************************
{
  // Quick Newton
  tk::real rho = m_rho0;
  std::size_t maxiter = 50;
  tk::real tol = 1.0e-04;
  tk::real err = tol + 1;
  for (std::size_t iter=0; iter<maxiter; ++iter)
  {
    
    auto rrho0a = std::pow(rho/m_rho0, m_alpha);
    auto p_coldcompr = m_K0/m_alpha * (rrho0a*rho/m_rho0) * (rrho0a-1.0);
    tk::real p = p_coldcompr - pr;
    auto dpdrho = DpccDrho(rho);
    auto delta = p/dpdrho;
    rho -= delta;
    err = std::sqrt(std::pow(p,2.0));
    if (err < tol) break;
  }
  return rho;
}

tk::real
GodunovRomenski::pressure(
  tk::real arho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real arhoE,
  tk::real alpha,
  std::size_t imat,
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad ) const
// *************************************************************************
//! \brief Calculate pressure from the material density, momentum, total energy
//!   and the inverse deformation gradient tensor using the GodunovRomenski
//!   equation of state
//! \param[in] arho Material partial density (alpha_k * rho_k)
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] arhoE Material total energy (alpha_k * rho_k * E_k)
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified
//!   by the calling code
//! \param[in] defgrad Material inverse deformation gradient tensor
//!   (g_k). Default is 0, so that for the single-material system,
//!   this argument can be left unspecified by the calling code
//! \return Material partial pressure (alpha_k * p_k) calculated using the
//!   GodunovRomenski EoS
//! \details The material pressure consists of the thermal component and the
//!   cold-compression component. Pressure-effects due to shear are stored in
//!   the elastic Cauchy stress tensor, not in the pressure.
// *************************************************************************
{
  if (!(alpha > 0.0)) return 0.0;

  // elastic part (α-weighted)
  std::array<std::array<tk::real,3>,3> devH{};
  const tk::real arhoEe = alpha * elasticEnergy(defgrad, devH);

  // cold-compression energy (α-weighted) with α-safe unweighting
  const tk::real aeff = alpha_eff(alpha);
  const tk::real rho  = std::max(arho/aeff, RHO_FLOOR);
  const tk::real arhoEc = alpha * coldcomprEnergy(rho);

  // thermal energy (α-weighted)
  const tk::real arhoEt = arhoE - arhoEe - arhoEc - 0.5*arho*(u*u + v*v + w*w);

  // cold-compression pressure (α-weighted) + Mie–Grüneisen thermal part
  const tk::real pcc_alpha = pressure_coldcompr(arho, alpha);
  tk::real partpressure = pcc_alpha + m_gamma * arhoEt;

  if (!std::isfinite(partpressure)) {
    // avoid throwing in hot path; return safe floor
    partpressure = 0.0;
  }
  return partpressure;
}

std::array< std::array< tk::real, 3 >, 3 >
GodunovRomenski::CauchyStress(
  tk::real,
  tk::real,
  tk::real,
  tk::real,
  tk::real,
  tk::real alpha,
  std::size_t /*imat*/,
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad ) const
// *************************************************************************
//! \brief Calculate the elastic Cauchy stress tensor from the material density,
//!   momentum, total energy, and inverse deformation gradient tensor using the
//!   GodunovRomenski equation of state
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
// //! \param[in] imat Material-id who's EoS is required. Default is 0, so that
// //!   for the single-material system, this argument can be left unspecified
// //!   by the calling code
//! \param[in] defgrad Material inverse deformation gradient tensor (g_k).
//! \return Elastic part of material Cauchy stress tensor (alpha_k * sigma_k)
//!   calculated using the GodunovRomenski EoS
// *************************************************************************
{
  std::array< std::array< tk::real, 3 >, 3 > asig{{{0,0,0}, {0,0,0}, {0,0,0}}};

  // obtain elastic contribution to energy and subtract it from pmean
  std::array< std::array< tk::real, 3 >, 3 > devH;

  // p_mean
  auto p_se = -elasticEnergy(defgrad, devH);
  auto pmean = alpha * p_se;

  // Pressure due to shear
  asig[0][0] = -pmean;
  asig[1][1] = -pmean;
  asig[2][2] = -pmean;

  // Add deviatoric component of Cauchy stress tensor
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      asig[i][j] += 2.0*m_mu*alpha*devH[i][j];
  }

  return asig;
}

tk::real
GodunovRomenski::soundspeed(
  tk::real arho,
  tk::real apr,
  tk::real alpha,
  std::size_t imat,
  const std::array< std::array< tk::real, 3 >, 3 >& ) const
// *************************************************************************
//! Calculate speed of sound from the material density and material pressure
//! \param[in] arho Material partial density (alpha_k * rho_k)
//! \param[in] apr Material partial pressure (alpha_k * p_k)
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified
//!   by the calling code
//!   (alpha * sigma_ij * n_j) projected onto the normal vector. Default is 0,
//!   so that for the single-material system, this argument can be left
//!   unspecified by the calling code
// //! \param[in] defgrad Material inverse deformation gradient tensor
// //!   (g_k) with the first dimension aligned to direction in which
// //!   wave speeds are required. Default is 0, so that for the single-material
// //!   system, this argument can be left unspecified by the calling code
//! \return Material speed of sound using the GodunovRomenski EoS
// *************************************************************************
{
  if (!(alpha > 0.0)) return 0.0;

  const tk::real aeff = alpha_eff(alpha);
  const tk::real rho  = std::max(arho/aeff, RHO_FLOOR);
  const tk::real pcc_alpha = pressure_coldcompr(arho, alpha);

  // hydro + thermal (all consistent with α-weighted apr)
  tk::real a2 = std::max<tk::real>(
                  1.0e-15,
                  DpccDrho(rho) + m_gamma * (apr - pcc_alpha) / std::max(arho, A_TRACE*rho)
                );

  // shear contribution
  a2 += (4.0/3.0) * m_mu / rho;

  tk::real a = std::sqrt(a2);
  return std::isfinite(a) ? a : 0.0;
}

tk::real
GodunovRomenski::shearspeed(
  tk::real arho,
  tk::real alpha,
  std::size_t imat ) const
// *************************************************************************
//! Calculate speed of sound from the material density and material pressure
//! \param[in] arho Material partial density (alpha_k * rho_k)
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified
//!   by the calling code
//! \return Material shear-wave speed speed using the SmallShearSolid EoS
// *************************************************************************
{
  // Approximate shear-wave speed. Ref. Barton, P. T. (2019).
  // An interface-capturing Godunov method for the simulation of compressible
  // solid-fluid problems. Journal of Computational Physics, 390, 25-50.
  if (!(alpha > 0.0)) return 0.0;
  const tk::real rho = std::max(arho/alpha_eff(alpha), RHO_FLOOR);
  tk::real a = std::sqrt(m_mu / rho);  // per-solid
  return std::isfinite(a) ? a : 0.0;
}

tk::real
GodunovRomenski::totalenergy(
  tk::real arho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real apr,
  tk::real alpha,
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad ) const
// *************************************************************************
//! \brief Calculate material specific total energy from the material
//!   density, momentum and material pressure
//! \param[in] arho Material partial density
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] apr Material partial pressure
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
//! \param[in] defgrad Material inverse deformation gradient tensor
//!   g_k. Default is 0, so that for the single-material system,
//!   this argument can be left unspecified by the calling code
//! \return Material specific total energy using the GodunovRomenski EoS
// *************************************************************************
{
  if (!(alpha > 0.0)) return 0.0;

  // alpha-weighted hydro + thermal
  const tk::real pcc_alpha = pressure_coldcompr(arho, alpha);
  const tk::real apt = apr - pcc_alpha;
  const tk::real arhoEh = apt/m_gamma + 0.5*arho*(u*u + v*v + w*w);

  // alpha-weighted elastic + cold-compression energies
  std::array<std::array<tk::real,3>,3> devH{};
  const tk::real arhoEe = alpha * elasticEnergy(defgrad, devH);
  const tk::real rho    = std::max(arho/alpha_eff(alpha), RHO_FLOOR);
  const tk::real arhoEc = alpha * coldcomprEnergy(rho);

  tk::real out = arhoEh + arhoEe + arhoEc;
  return std::isfinite(out) ? out : 0.0;
}

tk::real
GodunovRomenski::temperature(
  tk::real,
  tk::real,
  tk::real,
  tk::real,
  tk::real,
  tk::real,
  const std::array< std::array< tk::real, 3 >, 3 >& ) const
// *************************************************************************
//! \brief Calculate material temperature from the material density, and
//!   material specific total energy
// //! \param[in] arho Material partial density (alpha_k * rho_k)
// //! \param[in] u X-velocity
// //! \param[in] v Y-velocity
// //! \param[in] w Z-velocity
// //! \param[in] arhoE Material total energy (alpha_k * rho_k * E_k)
// //! \param[in] alpha Material volume fraction. Default is 1.0, so that for
// //!   the single-material system, this argument can be left unspecified by
// //!   the calling code
// //! \param[in] defgrad Material inverse deformation gradient tensor
// //!   (g_k). Default is 0, so that for the single-material system,
// //!   this argument can be left unspecified by the calling code
//! \return Material temperature using the GodunovRomenski EoS
// *************************************************************************
{
  // Temperature as a function of energy is not known.
  // So we just set a value.
  tk::real t = 300.0;

  return t;
}

tk::real
GodunovRomenski::min_eff_pressure(
  tk::real min,
  tk::real arho,
  tk::real alpha ) const
// *************************************************************************
//! Compute the minimum allowed pressure
//! \param[in] min Numerical threshold above which pressure needs to be limited
//! \param[in] arho Material partial density (alpha_k * rho_k)
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//! \return Minimum pressure allowed by physical constraints
// *************************************************************************
{
  if (!(alpha > 0.0)) return min;

  const tk::real aeff = alpha_eff(alpha);
  const tk::real rho  = std::max(arho/aeff, RHO_FLOOR);

  // p_cc per-solid
  const tk::real pcc_solid = pressure_coldcompr(aeff*rho, aeff) / aeff;

  // hydro constraint (non-negative c^2) in Mie–Grüneisen form
  tk::real pmin = min - (rho/m_gamma)*DpccDrho(rho) + pcc_solid;
  if (!std::isfinite(pmin)) pmin = min;
  return pmin;
}

tk::real
GodunovRomenski::coldcomprEnergy( tk::real rho ) const
// *************************************************************************
//! \brief Calculate cold-compression contribution to material energy from the
//!   material density
//! \param[in] rho Material density
//! \return Material cold compression energy using the GodunovRomenski EoS
// *************************************************************************
{
  auto rrho0a = std::pow(rho/m_rho0, m_alpha);
  return ( rho * m_K0/(2.0*m_rho0*m_alpha*m_alpha) * (rrho0a-1.0)*(rrho0a-1.0) );
}

tk::real
GodunovRomenski::pressure_coldcompr(
  tk::real arho,
  tk::real alpha ) const
// *************************************************************************
//! \brief Calculate cold-compression contribution to material pressure from the
//!   material density
//! \param[in] arho Material partial density (alpha_k * rho_k)
//! \param[in] alpha Material volume fraction. Default is 1.0.
//! \return Material cold compression pressure using the GodunovRomenski EoS
// *************************************************************************
{
  if (!(alpha > 0.0)) return 0.0;

  const tk::real aeff = alpha_eff(alpha);
  const tk::real rho  = std::max(arho/aeff, RHO_FLOOR);
  const tk::real rrho0a = std::pow(rho/m_rho0, m_alpha);
  return aeff * ( m_K0/m_alpha * (rrho0a * rho/m_rho0) * (rrho0a - 1.0) );
}

tk::real
GodunovRomenski::elasticEnergy(
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
  std::array< std::array< tk::real, 3 >, 3 >& devH ) const
// *************************************************************************
//! \brief Calculate elastic contribution to material energy from the material
//!   density, and deformation gradient tensor
//! \param[in] defgrad Material inverse deformation gradient tensor
//! \param[in/out] devH Deviatoric part of the Hensky tensor
//! \return Material elastic energy using the GodunovRomenski EoS
//! \details This function returns the material elastic energy, and also stores
//!   the elastic shear distortion for further use
// *************************************************************************
{
  // Symmetrize and ensure positive determinant
  auto g = defgrad;
  symmetrize(g);
  tk::real detg = tk::determinant(g);
  if (!(detg > DET_FLOOR)) {
    g = spherical_from_det(g);  // zero deviatoric if g is toxic
  }

  // Deviatoric Hencky (safe) and energy
  devH = tk::getDevHencky(g);
  tk::real rhoEe = 0.0;
  for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
    if (!std::isfinite(devH[i][j])) return 0.0;
    rhoEe += m_mu * devH[i][j] * devH[i][j];
  }
  return rhoEe;
}

tk::real
GodunovRomenski::DpccDrho( tk::real rho ) const
// *************************************************************************
//! Calculate the derivative of the cold compression pressure wrt. density
//! \param[in] rho Material partial density (alpha_k * rho_k)
//! \return Derivative of the cold compression pressure wrt. density
// *************************************************************************
{
  auto rrho0a = std::pow(rho/m_rho0, m_alpha);
  return m_K0/(m_rho0*m_alpha) *
    ((2.0*m_alpha+1.0)*(rrho0a*rrho0a) - (m_alpha+1.0)*rrho0a);
}
