// *****************************************************************************
/*!
  \file      src/PDE/EoS/LinearMieGruneisen.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Linear Mie-Gruneisen equation of state (a.k.a. shock-wave EOS)
  \details   This file defines functions for the LinearMieGruneisen equation of
             state for solids. The elastic contribution follows NeoHookeanSolid.
             The hydrodynamic contribution uses a linear Us-Up Hugoniot and a
             density-dependent Gruneisen coefficient. Ref. Shyue, K. M. (2001).
             A fluid-mixture type algorithm for compressible multicomponent flow
             with Mie–Grüneisen equation of state. J Comp Phys, 171(2), 678-707.
*/
// *****************************************************************************

#include <algorithm>
#include <cmath>
#include <iostream>

#include "Vector.hpp"
#include "EoS/LinearMieGruneisen.hpp"

using inciter::LinearMieGruneisen;

LinearMieGruneisen::LinearMieGruneisen(
  tk::real gamma0,
  tk::real rho0,
  tk::real alpha,
  tk::real c0,
  tk::real s1,
  tk::real cv,
  tk::real mu ) :
  m_gamma0(gamma0),
  m_rho0(rho0),
  m_alpha(alpha),
  m_c0(c0),
  m_s1(s1),
  m_cv(cv),
  m_mu(mu)
// *************************************************************************
//  Constructor
//! \param[in] gamma0 Reference Gruneisen coefficient
//! \param[in] rho0 Initial density
//! \param[in] alpha Density-dependence exponent for Gruneisen coefficient
//! \param[in] c0 Bulk sound speed coefficient in linear Us-Up relation
//! \param[in] s1 Slope coefficient in linear Us-Up relation
//! \param[in] cv Specific heat at constant volume
//! \param[in] mu Constant shear modulus
// *************************************************************************
{ }

tk::real
LinearMieGruneisen::density(
  tk::real pr,
  tk::real temp ) const
// *************************************************************************
//! \brief Calculate density from the material pressure and temperature
//!   using the LinearMieGruneisen equation of state
//! \param[in] pr Material pressure
//! \param[in] temp Material temperature
//! \return Material density calculated using the LinearMieGruneisen EoS
// *************************************************************************
{
  auto rho = m_rho0;
  const std::size_t maxiter = 50;
  const tk::real tol = 1.0e-10;

  for (std::size_t iter=0; iter<maxiter; ++iter) {
    const auto p = hugoniotPressure(rho) + gruneisen(rho)*rho*m_cv*temp;
    const auto dpdrho = dHugoniotPressureDrho(rho) +
      (gruneisen(rho) + rho*dGruneisenDrho(rho))*m_cv*temp;
    const auto rhoold = rho;
    const auto delta = (p - pr)/dpdrho;
    rho -= delta;
    if (rho <= 0.0) rho = 0.5*rhoold;
    if (std::abs(delta) <= tol*std::max(1.0, std::abs(rho))) break;
  }

  return rho;
}

tk::real
LinearMieGruneisen::pressure(
  tk::real arho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real arhoE,
  tk::real alpha,
  std::size_t /*imat*/,
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
  tk::real damage ) const
// *************************************************************************
//! \brief Calculate pressure from the material density, momentum, total energy
//!   and the inverse deformation gradient tensor using the LinearMieGruneisen
//!   equation of state
//! \param[in] arho Material partial density (alpha_k * rho_k)
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] arhoE Material total energy (alpha_k * rho_k * E_k)
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
// //! \param[in] imat Material-id who's EoS is required. Default is 0, so that
// //!   for the single-material system, this argument can be left unspecified
// //!   by the calling code
//! \param[in] defgrad Material inverse deformation gradient tensor
//!   (g_k). Default is 0, so that for the single-material system,
//!   this argument can be left unspecified by the calling code
//! \return Material partial pressure (alpha_k * p_k) calculated using the
//!   LinearMieGruneisen EoS
// *************************************************************************
{
  // obtain elastic contribution to energy
  tk::real eps2;
  const auto arhoEe = alpha*elasticEnergy(defgrad, eps2, damage);
  // obtain hydrodynamic internal energy
  const auto arhoEi = arhoE - arhoEe - 0.5*arho*(u*u + v*v + w*w);

  const auto rho = arho/alpha;
  auto partpressure = alpha*hugoniotPressure(rho) +
    gruneisen(rho)*(arhoEi - arho*hugoniotEnergy(rho));

  partpressure = std::max(min_eff_pressure(1e-10, arho, alpha), partpressure);

  return partpressure;
}

tk::real
LinearMieGruneisen::pressure_coldcompr(
  tk::real arho,
  tk::real alpha ) const
// *************************************************************************
//! \brief Calculate cold-compression component of pressure
//! \param[in] arho Material partial density
//! \param[in] alpha Material volume fraction
//! \return Material partial cold-compression pressure
// *************************************************************************
{
  return alpha*hugoniotPressure(arho/alpha);
}

std::array< std::array< tk::real, 3 >, 3 >
LinearMieGruneisen::CauchyStress(
  tk::real alpha,
  std::size_t /*imat*/,
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
  tk::real damage ) const
// *************************************************************************
//! \brief Calculate the elastic Cauchy stress tensor from the material
//!   inverse deformation gradient tensor using the LinearMieGruneisen EOS
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
// //! \param[in] imat Material-id who's EoS is required. Default is 0, so that
// //!   for the single-material system, this argument can be left unspecified
// //!   by the calling code
//! \param[in] defgrad Material inverse deformation gradient tensor (g_k).
//! \return Material Cauchy stress tensor (alpha_k * sigma_k) calculated using
//!   the LinearMieGruneisen EoS
// *************************************************************************
{
  std::array< std::array< tk::real, 3 >, 3 > asig{{{0,0,0}, {0,0,0}, {0,0,0}}};

  // obtain elastic contribution to energy
  tk::real eps2;
  elasticEnergy(defgrad, eps2, damage);

  // p_mean
  auto pmean = - alpha * m_mu * eps2;

  // Volumetric component of Cauchy stress tensor
  asig[0][0] = -pmean;
  asig[1][1] = -pmean;
  asig[2][2] = -pmean;

  // Deviatoric (trace-free) part of volume-preserving left Cauchy-Green tensor
  auto devbt = tk::getLeftCauchyGreen(defgrad);
  auto detb = std::pow(tk::determinant(devbt), 1.0/3.0);
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      devbt[i][j] /= detb;
  }
  auto trbt = (devbt[0][0]+devbt[1][1]+devbt[2][2])/3.0;
  devbt[0][0] -= trbt;
  devbt[1][1] -= trbt;
  devbt[2][2] -= trbt;

  // Add deviatoric component of Cauchy stress tensor
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      asig[i][j] += m_mu*alpha*devbt[i][j];
  }

  return asig;
}

tk::real
LinearMieGruneisen::soundspeed(
  tk::real arho,
  tk::real apr,
  tk::real alpha,
  std::size_t imat,
  const std::array< std::array< tk::real, 3 >, 3 >& /*defgrad*/,
  tk::real damage ) const
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
//  //! \param[in] defgrad Material inverse deformation gradient tensor
//  //!   (g_k) with the first dimension aligned to direction in which
//  //!   wave speeds are required. Default is 0, so that for the single-material
//  //!   system, this argument can be left unspecified by the calling code
//! \return Material speed of sound using the LinearMieGruneisen EoS
// *************************************************************************
{
  // Approximated elastic contribution, from Barton, P. T. (2019).
  // An interface-capturing Godunov method for the simulation of compressible
  // solid-fluid problems. Journal of Computational Physics, 390, 25-50
  auto al_eff = std::max( 1.0e-14, alpha );
  tk::real a = (4.0/3.0) * m_mu * al_eff / arho;

  // Hydrodynamic contribution from the Mie-Gruneisen derivatives.
  const auto rho = arho/al_eff;
  const auto gamma = gruneisen(rho);
  const auto e_thermal =
    (apr - al_eff*hugoniotPressure(rho))/(gamma*arho);
  a += dHugoniotPressureDrho(rho)
     + (gamma + rho*dGruneisenDrho(rho))*e_thermal
     - gamma*rho*dHugoniotEnergyDrho(rho)
     + gamma*apr/arho;

  // Compute square root
  a = std::sqrt(std::max(1.0e-15, a));

  // check sound speed divergence
  if (!std::isfinite(a)) {
    std::cout << "Material-id:      " << imat << std::endl;
    std::cout << "Volume-fraction:  " << alpha << std::endl;
    std::cout << "Partial density:  " << arho << std::endl;
    std::cout << "Partial pressure: " << apr << std::endl;
    Throw("Material-" + std::to_string(imat) + " has nan/inf sound speed: "
      + std::to_string(a) + ", material volume fraction: " +
      std::to_string(alpha));
  }

  return a;
}

tk::real
LinearMieGruneisen::shearspeed(
  tk::real arho,
  tk::real alpha,
  std::size_t imat,
  tk::real damage ) const
// *************************************************************************
//! Calculate speed of sound from the material density and material pressure
//! \param[in] arho Material partial density (alpha_k * rho_k)
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified
//!   by the calling code
//! \return Material shear-wave speed speed using the LinearMieGruneisen EoS
// *************************************************************************
{
  // Approximate shear-wave speed. Ref. Barton, P. T. (2019).
  // An interface-capturing Godunov method for the simulation of compressible
  // solid-fluid problems. Journal of Computational Physics, 390, 25-50.
  auto al_eff = std::max( 1e-14, alpha );
  tk::real a = std::sqrt(al_eff*m_mu/arho);

  // check shear-wave speed divergence
  if (!std::isfinite(a)) {
    std::cout << "Material-id:      " << imat << std::endl;
    std::cout << "Volume-fraction:  " << alpha << std::endl;
    std::cout << "Partial density:  " << arho << std::endl;
    Throw("Material-" + std::to_string(imat) + " has nan/inf shear-wave speed: "
      + std::to_string(a) + ", material volume fraction: " +
      std::to_string(alpha));
  }

  return a;
}

tk::real
LinearMieGruneisen::totalenergy(
  tk::real arho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real apr,
  tk::real alpha,
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
  tk::real damage ) const
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
//! \return Material specific total energy using the LinearMieGruneisen EoS
// *************************************************************************
{
  const auto rho = arho/alpha;

  // obtain hydrodynamic and kinetic energy
  tk::real arhoEh = arho*hugoniotEnergy(rho) +
    (apr - alpha*hugoniotPressure(rho))/gruneisen(rho) +
    0.5 * arho * (u*u + v*v + w*w);
  // obtain elastic contribution to energy
  tk::real eps2;
  tk::real arhoEe = alpha*elasticEnergy(defgrad, eps2, damage);

  return (arhoEh + arhoEe);
}

tk::real
LinearMieGruneisen::temperature(
  tk::real arho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real arhoE,
  tk::real alpha,
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
  tk::real damage ) const
// *************************************************************************
//! \brief Calculate material temperature from the material density, and
//!   material specific total energy
//! \param[in] arho Material partial density (alpha_k * rho_k)
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] arhoE Material total energy (alpha_k * rho_k * E_k)
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
//! \param[in] defgrad Material inverse deformation gradient tensor
//!   (g_k). Default is 0, so that for the single-material system,
//!   this argument can be left unspecified by the calling code
//! \return Material temperature using the LinearMieGruneisen EoS
// *************************************************************************
{
  // obtain elastic contribution to energy
  tk::real eps2;
  auto arhoEe = alpha*elasticEnergy(defgrad, eps2, damage);
  // obtain hydrodynamic internal energy
  auto arhoEi = arhoE - arhoEe - 0.5*arho*(u*u + v*v + w*w);

  const auto rho = arho/alpha;
  const auto t = (arhoEi/arho - hugoniotEnergy(rho))/m_cv;
  return t;
}

tk::real
LinearMieGruneisen::min_eff_pressure(
  tk::real min,
  tk::real arho,
  tk::real alpha ) const
// *************************************************************************
//! Compute the minimum allowed pressure
//! \param[in] min Numerical threshold above which pressure needs to be limited
//! \param[in] arho Material partial density (alpha_k * rho_k)
//! \param[in] alpha Material volume fraction
//! \return Minimum pressure allowed by physical constraints
// *************************************************************************
{
  // minimum pressure is constrained by zero soundspeed.
  auto al_eff = std::max( 1.0e-14, alpha );
  auto rho = arho/al_eff;
  auto gamma = gruneisen(rho);
  auto dgamma = dGruneisenDrho(rho);
  auto aph = al_eff*hugoniotPressure(rho);

  auto co = (gamma + rho*dgamma)/(gamma*arho);
  auto p_coeff = gamma/arho + co;
  auto p_const = dHugoniotPressureDrho(rho)
               - gamma*rho*dHugoniotEnergyDrho(rho)
               - co*aph;

  return -p_const/p_coeff + min;
}

tk::real
LinearMieGruneisen::elasticEnergy(
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
  tk::real& eps2,
  tk::real damage ) const
// *************************************************************************
//! \brief Calculate elastic contribution to material energy from the material
//!   density, and deformation gradient tensor
//! \param[in] defgrad Material inverse deformation gradient tensor
//! \param[in/out] eps2 Elastic shear distortion
//! \return Material elastic energy using the LinearMieGruneisen EoS
//! \details This function returns the material elastic energy, and also stores
//!   the elastic shear distortion for further use
// *************************************************************************
{
  // compute volume-preserving part of Right Cauchy-Green strain tensor
  auto Ct = tk::getIsochorRightCauchyGreen(defgrad);

  // compute elastic shear distortion
  eps2 = 0.5 * (Ct[0][0]+Ct[1][1]+Ct[2][2] - 3.0);

  // compute elastic energy
  auto rhoEe = m_mu * eps2;

  return rhoEe;
}

tk::real
LinearMieGruneisen::gruneisen( tk::real rho ) const
// *************************************************************************
//! Compute density-dependent Gruneisen coefficient
//! \param[in] rho Material density
//! \return Gruneisen coefficient
// *************************************************************************
{
  return m_gamma0 * std::pow(m_rho0/rho, m_alpha);
}

tk::real
LinearMieGruneisen::dGruneisenDrho( tk::real rho ) const
// *************************************************************************
//! Compute derivative of density-dependent Gruneisen coefficient wrt density
//! \param[in] rho Material density
//! \return Derivative of Gruneisen coefficient wrt density
// *************************************************************************
{
  return (- m_alpha * gruneisen(rho) / rho);
}

tk::real
LinearMieGruneisen::eta( tk::real rho ) const
// *************************************************************************
//! Compute compression measure
//! \param[in] rho Material density
//! \return Compression measure
// *************************************************************************
{
  return 1.0 - m_rho0/rho;
}

tk::real
LinearMieGruneisen::dEtaDrho( tk::real rho ) const
// *************************************************************************
//! Compute derivative of compression measure wrt density
//! \param[in] rho Material density
//! \return Derivative of compression measure wrt density
// *************************************************************************
{
  return m_rho0/(rho*rho);
}

tk::real
LinearMieGruneisen::hugoniotPressure( tk::real rho ) const
// *************************************************************************
//! Compute linear Us-Up Hugoniot pressure
//! \param[in] rho Material density
//! \return Hugoniot pressure
// *************************************************************************
{
  const auto et = eta(rho);
  const auto den = 1.0 - m_s1*et;
  return m_rho0*m_c0*m_c0*et/(den*den);
}

tk::real
LinearMieGruneisen::dHugoniotPressureDrho( tk::real rho ) const
// *************************************************************************
//! Compute derivative of Hugoniot pressure wrt density
//! \param[in] rho Material density
//! \return Derivative of Hugoniot pressure wrt density
// *************************************************************************
{
  const auto et = eta(rho);
  const auto den = 1.0 - m_s1*et;
  return m_rho0*m_c0*m_c0*(1.0 + m_s1*et)/(den*den*den)*
    dEtaDrho(rho);
}

tk::real
LinearMieGruneisen::hugoniotEnergy( tk::real rho ) const
// *************************************************************************
//! Compute Hugoniot internal energy
//! \param[in] rho Material density
//! \return Hugoniot internal energy
// *************************************************************************
{
  return 0.5*hugoniotPressure(rho)*eta(rho)/m_rho0;
}

tk::real
LinearMieGruneisen::dHugoniotEnergyDrho( tk::real rho ) const
// *************************************************************************
//! Compute derivative of Hugoniot internal energy wrt density
//! \param[in] rho Material density
//! \return Derivative of Hugoniot internal energy wrt density
// *************************************************************************
{
  return 0.5/m_rho0*( dHugoniotPressureDrho(rho)*eta(rho) +
    hugoniotPressure(rho)*dEtaDrho(rho) );
}
