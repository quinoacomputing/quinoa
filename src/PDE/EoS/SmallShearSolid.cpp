// *****************************************************************************
/*!
  \file      src/PDE/EoS/SmallShearSolid.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Small shear strain equation of state for solids
  \details   This file defines functions for the SmallShearSolid equation of
             state for the compressible flow equations. These functions are
             taken from Plohr, J. N., & Plohr, B. J. (2005). Linearized analysis
             of Richtmyer–Meshkov flow for elastic materials. Journal of Fluid
             Mechanics, 537, 55-89. The SmallShearSolid EOS uses a small-shear
             approximation for the elastic contribution, and a stiffened gas EOS
             for the hydrodynamic contribution of the internal energy.
*/
// *****************************************************************************

#include <cmath>
#include <iostream>
#include <lapacke.h>
#include "Vector.hpp"
#include "EoS/SmallShearSolid.hpp"

using inciter::SmallShearSolid;

SmallShearSolid::SmallShearSolid(
  tk::real gamma,
  tk::real pstiff,
  tk::real cv,
  tk::real mu ) :
  m_gamma(gamma),
  m_pstiff(pstiff),
  m_cv(cv),
  m_mu(mu)
// *************************************************************************
//  Constructor
//! \param[in] gamma Ratio of specific heats
//! \param[in] pstiff Stiffness pressure term
//! \param[in] cv Specific heat at constant volume
//! \param[in] mu Constant shear modulus
// *************************************************************************
{ }

tk::real
SmallShearSolid::density(
  tk::real pr,
  tk::real temp ) const
// *************************************************************************
//! \brief Calculate density from the material pressure and temperature 
//!   using the SmallShearSolid equation of state
//! \param[in] pr Material pressure
//! \param[in] temp Material temperature
//! \return Material density calculated using the SmallShearSolid EoS
// *************************************************************************
{
  tk::real g = m_gamma;
  tk::real p_c = m_pstiff;
  tk::real c_v = m_cv;

  tk::real rho = (pr + p_c) / ((g-1.0) * c_v * temp);
  return rho;
}

tk::real
SmallShearSolid::pressure(
  tk::real arho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real arhoE,
  tk::real alpha,
  std::size_t imat,
  const std::array< std::array< tk::real, 3 >, 3 >& adefgrad ) const
// *************************************************************************
//! \brief Calculate pressure from the material density, momentum, total energy
//!   and the inverse deformation gradient tensor using the SmallShearSolid
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
//! \param[in] adefgrad Material inverse deformation gradient tensor
//!   (alpha_k * g_k). Default is 0, so that for the single-material system,
//!   this argument can be left unspecified by the calling code
//! \return Material partial pressure (alpha_k * p_k) calculated using the
//!   SmallShearSolid EoS
// *************************************************************************
{
  // deformation gradient: 
  auto defgrad = adefgrad;
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      defgrad[i][j] /= alpha;
  }

  // obtain elastic contribution to energy
  tk::real eps2;
  auto arhoEe = alpha*elasticEnergy(arho/alpha, defgrad, eps2);
  // obtain hydro contribution to energy
  auto arhoEh = arhoE - arhoEe;

  // use stiffened gas eos to get pressure
  tk::real partpressure = (arhoEh - 0.5 * arho * (u*u + v*v + w*w) -
    alpha*m_pstiff) * (m_gamma-1.0) - alpha*m_pstiff;

  // check partial pressure divergence
  if (!std::isfinite(partpressure)) {
    std::cout << "Material-id:      " << imat << std::endl;
    std::cout << "Volume-fraction:  " << alpha << std::endl;
    std::cout << "Partial density:  " << arho << std::endl;
    std::cout << "Total energy:     " << arhoE << std::endl;
    std::cout << "Hydro energy:     " << arhoEh << std::endl;
    std::cout << "Velocity:         " << u << ", " << v << ", " << w
      << std::endl;
    Throw("Material-" + std::to_string(imat) +
      " has nan/inf partial pressure: " + std::to_string(partpressure) +
      ", material volume fraction: " + std::to_string(alpha));
  }

  return partpressure;
}

std::array< std::array< tk::real, 3 >, 3 >
SmallShearSolid::cauchyStressTensor(
  tk::real arho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real arhoE,
  tk::real alpha,
  std::size_t imat,
  const std::array< std::array< tk::real, 3 >, 3 >& adefgrad ) const
// *************************************************************************
//! \brief Calculate the Cauchy stress tensor from the material density,
//!   momentum, total energy, and inverse deformation gradient tensor using the
//!   SmallShearSolid equation of state
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
//! \param[in] adefgrad Material inverse deformation gradient tensor
//!   (alpha_k * g_k). Default is 0, so that for the single-material system,
//!   this argument can be left unspecified by the calling code
//! \return Material partial pressure (alpha_k * p_k) calculated using the
//!   SmallShearSolid EoS
// *************************************************************************
{
  std::array< std::array< tk::real, 3 >, 3 > sig{{{0,0,0}, {0,0,0}, {0,0,0}}};

  // deformation gradient: 
  auto defgrad = adefgrad;
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      defgrad[i][j] /= alpha;
  }

  // obtain elastic contribution to energy
  tk::real eps2;
  auto arhoEe = alpha*elasticEnergy(arho/alpha, defgrad, eps2);
  // obtain hydro contribution to energy
  auto arhoEh = arhoE - arhoEe;

  // use stiffened gas eos to get pressure
  tk::real partpressure = (arhoEh - 0.5 * arho * (u*u + v*v + w*w) -
    alpha*m_pstiff) * (m_gamma-1.0) - alpha*m_pstiff;

  // check partial pressure divergence
  if (!std::isfinite(partpressure)) {
    std::cout << "Material-id:      " << imat << std::endl;
    std::cout << "Volume-fraction:  " << alpha << std::endl;
    std::cout << "Partial density:  " << arho << std::endl;
    std::cout << "Total energy:     " << arhoE << std::endl;
    std::cout << "Hydro energy:     " << arhoEh << std::endl;
    std::cout << "Velocity:         " << u << ", " << v << ", " << w
      << std::endl;
    Throw("Material-" + std::to_string(imat) +
      " has nan/inf partial pressure: " + std::to_string(partpressure) +
      ", material volume fraction: " + std::to_string(alpha));
  }

  // p_mean
  auto pmean = partpressure/alpha - m_mu * eps2;

  // Volumetric component of Cauchy stress tensor
  sig[0][0] = -pmean;
  sig[1][1] = -pmean;
  sig[2][2] = -pmean;

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
      sig[i][j] += devbt[i][j];
  }

  return sig;
}

tk::real
SmallShearSolid::soundspeed(
  tk::real arho,
  tk::real apr,
  tk::real alpha,
  std::size_t imat ) const
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
//! \return Material speed of sound using the SmallShearSolid EoS
// *************************************************************************
{
  auto g = m_gamma;
  auto p_c = m_pstiff;

  auto p_eff = std::max( 1.0e-15, apr+(alpha*p_c) );

  tk::real a = std::sqrt( g * p_eff / arho );

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
SmallShearSolid::totalenergy(
  tk::real rho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real pr,
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad ) const
// *************************************************************************
//! \brief Calculate material specific total energy from the material
//!   density, momentum and material pressure
//! \param[in] rho Material density
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] pr Material pressure
//! \param[in] defgrad Material inverse deformation gradient tensor
//!   (alpha_k * g_k). Default is 0, so that for the single-material system,
//!   this argument can be left unspecified by the calling code
//! \return Material specific total energy using the SmallShearSolid EoS
// *************************************************************************
{
  // obtain hydro contribution to energy
  tk::real rhoEh = (pr + m_pstiff) / (m_gamma-1.0) + 0.5 * rho *
    (u*u + v*v + w*w) + m_pstiff;
  // obtain elastic contribution to energy
  tk::real eps2;
  tk::real rhoEe = elasticEnergy(rho, defgrad, eps2);

  return (rhoEh + rhoEe);
}

tk::real
SmallShearSolid::temperature(
  tk::real arho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real arhoE,
  tk::real alpha ) const
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
//! \return Material temperature using the SmallShearSolid EoS
// *************************************************************************
{
  auto c_v = m_cv;
  auto p_c = m_pstiff;

  tk::real t = (arhoE - 0.5 * arho * (u*u + v*v + w*w) - alpha*p_c) 
               / (arho*c_v);
  return t;
}

tk::real
SmallShearSolid::min_eff_pressure(
  tk::real min,
  tk::real,
  tk::real ) const
// *************************************************************************
//! Compute the minimum allowed pressure
//! \param[in] min Numerical threshold above which pressure needs to be limited
//! \return Minimum pressure allowed by physical constraints
// *************************************************************************
{
  // minimum pressure is constrained by zero soundspeed.
  return (min - m_pstiff);
}

tk::real
SmallShearSolid::elasticEnergy(
  tk::real rho,
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
  tk::real& eps2 ) const
// *************************************************************************
//! \brief Calculate elastic contribution to material energy from the material
//!   density, and deformation gradient tensor
//! \param[in] rho Material density
//! \param[in] defgrad Material inverse deformation gradient tensor
//! \param[in/out] eps2 Elastic shear distortion
//! \return Material elastic energy using the SmallShearSolid EoS
//! \details This function returns the material elastic energy, and also stores
//!   the elastic shear distortion for further use
// *************************************************************************
{
  // compute Right Cauchy-Green strain tensor
  auto Ct = tk::getRightCauchyGreen(defgrad);
  auto detC = std::pow(tk::determinant(Ct), 1.0/3.0);
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      Ct[i][j] /= detC;
  }

  // compute elastic shear distortion
  eps2 = 0.5 * (Ct[0][0]+Ct[1][1]+Ct[2][2] - 3.0);

  // compute elastic energy
  auto rhoEe = rho * m_mu * eps2;

  return rhoEe;
}
