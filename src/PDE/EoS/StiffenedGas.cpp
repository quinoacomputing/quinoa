// *****************************************************************************
/*!
  \file      src/PDE/EoS/StiffenedGas.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Stiffened-gas equation of state
  \details   This file defines functions for the stiffened gas equation of
             state for the compressible flow equations.
*/
// *****************************************************************************

#include <cmath>
#include <iostream>
#include "EoS/StiffenedGas.hpp"

using inciter::StiffenedGas;

StiffenedGas::StiffenedGas(
  tk::real gamma,
  tk::real pstiff,
  tk::real cv ) :
  m_gamma(gamma),
  m_pstiff(pstiff),
  m_cv(cv)
// *************************************************************************
//  Constructor
//! \param[in] gamma Ratio of specific heats
//! \param[in] pstiff Stiffened pressure term
//! \param[in] cv Specific heat at constant volume
// *************************************************************************
{ }

tk::real
StiffenedGas::density(
  tk::real pr,
  tk::real temp ) const
// *************************************************************************
//! \brief Calculate density from the material pressure and temperature 
//!   using the stiffened-gas equation of state
//! \param[in] pr Material pressure
//! \param[in] temp Material temperature
//! \return Material density calculated using the stiffened-gas EoS
// *************************************************************************
{
  tk::real g = m_gamma;
  tk::real p_c = m_pstiff;
  tk::real c_v = m_cv;

  tk::real rho = (pr + p_c) / ((g-1.0) * c_v * temp);
  return rho;
}

tk::real
StiffenedGas::pressure(
  tk::real arho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real arhoE,
  tk::real alpha,
  std::size_t imat,
  const std::array< std::array< tk::real, 3 >, 3 >& ) const
// *************************************************************************
//! \brief Calculate pressure from the material density, momentum and total
//!   energy using the stiffened-gas equation of state
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
//! \return Material partial pressure (alpha_k * p_k) calculated using the
//!   stiffened-gas EoS
// *************************************************************************
{
  tk::real g = m_gamma;
  tk::real p_c = m_pstiff;

  tk::real partpressure = (arhoE - 0.5 * arho * (u*u + v*v + w*w) -
    alpha*p_c) * (g-1.0) - alpha*p_c;

  // check partial pressure divergence
  if (!std::isfinite(partpressure)) {
    std::cout << "Material-id:      " << imat << std::endl;
    std::cout << "Volume-fraction:  " << alpha << std::endl;
    std::cout << "Partial density:  " << arho << std::endl;
    std::cout << "Total energy:     " << arhoE << std::endl;
    std::cout << "Velocity:         " << u << ", " << v << ", " << w
      << std::endl;
    Throw("Material-" + std::to_string(imat) +
      " has nan/inf partial pressure: " + std::to_string(partpressure) +
      ", material volume fraction: " + std::to_string(alpha));
  }

  return partpressure;
}

std::array< std::array< tk::real, 3 >, 3 >
StiffenedGas::CauchyStress(
  tk::real,
  tk::real,
  tk::real,
  tk::real,
  tk::real,
  tk::real,
  std::size_t,
  const std::array< std::array< tk::real, 3 >, 3 >& ) const
// *************************************************************************
//! \brief Calculate the Cauchy stress tensor from the material density,
//!   momentum, and total energy
//! \return Material Cauchy stress tensor (alpha_k * sigma_k)
// *************************************************************************
{
  std::array< std::array< tk::real, 3 >, 3 > asig{{{0,0,0}, {0,0,0}, {0,0,0}}};

  // No elastic contribution

  return asig;
}

tk::real
StiffenedGas::soundspeed(
  tk::real arho,
  tk::real apr,
  tk::real alpha,
  std::size_t imat,
  const std::array< std::array< tk::real, 3 >, 3 >&,
  const std::array< tk::real, 3 >&,
  const std::array< tk::real, 3 >& ) const
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
//! \return Material speed of sound using the stiffened-gas EoS
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
StiffenedGas::totalenergy(
  tk::real rho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real pr,
  const std::array< std::array< tk::real, 3 >, 3 >& ) const
// *************************************************************************
//! \brief Calculate material specific total energy from the material
//!   density, momentum and material pressure
//! \param[in] rho Material density
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] pr Material pressure
//! \return Material specific total energy using the stiffened-gas EoS
// *************************************************************************
{
  auto g = m_gamma;
  auto p_c = m_pstiff;

  tk::real rhoE = (pr + p_c) / (g-1.0) + 0.5 * rho * (u*u + v*v + w*w) + p_c;
  return rhoE;
}

tk::real
StiffenedGas::temperature(
  tk::real arho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real arhoE,
  tk::real alpha,
  const std::array< std::array< tk::real, 3 >, 3 >& ) const
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
//! \return Material temperature using the stiffened-gas EoS
// *************************************************************************
{
  auto c_v = m_cv;
  auto p_c = m_pstiff;

  tk::real t = (arhoE - 0.5 * arho * (u*u + v*v + w*w) - alpha*p_c) 
               / (arho*c_v);
  return t;
}

tk::real
StiffenedGas::min_eff_pressure(
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
