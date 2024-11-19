// *****************************************************************************
/*!
  \file      src/PDE/EoS/ThermallyPerfectGas.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Thermally perfect gas equation of state
  \details   This file defines functions for the thermally perfect gas equation
             of state for the compressible flow equations.
*/
// *****************************************************************************

#include <cmath>
#include <iostream>
#include "EoS/ThermallyPerfectGas.hpp"

using inciter::ThermallyPerfectGas;

ThermallyPerfectGas::ThermallyPerfectGas(
  tk::real gamma,
  tk::real R,
  std::vector< tk::real > cp_coeff) :
  m_gamma(gamma),
  m_R(R),
  m_cp_coeff(cp_coeff)
// *************************************************************************
//  Constructor
//! \param[in] gamma Ratio of specific heats
//! \param[in] R gas constant
//! \param[in] cp_coeff NASA Glenn polynomials coefficients for cp fit
// *************************************************************************
{ }

tk::real
ThermallyPerfectGas::density(
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
  tk::real R = m_R;

  tk::real rho = pr / (R * temp);
  return rho;
}

tk::real
ThermallyPerfectGas::pressure(
  tk::real rho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real rhoE,
  tk::real,
  std::size_t,
  const std::array< std::array< tk::real, 3 >, 3 >& ) const
// *************************************************************************
//! \brief Calculate pressure from the material density, momentum and total
//!   energy using the thermally perfect gas equation of state
//! \param[in] rho density
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] rhoE total energy
//! \return Pressure calculated using the thermally perfect gas EOS
// *************************************************************************
{
  tk::real R = m_R;

  tk::real temp = temperature(rho, u, v, w, rhoE);
  tk::real pres = rho * R * temp;

  return pres;
}

std::array< std::array< tk::real, 3 >, 3 >
ThermallyPerfectGas::CauchyStress(
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
ThermallyPerfectGas::soundspeed(
  tk::real rho,
  tk::real pr,
  tk::real,
  std::size_t,
  const std::array< std::array< tk::real, 3 >, 3 >&,
  const std::array< tk::real, 3 >& ) const
// *************************************************************************
//! Calculate speed of sound from the material density and material pressure
//! \param[in] rho density
//! \param[in] pr pressure
//! \return Material speed of sound using the ideal gas EoS
// *************************************************************************
{
  auto g = m_gamma;

  tk::real a = std::sqrt( g * pr / rho );

  return a;
}

tk::real
ThermallyPerfectGas::totalenergy(
  tk::real rho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real pr,
  const std::array< std::array< tk::real, 3 >, 3 >& ) const
// *************************************************************************
//! \brief Calculate material specific total energy from the material
//!   density, momentum and material pressure
//! \param[in] rho density
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] pr pressure
//! \return specific total energy using the thermally perfect gas EoS
// *************************************************************************
{
  auto R = m_R;

  tk::real temp = pr / (rho * R);
  tk::real e = R * (-m_cp_coeff[0] * pow(temp, -1) +
      m_cp_coeff[1] * log(temp) + (m_cp_coeff[2] - 1) * temp +
      m_cp_coeff[3] * pow(temp, 2) / 2 +
      m_cp_coeff[4] * pow(temp, 3) / 3 +
      m_cp_coeff[5] * pow(temp, 4) / 4 +
      m_cp_coeff[6] * pow(temp, 5) / 5 + m_cp_coeff[7]);

  return (rho * e + 0.5 * rho * (u*u + v*v + w*w));
}

tk::real
ThermallyPerfectGas::temperature(
  tk::real rho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real rhoE,
  tk::real,
  const std::array< std::array< tk::real, 3 >, 3 >& ) const
// *************************************************************************
//! \brief Calculate material temperature from the material density
//! \param[in] rho density
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] rhoE total energy
//! \return Material temperature using the thermally perfect gas EoS
// *************************************************************************
{
  auto R = m_R;

  // Solve for internal energy
  tk::real e = rhoE / rho - 0.5 * (u*u + v*v + w*w);

  // Solve for temperature
  tk::real temp = 1000; // Starting guess
  tk::real tol = 1e-8; // Stopping condition
  tk::real err = 1e8;
  while (err > tol) {
    tk::real f_T = R * (-m_cp_coeff[0] * pow(temp, -1) +
      m_cp_coeff[1] * log(temp) + (m_cp_coeff[2] - 1) * temp +
      m_cp_coeff[3] * pow(temp, 2) / 2 +
      m_cp_coeff[4] * pow(temp, 3) / 3 +
      m_cp_coeff[5] * pow(temp, 4) / 4 +
      m_cp_coeff[6] * pow(temp, 5) / 5 + m_cp_coeff[7]) - e;

    err = abs(f_T);

    tk::real fp_T = 0;
    tk::real power = -2;
    for (std::size_t k=0; k<m_cp_coeff.size()-1; ++k)
    {
      fp_T += m_cp_coeff[k] * pow(temp, power);
      if (k == 2) fp_T += -1;
      power += 1;
    }
    fp_T = fp_T * R;

    temp = temp - f_T / fp_T;
  }

  return temp;
}
