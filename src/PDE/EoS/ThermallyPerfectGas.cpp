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
  tk::real R,
  std::vector< std::vector< tk::real > > cp_coeff,
  std::vector< tk::real > t_range,
  tk::real dH_ref) :
  m_R(R),
  m_cp_coeff(cp_coeff),
  m_t_range(t_range),
  m_dH_ref(dH_ref)
// *************************************************************************
//  Constructor
//! \param[in] R gas constant
//! \param[in] cp_coeff NASA Glenn polynomials coefficients for cp fit
//! \param[in] t_range temperature range where polynomial coeffs are valid
//! \param[in] dH_ref reference enthalpy, h(t = 298.15 K) - h(t = 0 K)
// *************************************************************************
{ }

[[noreturn]] tk::real
ThermallyPerfectGas::density(
  tk::real ,
  tk::real ) const
// *************************************************************************
//! \brief Calculate density from the material pressure and temperature 
//!   using the stiffened-gas equation of state
//! \param[in] pr Material pressure
//! \param[in] temp Material temperature
//! \return Material density calculated using the stiffened-gas EoS
// *************************************************************************
{
  Throw("Direct call to TPG density should not occur. Use Mixture class.");
}

[[noreturn]] tk::real
ThermallyPerfectGas::pressure(
  tk::real ,
  tk::real ,
  tk::real ,
  tk::real ,
  tk::real ,
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
  Throw("Direct call to TPG pressure should not occur. Use Mixture class.");
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

[[noreturn]] tk::real
ThermallyPerfectGas::soundspeed(
  tk::real ,
  tk::real ,
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
  Throw("Direct call to TPG soundspeed should not occur. Use Mixture class.");
}

[[noreturn]] tk::real
ThermallyPerfectGas::totalenergy(
  tk::real ,
  tk::real ,
  tk::real ,
  tk::real ,
  tk::real ,
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
  Throw("Direct call to TPG totalenergy should not occur. Use Mixture class.");
}

[[noreturn]] tk::real
ThermallyPerfectGas::temperature(
  tk::real ,
  tk::real ,
  tk::real ,
  tk::real ,
  tk::real ,
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
  Throw("Direct call to TPG temperature should not occur. Use Mixture class.");
}

tk::real
ThermallyPerfectGas::internalenergy(tk::real temp) const
// *************************************************************************
//! \brief Calculate species internal energy
//! \param[in] temp Temperature
//! \return Species internal energy using the thermally perfect gas EoS
// *************************************************************************
{
  auto R = m_R;
  tk::real h = calc_h(temp) * R * temp + m_dH_ref;
  return h - R * temp;
}

tk::real
ThermallyPerfectGas::cv(tk::real temp) const
// *************************************************************************
//! \brief Calculate species specific heat (constant volume)
//! \param[in] temp Temperature
//! \return Species specific heat using the thermally perfect gas EoS
// *************************************************************************
{
  auto R = m_R;
  tk::real cp = calc_cp(temp) * R;
  return cp - R;
}
