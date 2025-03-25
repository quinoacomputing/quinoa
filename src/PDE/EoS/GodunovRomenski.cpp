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

// // Lapacke forward declarations
// extern "C" {

// using lapack_int = long;

// #define LAPACK_ROW_MAJOR 101

// lapack_int LAPACKE_dgeev(int, char, char, lapack_int, double*, lapack_int,
//   double*, double*, double*, lapack_int, double*, lapack_int );

// }

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
  tk::real ,
  tk::real ) const
// *************************************************************************
//! \brief Calculate density from the material pressure and temperature
//!   using the GodunovRomenski equation of state
//! \return Material density calculated using the density constaint
//!   rho = rho_0/det(F)
// *************************************************************************
{
  // punt by returning unstressed density for now
  return m_rho0;
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
// *************************************************************************
{
  // obtain elastic contribution to energy
  std::array< std::array< tk::real, 3 >, 3 > devH;
  auto arhoEe = alpha*elasticEnergy(defgrad, devH);
  // obtain cold compression contribution to energy
  auto rho = arho/alpha;
  auto arhoEc = alpha*coldcomprEnergy(rho);
  // obtain thermal contribution to energy
  auto arhoEt = arhoE - arhoEe - arhoEc - 0.5*arho*(u*u + v*v + w*w);

  // use Mie-Gruneisen form of Godunov-Romenski for pressure
  auto partpressure = alpha*coldcomprPressure(rho) - arhoEe + m_gamma*arhoEt;

  // check partial pressure divergence
  if (!std::isfinite(partpressure)) {
    std::cout << "Material-id:      " << imat << std::endl;
    std::cout << "Volume-fraction:  " << alpha << std::endl;
    std::cout << "Partial density:  " << arho << std::endl;
    Throw("Material-" + std::to_string(imat) +
      " has nan/inf partial pressure: " + std::to_string(partpressure) +
      ", material volume fraction: " + std::to_string(alpha));
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
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad ) const
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
//! \param[in] defgrad Material inverse deformation gradient tensor
//!   (g_k) with the first dimension aligned to direction in which
//!   wave speeds are required. Default is 0, so that for the single-material
//!   system, this argument can be left unspecified by the calling code
//! \return Material speed of sound using the GodunovRomenski EoS
// *************************************************************************
{
  tk::real a = 0.0;

  // Hydro contribution
  std::array< std::array< tk::real, 3 >, 3 > devH;
  tk::real rho = arho/alpha;
  auto p_se = -elasticEnergy(defgrad, devH);
  auto p_cc = coldcomprPressure(rho);
  auto rrho0a = std::pow(rho/m_rho0, m_alpha);
  a += std::max( 1.0e-15,
    m_K0/(m_rho0*m_alpha) * ((2.0*m_alpha+1.0)*(rrho0a*rrho0a) - (m_alpha+1.0)*rrho0a)
    + (m_gamma+1.0) * (apr - p_cc - p_se)/arho );

  // Shear contribution
  a += (4.0/3.0) * alpha * m_mu / arho;

  // Compute square root
  a = std::sqrt(a);

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
  tk::real a = std::sqrt(alpha*m_mu/arho);

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
GodunovRomenski::totalenergy(
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
//!   g_k. Default is 0, so that for the single-material system,
//!   this argument can be left unspecified by the calling code
//! \return Material specific total energy using the GodunovRomenski EoS
// *************************************************************************
{
  // obtain thermal and kinetic energy
  std::array< std::array< tk::real, 3 >, 3 > devH;
  auto p_se = -elasticEnergy(defgrad, devH);
  auto pt = pr - coldcomprPressure(rho) - p_se;
  auto rhoEh = pt/m_gamma + 0.5*rho*(u*u + v*v + w*w);
  // obtain elastic contribution to energy
  auto rhoEe = -p_se;
  auto rhoEc = coldcomprEnergy(rho);

  return (rhoEh + rhoEe + rhoEc);
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
  // minimum pressure is constrained by zero soundspeed.
  auto rho = arho/alpha;
  auto rrho0a = std::pow(rho/m_rho0, m_alpha);
  return min
    - rho/(m_gamma+1.0)
      * m_K0/(m_rho0*m_alpha) * ((2.0*m_alpha+1.0)*(rrho0a*rrho0a) - (m_alpha+1.0)*rrho0a)
    + coldcomprPressure(rho);
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
GodunovRomenski::coldcomprPressure( tk::real rho ) const
// *************************************************************************
//! \brief Calculate cold-compression contribution to material pressure from the
//!   material density
//! \param[in] rho Material density
//! \return Material cold compression pressure using the GodunovRomenski EoS
// *************************************************************************
{
  auto rrho0a = std::pow(rho/m_rho0, m_alpha);
  return (m_K0/m_alpha * (rrho0a*rho/m_rho0) * (rrho0a-1.0));
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
  // Compute deviatoric part of Hencky tensor
  devH = tk::getDevHencky(defgrad);

  // Compute elastic energy
  tk::real rhoEe = 0.0;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      rhoEe += m_mu*devH[i][j]*devH[i][j];

  return rhoEe;
}
