// *****************************************************************************
/*!
  \file      src/PDE/EoS/WilkinsAluminum.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Wilkins equation of state for aluminum
  \details   This file defines functions for the Wilkins equation of
             state for solids and a hydro EoS for aluminum. These functions were
             taken from Example 4 of Barton, Philip T. "An interface-capturing
             Godunov method for the simulation of compressible solid-fluid
             problems." Journal of Computational Physics 390 (2019): 25-50.
*/
// *****************************************************************************

#include <cmath>
#include <iostream>
#include "Vector.hpp"
#include "EoS/WilkinsAluminum.hpp"
#include "EoS/GetMatProp.hpp"

// // Lapacke forward declarations
// extern "C" {

// using lapack_int = long;

// #define LAPACK_ROW_MAJOR 101

// lapack_int LAPACKE_dgeev(int, char, char, lapack_int, double*, lapack_int,
//   double*, double*, double*, lapack_int, double*, lapack_int );

// }

static const tk::real e1 = -13.0e+09;
static const tk::real e2 = 20.0e+09;
static const tk::real e3 = 52.0e+09;
static const tk::real e4 = -59.0e+09;
static const tk::real e5 = 151.0e+09;

using inciter::WilkinsAluminum;

WilkinsAluminum::WilkinsAluminum(
  tk::real gamma,
  tk::real cv,
  tk::real mu ) :
  m_gamma(gamma),
  m_cv(cv),
  m_mu(mu)
// *************************************************************************
//  Constructor
//! \param[in] gamma Ratio of specific heats
//! \param[in] cv Specific heat at constant volume
//! \param[in] mu Constant shear modulus
// *************************************************************************
{
  // Since this is only for aluminum we hard set rho0
  m_rho0 = 2700.0;
}

void
WilkinsAluminum::setRho0( tk::real rho0 )
// *************************************************************************
//  Set rho0 EOS parameter; i.e. the initial density
//! \param[in] rho0 Initial material density that needs to be stored
// *************************************************************************
{
  m_rho0 = rho0;
}

tk::real
WilkinsAluminum::density(
  tk::real pr,
  tk::real ) const
// *************************************************************************
//! \brief Calculate density from the material pressure and temperature
//!   using the WilkinsAluminum equation of state
//! \param[in] pr Material pressure
// //! \param[in] temp Material temperature
//! \return Material density calculated using the WilkinsAluminum EoS
// *************************************************************************
{
  tk::real rho0 = m_rho0;
  // Quick Newton
  tk::real rho = rho0;
  std::size_t maxiter = 50;
  tk::real tol = 1.0e-04;
  tk::real err = tol + 1;
  for (std::size_t iter=0; iter<maxiter; ++iter)
  {
    tk::real p = 2*e2*std::pow(rho/rho0,3.0)
               + e3*std::pow(rho/rho0,2.0)
               - e5*rho/rho0 - e4 - pr;
    tk::real dpdrho = 6*e2*std::pow(rho/rho0,2.0)/rho0
                    + 2*e3*rho/(rho0*rho0) - e5/rho0;
    tk::real delta = p/dpdrho;
    rho -= delta;
    err = std::sqrt(std::pow(p,2.0));
    if (err < tol) break;
  }
  return rho;
}

tk::real
WilkinsAluminum::pressure(
  tk::real arho,
  tk::real,
  tk::real,
  tk::real,
  tk::real,
  tk::real alpha,
  std::size_t imat,
  const std::array< std::array< tk::real, 3 >, 3 >& ) const
// *************************************************************************
//! \brief Calculate pressure from the material density, momentum, total energy
//!   and the inverse deformation gradient tensor using the WilkinsAluminum
//!   equation of state
//! \param[in] arho Material partial density (alpha_k * rho_k)
// //! \param[in] u X-velocity
// //! \param[in] v Y-velocity
// //! \param[in] w Z-velocity
// //! \param[in] arhoE Material total energy (alpha_k * rho_k * E_k)
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified
//!   by the calling code
// //! \param[in] defgrad Material inverse deformation gradient tensor
// //!   (g_k). Default is 0, so that for the single-material system,
// //!   this argument can be left unspecified by the calling code
//! \return Material partial pressure (alpha_k * p_k) calculated using the
//!   WilkinsAluminum EoS
// *************************************************************************
{
  tk::real rho0 = m_rho0;
  tk::real rho = arho/alpha;
  tk::real partpressure = alpha*(2*e2*std::pow(rho/rho0,3.0)
                                 + e3*std::pow(rho/rho0,2.0)
                                 - e5*rho/rho0 - e4 );

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
WilkinsAluminum::CauchyStress(
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
//!   WilkinsAluminum equation of state
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
// //! \param[in] imat Material-id who's EoS is required. Default is 0, so that
// //!   for the single-material system, this argument can be left unspecified
// //!   by the calling code
//! \param[in] defgrad Material inverse deformation gradient tensor (g_k).
//! \return Material Cauchy stress tensor (alpha_k * sigma_k) calculated using
//!   the WilkinsAluminum EoS
// *************************************************************************
{
  std::array< std::array< tk::real, 3 >, 3 > asig{{{0,0,0}, {0,0,0}, {0,0,0}}};

  // obtain elastic contribution to energy and substract it from pmean
  std::array< std::array< tk::real, 3 >, 3 > devH;

  // p_mean
  auto pmean = - alpha * elasticEnergy(defgrad, devH);

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
WilkinsAluminum::soundspeed(
  tk::real arho,
  tk::real apr,
  tk::real alpha,
  std::size_t imat,
  const std::array< std::array< tk::real, 3 >, 3 >& /*defgrad*/ ) const
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
//! \return Material speed of sound using the WilkinsAluminum EoS
// *************************************************************************
{
  tk::real a = 0.0;

  // Hydro contribution
  tk::real rho0 = m_rho0;
  tk::real rho = arho/alpha;
  a += std::max( 1.0e-15, 6*e2*std::pow(rho/rho0,2.0)/rho0
                 + 2*e3*rho/(rho0*rho0) - e5/rho0 );

  // Shear contribution
  a += (4.0/3.0) * m_mu / (arho/alpha);

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
WilkinsAluminum::shearspeed(
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
WilkinsAluminum::totalenergy(
  tk::real arho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real,
  tk::real alpha,
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad ) const
// *************************************************************************
//! \brief Calculate material specific total energy from the material
//!   density, momentum and material pressure
//! \param[in] arho Material partial density
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
// //! \param[in] apr Material partial pressure
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
//! \param[in] defgrad Material inverse deformation gradient tensor
//!   g_k. Default is 0, so that for the single-material system,
//!   this argument can be left unspecified by the calling code
//! \return Material specific total energy using the WilkinsAluminum EoS
// *************************************************************************
{
  // obtain hydro contribution to energy
  tk::real rho0 = m_rho0;
  tk::real rho = arho/alpha;
  tk::real rhoEh = (e1+e2*std::pow(rho/rho0,2.0)+e3*(rho/rho0)
                    +e4*std::pow(rho/rho0,-1.0)-e5*std::log(rho/rho0))/rho0
                   + 0.5*rho*(u*u + v*v + w*w);
  // obtain elastic contribution to energy
  std::array< std::array< tk::real, 3 >, 3 > devH;
  tk::real rhoEe = elasticEnergy(defgrad, devH);

  return alpha*(rhoEh + rhoEe);
}

tk::real
WilkinsAluminum::temperature(
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
//! \return Material temperature using the WilkinsAluminum EoS
// *************************************************************************
{
  // Temperature does not directly contribute to energy
  // So we just set a value.
  tk::real t = 300.0;

  return t;
}

tk::real
WilkinsAluminum::min_eff_pressure(
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
  return min;
}

tk::real
WilkinsAluminum::elasticEnergy(
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
  std::array< std::array< tk::real, 3 >, 3 >& devH ) const
// *************************************************************************
//! \brief Calculate elastic contribution to material energy from the material
//!   density, and deformation gradient tensor
//! \param[in] defgrad Material inverse deformation gradient tensor
//! \param[in/out] devH Deviatoric part of the Hensky tensor
//! \return Material elastic energy using the WilkinsAluminum EoS
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
