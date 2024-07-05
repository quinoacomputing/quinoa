// *****************************************************************************
/*!
  \file      src/PDE/EoS/GodunovRomenskiSolid.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Godunov-Romenski equation of state for solids
  \details   This file defines functions for the Godunov-Romenski equation of
             state for solids. These function were mostly taken from Barton,
             Philip T. "An interface-capturing Godunov method for the simulation
             of compressible solid-fluid problems." Journal of Computational
             Physics 390 (2019): 25-50. The elastic energy and stress is
             obtained from a the deviatoric part of the Hencky strain, while the
             hydrodynamics contributions are obtained from a stiffened gas EOS.
*/
// *****************************************************************************

#include <cmath>
#include <iostream>
#include "Vector.hpp"
#include "EoS/GodunovRomenskiSolid.hpp"

// Lapacke forward declarations
extern "C" {

using lapack_int = long;

#define LAPACK_ROW_MAJOR 101

lapack_int LAPACKE_dgeev(int, char, char, lapack_int, double*, lapack_int,
  double*, double*, double*, lapack_int, double*, lapack_int );

}

using inciter::GodunovRomenskiSolid;

GodunovRomenskiSolid::GodunovRomenskiSolid(
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
GodunovRomenskiSolid::density(
  tk::real pr,
  tk::real temp ) const
// *************************************************************************
//! \brief Calculate density from the material pressure and temperature 
//!   using the GodunovRomenskiSolid equation of state
//! \param[in] pr Material pressure
//! \param[in] temp Material temperature
//! \return Material density calculated using the GodunovRomenskiSolid EoS
// *************************************************************************
{
  tk::real g = m_gamma;
  tk::real p_c = m_pstiff;
  tk::real c_v = m_cv;

  tk::real rho = (pr + p_c) / ((g-1.0) * c_v * temp);
  return rho;
}

tk::real
GodunovRomenskiSolid::pressure(
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
//!   and the inverse deformation gradient tensor using the GodunovRomenskiSolid
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
//!   GodunovRomenskiSolid EoS
// *************************************************************************
{
  // obtain elastic contribution to energy
  tk::real eps2;
  auto arhoEe = alpha*elasticEnergy(defgrad, eps2);
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
    std::cout << "det(defgrad):     " << tk::determinant(defgrad) << std::endl;
    std::cout << "Velocity:         " << u << ", " << v << ", " << w
      << std::endl;
    Throw("Material-" + std::to_string(imat) +
      " has nan/inf partial pressure: " + std::to_string(partpressure) +
      ", material volume fraction: " + std::to_string(alpha));
  }

  return partpressure;
}

std::array< std::array< tk::real, 3 >, 3 >
GodunovRomenskiSolid::CauchyStress(
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
//!   GodunovRomenskiSolid equation of state
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
// //! \param[in] imat Material-id who's EoS is required. Default is 0, so that
// //!   for the single-material system, this argument can be left unspecified
// //!   by the calling code
//! \param[in] defgrad Material inverse deformation gradient tensor (g_k).
//! \return Material Cauchy stress tensor (alpha_k * sigma_k) calculated using
//!   the GodunovRomenskiSolid EoS
// *************************************************************************
{
  std::array< std::array< tk::real, 3 >, 3 > asig{{{0,0,0}, {0,0,0}, {0,0,0}}};

  // obtain elastic contribution to energy
  tk::real eps2;
  elasticEnergy(defgrad, eps2);

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
  auto devH = tk::getDevHencky(defgrad);
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      asig[i][j] += 2.0*m_mu*alpha*devH[i][j];
  }

  return asig;
}

tk::real
GodunovRomenskiSolid::soundspeed(
  tk::real arho,
  tk::real apr,
  tk::real alpha,
  std::size_t imat,
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
  const std::array< tk::real, 3 >& /*adefgradn*/,
  const std::array< tk::real, 3 >& /*asigman*/ ) const
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
//  //! \param[in] adefgradn Material inverse deformation gradient tensor in
//  //!   direction of vector n (alpha_k * g_ij * n_j). Default is 0, so that for
//  //!   the single-material system, this argument can be left unspecified by the
//  //!   calling code
//  //! \param[in] asigman Material traction vector in normal direction
//  //!   (alpha * sigma_ij * n_j ). Default is 0, so that for the single-material
//  //!   system, this argument can be left unspecified by the calling code
//! \return Material speed of sound using the GodunovRomenskiSolid EoS
// *************************************************************************
{
  // deformation gradient
  tk::real g__11 = defgrad[0][0];
  tk::real g__12 = defgrad[0][1];
  tk::real g__13 = defgrad[0][2];
  tk::real g__21 = defgrad[1][0];
  tk::real g__22 = defgrad[1][1];
  tk::real g__23 = defgrad[1][2];
  tk::real g__31 = defgrad[2][0];
  tk::real g__32 = defgrad[2][1];
  tk::real g__33 = defgrad[2][2];

  tk::real a = 0.0;
  
  // hydrodynamic contribution
  auto p_eff = std::max( 1.0e-15, apr+(alpha*m_pstiff) );
  a += m_gamma * p_eff / arho;

  // Other pressure contributions
  a += 

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
GodunovRomenskiSolid::totalenergy(
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
//! \return Material specific total energy using the GodunovRomenskiSolid EoS
// *************************************************************************
{
  // obtain hydro contribution to energy
  tk::real rhoEh = (pr + m_pstiff) / (m_gamma-1.0) + 0.5 * rho *
    (u*u + v*v + w*w) + m_pstiff;
  // obtain elastic contribution to energy
  tk::real eps2;
  tk::real rhoEe = elasticEnergy(defgrad, eps2);

  return (rhoEh + rhoEe);
}

tk::real
GodunovRomenskiSolid::temperature(
  tk::real arho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real arhoE,
  tk::real alpha,
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad ) const
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
//! \return Material temperature using the GodunovRomenskiSolid EoS
// *************************************************************************
{
  // obtain elastic contribution to energy
  tk::real eps2;
  auto arhoEe = alpha*elasticEnergy(defgrad, eps2);
  // obtain hydro contribution to energy
  auto arhoEh = arhoE - arhoEe;

  tk::real t = (arhoEh - 0.5 * arho * (u*u + v*v + w*w) - alpha*m_pstiff)
               / (arho*m_cv);
  return t;
}

tk::real
GodunovRomenskiSolid::min_eff_pressure(
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
GodunovRomenskiSolid::elasticEnergy(
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
  tk::real& eps2 ) const
// *************************************************************************
//! \brief Calculate elastic contribution to material energy from the material
//!   density, and deformation gradient tensor
//! \param[in] defgrad Material inverse deformation gradient tensor
//! \param[in/out] eps2 Elastic shear distortion
//! \return Material elastic energy using the GodunovRomenskiSolid EoS
//! \details This function returns the material elastic energy, and also stores
//!   the elastic shear distortion for further use
// *************************************************************************
{
  // Compute deviatoric part of Hencky tensor
  auto devH = tk::getDevHencky(defgrad);

  // Compute elastic shear distortion
  eps2 = 0.0;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      eps2 += devH[i][j]*devH[i][j];

  // compute elastic energy
  auto rhoEe = m_mu * eps2;

  return rhoEe;
}
