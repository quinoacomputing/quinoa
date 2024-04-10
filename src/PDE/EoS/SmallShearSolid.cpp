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
#include "Vector.hpp"
#include "EoS/SmallShearSolid.hpp"

// Lapacke forward declarations
extern "C" {

using lapack_int = long;

#define LAPACK_ROW_MAJOR 101

lapack_int LAPACKE_dgeev(int, char, char, lapack_int, double*, lapack_int,
  double*, double*, double*, lapack_int, double*, lapack_int );

}

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
  // deformation gradient
  auto defgrad = adefgrad;
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      defgrad[i][j] /= alpha;
  }

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
SmallShearSolid::CauchyStress(
  tk::real,
  tk::real,
  tk::real,
  tk::real,
  tk::real,
  tk::real alpha,
  std::size_t /*imat*/,
  const std::array< std::array< tk::real, 3 >, 3 >& adefgrad ) const
// *************************************************************************
//! \brief Calculate the elastic Cauchy stress tensor from the material density,
//!   momentum, total energy, and inverse deformation gradient tensor using the
//!   SmallShearSolid equation of state
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for
//!   the single-material system, this argument can be left unspecified by
//!   the calling code
// //! \param[in] imat Material-id who's EoS is required. Default is 0, so that
// //!   for the single-material system, this argument can be left unspecified
// //!   by the calling code
//! \param[in] adefgrad Material inverse deformation gradient tensor
//!   (alpha_k * g_k).
//! \return Material Cauchy stress tensor (alpha_k * sigma_k) calculated using
//!   the SmallShearSolid EoS
// *************************************************************************
{
  std::array< std::array< tk::real, 3 >, 3 > asig{{{0,0,0}, {0,0,0}, {0,0,0}}};

  // deformation gradient
  auto defgrad = adefgrad;
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      defgrad[i][j] /= alpha;
  }

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
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      asig[i][j] += m_mu*alpha*devbt[i][j];
  }

  return asig;
}

tk::real
SmallShearSolid::soundspeed(
  tk::real arho,
  tk::real apr,
  tk::real alpha,
  std::size_t imat,
  const std::array< std::array< tk::real, 3 >, 3 >& adefgrad,
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
//! \param[in] adefgrad Material inverse deformation gradient tensor
//!   (alpha_k * g_k) with the first dimension aligned to direction in which
//!   wave speeds are required. Default is 0, so that for the single-material
//!   system, this argument can be left unspecified by the calling code
//  //! \param[in] adefgradn Material inverse deformation gradient tensor in
//  //!   direction of vector n (alpha_k * g_ij * n_j). Default is 0, so that for
//  //!   the single-material system, this argument can be left unspecified by the
//  //!   calling code
//  //! \param[in] asigman Material traction vector in normal direction
//  //!   (alpha * sigma_ij * n_j ). Default is 0, so that for the single-material
//  //!   system, this argument can be left unspecified by the calling code
//! \return Material speed of sound using the SmallShearSolid EoS
// *************************************************************************
{
  // deformation gradient
  auto g = adefgrad;
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      g[i][j] /= alpha;
  }

  tk::real g__11 = g[0][0];
  tk::real g__12 = g[0][1];
  tk::real g__13 = g[0][2];
  tk::real g__21 = g[1][0];
  tk::real g__22 = g[1][1];
  tk::real g__23 = g[1][2];
  tk::real g__31 = g[2][0];
  tk::real g__32 = g[2][1];
  tk::real g__33 = g[2][2];

  std::array< std::array< tk::real, 3 >, 3> dsigdg;

  // d(sigma_11)/d(g_11)
  dsigdg[0][0] = ( (((g__13 * g__13 * g__22 - 3 * g__12 * g__13 * g__23
    - 2 * g__22 * (g__12 * g__12 + g__22 * g__22 + g__23 * g__23)) * g__33
    - (-2 * g__13 * g__13 * g__23 - 3 * g__12 * g__13 * g__22
    + g__23 * (g__12 * g__12 - 2 * g__22 * g__22 - 2 * g__23 * g__23)) * g__32)
    * g__31 * g__31) + (( ((3 * g__12 * g__23 - 2 * g__13 * g__22) * g__33
    * g__33) +  (g__32 * (g__12 * g__22 - g__13 * g__23) * g__33) -  (3 * g__22
    * (g__22 * g__22 + g__23 * g__23 + g__32 * g__32) * g__13) + 0.3e1 *  g__12
    *  g__23 * ( (g__22 * g__22) +  (g__23 * g__23) + 0.2e1 / 0.3e1 *  g__32
    *  g__32)) * g__11 + 0.3e1 * (( (g__12 * g__13) + 0.4e1 / 0.3e1 *  g__22
    *  g__23) *  (g__33 * g__33) +  g__32 * ( (g__12 * g__12) -  (g__13 * g__13)
    + 0.4e1 / 0.3e1 *  g__22 *  g__22 - 0.4e1 / 0.3e1 *  g__23 *  g__23) * g__33
    + (g__13 * g__13 * g__22 * g__23) +  (g__12 * (g__22 * g__22 - g__23 * g__23
    - g__32 * g__32) * g__13) -  g__22 *  g__23 * ( (g__12 * g__12)
    + 0.4e1 / 0.3e1 *  g__32 *  g__32)) * g__21) *  g__31 +  (g__33 * g__22
    - g__23 * g__32) *  (g__22 * g__22 + g__23 * g__23 + g__32 * g__32 + g__33
    * g__33) * g__11 * g__11 - 0.3e1 * ( ( std::pow( g__33,  3) * g__12)
    -  (g__33 * g__33 * g__13 * g__32) + (- (g__13 * g__22 * g__23) / 0.3e1
    + 0.2e1 / 0.3e1 *  g__12 * ( (g__22 * g__22) + 0.3e1 / 0.2e1 *  g__23
    * g__23 + 0.3e1 / 0.2e1 *  g__32 *  g__32)) *  g__33 +  (((-3 * g__22
    * g__22 - 2 * g__23 * g__23 - 3 * g__32 * g__32) * g__13 + g__12 * g__22
    * g__23) * g__32) / 0.3e1) * g__21 * g__11 + (- (14 * g__12 * g__12 * g__22)
    - 0.2e1 * g__21 * g__21 *  g__22 -  (14 *  std::pow( g__22,  3)))
    * std::pow( g__33,  3) + 0.14e2 * ( (2 * g__12 * g__13 * g__22)
    + g__23 * ( (g__12 * g__12) + g__21 * g__21 / 0.7e1
    + (3 * g__22 * g__22))) *  g__32 *  (g__33 * g__33)
    + (-0.2e1 *  g__22 * (g__21 * g__21 +  (7 * g__22 * g__22)
    +  (7 * g__32 * g__32)) *  (g__13 * g__13) + 0.3e1 *  g__12
    * g__23 * (g__21 * g__21 + 0.28e2 / 0.3e1 *  g__22 *  g__22
    - 0.28e2 / 0.3e1 *  g__32 *  g__32) *  g__13 +  g__22
    * ((g__21 * g__21 -  (14 * g__23 * g__23)) *  (g__12 * g__12)
    - 0.2e1 *  (g__32 * g__32) * (g__21 * g__21 +  (21 * g__23 * g__23))))
    *  g__33 + 0.2e1 *  g__32 * (- g__23 * (g__21 * g__21 -  (14 * g__22
    * g__22) -  (14 * g__32 * g__32)) *  (g__13 * g__13) / 0.2e1
    - 0.3e1 / 0.2e1 *  g__22 *  g__12 * (g__21 * g__21 + 0.28e2 / 0.3e1
    *  g__23 *  g__23) *  g__13 +  g__23 * (g__21 * g__21
    + (7 * g__23 * g__23)) *  (g__12 * g__12 + g__32 * g__32)))
    * m_mu * std::pow(std::pow( g__33 * (g__11 *  g__22 -  g__12 * g__21)
    +  (g__12 * g__23 * g__31) +  g__13 * (g__21 *  g__32 -  (g__22 * g__31))
    - g__11 *  g__23 *  g__32, 0.2e1), -0.2e1 / 0.3e1) / ( ((g__12 * g__23
    - g__13 * g__22) * g__31) +  (g__33 * g__22 - g__23 * g__32) * g__11
    - g__21 *  (g__33 * g__12 - g__13 * g__32)) / 0.9e1;

  // d(sigma_11)/d(g_21)
  dsigdg[0][1] = 0.2e1 / 0.9e1 * m_mu * (((-g__33 * g__12/0.2e1 - g__13 * g__32)
    * g__31 * g__31 + 0.3e1 / 0.2e1 * g__11 * (-g__12 * g__13 + g__32 * g__33)
    * g__31 + g__12 * (g__11 * g__11 + 0.7e1 * g__12 * g__12 + 0.7e1 * g__32
    * g__32) * g__33 + g__13 * g__32 * (g__11 * g__11 - 0.14e2 * g__12 * g__12
    - 0.14e2 * g__32 * g__32) / 0.2e1) * g__23 * g__23 + (0.3e1 / 0.2e1 * g__22
    * (-g__12 * g__32 + g__13 * g__33) * g__31 * g__31 + ((g__33 * g__33 * g__12
    + g__33 * g__13 * g__32 / 0.2e1 + 0.3e1 / 0.2e1 * g__12 * (g__12 * g__12
    + g__13 * g__13 + g__32 * g__32)) * g__21 - 0.3e1 / 0.2e1 * g__11 * g__22
    * (g__12 * g__12 - g__13 * g__13 - g__32 * g__32 + g__33 * g__33)) * g__31
    - 0.3e1 / 0.2e1 * g__11 * (g__33 * g__33 * g__32 + g__12 * g__13
    * g__33 / 0.3e1 + g__32 * (g__12 * g__12 + 0.2e1 / 0.3e1 * g__13 * g__13
    + g__32 * g__32)) * g__21 - 0.14e2 * (g__33 * g__33 * g__12 * g__32
    + (0.3e1 / 0.28e2 * g__11 * g__11 * g__13 + g__12 * g__12 * g__13 - g__13
    * g__32 * g__32) * g__33 - 0.3e1 / 0.28e2 * (g__11 * g__11 + 0.28e2 / 0.3e1
    * g__13 * g__13) * g__12 * g__32) * g__22) * g__23
    + (g__12 * (g__12 * g__12 + g__13 * g__13 + g__22 * g__22) * g__33
    - g__13 * (g__12 * g__12 + g__13 * g__13 - g__22 * g__22 / 0.2e1) * g__32)
    * g__31 * g__31 + (-0.3e1 / 0.2e1 * g__22 * (g__33 * g__33 * g__13 + g__33
    * g__12 * g__32 / 0.3e1 + g__13 * (g__12 * g__12 + g__13 * g__13
    + 0.2e1 / 0.3e1 * g__32 * g__32)) * g__21 - 0.2e1 * g__11 * (g__33 * g__33
    * g__12 * g__13 + g__32 * (g__12 * g__12 - g__13 * g__13 + 0.3e1 / 0.4e1
    * g__22 * g__22) * g__33 - 0.3e1 / 0.4e1 * g__13 * g__12 * (g__22 * g__22
    + 0.4e1 / 0.3e1 * g__32 * g__32))) * g__31 - (g__33 * g__12 - g__13 * g__32)
    * (g__12 * g__12 + g__13 * g__13 + g__32 * g__32 + g__33 * g__33) * g__21
    * g__21 / 0.2e1 + 0.3e1 / 0.2e1 * g__22 * g__11 * (std::pow(g__33, 0.3e1)
    + (0.2e1 / 0.3e1 * g__12 * g__12 + g__13 * g__13 + g__32 * g__32) * g__33
    + g__12 * g__13 * g__32 / 0.3e1) * g__21 + g__12 * (g__11 * g__11 + 0.7e1
    * g__12 * g__12 + 0.7e1 * g__22 * g__22) * std::pow(g__33, 0.3e1) - g__13
    * g__32 * (g__11 * g__11 + 0.21e2 * g__12 * g__12 + 0.7e1 * g__22 * g__22)
    * g__33 * g__33 - ((g__11 * g__11 - 0.14e2 * g__13 * g__13) * g__22 * g__22
    - 0.2e1 * g__32 * g__32 * (g__11 * g__11 + 0.21e2 * g__13 * g__13)) * g__12
    * g__33 / 0.2e1 - g__13 * g__32 * (g__22 * g__22 + g__32 * g__32) * (g__11
    * g__11 + 0.7e1 * g__13 * g__13)) * std::pow(std::pow(g__33 * (g__11 * g__22
    - g__12 * g__21) + g__12 * g__23 * g__31 + g__13 * (g__21 * g__32 - g__22
    * g__31) - g__11 * g__23 * g__32, 0.2e1), -0.2e1 / 0.3e1) / ((-g__11
    * g__32 + g__12 * g__31) * g__23 - g__13 * g__22 * g__31
    + (-g__33 * g__12 + g__13 * g__32) * g__21 + g__22 * g__33 * g__11);

  // d(sigma_11)/d(g_31)
  dsigdg[0][2] = -0.2e1 / 0.9e1 * m_mu * ((-g__12 * g__31 * g__31 / 0.2e1
    + 0.3e1 / 0.2e1 * g__11 * g__31 * g__32 + g__12 * (g__11 * g__11
    + 0.7e1 * g__12 * g__12 + 0.7e1 * g__32 * g__32)) * std::pow(g__23, 0.3e1)
    + (((-0.3e1 / 0.2e1 * g__32 * g__11 + g__12 * g__31) * g__33
    - 0.2e1 * (g__11 * g__12 + 0.3e1 / 0.4e1 * g__31 * g__32) * g__13) * g__21
    - 0.3e1 / 0.2e1 * ((g__11 * g__31 + 0.28e2 / 0.3e1 * g__12 * g__32) * g__33
    + 0.2e1 / 0.3e1 * g__13 * (g__11 * g__11 + 0.21e2 * g__12 * g__12
    - g__31 * g__31 / 0.2e1 + 0.7e1 * g__32 * g__32)) * g__22) * g__23 * g__23
    + ((-g__33 * g__33 * g__12 / 0.2e1 + 0.3e1 / 0.2e1 * g__33 * g__13 * g__32
    + g__12 * (g__12 * g__12 + g__13 * g__13 + g__32 * g__32)) * g__21 * g__21
    + 0.3e1 / 0.2e1 * (g__33 * g__33 * g__11 + g__33 * g__13 * g__31 / 0.3e1
    - g__12 * g__31 * g__32 / 0.3e1 - 0.4e1 / 0.3e1 * g__11 * (g__12 * g__12
    - g__13 * g__13 + 0.3e1 / 0.4e1 * g__32 * g__32)) * g__22 * g__21
    + g__12 * (g__11 * g__11 + 0.7e1 * g__12 * g__12 + 0.7e1 * g__22 * g__22)
    * g__33 * g__33 - 0.3e1 / 0.2e1 * g__13 * (g__11 * g__12 * g__31 / 0.3e1
    + g__32 * (g__11 * g__11 + 0.28e2 / 0.3e1 * g__12 * g__12 - 0.28e2 / 0.3e1
    * g__22 * g__22)) * g__33 - g__12 * (g__12 * g__12 + g__13 * g__13
    + g__22 * g__22) * g__31 * g__31 / 0.2e1 + g__32 * g__11 * (g__12 * g__12
    + 0.3e1 / 0.2e1 * g__13 * g__13 + 0.3e1 / 0.2e1 * g__22 * g__22) * g__31
    + g__12 * ((g__22 * g__22 - g__32 * g__32 / 0.2e1) * g__11 * g__11
    + (0.21e2 * g__22 * g__22 + 0.7e1 * g__32 * g__32) * g__13 * g__13))
    * g__23 - g__22 * (g__33 * g__33 * g__13 + 0.3e1 / 0.2e1 * g__33
    * g__12 * g__32 + g__13 * (g__12 * g__12 + g__13 * g__13 - g__32
    * g__32 / 0.2e1)) * g__21 * g__21 + (-0.3e1 / 0.2e1 * g__33 * g__33
    * g__11 * g__12 * g__13 + (0.3e1 / 0.2e1 * g__12 * (g__12 * g__12
    + g__13 * g__13 + g__22 * g__22) * g__31 - 0.3e1 / 0.2e1 * g__32
    * g__11 * (g__12 * g__12 - g__13 * g__13 - g__22 * g__22)) * g__33
    + 0.2e1 * g__13 * (-0.3e1 / 0.4e1 * (g__12 * g__12 + g__13 * g__13
    + 0.2e1 / 0.3e1 * g__22 * g__22) * g__32 * g__31 + g__11 * g__12
    * (g__22 * g__22 + 0.3e1 / 0.4e1 * g__32 * g__32))) * g__21 + g__22
    * (g__13 * (g__11 * g__11 - 0.14e2 * g__12 * g__12 - 0.14e2 * g__22 * g__22)
    * g__33 * g__33 + (-0.3e1 * g__11 * (g__12 * g__12 + 0.2e1 / 0.3e1 * g__13
    * g__13 + g__22 * g__22) * g__31 + 0.3e1 * (g__11 * g__11 + 0.28e2 / 0.3e1
    * g__13 * g__13) * g__12 * g__32) * g__33 - 0.2e1 * g__13 * ((-g__12
    * g__12 / 0.2e1 - g__13 * g__13 / 0.2e1 - g__22 * g__22 / 0.2e1)
    * g__31 * g__31 - g__11 * g__12 * g__31 * g__32 / 0.2e1 + (g__22 * g__22
    + g__32 * g__32) * (g__11 * g__11 + 0.7e1 * g__13 * g__13))) / 0.2e1)
    * std::pow(std::pow(g__33 * (g__11 * g__22 - g__12 * g__21) + g__12
    * g__23 * g__31 + g__13 * (g__21 * g__32 - g__22 * g__31) - g__11
    * g__23 * g__32, 0.2e1), -0.2e1 / 0.3e1) / ((-g__32 * g__11 + g__12 * g__31)
    * g__23 + (-g__33 * g__12 + g__13 * g__32) * g__21 + (g__11 * g__33 - g__13
    * g__31) * g__22);

  // d(sigma_21)/d(g_11)
  dsigdg[1][0] = ( ((g__33 * (3 * g__13 * g__13 + 4 * g__23 * g__23) * g__21
    + ((-4 * g__13 * g__13 - 4 * g__23 * g__23) * g__31 + g__33 * g__11
    * g__13) * g__23) * g__32 * g__32) +  ((((g__13 * g__23 * g__23
    - 6 * g__33 * g__33 * g__13) * g__21 - g__23 * (-7 * g__33 * g__13
    * g__31 + g__11 * (g__23 * g__23 + g__33 * g__33))) * g__12 - ((g__13
    * g__13 * g__23 + 8 * g__23 * g__33 * g__33) * g__21 - g__33 * (g__13
    * g__13 + 8 * g__23 * g__23) * g__31 + g__13 * g__11 * (g__33 - g__23)
    * (g__33 + g__23)) * g__22) * g__32) +  (3 * (g__23 * g__23 + g__33 * g__33)
    * (g__21 * g__33 - g__23 * g__31) * g__12 * g__12) +  (g__22 * (-7 * g__33
    * g__13 * g__23 * g__21 - g__13 * (-6 * g__23 * g__23 + g__33 * g__33)
    * g__31 + g__33 * g__11 * (g__23 * g__23 + g__33 * g__33)) * g__12)
    + 0.4e1 *  (g__22 * g__22) * ( ((g__13 * g__13 * g__33
    + std::pow( g__33,  3)) * g__21) - (( (g__33 * g__33) + 0.3e1 / 0.4e1
    * g__13 *  g__13) *  g__31 +  (g__33 * g__11 * g__13) / 0.4e1)
    * g__23)) * m_mu * std::pow( std::pow( (g__33 * (g__11 * g__22
    - g__12 * g__21) + g__12 * g__23 * g__31 + g__13 * (g__21 * g__32
    - g__22 * g__31) - g__11 * g__23 * g__32),  2), -0.2e1 / 0.3e1)
    / ((-g__11 * g__23 + g__13 * g__21) * g__32 + (-g__21 * g__33
    + g__23 * g__31) * g__12 + (g__11 * g__33 - g__13 * g__31)
    * g__22) / 0.3e1;

  // d(sigma_21)/d(g_21)
  dsigdg[1][1] = -0.4e1 / 0.3e1 * ((g__33 * (g__13 * g__13
    + 0.3e1 / 0.4e1 * g__23 * g__23) * g__11
    + ((-0.4e1 * g__13 * g__13 - 0.4e1 * g__23 * g__23)
    * g__31 + g__33 * g__21 * g__23) * g__13 / 0.4e1) * g__32
    * g__32 + ((-0.3e1 / 0.2e1 * (g__33 * g__33 - g__13 * g__13 / 0.6e1)
    * g__23 * g__11 - g__13 * (-0.7e1 * g__33 * g__23 * g__31 + g__21
    * (g__13 * g__13 + g__33 * g__33)) / 0.4e1) * g__22 - 0.2e1 * (g__13
    * (g__33 * g__33 + g__23 * g__23 / 0.8e1) * g__11 - g__33 * (g__13
    * g__13 + g__23 * g__23 / 0.8e1) * g__31 + g__21 * g__23 * (g__33 - g__13)
    * (g__33 + g__13) / 0.8e1) * g__12) * g__32 + 0.3e1 / 0.4e1
    * (g__13 * g__13 + g__33 * g__33) * (g__11 * g__33 - g__13 * g__31)
    * g__22 * g__22 + g__12 * (-0.7e1 * g__33 * g__11 * g__13 * g__23 - g__23
    * (-0.6e1 * g__13 * g__13 + g__33 * g__33) * g__31 + g__33 * g__21
    * (g__13 * g__13 + g__33 * g__33)) * g__22 / 0.4e1 + ((g__23 * g__23
    * g__33 + std::pow(g__33, 0.3e1)) * g__11 - g__13 * ((g__33 * g__33
    + 0.3e1 / 0.4e1 * g__23 * g__23) * g__31 + g__33 * g__21 * g__23 / 0.4e1))
    * g__12 * g__12) * m_mu * std::pow(std::pow(g__33 * (g__11 * g__22 - g__12
    * g__21) + g__12 * g__23 * g__31 + g__13 * (g__21 * g__32 - g__22 * g__31)
    - g__11 * g__23 * g__32, 0.2e1), -0.2e1 / 0.3e1) / ((-g__11 * g__23 + g__13
    * g__21) * g__32 + (g__11 * g__33 - g__13 * g__31) * g__22 - g__12
    * (g__21 * g__33 - g__23 * g__31));

  // d(sigma_21)/d(g_31)
  dsigdg[1][2] = 0.4e1 / 0.3e1 * m_mu * (((-std::pow(g__13, 0.3e1)
    - g__33 * g__33 * g__13) * g__21 + 0.3e1 / 0.4e1 * ((g__33 * g__33
    + 0.4e1 / 0.3e1 * g__13 * g__13) * g__11 + g__33 * g__13 * g__31 / 0.3e1)
    * g__23) * g__22 * g__22 + ((0.7e1 / 0.4e1 * g__33 * g__13 * g__23 * g__21
    + g__33 * (g__13 * g__13 - 0.6e1 * g__23 * g__23) * g__11 / 0.4e1
    - g__13 * g__31 * (g__13 * g__13 + g__23 * g__23) / 0.4e1) * g__32
    - (-g__23 * (0.8e1 * g__13 * g__13 + g__33 * g__33) * g__21
    + g__13 * (0.8e1 * g__23 * g__23 + g__33 * g__33) * g__11
    - g__33 * g__31 * (g__13 - g__23) * (g__13 + g__23)) * g__12 / 0.4e1)
    * g__22 + 0.3e1 / 0.4e1 * (g__13 * g__13 + g__23 * g__23)
    * (g__11 * g__23 - g__13 * g__21) * g__32 * g__32 - 0.7e1 / 0.4e1
    * (-0.6e1 / 0.7e1 * g__33 * (g__13 * g__13 - g__23 * g__23 / 0.6e1)
    * g__21 + (g__33 * g__11 * g__13 - g__31 * (g__13 * g__13 + g__23
    * g__23) / 0.7e1) * g__23) * g__12 * g__32 + g__12 * g__12
    * ((-0.3e1 / 0.4e1 * g__33 * g__33 * g__13 - g__13 * g__23 * g__23)
    * g__21 + g__23 * (g__11 * (g__23 * g__23 + g__33 * g__33)
    - g__33 * g__13 * g__31 / 0.4e1))) * std::pow(std::pow(g__33
    * (g__11 * g__22 - g__12 * g__21) + g__12 * g__23 * g__31 + g__13
    * (g__21 * g__32 - g__22 * g__31) - g__11 * g__23 * g__32, 0.2e1),
    -0.2e1 / 0.3e1) / ((-g__11 * g__23 + g__13 * g__21) * g__32
    + (g__11 * g__33 - g__13 * g__31) * g__22 - g__12 * (g__21
    * g__33 - g__23 * g__31));

  // d(sigma_31)/d(g_11)
  dsigdg[2][0] = - ((g__32 * (3 * g__12 * g__12 + 4 * g__22 * g__22)
    * g__21 + g__22 * ((-4 * g__12 * g__12 - 4 * g__22 * g__22) * g__31
    + g__11 * g__12 * g__32)) * g__33 * g__33 + ((g__12 * (g__22 * g__22
    - 6 * g__32 * g__32) * g__21 - g__22 * (-7 * g__12 * g__31 * g__32
    + g__11 * (g__22 * g__22 + g__32 * g__32))) * g__13
    + (-g__22 * (g__12 * g__12 + 8 * g__32 * g__32) * g__21 + g__32
    * (g__12 * g__12 + 8 * g__22 * g__22) * g__31 + g__11 * g__12
    * (g__22 - g__32) * (g__22 + g__32)) * g__23) * g__33
    + 3 * (g__22 * g__22 + g__32 * g__32) * (g__21 * g__32 - g__22 * g__31)
    * g__13 * g__13 + g__23 * (-7 * g__12 * g__22 * g__32 * g__21
    + (6 * g__12 * g__22 * g__22 - g__12 * g__32 * g__32) * g__31
    + g__32 * g__11 * (g__22 * g__22 + g__32 * g__32)) * g__13 - g__23
    * g__23 * ((-4 * g__12 * g__12 * g__32 - 4 *  std::pow( g__32,  3))
    * g__21 + g__22 * (g__31 * (3 * g__12 * g__12 + 4 * g__32 * g__32)
    + g__11 * g__12 * g__32))) * m_mu * std::pow(  std::pow( (g__33
    * (g__11 * g__22 - g__12 * g__21) + g__12 * g__23 * g__31 + g__13
    * (g__21 * g__32 - g__22 * g__31) - g__11 * g__23 * g__32),  2),
    -0.2e1 / 0.3e1) / (g__33 * (g__11 * g__22 - g__12 * g__21) + g__13
    * (g__21 * g__32 - g__22 * g__31) - g__23 * (g__32 * g__11 - g__12
    * g__31)) / 0.3e1;

  // d(sigma_31)/d(g_21)
  dsigdg[2][1] = 0.4e1 / 0.3e1 * (((-std::pow(g__12, 0.3e1)
    - g__12 * g__22 * g__22) * g__31 + ((g__12 * g__12 + 0.3e1 / 0.4e1
    * g__22 * g__22) * g__11 + g__12 * g__21 * g__22 / 0.4e1) * g__32)
    * g__33 * g__33 + ((0.7e1 / 0.4e1 * g__12 * g__22 * g__31 * g__32
    + g__22 * (g__12 * g__12 - 0.6e1 * g__32 * g__32) * g__11 / 0.4e1
    - g__12 * g__21 * (g__12 * g__12 + g__32 * g__32) / 0.4e1) * g__23
    - g__13 * ((-0.8e1 * g__12 * g__12 - g__22 * g__22) * g__32 * g__31
    + g__12 * (g__22 * g__22 + 0.8e1 * g__32 * g__32) * g__11 - g__21
    * g__22 * (g__12 - g__32) * (g__12 + g__32)) / 0.4e1) * g__33
    + 0.3e1 / 0.4e1 * (g__12 * g__12 + g__32 * g__32) * (g__32 * g__11
    - g__12 * g__31) * g__23 * g__23 - 0.7e1 / 0.4e1 * g__13
    * (-0.6e1 / 0.7e1 * g__22 * (g__12 * g__12 - g__32 * g__32 / 0.6e1)
    * g__31 + g__32 * (g__11 * g__12 * g__22 - g__21 * (g__12 * g__12
    + g__32 * g__32) / 0.7e1)) * g__23 + g__13 * g__13 * ((-0.3e1 / 0.4e1
    * g__12 * g__22 * g__22 - g__12 * g__32 * g__32) * g__31
    + (g__11 * (g__22 * g__22 + g__32 * g__32) - g__12 * g__21
    * g__22 / 0.4e1) * g__32)) * m_mu * std::pow(std::pow(g__33
    * (g__11 * g__22 - g__12 * g__21) + g__12 * g__23 * g__31 + g__13
    * (g__21 * g__32 - g__22 * g__31) - g__11 * g__23 * g__32, 0.2e1),
    -0.2e1 / 0.3e1) / (g__13 * (g__21 * g__32 - g__22 * g__31)
    + (-g__32 * g__11 + g__12 * g__31) * g__23 + g__33
    * (g__11 * g__22 - g__12 * g__21));

  // d(sigma_31)/d(g_31)
  dsigdg[2][2] = -(((-0.4e1 / 0.3e1 * g__12 * g__32 * g__32
    - 0.4e1 / 0.3e1 * std::pow(g__12, 0.3e1)) * g__21 + 0.4e1 / 0.3e1
    * g__22 * ((g__12 * g__12 + 0.3e1 / 0.4e1 * g__32 * g__32) * g__11
    + g__12 * g__31 * g__32 / 0.4e1)) * g__23 * g__23 + ((0.7e1 / 0.3e1
    * g__12 * g__22 * g__32 * g__21 + g__32 * (g__12 * g__12 - 0.6e1
    * g__22 * g__22) * g__11 / 0.3e1 - g__12 * g__31 * (g__12 * g__12
    + g__22 * g__22) / 0.3e1) * g__33 - 0.8e1 / 0.3e1 * g__13 * (-g__22
    * (g__12 * g__12 + g__32 * g__32 / 0.8e1) * g__21 + g__12 * (g__22
    * g__22 + g__32 * g__32 / 0.8e1) * g__11 - g__31 * g__32 * (g__12
    - g__22) * (g__12 + g__22) / 0.8e1)) * g__23 + (g__12 * g__12
    + g__22 * g__22) * (g__11 * g__22 - g__12 * g__21) * g__33 * g__33
    - 0.7e1 / 0.3e1 * (-0.6e1 / 0.7e1 * g__32 * (g__12 * g__12 - g__22
    * g__22 / 0.6e1) * g__21 + (g__11 * g__12 * g__32 - g__31 * (g__12
    * g__12 + g__22 * g__22) / 0.7e1) * g__22) * g__13 * g__33
    + 0.4e1 / 0.3e1 * (-g__12 * (g__22 * g__22 + 0.3e1 / 0.4e1
    * g__32 * g__32) * g__21 + g__22 * (g__11 * (g__22 * g__22
    + g__32 * g__32) - g__12 * g__31 * g__32 / 0.4e1)) * g__13 * g__13)
    * m_mu * std::pow(std::pow(g__33 * (g__11 * g__22 - g__12 * g__21)
    + g__12 * g__23 * g__31 + g__13 * (g__21 * g__32 - g__22 * g__31)
    - g__11 * g__23 * g__32, 0.2e1), -0.2e1 / 0.3e1) / (g__13 * (g__21
    * g__32 - g__22 * g__31) + (-g__32 * g__11 + g__12 * g__31) * g__23
    + g__33 * (g__11 * g__22 - g__12 * g__21));

  // Define amat
  double amat[9];
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
    {
      amat[i*3+j] = 0.0;
      for (std::size_t k=0; k<3; ++k)
      {
        amat[i*3+j] -= adefgrad[k][j]*dsigdg[i][k];
      }
      amat[i*3+j] /= arho;
    }

  double eig_real[3], eig_imag[3];
  double vl[3], vr[3];
  #ifndef NDEBUG
  lapack_int ierr =
  #endif
    LAPACKE_dgeev (LAPACK_ROW_MAJOR, 'N', 'N', 3, amat, 3,
      eig_real, eig_imag, vl, 3, vr, 3);
  Assert(ierr==0, "Lapack failed to compute eigenvalues");

  // Manually find max
  tk::real eig_max = eig_real[0];
  for (std::size_t i=1; i<3; i++)
    if (eig_real[i] > eig_max)
      eig_max = eig_real[i];
  tk::real a = std::sqrt(eig_max);

  // hydrodynamic contribution
  auto p_eff = std::max( 1.0e-15, apr+(alpha*m_pstiff) );
  a += std::sqrt( m_gamma * p_eff / arho );

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
//!   g_k. Default is 0, so that for the single-material system,
//!   this argument can be left unspecified by the calling code
//! \return Material specific total energy using the SmallShearSolid EoS
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
SmallShearSolid::temperature(
  tk::real arho,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real arhoE,
  tk::real alpha,
  const std::array< std::array< tk::real, 3 >, 3 >& adefgrad ) const
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
//!   (alpha_k * g_k). Default is 0, so that for the single-material system,
//!   this argument can be left unspecified by the calling code
//! \return Material temperature using the SmallShearSolid EoS
// *************************************************************************
{
  // deformation gradient
  auto defgrad = adefgrad;
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      defgrad[i][j] /= alpha;
  }

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
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
  tk::real& eps2 ) const
// *************************************************************************
//! \brief Calculate elastic contribution to material energy from the material
//!   density, and deformation gradient tensor
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
  auto rhoEe = m_mu * eps2;

  return rhoEe;
}
