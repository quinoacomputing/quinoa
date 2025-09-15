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

void
SmallShearSolid::setRho0( tk::real rho0 )
// *************************************************************************
//  Set rho0 EOS parameter; i.e. the initial density
//! \param[in] rho0 Initial material density that needs to be stored
// *************************************************************************
{
  m_rho0 = rho0;
}

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
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad ) const
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
//! \param[in] defgrad Material inverse deformation gradient tensor
//!   (g_k). Default is 0, so that for the single-material system,
//!   this argument can be left unspecified by the calling code
//! \return Material partial pressure (alpha_k * p_k) calculated using the
//!   SmallShearSolid EoS
// *************************************************************************
{
  if (!(alpha > 0.0)) return 0.0;

  // elastic with robust elasticEnergy()
  tk::real eps2;
  const tk::real arhoEe = alpha * elasticEnergy(defgrad, eps2);

  // hydro part (α-weighted)
  const tk::real arhoEh = arhoE - arhoEe;
  tk::real partpressure = (arhoEh - 0.5*arho*(u*u+v*v+w*w))*(m_gamma-1.0)
                          - alpha*m_gamma*m_pstiff;

  if (!std::isfinite(partpressure)) partpressure = 0.0;      // finite fallback
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
  const std::array< std::array< tk::real, 3 >, 3 >& defgrad ) const
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
//! \param[in] defgrad Material inverse deformation gradient tensor (g_k).
//! \return Material Cauchy stress tensor (alpha_k * sigma_k) calculated using
//!   the SmallShearSolid EoS
// *************************************************************************
{
  std::array<std::array<tk::real,3>,3> asig{{{0,0,0},{0,0,0},{0,0,0}}};
  if (!(alpha > 0.0)) return asig;

  // volumetric shear (spherical) via eps2
  tk::real eps2;
  (void) elasticEnergy(defgrad, eps2);
  const tk::real pmean = - alpha * m_mu * eps2;

  asig[0][0] = asig[1][1] = asig[2][2] = -pmean;

  // deviatoric with det guard
  auto B = tk::getLeftCauchyGreen(defgrad);
  tk::real detB = safe_pos(tk::determinant(B), DET_FLOOR);
  const tk::real Jm23 = std::pow(detB, 1.0/3.0); // = J^{-2/3}
  for (int i=0;i<3;++i) for (int j=0;j<3;++j) B[i][j] /= Jm23;
  const tk::real trB3 = (B[0][0]+B[1][1]+B[2][2]) / 3.0;
  B[0][0]-=trB3; B[1][1]-=trB3; B[2][2]-=trB3;

  for (int i=0;i<3;++i) for (int j=0;j<3;++j)
    asig[i][j] += m_mu * alpha * B[i][j];

  return asig;
}

tk::real
SmallShearSolid::soundspeed(
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
//  //! \param[in] defgrad Material inverse deformation gradient tensor
//  //!   (g_k) with the first dimension aligned to direction in which
//  //!   wave speeds are required. Default is 0, so that for the single-material
//  //!   system, this argument can be left unspecified by the calling code
//! \return Material speed of sound using the SmallShearSolid EoS
// *************************************************************************
{
  /* Rigorous approach: Eigenvalues of full elastic tensor

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
        amat[i*3+j] -= defgrad[k][j]*dsigdg[i][k];
      }
      amat[i*3+j] /= (arho/alpha);
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

  */

  // Approximated elastic contribution, from Barton, P. T. (2019).
  // An interface-capturing Godunov method for the simulation of compressible
  // solid-fluid problems. Journal of Computational Physics, 390, 25-50
  if (!(alpha > 0.0)) return 0.0;

  // per-solid variables with α floor (only for unweighting)
  const tk::real aeff   = alpha_eff(alpha);
  const tk::real rho    = std::max(arho/aeff, RHO_FLOOR);
  const tk::real p_eff  = std::max<tk::real>(apr/aeff + m_pstiff, 1.0e-15);

  // Barton-style approx
  tk::real a2 = (4.0/3.0) * m_mu / rho + m_gamma * p_eff / rho;
  tk::real a  = std::sqrt(std::max<tk::real>(a2, 0.0));
  return std::isfinite(a) ? a : 0.0;
}

tk::real
SmallShearSolid::shearspeed(
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
  if (!(alpha > 0.0)) return 0.0;
  const tk::real rho = std::max(arho/alpha_eff(alpha), RHO_FLOOR);
  tk::real a = std::sqrt(m_mu / rho);
  return std::isfinite(a) ? a : 0.0;
}

tk::real
SmallShearSolid::totalenergy(
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
//! \return Material specific total energy using the SmallShearSolid EoS
// *************************************************************************
{
  if (!(alpha > 0.0)) return 0.0;

  // hydro for stiffened gas
  const tk::real arhoEh = (apr + alpha*m_gamma*m_pstiff)/(m_gamma-1.0)
                          + 0.5*arho*(u*u+v*v+w*w);

  // elastic
  tk::real eps2;
  const tk::real arhoEe = alpha * elasticEnergy(defgrad, eps2);

  const tk::real out = arhoEh + arhoEe;
  return std::isfinite(out) ? out : 0.0;
}

tk::real
SmallShearSolid::temperature(
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
//! \return Material temperature using the SmallShearSolid EoS
// *************************************************************************
{
  if (!(alpha > 0.0)) return 0.0;

  tk::real eps2;
  const tk::real arhoEe = alpha * elasticEnergy(defgrad, eps2);
  const tk::real arhoEh = arhoE - arhoEe;

  const tk::real arho_eff = std::max<tk::real>(arho, A_TRACE*RHO_FLOOR);
  tk::real T = (arhoEh - 0.5*arho*(u*u+v*v+w*w) - alpha*m_pstiff) / (arho_eff*m_cv);
  return std::isfinite(T) ? T : 0.0;
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
  // Precondition g: symmetrize and fix determinant
  auto g = defgrad;
  symmetrize(g);
  tk::real detg = tk::determinant(g);
  if (!(detg > DET_FLOOR)) {
    g = spherical_from_det(g);
  }

  // Isochoric Right Cauchy–Green from the *fixed* g
  auto Ct = tk::getIsochorRightCauchyGreen(g);

  // eps2
  eps2 = 0.5 * (Ct[0][0] + Ct[1][1] + Ct[2][2] - 3.0);

  // elastic energy per-solid (distortional)
  tk::real rhoEe = m_mu * eps2;
  if (!std::isfinite(rhoEe)) { eps2 = 0.0; rhoEe = 0.0; }
  return rhoEe;
}
