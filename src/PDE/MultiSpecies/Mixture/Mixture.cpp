// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Mixture.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Multispecies mixture function
  \details   This file declares functions for computing mixture flow quantities
*/
// *****************************************************************************

#include "MultiSpecies/Mixture/Mixture.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"

using inciter::Mixture;

Mixture::Mixture(
  std::size_t nspec,
  const std::vector< tk::real >& ugp,
  const std::vector< EOS >& mat_blk) :
  m_nspec(nspec),
  m_mix_density(0),
  m_mix_R(0),
  m_Ys(nspec, 0)
// *************************************************************************
//  Constructor (use during timestepping)
//! \brief Initialize a mixture class using the state vector
//! \param[in] nspec Number of species in mixture
//! \param[in] ugp State vector
//! \param[in] mat_blk EOS material block
// *************************************************************************
{
  // Compute total density
  m_mix_density = 0.;
  for (std::size_t k=0; k<m_nspec; ++k)
    m_mix_density += ugp[multispecies::densityIdx(m_nspec, k)];

  // Compute mass fractions
  for (std::size_t k=0; k<m_nspec; ++k)
    m_Ys[k] = ugp[multispecies::densityIdx(m_nspec, k)] / m_mix_density;

  // Compute mixture gas constant
  m_mix_R = 0.;
  for (std::size_t k = 0; k < m_nspec; k++)
    m_mix_R += m_Ys[k] * mat_blk[k].compute< EOS::gas_constant >();
}

Mixture::Mixture(
  std::size_t nspec,
  const std::vector< tk::real >& Ys,
  tk::real mix_pressure,
  tk::real temperature,
  const std::vector< EOS >& mat_blk) :
  m_nspec(nspec),
  m_mix_density(0),
  m_mix_R(0),
  m_Ys(Ys)
// *************************************************************************
//  Constructor (use during initialization)
//! \brief Initialize a mixture class using the mixture thermodynamics and
//!   known mass fractions.
//! \param[in] nspec Number of species in mixture
//! \param[in] Ys Mass fractions
//! \param[in] mix_pressure Mixture pressure
//! \param[in] temperature Temperature
//! \param[in] mat_blk EOS material block
// *************************************************************************
{
  // Compute mixture gas constant
  m_mix_R = 0.;
  for (std::size_t k = 0; k < m_nspec; k++)
    m_mix_R += m_Ys[k] * mat_blk[k].compute< EOS::gas_constant >();

  // Compute total density (via ideal gas EOS)
  m_mix_density = mix_pressure / (m_mix_R * temperature);
}

tk::real
Mixture::frozen_soundspeed(
  tk::real mix_density,
  tk::real mix_temp,
  const std::vector< EOS >& mat_blk) const
// *************************************************************************
//! \brief Calculate frozen speed of sound based on the mixture composition
//!   and species parameters.
//! \param[in] mix_density Mixture density (sum of species density)
//! \param[in] mix_temp Mixture temperature (provided at call-site, since
//!   it is reconstructed separately)
//! \param[in] mat_blk EOS material block
//! \return Mixture speed of sound using the ideal gas EoS
// *************************************************************************
{
  // Clip pressure
  auto mix_peff = std::max( 1.0e-15, pressure(mix_density, mix_temp) );

  // Compute beta, mixture parameters for sound speed calc.
  tk::real mix_Cv = 0.;
  for (std::size_t k = 0; k < m_nspec; k++) {
    mix_Cv += mat_blk[k].compute< EOS::cv >(mix_temp) * m_Ys[k];
  }
  tk::real beta = m_mix_R / mix_Cv;

  // Compute speed of sound
  tk::real a_sq = (1. + beta) * mix_peff / mix_density;
  return std::sqrt(a_sq);
}

tk::real
Mixture::totalenergy(
  tk::real mix_density,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real mix_temp,
  const std::vector< EOS >& mat_blk) const
// *************************************************************************
//! \brief Calculate total energy based on the mixture composition
//!   and species parameters.
//! \param[in] mix_density Mixture density (sum of species density)
//! \param[in] u Velocity component
//! \param[in] v Velocity component
//! \param[in] w Velocity component
//! \param[in] mix_temp Mixture temperature (provided at call-site, since
//!   it is reconstructed separately
//! \param[in] mat_blk EOS material block
//! \return Total energy
// *************************************************************************
{
  // Compute mixture internal energy
  tk::real mix_e = 0.;
  for (std::size_t k = 0; k < m_nspec; k++) {
    mix_e += m_Ys[k] * mat_blk[k].compute< EOS::internalenergy >(mix_temp);
  }

  // Compute total energy
  return mix_density * (mix_e + 0.5 * (u*u + v*v + w*w));
}

tk::real
Mixture::pressure(
  tk::real mix_density,
  tk::real mix_temp ) const
// *************************************************************************
//! \brief Calculate mixture pressure based on the mixture composition
//!   and species parameters.
//! \param[in] mix_density Mixture density (sum of species density)
//! \param[in] mix_temp Mixture temperature
//! \return Mixture pressure
// *************************************************************************
{
  // Compute pressure based on the ideal gas EOS
  return mix_density * m_mix_R * mix_temp;
}

tk::real
Mixture::Cp(
  tk::real mix_temp,
  const std::vector< EOS >& mat_blk) const
// *************************************************************************
//! \brief Calculate mixture specific heat at constant pressure based on the
//!   mixture composition and species parameters.
//! \param[in] mix_temp Mixture temperature
//! \param[in] mat_blk EOS material block
//! \return Mass-fraction averaged mixture specific heat at constant pressure
// *************************************************************************
{
  tk::real mix_Cp = 0.;
  for (std::size_t k = 0; k < m_nspec; k++)
    mix_Cp += m_Ys[k] * mat_blk[k].compute< EOS::cp >(mix_temp);

  return mix_Cp;
}

tk::real
Mixture::viscCoeff(
  const std::vector< EOS >& mat_blk) const
// *************************************************************************
//! \brief Calculate mixture dynamic viscosity coefficient based on the mixture
//!   composition and species parameters.
//! \param[in] mat_blk EOS material block
//! \return Mass-fraction averaged mixture dynamic viscosity coefficient
// *************************************************************************
{
  tk::real mix_visc = 0.;
  for (std::size_t k = 0; k < m_nspec; k++)
    mix_visc += m_Ys[k] * mat_blk[k].compute< EOS::viscCoeff >();

  return mix_visc;
}


tk::real
Mixture::compute_f(
  tk::real e,
  const std::vector< EOS >& mat_blk,
  tk::real x ) const
{
  // get temp from x
  tk::real temp = std::exp(x);
  tk::real f = 0.0;
  for (std::size_t k = 0; k < m_nspec; k++) {
    f += m_Ys[k] * mat_blk[k].compute< EOS::internalenergy >(temp);
  }
  f -= e;
  return f;
}

tk::real
Mixture::compute_fp(
  const std::vector< EOS >& mat_blk,
  tk::real x ) const
{
  // get temp from x
  tk::real temp = std::exp(x);
  tk::real fp = 0.0;
  for (std::size_t k = 0; k < m_nspec; k++) {
    fp += m_Ys[k] * mat_blk[k].compute< EOS::cv >(temp);
  }
  return fp*temp;
}

tk::real
Mixture::temperature(
  tk::real mix_density,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real rhoE,
  const std::vector< EOS >& mat_blk,
  int& converged,
  tk::real T_init ) const
// *************************************************************************
//! \brief Calculate temperature based on the mixture composition
//!   and species parameters.
//! \param[in] mix_density Mixture density (sum of species density)
//! \param[in] u Velocity component
//! \param[in] v Velocity component
//! \param[in] w Velocity component
//! \param[in] rhoE Total energy of the mixture
//! \param[in] mat_blk EOS material block
//! \param[in,out] converged Indicator of Newton method convergence
//! \param[in] T_init Initial temperature guess; default is 1500.
//! \return Mixture pressure
// *************************************************************************
{
  // Compute internal energy
  tk::real e = rhoE / mix_density - 0.5 * (u*u + v*v + w*w);

  // Parameters
  tk::real temp = std::max( 10.0, T_init); // Starting guess
  tk::real tol = std::max(1e-6, 1e-10 * std::abs(e)); // Stopping condition
  tk::real err;
  std::size_t maxiter = 100;
  std::size_t i(0);
  std::array< tk::real, 100 > debug_temp, debug_f;

  // Solve for x, with temp = exp(x) using Newton's method.
  tk::real x = std::log(temp);

  // Compute f
  tk::real f = compute_f(e, mat_blk, x);

  // Convergence flag
  converged = 0;
  if (std::abs(f) <= tol) {
    converged = 1;
    return temp;
  } 

  while (i < maxiter && std::abs(f) > tol) {

    // Compute fp
    tk::real fp = compute_fp(mat_blk, x);

    bool accepted = false;
    std::size_t nline = 10;
    tk::real alpha = 1.0e-00;
    tk::real x_test = x;
    tk::real f_test = f;
    for (std::size_t iline=0; iline<nline; ++iline) {
      // Compute the step, bound it
      tk::real step = std::clamp(f/fp, -5.0, 5.0);
      // Compute new solution
      x_test = x - alpha * step;

      f_test = compute_f(e, mat_blk, x_test);

      // Check if error decreases with current alpha
      if (std::abs(f_test) >= std::abs(f))
        alpha *= 0.5;
      else {
        accepted = true;
        break;
      }
    }

    if (!accepted) {
      converged = 0;
      debug_temp[i] = std::exp(x);
      debug_f[i] = f;
      break;
    }

    f = f_test;
    x = x_test;
    debug_temp[i] = std::exp(x);
    debug_f[i] = f;

    // Check stopping conditions
    err = std::abs(f);
    if (err <= tol) {
      converged = 1;
      break;
    }
    i++;
    if ( i == maxiter ) {
      converged = 0;
      break;
    }
  }
  temp = std::exp(x);
  if (!converged) {
    temp = std::exp(x);
    printf("e, rhoE, mix_density = %e, %e, %e\n", e, rhoE, mix_density);
    printf("u, v, w = %e, %e, %e\n", u, v, w);
    for (std::size_t j=0; j<i; ++j)
      printf("T[%lu], f[%lu] = %e, %e\n", j, j, debug_temp[j], debug_f[j]);
    Throw("Mixture Newton's Method for temperature failed to converge after iterations "
          + std::to_string(i)  + " with temperature " + std::to_string(temp) +
          " at final iteration" );
  }
  return temp;
}

// tk::real
// Mixture::temperature(
//   tk::real mix_density,
//   tk::real u,
//   tk::real v,
//   tk::real w,
//   tk::real rhoE,
//   const std::vector< EOS >& mat_blk,
//   int& converged,
//   tk::real T_init ) const
// // *************************************************************************
// //! \brief Calculate temperature based on the mixture composition
// //!   and species parameters.
// //! \param[in] mix_density Mixture density (sum of species density)
// //! \param[in] u Velocity component
// //! \param[in] v Velocity component
// //! \param[in] w Velocity component
// //! \param[in] rhoE Total energy of the mixture
// //! \param[in] mat_blk EOS material block
// //! \param[in,out] converged Indicator of Newton method convergence
// //! \param[in] T_init Initial temperature guess; default is 1500.
// //! \return Mixture pressure
// // *************************************************************************
// {
//   // Compute internal energy
//   tk::real e = rhoE / mix_density - 0.5 * (u*u + v*v + w*w);

//   // Solve for temperature -- Newton's method
//   tk::real temp = std::max( 10.0, T_init); // Starting guess
//   tk::real tol = std::max(1e-8, 1e-8 * e); // Stopping condition
//   tk::real err;
//   std::size_t maxiter = 10;
//   std::size_t i(0);
//   while (i < maxiter) {
//     // Construct f(T) = e(temp) - e
//     tk::real f_T = 0.;
//     for (std::size_t k = 0; k < m_nspec; k++) {
//       f_T += m_Ys[k] * mat_blk[k].compute< EOS::internalenergy >(temp);
//     }
//     f_T -= e;

//     // Construct f'(T) = cv(temp)
//     tk::real fp_T = 0.;
//     for (std::size_t k = 0; k < m_nspec; k++) {
//       fp_T += m_Ys[k] * mat_blk[k].compute< EOS::cv >(temp);
//     }

//     // Calculate next guess
//     temp = temp - f_T / fp_T;

//     // Check stopping conditions
//     err = abs(f_T);
//     if (err <= tol) {
//       converged = 1;
//       break;
//     }
//     i++;
//     if ( i == maxiter ) {
//       converged = 0;
//     }
//   }

//   return temp;
// }
