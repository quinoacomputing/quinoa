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
  std::size_t nspec) :
  m_nspec(nspec)
// *************************************************************************
//  Constructor (use during timestepping)
//! \brief Initialize a mixture class
//! \param[in] nspec Number of species in mixture
// *************************************************************************
{}

void
Mixture::set_state(std::vector< tk::real > ugp,
  const std::vector< EOS >& mat_blk)
// *************************************************************************
//! \brief Set commonly used mixture properties based off a given state
//!   vector (rho_s, rhou, rhov, rhow, rhoE). Used during timestepping
//! \param[in] ugp State vector
// *************************************************************************
{
  // Compute total density
  m_mix_density = 0.;
  for (std::size_t k=0; k<m_nspec; ++k)
    m_mix_density += ugp[multispecies::densityIdx(m_nspec, k)];

  // Compute mass fractions
  for (std::size_t k=0; k<m_nspec; ++k)
    m_Ys.push_back(ugp[multispecies::densityIdx(m_nspec, k)] /
      m_mix_density);

  // Compute mixture gas constant
  m_mix_R = 0.;
  for (std::size_t k = 0; k < m_nspec; k++)
    m_mix_R += m_Ys[k] * mat_blk[k].compute< EOS::gas_constant >();
}

void
Mixture::set_massfrac(std::vector< tk::real > Ys,
                  tk::real mix_pressure,
                  tk::real temperature,
                  const std::vector< EOS >& mat_blk)
// *************************************************************************
//! \brief Set commonly used mixture properties based off a given set of
//!   mass fractions, a mixture pressure, and a temperature. Used during
//!   initialization, or when setting boundary conditions.
//! \param[in] Ys Mass fractions
//! \param[in] mix_pressure Mixture pressure
//! \param[in] temperature Temperature
// *************************************************************************
{
  // Compute mixture gas constant
  m_mix_R = 0.;
  for (std::size_t k = 0; k < m_nspec; k++)
    m_mix_R += Ys[k] * mat_blk[k].compute< EOS::gas_constant >();

  // Compute total density (via ideal gas EOS)
  m_mix_density = mix_pressure / (m_mix_R * temperature);

  // Store mass fractions
  for (std::size_t k=0; k<m_nspec; ++k)
    m_Ys.push_back(Ys[k]);
}

tk::real
Mixture::frozen_soundspeed(
  tk::real mix_density,
  tk::real mix_pressure,
  const std::vector< EOS >& mat_blk) const
// *************************************************************************
//! \brief Calculate frozen speed of sound based on the mixture composition
//!   and species parameters.
//! \param[in] mix_density Mixture density (sum of species density)
//! \param[in] mix_pressure Mixture pressure (sum of species pressure)
//! \return Mixture speed of sound using the ideal gas EoS
// *************************************************************************
{
  // Clip pressure
  auto mix_peff = std::max( 1.0e-15, mix_pressure );

  // Compute temperature based on equation of state
  tk::real temp = mix_peff / (m_mix_R * mix_density);

  // Compute beta, mixture parameters for sound speed calc.
  tk::real mix_Cv = 0.;
  for (std::size_t k = 0; k < m_nspec; k++) {
    mix_Cv += mat_blk[k].compute< EOS::calc_cv >(temp) * m_Ys[k];
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
  tk::real mix_pressure,
  const std::vector< EOS >& mat_blk) const
// *************************************************************************
//! \brief Calculate total energy based on the mixture composition
//!   and species parameters.
//! \param[in] mix_density Mixture density (sum of species density)
//! \param[in] u Velocity component
//! \param[in] v Velocity component
//! \param[in] w Velocity component
//! \param[in] mix_pressure Mixture pressure
//! \return Total energy
// *************************************************************************
{
  // Calculate temperature from ideal gas EOS
  tk::real temp = mix_pressure / (m_mix_R * mix_density);

  // Compute mixture internal energy
  tk::real mix_e = 0.;
  for (std::size_t k = 0; k < m_nspec; k++) {
    mix_e += m_Ys[k] * mat_blk[k].compute< EOS::calc_e >(temp);
  }

  // Compute total energy
  return mix_density * (mix_e + 0.5 * (u*u + v*v + w*w));
}

tk::real
Mixture::pressure(
  tk::real mix_density,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real rhoE,
  const std::vector< EOS >& mat_blk) const
// *************************************************************************
//! \brief Calculate mixture pressure based on the mixture composition
//!   and species parameters.
//! \param[in] mix_density Mixture density (sum of species density)
//! \param[in] u Velocity component
//! \param[in] v Velocity component
//! \param[in] w Velocity component
//! \param[in] rhoE Total energy of the mixture
//! \return Mixture pressure
// *************************************************************************
{
  // Compute temperature from the state vector
  tk::real temp = temperature(mix_density, u, v, w, rhoE, mat_blk);

  // Compute pressure based on the ideal gas EOS
  return mix_density * m_mix_R * temp;
}

tk::real
Mixture::temperature(
  tk::real mix_density,
  tk::real u,
  tk::real v,
  tk::real w,
  tk::real rhoE,
  const std::vector< EOS >& mat_blk) const
// *************************************************************************
//! \brief Calculate temperature based on the mixture composition
//!   and species parameters.
//! \param[in] mix_density Mixture density (sum of species density)
//! \param[in] u Velocity component
//! \param[in] v Velocity component
//! \param[in] w Velocity component
//! \param[in] rhoE Total energy of the mixture
//! \return Mixture pressure
// *************************************************************************
{
  // Compute internal energy
  tk::real e = rhoE / mix_density - 0.5 * (u*u + v*v + w*w);
  if (e < 1e-8) e = 1e-8; // TODO: standin until positivity is implemented

  // Solve for temperature -- Newton's method
  tk::real temp = 1500; // Starting guess
  tk::real tol = std::max(1e-8, 1e-8 * e); // Stopping condition
  tk::real err;
  std::size_t maxiter = 10;
  std::size_t i(0);
  while (i < maxiter) {
    // Construct f(T) = e(temp) - e
    tk::real f_T = 0.;
    for (std::size_t k = 0; k < m_nspec; k++) {
      f_T += m_Ys[k] * mat_blk[k].compute< EOS::calc_e >(temp);
    }
    f_T -= e;

    // Construct f'(T) = cv(temp)
    tk::real fp_T = 0.;
    for (std::size_t k = 0; k < m_nspec; k++) {
      fp_T += m_Ys[k] * mat_blk[k].compute< EOS::calc_cv >(temp);
    }

    // Calculate next guess
    temp = temp - f_T / fp_T;

    // Check stopping conditions
    err = abs(f_T);
    if (err <= tol) break;
    i++;
    if ( i == maxiter ) {
      Throw("Mixture Newton's Method for temperature failed to converge after iterations "
      + std::to_string(i));
    }
  }

  return temp;
}

