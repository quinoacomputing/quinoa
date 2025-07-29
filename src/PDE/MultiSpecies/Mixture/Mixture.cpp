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
Mixture::mix_Cv(
  tk::real mix_temp,
  const std::vector< EOS >& mat_blk) const
// *************************************************************************
//! \brief Calculate the mixture Cv
//! \param[in] mix_temp Mixture temperature (provided at call-site, since
//!   it is reconstructed separately
//! \param[in] mat_blk EOS material block
//! \return Mixture Cv using the ideal gas EoS
// *************************************************************************
{
  tk::real mix_Cv = 0.;
  for (std::size_t k = 0; k < m_nspec; k++) {
    mix_Cv += mat_blk[k].compute< EOS::cv >(mix_temp) * m_Ys[k];
  }
  return mix_Cv;
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
//! \param[in] mat_blk EOS material block
//! \return Mixture pressure
// *************************************************************************
{
  // Compute internal energy
  tk::real e = rhoE / mix_density - 0.5 * (u*u + v*v + w*w);

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
      f_T += m_Ys[k] * mat_blk[k].compute< EOS::internalenergy >(temp);
    }
    f_T -= e;

    // Construct f'(T) = cv(temp)
    tk::real fp_T = 0.;
    for (std::size_t k = 0; k < m_nspec; k++) {
      fp_T += m_Ys[k] * mat_blk[k].compute< EOS::cv >(temp);
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

std::vector < tk::real >
Mixture::pressure_prim_partials(
  tk::real mix_density,
  tk::real mix_temp,
  const std::vector< EOS >& mat_blk ) const
// *************************************************************************
//! \brief Calculate mixture pressure partial derivatives with respect to the
//!   primitives.
//! \param[in] mix_density Mixture density (sum of species density)
//! \param[in] mix_temp Mixture temperature
//! \param[in] mat_blk EOS material block
//! \return Mixture pressure
// *************************************************************************
{
  std::vector< tk::real > dpdP(m_nspec + 3 + 1, 0.0);
  std::vector< tk::real > dRdP = mix_R_prim_partials(mix_density, mat_blk);
  for (std::size_t k = 0; k < m_nspec; k++) {
    dpdP[k] = m_mix_R * mix_temp + mix_density * dRdP[k] * mix_temp;
  }
  dpdP[m_nspec + 3] = mix_density*m_mix_R;
  return dpdP;
}

std::vector < tk::real >
Mixture::mix_R_prim_partials(
  tk::real mix_density,
  const std::vector< EOS >& mat_blk) const
// *************************************************************************
//! \brief Calculate mixture gas constant partial derivatives with respect to
//!   the primitives.
//! \param[in] mix_density Mixture density (sum of species density)
//! \param[in] mat_blk EOS material block
//! \return Mixture pressure
// *************************************************************************
{
  std::vector< tk::real > dRdP(m_nspec + 3 + 1, 0.0);
  for (std::size_t k = 0; k < m_nspec; k++) {
    dRdP[k] = mat_blk[k].compute< EOS::gas_constant >() / mix_density
            - m_mix_R / mix_density;
  }
  return dRdP;
}

std::vector < tk::real >
Mixture::mix_Cv_prim_partials(
  tk::real mix_density,
  tk::real mix_temp,
  const std::vector< EOS >& mat_blk) const
// *************************************************************************
//! \brief Calculate mixture specific heat partial derivatives with respect to
//!   the primitives.
//! \param[in] mix_density Mixture density (sum of species density)
//! \param[in] mix_temp Mixture temperature
//! \param[in] mat_blk EOS material block
//! \return Mixture pressure
// *************************************************************************
{
  std::vector< tk::real > dCvdP(m_nspec + 3 + 1, 0.0);
  tk::real mix_Cv = 0.;

  for (std::size_t k = 0; k < m_nspec; k++) {
    mix_Cv += mat_blk[k].compute< EOS::cv >(mix_temp) * m_Ys[k];
  }

  for (std::size_t k = 0; k < m_nspec; k++) {
    dCvdP[k] = mat_blk[k].compute< EOS::cv >(mix_temp) / mix_density
             - mix_Cv / mix_density;
  }
  return dCvdP;
}

std::vector < tk::real >
Mixture::soundspeed_prim_partials(
  tk::real mix_density,
  tk::real mix_temp,
  const std::vector< EOS >& mat_blk ) const
// *************************************************************************
//! \brief Calculate mixture sound speed partial derivatives with respect to
//!   the primitives.
//! \param[in] mix_density Mixture density (sum of species density)
//! \param[in] mix_temp Mixture temperature
//! \param[in] mat_blk EOS material block
//! \return Mixture pressure
// *************************************************************************
{
  std::vector< tk::real > dadP(m_nspec + 3 + 1, 0.0),
                          drhodP(m_nspec + 3 + 1, 0.0);
  auto dpdP = pressure_prim_partials(mix_density, mix_temp, mat_blk);
  auto dCvdP = mix_Cv_prim_partials(mix_density, mix_temp, mat_blk);
  auto dRdP = mix_R_prim_partials(mix_density, mat_blk);
  auto a = frozen_soundspeed(mix_density, mix_temp, mat_blk);
  auto p = pressure(mix_density, mix_temp);
  tk::real mix_Cv(0.), dbetadP(0.);
  for (std::size_t k = 0; k < m_nspec; k++) {
    mix_Cv += mat_blk[k].compute< EOS::cv >(mix_temp) * m_Ys[k];
    drhodP[k] = 1;
  }
  tk::real beta = m_mix_R / mix_Cv;

  // Add constituent partials together
  for (std::size_t k = 0; k < m_nspec + 3 + 1; k++) {
    dbetadP = dRdP[k] / mix_Cv - m_mix_R / ( mix_Cv * mix_Cv ) * dCvdP[k];
    dadP[k] = 0.5 / a * ( (1 + beta) / mix_density * dpdP[k]
                      + p / mix_density * dbetadP
                      - (1 + beta) * p / (mix_density*mix_density) * drhodP[k]
                      );
  }
  return dadP;
}
