// *****************************************************************************
/*!
  \file      src/DiffEq/Velocity/VelocityCoeffPolicy.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Velocity coefficients policies
  \details   This file defines coefficients policy classes for the velocity
             SDE, defined in DiffEq/Velocity/Velocity.h. For general 
             requirements on velocity SDE coefficients policy classes see the
             header file.
*/
// *****************************************************************************

#include "VelocityCoeffPolicy.hpp"
#include "Table.hpp"

walker::VelocityCoeffConstShear::VelocityCoeffConstShear(
  kw::sde_c0::info::expect::type C0_,
  kw::sde_c0::info::expect::type& C0,
  std::array< tk::real, 9 >& dU ) : m_dU( {{ 0.0, 1.0, 0.0,
                                             0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0 }} )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] C0_ Value of C0 parameter in the Langevin model
//! \param[in,out] C0 Value of to set the C0 parameter in the Langevin model
//! \param[in,out] dU Prescribed mean velocity gradient
// *****************************************************************************
{
  C0 = C0_;
  dU = m_dU;
}

void
walker::VelocityCoeffConstShear::update(
  char depvar,
  char dissipation_depvar,
  const std::map< tk::ctr::Product, tk::real >& moments,
  const tk::Table<1>&,
  ctr::DepvarType solve,
  ctr::VelocityVariantType variant,
  kw::sde_c0::info::expect::type C0,
  tk::real,
  tk::real& eps,
  std::array< tk::real, 9 >& G ) const
// *****************************************************************************
//  Update the model coefficients prescribing shear
//! \param[in] depvar Dependent variable for of this SDE
//! \param[in] dissipation_depvar Dependent variable for coupled dissipation eq
//! \param[in] moments Map of statistical moments
//! \param[in] solve Configured dependent variable to solve for
//! \param[in] variant Velocity model variant configured
//! \param[in] C0 Coefficient C0 in the Langevin model
//! \param[in,out] eps Dissipation rate of turbulent kinetic energy
//! \param[in,out] G Coefficient tensor (3x3) in the Langevin equation
//! \details Update the dissipation rate (eps) and G_{ij} based on the
//!   turbulent kinetic energy (k) for a prescribed honmogeneous shear flow.
// *****************************************************************************
{
  using tk::ctr::lookup;
  using tk::ctr::mean;

  // Compute turbulent kinetic energy
  auto rs = reynoldsStress( depvar, solve, moments );

  // Compute turbulent kinetic energy
  auto k = (rs[0] + rs[1] + rs[2]) / 2.0;

  // Access mean turbulence frequency
  tk::real O = lookup( mean(dissipation_depvar,0), moments );

  // compute turbulent kinetic energy dissipation rate
  eps = O*k;

  // update drift tensor based on the Langevin model variant configured
  if (variant == ctr::VelocityVariantType::SLM)     // simplified
    G = slm( O, C0 );
  else if (variant == ctr::VelocityVariantType::GLM)// generalized
    G = glm( O, C0, rs, m_dU );
  else Throw( "Velocity variant type not implemented" );
}

walker::VelocityCoeffStationary::VelocityCoeffStationary(
  kw::sde_c0::info::expect::type C0_,
  kw::sde_c0::info::expect::type& C0,
  std::array< tk::real, 9 >& dU )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] C0_ Value of C0 parameter in the Langevin model
//! \param[in,out] C0 Value of to set the C0 parameter in the Langevin model
//! \param[in,out] dU Prescribed mean velocity gradient
//! \details Prescribe no shear. The value of C0 is insignificant for a forced
//!   stationary velocity PDF because drift and diffusion are in balance, so
//!   that dk/dt = 0.
// *****************************************************************************
{
  C0 = C0_;
  dU = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
}

void
walker::VelocityCoeffStationary::update(
  char,
  char,
  const std::map< tk::ctr::Product, tk::real >&,
  const tk::Table<1>&,
  ctr::DepvarType,
  ctr::VelocityVariantType,
  kw::sde_c0::info::expect::type C0,
  tk::real,
  tk::real& eps,
  std::array< tk::real, 9 >& G ) const
// *****************************************************************************
//  Update the model coefficients forcing a statistically stationary PDF
//! \param[in] C0 Coefficient C0 in the Langevin model, should not affect the
//!   solution for forced velocity PDF
//! \param[in,out] eps Dissipation rate of turbulent kinetic energy, force = 1
//! \param[in,out] G Coefficient tensor (3x3) in the Langevin equation
//! \details Update the dissipation rate (eps) and G_{ij} so that the velocity
//!   PDF is stationary. The value of C0 is insignificant for a forced
//!   stationary velocity PDF because drift and diffusion are in balance, so
//!   that dk/dt = 0.
// *****************************************************************************
{
  // Override turbulent kinetic energy to keep the velocity PDF exactly
  // stationary
  tk::real k = 1.0;

  // Do not couple a dissipation eq for forced velocity PDF, but set to unity
  // and keep the PDF stationary.
  tk::real O = 1.0;

  // Compute turbulent kinetic energy dissipation rate
  eps = O*k;

  // Update drift tensor to force the velocity PDF stationary. Note that his is
  // NOT the simplified or generalized Langevin model, but a modification to
  // keep the PDF stationary, see Pope, Turbulent Flows, 2000, Eq.12.100.
  G.fill( 0.0 );
  G[0] = G[4] = G[8] = -0.75*C0*O;
}

walker::VelocityCoeffHydroTimeScale::VelocityCoeffHydroTimeScale(
  kw::sde_c0::info::expect::type C0_,
  kw::sde_c0::info::expect::type& C0,
  std::array< tk::real, 9 >& dU )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] C0_ Value of C0 parameter in the Langevin model
//! \param[in,out] C0 Value of to set the C0 parameter in the Langevin model
//! \param[in,out] dU Prescribed mean velocity gradient
// *****************************************************************************
{
  C0 = C0_;
  dU = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
}

void
walker::VelocityCoeffHydroTimeScale::update(
  char depvar,
  char,
  const std::map< tk::ctr::Product, tk::real >& moments,
  const tk::Table<1>& hts,
  ctr::DepvarType solve,
  ctr::VelocityVariantType,
  kw::sde_c0::info::expect::type C0,
  tk::real t,
  tk::real& eps,
  std::array< tk::real, 9 >& G ) const
// *****************************************************************************
//  Update the model coefficients sampling the hydrodynamics time scale from a
//  prescribed function table
//! \param[in] depvar Dependent variable for of this SDE
//! \param[in] moments Map of statistical moments
//! \param[in] hts Table to take hydrodynamics time scale from
//! \param[in] solve Configured dependent variable to solve for
//! \param[in] C0 Coefficient C0 in the Langevin model
//! \param[in] t Physical time to sample hydrodynamics time scale at
//! \param[in,out] eps Dissipation rate of turbulent kinetic energy
//! \param[in,out] G Coefficient tensor (3x3) in the Langevin equation
//! \details Update the dissipation rate (eps) based on eps/k (from DNS) and the
//!   turbulent kinetic energy (k) (from the SDE)
// *****************************************************************************
{
  // Compute turbulent kinetic energy
  auto k = tke( depvar, solve, moments );

  // Sample the inverse hydrodynamics timescale at time t
  auto ts = tk::sample<1>( t, hts )[ 0 ];  // eps/k

  // compute turbulent kinetic energy dissipation rate
  eps = ts * k;

  // update drift tensor based on the simplified Langevin model
  G.fill( 0.0 );
  G[0] = G[4] = G[8] = -(0.5+0.75*C0) * ts;
}
