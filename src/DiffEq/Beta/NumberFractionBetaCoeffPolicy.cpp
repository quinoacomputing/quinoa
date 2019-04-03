// *****************************************************************************
/*!
  \file      src/DiffEq/Beta/NumberFractionBetaCoeffPolicy.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mass-fraction beta SDE coefficients policies
  \details   This file defines coefficients policy classes for the mass-fraction
             beta SDE, defined in DiffEq/Beta/NumberFractionBeta.h. For
             general requirements on number-fraction beta SDE coefficients
             policy classes see the header file.
*/
// *****************************************************************************

#include "NumberFractionBetaCoeffPolicy.hpp"

walker::NumFracBetaCoeffConst::NumFracBetaCoeffConst(
  tk::ctr::ncomp_type ncomp,
  const std::vector< kw::sde_b::info::expect::type >& b_,
  const std::vector< kw::sde_S::info::expect::type >& S_,
  const std::vector< kw::sde_kappa::info::expect::type >& k_,
  const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
  const std::vector< kw::sde_rcomma::info::expect::type >& rcomma_,
  std::vector< kw::sde_b::info::expect::type >& b,
  std::vector< kw::sde_S::info::expect::type >& S,
  std::vector< kw::sde_kappa::info::expect::type >& k,
  std::vector< kw::sde_rho2::info::expect::type >& rho2,
  std::vector< kw::sde_rcomma::info::expect::type >& rcomma )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] b_ Vector used to initialize coefficient vector b
//! \param[in] S_ Vector used to initialize coefficient vector S
//! \param[in] k_ Vector used to initialize coefficient vector k
//! \param[in] rho2_ Vector used to initialize coefficient vector rho2
//! \param[in] rcomma_ Vector used to initialize coefficient vector rcomma
//! \param[in,out] b Coefficient vector to be initialized
//! \param[in,out] S Coefficient vector to be initialized
//! \param[in,out] k Coefficient vector to be initialized
//! \param[in,out] rho2 Coefficient vector to be initialized
//! \param[in,out] rcomma Coefficient vector to be initialized
// *****************************************************************************
{
  ErrChk( b_.size() == ncomp,
          "Wrong number of number-fraction beta SDE parameters 'b'");
  ErrChk( S_.size() == ncomp,
          "Wrong number of number-fraction beta SDE parameters 'S'");
  ErrChk( k_.size() == ncomp,
          "Wrong number of number-fraction beta SDE parameters 'kappa'");
  ErrChk( rho2_.size() == ncomp,
          "Wrong number of number-fraction beta SDE parameters 'rho2'");
  ErrChk( rcomma_.size() == ncomp,
          "Wrong number of number-fraction beta SDE parameters 'rcomma'");

  b = b_;
  S = S_;
  k = k_;
  rho2 = rho2_;
  rcomma = rcomma_;
}
