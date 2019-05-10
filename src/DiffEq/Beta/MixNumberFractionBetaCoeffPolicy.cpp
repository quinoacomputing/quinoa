// *****************************************************************************
/*!
  \file      src/DiffEq/Beta/MixNumberFractionBetaCoeffPolicy.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mass-fraction beta SDE coefficients policies
  \details   This file defines coefficients policy classes for the mass-fraction
             beta SDE, defined in DiffEq/Beta/MixNumberFractionBeta.h. For
             general requirements on mixture number-fraction beta SDE
             coefficients policy classes see the header file.
*/
// *****************************************************************************

#include "MixNumberFractionBetaCoeffPolicy.hpp"

walker::MixNumFracBetaCoeffDecay::MixNumFracBetaCoeffDecay(
  ncomp_t ncomp,
  const std::vector< kw::sde_bprime::info::expect::type >& bprime_,
  const std::vector< kw::sde_S::info::expect::type >& S_,
  const std::vector< kw::sde_kappaprime::info::expect::type >& kprime_,
  const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
  const std::vector< kw::sde_rcomma::info::expect::type >& rcomma_,
  std::vector< kw::sde_bprime::info::expect::type  >& bprime,
  std::vector< kw::sde_S::info::expect::type >& S,
  std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
  std::vector< kw::sde_rho2::info::expect::type >& rho2,
  std::vector< kw::sde_rcomma::info::expect::type >& rcomma,
  std::vector< kw::sde_b::info::expect::type  >& b,
  std::vector< kw::sde_kappa::info::expect::type >& k )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] bprime_ Vector used to initialize coefficient vector bprime
//! \param[in] S_ Vector used to initialize coefficient vector S
//! \param[in] kprime_ Vector used to initialize coefficient vector kprime
//! \param[in] rho2_ Vector used to initialize coefficient vector rho2
//! \param[in] rcomma_ Vector used to initialize coefficient vector rcomma
//! \param[in,out] bprime Coefficient vector to be initialized
//! \param[in,out] S Coefficient vector to be initialized
//! \param[in,out] kprime Coefficient vector to be initialized
//! \param[in,out] rho2 Coefficient vector to be initialized
//! \param[in,out] rcomma Coefficient vector to be initialized
//! \param[in,out] b Coefficient vector to be initialized
//! \param[in,out] k Coefficient vector to be initialized
// *****************************************************************************
{
  ErrChk( bprime_.size() == ncomp,
    "Wrong number of mix number-fraction beta SDE parameters 'b''");
  ErrChk( S_.size() == ncomp,
    "Wrong number of mix number-fraction beta SDE parameters 'S'");
  ErrChk( kprime_.size() == ncomp,
    "Wrong number of mix number-fraction beta SDE parameters 'kappa''");
  ErrChk( rho2_.size() == ncomp,
    "Wrong number of mix number-fraction beta SDE parameters 'rho2'");
  ErrChk( rcomma_.size() == ncomp,
    "Wrong number of mix number-fraction beta SDE parameters 'rcomma'");

  bprime = bprime_;
  S = S_;
  kprime = kprime_;
  rho2 = rho2_;
  rcomma = rcomma_;

  b.resize( bprime.size() );
  k.resize( kprime.size() );
}

void
walker::MixNumFracBetaCoeffDecay::update(
  char depvar,
  ncomp_t ncomp,
  const std::map< tk::ctr::Product, tk::real >& moments,
  const std::vector< kw::sde_bprime::info::expect::type  >& bprime,
  const std::vector< kw::sde_kappaprime::info::expect::type >& kprime,
  std::vector< kw::sde_b::info::expect::type  >& b,
  std::vector< kw::sde_kappa::info::expect::type >& k ) const
// *****************************************************************************
//  Update coefficients
//! \details This where the mix number-fraction beta SDE is made consistent
//!   with the no-mix and fully mixed limits by specifying the SDE
//!   coefficients, b and kappa as functions of b' and kappa'.
// *****************************************************************************
{
  for (ncomp_t c=0; c<ncomp; ++c) {
    tk::real m = tk::ctr::lookup( tk::ctr::mean(depvar,c), moments );
    tk::real v = tk::ctr::lookup( tk::ctr::variance(depvar,c), moments );

    if (m<1.0e-8 || m>1.0-1.0e-8) m = 0.5;
    if (v<1.0e-8 || v>1.0-1.0e-8) v = 0.5;

    b[c] = bprime[c] * (1.0 - v / m / ( 1.0 - m ));
    k[c] = kprime[c] * v;
  }
}

