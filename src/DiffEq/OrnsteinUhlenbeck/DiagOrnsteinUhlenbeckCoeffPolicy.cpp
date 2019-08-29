// *****************************************************************************
/*!
  \file      src/DiffEq/OrnsteinUhlenbeck/DiagOrnsteinUhlenbeckCoeffPolicy.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Ornstein-Uhlenbeck coefficients policies
  \details   This file defines coefficients policy classes for the diagonal
             Ornstein-Uhlenbeck SDE, defined in
             DiffEq/OrnsteinUhlenbeck/DiagOrnsteinUhlenbeck.h. For general
             requirements on the diagonal Ornstein-Uhlenbeck SDE coefficients
             policy classes see the header file.
*/
// *****************************************************************************

#include "DiagOrnsteinUhlenbeckCoeffPolicy.hpp"

walker::DiagOrnsteinUhlenbeckCoeffConst::DiagOrnsteinUhlenbeckCoeffConst(
  tk::ctr::ncomp_t ncomp,
  const std::vector< kw::sde_sigmasq::info::expect::type >& sigmasq_,
  const std::vector< kw::sde_theta::info::expect::type >& theta_,
  const std::vector< kw::sde_mu::info::expect::type >& mu_,
  std::vector< kw::sde_sigmasq::info::expect::type >& sigmasq,
  std::vector< kw::sde_theta::info::expect::type >& theta,
  std::vector< kw::sde_mu::info::expect::type >& mu )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] sigmasq_ Vector used to initialize coefficient vector sigmasq
//! \param[in] theta_ Vector used to initialize coefficient vector theta
//! \param[in] mu_ Vector used to initialize coefficient vector mu
//! \param[in,out] sigmasq Coefficient vector to be initialized
//! \param[in,out] theta Coefficient vector to be initialized
//! \param[in,out] mu Coefficient vector to be initialized
// *****************************************************************************
{
  ErrChk( sigmasq_.size() == ncomp,
   "Wrong number of diagonal Ornstein-Uhlenbeck SDE parameters 'sigmasq'");
  ErrChk( theta_.size() == ncomp,
   "Wrong number of diagonal Ornstein_uhlenbeck SDE parameters 'theta'");
  ErrChk( mu_.size() == ncomp,
   "Wrong number of diagonal Ornstein_uhlenbeck SDE parameters 'mu'");

  sigmasq = sigmasq_;
  theta = theta_;
  mu = mu_;
}
