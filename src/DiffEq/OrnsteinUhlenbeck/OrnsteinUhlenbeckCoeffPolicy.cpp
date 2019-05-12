// *****************************************************************************
/*!
  \file      src/DiffEq/OrnsteinUhlenbeck/OrnsteinUhlenbeckCoeffPolicy.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Ornstein-Uhlenbeck coefficients policies
  \details   This file defines coefficients policy classes for the diagonal
             Ornstein-Uhlenbeck SDE, defined in
             DiffEq/OrnsteinUhlenbeck/OrnsteinUhlenbeck.h. For general
             requirements on the diagonal Ornstein-Uhlenbeck SDE coefficients
             policy classes see the header file.
*/
// *****************************************************************************

#include "OrnsteinUhlenbeckCoeffPolicy.hpp"

walker::OrnsteinUhlenbeckCoeffConst::OrnsteinUhlenbeckCoeffConst(
  tk::ctr::ncomp_type ncomp,
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
  ErrChk( sigmasq_.size() == ncomp*(ncomp+1)/2,
          "Wrong number of Ornstein-Uhlenbeck SDE parameters 'sigmasq'");
  ErrChk( theta_.size() == ncomp,
          "Wrong number of Ornstein-Uhlenbeck SDE parameters 'theta'");
  ErrChk( mu_.size() == ncomp,
          "Wrong number of Ornstein-Uhlenbeck SDE parameters 'mu'");

  // Prepare upper triangle for Cholesky-decomposition using LAPACK
  sigmasq.resize( ncomp * ncomp );
  std::size_t c = 0;
  for (tk::ctr::ncomp_type i=0; i<ncomp; ++i)
    for (tk::ctr::ncomp_type j=0; j<ncomp; ++j)
      if (i<=j)
        sigmasq[ i*ncomp+j ] = sigmasq_[ c++ ];
      else
        sigmasq[ i*ncomp+j ] = 0.0;

  theta = theta_;
  mu = mu_;
}
