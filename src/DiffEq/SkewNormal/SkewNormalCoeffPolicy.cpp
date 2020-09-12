// *****************************************************************************
/*!
  \file      src/DiffEq/SkewNormal/SkewNormalCoeffPolicy.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     SkewNormal coefficients policies
  \details   This file defines coefficients policy classes for the skew-normal
             SDE, defined in DiffEq/SkewNormal/SkewNormal.h. For general
             requirements on skew-normal SDE coefficients policy classes see the
             header file.
*/
// *****************************************************************************

#include "SkewNormalCoeffPolicy.hpp"

walker::SkewNormalCoeffConst::SkewNormalCoeffConst(
  tk::ctr::ncomp_t ncomp,
  const std::vector< kw::sde_T::info::expect::type >& timescale_,
  const std::vector< kw::sde_sigmasq::info::expect::type >& sigmasq_,
  const std::vector< kw::sde_lambda::info::expect::type >& lambda_,
  std::vector< kw::sde_T::info::expect::type >& timescale,
  std::vector< kw::sde_sigmasq::info::expect::type >& sigmasq,
  std::vector< kw::sde_lambda::info::expect::type >& lambda )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] timescale_ Vector used to initialize coefficient vector timescale
//! \param[in] sigmasq_ Vector used to initialize coefficient vector sigmasq
//! \param[in] lambda_ Vector used to initialize coefficient vector lambda
//! \param[in,out] timescale Coefficient vector to be initialized
//! \param[in,out] sigmasq Coefficient vector to be initialized
//! \param[in,out] lambda Coefficient vector to be initialized
// *****************************************************************************
{
  ErrChk( timescale_.size() == ncomp,
    "Wrong number of diagonal Skew-normal SDE parameters 'timescale'");
  ErrChk( sigmasq_.size() == ncomp,
    "Wrong number of diagonal Skew-normal SDE parameters 'sigmasq'");
  ErrChk( lambda_.size() == ncomp,
    "Wrong number of diagonal Skew-normal SDE parameters 'lambda'");

  timescale = timescale_;
  sigmasq = sigmasq_;
  lambda = lambda_;
}
