// *****************************************************************************
/*!
  \file      src/DiffEq/WrightFisher/WrightFisherCoeffPolicy.C
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Wright-Fisher coefficients policies
  \details   This file defines coefficients policy classes for the Wright-Fisher
             SDE, defined in DiffEq/WrightFisher/WrightFisher.h. For general
             requirements on Wright-Fisher SDE coefficients policy classes see
             the header file.
*/
// *****************************************************************************

#include "WrightFisherCoeffPolicy.h"

walker::WrightFisherCoeffConst::WrightFisherCoeffConst(
  tk::ctr::ncomp_type ncomp,
  const std::vector< kw::sde_omega::info::expect::type >& omega_,
  std::vector< kw::sde_omega::info::expect::type >& omega )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] ncomp Number of scalar components in this SDE system
//! \param[in] omega_ Vector used to initialize coefficient vector omega
//! \param[in,out] omega Coefficient vector to be initialized
// *****************************************************************************
{
  ErrChk( omega_.size() == ncomp,
          "Wrong number of Wright-Fisher SDE parameters 'omega'");

  omega = omega_;
}
