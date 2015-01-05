//******************************************************************************
/*!
  \file      src/DiffEq/DiagOUCoeffPolicy.h
  \author    J. Bakosi
  \date      Tue 13 Jan 2015 10:54:30 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Diagonal Ornstein-Uhlenbeck coefficients policies
  \details   Diagonal Ornstein-Uhlenbeck coefficients policies
*/
//******************************************************************************
#ifndef DiagOUCoeffPolicy_h
#define DiagOUCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! Diagonal Ornstein-Uhlenbeck constant coefficients policity: constants in time
struct DiagOUCoeffConst {

  //! Constructor: default for accessing policy name, type, etc.
  DiagOUCoeffConst() = default;
  //! Constructor: initialize coefficients
  DiagOUCoeffConst( tk::ctr::ncomp_type ncomp,
                    const std::vector< tk::real >& sigma_,
                    const std::vector< tk::real >& theta_,
                    const std::vector< tk::real >& mu_,
                    std::vector< tk::real >& sigma,
                    std::vector< tk::real >& theta,
                    std::vector< tk::real >& mu )
  {
    sigma = sigma_;
    theta = theta_;
    mu = mu_;
    ErrChk( sigma.size() == ncomp,
            "Wrong number of diagonal OU SDE parameters 'sigma'");
    ErrChk( theta.size() == ncomp,
            "Wrong number of diagonal OU SDE parameters 'theta'");
    ErrChk( mu.size() == ncomp,
            "Wrong number of diagonal OU SDE parameters 'mu'");
  }

  std::string policy() const noexcept
  { return tk::ctr::CoeffPolicy().name( tk::ctr::CoeffPolicyType::CONSTANT ); }

  tk::ctr::CoeffPolicyType type() const noexcept
  { return tk::ctr::CoeffPolicyType::CONSTANT; }

  //! Function call: no-op for constant coefficients
  void operator()( const tk::real&,
                   std::vector< tk::real >&,
                   std::vector< tk::real >&,
                   std::vector< tk::real >& ) {}
};

//! List of all Ornstein-Uhlenbeck's coefficients policies
using DiagOUCoeffPolicies = boost::mpl::vector< DiagOUCoeffConst >;

} // walker::

#endif // DiagOUCoeffPolicy_h
