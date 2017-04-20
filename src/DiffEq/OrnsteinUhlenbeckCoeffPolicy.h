// *****************************************************************************
/*!
  \file      src/DiffEq/OrnsteinUhlenbeckCoeffPolicy.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Ornstein-Uhlenbeck coefficients policies
  \details   This file defines coefficients policy classes for the
    Ornstein-Uhlenbeck SDE, defined in DiffEq/OrnsteinUhlenbeck.h.

    General requirements on the Ornstein-Uhlenbeck SDE coefficients policy
    classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients, sigmasq, theta, and mu. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_type ncomp,
          const std::vector< kw::sde_sigmasq::info::expect::type >& sigmasq_,
          const std::vector< kw::sde_theta::info::expect::type >& theta_,
          const std::vector< kw::sde_mu::info::expect::type >& mu_,
          std::vector< kw::sde_sigmasq::info::expect::type >& sigmasq,
          std::vector< kw::sde_theta::info::expect::type >& theta,
          std::vector< kw::sde_mu::info::expect::type >& mu )
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of the
        Ornstein-Uhlenbeck SDEs.
      - Constant references to sigmasq_, theta_, and mu_, which denote three
        vectors of real values used to initialize the parameter vectors of the
        system of Ornstein-Uhlenbeck SDEs. The length of the vector sigmasq_
        must be equal to ncomp*(ncomp+1)/2, while the number of components of
        the vectors theta_, and mu_, must be equal to ncomp.
      - References to sigmasq, theta, and mu, which denote the parameter vectors
        to be initialized based on sigmasq_, theta_, and mu_.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::CoeffPolicyType type() noexcept {
          return ctr::CoeffPolicyType::CONSTANT;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef OrnsteinUhlenbeckCoeffPolicy_h
#define OrnsteinUhlenbeckCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Walker/Options/CoeffPolicy.h"
#include "SystemComponents.h"

namespace walker {

//! Ornstein-Uhlenbeck constant coefficients policity: constants in time
class OrnsteinUhlenbeckCoeffConst {

  public:
    //! Constructor: initialize coefficients
    OrnsteinUhlenbeckCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_sigmasq::info::expect::type >& sigmasq_,
      const std::vector< kw::sde_theta::info::expect::type >& theta_,
      const std::vector< kw::sde_mu::info::expect::type >& mu_,
      std::vector< kw::sde_sigmasq::info::expect::type >& sigmasq,
      std::vector< kw::sde_theta::info::expect::type >& theta,
      std::vector< kw::sde_mu::info::expect::type >& mu )
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

    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONSTANT; }
};

//! List of all Ornstein-Uhlenbeck's coefficients policies
using OrnsteinUhlenbeckCoeffPolicies =
  boost::mpl::vector< OrnsteinUhlenbeckCoeffConst >;

} // walker::

#endif // OrnsteinUhlenbeckCoeffPolicy_h
