//******************************************************************************
/*!
  \file      src/DiffEq/DiagOrnsteinUhlenbeckCoeffPolicy.h
  \author    J. Bakosi
  \date      Sat 07 Feb 2015 07:45:00 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Diagonal Ornstein-Uhlenbeck coefficients policies
  \details   This file defines coefficients policy classes for the diagonal
    Ornstein-Uhlenbeck SDE, defined in DiffEq/DiagOrnsteinUhlenbeck.h.

    General requirements on the diagonal Ornstein-Uhlenbeck SDE coefficients
    policy classes:

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
        diagonal Ornstein-Uhlenbeck SDEs.
      - Constant references to sigmasq_, theta_, and mu_, which denote three
        vectors of real values used to initialize the parameter vectors of the
        system of diagonal Ornstein-Uhlenbeck SDEs. The length of the vectors
        must be equal to the number of components given by ncomp.
      - References to sigmasq, theta, and mu, which denote the parameter vectors
        to be initialized based on sigmasq_, theta_, and mu_.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static tk::ctr::CoeffPolicyType type() noexcept {
          return tk::ctr::CoeffPolicyType::CONSTANT;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.

    - Must define the function _lookup()_, called from
      DiagOrnsteinUhlenbeck::initialize(), performing pre-lookup of the
      locations of the statistical moments required by the given model. Required
      signature:
      \code{.cpp}
        void lookup( const tk::Statistics& stat, char depvar )
      \endcode
      where _stat_ is the Statistics object, allowing access to the location of
      the various moments in memory, and _depvar_ is the dependent variable
      associated with the diagonal Ornstein-Uhlenbeck SDE, given in the control
      file by the user.
*/
//******************************************************************************
#ifndef DiagOrnsteinUhlenbeckCoeffPolicy_h
#define DiagOrnsteinUhlenbeckCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! Diagonal Ornstein-Uhlenbeck constant coefficients policity: constants in time
class DiagOrnsteinUhlenbeckCoeffConst {

  public:
    //! Constructor: initialize coefficients
    DiagOrnsteinUhlenbeckCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_sigmasq::info::expect::type >& sigmasq_,
      const std::vector< kw::sde_theta::info::expect::type >& theta_,
      const std::vector< kw::sde_mu::info::expect::type >& mu_,
      std::vector< kw::sde_sigmasq::info::expect::type >& sigmasq,
      std::vector< kw::sde_theta::info::expect::type >& theta,
      std::vector< kw::sde_mu::info::expect::type >& mu )
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

    //! Coefficients policy type accessor
    static tk::ctr::CoeffPolicyType type() noexcept
    { return tk::ctr::CoeffPolicyType::CONSTANT; }

    //! Lookup statistical moments required: no-op for constant coefficients
    void lookup( const tk::Statistics&, char ) {}

    //! Function call: no-op for constant coefficients
    void operator()( const tk::real&,
                     std::vector< tk::real >&,
                     std::vector< tk::real >&,
                     std::vector< tk::real >& ) {}
};

//! List of all Ornstein-Uhlenbeck's coefficients policies
using DiagOrnsteinUhlenbeckCoeffPolicies =
  boost::mpl::vector< DiagOrnsteinUhlenbeckCoeffConst >;

} // walker::

#endif // DiagOrnsteinUhlenbeckCoeffPolicy_h
