//******************************************************************************
/*!
  \file      src/DiffEq/OrnsteinUhlenbeckCoeffPolicy.h
  \author    J. Bakosi
  \date      Mon 26 Jan 2015 11:52:59 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Ornstein-Uhlenbeck coefficients policies
  \details   This file defines coefficients policy classes for the
    Ornstein-Uhlenbeck SDE, defined in DiffEq/OrnsteinUhlenbeck.h.

    General requirements on the Ornstein-Uhlenbeck SDE coefficients policy
    classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients, sigma, theta, and mu. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_type ncomp,
          const std::vector< kw::sde_sigma::info::expect::type >& sigma_,
          const std::vector< kw::sde_theta::info::expect::type >& theta_,
          const std::vector< kw::sde_mu::info::expect::type >& mu_,
          std::vector< kw::sde_sigma::info::expect::type >& sigma,
          std::vector< kw::sde_theta::info::expect::type >& theta,
          std::vector< kw::sde_mu::info::expect::type >& mu )
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of the
        Ornstein-Uhlenbeck SDEs.
      - Constant references to sigma_, theta_, and mu_, which denote three
        vectors of real values used to initialize the parameter vectors of the
        system of Ornstein-Uhlenbeck SDEs. The length of the vector sigma_ must
        be equal to ncomp*(ncomp+1)/2, while the number of components of the
        vectors theta_, and mu_, must be equal to ncomp.
      - References to sigma, theta, and mu, which denote the parameter vectors
        to be initialized based on sigma_, theta_, and mu_.

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
      OrnsteinUhlenbeck::initialize(), performing pre-lookup of the locations of
     the statistical moments required by the given model. Required
      signature:
      \code{.cpp}
        void lookup( const tk::Statistics& stat, char depvar )
      \endcode
      where _stat_ is the Statistics object, allowing access to the location of
      the various moments in memory, and _depvar_ is the dependent variable
      associated with the Ornstein-Uhlenbeck SDE, given in the control file by
      the user.
*/
//******************************************************************************
#ifndef OrnsteinUhlenbeckCoeffPolicy_h
#define OrnsteinUhlenbeckCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Options/CoeffPolicy.h>

namespace walker {

//! Ornstein-Uhlenbeck constant coefficients policity: constants in time
class OrnsteinUhlenbeckCoeffConst {

  public:
    //! Constructor: initialize coefficients
    OrnsteinUhlenbeckCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_sigma::info::expect::type >& sigma_,
      const std::vector< kw::sde_theta::info::expect::type >& theta_,
      const std::vector< kw::sde_mu::info::expect::type >& mu_,
      std::vector< kw::sde_sigma::info::expect::type >& sigma,
      std::vector< kw::sde_theta::info::expect::type >& theta,
      std::vector< kw::sde_mu::info::expect::type >& mu )
    {
      ErrChk( sigma_.size() == ncomp*(ncomp+1)/2,
              "Wrong number of Ornstein-Uhlenbeck SDE parameters 'sigma'");
      sigma.resize( ncomp * ncomp );
      std::size_t c = 0;
      for (tk::ctr::ncomp_type i=0; i<ncomp; ++i)
        for (tk::ctr::ncomp_type j=0; j<ncomp; ++j)
          if (i<=j) sigma[ i*ncomp+j ] = sigma_[ c++ ];

      theta = theta_;
      mu = mu_;
      ErrChk( theta.size() == ncomp,
              "Wrong number of Ornstein-Uhlenbeck SDE parameters 'theta'");
      ErrChk( mu.size() == ncomp,
              "Wrong number of Ornstein-Uhlenbeck SDE parameters 'mu'");
    }

    static tk::ctr::CoeffPolicyType type() noexcept
    { return tk::ctr::CoeffPolicyType::CONSTANT; }

    //! Lookup statistical moments required: no-op for constant coefficients
    void lookup( const tk::Statistics& stat, char ) {}

    //! Function call: no-op for constant coefficients
    void operator()( const tk::real&,
                     std::vector< tk::real >&,
                     std::vector< tk::real >&,
                     std::vector< tk::real >& ) {}
};

//! List of all Ornstein-Uhlenbeck's coefficients policies
using OrnsteinUhlenbeckCoeffPolicies =
  boost::mpl::vector< OrnsteinUhlenbeckCoeffConst >;

} // walker::

#endif // OrnsteinUhlenbeckCoeffPolicy_h
