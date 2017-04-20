// *****************************************************************************
/*!
  \file      src/DiffEq/BetaCoeffPolicy.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Beta coefficients policies
  \details   This file defines coefficients policy classes for the beta SDE,
    defined in DiffEq/Beta.h.

    General requirements on beta SDE coefficients policy classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients, b, S, and kappa. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_type ncomp,
          const std::vector< kw::sde_b::info::expect::type >& b_,
          const std::vector< kw::sde_S::info::expect::type >& S_,
          const std::vector< kw::sde_kappa::info::expect::type >& k_,
          std::vector< kw::sde_b::info::expect::type  >& b,
          std::vector< kw::sde_S::info::expect::type >& S,
          std::vector< kw::sde_kappa::info::expect::type >& k )
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of beta
        SDEs.
      - Constant references to b_, S_, and k_, which denote three vectors of
        real values used to initialize the parameter vectors of the system of
        beta SDEs. The length of the vectors must be equal to the number of
        components given by ncomp.
      - References to b, S, and k, which denote the parameter vectors to be
        initialized based on b_, S_, and k_.

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
#ifndef BetaCoeffPolicy_h
#define BetaCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "SystemComponents.h"
#include "Walker/Options/CoeffPolicy.h"

namespace walker {

//! \brief Beta constant coefficients policity: constants in time
class BetaCoeffConst {

  public:
    //! Constructor: initialize coefficients
    BetaCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_b::info::expect::type >& b_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappa::info::expect::type >& k_,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappa::info::expect::type >& k )
    {
      ErrChk( b_.size() == ncomp,
              "Wrong number of beta SDE parameters 'b'");
      ErrChk( S_.size() == ncomp,
              "Wrong number of beta SDE parameters 'S'");
      ErrChk( k_.size() == ncomp,
              "Wrong number of beta SDE parameters 'kappa'");

      b = b_;
      S = S_;
      k = k_;
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONSTANT; }
};

//! List of all beta's coefficients policies
using BetaCoeffPolicies = boost::mpl::vector< BetaCoeffConst >;

} // walker::

#endif // BetaCoeffPolicy_h
