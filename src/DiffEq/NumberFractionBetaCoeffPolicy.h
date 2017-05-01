// *****************************************************************************
/*!
  \file      src/DiffEq/NumberFractionBetaCoeffPolicy.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Number-fraction beta SDE coefficients policies
  \details   This file defines coefficients policy classes for the
    number-fraction beta SDE, defined in DiffEq/NumberFractionBeta.h.

    General requirements on number-fraction beta SDE coefficients policy
    classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients, b, S, kappa, rho2, and rcomma. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_type ncomp,
          const std::vector< kw::sde_b::info::expect::type >& b_,
          const std::vector< kw::sde_S::info::expect::type >& S_,
          const std::vector< kw::sde_kappa::info::expect::type >& k_,
          const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
          const std::vector< kw::sde_rcomma::info::expect::type >& rcomma_,
          std::vector< kw::sde_b::info::expect::type  >& b,
          std::vector< kw::sde_S::info::expect::type >& S,
          std::vector< kw::sde_kappa::info::expect::type >& k,
          std::vector< kw::sde_rho2::info::expect::type >& rho2_,
          std::vector< kw::sde_rcomma::info::expect::type >& rcomma_ );
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of
        number-fraction beta SDEs.
      - Constant references to b_, S_, k_, rho2_, and rcomma_, which denote five
        vectors of real values used to initialize the parameter vectors of the
        system of number-fraction beta SDEs. The length of the vectors must be
        equal to the number of components given by ncomp.
      - References to b, S, k, rho2_, and rcomma, which denote the parameter
        vectors to be initialized based on b_, S_, k_, rho2_, and rcomma_.

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
#ifndef NumberFractionBetaCoeffPolicy_h
#define NumberFractionBetaCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Walker/Options/CoeffPolicy.h"
#include "SystemComponents.h"

namespace walker {

//! \brief Number-fraction beta SDE constant coefficients policity: constants in
//!   time
class NumberFractionBetaCoeffConst {

  public:
    //! Constructor: initialize coefficients
    NumberFractionBetaCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_b::info::expect::type >& b_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappa::info::expect::type >& k_,
      const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
      const std::vector< kw::sde_rcomma::info::expect::type >& rcomma_,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappa::info::expect::type >& k,
      std::vector< kw::sde_rho2::info::expect::type >& rho2,
      std::vector< kw::sde_rcomma::info::expect::type >& rcomma )
    {
      ErrChk( b_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'b'");
      ErrChk( S_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'S'");
      ErrChk( k_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'kappa'");
      ErrChk( rho2_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'rho2'");
      ErrChk( rcomma_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'rcomma'");

      b = b_;
      S = S_;
      k = k_;
      rho2 = rho2_;
      rcomma = rcomma_;
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONSTANT; }
};

//! List of all number-fraction beta's coefficients policies
using NumberFractionBetaCoeffPolicies =
  boost::mpl::vector< NumberFractionBetaCoeffConst >;

} // walker::

#endif // NumberFractionBetaCoeffPolicy_h
