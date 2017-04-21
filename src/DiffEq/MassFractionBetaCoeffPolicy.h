// *****************************************************************************
/*!
  \file      src/DiffEq/MassFractionBetaCoeffPolicy.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Mass-fraction beta SDE coefficients policies
  \details   This file defines coefficients policy classes for the
    number-fraction beta SDE, defined in DiffEq/MassFractionBeta.h.

    General requirements on number-fraction beta SDE coefficients policy
    classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients, b, S, kappa, rho2, and r. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_type ncomp,
          const std::vector< kw::sde_b::info::expect::type >& b_,
          const std::vector< kw::sde_S::info::expect::type >& S_,
          const std::vector< kw::sde_kappa::info::expect::type >& k_,
          const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
          const std::vector< kw::sde_r::info::expect::type >& r_,
          std::vector< kw::sde_b::info::expect::type  >& b,
          std::vector< kw::sde_S::info::expect::type >& S,
          std::vector< kw::sde_kappa::info::expect::type >& k,
          std::vector< kw::sde_rho2::info::expect::type >& rho2_,
          std::vector< kw::sde_r::info::expect::type >& r_ );
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of
        number-fraction beta SDEs.
      - Constant references to b_, S_, k_, rho2_, and r_, which denote five
        vectors of real values used to initialize the parameter vectors of the
        system of number-fraction beta SDEs. The length of the vectors must be
        equal to the number of components given by ncomp.
      - References to b, S, k, rho2_, and r, which denote the parameter vectors
        to be initialized based on b_, S_, k_, rho2_, and r_.

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
#ifndef MassFractionBetaCoeffPolicy_h
#define MassFractionBetaCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Walker/Options/CoeffPolicy.h"
#include "SystemComponents.h"

namespace walker {

//! \brief Mass-fraction beta SDE constant coefficients policity: constants in
//!   time
class MassFractionBetaCoeffConst {

  public:
    //! Constructor: initialize coefficients
    MassFractionBetaCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_b::info::expect::type >& b_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappa::info::expect::type >& k_,
      const std::vector< kw::sde_rho2::info::expect::type >& rho2_,
      const std::vector< kw::sde_r::info::expect::type >& r_,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappa::info::expect::type >& k,
      std::vector< kw::sde_rho2::info::expect::type >& rho2,
      std::vector< kw::sde_r::info::expect::type >& r )
    {
      ErrChk( b_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'b'");
      ErrChk( S_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'S'");
      ErrChk( k_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'kappa'");
      ErrChk( rho2_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'rho2'");
      ErrChk( r_.size() == ncomp,
              "Wrong number of number-fraction beta SDE parameters 'r'");

      b = b_;
      S = S_;
      k = k_;
      rho2 = rho2_;
      r = r_;
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONSTANT; }
};

//! List of all mass-fraction beta's coefficients policies
using MassFractionBetaCoeffPolicies =
  boost::mpl::vector< MassFractionBetaCoeffConst >;

} // walker::

#endif // MassFractionBetaCoeffPolicy_h
