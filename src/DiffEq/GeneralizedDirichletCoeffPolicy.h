// *****************************************************************************
/*!
  \file      src/DiffEq/GeneralizedDirichletCoeffPolicy.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Lochner's generalized Dirichlet coefficients policies
  \details   This file defines coefficients policy classes for the generalized
    Dirichlet SDE, defined in DiffEq/GeneralizedDirichlet.h.

    General requirements on generalized Dirichlet SDE coefficients policy
    classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients, b, S, kappa, and c. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_type ncomp,
          const std::vector< kw::sde_b::info::expect::type >& b_,
          const std::vector< kw::sde_S::info::expect::type >& S_,
          const std::vector< kw::sde_kappa::info::expect::type >& k_,
          const std::vector< kw::sde_c::info::expect::type >& c_,
          std::vector< kw::sde_b::info::expect::type  >& b,
          std::vector< kw::sde_S::info::expect::type >& S,
          std::vector< kw::sde_kappa::info::expect::type >& k,
          std::vector< kw::sde_c::info::expect::type >& c )
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of beta
        SDEs.
      - Constant references to b_, S_, k_, and c_, which denote four vectors of
        real values used to initialize the parameter vectors of the generalized
        Dirichlet SDEs. The length of the vectors b_, S_, and kappa_, must be
        equal to the number of components given by ncomp, while the length of
        vector c_ must be ncomp*(ncomp-1)/2.
      - References to b, S, k, and c, which denote the parameter vectors to be
        initialized based on b_, S_, k_, and c_.

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
#ifndef GeneralizedDirichletCoeffPolicy_h
#define GeneralizedDirichletCoeffPolicy_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Walker/Options/CoeffPolicy.h"
#include "SystemComponents.h"

namespace walker {

//! Generalized Dirichlet constant coefficients policity: constants in time
class GeneralizedDirichletCoeffConst {

  public:
    //! Constructor: initialize coefficients
    GeneralizedDirichletCoeffConst(
      tk::ctr::ncomp_type ncomp,
      const std::vector< kw::sde_b::info::expect::type >& b_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappa::info::expect::type >& k_,
      const std::vector< kw::sde_c::info::expect::type >& c_,
      std::vector< kw::sde_b::info::expect::type >& b,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappa::info::expect::type >& k,
      std::vector< kw::sde_c::info::expect::type >& c )
    {
      ErrChk( b_.size() == ncomp,
              "Wrong number of generalized Dirichlet SDE parameters 'b'");
      ErrChk( S_.size() == ncomp,
              "Wrong number of generalized Dirichlet SDE parameters 'S'");
      ErrChk( k_.size() == ncomp,
              "Wrong number of generalized Dirichlet SDE parameters 'kappa'");
      ErrChk( c_.size() == ncomp*(ncomp-1)/2,
              "Wrong number of generalized Dirichlet SDE parameters 'c'");

      b = b_;
      S = S_;
      k = k_;
      c = c_;
    }

    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONSTANT; }
};

//! List of all generalized Dirichlet's coefficients policies
using GeneralizedDirichletCoeffPolicies =
  boost::mpl::vector< GeneralizedDirichletCoeffConst >;

} // walker::

#endif // GeneralizedDirichletCoeffPolicy_h
