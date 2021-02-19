// *****************************************************************************
/*!
  \file      src/DiffEq/Dirichlet/DirichletCoeffPolicy.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Dirichlet coefficients policies
  \details   This file defines coefficients policy classes for the Dirichlet
    SDE, defined in DiffEq/Dirichlet.h.

    General requirements on Dirichlet SDE coefficients policy classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients, b, S, and kappa. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_t ncomp,
          const std::vector< kw::sde_b::info::expect::type >& b_,
          const std::vector< kw::sde_S::info::expect::type >& S_,
          const std::vector< kw::sde_kappa::info::expect::type >& k_,
          std::vector< kw::sde_b::info::expect::type  >& b,
          std::vector< kw::sde_S::info::expect::type >& S,
          std::vector< kw::sde_kappa::info::expect::type >& k )
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of Dirichlet
        SDEs.
      - Constant references to b_, S_, and k_, which denote three vectors of
        real values used to initialize the parameter vectors of the Dirichlet
        SDEs. The length of the vectors must be equal to the number of
        components given by ncomp.
      - References to b, S, and k, which denote the parameter vectors to be
        initialized based on b_, S_, and k_.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::CoeffPolicyType type() noexcept {
          return ctr::CoeffPolicyType::CONST_COEFF;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef DirichletCoeffPolicy_h
#define DirichletCoeffPolicy_h

#include <brigand/sequences/list.hpp>

#include "Types.hpp"
#include "Walker/Options/CoeffPolicy.hpp"
#include "SystemComponents.hpp"

namespace walker {

//! Dirichlet constant coefficients policity: constants in time
class DirichletCoeffConst {

  public:
    //! Constructor: initialize coefficients
    DirichletCoeffConst(
      tk::ctr::ncomp_t ncomp,
      const std::vector< kw::sde_b::info::expect::type >& b_,
      const std::vector< kw::sde_S::info::expect::type >& S_,
      const std::vector< kw::sde_kappa::info::expect::type >& k_,
      std::vector< kw::sde_b::info::expect::type  >& b,
      std::vector< kw::sde_S::info::expect::type >& S,
      std::vector< kw::sde_kappa::info::expect::type >& k );

    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONST_COEFF; }
};

//! List of all Dirichlet's coefficients policies
using DirichletCoeffPolicies = brigand::list< DirichletCoeffConst >;

} // walker::

#endif // DirichletCoeffPolicy_h
