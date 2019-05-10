// *****************************************************************************
/*!
  \file      src/DiffEq/WrightFisher/WrightFisherCoeffPolicy.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Wright-Fisher coefficients policies
  \details   This file defines coefficients policy classes for the Wright-Fisher
    SDE, defined in DiffEq/WrightFisher.h.

    General requirements on the Wright-Fisher SDE coefficients policy classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficient vector, omega. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          tk::ctr::ncomp_t ncomp,
          const std::vector< kw::sde_omega::info::expect::type >& omega_,
          std::vector< kw::sde_omega::info::expect::type >& omega )
      \endcode
      where
      - ncomp denotes the number of scalar components of the system of the
        Wright-Fisher SDE.
      - Constant references to omega_ which denote a vector of real values used
        to initialize the parameter vector of the Wright-Fisher SDE. The length
        of the vector must be equal to the number of components given by ncomp.
      - Reference to omega which denote the parameter vector to be initialized
        based on omega_.

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
#ifndef WrightFisherCoeffPolicy_h
#define WrightFisherCoeffPolicy_h

#include <brigand/sequences/list.hpp>

#include "Types.hpp"
#include "Walker/Options/CoeffPolicy.hpp"
#include "SystemComponents.hpp"

namespace walker {

//! Wright-Fisher constant coefficients policity: constants in time
class WrightFisherCoeffConst {

  public:
    //! Constructor: initialize coefficients
    WrightFisherCoeffConst(
      tk::ctr::ncomp_t ncomp,
      const std::vector< kw::sde_omega::info::expect::type >& omega_,
      std::vector< kw::sde_omega::info::expect::type >& omega );

    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONST_COEFF; }
};

//! List of all Wright-Fisher's coefficients policies
using WrightFisherCoeffPolicies = brigand::list< WrightFisherCoeffConst >;

} // walker::

#endif // WrightFisherCoeffPolicy_h
