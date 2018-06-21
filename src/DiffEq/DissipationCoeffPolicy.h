// *****************************************************************************
/*!
  \file      src/DiffEq/DissipationCoeffPolicy.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Particle dissipation equation coefficients policies
  \details   This file defines coefficients policy classes for the Lagrangian
    particle dissipation equation defined in DiffEq/Dissipation.h.

    General requirements on dissipation equation coefficients policy classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficient, C0. Required signature:
      \code{.cpp}
        CoeffPolicyName();
      \endcode

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
#ifndef DissipationCoeffPolicy_h
#define DissipationCoeffPolicy_h

#include <array>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "SystemComponents.h"
#include "Walker/Options/CoeffPolicy.h"

namespace walker {

//! Dissipation equation coefficients policy
class DissipationCoeffConstant {

  public:
    //! Constructor: initialize coefficients
    DissipationCoeffConstant() {}

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONSTANT; }
};

//! List of all dissipation eq coefficients policies
using DissipationCoeffPolicies = boost::mpl::vector< DissipationCoeffConstant
                                                   >;

} // walker::

#endif // DissipationCoeffPolicy_h
