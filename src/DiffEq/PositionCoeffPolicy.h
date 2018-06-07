// *****************************************************************************
/*!
  \file      src/DiffEq/PositionCoeffPolicy.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Particle position equation coefficients policies
  \details   This file defines coefficients policy classes for the Lagrangian
    particle position equation defined in DiffEq/Position.h.

    General requirements on position equation coefficients policy classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficient, C0. Required signature:
      \code{.cpp}
        CoeffPolicyName();
      \endcode

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::CoeffPolicyType type() noexcept {
          return ctr::CoeffPolicyType::INSTANTANEOUS_VELOCITY;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef PositionCoeffPolicy_h
#define PositionCoeffPolicy_h

#include <array>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "SystemComponents.h"
#include "Walker/Options/CoeffPolicy.h"

namespace walker {

//! Position equation coefficients policy
class InstantaneousVelocity {

  public:
    //! Constructor: initialize coefficients
    InstantaneousVelocity() {}

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::INSTANTANEOUS_VELOCITY; }
};

//! List of all position eq coefficients policies
using PositionCoeffPolicies = boost::mpl::vector< InstantaneousVelocity
                                                >;

} // walker::

#endif // PositionCoeffPolicy_h
