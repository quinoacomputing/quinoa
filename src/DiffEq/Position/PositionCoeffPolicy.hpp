// *****************************************************************************
/*!
  \file      src/DiffEq/Position/PositionCoeffPolicy.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Particle position equation coefficients policies
  \details   This file defines coefficients policy classes for the Lagrangian
    particle position equation defined in DiffEq/Position.h.

    General requirements on position equation coefficients policy classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficients. Required signature:
      \code{.cpp}
        CoeffPolicyName( std::array< tk::real, 9 >& dU );
      \endcode
      where _dU_ is an optionally prescribed mean velocity gradient.

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

#include <brigand/sequences/list.hpp>

#include "Types.hpp"
#include "SystemComponents.hpp"
#include "Walker/Options/CoeffPolicy.hpp"

namespace walker {

//! Position equation coefficients policy given by the instantaneous velocity
class PositionInstVel {

  public:
    //! Constructor
    explicit PositionInstVel( std::array< tk::real, 9 >& );

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::INSTANTANEOUS_VELOCITY; }
};

//! \brief Position equation coefficients policy using a prescribed constant
//!   mean velocity gradient for homogeneous shear flow
class PositionConstShear {

  public:
    //! Constructor: prescribe mean shear
    explicit PositionConstShear( std::array< tk::real, 9 >& dU );

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONST_SHEAR; }
};

//! List of all position eq coefficients policies
using PositionCoeffPolicies = brigand::list< PositionInstVel
                                           , PositionConstShear
                                           >;

} // walker::

#endif // PositionCoeffPolicy_h
