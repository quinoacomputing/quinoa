// *****************************************************************************
/*!
  \file      src/PDE/Transport/Physics.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     All physics configurations for the scalar transport equations
  \details   This file collects all Physics policy classes for the scalar
    transport equations, defined in PDE/Transport/Transport.h.

    General requirements on Transport Physics policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::PhysicsType type() noexcept {
          return ctr::PhysicsType::ADVECTION;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for Physics policies.

    - Must define the static function _diffusionRhs()_, adding diffusion terms
      to the right hand side.

    - Must define the static function _diffusion_dt()_, computing the minumum
      time step size based on the diffusion term.
*/
// *****************************************************************************
#ifndef TransportPhysics_h
#define TransportPhysics_h

#include <boost/mpl/vector.hpp>

#include "Physics/Advection.h"
#include "Physics/AdvDiff.h"

namespace inciter {

//! List of all Transport Physics policies (defined in the includes above)
using TransportPhysics = boost::mpl::vector< TransportPhysicsAdvection
                                           , TransportPhysicsAdvDiff >;

} // inciter::

#endif // TransportPhysics_h
