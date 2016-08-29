// *****************************************************************************
/*!
  \file      src/PDE/TransportPhysics.h
  \author    J. Bakosi
  \date      Mon 29 Aug 2016 01:13:31 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Physics configurations for a system of transport equations
  \details   This file defines policy classes for transport equations,
    defined in PDE/Transport.h.

    General requirements on transport equation physics policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::PhysicsType type() noexcept {
          return ctr::PhysicsType::Advection;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef TransportPhysics_h
#define TransportPhysics_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Physics.h"

namespace inciter {

//! Transport equation system of PDEs problem: advection
class TransportPhysicsAdvection {
  public:


    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::ADVECTION; }
};

//! List of all Transport equation problem policies
using TransportPhysics = boost::mpl::vector< TransportPhysicsAdvection >;

} // inciter::

#endif // AdvDifTransport_h
