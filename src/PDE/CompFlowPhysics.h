// *****************************************************************************
/*!
  \file      src/PDE/CompFlowPhysics.h
  \author    J. Bakosi
  \date      Fri 22 Jul 2016 02:49:02 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Physics configurations for the compressible flow equations
  \details   This file defines policy classes for the compressible flow
    equations, defined in PDE/CompFlow.h.

    General requirements on flow equations problem policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::PhysicsType type() noexcept {
          return ctr::PhysicsType::NAVIERSTOKES;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef CompFlowPhysics_h
#define CompFlowPhysics_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Physics.h"

namespace inciter {

//! CompFlow system of PDEs problem: Navier-Stokes (viscous flow)
class CompFlowPhysicsNavierStokes {
  public:


    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::NAVIERSTOKES; }
};

//! CompFlow system of PDEs problem: Euler (inviscid flow)
class CompFlowPhysicsEuler {
  public:


    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::EULER; }
};

//! List of all CompFlow problem policies
using CompFlowPhysics = boost::mpl::vector< CompFlowPhysicsNavierStokes,
                                            CompFlowPhysicsEuler >;

} // inciter::

#endif // CompFlowPhysics_h
