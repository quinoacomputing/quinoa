// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Physics.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     All physics configurations for the compressible flow equations
  \details   This file collects all Physics policy classes for the compressible
    flow equations, defined in PDE/CompFlow/CompFlow.h.

    General requirements on CompFlow Physics policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::PhysicsType type() noexcept {
          return ctr::PhysicsType::NAVIERSTOKES;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for Physics policies.

    - Must define the static function _viscousRhs()_, adding the viscous terms
      to the right hand side.

    - Must define the static function _viscous_dt()_, computing the minumum time
      step size based on the viscous term.

    - Must define the static function _conductRhs()_, adding the heat conduction
      terms to the right hand side.

    - Must define the static function _conduct_dt()_, computing the minumum time
      step size based on the heat diffusion term.
*/
// *****************************************************************************
#ifndef CompFlowPhysics_h
#define CompFlowPhysics_h

#include <boost/mpl/vector.hpp>

#include "Physics/Euler.h"
#include "Physics/NavierStokes.h"

namespace inciter {

//! List of all CompFlow Physics policies (defined in the includes above)
using CompFlowPhysics = boost::mpl::vector< CompFlowPhysicsEuler
                                          , CompFlowPhysicsNavierStokes >;

} // inciter::

#endif // CompFlowPhysics_h
