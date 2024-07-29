// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Physics/CG.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Physics configurations for compressible single-material flow using
    continuous Galerkin discretization
  \details   This file configures all Physics policy classes for compressible
    single-material flow implementied using continuous Galerkin discretization,
    defined in PDE/CompFlow/CGCompFlow.h.

    General requirements on CompFlow Physics policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::PhysicsType type() noexcept {
          return ctr::PhysicsType::EULER;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for Physics policies.

    - Must define the function _viscousRhs()_, adding the viscous terms
      to the right hand side.

    - Must define the function _viscous_dt()_, computing the minumum time
      step size based on the viscous term.

    - Must define the function _conductRhs()_, adding the heat conduction
      terms to the right hand side.

    - Must define the function _conduct_dt()_, computing the minumum time
      step size based on the heat diffusion term.
*/
// *****************************************************************************
#ifndef CompFlowPhysicsCG_h
#define CompFlowPhysicsCG_h

#include <brigand/sequences/list.hpp>

#include "CGEuler.hpp"

namespace inciter {
namespace cg {

//! CompFlow Physics policies implemented using continuous Galerkin
using CompFlowPhysics = brigand::list< CompFlowPhysicsEuler >;

} // cg::
} // inciter::

#endif // CompFlowPhysicsCG_h
