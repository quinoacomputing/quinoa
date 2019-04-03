// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Physics/DG.hpppp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Physics configurations for compressible flow using discontinuous
    Galerkin finite element methods
  \details   This file configures all Physics policy classes for compressible
    flow implementied using discontinuous Galerkin finite element
    discretizations, defined in PDE/CompFlow/DGCompFlow.h.

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
*/
// *****************************************************************************
#ifndef CompFlowPhysicsDG_h
#define CompFlowPhysicsDG_h

#include <brigand/sequences/list.hpp>

#include "DGEuler.hpp"

namespace inciter {
namespace dg {

//! CompFlow Physics policies implemented using discontinuous Galerkin
using CompFlowPhysics = brigand::list< CompFlowPhysicsEuler >;

} // dg::
} // inciter::

#endif // CompFlowPhysicsDG_h
