// *****************************************************************************
/*!
  \file      src/PDE/MultiMatCompFlow/Physics/DG.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Physics configurations for multi-material compressible flow using
    discontinuous Galerkin finite element methods
  \details   This file configures all Physics policy classes for multi-material
    compressible flow implementied using discontinuous Galerkin finite element
    discretizations, defined in PDE/MultiMatCompFlow/DGMultiMatCompFlow.h.

    General requirements on MultiMatCompFlow Physics policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::PhysicsType type() noexcept {
          return ctr::PhysicsType::MULTIMAT_VELEQ;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for Physics policies.
*/
// *****************************************************************************
#ifndef MultiMatCompFlowPhysicsDG_h
#define MultiMatCompFlowPhysicsDG_h

#include <brigand/sequences/list.hpp>

#include "DGVelEq.h"

namespace inciter {
namespace dg {

//! MultiMatCompFlow Physics policies implemented using discontinuous Galerkin
using MultiMatCompFlowPhysics = brigand::list< MultiMatCompFlowPhysicsVelEq >;

} // dg::
} // inciter::

#endif // MultiMatCompFlowPhysicsDG_h
