// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Physics/DG.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Physics configurations for compressible single-material flow using
    discontinuous Galerkin
  \details   This file configures all Physics policy classes for compressible
    single-material flow implementied using discontinuous Galerkin
    discretization, defined in PDE/CompFlow/DGCompFlow.h.

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
#ifndef CompFlowPhysicsDG_h
#define CompFlowPhysicsDG_h

#include <boost/mpl/vector.hpp>

#include "DGEuler.h"
#include "DGNavierStokes.h"

namespace inciter {
namespace dg {

//! CompFlow Physics policies implemented using discontinuous Galerkin
using CompFlowPhysics = boost::mpl::vector< CompFlowPhysicsEuler
                                          , CompFlowPhysicsNavierStokes >;

} // dg::
} // inciter::

#endif // CompFlowPhysicsDG_h
