// *****************************************************************************
/*!
  \file      src/PDE/Transport/Physics/CG.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Physics configurations for scalar transport using continuous
             Galerkin
  \details   This file configures all Physics policy classes for the scalar
    transport equations implemented using continuous Galerkin discretization,
    defined in PDE/Transport/CGTransport.h.

    General requirements on CGTransport Physics policy classes:

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
#ifndef TransportCGPhysics_h
#define TransportCGPhysics_h

#include <boost/mpl/vector.hpp>

#include "CGAdvection.h"
#include "CGAdvDiff.h"

namespace inciter {
namespace cg {

//! Transport Physics policies implemented using continuous Galerkin
using TransportPhysics = boost::mpl::vector< TransportPhysicsAdvection
                                           , TransportPhysicsAdvDiff >;

} // cg::
} // inciter::

#endif // TransportPhysics_h
