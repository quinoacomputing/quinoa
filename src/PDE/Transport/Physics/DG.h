// *****************************************************************************
/*!
  \file      src/PDE/Transport/Physics/DG.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Physics configurations for scalar transport using discontinuous
             Galerkin
  \details   This file configures all Physics policy classes for the scalar
    transport equations implemented using discontinuous Galerkin discretization,
    defined in PDE/Transport/DGTransport.h.

    General requirements on DGTransport Physics policy classes:

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
#ifndef TransportDGPhysics_h
#define TransportDGPhysics_h

#include <boost/mpl/vector.hpp>

#include "DGAdvection.h"

namespace inciter {
namespace dg {

//! Transport Physics policies implemented using discontinuous Galerkin
using TransportPhysics = boost::mpl::vector< TransportPhysicsAdvection >;

} // dg::
} // inciter::

#endif // TransportDGPhysics_h
