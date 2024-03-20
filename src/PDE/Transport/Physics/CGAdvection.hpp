// *****************************************************************************
/*!
  \file      src/PDE/Transport/Physics/CGAdvection.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Physics policy for advection using continuous Galerkin (no-op)
  \details   This file defines a Physics policy class for scalar transport using
     continuous Galerkin discretization, defined in PDE/Transport/CGTransport.h.
     The class defined here is used to configure the behavior of CGTransport.
     See PDE/Transport/Physics/CG.h for general requirements on Physics policy
     classes for CGTransport.
*/
// *****************************************************************************
#ifndef TransportPhysicsCGAdvection_h
#define TransportPhysicsCGAdvection_h

#include <limits>

#include "Inciter/Options/Physics.hpp"

namespace inciter {
namespace cg {

//! Physics policy for advection using continuous Galerkin (no-op)
//! \details This class is a no-op, consistent with no additional physics needed
//!   to make the basic implementation in Transport the advection equation.
class TransportPhysicsAdvection {

  private:
    using ncomp_t = tk::ctr::ncomp_t;

  public:
    //! Add diffusion contribution to rhs (no-op for advection only)
    void
    diffusionRhs( ncomp_t,
                  tk::real,
                  tk::real,
                  const std::array< std::array< tk::real, 3 >, 4 >&,
                  const std::array< std::size_t, 4 >&,
                  const std::vector< std::array< tk::real, 4 > >&,
                  const std::vector< const tk::real* >&,
                  tk::Fields& ) const {}

    //! Compute the minimum time step size based on the diffusion
    //! \return A large time step size, i.e., ignore
    tk::real
    diffusion_dt( tk::ctr::ncomp_t,
                  tk::real,
                  const std::vector< std::array< tk::real, 4 > >& ) const
    { return std::numeric_limits< tk::real >::max(); }

    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::ADVECTION; }
};

} // cg::
} // inciter::

#endif // TransportPhysicsCGAdvection_h
