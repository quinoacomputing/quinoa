// *****************************************************************************
/*!
  \file      src/PDE/Transport/Physics/DGAdvection.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Physics policy for advection using continuous Galerkin (no-op)
  \details   This file defines a Physics policy class for scalar transport using
     discontinuous Galerkin discretization, defined in
     PDE/Transport/DGTransport.h. The class defined here is used to configure
     the behavior of DGTransport. See PDE/Transport/Physics/DG.h for general
     requirements on Physics policy classes for DGTransport.
*/
// *****************************************************************************
#ifndef TransportPhysicsDGAdvection_h
#define TransportPhysicsDGAdvection_h

#include <limits>

#include "Types.h"
#include "Inciter/Options/Physics.h"

namespace inciter {
namespace dg {

//! Physics policy for advection using continuous Galerkin (no-op)
//! \details This class is a no-op, consistent with no additional physics needed
//!   to make the basic implementation in Transport the advection equation.
class TransportPhysicsAdvection {

  public:

    //! Add diffusion contribution to rhs at 2nd step stage (no-op)
    static void
    diffusionRhs( tk::ctr::ncomp_type,
                  tk::ctr::ncomp_type,
                  tk::real,
                  tk::real,
                  const std::array< std::array< tk::real, 3 >, 4 >&,
                  const std::array< std::size_t, 4 >&,
                  const std::vector< std::array< tk::real, 4 > >&,
                  const std::vector< const tk::real* >&,
                  tk::Fields& )
    {}

    //! Compute the minimum time step size based on the diffusion
    //! \return A large time step size, i.e., ignore
    static tk::real
    diffusion_dt( tk::ctr::ncomp_type,
                  tk::ctr::ncomp_type,
                  tk::real,
                  const std::vector< std::array< tk::real, 4 > >& )
    { return std::numeric_limits< tk::real >::max(); }

    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::ADVECTION; }
};

} // dg::
} // inciter::

#endif // TransportPhysicsDGAdvection_h
