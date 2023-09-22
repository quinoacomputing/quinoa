// *****************************************************************************
/*!
  \file      src/PDE/Transport/Physics/CGAdvDiff.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Physics policy for advection-diffusion using continuous Galerkin
  \details   This file declares a Physics policy class for the transport
    equations, defined in PDE/Transport/CGTransport.h implementing
    node-centered continuous Galerkin (CG) discretizations.
    See PDE/Transport/Physics/CG.h for general requirements on Physics policy
    classes for cg::Transport.
*/
// *****************************************************************************
#ifndef TransportPhysicsCGAdvDiff_h
#define TransportPhysicsCGAdvDiff_h

#include "Types.hpp"
#include "Fields.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/Physics.hpp"

namespace inciter {
namespace cg {

//! Physics policy for advection-diffusion using continuous Galerkin
class TransportPhysicsAdvDiff {

  private:
    using ncomp_t = tk::ctr::ncomp_t;
    using eq = tag::transport;

  public:
    //! Add diffusion contribution to rhs
    void
    diffusionRhs( ncomp_t ncomp,
                  tk::real deltat,
                  tk::real J,
                  const std::array< std::array< tk::real, 3 >, 4 >& grad,
                  const std::array< std::size_t, 4 >& N,
                  const std::vector< std::array< tk::real, 4 > >& u,
                  const std::vector< const tk::real* >& r,
                  tk::Fields& R ) const;

    //! Compute the minimum time step size based on the diffusion
    tk::real
    diffusion_dt( ncomp_t ncomp,
                  tk::real L,
                  const std::vector< std::array< tk::real, 4 > >& ) const;

    //! Return physics type
    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::ADVDIFF; }
};

} // cg::
} // inciter::

#endif // TransportPhysicsCGAdvDiff_h
