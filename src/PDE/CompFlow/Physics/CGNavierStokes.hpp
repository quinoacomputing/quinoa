// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Physics/CGNavierStokes.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Physics policy for the Navier-Stokes equation using continuous
    Galerkin
  \details   This file declares a Physics policy class for the compressible
    single-material viscous flow equations using continuous Galerkin
    discretization, defined in PDE/CompFlow/CGCompFlow.h. The class defined here
    is used to configure the behavior of CGCompFlow. See
    PDE/CompFlow/Physics/CG.h for general requirements on Physics policy classes
    for CGCompFlow.
*/
// *****************************************************************************
#ifndef CompFlowPhysicsCGNavierStokes_h
#define CompFlowPhysicsCGNavierStokes_h

#include <array>
#include <limits>

#include "Fields.hpp"
#include "Inciter/Options/Physics.hpp"

namespace inciter {

namespace cg {

//! CompFlow system of PDEs problem: Navier-Stokes (viscous flow)
//! \details This class adds the viscous force contributions to the momentum and
//!    energy conservation equations governing compressible flow.
class CompFlowPhysicsNavierStokes {

  public:
    //! Add viscous stress contribution to momentum and energy rhs
    void
    viscousRhs( tk::real dt,
                tk::real J,
                const std::array< std::size_t, 4 >& N,
                const std::array< std::array< tk::real, 3 >, 4 >& grad,
                const std::array< std::array< tk::real, 4 >, 5 >& u,
                const std::array< const tk::real*, 5 >& r,
                tk::Fields& R ) const;

    //! Compute the minimum time step size based on the viscous force
    tk::real
    viscous_dt( tk::real L,
                const std::array< std::array< tk::real, 4 >, 5 >& u ) const;

    //! Add heat conduction contribution to the energy rhs
    void
    conductRhs( tk::real dt,
                tk::real J,
                const std::array< std::size_t, 4 >& N,
                const std::array< std::array< tk::real, 3 >, 4 >& grad,
                const std::array< std::array< tk::real, 4 >, 5 >& u,
                const std::array< const tk::real*, 5 >& r,
                tk::Fields& R ) const;

    //! Compute the minimum time step size based on thermal diffusion
    tk::real
    conduct_dt( tk::real L,
                tk::real g,
                const std::array< std::array< tk::real, 4 >, 5 >& u ) const;

    //! Return phsyics type
    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::NAVIERSTOKES; }
};

} // cg::

} // inciter::

#endif // CompFlowPhysicsCGNavierStokes_h
