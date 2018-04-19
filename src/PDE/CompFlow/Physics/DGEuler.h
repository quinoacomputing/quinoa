// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Physics/DGEuler.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Physics policy for the Euler equation using continuous Galerkin
  \details   This file defines a Physics policy class for the compressible
    single-material inviscid flow equations using continuous Galerkin
    discretization, defined in PDE/CompFlow/DGCompFlow.h. The class defined here
    is used to configure the behavior of DGCompFlow. See
    PDE/CompFlow/Physics/DG.h for general requirements on Physics policy classes
    for DGCompFlow.
*/
// *****************************************************************************
#ifndef CompFlowPhysicsDGEuler_h
#define CompFlowPhysicsDGEuler_h

#include <array>
#include <limits>

#include "Types.h"
#include "Inciter/Options/Physics.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! CompFlow system of PDEs problem: Euler (inviscid flow)
//! \details This class is a no-op, consistent with no additional physics needed
//!   to make the basic implementation in CompFlow the Euler equations governing
//!   compressible flow.
class CompFlowPhysicsEuler {

  public:
    //! Add viscous stress contribution to momentum and energy rhs (no-op)
    static void
    viscousRhs( tk::real,
                tk::real,
                const std::array< std::size_t, 4 >&,
                const std::array< std::array< tk::real, 3 >, 4 >&,
                const std::array< std::array< tk::real, 4 >, 5 >&,
                const std::array< const tk::real*, 5 >&,
                tk::Fields& ) {}

    //! Compute the minimum time step size based on the viscous force
    //! \return A large time step size, i.e., ignore
    static tk::real
    viscous_dt( tk::real, const std::array< std::array< tk::real, 4 >, 5 >& )
    { return std::numeric_limits< tk::real >::max(); }

    //! Add heat conduction contribution to energy rhs (no-op)
    static void
    conductRhs( tk::real,
                tk::real,
                const std::array< std::size_t, 4 >&,
                const std::array< std::array< tk::real, 3 >, 4 >&,
                const std::array< std::array< tk::real, 4 >, 5 >&,
                const std::array< const tk::real*, 5 >&,
                tk::Fields& ) {}

    //! Compute the minimum time step size based on thermal diffusion
    //! \return A large time step size, i.e., ignore
    static tk::real
    conduct_dt( tk::real,  const std::array< std::array< tk::real, 4 >, 5 >& )
    { return std::numeric_limits< tk::real >::max(); }

    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::EULER; }
};

} // dg::

} // inciter::

#endif // CompFlowPhysicsDGEuler_h
