// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Physics/DGEuler.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Physics policy for the Euler equation governing multi-material flow
    using a DG method
  \details   This file defines a Physics policy class for the compressible
    flow equations class dg::MultiMat, defined in PDE/MultiMat/DGMultiMat.h.
    This specific algorithm solves the Euler (inviscid) equations of
    multi-material flow and uses a finite volume discretization scheme. See
    PDE/MultiMat/Physics/DG.h for general requirements on Physics policy classes
    for dg::MultiMat.
*/
// *****************************************************************************
#ifndef MultiMatPhysicsDGEuler_h
#define MultiMatPhysicsDGEuler_h

#include "Types.hpp"
#include "Exception.hpp"
#include "Inciter/Options/Physics.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! MultiMat system of PDEs problem: Euler (inviscid)
//! \details This class is a no-op, consistent with no additional physics needed
//!   to make the basic implementation in MultiMat the Euler equations
//!   governing multi-material compressible flow.
class MultiMatPhysicsEuler {

  public:
    //! Compute the time step size restriction based on this physics
    //! \return A large time step size, i.e., ignore
    tk::real dtRestriction( std::size_t,
      const tk::Fields&,
      std::size_t,
      const int ) const
    { return std::numeric_limits< tk::real >::max(); }

    //! Compute sources corresponding to this physics
    void physSrc( std::size_t,
      std::size_t,
      tk::real,
      const tk::Fields&,
      const std::unordered_map< std::size_t, std::set< std::size_t > >&,
      tk::Fields&,
      int& ) const {}

    //! Return enum denoting physics policy
    //! \return Enum denoting physics policy.
    static ctr::PhysicsType type() noexcept { return ctr::PhysicsType::EULER; }
};

} // dg::

} // inciter::

#endif // CompFlowPhysicsDGMultiMatEuler_h
