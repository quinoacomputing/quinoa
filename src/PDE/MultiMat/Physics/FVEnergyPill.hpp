// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Physics/FVEnergyPill.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Physics policy for the Euler equation governing multi-material flow
    using a finite volume method
  \details   This file defines a Physics policy class for the compressible
    flow equations class fv::MultiMat, defined in PDE/MultiMat/FVMultiMat.h.
    This specific class allows energy pill initialization of a user defined
    box/block for multi-material flow. See PDE/MultiMat/Physics/FV.h for general
    requirements on Physics policy classes for fv::MultiMat.
*/
// *****************************************************************************
#ifndef MultiMatPhysicsFVEnergyPill_h
#define MultiMatPhysicsFVEnergyPill_h

#include "Inciter/Options/Physics.hpp"
#include "Fields.hpp"

namespace inciter {

namespace fv {

//! MultiMat system of PDEs problem: EnergyPill (velocity equilibrium)
//! \details This class is a no-op, consistent with no additional physics needed
//!   to make the basic implementation in MultiMat the Euler equations
//!   governing multi-material compressible flow.
class MultiMatPhysicsEnergyPill {

  public:
    //! Compute the time step size restriction based on this physics
    tk::real dtRestriction( std::size_t system,
      const tk::Fields& geoElem,
      std::size_t nelem,
      const int engSrcAd ) const;

    //! Return enum denoting physics policy
    //! \return Enum denoting physics policy.
    static ctr::PhysicsType type() noexcept {
      return ctr::PhysicsType::ENERGYPILL; }
};

} // fv::

} // inciter::

#endif // CompFlowPhysicsFVMultiMatEnergyPill_h
