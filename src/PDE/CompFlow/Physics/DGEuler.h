// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Physics/DGEuler.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Physics policy for the Euler equation governing single-material
    flow using a continuous Galerkin finite element method
  \details   This file defines a Physics policy class for the compressible
    flow equations class dg::CompFlow, defined in PDE/CompFlow/DGCompFlow.h.
    This specific algorithm assumes single-material flow and uses a
    discontinuous Galerkin finite element discretization scheme. See
    PDE/CompFlow/Physics/DG.h for general requirements on Physics policy classes
    for dg::CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowPhysicsDGEuler_h
#define CompFlowPhysicsDGEuler_h

#include "Inciter/Options/Physics.h"

namespace inciter {
namespace dg {

//! CompFlow system of PDEs problem: Euler (inviscid flow)
//! \details This class is a no-op, consistent with no additional physics needed
//!   to make the basic implementation in CompFlow the Euler equations governing
//!   compressible flow.
class CompFlowPhysicsEuler {

  private:
    using ncomp_t = tk::ctr::ncomp_type;

  public:
    //! Return physics type
    //! \return Physics type
    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::EULER; }
};

} // dg::
} // inciter::

#endif // CompFlowPhysicsDGEuler_h
