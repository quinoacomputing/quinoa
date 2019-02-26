// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Physics/DGVelEq.h
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Physics policy for the Euler equation governing multi-material flow
    using a finite volume method
  \details   This file defines a Physics policy class for the compressible
    flow equations class dg::MultiMat, defined in PDE/MultiMat/DGMultiMat.h.
    This specific algorithm assumes multi-material flow with a single velocity
    (velocity equilibirum) and uses a finite volume discretization scheme. See
    PDE/MultiMat/Physics/DG.h for general requirements on Physics policy classes
    for dg::MultiMat.
*/
// *****************************************************************************
#ifndef MultiMatPhysicsDGVelEq_h
#define MultiMatPhysicsDGVelEq_h

#include "Types.h"
#include "Exception.h"
#include "Inciter/Options/Physics.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! MultiMat system of PDEs problem: VelEq (velocity equilibrium)
//! \details This class is a no-op, consistent with no additional physics needed
//!   to make the basic implementation in MultiMat the Euler equations
//!   governing multi-material compressible flow.
class MultiMatPhysicsVelEq {

  public:
    //! Return enum denoting physics policy
    //! \return Enum denoting physics policy.
    static ctr::PhysicsType type() noexcept { return ctr::PhysicsType::VELEQ; }
};

} // dg::

} // inciter::

#endif // CompFlowPhysicsDGMultiMatVelEq_h
