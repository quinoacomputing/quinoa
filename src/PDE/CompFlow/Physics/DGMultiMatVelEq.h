// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Physics/DGMultiMatVelEq.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Physics policy for the Euler equation governing multi-material flow
    using a finite volume method
  \details   This file defines a Physics policy class for the compressible
    flow equations class dg::CompFlow, defined in PDE/CompFlow/DGCompFlow.h.
    This specific algorithm assumes multi-material flow with a single velocity
    (velocity equilibirum) and uses a finite volume discretization scheme. See
    PDE/CompFlow/Physics/DG.h for general requirements on Physics policy classes
    for dg::CompFlow.
*/
// *****************************************************************************
#ifndef CompFlowPhysicsDGMultiMatVelEq_h
#define CompFlowPhysicsDGMultiMatVelEq_h

#include "Types.h"
#include "Exception.h"
#include "Inciter/Options/Physics.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! CompFlow system of PDEs problem: Euler (inviscid flow)
//! \details This class is a no-op, consistent with no additional physics needed
//!   to make the basic implementation in CompFlow the Euler equations governing
//!   compressible flow.
class CompFlowPhysicsMultiMatVelEq {

  private:
    using ncomp_t = tk::ctr::ncomp_type;

  public:
    //! ...
    static ctr::PhysicsType type() noexcept
    { return ctr::PhysicsType::MULTIMAT_VELEQ; }
};

} // dg::

} // inciter::

#endif // CompFlowPhysicsDGMultiMatVelEq_h
