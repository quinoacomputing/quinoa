//******************************************************************************
/*!
  \file      src/PDE/AdvDiffProblem.h
  \author    J. Bakosi
  \date      Wed 03 Feb 2016 04:17:50 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Problem configurations for the advection-diffusion equation
  \details   This file defines policy classes for the advection-diffusion
    partial differential equation, defined in PDE/AdvDiff.h.

    General requirements on advection-diffusion partial differential equation
    problem policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::ProblemType type() noexcept {
          return ctr::ProblemType::SHEAR_DIFF;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
//******************************************************************************
#ifndef AdvDiffProblem_h
#define AdvDiffProblem_h

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! Advection-diffusion PDE problem: user-defined
class AdvDiffProblemUserDefined {
  public:
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::USER_DEFINED; }
};

//! Advection-diffusion PDE problem: diffusion of a shear layer
class AdvDiffProblemShearDiff {
  public:
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::SHEAR_DIFF; }
};

//! Advection-diffusion PDE problem: Zalesak's slotted cylinder
class AdvDiffProblemSlotCyl {
  public:
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::SLOT_CYL; }
};


//! List of all advection-diffusion PDE's problem policies
using AdvDiffProblems = boost::mpl::vector< AdvDiffProblemUserDefined
                                          , AdvDiffProblemShearDiff
                                          , AdvDiffProblemSlotCyl 
                                          >;

} // inciter::

#endif // AdvDiffProblem_h
