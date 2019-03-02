// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     All problem configurations for the scalar transport equations
  \details   This file collects all Problem policy classes for the scalar
    transport equations, defined in PDE/Transport/Transport.h.

    General requirements on Transport Problem policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::ProblemType type() noexcept {
          return ctr::ProblemType::SHEAR_DIFF;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.

    - Must define the static function _errchk()_, doing general error checks.

    - Must define the static function _solution()_, used to evaluate the
      analytic solution (if defined) and for initialization of the computed
      fields at time _t_.

    - Must define the static function _solinc()_, used to evaluate the
      increment from t to t+dt of the analytic solution (if defined).

    - Must define the static function _side()_,  used to query all side set IDs
      the user has configured for all components.
  
    - Must define the static function _dirbc()_,  used to query Dirichlet
      boundary condition value on a given side set for all components in the PDE
      system.
  
    - Must define the static function _prescribedVelocity()_, used to query the
      prescribed velocity at a point.
*/
// *****************************************************************************
#ifndef TransportProblem_h
#define TransportProblem_h

#include <brigand/sequences/list.hpp>

#include "Problem/ShearDiff.h"
#include "Problem/SlotCyl.h"
#include "Problem/GaussHump.h"
#include "Problem/CylAdvect.h"

namespace inciter {

//! List of all Transport Problem policies (defined in the includes above)
using TransportProblems = brigand::list< TransportProblemShearDiff
                                       , TransportProblemSlotCyl
                                       , TransportProblemGaussHump
                                       , TransportProblemCylAdvect
                                       >;

} // inciter::

#endif // TransportProblem_h
