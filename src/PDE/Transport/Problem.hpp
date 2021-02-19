// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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

    - Must define the function _errchk()_, doing general error checks.

    - Must define the static function _initialize()_, used for initialization of
      the computed fields at time t.

    - Must define the static function _analyticSolution()_, used for
      sampling the analytical solution if exist) at time t.

    - Must define the static function _prescribedVelocity()_, used to query the
      prescribed velocity at a point.
*/
// *****************************************************************************
#ifndef TransportProblem_h
#define TransportProblem_h

#include <brigand/sequences/list.hpp>

#include "Problem/ShearDiff.hpp"
#include "Problem/SlotCyl.hpp"
#include "Problem/GaussHump.hpp"
#include "Problem/CylAdvect.hpp"
#include "Problem/CylVortex.hpp"

namespace inciter {

//! List of all Transport Problem policies (defined in the includes above)
using TransportProblems = brigand::list< TransportProblemShearDiff
                                       , TransportProblemSlotCyl
                                       , TransportProblemGaussHump
                                       , TransportProblemCylAdvect
                                       , TransportProblemCylVortex
                                       >;

} // inciter::

#endif // TransportProblem_h
