// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
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

   - Must define the static function _init()_, used for initialization of the
      computed fields as well as sampling the analytical solution (if exist) at
      time t.

    - Must define the static function _side()_,  used to query all side set IDs
      the user has configured for all components.

    - Must define the static function _dirbc()_,  used to query Dirichlet
      boundary condition value on a given side set for all components in the PDE
      system.

    - Must define the static function _velocity()_, used to query the velocity
      field at cell nodes.

    - Must define the static function _prescribedVelocity()_, used to query the
      prescribed velocity at a point.
*/
// *****************************************************************************
#ifndef TransportProblem_h
#define TransportProblem_h

#include <boost/mpl/vector.hpp>

#include "Problem/ShearDiff.h"
#include "Problem/SlotCyl.h"

namespace inciter {

//! List of all Transport Problem policies (defined in the includes above)
using TransportProblems = boost::mpl::vector< TransportProblemShearDiff
                                            , TransportProblemSlotCyl
                                            >;

} // inciter::

#endif // TransportProblem_h
