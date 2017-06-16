// *****************************************************************************
/*!
  \file      src/PDE/CompFlowProblem.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem configurations for the compressible flow equations
  \details   This file includes policy classes for the compressible flow
    equations, defined in PDE/CompFlow.h.

    General requirements on flow equations problem policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::ProblemType type() noexcept {
          return ctr::ProblemType::USER_DEFINED;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for problem policies.

    - Must define the static function _names()_, returning the names of integral
      variables to be output to diagnostics file.

    - Must define the static function _init()_, used for initialization of the
      computed fields as well as sampling the analytical solution (if exist) at
      time t.

    - Must define the static function _sourceRhs()_, used for adding source
      terms to the righ hand side.

    - Must define the static function _side()_,  used to uery all side set IDs
      the user has configured for all components.

    - Must define the static function _dirbc()_,  used to query Dirichlet
      boundary condition value on a given side set for all components in the PDE
      system.

    - Must define the static function _fieldNames()_, used to provide the field
      names to be output to file.

    - Must define the static function _fieldOutput()_, used to provide the field
      output.
*/
// *****************************************************************************
#ifndef CompFlowProblem_h
#define CompFlowProblem_h

#include <boost/mpl/vector.hpp>

#include "CompFlowProblem/UserDefined.h"
#include "CompFlowProblem/VorticalFlow.h"
#include "CompFlowProblem/NLEnergyGrowth.h"
#include "CompFlowProblem/RayleighTaylor.h"
#include "CompFlowProblem/TaylorGreen.h"

namespace inciter {

//! List of all CompFlow problem policies (defined in the includes above)
using CompFlowProblems = boost::mpl::vector< CompFlowProblemUserDefined
                                           , CompFlowProblemVorticalFlow
                                           , CompFlowProblemNLEnergyGrowth
                                           , CompFlowProblemRayleighTaylor
                                           , CompFlowProblemTaylorGreen >;

} // inciter::

#endif // CompFlowProblem_h
