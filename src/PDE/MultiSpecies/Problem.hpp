// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Problem.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     All problem configurations for the compressible flow equations
  \details   This file collects all Problem policy classes for the compressible
    flow equations, defined in PDE/MultiSpecies/MultiSpecies.h.

    General requirements on MultiSpecies Problem policy classes:

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

    - Must define the static function _initialize()_, used for initialization of
      the computed fields at time t.

    - Must define the static function _analyticSolution()_, used for
      sampling the analytical solution if exist) at time t.

    - Must define the static function _src()_, used for adding source terms to
      the righ hand side.

    - Must define the static function _dirbc()_,  used to query Dirichlet
      boundary condition value on a given side set for all components in the PDE
      system.
*/
// *****************************************************************************
#ifndef MultiSpeciesProblem_h
#define MultiSpeciesProblem_h

#include <brigand/sequences/list.hpp>

#include "Problem/UserDefined.hpp"

namespace inciter {

//! List of all MultiSpecies Problem policies (defined in the includes above)
using MultiSpeciesProblems =
  brigand::list< MultiSpeciesProblemUserDefined
               >;

} // inciter::

#endif // MultiSpeciesProblem_h
