// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     All problem configurations for the compressible flow equations
  \details   This file collects all Problem policy classes for the compressible
    flow equations, defined in PDE/MultiMat/MultiMat.h.

    General requirements on MuliMatCompFlow Problem policy classes:

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
#ifndef MultiMatProblem_h
#define MultiMatProblem_h

#include <brigand/sequences/list.hpp>

#include "Problem/UserDefined.hpp"
#include "Problem/InterfaceAdvection.hpp"
#include "Problem/SodShocktube.hpp"
#include "Problem/WaterAirShocktube.hpp"
#include "Problem/ShockHeBubble.hpp"
#include "Problem/UnderwaterEx.hpp"
#include "Problem/ShockDensityWave.hpp"
#include "Problem/EquilInterfaceAdvect.hpp"
#include "Problem/RichtmyerMeshkov.hpp"
#include "Problem/SinewavePacket.hpp"
#include "Problem/MixedCell.hpp"

namespace inciter {

//! List of all MultiMat Problem policies (defined in the includes above)
using MultiMatProblems =
  brigand::list< MultiMatProblemUserDefined
               , MultiMatProblemSodShocktube
               , MultiMatProblemInterfaceAdvection
               , MultiMatProblemWaterAirShocktube
               , MultiMatProblemShockHeBubble
               , MultiMatProblemUnderwaterEx
               , MultiMatProblemShockDensityWave
               , MultiMatProblemEquilInterfaceAdvect
               , MultiMatProblemRichtmyerMeshkov
               , MultiMatProblemSinewavePacket
               , MultiMatProblemMixedCell >;

} // inciter::

#endif // MultiMatProblem_h
