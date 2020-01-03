// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/RiemannFactory.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register available Riemann solvers for multimaterial compressible
             hydrodynamics into a factory
  \details   Register available Riemann solvers for multimaterial compressible
             hydrodynamics into a factory.
*/
// *****************************************************************************

#include <brigand/sequences/list.hpp>
#include <brigand/algorithms/for_each.hpp>

#include "RiemannFactory.hpp"
#include "Riemann/HLL.hpp"
#include "Riemann/AUSM.hpp"

inciter::multimatRiemannFactory
inciter::multimatRiemannSolvers()
// *****************************************************************************
// \brief Register available Riemann solvers for multimaterial compressible
//    hydrodynamics into a factory
//! \return Riemann solver factory
// *****************************************************************************
{
  using RiemannSolverList = brigand::list< AUSM, HLL >;
  multimatRiemannFactory r;
  brigand::for_each< RiemannSolverList >( registerRiemannSolver( r ) );
  return r;
}
