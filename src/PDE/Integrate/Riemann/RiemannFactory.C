// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Riemann/RiemannFactory.C
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register available Riemann solvers into a factory
  \details   Register available Riemann solvers into a factory.
*/
// *****************************************************************************

#include <brigand/sequences/list.hpp>
#include <brigand/algorithms/for_each.hpp>

#include "RiemannFactory.h"
#include "HLLC.h"
#include "LaxFriedrichs.h"

inciter::RiemannFactory
inciter::RiemannSolvers()
// *****************************************************************************
// Register available Riemann solvers into a factory
//! \return Riemann solver factory
// *****************************************************************************
{
  using RiemannSolverList = brigand::list< HLLC, LaxFriedrichs >;
  RiemannFactory r;
  brigand::for_each< RiemannSolverList >( registerRiemannSolver( r ) );
  return r;
}
