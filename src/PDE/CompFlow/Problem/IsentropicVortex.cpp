// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/IsentropicVortex.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a Problem policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problem.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************

#include "IsentropicVortex.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"
#include "EoS/GetMatProp.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemIsentropicVortex;

tk::InitializeFn::result_type
CompFlowProblemIsentropicVortex::initialize(
  [[maybe_unused]] ncomp_t ncomp,
  const std::vector< EOS >& mat_blk,
  tk::real x,
  tk::real y,
  tk::real z,
  tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  return analyticSolution(ncomp, mat_blk, x, y, z, t);
}

tk::InitializeFn::result_type
CompFlowProblemIsentropicVortex::analyticSolution(
  [[maybe_unused]] ncomp_t ncomp,
  const std::vector< EOS >& mat_blk,
  tk::real x,
  tk::real y,
  tk::real,
  tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  Assert( ncomp == 5, "Number of scalar components must be 5" );

  using tag::param;

  // velocity
  std::array< tk::real, 3 > vel{{1.0, 1.0, 0.0}};

  // center of the vortex
  std::array< tk::real, 3 > x0{{5.0, 5.0, 0.0}};
  for (std::size_t i=0; i<3; ++i)
    x0[i] += vel[i]*t;

  // distance from the center
  tk::real r2 = (x-x0[0])*(x-x0[0]) + (y-x0[1])*(y-x0[1]);

  tk::real r, p, rE;
  tk::real estr(5.0);
  tk::real gam = getmatprop< tag::compflow, tag::gamma >();

  // perturbed density
  r = std::pow( (1.0 - ((gam-1.0)*estr*estr)/(8.0*gam*M_PI*M_PI) * exp(1.0-r2)),
    (1.0/(gam-1.0)) );
  // perturbed velocity
  vel[0] -= 0.5*estr/M_PI * exp((1.0-r2)/0.5) * (y-x0[1]);
  vel[1] += 0.5*estr/M_PI * exp((1.0-r2)/0.5) * (x-x0[0]);
  // pressure
  p = 1.0;
  // total specific energy
  rE = mat_blk[0].compute< EOS::totalenergy >( r, vel[0], vel[1], vel[2], p );

  return {{ r, r*vel[0], r*vel[1], r*vel[2], rE }};
}

std::vector< std::string >
CompFlowProblemIsentropicVortex::analyticFieldNames( ncomp_t ) const
// *****************************************************************************
// Return analytic field names to be output to file
//! \return Vector of strings labelling fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;
  n.push_back( "density_analytical" );
  n.push_back( "x-momentum_analytical" );
  n.push_back( "y-momentum_analytical" );
  n.push_back( "z-momentum_analytical" );
  n.push_back( "total_energy_analytical" );

  return n;
}

std::vector< std::string >
CompFlowProblemIsentropicVortex::names( ncomp_t ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
