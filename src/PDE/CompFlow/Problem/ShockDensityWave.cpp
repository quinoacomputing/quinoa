// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/ShockDensityWave.cpp
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

#include "ShockDensityWave.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemShockDensityWave;

tk::InitializeFn::result_type
CompFlowProblemShockDensityWave::initialize( ncomp_t,
                                             const std::vector< EOS >& mat_blk,
                                             tk::real x,
                                             tk::real,
                                             tk::real,
                                             tk::real )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] x X coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
//! \details This function only initializes the Shock-density wave problem, but
//!   does not actually give the analytical solution at time greater than 0.
//!   This problem does not have an analytical solution.
// *****************************************************************************
{
  using tag::param;

  tk::real r, p, u, v, w, rE;
  if (x > -4.0) {
    // density
    r = 1.0 + 0.2 * sin(5.0 * x);
    // pressure
    p = 1.0;
    // velocity
    u = 0.0;
    v = 0.0;
    w = 0.0;
  }
  else {
    // density
    r = 3.8571;
    // pressure
    p = 10.3333;
    // velocity
    u = 2.6294;
    v = 0.0;
    w = 0.0;
  }
  // total specific energy
  rE = mat_blk[0].compute< EOS::totalenergy >( r, u, v, w, p );

  return {{ r, r*u, r*v, r*w, rE }};
}

tk::InitializeFn::result_type
CompFlowProblemShockDensityWave::analyticSolution(
 ncomp_t,
 const std::vector< EOS >& mat_blk,
 tk::real x,
 tk::real,
 tk::real,
 tk::real )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] x X coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
//! \warning This is NOT the analytic solution at all times, only at t=0
// *****************************************************************************
{
  return initialize( 0, mat_blk, x, 0, 0, 0 );
}

std::vector< std::string >
CompFlowProblemShockDensityWave::analyticFieldNames( ncomp_t ) const
// *****************************************************************************
// Return analytic field names to be output to file
//! \return Vector of strings labelling analytic fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;
  n.push_back( "density_analytical" );
  n.push_back( "x-velocity_analytical" );
  n.push_back( "y-velocity_analytical" );
  n.push_back( "z-velocity_analytical" );
  n.push_back( "specific_total_energy_analytical" );
  n.push_back( "pressure_analytical" );

  return n;
}

std::vector< std::string >
CompFlowProblemShockDensityWave::names( ncomp_t ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
