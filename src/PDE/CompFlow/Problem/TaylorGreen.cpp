// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/TaylorGreen.cpp
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

#include "TaylorGreen.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemTaylorGreen;

tk::InitializeFn::result_type
CompFlowProblemTaylorGreen::initialize( ncomp_t system,
                                        ncomp_t,
                                        tk::real x,
                                        tk::real y,
                                        tk::real,
                                        tk::real )
// *****************************************************************************
//! Initialize numerical solution
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  using tag::param; using std::sin; using std::cos;

  // density
  auto r = 1.0;
  // pressure
  auto p = 10.0 + r/4.0*(cos(2.0*M_PI*x) + cos(2.0*M_PI*y));
  // velocity
  auto u =  sin(M_PI*x) * cos(M_PI*y);
  auto v = -cos(M_PI*x) * sin(M_PI*y);
  auto w = 0.0;
  // total specific energy
  auto rE = eos_totalenergy< eq >( system, r, u, v, w, p );

  return {{ r, r*u, r*v, r*w, rE }};
}

tk::InitializeFn::result_type
CompFlowProblemTaylorGreen::analyticSolution( ncomp_t system,
                                              ncomp_t,
                                              std::vector< EoS_Base* >,
                                              tk::real x,
                                              tk::real y,
                                              tk::real,
                                              tk::real )
// *****************************************************************************
//  Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  using tag::param; using std::sin; using std::cos;

  // density
  auto r = 1.0;
  // pressure
  auto p = 10.0 + r/4.0*(cos(2.0*M_PI*x) + cos(2.0*M_PI*y));
  // velocity
  auto u =  sin(M_PI*x) * cos(M_PI*y);
  auto v = -cos(M_PI*x) * sin(M_PI*y);
  auto w = 0.0;
  // total specific energy
  auto E = eos_totalenergy< eq >( system, r, u, v, w, p );

  return {{ r, u, v, w, E, p }};
}

std::vector< std::string >
CompFlowProblemTaylorGreen::analyticFieldNames( ncomp_t ) const
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
CompFlowProblemTaylorGreen::names( ncomp_t ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
