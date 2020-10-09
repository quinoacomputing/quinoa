// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/SodShocktube.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a Problem policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problem.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************

#include "SodShocktube.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemSodShocktube;

tk::SolutionFn::result_type
CompFlowProblemSodShocktube::solution( ncomp_t system,
                                       [[maybe_unused]] ncomp_t ncomp,
                                       tk::real x,
                                       tk::real,
                                       tk::real,
                                       tk::real,
                                       int& )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::SolutionFn
//! \details This function only initializes the Sod shock tube problem, but does
//!   not actually give the analytical solution at time greater than 0. The
//!   analytical solution would require an exact Riemann solver, which has not
//!   been implemented yet.
// *****************************************************************************
{
  Assert( ncomp == ncomp, "Number of scalar components must be " +
                          std::to_string(ncomp) );
  using tag::param;

  tk::real r, p, u, v, w, rE;
  if (x<0.5) {
    // density
    r = 1.0;
    // pressure
    p = 1.0;
    // velocity
    u = 0.0;
    v = 0.0;
    w = 0.0;
  }
  else {
    // density
    r = 0.125;
    // pressure
    p = 0.1;
    // velocity
    u = 0.0;
    v = 0.0;
    w = 0.0;
  }
  // total specific energy
  rE = eos_totalenergy< eq >( system, r, u, v, w, p );

  return {{ r, r*u, r*v, r*w, rE }};
}

std::vector< std::string >
CompFlowProblemSodShocktube::fieldNames( ncomp_t ) const
// *****************************************************************************
// Return field names to be output to file
//! \return Vector of strings labelling fields output in file
// *****************************************************************************
{
  auto n = CompFlowFieldNames();

  const auto pref = inciter::g_inputdeck.get< tag::pref, tag::pref >();
  if (pref) n.push_back( "number of degrees of freedom" );

  return n;
}

std::vector< std::vector< tk::real > >
CompFlowProblemSodShocktube::fieldOutput(
  ncomp_t system,
  ncomp_t,
  ncomp_t offset,
  std::size_t nunk,
  tk::real,
  tk::real,
  const std::vector< tk::real >&,
  const std::array< std::vector< tk::real >, 3 >&,
  const tk::Fields& U ) const
// *****************************************************************************
//  Return field output going to file
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] offset System offset specifying the position of the system of
//!   PDEs among other systems
//! \param[in] nunk Number of unknowns to extract
//! \param[in] U Solution vector at recent time step
//! \return Vector of vectors to be output to file
// *****************************************************************************
{
  return CompFlowFieldOutput( system, offset, nunk, U );
}

std::vector< std::string >
CompFlowProblemSodShocktube::names( ncomp_t ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
