// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/VorticalFlow.cpp
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

#include "VorticalFlow.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemVorticalFlow;

tk::InitializeFn::result_type
CompFlowProblemVorticalFlow::initialize( ncomp_t system,
                                         ncomp_t,
                                         tk::real x,
                                         tk::real y,
                                         tk::real z,
                                         tk::real )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  using tag::param; using tag::compflow;

  // manufactured solution parameters
  auto a = g_inputdeck.get< param, compflow, tag::alpha >()[ system ];
  auto b = g_inputdeck.get< param, compflow, tag::beta >()[ system ];
  auto p0 = g_inputdeck.get< param, compflow, tag::p0 >()[ system ];
  // ratio of specific heats
  const auto& matprop =
    g_inputdeck.get< tag::param, tag::compflow, tag::material >()[system];
  const auto& meos = g_inputdeck.get< tag::param, tag::compflow,
    tag::matidxmap >().get< tag::eosidx >()[0];
  auto g = matprop[meos].get< tag::gamma >()[0];
  // velocity
  auto ru = a*x - b*y;
  auto rv = b*x + a*y;
  auto rw = -2.0*a*z;
  // total specific energy
  auto rE = (ru*ru+rv*rv+rw*rw)/2.0 + (p0-2.0*a*a*z*z)/(g-1.0);

  return {{ 1.0, ru, rv, rw, rE }};
}

tk::InitializeFn::result_type
CompFlowProblemVorticalFlow::analyticSolution( ncomp_t system,
                                               ncomp_t,
                                               tk::real x,
                                               tk::real y,
                                               tk::real z,
                                               tk::real )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  using tag::param; using tag::compflow;

  // manufactured solution parameters
  auto a = g_inputdeck.get< param, compflow, tag::alpha >()[ system ];
  auto b = g_inputdeck.get< param, compflow, tag::beta >()[ system ];
  auto p0 = g_inputdeck.get< param, compflow, tag::p0 >()[ system ];
  // ratio of specific heats
  const auto& matprop =
    g_inputdeck.get< tag::param, tag::compflow, tag::material >()[system];
  const auto& meos = g_inputdeck.get< tag::param, tag::compflow,
    tag::matidxmap >().get< tag::eosidx >()[0];
  auto g = matprop[meos].get< tag::gamma >()[0];
  // velocity
  auto ru = a*x - b*y;
  auto rv = b*x + a*y;
  auto rw = -2.0*a*z;
  // total specific energy
  auto rE = (ru*ru+rv*rv+rw*rw)/2.0 + (p0-2.0*a*a*z*z)/(g-1.0);
  // pressure
  auto p = p0 - 2.0*a*a*z*z;

  return {{ 1.0, ru, rv, rw, rE, p }};
}

std::vector< std::string >
CompFlowProblemVorticalFlow::analyticFieldNames( ncomp_t ) const
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
CompFlowProblemVorticalFlow::names( ncomp_t ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
