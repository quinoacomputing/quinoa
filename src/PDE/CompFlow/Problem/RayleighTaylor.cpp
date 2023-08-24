// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/RayleighTaylor.cpp
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

#include "RayleighTaylor.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemRayleighTaylor;

tk::InitializeFn::result_type
CompFlowProblemRayleighTaylor::initialize( ncomp_t,
                                           const std::vector< EOS >& mat_blk,
                                           tk::real x,
                                           tk::real y,
                                           tk::real z,
                                           tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  using tag::param; using std::sin; using std::cos;

  // manufactured solution parameters
  const auto a = g_inputdeck.get< param, eq, tag::alpha >()[0];
  const auto bx = g_inputdeck.get< param, eq, tag::betax >()[0];
  const auto by = g_inputdeck.get< param, eq, tag::betay >()[0];
  const auto bz = g_inputdeck.get< param, eq, tag::betaz >()[0];
  const auto p0 = g_inputdeck.get< param, eq, tag::p0 >()[0];
  const auto r0 = g_inputdeck.get< param, eq, tag::r0 >()[0];
  const auto k = g_inputdeck.get< param, eq, tag::kappa >()[0];
  // spatial component of density and pressure fields
  const tk::real gx = bx*x*x + by*y*y + bz*z*z;
  // density
  const tk::real r = r0 - gx;
  // pressure
  const tk::real p = p0 + a*gx;
  // velocity
  const tk::real ft = cos(k*M_PI*t);
  const tk::real u = ft*z*sin(M_PI*x);
  const tk::real v = ft*z*cos(M_PI*y);
  const tk::real w = ft*(-0.5*M_PI*z*z*(cos(M_PI*x)-sin(M_PI*y)));
  // total specific energy
  const tk::real rE = mat_blk[0].compute< EOS::totalenergy >( r, u, v, w, p );

  return {{ r, r*u, r*v, r*w, rE }};
}

tk::InitializeFn::result_type
CompFlowProblemRayleighTaylor::analyticSolution(
  ncomp_t,
  const std::vector< EOS >& mat_blk,
  tk::real x,
  tk::real y,
  tk::real z,
  tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  using tag::param; using std::sin; using std::cos;

  // manufactured solution parameters
  auto a = g_inputdeck.get< param, eq, tag::alpha >()[0];
  auto bx = g_inputdeck.get< param, eq, tag::betax >()[0];
  auto by = g_inputdeck.get< param, eq, tag::betay >()[0];
  auto bz = g_inputdeck.get< param, eq, tag::betaz >()[0];
  auto p0 = g_inputdeck.get< param, eq, tag::p0 >()[0];
  auto r0 = g_inputdeck.get< param, eq, tag::r0 >()[0];
  auto k = g_inputdeck.get< param, eq, tag::kappa >()[0];
  // spatial component of density and pressure fields
  auto gx = bx*x*x + by*y*y + bz*z*z;
  // density
  auto r = r0 - gx;
  // pressure
  auto p = p0 + a*gx;
  // velocity
  auto ft = cos(k*M_PI*t);
  auto u = ft*z*sin(M_PI*x);
  auto v = ft*z*cos(M_PI*y);
  auto w = ft*(-0.5*M_PI*z*z*(cos(M_PI*x)-sin(M_PI*y)));
  // total specific energy
  auto E = mat_blk[0].compute< EOS::totalenergy >( r, u, v, w, p ) / r;

  return {{ r, u, v, w, E, p }};
}

std::vector< std::string >
CompFlowProblemRayleighTaylor::analyticFieldNames( ncomp_t ) const
// *****************************************************************************
// Return analytic field names to be output to file
//! \return Vector of strings labelling fields output in file
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
CompFlowProblemRayleighTaylor::names( ncomp_t /*ncomp*/ ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
