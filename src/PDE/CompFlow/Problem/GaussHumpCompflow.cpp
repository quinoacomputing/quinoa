// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/GaussHumpCompflow.cpp
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

#include "GaussHumpCompflow.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemGaussHump;

tk::InitializeFn::result_type
CompFlowProblemGaussHump::initialize( ncomp_t system,
                                      ncomp_t ncomp,
                                      const std::vector< EOS >& mat_blk,
                                      tk::real x,
                                      tk::real y,
                                      tk::real,
                                      tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
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

  const auto vel = prescribedVelocity( system, ncomp, x, y, 0.0 );

  tk::real r, p, u, v, w, rE;

  // center of the hump
  auto x0 = 0.25 + vel[0][0]*t;
  auto y0 = 0.25 + vel[0][1]*t;

  // density
  r = 1.0 + exp( -((x-x0)*(x-x0)
                 + (y-y0)*(y-y0))/(2.0 * 0.005) );
  // pressure
  p = 1.0;
  // velocity
  u = 1;
  v = 1;
  w = 0;
  // total specific energy
  rE = mat_blk[0].eosCall< EOS::totalenergy >( r, u, v, w, p );

  return {{ r, r*u, r*v, r*w, rE }};
}

tk::InitializeFn::result_type
CompFlowProblemGaussHump::analyticSolution( ncomp_t system,
                                            ncomp_t ncomp,
                                            const std::vector< EOS >& mat_blk,
                                            tk::real x,
                                            tk::real y,
                                            tk::real,
                                            tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
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

  const auto vel = prescribedVelocity( system, ncomp, x, y, 0.0 );

  // center of the hump
  auto x0 = 0.25 + vel[0][0]*t;
  auto y0 = 0.25 + vel[0][1]*t;

  // density
  auto r = 1.0 + exp( -((x-x0)*(x-x0) + (y-y0)*(y-y0))/(2.0 * 0.005) );
  // pressure
  auto p = 1.0;
  // velocity
  auto u = 1.0;
  auto v = 1.0;
  auto w = 0.0;
  // total specific energy
  auto E = mat_blk[0].eosCall< EOS::totalenergy >( r, u, v, w, p ) / r;

  return {{ r, u, v, w, E, p }};
}

std::vector< std::string >
CompFlowProblemGaussHump::analyticFieldNames( ncomp_t ) const
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
CompFlowProblemGaussHump::names( ncomp_t ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}

std::vector< std::array< tk::real, 3 > >
CompFlowProblemGaussHump::prescribedVelocity( ncomp_t, ncomp_t ncomp, tk::real,
                                             tk::real, tk::real )
// *****************************************************************************
//! Assign prescribed velocity at a point
//! \param[in] ncomp Number of components in this transport equation
//! \return Velocity assigned to all vertices of a tetrehedron, size:
//!   ncomp * ndim = [ncomp][3]
// *****************************************************************************
{
  std::vector< std::array< tk::real, 3 > > vel( ncomp );

  for (ncomp_t c=0; c<ncomp; ++c)
    vel[c] = {{ 1, 1, 0.0 }};

  return vel;
}
