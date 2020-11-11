// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/RayleighTaylor.cpp
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

#include "RayleighTaylor.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemRayleighTaylor;

tk::SolutionFn::result_type
CompFlowProblemRayleighTaylor::solution( ncomp_t system,
                                         ncomp_t,
                                         tk::real x,
                                         tk::real y,
                                         tk::real z,
                                         tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
//! \note The function signature must follow tk::SolutionFn
// *****************************************************************************
{
  using tag::param; using std::sin; using std::cos;

  // manufactured solution parameters
  const auto a = g_inputdeck.get< param, eq, tag::alpha >()[system];
  const auto bx = g_inputdeck.get< param, eq, tag::betax >()[system];
  const auto by = g_inputdeck.get< param, eq, tag::betay >()[system];
  const auto bz = g_inputdeck.get< param, eq, tag::betaz >()[system];
  const auto p0 = g_inputdeck.get< param, eq, tag::p0 >()[system];
  const auto r0 = g_inputdeck.get< param, eq, tag::r0 >()[system];
  const auto k = g_inputdeck.get< param, eq, tag::kappa >()[system];
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
  const tk::real rE = eos_totalenergy< eq >( system, r, u, v, w, p );

  return {{ r, r*u, r*v, r*w, rE }};
}

std::vector< std::string >
CompFlowProblemRayleighTaylor::fieldNames( ncomp_t ) const
// *****************************************************************************
// Return field names to be output to file
//! \return Vector of strings labelling fields output in file
// *****************************************************************************
{
  auto n = CompFlowFieldNames();

  n.push_back( "density_analytical" );
  n.push_back( "x-velocity_analytical" );
  n.push_back( "y-velocity_analytical" );
  n.push_back( "z-velocity_analytical" );
  n.push_back( "specific_total_energy_analytical" );
  n.push_back( "pressure_analytical" );
  n.push_back( "err(rho)" );
  n.push_back( "err(e)" );
  n.push_back( "err(p)" );
  n.push_back( "err(u)" );
  n.push_back( "err(v)" );
  n.push_back( "err(w)" );

  return n;
}

std::vector< std::vector< tk::real > >
CompFlowProblemRayleighTaylor::fieldOutput(
  ncomp_t system,
  ncomp_t ncomp,
  ncomp_t offset,
  std::size_t nunk,
  std::size_t rdof,
  tk::real t,
  tk::real V,
  const std::vector< tk::real >& vol,
  const std::array< std::vector< tk::real >, 3 >& coord,
  const tk::Fields& U ) const
// *****************************************************************************
//  Return field output going to file
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset System offset specifying the position of the system of
//!   PDEs among other systems
//! \param[in] nunk Number of unknowns to extract
//! \param[in] rdof Number of reconstructed degrees of freedom. This is used as
//!   the number of scalar components to shift when extracting scalar
//!   components.
//! \param[in] t Physical time
//! \param[in] V Total mesh volume (across the whole problem)
//! \param[in] vol Nodal mesh volumes
//! \param[in] coord Mesh node coordinates
//! \param[in] U Solution vector at recent time step
//! \return Vector of vectors to be output to file
// *****************************************************************************
{
  auto out = CompFlowFieldOutput( system, offset, nunk, rdof, U );

  auto r = U.extract( 0*rdof, offset );
  auto u = U.extract( 1*rdof, offset );
  auto v = U.extract( 2*rdof, offset );
  auto w = U.extract( 3*rdof, offset );
  auto E = U.extract( 4*rdof, offset );

  // mesh node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  auto er = r, ee = r, ep = r, eu = r, ev = r, ew = r, p = r;
  for (std::size_t i=0; i<nunk; ++i) {
    auto s = solution( system, ncomp, x[i], y[i], z[i], t );
    er[i] = std::pow( r[i] - s[0], 2.0 ) * vol[i] / V;
    ee[i] = std::pow( E[i] - s[4]/s[0], 2.0 ) * vol[i] / V;
    eu[i] = std::pow( u[i] - s[1]/s[0], 2.0 ) * vol[i] / V;
    ev[i] = std::pow( v[i] - s[2]/s[0], 2.0 ) * vol[i] / V;
    ew[i] = std::pow( w[i] - s[3]/s[0], 2.0 ) * vol[i] / V;
    auto ap = eos_pressure< eq >( system, s[0], s[1]/s[0], s[2]/s[0], s[3]/s[0],
                                  s[4] );
    r[i] = s[0];
    u[i] = s[1]/s[0];
    v[i] = s[2]/s[0];
    w[i] = s[3]/s[0];
    E[i] = s[4]/s[0];
    p[i] = eos_pressure< eq >( system, r[i], u[i], v[i], w[i], r[i]*E[i] );
    ep[i] = std::pow( ap - p[i], 2.0 ) * vol[i] / V;
  }

  out.push_back( r );
  out.push_back( u );
  out.push_back( v );
  out.push_back( w );
  out.push_back( E );
  out.push_back( p );

  out.push_back( er );
  out.push_back( ee );
  out.push_back( ep );
  out.push_back( eu );
  out.push_back( ev );
  out.push_back( ew );

  return out;
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
