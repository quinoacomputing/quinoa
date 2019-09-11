// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/RayleighTaylor.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a Problem policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problem.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************

#include "RayleighTaylor.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "EoS/EoS.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemRayleighTaylor;

tk::SolutionFn::result_type
CompFlowProblemRayleighTaylor::solution( ncomp_t system,
                                         [[maybe_unused]] ncomp_t ncomp,
                                         tk::real x,
                                         tk::real y,
                                         tk::real z,
                                         tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
//! \note The function signature must follow tk::SolutionFn
// *****************************************************************************
{
  Assert( ncomp == ncomp, "Number of scalar components must be " +
                          std::to_string(ncomp) );
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

std::vector< tk::real >
CompFlowProblemRayleighTaylor::solinc( ncomp_t system, ncomp_t ncomp,
  tk::real x, tk::real y, tk::real z, tk::real t, tk::real dt ) const
// *****************************************************************************
// Evaluate the increment from t to t+dt of the analytical solution at (x,y,z)
// for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution increment starting from
//! \param[in] dt Time increment at which evaluate the solution increment to
//! \return Increment in values of all components evaluated at (x,y,z,t+dt)
// *****************************************************************************
{
  auto st1 = solution( system, ncomp, x, y, z, t );
  auto st2 = solution( system, ncomp, x, y, z, t+dt );
  std::transform( begin(st1), end(st1), begin(st2), begin(st2),
                  []( tk::real s, tk::real& d ){ return d -= s; } );
  return st2;
}

tk::SrcFn::result_type
CompFlowProblemRayleighTaylor::src( ncomp_t system, ncomp_t ncomp, tk::real x,
                                    tk::real y, tk::real z, tk::real t )
// *****************************************************************************
//  Compute and return source term for manufactured solution
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Physical time at which to evaluate the source
//! \return Array of reals containing the source for all components
//! \note The function signature must follow tk::SrcFn
// *****************************************************************************
{
  using tag::param; using std::sin; using std::cos;

  // manufactured solution parameters
  auto a = g_inputdeck.get< param, eq, tag::alpha >()[system];
  auto bx = g_inputdeck.get< param, eq, tag::betax >()[system];
  auto by = g_inputdeck.get< param, eq, tag::betay >()[system];
  auto bz = g_inputdeck.get< param, eq, tag::betaz >()[system];
  auto k = g_inputdeck.get< param, eq, tag::kappa >()[system];
  auto p0 = g_inputdeck.get< param, eq, tag::p0 >()[system];
  // ratio of specific heats
  tk::real g = g_inputdeck.get< param, eq, tag::gamma >()[system][0];

  // evaluate solution at x,y,z,t
  auto s = solution( system, ncomp, x, y, z, t );

  // density, velocity, energy, pressure
  auto rho = s[0];
  auto u = s[1]/s[0];
  auto v = s[2]/s[0];
  auto w = s[3]/s[0];
  auto E = s[4]/s[0];
  auto p = p0 + a*(bx*x*x + by*y*y + bz*z*z);

  // spatial gradients
  std::array< tk::real, 3 > drdx{{ -2.0*bx*x, -2.0*by*y, -2.0*bz*z }};
  std::array< tk::real, 3 > dpdx{{ 2.0*a*bx*x, 2.0*a*by*y, 2.0*a*bz*z }};
  tk::real ft = cos(k*M_PI*t);
  std::array< tk::real, 3 > dudx{{ ft*M_PI*z*cos(M_PI*x),
                                   0.0,
                                   ft*sin(M_PI*x) }};
  std::array< tk::real, 3 > dvdx{{ 0.0,
                                   -ft*M_PI*z*sin(M_PI*y),
                                   ft*cos(M_PI*y) }};
  std::array< tk::real, 3 > dwdx{{ ft*M_PI*0.5*M_PI*z*z*sin(M_PI*x),
                                   ft*M_PI*0.5*M_PI*z*z*cos(M_PI*y),
                                  -ft*M_PI*z*(cos(M_PI*x) - sin(M_PI*y)) }};
  std::array< tk::real, 3 > dedx{{
    dpdx[0]/rho/(g-1.0) - p/(g-1.0)/rho/rho*drdx[0]
    + u*dudx[0] + v*dvdx[0] + w*dwdx[0],
    dpdx[1]/rho/(g-1.0) - p/(g-1.0)/rho/rho*drdx[1]
    + u*dudx[1] + v*dvdx[1] + w*dwdx[1],
    dpdx[2]/rho/(g-1.0) - p/(g-1.0)/rho/rho*drdx[2]
    + u*dudx[2] + v*dvdx[2] + w*dwdx[2] }};
  
  // time derivatives
  auto dudt = -k*M_PI*sin(k*M_PI*t)*z*sin(M_PI*x);
  auto dvdt = -k*M_PI*sin(k*M_PI*t)*z*cos(M_PI*y);
  auto dwdt =  k*M_PI*sin(k*M_PI*t)/2*M_PI*z*z*(cos(M_PI*x) - sin(M_PI*y));
  auto dedt = u*dudt + v*dvdt + w*dwdt;

  std::vector< tk::real > r( ncomp );
  // density source
  r[0] = u*drdx[0] + v*drdx[1] + w*drdx[2];
  // momentum source
  r[1] = rho*dudt+u*r[0]+dpdx[0] + s[1]*dudx[0]+s[2]*dudx[1]+s[3]*dudx[2];
  r[2] = rho*dvdt+v*r[0]+dpdx[1] + s[1]*dvdx[0]+s[2]*dvdx[1]+s[3]*dvdx[2];
  r[3] = rho*dwdt+w*r[0]+dpdx[2] + s[1]*dwdx[0]+s[2]*dwdx[1]+s[3]*dwdx[2];
  // energy source
  r[4] = rho*dedt + E*r[0] + s[1]*dedx[0]+s[2]*dedx[1]+s[3]*dedx[2]
         + u*dpdx[0]+v*dpdx[1]+w*dpdx[2];

  return r;
}

void
CompFlowProblemRayleighTaylor::side( std::unordered_set< int >& conf ) const
// *****************************************************************************
//  Query all side set IDs the user has configured for all components in this
//  PDE system
//! \param[in,out] conf Set of unique side set IDs to add to
// *****************************************************************************
{
  using tag::param; using tag::bcdir;

  for (const auto& s : g_inputdeck.get< param, eq, bcdir >())
    for (const auto& i : s)
      conf.insert( std::stoi(i) );
}

std::vector< std::string >
CompFlowProblemRayleighTaylor::fieldNames( ncomp_t ) const
// *****************************************************************************
// Return field names to be output to file
//! \return Vector of strings labelling fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;

  n.push_back( "density_numerical" );
  n.push_back( "x-velocity_numerical" );
  n.push_back( "y-velocity_numerical" );
  n.push_back( "z-velocity_numerical" );
  n.push_back( "specific_total_energy_numerical" );
  n.push_back( "pressure_numerical" );
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
  tk::real t,
  tk::real V,
  const std::vector< tk::real >& vol,
  const std::array< std::vector< tk::real >, 3 >& coord,
  tk::Fields& U ) const
// *****************************************************************************
//  Return field output going to file
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset System offset specifying the position of the system of
//!   PDEs among other systems
//! \param[in] t Physical time
//! \param[in] V Total mesh volume (across the whole problem)
//! \param[in] vol Nodal mesh volumes
//! \param[in] coord Mesh node coordinates
//! \param[in] U Solution vector at recent time step
//! \return Vector of vectors to be output to file
// *****************************************************************************
{
  std::vector< std::vector< tk::real > > out;
  auto r = U.extract( 0, offset );
  auto u = U.extract( 1, offset );
  auto v = U.extract( 2, offset );
  auto w = U.extract( 3, offset );
  auto E = U.extract( 4, offset );

  // mesh node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  out.push_back( r );
  std::transform( r.begin(), r.end(), u.begin(), u.begin(),
                  []( tk::real s, tk::real& d ){ return d /= s; } );
  out.push_back( u );
  std::transform( r.begin(), r.end(), v.begin(), v.begin(),
                  []( tk::real s, tk::real& d ){ return d /= s; } );
  out.push_back( v );
  std::transform( r.begin(), r.end(), w.begin(), w.begin(),
                  []( tk::real s, tk::real& d ){ return d /= s; } );
  out.push_back( w );
  std::transform( r.begin(), r.end(), E.begin(), E.begin(),
                  []( tk::real s, tk::real& d ){ return d /= s; } );
  out.push_back( E );

  auto p = r;
  for (std::size_t i=0; i<r.size(); ++i)
    p[i] = eos_pressure< eq >( system, r[i], u[i], v[i], w[i], r[i]*E[i] );
  out.push_back( p );

  auto er = r, ee = r, ep = r, eu = r, ev = r, ew = r;
  for (std::size_t i=0; i<r.size(); ++i) {
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
