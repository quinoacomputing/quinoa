// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/GaussHumpCompflow.cpp
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

#include "GaussHumpCompflow.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "EoS/EoS.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemGaussHump;

tk::SolutionFn::result_type
CompFlowProblemGaussHump::solution( ncomp_t system,
                                    ncomp_t ncomp,
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
//! \note The function signature must follow tk::SolutionFn
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
  rE = eos_totalenergy< eq >( system, r, u, v, w, p );

  return {{ r, r*u, r*v, r*w, rE }};
}

std::vector< tk::real >
CompFlowProblemGaussHump::solinc( ncomp_t system, ncomp_t ncomp, tk::real x,
        tk::real y, tk::real z, tk::real t, tk::real dt ) const
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
CompFlowProblemGaussHump::src( ncomp_t, ncomp_t, tk::real,
                               tk::real, tk::real, tk::real )
// *****************************************************************************
//  Compute and return source term for manufactured solution
//! \return Array of reals containing the source for all components
//! \note The function signature must follow tk::SrcFn
// *****************************************************************************
{
  return {{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
}

void
CompFlowProblemGaussHump::side( std::unordered_set< int >& conf ) const
// *****************************************************************************
//  Query all side set IDs the user has configured for all components in this
//  PDE system
//! \param[in,out] conf Set of unique side set IDs to add to
// *****************************************************************************
{
  using tag::param;

  for (const auto& s : g_inputdeck.get< param, eq, tag::bcinlet >())
    for (const auto& i : s)
      conf.insert( std::stoi(i) );

  for (const auto& s : g_inputdeck.get< param, eq, tag::bcsubsonicoutlet >())
    for (const auto& i : s)
      conf.insert( std::stoi(i) );

  for (const auto& s : g_inputdeck.get< param, eq, tag::bcextrapolate >())
    for (const auto& i : s)
      conf.insert( std::stoi(i) );

  for (const auto& s : g_inputdeck.get< param, eq, tag::bcdir >())
    for (const auto& i : s)
      conf.insert( std::stoi(i) );
}

std::vector< std::string >
CompFlowProblemGaussHump::fieldNames( ncomp_t ) const
// *****************************************************************************
// Return field names to be output to file
//! \return Vector of strings labelling fields output in file
// *****************************************************************************
{
  const auto pref = inciter::g_inputdeck.get< tag::pref, tag::pref >();

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

  if(pref)
    n.push_back( "number of degree of freedom" );

  return n;
}

std::vector< std::vector< tk::real > >
CompFlowProblemGaussHump::fieldOutput(
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
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();

  std::vector< std::vector< tk::real > > out;
  auto r = U.extract( 0*rdof, offset );
  auto u = U.extract( 1*rdof, offset );
  auto v = U.extract( 2*rdof, offset );
  auto w = U.extract( 3*rdof, offset );
  auto e = U.extract( 4*rdof, offset );

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
  std::transform( r.begin(), r.end(), e.begin(), e.begin(),
                  []( tk::real s, tk::real& d ){ return d /= s; } );
  out.push_back( e );

  auto p = r;
  for (std::size_t i=0; i<r.size(); ++i)
    p[i] = eos_pressure< eq >( system, r[i], u[i], v[i], w[i], r[i]*e[i] );
  out.push_back( p );

  auto er = r;
  for (std::size_t i=0; i<r.size(); ++i) {
    auto s = solution( system, ncomp, x[i], y[i], z[i], t );
    er[i] = std::pow( r[i] - s[0], 2.0 ) * vol[i] / V;
    r[i] = s[0];
    u[i] = s[1]/s[0];
    v[i] = s[2]/s[0];
    w[i] = s[3]/s[0];
    e[i] = s[4]/s[0];
    p[i] = eos_pressure< eq >( system, r[i], u[i], v[i], w[i], r[i]*e[i] );
  }

  out.push_back( r );
  out.push_back( u );
  out.push_back( v );
  out.push_back( w );
  out.push_back( e );
  out.push_back( p );

  out.push_back( er );

  return out;
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
