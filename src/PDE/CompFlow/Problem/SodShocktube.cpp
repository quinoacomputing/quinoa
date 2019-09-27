// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/SodShocktube.cpp
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

#include "SodShocktube.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "EoS/EoS.hpp"

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
                                       tk::real )
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

std::vector< tk::real >
CompFlowProblemSodShocktube::solinc( ncomp_t system, ncomp_t ncomp, tk::real x,
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
CompFlowProblemSodShocktube::src( ncomp_t, ncomp_t, tk::real,
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
CompFlowProblemSodShocktube::side( std::unordered_set< int >& conf ) const
// *****************************************************************************
//  Query all side set IDs the user has configured for all components in this
//  PDE system
//! \param[in,out] conf Set of unique side set IDs to add to
// *****************************************************************************
{
  using tag::param;

  for (const auto& s : g_inputdeck.get< param, eq, tag::bcextrapolate >())
    for (const auto& i : s) conf.insert( std::stoi(i) );

  for (const auto& s : g_inputdeck.get< param, eq, tag::bcsym >())
    for (const auto& i : s) conf.insert( std::stoi(i) );
}

std::vector< std::string >
CompFlowProblemSodShocktube::fieldNames( ncomp_t ) const
// *****************************************************************************
// Return field names to be output to file
//! \return Vector of strings labelling fields output in file
// *****************************************************************************
{
  const auto pref = inciter::g_inputdeck.get< tag::pref, tag::pref >();

  std::vector< std::string > n;

  n.push_back( "density_numerical" );
  //n.push_back( "density_analytical" );
  n.push_back( "x-velocity_numerical" );
  //n.push_back( "x-velocity_analytical" );
  //n.push_back( "err(u)" );
  n.push_back( "y-velocity_numerical" );
  //n.push_back( "y-velocity_analytical" );
  n.push_back( "z-velocity_numerical" );
  //n.push_back( "z-velocity_analytical" );
  n.push_back( "specific_total_energy_numerical" );
  //n.push_back( "specific_total_energy_analytical" );
  //n.push_back( "err(E)" );
  n.push_back( "pressure_numerical" );
  //n.push_back( "pressure_analytical" );

  if(pref)
    n.push_back( "number of degree of freedom" );

  return n;
}

std::vector< std::vector< tk::real > >
CompFlowProblemSodShocktube::fieldOutput(
  ncomp_t system,
  ncomp_t,
  ncomp_t offset,
  tk::real,
  tk::real,
  const std::vector< tk::real >&,
  const std::array< std::vector< tk::real >, 3 >&,
  tk::Fields& U ) const
// *****************************************************************************
//  Return field output going to file
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] offset System offset specifying the position of the system of
//!   PDEs among other systems
//! \param[in] U Solution vector at recent time step
//! \return Vector of vectors to be output to file
// *****************************************************************************
{
  // number of degree of freedom
  const std::size_t rdof =
    g_inputdeck.get< tag::discr, tag::rdof >();

  std::vector< std::vector< tk::real > > out;
  const auto r  = U.extract( 0*rdof, offset );
  const auto ru = U.extract( 1*rdof, offset );
  const auto rv = U.extract( 2*rdof, offset );
  const auto rw = U.extract( 3*rdof, offset );
  const auto re = U.extract( 4*rdof, offset );

  // mesh node coordinates
  //const auto& x = coord[0];
  //const auto& y = coord[1];

  out.push_back( r );
  //out.push_back( std::vector< tk::real >( r.size(), 1.0 ) );

  std::vector< tk::real > u = ru;
  std::transform( r.begin(), r.end(), u.begin(), u.begin(),
                  []( tk::real s, tk::real& d ){ return d /= s; } );
  out.push_back( u );
  //std::vector< tk::real > ua = ru;
  //for (std::size_t i=0; i<ua.size(); ++i)
  //  ua[i] = std::sin(M_PI*x[i]) * std::cos(M_PI*y[i]);
  //out.push_back( ua );

  //// error in x-velocity
  //auto err = u;
  //for (std::size_t i=0; i<u.size(); ++i)
  //   err[i] = std::pow( ua[i] - u[i], 2.0 ) * vol[i] / V;
  // out.push_back( err );

  std::vector< tk::real > v = rv;
  //std::vector< tk::real > va = rv;
  std::transform( r.begin(), r.end(), v.begin(), v.begin(),
                  []( tk::real s, tk::real& d ){ return d /= s; } );
  out.push_back( v );
  //for (std::size_t i=0; i<va.size(); ++i)
  //  va[i] = -std::cos(M_PI*x[i]) * std::sin(M_PI*y[i]);
  //out.push_back( va );

  std::vector< tk::real > w = rw;
  //std::vector< tk::real > wa = rw;
  std::transform( r.begin(), r.end(), w.begin(), w.begin(),
                  []( tk::real s, tk::real& d ){ return d /= s; } );
  out.push_back( w );
  //for (std::size_t i=0; i<wa.size(); ++i)
  //  wa[i] = 0.0;
  //out.push_back( wa );

  std::vector< tk::real > E = re;
  //std::vector< tk::real > Ea = re;
  //std::vector< tk::real > Pa( r.size(), 0.0 );
  std::transform( r.begin(), r.end(), E.begin(), E.begin(),
                  []( tk::real s, tk::real& d ){ return d /= s; } );
  out.push_back( E );
  //for (std::size_t i=0; i<Ea.size(); ++i) {
  //  Pa[i] = 10.0 +
  //    r[i]/4.0*(std::cos(2.0*M_PI*x[i]) + std::cos(2.0*M_PI*y[i]));
  //  Ea[i] = Pa[i]/(g-1.0)/r[i] +
  //          0.5*(ua[i]*ua[i] + va[i]*va[i] + wa[i]*wa[i])/r[i];
  //}
  //out.push_back( Ea );

  //// error in total specific energy
  //for (std::size_t i=0; i<v.size(); ++i)
  //  err[i] = std::pow( Ea[i] - E[i], 2.0 ) * vol[i] / V;
  //out.push_back( err );

  std::vector< tk::real > P( r.size(), 0.0 );
  for (std::size_t i=0; i<P.size(); ++i)
    P[i] = eos_pressure< eq >( system, r[i], u[i], v[i], w[i], r[i]*E[i] );
  out.push_back( P );
  //out.push_back( Pa );

  return out;
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
