// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/NLEnergyGrowth.cpp
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

#include "NLEnergyGrowth.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "EoS/EoS.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemNLEnergyGrowth;

tk::real
CompFlowProblemNLEnergyGrowth::hx( tk::real bx, tk::real by, tk::real bz,
                                   tk::real x, tk::real y, tk::real z )
// *****************************************************************************
//  Compute internal energy parameter
//! \param[in] bx Parameter betax
//! \param[in] by Parameter betay
//! \param[in] bz Parameter betaz
//! \param[in] x X coordinate to evaluate at
//! \param[in] y Y coordinate to evaluate at
//! \param[in] z Z coordinate to evaluate at
//! \return Internal energy parameter
// *****************************************************************************
{
  return std::cos(bx*M_PI*x) * std::cos(by*M_PI*y) * std::cos(bz*M_PI*z);
}

tk::real
CompFlowProblemNLEnergyGrowth::ec( tk::real ce, tk::real kappa, tk::real t,
                                   tk::real h, tk::real p )
// *****************************************************************************
//  Compute a power of the internal energy
//! \param[in] ce Internal energy parameter
//! \param[in] kappa Internal energy parameter
//! \param[in] t Physical time
//! \param[in] h Internal energy parameter
//! \param[in] p Power
//! \return Internal energy raised to power p
// *****************************************************************************
{
  return std::pow( -3.0*(ce + kappa*h*h*t), p );
}

tk::SolutionFn::result_type
CompFlowProblemNLEnergyGrowth::solution( ncomp_t system,
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
  using tag::param;

  // manufactured solution parameters
  const auto ce = g_inputdeck.get< param, eq, tag::ce >()[system];
  const auto r0 = g_inputdeck.get< param, eq, tag::r0 >()[system];
  const auto a = g_inputdeck.get< param, eq, tag::alpha >()[system];
  const auto k = g_inputdeck.get< param, eq, tag::kappa >()[system];
  const auto bx = g_inputdeck.get< param, eq, tag::betax >()[system];
  const auto by = g_inputdeck.get< param, eq, tag::betay >()[system];
  const auto bz = g_inputdeck.get< param, eq, tag::betaz >()[system];
  // spatial component of density field
  const tk::real gx = 1.0 - x*x - y*y - z*z;
  // internal energy parameter
  const auto h = hx( bx, by, bz, x, y, z );
  // temporal component of the density field
  tk::real ft = std::exp( -a*t );
  // solution at t
  auto r = r0 + ft*gx;

  return {{ r, 0.0, 0.0, 0.0, r*ec(ce,k,t,h,-1.0/3.0) }};
}

std::vector< tk::real >
CompFlowProblemNLEnergyGrowth::solinc( ncomp_t system, ncomp_t ncomp,
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
CompFlowProblemNLEnergyGrowth::src( ncomp_t system, ncomp_t ncomp, tk::real x,
                                    tk::real y, tk::real z, tk::real t )
// *****************************************************************************
//  Compute and return source term for NLEG manufactured solution
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
  const auto a = g_inputdeck.get< param, eq, tag::alpha >()[system];
  const auto bx = g_inputdeck.get< param, eq, tag::betax >()[system];
  const auto by = g_inputdeck.get< param, eq, tag::betay >()[system];
  const auto bz = g_inputdeck.get< param, eq, tag::betaz >()[system];
  const auto ce = g_inputdeck.get< param, eq, tag::ce >()[system];
  const auto kappa = g_inputdeck.get< param, eq, tag::kappa >()[system];
  const auto r0 = g_inputdeck.get< param, eq, tag::r0 >()[system];
  // ratio of specific heats
  const auto g = g_inputdeck.get< param, eq, tag::gamma >()[system][0];
  // spatial component of density field
  const auto gx = 1.0 - x*x - y*y - z*z;
  // derivative of spatial component of density field
  const std::array< tk::real, 3 > dg{{ -2.0*x, -2.0*y, -2.0*z }};
  // spatial component of energy field
  const auto h = hx( bx, by, bz, x, y, z );
  // derivative of spatial component of energy field
  std::array< tk::real, 3 >
    dh{{ -bx*M_PI*sin(bx*M_PI*x)*cos(by*M_PI*y)*cos(bz*M_PI*z),
         -by*M_PI*cos(bx*M_PI*x)*sin(by*M_PI*y)*cos(bz*M_PI*z),
         -bz*M_PI*cos(bx*M_PI*x)*cos(by*M_PI*y)*sin(bz*M_PI*z) }};
  // temporal function f and its derivative
  const auto ft = std::exp(-a*t);
  const auto dfdt = -a*ft;
  // density and its derivatives
  const auto rho = r0 + ft*gx;
  const std::array< tk::real, 3 > drdx{{ ft*dg[0], ft*dg[1], ft*dg[2] }};
  const auto drdt = gx*dfdt;
  // internal energy and its derivatives
  const auto ie = ec( ce, kappa, t, h, -1.0/3.0 );
  const std::array< tk::real, 3 > dedx{{
    2.0*std::pow(ie,4.0)*kappa*h*dh[0]*t,
    2.0*std::pow(ie,4.0)*kappa*h*dh[1]*t,
    2.0*std::pow(ie,4.0)*kappa*h*dh[2]*t }};
  const auto dedt = kappa*h*h*std::pow(ie,4.0);
  // sources
  std::vector< tk::real > r( ncomp );
  // density source
  r[0] = drdt;
  // momentum source
  r[1] = (g-1.0)*(rho*dedx[0] + ie*drdx[0]);
  r[2] = (g-1.0)*(rho*dedx[1] + ie*drdx[1]);
  r[3] = (g-1.0)*(rho*dedx[2] + ie*drdx[2]);
  // energy source
  r[4] = rho*dedt + ie*drdt;
  return r;
}

void
CompFlowProblemNLEnergyGrowth::side( std::unordered_set< int >& conf ) const
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
CompFlowProblemNLEnergyGrowth::fieldNames( ncomp_t ) const
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
  n.push_back( "err(e)" );

  if(pref)
    n.push_back( "number of degree of freedom" );

  return n;
}

std::vector< std::vector< tk::real > >
CompFlowProblemNLEnergyGrowth::fieldOutput(
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
  // number of degree of freedom
  const std::size_t rdof =
    g_inputdeck.get< tag::discr, tag::rdof >();

  std::vector< std::vector< tk::real > > out;
  auto r = U.extract( 0*rdof, offset );
  auto u = U.extract( 1*rdof, offset );
  auto v = U.extract( 2*rdof, offset );
  auto w = U.extract( 3*rdof, offset );
  auto E = U.extract( 4*rdof, offset );

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

  auto er = r, ee = r;
  for (std::size_t i=0; i<r.size(); ++i) {
    auto s = solution( system, ncomp, x[i], y[i], z[i], t );
    er[i] = std::pow( r[i] - s[0], 2.0 ) * vol[i] / V;
    ee[i] = std::pow( E[i] - s[4]/s[0], 2.0 ) * vol[i] / V;
    r[i] = s[0];
    u[i] = s[1]/s[0];
    v[i] = s[2]/s[0];
    w[i] = s[3]/s[0];
    E[i] = s[4]/s[0];
    p[i] = eos_pressure< eq >( system, r[i], u[i], v[i], w[i], r[i]*E[i] );
  }

  out.push_back( r );
  out.push_back( u );
  out.push_back( v );
  out.push_back( w );
  out.push_back( E );
  out.push_back( p );

  out.push_back( er );
  out.push_back( ee );

  return out;
}

std::vector< std::string >
CompFlowProblemNLEnergyGrowth::names( ncomp_t ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
