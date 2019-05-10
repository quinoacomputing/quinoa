// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem/CylAdvect.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for transport equations
  \details   This file defines a Problem policy class for the transport
    equations, defined in PDE/Transport/CGTransport.h implementing
    node-centered continuous Galerkin (CG) and PDE/Transport/DGTransport.h
    implementing cell-centered discontinuous Galerkin (DG) discretizations.
    See PDE/Transport/Problem.h for general requirements on Problem policy
    classes for cg::Transport and dg::Transport.
*/
// *****************************************************************************

#include "CylAdvect.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::TransportProblemCylAdvect;

std::vector< tk::real >
TransportProblemCylAdvect::solution( ncomp_t system, ncomp_t ncomp,
          tk::real x, tk::real y, tk::real, tk::real t )
// *****************************************************************************
//  Evaluate analytical solution at (x,y,t) for all components
//! \param[in] system Equation system index
//! \param[in] ncomp Number of components in this transport equation system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x,y,t)
// *****************************************************************************
{
  const auto vel = prescribedVelocity( system, ncomp, x, y, 0.0 );

  std::vector< tk::real > s( ncomp, 0.0 );
  for (ncomp_t c=0; c<ncomp; ++c)
  {
    // center of the cylinder
    auto x0 = 0.25 + vel[c][0]*t;
    auto y0 = 0.25 + vel[c][1]*t;

    // square wave
    auto r = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
    if (r<0.2)
      s[c] = 1.0;
    else
      s[c] = 0.0;
  }

  return s;
}

std::vector< tk::real >
TransportProblemCylAdvect::solinc( ncomp_t, ncomp_t ncomp, tk::real x,
                                tk::real y, tk::real, tk::real t, tk::real dt )
const
// *****************************************************************************
//  Evaluate the increment from t to t+dt of the analytical solution at (x,y,z)
//  for all components
//! \param[in] ncomp Number of components in this transport equation system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution increment starting from
//! \param[in] dt Time increment at which evaluate the solution increment to
//! \return Increment in values of all components evaluated at (x,y,t+dt)
// *****************************************************************************
{
  auto st1 = solution( 0, ncomp, x, y, 0.0, t );
  auto st2 = solution( 0, ncomp, x, y, 0.0, t+dt );

  std::transform( begin(st1), end(st1), begin(st2), begin(st2),
                  []( tk::real s, tk::real& d ){ return d -= s; } );

  return st2;
}

void
TransportProblemCylAdvect::side( std::unordered_set< int >& conf ) const
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

  for (const auto& s : g_inputdeck.get< param, eq, tag::bcoutlet >())
    for (const auto& i : s)
      conf.insert( std::stoi(i) );

  for (const auto& s : g_inputdeck.get< param, eq, tag::bcextrapolate >())
    for (const auto& i : s)
      conf.insert( std::stoi(i) );

  for (const auto& s : g_inputdeck.get< param, eq, tag::bcdir >())
    for (const auto& i : s)
      conf.insert( std::stoi(i) );
}

std::vector< std::array< tk::real, 3 > >
TransportProblemCylAdvect::prescribedVelocity( ncomp_t, ncomp_t ncomp, tk::real,
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
    vel[c] = {{ 0.1, 0.1, 0.0 }};

  return vel;
}
