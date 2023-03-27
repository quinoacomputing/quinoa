// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem/CylVortex.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for transport equations
  \details   This file defines a Problem policy class for the transport
    equations, defined in PDE/Transport/CGTransport.hpp implementing
    node-centered continuous Galerkin (CG) and PDE/Transport/DGTransport.hpp
    implementing cell-centered discontinuous Galerkin (DG) discretizations.
    See PDE/Transport/Problem.hpp for general requirements on Problem policy
    classes for cg::Transport and dg::Transport.
*/
// *****************************************************************************

#include "CylVortex.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::TransportProblemCylVortex;

std::vector< tk::real >
TransportProblemCylVortex::initialize( ncomp_t system, ncomp_t ncomp,
                                       const std::vector< EOS >&,
                                       tk::real x, tk::real y, tk::real,
                                       tk::real t )
// *****************************************************************************
//  Evaluate initial solution at (x,y,t) for all components
//! \param[in] system Equation system index
//! \param[in] ncomp Number of components in this transport equation system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x,y,t=0)
//! \details This function only gives the initial condition for the cylinder,
//!   and not the solution at any time t>0.
// *****************************************************************************
{
  const auto vel = prescribedVelocity( system, ncomp, x, y, 0.0, t );

  if (ncomp != 4) Throw("Cylinder deformation in vortex is only set up for 4 "
    "components");

  std::vector< tk::real > s( ncomp, 0.0 );

  // center of the cylinder
  auto x0 = 0.5;
  auto y0 = 0.75;
  auto r = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));

  if (r<=0.15) {
    if (x<x0 && y>=y0) s[0] = 1.0;
    else if (x>=x0 && y>=y0) s[1] = 1.0;
    else if (x>=x0 && y<y0) s[2] = 1.0;
    else if (x<x0 && y<y0) s[3] = 1.0;
  }

  return s;
}

std::vector< std::array< tk::real, 3 > >
TransportProblemCylVortex::prescribedVelocity( ncomp_t, ncomp_t ncomp,
  tk::real x, tk::real y, tk::real, tk::real t )
// *****************************************************************************
//! Assign prescribed velocity at a point
//! \param[in] ncomp Number of components in this transport equation
//! \param[in] x x coordinate at which to assign velocity
//! \param[in] y y coordinate at which to assign velocity
//! \param[in] t time at which to assign velocity
//! \return Velocity assigned to all vertices of a tetrehedron, size:
//!   ncomp * ndim = [ncomp][3]
// *****************************************************************************
{
  std::vector< std::array< tk::real, 3 > > vel( ncomp );

  auto pi = 4.0 * std::atan(1.0);
  for (ncomp_t c=0; c<ncomp; ++c) {
    vel[c][0] = - 2.0*std::cos(t*pi/4.0) * std::pow(std::sin(pi*x), 2)
      * std::sin(pi*y) * std::cos(pi*y);
    vel[c][1] = 2.0*std::cos(t*pi/4.0) * std::pow(std::sin(pi*y), 2)
      * std::sin(pi*x) * std::cos(pi*x);
    vel[c][2] = 0.0;
  }

  return vel;
}
