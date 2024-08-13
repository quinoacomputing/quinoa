// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem/GaussHump.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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

#include "GaussHump.hpp"

using inciter::TransportProblemGaussHump;

std::vector< tk::real >
TransportProblemGaussHump::initialize( ncomp_t ncomp,
  const std::vector< EOS >&, tk::real x, tk::real y, tk::real,
  tk::real t )
// *****************************************************************************
//  Evaluate analytical solution at (x,y,t) for all components
//! \param[in] ncomp Number of components in this transport equation system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x,y,t)
// *****************************************************************************
{
  const auto vel = prescribedVelocity( ncomp, x, y, 0.0, t );

  std::vector< tk::real > s( ncomp, 0.0 );
  for (ncomp_t c=0; c<ncomp; ++c)
  {
    // center of the hump
    auto x0 = 0.25 + vel[c][0]*t;
    auto y0 = 0.25 + vel[c][1]*t;

    // hump
    s[c] = 1.0 * exp( -((x-x0)*(x-x0) + (y-y0)*(y-y0))/(2.0 * 0.005) );
  }
  return s;
}

std::vector< std::array< tk::real, 3 > >
TransportProblemGaussHump::prescribedVelocity( ncomp_t ncomp, tk::real,
                                               tk::real, tk::real, tk::real )
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
