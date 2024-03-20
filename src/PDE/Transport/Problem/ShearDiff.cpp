// *****************************************************************************
/*!
  \file      src/PDE/Transport/Problem/ShearDiff.cpp
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

#include "ShearDiff.hpp"
#include "Inciter/InputDeck/New2InputDeck.hpp"

namespace inciter {

extern ctr::New2InputDeck g_inputdeck;

} // ::inciter

using inciter::TransportProblemShearDiff;

std::vector< tk::real >
TransportProblemShearDiff::initialize( ncomp_t ncomp,
  const std::vector< EOS >&, tk::real x, tk::real y, tk::real z,
  tk::real t )
// *****************************************************************************
//  Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] ncomp Number of components in this transport equation system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x,y,t)
// *****************************************************************************
{
  const auto& u0 = g_inputdeck.get< eq, newtag::u0 >();
  const auto& d = g_inputdeck.get< eq, newtag::diffusivity >();
  const auto& l = g_inputdeck.get< eq, newtag::lambda >();

  std::vector< tk::real > r( ncomp );
  for (ncomp_t c=0; c<ncomp; ++c) {
    const auto li = 2*c;
    const auto di = 3*c;
    const auto phi3s = (l[li+0]*l[li+0]*d[di+1]/d[di+0] +
                        l[li+1]*l[li+1]*d[di+2]/d[di+0]) / 12.0;
    r[c] =
        1.0 / ( 8.0 * std::pow(M_PI,3.0/2.0) *
                std::sqrt(d[di+0]*d[di+1]*d[di+2]) *
                std::pow(t,3.0/2.0) * std::sqrt(1.0+phi3s*t*t) ) *
        exp( -std::pow( x - u0[c]*t -
                        0.5*(l[li+0]*y + l[li+1]*z)*t, 2.0 ) /
              ( 4.0 * d[di+0] * t * (1.0 + phi3s*t*t) )
             -y*y / ( 4.0 * d[di+1] * t )
             -z*z / ( 4.0 * d[di+2] * t ) );
  }

  return r;
}

void
TransportProblemShearDiff::errchk( ncomp_t ncomp ) const
// *****************************************************************************
//  Do error checking on PDE parameters
//! \param[in] ncomp Number of components in this transport equation
// *****************************************************************************
{
  auto u0 = g_inputdeck.get< eq, newtag::u0 >();
  ErrChk( ncomp == u0.size(),
    "Wrong number of advection-diffusion PDE parameters 'u0'" );

  auto lambda = g_inputdeck.get< eq, newtag::lambda >();
  ErrChk( 2*ncomp == lambda.size(),
    "Wrong number of advection-diffusion PDE parameters 'lambda'" );

  auto& d = g_inputdeck.get< eq, newtag::diffusivity >();
  ErrChk( 3*ncomp == d.size(),
    "Wrong number of advection-diffusion PDE parameters 'diffusivity'" );
}

std::vector< std::array< tk::real, 3 > >
TransportProblemShearDiff::prescribedVelocity( ncomp_t ncomp,
                                               tk::real, tk::real y, tk::real z,
                                               tk::real )
// *****************************************************************************
//  Assign prescribed shear velocity at a point
//! \param[in] ncomp Number of components in this transport equation
//! \param[in] y y coordinate at which to assign velocity
//! \param[in] z Z coordinate at which to assign velocity
//! \return Velocity assigned to all vertices of a tetrehedron, size:
//!   ncomp * ndim = [ncomp][3]
// *****************************************************************************
{
  auto u0 = g_inputdeck.get< eq, newtag::u0 >();
  auto l = g_inputdeck.get< eq, newtag::lambda >();

  std::vector< std::array< tk::real, 3 > > vel( ncomp );
  for (ncomp_t c=0; c<ncomp; ++c)
    vel[c] = {{ u0[c] + l[2*c+0]*y + l[2*c+1]*z, 0.0, 0.0 }};

  return vel;
}
