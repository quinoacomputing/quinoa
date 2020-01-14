// *****************************************************************************
/*!
  \file      src/PDE/Transport/Physics/CGAdvDiff.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Physics policy for the transport equations using continuous
             Galerkin (CG)
  \details   This file defines a Physics policy class for the transport
    equations, defined in PDE/Transport/CGTransport.h implementing
    node-centered continuous Galerkin (CG) discretizations.
    See PDE/Transport/Physics/CG.h for general requirements on Physics policy
    classes for cg::Transport.
*/
// *****************************************************************************

#include "CGAdvDiff.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::cg::TransportPhysicsAdvDiff;

tk::real
TransportPhysicsAdvDiff::diffusion_dt(
  ncomp_t e,
  ncomp_t ncomp,
  tk::real L,
  const std::vector< std::array< tk::real, 4 > >& ) const
// *****************************************************************************
//! Compute the minimum time step size based on the diffusion
//! \param[in] e Equation system index, i.e., which transport equation
//!   system we operate on among the systems of PDEs
//! \param[in] ncomp Number of components in this PDE
//! \param[in] L Characteristic length scale
//! \return Minimum time step size based on diffusion
// *****************************************************************************
{
  // diffusivities for all components
  const auto& df = g_inputdeck.get< tag::param, eq, tag::diffusivity >()[e];

  // compute the minimum diffusion time step size across the four nodes
  tk::real mindt = std::numeric_limits< tk::real >::max();
  for (ncomp_t c=0; c<ncomp; ++c) {
    const auto di = 3*c;
    const auto d = std::max( df[di+2], std::max( df[di+0], df[di+1] ) );
    const auto dt = L * L / (2.0*d);  // dt ~ dx^2/(2D)
    if (dt < mindt) mindt = dt;
  }

  return mindt;
}

