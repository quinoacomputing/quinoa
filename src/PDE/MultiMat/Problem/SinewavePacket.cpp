// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/SinewavePacket.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a Problem policy class for the compressible flow
    equations, defined in PDE/MultiMat/MultiMat.h. See PDE/MultiMat/Problem.h
    for general requirements on Problem policy classes for MultiMat.
*/
// *****************************************************************************

#include "SinewavePacket.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

}

using inciter::MultiMatProblemSinewavePacket;

tk::InitializeFn::result_type
MultiMatProblemSinewavePacket::initialize( ncomp_t ncomp,
                                           const std::vector< EOS >& mat_blk,
                                           tk::real x,
                                           tk::real /*y*/,
                                           tk::real /*z*/,
                                           tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x,t)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  auto nmat = g_inputdeck.get< tag::param, eq, tag::nmat >()[0];

  Assert(nmat == 1, "Sinewave packet advection not set up for more than one "
    "material");

  // see also Control/Inciter/InputDeck/Grammar.hpp
  Assert( ncomp == 3*nmat+3, "Incorrect number of components in multi-material "
          "system" );

  std::vector< tk::real > s( ncomp, 0.0 );
  auto u = 1.0;
  auto v = 0.0;
  auto w = 0.0;
  auto pi = 4.0 * std::atan(1.0);

  // location
  auto xi = x - u*t;

  // exact solution
  auto a = 20.0;
  auto c = 20.0;
  auto rhob = std::exp(-a*(xi-0.5)*(xi-0.5)) * std::sin(c*pi*xi) + 2.0;
  s[volfracIdx(nmat, 0)] = 1.0;
  s[densityIdx(nmat, 0)] = s[volfracIdx(nmat, 0)] * rhob;
  s[energyIdx(nmat, 0)] = s[volfracIdx(nmat, 0)]
    * mat_blk[0].compute< EOS::totalenergy >( rhob, u, v, w, 1.0 );
  s[momentumIdx(nmat, 0)] = rhob * u;
  s[momentumIdx(nmat, 1)] = rhob * v;
  s[momentumIdx(nmat, 2)] = rhob * w;

  return s;
}
