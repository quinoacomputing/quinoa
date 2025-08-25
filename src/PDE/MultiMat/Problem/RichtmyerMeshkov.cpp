// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/RichtmyerMeshkov.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a Problem policy class for the multi-material
    compressible flow equations, defined in PDE/MultiMat/MultiMat.hpp. See
    PDE/MultiMat/Problem.hpp for general requirements on Problem policy classes
    for MultiMat.
*/
// *****************************************************************************

#include "RichtmyerMeshkov.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiMatProblemRichtmyerMeshkov;

tk::InitializeFn::result_type
MultiMatProblemRichtmyerMeshkov::initialize( ncomp_t ncomp,
  const std::vector< EOS >& mat_blk,
  tk::real x, tk::real y, tk::real,
  tk::real )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
//! \note The function signature must follow tk::InitializeFn
//! \details This function only initializes the Richtmyer-Meshkov instability
//!   problem, but does not actually give the analytical solution at time
//!   greater than 0.
// *****************************************************************************
{
  // see also Control/Inciter/InputDeck/Grammar.hpp
  Assert( ncomp == 9, "Number of scalar components must be 9" );

  auto nmat = g_inputdeck.get< eq, tag::nmat >();
  auto alphamin = g_inputdeck.get< eq, tag::min_volumefrac >();

  std::vector< tk::real > s( ncomp, 0.0 );
  tk::real p, T, u, v, w;

  // unshocked state
  p = 95600.0;
  T = 296.0;
  u = 0.0;
  v = 0.0;
  w = 0.0;
  s[volfracIdx(nmat,0)] = 1.0-alphamin;

  // shocked state
  if (x < 0.01) {
    p = 145347.8785;
    T = 324.6772971;
    u = 101.2761625;
  }

  // location of interface
  auto pi = 4.0 * std::atan(1.0);
  auto wvln = 0.059333;
  auto ampl = 0.002;
  if (x>0.03+ampl*std::sin(2.0*pi*y/wvln + pi/2.0)) {
    s[volfracIdx(nmat,0)] = alphamin;
  }

  // volume-fraction of 2nd material
  s[volfracIdx(nmat, 1)] = 1.0 - s[volfracIdx(nmat, 0)];

  auto rb = 0.0;
  for (std::size_t k=0; k<nmat; ++k) {
    // density
    auto r = mat_blk[k].compute< EOS::density >(p, T);
    s[densityIdx(nmat, k)] = s[volfracIdx(nmat, k)]*r;
    rb += s[densityIdx(nmat, k)];
    // total specific energy
    s[energyIdx(nmat, k)] =
      mat_blk[k].compute< EOS::totalenergy >( s[volfracIdx(nmat, k)]*r, u, v, w,
      s[volfracIdx(nmat, k)]*p, s[volfracIdx(nmat, k)] );
  }
  // momentum
  s[momentumIdx(nmat, 0)] = rb*u;
  s[momentumIdx(nmat, 1)] = rb*v;
  s[momentumIdx(nmat, 2)] = rb*w;

  return s;
}
