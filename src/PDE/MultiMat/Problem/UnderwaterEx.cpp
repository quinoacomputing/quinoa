// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/UnderwaterEx.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the multi-material flow equations
  \details   This file defines a Problem policy class for the multi-material
    compressible flow equations, defined in PDE/MultiMat/DGMultiMat.hpp. See
    PDE/MultiMat/Problem.hpp for general requirements on Problem policy classes
    for MultiMat.
*/
// *****************************************************************************

#include "UnderwaterEx.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiMatProblemUnderwaterEx;

tk::InitializeFn::result_type
MultiMatProblemUnderwaterEx::initialize( ncomp_t ncomp,
                                         const std::vector< EOS >& mat_blk,
                                         tk::real x,
                                         tk::real y,
                                         tk::real z,
                                         tk::real )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
//! \details This function only initializes the underwater explosion problem,
//!   but does not actually give the analytical solution at time greater than 0.
// *****************************************************************************
{
  // see also Control/Inciter/InputDeck/Grammar.hpp
  Assert( ncomp == 12, "Number of scalar components must be 12" );

  auto nmat = g_inputdeck.get< eq, tag::nmat >();
  auto alphamin = g_inputdeck.get< eq, tag::min_volumefrac >();

  std::vector< tk::real > s(ncomp, 0.0), r(nmat, 0.0);
  tk::real p, u, v, w, temp;

  // velocity
  u = 0.0;
  v = 0.0;
  w = 0.0;

  // background state (air)
  // volume-fraction
  s[volfracIdx(nmat, 0)] = alphamin;
  s[volfracIdx(nmat, 1)] = alphamin;
  s[volfracIdx(nmat, 2)] = 1.0-2.0*alphamin;
  // pressure
  p = 1.01325e5;
  // temperature
  temp = 288.2;

  auto radb = std::sqrt((y-1.0)*(y-1.0) + x*x + z*z);

  // high-pressure gas bubble
  if (radb <= 0.3) {
    // volume-fraction
    s[volfracIdx(nmat, 0)] = alphamin;
    s[volfracIdx(nmat, 1)] = 1.0-2.0*alphamin;
    s[volfracIdx(nmat, 2)] = alphamin;
    // pressure
    p = 1.0e9;
    // temperature
    temp = 2000.0;
  }
  // water level
  else if (y <= 1.5) {
    // volume-fraction
    s[volfracIdx(nmat, 0)] = 1.0-2.0*alphamin;
    s[volfracIdx(nmat, 1)] = alphamin;
    s[volfracIdx(nmat, 2)] = alphamin;
    // temperature
    temp = 185.52;
  }

  auto rb(0.0);
  for (std::size_t k=0; k<nmat; ++k)
  {
    // densities
    r[k] = mat_blk[k].compute< EOS::density >( p, temp );
    // partial density
    s[densityIdx(nmat, k)] = s[volfracIdx(nmat, k)]*r[k];
    // total specific energy
    s[energyIdx(nmat, k)] =
      mat_blk[k].compute< EOS::totalenergy >( s[volfracIdx(nmat, k)]*r[k],
      u, v, w, s[volfracIdx(nmat, k)]*p, s[volfracIdx(nmat, k)] );
    rb += s[densityIdx(nmat, k)];
  }

  s[momentumIdx(nmat, 0)] = rb * u;
  s[momentumIdx(nmat, 1)] = rb * v;
  s[momentumIdx(nmat, 2)] = rb * w;

  return s;
}
