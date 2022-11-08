// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/ShockHeBubble.cpp
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

#include "ShockHeBubble.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiMatProblemShockHeBubble;

tk::InitializeFn::result_type
MultiMatProblemShockHeBubble::initialize( ncomp_t system,
                                          ncomp_t ncomp,
                                          const std::vector< EOS >& mat_blk,
                                          tk::real x,
                                          tk::real y,
                                          tk::real z,
                                          tk::real )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which multi-material
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
//! \details This function only initializes the shock He-bubble interaction
//!   problem, but does not actually give the analytical solution at time
//!   greater than 0.
// *****************************************************************************
{
  // see also Control/Inciter/InputDeck/Grammar.hpp
  Assert( ncomp == 9, "Number of scalar components must be 9" );

  auto nmat =
    g_inputdeck.get< tag::param, eq, tag::nmat >()[system];

  std::vector< tk::real > s(ncomp, 0.0), r(nmat, 0.0);
  tk::real p, u, v, w, temp;
  auto alphamin = 1.0e-12;

  // pre-shock state
  s[volfracIdx(nmat, 0)] = alphamin;
  s[volfracIdx(nmat, 1)] = 1.0-alphamin;
  p = 1.0e5;
  u = 0.0;
  v = 0.0;
  w = 0.0;
  temp = 248.88;

  // post-shock state
  if (x >= 0.2125) {
    // volume-fraction
    s[volfracIdx(nmat, 0)] = alphamin;
    s[volfracIdx(nmat, 1)] = 1.0-alphamin;
    // pressure
    p = 1.5698e5;
    // velocity
    u = -113.5;
    // temperature
    temp = 283.86;
  }

  // bubble
  auto radb = std::sqrt((x-0.1725)*(x-0.1725) + y*y + z*z);
  if (radb <= 0.025) {
    // volume-fraction
    s[volfracIdx(nmat, 0)] = 1.0-alphamin;
    s[volfracIdx(nmat, 1)] = alphamin;
  }

  auto rb(0.0);
  for (std::size_t k=0; k<nmat; ++k)
  {
    // densities
    r[k] = mat_blk[k].compute< EOS::density >( p, temp );
    // partial density
    s[densityIdx(nmat, k)] = s[volfracIdx(nmat, k)]*r[k];
    // total specific energy
    s[energyIdx(nmat, k)] = s[volfracIdx(nmat, k)]*
      mat_blk[k].compute< EOS::totalenergy >( r[k], u, v, w, p );
    rb += s[densityIdx(nmat, k)];
  }

  s[momentumIdx(nmat, 0)] = rb * u;
  s[momentumIdx(nmat, 1)] = rb * v;
  s[momentumIdx(nmat, 2)] = rb * w;

  return s;
}
