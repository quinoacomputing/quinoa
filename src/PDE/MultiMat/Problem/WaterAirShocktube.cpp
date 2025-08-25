// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/WaterAirShocktube.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the multi-material flow equations
  \details   This file defines a Problem policy class for the multi-material
    compressible flow equations, defined in PDE/MultiMat/MultiMat.hpp. See
    PDE/MultiMat/Problem.hpp for general requirements on Problem policy classes
    for MultiMat.
*/
// *****************************************************************************

#include "WaterAirShocktube.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiMatProblemWaterAirShocktube;

tk::InitializeFn::result_type
MultiMatProblemWaterAirShocktube::initialize( ncomp_t ncomp,
                                              const std::vector< EOS >& mat_blk,
                                              tk::real x,
                                              tk::real,
                                              tk::real,
                                              tk::real )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
//! \details This function only initializes the Water-Air shock tube problem,
//!   but does not actually give the analytical solution at time greater than 0.
//!   The analytical solution would require an exact Riemann solver for
//!   stiffened gas EoS, which has not been implemented yet.
// *****************************************************************************
{
  // see also Control/Inciter/InputDeck/Grammar.hpp
  Assert( ncomp == 9, "Number of scalar components must be 9" );

  auto nmat = g_inputdeck.get< eq, tag::nmat >();
  auto alphamin = g_inputdeck.get< eq, tag::min_volumefrac >();

  std::vector< tk::real > s(ncomp, 0.0), r(nmat, 0.0);
  tk::real p, u, v, w;

  if (x<0.75) {
    // volume-fraction
    s[volfracIdx(nmat, 0)] = 1.0-alphamin;
    s[volfracIdx(nmat, 1)] = alphamin;
    // pressure
    p = 1.0e9;
    // densities
    for (std::size_t k=0; k<nmat; ++k)
      r[k] = mat_blk[k].compute< EOS::density >( p, 494.646 );
    // velocity
    u = 0.0;
    v = 0.0;
    w = 0.0;
  }
  else {
    // volume-fraction
    s[volfracIdx(nmat, 0)] = alphamin;
    s[volfracIdx(nmat, 1)] = 1.0-alphamin;
    // pressure
    p = 1.0e5;
    // densities
    for (std::size_t k=0; k<nmat; ++k)
      r[k] = mat_blk[k].compute< EOS::density >( p, 34.844 );
    // velocity
    u = 0.0;
    v = 0.0;
    w = 0.0;
  }
  for (std::size_t k=0; k<nmat; ++k)
  {
    // partial density
    s[densityIdx(nmat, k)] = s[volfracIdx(nmat, k)]*r[k];
    // total specific energy
    s[energyIdx(nmat, k)] =
      mat_blk[k].compute< EOS::totalenergy >( s[volfracIdx(nmat, k)]*r[k],
      u, v, w, s[volfracIdx(nmat, k)]*p, s[volfracIdx(nmat, k)] );
  }

  return s;
}
