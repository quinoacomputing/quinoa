// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/ShockDensityWave.cpp
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

#include "ShockDensityWave.hpp"
#include "Inciter/InputDeck/New2InputDeck.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::New2InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiMatProblemShockDensityWave;

tk::InitializeFn::result_type
MultiMatProblemShockDensityWave::initialize( ncomp_t ncomp,
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
//! \details This function only initializes the Shock-density wave problem, but
//!   does not actually give the analytical solution at time greater than 0.
//!   This problem does not have an analytical solution.
// *****************************************************************************
{
  auto nmat = g_inputdeck.get< eq, newtag::nmat >();

  // see also Control/Inciter/InputDeck/Grammar.hpp
  Assert( ncomp == 3*nmat+3, "Number of scalar components must be 6 or 9" );

  std::vector< tk::real > s( ncomp, 0.0 );
  tk::real r, p, u, v, w;
  auto alphamin = 1.0e-12;

  if(nmat > 1) {                  // If this is multi-material test
    if(x > -4.0) {
      s[volfracIdx(nmat, 0)] = 1.0-alphamin;
      s[volfracIdx(nmat, 1)] = alphamin;
    } else {
      s[volfracIdx(nmat, 0)] = alphamin;
      s[volfracIdx(nmat, 1)] = 1.0-alphamin;
    }
  } else if(nmat == 1){           // If this is a single-material test
      s[volfracIdx(nmat, 0)] = 1.0;
  } else {
    Throw("The test is not configured for nmat > 2");
  }

  if (x > -4.0) {
    // density
    r = 1.0 + 0.2 * sin(5.0 * x);
    // pressure
    p = 1.0;
    // velocity
    u = 0.0;
    v = 0.0;
    w = 0.0;
  }
  else {
    // density
    r = 3.8571;
    // pressure
    p = 10.3333;
    // velocity
    u = 2.6294;
    v = 0.0;
    w = 0.0;
  }

  // bulk density
  tk::real rb(0.0);

  for(std::size_t imat = 0; imat < nmat; imat++) {
    // partial density
    s[densityIdx(nmat, imat)] = s[volfracIdx(nmat, imat)]*r;

    // total specific energy
    s[energyIdx(nmat, imat)] = s[volfracIdx(nmat, imat)]*
      mat_blk[imat].compute< EOS::totalenergy >( r, u, v, w, p );

    rb += s[densityIdx(nmat, imat)];
  }

  // momentum
  s[momentumIdx(nmat, 0)] = rb * u;
  s[momentumIdx(nmat, 1)] = rb * v;
  s[momentumIdx(nmat, 2)] = rb * w;

  return s;
}
