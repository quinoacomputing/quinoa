// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/SodShocktube.cpp
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

#include "SodShocktube.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiMatProblemSodShocktube;

tk::InitializeFn::result_type
MultiMatProblemSodShocktube::initialize( ncomp_t ncomp,
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
//! \details This function only initializes the Sod shock tube problem, but does
//!   not actually give the analytical solution at time greater than 0. The
//!   analytical solution would require an exact Riemann solver, which has not
//!   been implemented yet.
// *****************************************************************************
{
  // see also Control/Inciter/InputDeck/Grammar.hpp
  Assert( ncomp == 9, "Number of scalar components must be 9" );

  auto nmat = g_inputdeck.get< tag::param, eq, tag::nmat >();

  std::vector< tk::real > s( ncomp, 0.0 );
  tk::real r, p, u, v, w;
  auto alphamin = 1.0e-12;

  if (x<0.5) {
    // volume-fraction
    s[volfracIdx(nmat, 0)] = 1.0-alphamin;
    s[volfracIdx(nmat, 1)] = alphamin;
    // density
    r = 1.0;
    // pressure
    p = 1.0;
    // velocity
    u = 0.0;
    v = 0.0;
    w = 0.0;
  }
  else {
    // volume-fraction
    s[volfracIdx(nmat, 0)] = alphamin;
    s[volfracIdx(nmat, 1)] = 1.0-alphamin;
    // density
    r = 0.125;
    // pressure
    p = 0.1;
    // velocity
    u = 0.0;
    v = 0.0;
    w = 0.0;
  }
  s[densityIdx(nmat, 0)] = s[volfracIdx(nmat, 0)]*r;
  s[densityIdx(nmat, 1)] = s[volfracIdx(nmat, 1)]*r;
  // total specific energy
  s[energyIdx(nmat, 0)] = s[volfracIdx(nmat, 0)]*
    mat_blk[0].compute< EOS::totalenergy >( r, u, v, w, p );
  s[energyIdx(nmat, 1)] = s[volfracIdx(nmat, 1)]*
    mat_blk[1].compute< EOS::totalenergy >( r, u, v, w, p );

  return s;
}
