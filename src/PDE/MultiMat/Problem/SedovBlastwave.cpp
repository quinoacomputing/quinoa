// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/SedovBlastwave.cpp
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

#include "SedovBlastwave.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiMatProblemSedovBlastwave;

tk::InitializeFn::result_type
MultiMatProblemSedovBlastwave::initialize( ncomp_t ncomp,
                                           const std::vector< EOS >& mat_blk,
                                           tk::real x,
                                           tk::real y,
                                           tk::real z,
                                           tk::real )
// *****************************************************************************
//  Set initial conditions for the Sedov blast-wave problem
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x,y)
//! \note The function signature must follow tk::InitializeFn
//! \details This function only configures for nmat=1, cloning the Sedov
//!   blast-wave setup from CompFlowProblemSedovBlastwave.
// *****************************************************************************
{
  auto nmat = g_inputdeck.get< eq, tag::nmat >();

  Assert( nmat == 1, "Sedov blast-wave problem configured for nmat=1 only" );
  Assert( ncomp == 3*nmat+3, "Number of scalar components must be 6" );

  std::vector< tk::real > s( ncomp, 0.0 );

  // pressure
  tk::real p = ( x < 0.05 && y < 0.05 && z < 0.05 ) ? 783.4112 : 1.0e-6;
  // density
  tk::real r = 1.0;
  // velocity
  tk::real u = 0.0, v = 0.0, w = 0.0;

  // volume fraction (single material)
  s[volfracIdx(nmat, 0)] = 1.0;
  // partial density
  s[densityIdx(nmat, 0)] = r;
  // total specific energy
  s[energyIdx(nmat, 0)] =
    mat_blk[0].compute< EOS::totalenergy >( r, u, v, w, p, 1.0 );

  // momentum
  s[momentumIdx(nmat, 0)] = r * u;
  s[momentumIdx(nmat, 1)] = r * v;
  s[momentumIdx(nmat, 2)] = r * w;

  return s;
}
