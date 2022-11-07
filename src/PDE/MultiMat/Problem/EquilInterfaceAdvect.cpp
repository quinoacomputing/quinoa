// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/EquilInterfaceAdvect.cpp
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

#include "EquilInterfaceAdvect.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::MultiMatProblemEquilInterfaceAdvect;

tk::InitializeFn::result_type
MultiMatProblemEquilInterfaceAdvect::initialize( ncomp_t system,
                                                 ncomp_t ncomp,
                                              const std::vector< EOS >& mat_blk,
                                                 tk::real x,
                                                 tk::real y,
                                                 tk::real z,
                                                 tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time at which to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
//! \note The function signature must follow tk::InitializeFn
//! \details This function initializes the equilibrium interface advection
//!   problem, and gives the analytical solution at time greater than 0.
// *****************************************************************************
{
  // see also Control/Inciter/InputDeck/Grammar.hpp
  Assert( ncomp == 9, "Number of scalar components must be 9" );

  auto nmat = g_inputdeck.get< tag::param, eq, tag::nmat >()[system];

  std::vector< tk::real > s( ncomp, 0.0 );
  tk::real r, p, u, v, w;
  auto alphamin = 1.0e-12;

  // pressure
  p = 0.4;
  // velocity
  u = 3.0;
  v = 2.0;
  w = 1.0;
  // location of interface
  auto xc = 0.45 + u*t;
  auto yc = 0.45 + v*t;
  auto zc = 0.45 + w*t;
  // volume-fraction
  s[volfracIdx(nmat, 0)] = (1.0-2.0*alphamin) * 0.5 *
    (1.0-std::tanh(5.0*((x-xc)+(y-yc)+(z-zc)))) + alphamin;
  s[volfracIdx(nmat, 1)] = 1.0 - s[volfracIdx(nmat, 0)];
  // density
  r = 5.0 + x + y + z;
  s[densityIdx(nmat, 0)] = s[volfracIdx(nmat, 0)]*r;
  s[densityIdx(nmat, 1)] = s[volfracIdx(nmat, 1)]*r;
  // total specific energy
  s[energyIdx(nmat, 0)] = s[volfracIdx(nmat, 0)]*
    mat_blk[0].eosCall< EOS::totalenergy >( r, u, v, w, p );
  s[energyIdx(nmat, 1)] = s[volfracIdx(nmat, 1)]*
    mat_blk[1].eosCall< EOS::totalenergy >( r, u, v, w, p );
  // momentum
  s[momentumIdx(nmat, 0)] = r*u;
  s[momentumIdx(nmat, 1)] = r*v;
  s[momentumIdx(nmat, 2)] = r*w;

  return s;
}
