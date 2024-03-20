// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/SheddingFlow.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a Problem policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problem.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************

#include "SheddingFlow.hpp"
#include "Inciter/InputDeck/New2InputDeck.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::New2InputDeck g_newinputdeck;

} // ::inciter

using inciter::CompFlowProblemSheddingFlow;

tk::InitializeFn::result_type
CompFlowProblemSheddingFlow::initialize( ncomp_t,
                                         const std::vector< EOS >& mat_blk,
                                         tk::real,
                                         tk::real,
                                         tk::real,
                                         tk::real )
// *****************************************************************************
//! Evaluate initial solution at (x,y,z,t) for all components
//! \param[in] x X coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  // Assign uniform initial condition according to the farfield state
  auto r = g_newinputdeck.get< newtag::bc >()[0].get< newtag::density >();
  auto p = g_newinputdeck.get< newtag::bc >()[0].get< newtag::pressure >();
  const auto& u = g_newinputdeck.get< newtag::bc >()[0].get< newtag::velocity >();
  auto rE = mat_blk[0].compute< EOS::totalenergy >( r, u[0], u[1], u[2], p );

  return {{ r, r*u[0], r*u[1], r*u[2], rE }};
}

std::vector< std::string >
CompFlowProblemSheddingFlow::names( ncomp_t ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
