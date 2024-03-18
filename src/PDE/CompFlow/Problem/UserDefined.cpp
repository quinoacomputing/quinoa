// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/UserDefined.cpp
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

#include <limits>

#include "UserDefined.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"
#include "EoS/GetMatProp.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemUserDefined;

tk::InitializeFn::result_type
CompFlowProblemUserDefined::initialize( ncomp_t ncomp,
                                        const std::vector< EOS >& mat_blk,
                                        tk::real,
                                        tk::real,
                                        tk::real,
                                        tk::real )
// *****************************************************************************
//! Set initial conditions
//! \param[in] ncomp Number of scalar components in this PDE system
//! \return Values of all components
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  tk::InitializeFn::result_type u( ncomp, 0.0 );

  // Set background ICs
  const auto& ic = g_inputdeck.get< tag::param, eq, tag::ic >();
  const auto& bgrhoic = ic.get< tag::density >();
  const auto& bgvelic = ic.get< tag::velocity >();
  const auto& bgpreic = ic.get< tag::pressure >();
  const auto& bgenic = ic.get< tag::energy >();
  const auto& bgtempic = ic.get< tag::temperature >();

  u[0] = bgrhoic.at(0);
  u[1] = u[0] * bgvelic.at(0);
  u[2] = u[0] * bgvelic.at(1);
  u[3] = u[0] * bgvelic.at(2);

  if (!bgpreic.empty()) {
    u[4] = mat_blk[0].compute< EOS::totalenergy >( u[0], u[1]/u[0], u[2]/u[0],
      u[3]/u[0], bgpreic.at(0) );
  } else if (!bgenic.empty()) {
    u[4] = u[0] * bgenic[0];
  } else if (!bgtempic.empty()) {
    const auto& c_v = getmatprop< newtag::cv >();
    u[4] = u[0] * bgtempic[0] * c_v;
  } else Throw( "IC background energy cannot be computed. User must specify "
                "one of background pressure, energy, or temperature." );

  return u;
}

std::vector< std::string >
CompFlowProblemUserDefined::names( ncomp_t ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
