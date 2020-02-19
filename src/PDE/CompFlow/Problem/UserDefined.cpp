// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/UserDefined.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
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

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemUserDefined;

tk::SolutionFn::result_type
CompFlowProblemUserDefined::solution( ncomp_t system,
                                      [[maybe_unused]] ncomp_t ncomp,
                                      [[maybe_unused]] tk::real x,
                                      [[maybe_unused]] tk::real y,
                                      [[maybe_unused]] tk::real z,
                                      [[maybe_unused]] tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \return Values of all components
//! \note The function signature must follow tk::SolutionFn
// *****************************************************************************
{
  Assert( ncomp == ncomp, "Number of scalar components must be " +
                          std::to_string(ncomp) );

  tk::SolutionFn::result_type u( ncomp, 0.0 );

  // Set background ICs
  const auto& ic = g_inputdeck.get< tag::param, eq, tag::ic >();
  const auto& bgdensityic = ic.get< tag::density >();
  const auto& bgvelocityic = ic.get< tag::velocity >();
  const auto& bgpressureic = ic.get< tag::pressure >();
  const auto& bgenergyic = ic.get< tag::energy >();
  const auto& bgtemperatureic = ic.get< tag::temperature >();

  Assert( bgdensityic.size() > system, "No background density IC" );
  Assert( bgvelocityic.size() > 3*system, "No background velocity IC" );

  u[0] = bgdensityic.at(system).at(0);
  u[1] = u[0] * bgvelocityic.at(system).at(0);
  u[2] = u[0] * bgvelocityic.at(system).at(1);
  u[3] = u[0] * bgvelocityic.at(system).at(2);

  if (bgpressureic.size() > system && !bgpressureic[system].empty()) {
    u[4] = eos_totalenergy< eq >( system, u[0], u[1], u[2], u[3],
                                  bgpressureic.at(system).at(0) );
  } else if (bgenergyic.size() > system && !bgenergyic[system].empty()) {
    u[4] = u[0] * bgenergyic[system][0];
  } else
    if (bgtemperatureic.size() > system && !bgtemperatureic[system].empty())
  {
    const auto& cv = g_inputdeck.get< tag::param, eq, tag::cv >();
    u[4] = u[0] * bgtemperatureic[system][0] * cv.at(system).at(0);
  }

  // Apply optional box ICs on top of background ICs
  const auto& icbox = ic.get< tag::box >();
  std::vector< tk::real > box{ icbox.get< tag::xmin >(),
                               icbox.get< tag::xmax >(),
                               icbox.get< tag::ymin >(),
                               icbox.get< tag::ymax >(),
                               icbox.get< tag::zmin >(),
                               icbox.get< tag::zmax >() };
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  if (std::any_of( begin(box), end(box),
        [=]( tk::real p ){ return std::abs(p) > eps; }))
  {
    const auto& boxdensityic = icbox.get< tag::density >();
    const auto& boxvelocityic = icbox.get< tag::velocity >();
    const auto& boxpressureic = icbox.get< tag::pressure >();
    const auto& boxenergyic = icbox.get< tag::energy >();
    const auto& boxtemperatureic = icbox.get< tag::temperature >();

    if (x>box[0] && x<box[1] && y>box[2] && y<box[3] && z>box[4] && z<box[5]) {
      if (boxdensityic.size() > system && !boxdensityic[system].empty()) {
        u[0] = boxdensityic[system][0];
      }
      if (boxvelocityic.size() > system && boxvelocityic[system].size() > 2) {
        u[1] = u[0] * boxvelocityic[system][0];
        u[2] = u[0] * boxvelocityic[system][1];
        u[3] = u[0] * boxvelocityic[system][2];
      }
      if (boxpressureic.size() > system && !boxpressureic[system].empty()) {
        u[4] = eos_totalenergy< eq >( system, u[0], u[1], u[2], u[3],
                                      boxpressureic[system][0] );
      }
      if (boxenergyic.size() > system && !boxenergyic[system].empty()) {
        u[4] = u[0] * boxenergyic[system][0];
      }
      if (boxtemperatureic.size() > system &&
         !boxtemperatureic[system].empty())
      {
        const auto& cv = g_inputdeck.get< tag::param, eq, tag::cv >();
        u[4] = u[0] * boxtemperatureic[system][0] * cv.at(system).at(0);
      }
    }
  }

  return u;
}

tk::SrcFn::result_type
CompFlowProblemUserDefined::src( ncomp_t, ncomp_t, tk::real,
                                 tk::real, tk::real, tk::real )
// *****************************************************************************
//  Compute and return source term for manufactured solution
//! \details No-op for user-defined problems
//! \return Array of reals containing the source for all components
//! \note The function signature must follow tk::SrcFn
// *****************************************************************************
{
  return {{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
}

std::vector< std::string >
CompFlowProblemUserDefined::fieldNames( ncomp_t ) const
// *****************************************************************************
// Return field names to be output to file
//! \return Vector of strings labelling fields output in file
// *****************************************************************************
{
  return CompFlowFieldNames();
}

std::vector< std::vector< tk::real > >
CompFlowProblemUserDefined::fieldOutput(
  ncomp_t system,
  ncomp_t,
  ncomp_t offset,
  tk::real,
  tk::real,
  const std::vector< tk::real >&,
  const std::array< std::vector< tk::real >, 3 >&,
  tk::Fields& U ) const
// *****************************************************************************
//  Return field output going to file
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] offset System offset specifying the position of the system of
//!   PDEs among other systems
//! \param[in] U Solution vector at recent time step
//! \return Vector of vectors to be output to file
// *****************************************************************************
{
  return CompFlowFieldOutput( system, offset, U );
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
