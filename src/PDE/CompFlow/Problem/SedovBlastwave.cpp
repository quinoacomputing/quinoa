// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/SedovBlastwave.cpp
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

#include "SedovBlastwave.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemSedovBlastwave;

tk::InitializeFn::result_type
CompFlowProblemSedovBlastwave::initialize( ncomp_t system,
                                           ncomp_t,
                                           tk::real x,
                                           tk::real y,
                                           tk::real z,
                                           tk::real )
// *****************************************************************************
//  Initialize numerical solution
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  tk::real r=0, p=0, u=0, v=0, w=0, rE=0;

  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
  const auto centering = ctr::Scheme().centering( scheme );

  // pressure
  if (centering == tk::Centering::ELEM) {

    if ( (x<0.05) && (y<0.05) ) p = 783.4112; else p = 1.0e-6;

  } else if (centering == tk::Centering::NODE) {

    auto eps = std::numeric_limits< tk::real >::epsilon();
    if (std::abs(x) < eps && std::abs(y) < eps && std::abs(z) < eps)
      p = g_inputdeck.get< tag::param, tag::compflow, tag::p0 >()[ system ];
    else
      p = 0.67e-4;

  }

  // density
  r = 1.0;
  // velocity
  u = v = w = 0.0;
  // total specific energy
  rE = eos_totalenergy< eq >( system, r, u, v, w, p );

  return {{ r, r*u, r*v, r*w, rE }};
}

tk::InitializeFn::result_type
CompFlowProblemSedovBlastwave::analyticSolution( ncomp_t system,
                                                 ncomp_t,
                                                 std::vector< EoS_Base* >,
                                                 tk::real x,
                                                 tk::real y,
                                                 tk::real z,
                                                 tk::real )
// *****************************************************************************
//  Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
//! \warning This is NOT the analytic solution at all times, only at t=0
// *****************************************************************************
{
  return initialize( system, 0, x, y, z, 0 );
}

std::vector< std::string >
CompFlowProblemSedovBlastwave::analyticFieldNames( ncomp_t ) const
// *****************************************************************************
// Return analytic field names to be output to file
//! \return Vector of strings labelling analytic fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;
  n.push_back( "density_analytical" );
  n.push_back( "x-velocity_analytical" );
  n.push_back( "y-velocity_analytical" );
  n.push_back( "z-velocity_analytical" );
  n.push_back( "specific_total_energy_analytical" );
  n.push_back( "pressure_analytical" );

  return n;
}
std::vector< std::string >
CompFlowProblemSedovBlastwave::names( ncomp_t ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
