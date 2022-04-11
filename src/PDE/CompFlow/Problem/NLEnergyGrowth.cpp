// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/NLEnergyGrowth.cpp
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

#include "NLEnergyGrowth.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemNLEnergyGrowth;

tk::real
CompFlowProblemNLEnergyGrowth::hx( tk::real bx, tk::real by, tk::real bz,
                                   tk::real x, tk::real y, tk::real z )
// *****************************************************************************
//  Compute internal energy parameter
//! \param[in] bx Parameter betax
//! \param[in] by Parameter betay
//! \param[in] bz Parameter betaz
//! \param[in] x X coordinate to evaluate at
//! \param[in] y Y coordinate to evaluate at
//! \param[in] z Z coordinate to evaluate at
//! \return Internal energy parameter
// *****************************************************************************
{
  return std::cos(bx*M_PI*x) * std::cos(by*M_PI*y) * std::cos(bz*M_PI*z);
}

tk::real
CompFlowProblemNLEnergyGrowth::ec( tk::real ce, tk::real kappa, tk::real t,
                                   tk::real h, tk::real p )
// *****************************************************************************
//  Compute a power of the internal energy
//! \param[in] ce Internal energy parameter
//! \param[in] kappa Internal energy parameter
//! \param[in] t Physical time
//! \param[in] h Internal energy parameter
//! \param[in] p Power
//! \return Internal energy raised to power p
// *****************************************************************************
{
  return std::pow( -3.0*(ce + kappa*h*h*t), p );
}

tk::InitializeFn::result_type
CompFlowProblemNLEnergyGrowth::initialize( ncomp_t system,
                                           ncomp_t,
                                           tk::real x,
                                           tk::real y,
                                           tk::real z,
                                           tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  using tag::param;

  // manufactured solution parameters
  auto ce = g_inputdeck.get< param, eq, tag::ce >()[system];
  auto r0 = g_inputdeck.get< param, eq, tag::r0 >()[system];
  auto a = g_inputdeck.get< param, eq, tag::alpha >()[system];
  auto k = g_inputdeck.get< param, eq, tag::kappa >()[system];
  auto bx = g_inputdeck.get< param, eq, tag::betax >()[system];
  auto by = g_inputdeck.get< param, eq, tag::betay >()[system];
  auto bz = g_inputdeck.get< param, eq, tag::betaz >()[system];
  // spatial component of density field
  auto gx = 1.0 - x*x - y*y - z*z;
  // internal energy parameter
  auto h = hx( bx, by, bz, x, y, z );
  // temporal component of the density field
  tk::real ft = std::exp( -a*t );
  // density
  auto r = r0 + ft*gx;
  // energy
  auto re = r*ec(ce,k,t,h,-1.0/3.0);

  return {{ r, 0.0, 0.0, 0.0, re }};
}

tk::InitializeFn::result_type
CompFlowProblemNLEnergyGrowth::analyticSolution( ncomp_t system,
                                                 ncomp_t,
                                               std::vector< EoS_Base* > mat_blk,
                                                 tk::real x,
                                                 tk::real y,
                                                 tk::real z,
                                                 tk::real t )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \param[in] t Time where to evaluate the solution
//! \return Values of all components evaluated at (x,y,z,t)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  using tag::param;

  // manufactured solution parameters
  auto ce = g_inputdeck.get< param, eq, tag::ce >()[system];
  auto r0 = g_inputdeck.get< param, eq, tag::r0 >()[system];
  auto a = g_inputdeck.get< param, eq, tag::alpha >()[system];
  auto k = g_inputdeck.get< param, eq, tag::kappa >()[system];
  auto bx = g_inputdeck.get< param, eq, tag::betax >()[system];
  auto by = g_inputdeck.get< param, eq, tag::betay >()[system];
  auto bz = g_inputdeck.get< param, eq, tag::betaz >()[system];
  // spatial component of density field
  auto gx = 1.0 - x*x - y*y - z*z;
  // internal energy parameter
  auto h = hx( bx, by, bz, x, y, z );
  // temporal component of the density field
  tk::real ft = std::exp( -a*t );
  // density
  auto r = r0 + ft*gx;
  // energy
  auto re = r*ec(ce,k,t,h,-1.0/3.0);
  // pressure
  auto p = mat_blk[0]->eos_pressure( system, r, 0.0, 0.0, 0.0, re );

  return {{ r, 0.0, 0.0, 0.0, re/r, p }};
}

std::vector< std::string >
CompFlowProblemNLEnergyGrowth::analyticFieldNames( ncomp_t ) const
// *****************************************************************************
// Return analytic field names to be output to file
//! \return Vector of strings labelling fields output in file
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
CompFlowProblemNLEnergyGrowth::names( ncomp_t ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
