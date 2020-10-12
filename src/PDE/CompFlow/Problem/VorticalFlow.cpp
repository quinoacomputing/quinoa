// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/VorticalFlow.cpp
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

#include "VorticalFlow.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemVorticalFlow;

tk::SolutionFn::result_type
CompFlowProblemVorticalFlow::solution( ncomp_t system,
                                       [[maybe_unused]] ncomp_t ncomp,
                                       tk::real x,
                                       tk::real y,
                                       tk::real z,
                                       tk::real,
                                       int& )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \param[in] z Z coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::SolutionFn
// *****************************************************************************
{
  Assert( ncomp == ncomp, "Number of scalar components must be " +
                          std::to_string(ncomp) );
  using tag::param; using tag::compflow;

  // manufactured solution parameters
  const auto& a =
    g_inputdeck.get< param, compflow, tag::alpha >()[ system ];
  const auto& b = g_inputdeck.get< param, compflow, tag::beta >()[ system ];
  const auto& p0 = g_inputdeck.get< param, compflow, tag::p0 >()[ system ];
  // ratio of specific heats
  tk::real g = g_inputdeck.get< param, compflow, tag::gamma >()[ system ][0];
  // velocity
  const tk::real ru = a*x - b*y;
  const tk::real rv = b*x + a*y;
  const tk::real rw = -2.0*a*z;
  // total specific energy
  const tk::real rE = (ru*ru+rv*rv+rw*rw)/2.0 + (p0-2.0*a*a*z*z)/(g-1.0);

  return {{ 1.0, ru, rv, rw, rE }};
}

std::vector< std::string >
CompFlowProblemVorticalFlow::fieldNames( ncomp_t ) const
// *****************************************************************************
// Return field names to be output to file
//! \return Vector of strings labelling fields output in file
// *****************************************************************************
{
  auto n = CompFlowFieldNames();

  n.push_back( "density_analytical" );
  n.push_back( "x-velocity_analytical" );
  n.push_back( "y-velocity_analytical" );
  n.push_back( "z-velocity_analytical" );
  n.push_back( "specific_total_energy_analytical" );
  n.push_back( "pressure_analytical" );

  return n;
}

std::vector< std::vector< tk::real > >
CompFlowProblemVorticalFlow::fieldOutput(
  ncomp_t system,
  ncomp_t,
  ncomp_t offset,
  std::size_t nunk,
  std::size_t rdof,
  tk::real,
  tk::real,
  const std::vector< tk::real >&,
  const std::array< std::vector< tk::real >, 3 >& coord,
  const tk::Fields& U ) const
// *****************************************************************************
//  Return field output going to file
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] offset System offset specifying the position of the system of
//!   PDEs among other systems
//! \param[in] nunk Number of unknowns to extract
//! \param[in] coord Mesh node coordinates
//! \param[in] U Solution vector at recent time step
//! \return Vector of vectors to be output to file
// *****************************************************************************
{
   // manufactured solution parameters
   const auto& a =
     g_inputdeck.get< tag::param, tag::compflow, tag::alpha >()[system];
   const auto& b =
     g_inputdeck.get< tag::param, tag::compflow, tag::beta >()[system];
   const auto& p0 =
     g_inputdeck.get< tag::param, tag::compflow, tag::p0 >()[system];
   // ratio of specific heats
   tk::real g =
     g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[system][0];

   auto out = CompFlowFieldOutput( system, offset, nunk, rdof, U );

   const auto r  = U.extract( 0*rdof, offset );
   const auto ru = U.extract( 1*rdof, offset );
   const auto rv = U.extract( 2*rdof, offset );
   const auto rw = U.extract( 3*rdof, offset );
   const auto re = U.extract( 4*rdof, offset );

   // mesh node coordinates
   const auto& x = coord[0];
   const auto& y = coord[1];
   const auto& z = coord[2];

   out.push_back( std::vector< tk::real >( nunk, 1.0 ) );

   std::vector< tk::real > u = ru;
   for (std::size_t i=0; i<nunk; ++i) u[i] = a*x[i] - b*y[i];
   out.push_back( u );

   std::vector< tk::real > v = rv;
   for (std::size_t i=0; i<nunk; ++i) v[i] = b*x[i] + a*y[i];
   out.push_back( v );

   std::vector< tk::real > w = rw;
   for (std::size_t i=0; i<nunk; ++i) w[i] = -2.0*a*z[i];
   out.push_back( w );

   std::vector< tk::real > E = re;
   for (std::size_t i=0; i<nunk; ++i)
     E[i] = 0.5*(u[i]*u[i] + v[i]*v[i] + w[i]*w[i]) +
            (p0 - 2.0*a*a*z[i]*z[i])/(g-1.0);
   out.push_back( E );

   std::vector< tk::real > P( nunk, 0.0 );
   for (std::size_t i=0; i<nunk; ++i)
     P[i] = p0 - 2.0*a*a*z[i]*z[i];
   out.push_back( P );

   return out;
}

std::vector< std::string >
CompFlowProblemVorticalFlow::names( ncomp_t ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
