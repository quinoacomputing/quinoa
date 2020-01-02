// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/UserDefined.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Problem configuration for the compressible flow equations
  \details   This file defines a Problem policy class for the compressible flow
    equations, defined in PDE/CompFlow/CompFlow.h. See PDE/CompFlow/Problem.h
    for general requirements on Problem policy classes for CompFlow.
*/
// *****************************************************************************

#include "UserDefined.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemUserDefined;

tk::SolutionFn::result_type
CompFlowProblemUserDefined::solution( ncomp_t,
                                      [[maybe_unused]] ncomp_t ncomp,
                                      tk::real,
                                      tk::real,
                                      tk::real,
                                      tk::real )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] ncomp Number of scalar components in this PDE system
//! \return Values of all components
//! \note The function signature must follow tk::SolutionFn
// *****************************************************************************
{
  Assert( ncomp == ncomp, "Number of scalar components must be " +
                          std::to_string(ncomp) );
  return {{ 1.0, 0.0, 0.0, 1.0, 293.0 }};
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
  auto n = CompFlowFieldNames();
  n.push_back( "temperature" );

  return n;
}

std::vector< std::vector< tk::real > >
CompFlowProblemUserDefined::fieldOutput(
  ncomp_t,
  ncomp_t,
  ncomp_t offset,
  tk::real,
  tk::real,
  const std::vector< tk::real >&,
  const std::array< std::vector< tk::real >, 3 >&,
  tk::Fields& U ) const
// *****************************************************************************
//  Return field output going to file
//! \param[in] offset System offset specifying the position of the system of
//!   PDEs among other systems
//! \param[in] U Solution vector at recent time step
//! \return Vector of vectors to be output to file
// *****************************************************************************
{
  // number of degrees of freedom
  const std::size_t rdof =
    g_inputdeck.get< tag::discr, tag::rdof >();

  auto out = CompFlowFieldOutput(0, offset, U);

  const auto r  = U.extract( 0*rdof, offset );
  const auto ru = U.extract( 1*rdof, offset );
  const auto rv = U.extract( 2*rdof, offset );
  const auto rw = U.extract( 3*rdof, offset );
  const auto re = U.extract( 4*rdof, offset );

  std::vector< tk::real > p = r;
  for (std::size_t i=0; i<p.size(); ++i)
    p[i] = eos_pressure< eq >( 0, r[i], ru[i], rv[i], rw[i], re[i] );
  out.push_back( p );

  std::vector< tk::real > T = r;
  tk::real cv = g_inputdeck.get< tag::param, eq, tag::cv >()[0][0];
  for (std::size_t i=0; i<T.size(); ++i)
    T[i] = cv*(re[i] - (ru[i]*ru[i] + rv[i]*rv[i] + rw[i]*rw[i])/2.0);
  out.push_back( T );

  return out;
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
