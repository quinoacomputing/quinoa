// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/TaylorGreen.cpp
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

#include "TaylorGreen.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FieldOutput.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemTaylorGreen;

tk::InitializeFn::result_type
CompFlowProblemTaylorGreen::initialize( ncomp_t system,
                                        ncomp_t,
                                        tk::real x,
                                        tk::real y,
                                        tk::real,
                                        tk::real )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  using tag::param; using std::sin; using std::cos;

  // density
  auto r = 1.0;
  // pressure
  auto p = 10.0 + r/4.0*(cos(2.0*M_PI*x) + cos(2.0*M_PI*y));
  // velocity
  auto u =  sin(M_PI*x) * cos(M_PI*y);
  auto v = -cos(M_PI*x) * sin(M_PI*y);
  auto w = 0.0;
  // total specific energy
  auto rE = eos_totalenergy< eq >( system, r, u, v, w, p );

  return {{ r, r*u, r*v, r*w, rE }};
}

tk::InitializeFn::result_type
CompFlowProblemTaylorGreen::analyticSolution( ncomp_t system,
                                              ncomp_t,
                                              tk::real x,
                                              tk::real y,
                                              tk::real,
                                              tk::real )
// *****************************************************************************
//! Evaluate analytical solution at (x,y,z,t) for all components
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] x X coordinate where to evaluate the solution
//! \param[in] y Y coordinate where to evaluate the solution
//! \return Values of all components evaluated at (x)
//! \note The function signature must follow tk::InitializeFn
// *****************************************************************************
{
  using tag::param; using std::sin; using std::cos;

  // density
  auto r = 1.0;
  // pressure
  auto p = 10.0 + r/4.0*(cos(2.0*M_PI*x) + cos(2.0*M_PI*y));
  // velocity
  auto u =  sin(M_PI*x) * cos(M_PI*y);
  auto v = -cos(M_PI*x) * sin(M_PI*y);
  auto w = 0.0;
  // total specific energy
  auto E = eos_totalenergy< eq >( system, r, u, v, w, p ) / r;

  return {{ r, u, v, w, E }};
}

std::vector< std::string >
CompFlowProblemTaylorGreen::fieldNames( ncomp_t ) const
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
CompFlowProblemTaylorGreen::fieldOutput(
  ncomp_t system,
  ncomp_t,
  ncomp_t offset,
  std::size_t nunk,
  std::size_t rdof,
  tk::real,
  tk::real V,
  const std::vector< tk::real >& vol,
  const std::array< std::vector< tk::real >, 3 >& coord,
  const tk::Fields& U ) const
// *****************************************************************************
//  Return field output going to file
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] offset System offset specifying the position of the system of
//!   PDEs among other systems
//! \param[in] nunk Number of unknowns to extract
//! \param[in] rdof Number of reconstructed degrees of freedom. This is used as
//!   the number of scalar components to shift when extracting scalar
//!   components.
//! \param[in] V Total mesh volume (across the whole problem)
//! \param[in] vol Nodal mesh volumes
//! \param[in] coord Mesh point coordinates
//! \param[in] U Solution vector at recent time step
//! \return Vector of vectors to be output to file
// *****************************************************************************
{
  auto out = CompFlowFieldOutput( system, offset, nunk, rdof, U );

  const auto r  = U.extract( 0*rdof, offset );
  const auto ru = U.extract( 1*rdof, offset );
  const auto rv = U.extract( 2*rdof, offset );
  const auto re = U.extract( 4*rdof, offset );

  // mesh node coordinates
  const auto& x = coord[0];
  const auto& y = coord[1];

  out.push_back( std::vector< tk::real >( nunk, 1.0 ) );

  std::vector< tk::real > ua( nunk );
  for (std::size_t i=0; i<nunk; ++i)
    ua[i] = std::sin(M_PI*x[i]) * std::cos(M_PI*y[i]);
  out.push_back( ua );

  // error in x-velocity
  std::vector< tk::real > err( nunk );
  for (std::size_t i=0; i<nunk; ++i)
    err[i] = std::pow( ua[i] - ru[i], 2.0 ) * vol[i] / V;
  out.push_back( err );

  std::vector< tk::real > va( nunk );
  for (std::size_t i=0; i<nunk; ++i)
    va[i] = -std::cos(M_PI*x[i]) * std::sin(M_PI*y[i]);
  out.push_back( va );

  // error in v-velocity
  for (std::size_t i=0; i<nunk; ++i)
    err[i] = std::pow( va[i] - rv[i], 2.0 ) * vol[i] / V;
  out.push_back( err );

  out.push_back( std::vector< tk::real >( nunk, 0.0 ) );

  std::vector< tk::real > Ea( nunk ), Pa( nunk );
  for (std::size_t i=0; i<nunk; ++i) {
    Pa[i] = 10.0 + r[i]/4.0*(std::cos(2.0*M_PI*x[i]) + std::cos(2.0*M_PI*y[i]));
    Ea[i] = eos_totalenergy< eq >( system, r[i], ua[i]/r[i], va[i]/r[i],
                                   0.0, Pa[i]/r[i] );
  }
  out.push_back( Ea );

  // error in total specific energy
  for (std::size_t i=0; i<nunk; ++i)
    err[i] = std::pow( Ea[i] - re[i], 2.0 ) * vol[i] / V;
  out.push_back( err );

  out.push_back( Pa );

  return out;
}

std::vector< std::string >
CompFlowProblemTaylorGreen::names( ncomp_t ) const
// *****************************************************************************
//  Return names of integral variables to be output to diagnostics file
//! \return Vector of strings labelling integral variables output
// *****************************************************************************
{
  return { "r", "ru", "rv", "rw", "re" };
}
