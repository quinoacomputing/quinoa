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

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // ::inciter

using inciter::CompFlowProblemUserDefined;

tk::SolutionFn::result_type
CompFlowProblemUserDefined::solution( ncomp_t,
                                      ncomp_t ncomp,
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
  Assert( ncomp == m_ncomp, "Number of scalar components must be " +
                            std::to_string(m_ncomp) );
  IGNORE(ncomp);

  return {{ 1.0, 0.0, 0.0, 1.0, 293.0 }};
}

std::array< tk::real, 5 >
CompFlowProblemUserDefined::solinc( ncomp_t, tk::real, tk::real,
                                    tk::real, tk::real, tk::real ) const
// *****************************************************************************
// Evaluate the increment from t to t+dt of the analytical solution at (x,y,z)
// for all components
//! \return Increment in values of all components evaluated at (x,y,z,t+dt)
// *****************************************************************************
{
  return {{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
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

void
CompFlowProblemUserDefined::side( std::unordered_set< int >& conf ) const
// *****************************************************************************
//  Query all side set IDs the user has configured for all components in this
//  PDE system
//! \param[in,out] conf Set of unique side set IDs to add to
// *****************************************************************************
{
  using tag::param; using tag::bcdir;

  for (const auto& s : g_inputdeck.get< param, eq, bcdir >())
    conf.insert( std::stoi(s[0]) );
}

std::vector< std::string >
CompFlowProblemUserDefined::fieldNames( ncomp_t ) const
// *****************************************************************************
// Return field names to be output to file
//! \return Vector of strings labelling fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;

  n.push_back( "density" );
  n.push_back( "x-velocity" );
  n.push_back( "y-velocity" );
  n.push_back( "z-velocity" );
  n.push_back( "specific total energy" );
  n.push_back( "pressure" );
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
  std::vector< std::vector< tk::real > > out;

  const auto r = U.extract( 0, offset );
  const auto ru = U.extract( 1, offset );
  const auto rv = U.extract( 2, offset );
  const auto rw = U.extract( 3, offset );
  const auto re = U.extract( 4, offset );

  out.push_back( r );

  std::vector< tk::real > u = ru;
  std::transform( r.begin(), r.end(), u.begin(), u.begin(),
                  []( tk::real s, tk::real& d ){ return d /= s; } );
  out.push_back( u );

  std::vector< tk::real > v = rv;
  std::transform( r.begin(), r.end(), v.begin(), v.begin(),
                  []( tk::real s, tk::real& d ){ return d /= s; } );
  out.push_back( v );

  std::vector< tk::real > w = rw;
  std::transform( r.begin(), r.end(), w.begin(), w.begin(),
                  []( tk::real s, tk::real& d ){ return d /= s; } );
  out.push_back( w );

  std::vector< tk::real > E = re;
  std::transform( r.begin(), r.end(), E.begin(), E.begin(),
                  []( tk::real s, tk::real& d ){ return d /= s; } );
  out.push_back( E );

  std::vector< tk::real > p = r;
  tk::real g = g_inputdeck.get< tag::param, eq, tag::gamma >()[0];
  for (std::size_t i=0; i<p.size(); ++i)
    p[i] = (g-1.0)*r[i]*(E[i] - (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0);
  out.push_back( p );

  std::vector< tk::real > T = r;
  tk::real cv = g_inputdeck.get< tag::param, eq, tag::cv >()[0];
  for (std::size_t i=0; i<T.size(); ++i)
    T[i] = cv*(E[i] - (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0);
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
