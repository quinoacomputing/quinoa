// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/FieldOutput.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field outputs for single-material equation solver
  \details   This file defines functions for field quantites to be output to
    files for compressible single-material equations.
*/
// *****************************************************************************

#include "FieldOutput.hpp"
#include "EoS/EoS.hpp"
#include "ContainerUtil.hpp"

namespace inciter {

std::vector< std::string > CompFlowFieldNames()
// *****************************************************************************
// Return field names to be output to file
//! \return Vector of strings labelling fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;

  n.push_back( "density_numerical" );
  n.push_back( "x-velocity_numerical" );
  n.push_back( "y-velocity_numerical" );
  n.push_back( "z-velocity_numerical" );
  n.push_back( "specific_total_energy_numerical" );
  n.push_back( "pressure_numerical" );

  return n;
}

std::vector< std::vector< tk::real > > 
CompFlowFieldOutput( ncomp_t system,
                     ncomp_t offset,
                     tk::Fields& U )
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
  // number of degree of freedom
  const std::size_t rdof = g_inputdeck.get< tag::discr, tag::rdof >();

  std::vector< std::vector< tk::real > > out;
  const auto r  = U.extract( 0*rdof, offset );
  const auto ru = U.extract( 1*rdof, offset );
  const auto rv = U.extract( 2*rdof, offset );
  const auto rw = U.extract( 3*rdof, offset );
  const auto re = U.extract( 4*rdof, offset );

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

  std::vector< tk::real > P( r.size(), 0.0 );
  for (std::size_t i=0; i<P.size(); ++i) {
    P[i] = eos_pressure< tag::compflow >
             ( system, r[i], u[i], v[i], w[i], r[i]*E[i] );
  }
  out.push_back( P );

  return out;
}

std::vector< std::string > CompFlowSurfNames()
// *****************************************************************************
// Return surface field names to be output to file
//! \note Every surface will output these fields.
//! \return Vector of strings labelling surface fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;

  n.push_back( "density_numerical" );
  n.push_back( "x-velocity_numerical" );
  n.push_back( "y-velocity_numerical" );
  n.push_back( "z-velocity_numerical" );
  n.push_back( "specific_total_energy_numerical" );
  n.push_back( "pressure_numerical" );

  return n;
}

std::vector< std::vector< tk::real > >
CompFlowSurfOutput( ncomp_t system,
                    const std::map< int, std::vector< std::size_t > >& bnd,
                    tk::Fields& U )
// *****************************************************************************
//  Return surface field output going to file
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] bnd Boundary node/elem lists mapped to side set ids
//! \param[in] U Solution vector at recent time step
//! \return Vector of vectors of solution along side sets to be output to file
// *****************************************************************************
{
  std::vector< std::vector< tk::real > > out;

  // extract field output along side sets requested
  for (auto s : g_inputdeck.outsets()) {
    // get node list for side set requested
    auto b = bnd.find(s);
    if (b == end(bnd)) continue;
    const auto& nodes = b->second;
    std::vector< tk::real > surfaceSol( nodes.size() );
    auto i = out.size();
    out.insert( end(out), 6, surfaceSol );
    std::size_t j = 0;
    for (auto n : nodes) {
      const auto u = U.extract( n );
      Assert( u.size() == 5, "Size mismatch" );
      out[i+0][j] = u[0];
      out[i+1][j] = u[1];
      out[i+2][j] = u[2];
      out[i+3][j] = u[3];
      out[i+4][j] = u[4];
      out[i+5][j] = eos_pressure< tag::compflow >
                      ( system, u[0], u[1], u[2], u[3], u[4] );
      ++j;
    }
  }

  return out;
}

} //inciter::
