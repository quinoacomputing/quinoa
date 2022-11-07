// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/Problem/FieldOutput.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field outputs for single-material equation solver
  \details   This file defines functions for field quantites to be output to
    files for compressible single-material equations.
*/
// *****************************************************************************

#include "FieldOutput.hpp"
#include "EoS/EoS.hpp"
#include "EoS/EosVariant.hpp"
#include "ContainerUtil.hpp"
#include "History.hpp"

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
CompFlowFieldOutput( ncomp_t,
                     const std::vector< EOS >& mat_blk,
                     std::size_t nunk,
                     std::size_t rdof,
                     const tk::Fields& U )
// *****************************************************************************
//  Return field output going to file
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] nunk Number of unknowns to extract
//! \param[in] rdof Number of reconstructed degrees of freedom. This is used as
//!   the number of scalar components to shift when extracting scalar
//!   components.
//! \param[in] U Solution vector at recent time step
//! \return Vector of vectors to be output to file
// *****************************************************************************
{
  std::vector< std::vector< tk::real > > out;
  const auto r  = U.extract_comp( 0*rdof );
  const auto ru = U.extract_comp( 1*rdof );
  const auto rv = U.extract_comp( 2*rdof );
  const auto rw = U.extract_comp( 3*rdof );
  const auto re = U.extract_comp( 4*rdof );

  Assert( r.size() >= nunk, "Size mismatch" );
  Assert( ru.size() >= nunk, "Size mismatch" );
  Assert( rv.size() >= nunk, "Size mismatch" );
  Assert( rw.size() >= nunk, "Size mismatch" );
  Assert( re.size() >= nunk, "Size mismatch" );

  out.push_back( r );

  std::vector< tk::real > u = ru;
  for (std::size_t i=0; i<nunk; ++i) u[i] /= r[i];
  out.push_back( u );

  std::vector< tk::real > v = rv;
  for (std::size_t i=0; i<nunk; ++i) v[i] /= r[i];
  out.push_back( v );

  std::vector< tk::real > w = rw;
  for (std::size_t i=0; i<nunk; ++i) w[i] /= r[i];
  out.push_back( w );

  std::vector< tk::real > E = re;
  for (std::size_t i=0; i<nunk; ++i) E[i] /= r[i];
  out.push_back( E );

  std::vector< tk::real > P( nunk, 0.0 );
  for (std::size_t i=0; i<nunk; ++i) {
    P[i] = mat_blk[0].eosCall< EOS::pressure >( r[i], u[i], v[i], w[i],
      r[i]*E[i] );
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
CompFlowSurfOutput( ncomp_t,
                    const std::vector< EOS >& mat_blk,
                    const std::map< int, std::vector< std::size_t > >& bnd,
                    const tk::Fields& U )
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
      out[i+1][j] = u[1]/u[0];
      out[i+2][j] = u[2]/u[0];
      out[i+3][j] = u[3]/u[0];
      out[i+4][j] = u[4]/u[0];
      out[i+5][j] = mat_blk[0].eosCall< EOS::pressure >( u[0], u[1]/u[0],
        u[2]/u[0], u[3]/u[0], u[4] );
      ++j;
    }
  }

  return out;
}

std::vector< std::string > CompFlowHistNames()
// *****************************************************************************
// Return time history field names to be output to file
//! \note Every time history point will output these fields.
//! \return Vector of strings labelling time history fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;

  n.push_back( "density" );
  n.push_back( "x-velocity" );
  n.push_back( "y-velocity" );
  n.push_back( "z-velocity" );
  n.push_back( "energy" );
  n.push_back( "pressure" );

  return n;
}

std::vector< std::vector< tk::real > >
CompFlowHistOutput( ncomp_t,
                    const std::vector< EOS >& mat_blk,
                    const std::vector< HistData >& h,
                    const std::vector< std::size_t >& inpoel,
                    const tk::Fields& U )
// *****************************************************************************
//  Return time history field output evaluated at time history points
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] h History point data
//! \param[in] inpoel Mesh element connectivity
//! \param[in] U Solution vector at recent time step
//! \return Vector of vectors of solution variables evaluated in all history
//!   points. Inner vector: variables, outer vector: points.
// *****************************************************************************
{
  std::vector< std::vector< tk::real > > out( h.size() );

  std::size_t j = 0;
  for (const auto& p : h) {
    auto e = p.get< tag::elem >();        // host element id
    const auto& n = p.get< tag::fn >();   // shapefunctions evaluated at point
    out[j].resize( 6, 0.0 );
    for (std::size_t i=0; i<4; ++i) {
      const auto u = U.extract( inpoel[e*4+i] );
      Assert( u.size() == 5, "Size mismatch" );
      out[j][0] += n[i] * u[0];
      out[j][1] += n[i] * u[1]/u[0];
      out[j][2] += n[i] * u[2]/u[0];
      out[j][3] += n[i] * u[3]/u[0];
      out[j][4] += n[i] * u[4]/u[0];
      out[j][5] += n[i] * mat_blk[0].eosCall< EOS::pressure >( u[0], u[1]/u[0],
        u[2]/u[0], u[3]/u[0], u[4] );
    }
    ++j;
  }

  return out;
}

} //inciter::
