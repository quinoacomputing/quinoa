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
#include "EoS/GetMatProp.hpp"
#include "EoS/EOS.hpp"
#include "ContainerUtil.hpp"
#include "History.hpp"

namespace inciter {

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
CompFlowSurfOutput( const std::vector< EOS >& mat_blk,
                    const std::map< int, std::vector< std::size_t > >& bnd,
                    const tk::Fields& U )
// *****************************************************************************
//  Return surface field output going to file
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
      out[i+5][j] = mat_blk[0].compute< EOS::pressure >( u[0], u[1]/u[0],
        u[2]/u[0], u[3]/u[0], u[4] );
      ++j;
    }
  }

  return out;
}

std::vector< std::vector< tk::real > >
CompFlowElemSurfOutput(
  const std::vector< EOS >& mat_blk,
  const std::map< int, std::vector< std::size_t > >& bface,
  const std::vector< std::size_t >& triinpoel,
  const tk::Fields& U )
// *****************************************************************************
//  Return element surface field output (on triangle faces) going to file
//! \param[in] mat_blk Material EOS block
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] triinpoel Boundary triangle face connecitivity with local ids
//! \param[in] U Solution vector at recent time step
//! \return Vector of vectors of solution on side set faces to be output to file
// *****************************************************************************
{
  std::vector< std::vector< tk::real > > out;

  // extract field output along side sets requested
  for (auto s : g_inputdeck.outsets()) {
    // get face list for side set requested
    auto b = bface.find(s);
    if (b == end(bface)) continue;
    const auto& faces = b->second;
    std::vector< tk::real > surfaceSol( faces.size() );
    auto i = out.size();
    out.insert( end(out), 6, surfaceSol );
    std::size_t j = 0;
    for (auto f : faces) {
      // access solutions at nodes
      auto UA = U.extract(triinpoel[f*3+0]);
      auto UB = U.extract(triinpoel[f*3+1]);
      auto UC = U.extract(triinpoel[f*3+2]);

      auto rho = (UA[0] + UB[0] + UC[0]) / 3.0;
      auto u = (UA[1]/UA[0] + UB[1]/UB[0] + UC[1]/UC[0]) / 3.0;
      auto v = (UA[2]/UA[0] + UB[2]/UB[0] + UC[2]/UC[0]) / 3.0;
      auto w = (UA[3]/UA[0] + UB[3]/UB[0] + UC[3]/UC[0]) / 3.0;
      auto E = (UA[4]/UA[0] + UB[4]/UB[0] + UC[4]/UC[0]) / 3.0;

      out[i+0][j] = rho;
      out[i+1][j] = u;
      out[i+2][j] = v;
      out[i+3][j] = w;
      out[i+4][j] = E;
      out[i+5][j] = mat_blk[0].compute< EOS::pressure >( rho, u, v, w, rho*E );
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
CompFlowHistOutput( const std::vector< EOS >& mat_blk,
                    const std::vector< HistData >& h,
                    const std::vector< std::size_t >& inpoel,
                    const tk::Fields& U )
// *****************************************************************************
//  Return time history field output evaluated at time history points
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
      out[j][5] += n[i] * mat_blk[0].compute< EOS::pressure >( u[0], u[1]/u[0],
        u[2]/u[0], u[3]/u[0], u[4] );
    }
    ++j;
  }

  return out;
}

} //inciter::
