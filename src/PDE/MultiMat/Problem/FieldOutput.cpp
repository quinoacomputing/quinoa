// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/FieldOutput.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field outputs for multi-material equation solver
  \details   This file defines functions for field quantites to be output to
    files for compressible multi-material equations.
*/
// *****************************************************************************
#include "FieldOutput.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "EoS/EoS.hpp"
#include "Vector.hpp"

namespace inciter {

std::vector< std::string >
MultiMatFieldNames( std::size_t nmat )
// *****************************************************************************
// Return multi-material field names to be output to file
//! \param[in] nmat Number of materials in system
//! \return Vector of strings labelling fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;

  for (std::size_t k=0; k<nmat; ++k)
    n.push_back( "volfrac"+std::to_string(k+1)+"_numerical" );
  for (std::size_t k=0; k<nmat; ++k)
    n.push_back( "density"+std::to_string(k+1)+"_numerical" );
  n.push_back( "density_numerical" );
  n.push_back( "x-velocity_numerical" );
  n.push_back( "y-velocity_numerical" );
  n.push_back( "z-velocity_numerical" );
  for (std::size_t k=0; k<nmat; ++k)
    n.push_back( "pressure"+std::to_string(k+1)+"_numerical" );
  n.push_back( "pressure_numerical" );
  for (std::size_t k=0; k<nmat; ++k)
    n.push_back( "soundspeed"+std::to_string(k+1) );
  n.push_back( "total_energy_density_numerical" );
  n.push_back( "material_indicator" );
  n.push_back( "timestep" );

  return n;
}

std::vector< std::vector< tk::real > >
MultiMatFieldOutput(
  ncomp_t,
  std::size_t nmat,
  const std::vector< EoS_Base* >& mat_blk,
  std::size_t nunk,
  std::size_t rdof,
  const std::vector< tk::real >&,
  const std::array< std::vector< tk::real >, 3 >&,
  const tk::Fields& U,
  const tk::Fields& P )
// *****************************************************************************
//  Return field output going to file
//! \param[in] nmat Number of materials in systen
//! \param[in] nunk Number of unknowns to extract
//! \param[in] rdof Number of reconstructed degrees of freedom
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitive quantities at recent time step
//! \return Vector of vectors to be output to file
// *****************************************************************************
{
  // field-output vector with:
  // - nmat volume fractions
  // - nmat material densities
  // - 1 bulk density
  // - 3 velocity components
  // - nmat material pressures
  // - 1 bulk pressure
  // - nmat material sound-speeds
  // - 1 bulk total energy density
  // - 1 material indicator
  // - 1 time-step
  // leading to a size of 4*nmat+8
  std::vector< std::vector< tk::real > >
    out( 4*nmat+8, std::vector< tk::real >( nunk ) );

  //// mesh node coordinates
  //const auto& x = coord[0];
  //const auto& y = coord[1];
  //const auto& z = coord[2];

  // material volume-fractions
  for (std::size_t k=0; k<nmat; ++k) {
    for (std::size_t i=0; i<nunk; ++i)
      out[k][i] = U(i, volfracDofIdx(nmat, k, rdof, 0));
  }

  // material densities
  for (std::size_t k=0; k<nmat; ++k) {
    for (std::size_t i=0; i<nunk; ++i) {
      out[nmat+k][i] = U(i, densityDofIdx(nmat, k, rdof, 0)) /
        std::max(1e-16, U(i, volfracDofIdx(nmat, k, rdof, 0)));
    }
  }

  // bulk density
  for (std::size_t i=0; i<nunk; ++i) {
    for (std::size_t k=0; k<nmat; ++k)
      out[2*nmat][i] += U(i, densityDofIdx(nmat, k, rdof, 0));
  }

  // velocity components
  for (std::size_t i=0; i<nunk; ++i) {
    out[2*nmat+1][i] = P(i, velocityDofIdx(nmat, 0, rdof, 0));
    out[2*nmat+2][i] = P(i, velocityDofIdx(nmat, 1, rdof, 0));
    out[2*nmat+3][i] = P(i, velocityDofIdx(nmat, 2, rdof, 0));
  }

  // material pressures
  for (std::size_t k=0; k<nmat; ++k) {
    for (std::size_t i=0; i<nunk; ++i) {
      out[2*nmat+4+k][i] =
        P(i, pressureDofIdx(nmat, k, rdof, 0)) /
        std::max(1e-16, U(i, volfracDofIdx(nmat, k, rdof, 0)));
    }
  }

  // bulk pressure
  for (std::size_t i=0; i<nunk; ++i) {
    for (std::size_t k=0; k<nmat; ++k)
      out[3*nmat+4][i] += P(i, pressureDofIdx(nmat, k, rdof, 0));
  }

  // material sound speeds
  for (std::size_t k=0; k<nmat; ++k) {
    for (std::size_t i=0; i<nunk; ++i) {
      out[3*nmat+5+k][i] =
      mat_blk[k]->eos_soundspeed(
        std::max(1e-16, U(i, densityDofIdx(nmat,k,rdof,0))),
        P(i, pressureDofIdx(nmat,k,rdof,0)),
        U(i, volfracDofIdx(nmat,k,rdof,0)) );
    }
  }

  // bulk total energy density
  for (std::size_t i=0; i<nunk; ++i) {
    for (std::size_t k=0; k<nmat; ++k)
      out[4*nmat+5][i] += U(i, energyDofIdx(nmat, k, rdof, 0));
  }

  // material indicator
  for (std::size_t i=0; i<nunk; ++i) {
    out[4*nmat+6][i] = 0.0;
    for (std::size_t k=0; k<nmat; ++k) {
      out[4*nmat+6][i] += U(i, volfracDofIdx(nmat,k,rdof,0))
        * static_cast< tk::real >(k+1);
    }
  }

  // time-step
  //for (std::size_t i=0; i<nunk; ++i) {
  //  // advection velocity
  //  auto u = U(i, velocityDofIdx(nmat,0,rdof,0));
  //  auto v = U(i, velocityDofIdx(nmat,1,rdof,0));
  //  auto w = U(i, velocityDofIdx(nmat,2,rdof,0));

  //  auto vn = std::sqrt(tk::dot({{u, v, w}}, {{u, v, w}}));

  //  // acoustic speed
  //  auto a = 0.0;
  //  for (std::size_t k=0; k<nmat; ++k)
  //  {
  //    if (U(i, volfracDofIdx(nmat,k,rdof,0)) > 1.0e-04) {
  //      a = std::max( a, eos_soundspeed< tag::multimat >( 0,
  //        U(i, densityDofIdx(nmat,k,rdof,0)),
  //        P(i, pressureDofIdx(nmat,k,rdof,0)),
  //        U(i, volfracDofIdx(nmat,k,rdof,0)), k ) );
  //    }
  //  }

  //  out[4*nmat+7][i] = geoElem(i,4) / (std::fabs(vn) + a);
  //}

  return out;
}

std::vector< std::string > MultiMatHistNames()
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

std::vector< std::string > MultiMatDiagNames(std::size_t nmat)
// *****************************************************************************
// Return diagnostic var names to be output to file
//! \param[in] nmat Number of materials in systen
//! \return Vector of strings labelling diagnostic fields output in file
// *****************************************************************************
{
  std::vector< std::string > n;

  for (std::size_t k=0; k<nmat; ++k)
    n.push_back( "f"+std::to_string(k+1) );
  for (std::size_t k=0; k<nmat; ++k)
    n.push_back( "fr"+std::to_string(k+1) );
  n.push_back( "fru" );
  n.push_back( "frv" );
  n.push_back( "frw" );
  for (std::size_t k=0; k<nmat; ++k)
    n.push_back( "fre"+std::to_string(k+1) );

  return n;
}

} //inciter::
