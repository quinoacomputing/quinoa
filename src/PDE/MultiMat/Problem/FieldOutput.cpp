// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/Problem/FieldOutput.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field outputs for multi-material equation solver
  \details   This file defines functions for field quantites to be output to
    files for compressible multi-material equations.
*/
// *****************************************************************************
#include "FieldOutput.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "EoS/EoS.hpp"

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
  n.push_back( "total_energy_density_numerical" );

  return n;
}

std::vector< std::vector< tk::real > >
MultiMatFieldOutput(
  ncomp_t,
  std::size_t nmat,
  ncomp_t offset,
  std::size_t rdof,
  tk::Fields& U,
  const tk::Fields& P )
// *****************************************************************************
//  Return field output going to file
//! \param[in] nmat Number of materials in systen
//! \param[in] offset System offset specifying the position of the system of
//!   PDEs among other systems
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
  // - 1 bulk total energy density
  // leading to a size of 3*nmat+6
  std::vector< std::vector< tk::real > >
    out( 3*nmat+6, std::vector< tk::real >( U.nunk(), 0.0 ) );

  //// mesh node coordinates
  //const auto& x = coord[0];
  //const auto& y = coord[1];
  //const auto& z = coord[2];

  // material volume-fractions
  for (std::size_t k=0; k<nmat; ++k) {
    for (std::size_t i=0; i<U.nunk(); ++i)
      out[k][i] = U(i, volfracDofIdx(nmat, k, rdof, 0), offset);
  }

  // material densities
  for (std::size_t k=0; k<nmat; ++k) {
    for (std::size_t i=0; i<U.nunk(); ++i) {
      out[nmat+k][i] = U(i, densityDofIdx(nmat, k, rdof, 0), offset)
        /U(i, volfracDofIdx(nmat, k, rdof, 0), offset);
    }
  }

  // bulk density
  for (std::size_t i=0; i<U.nunk(); ++i) {
    for (std::size_t k=0; k<nmat; ++k)
      out[2*nmat][i] += U(i, densityDofIdx(nmat, k, rdof, 0), offset);
  }

  // velocity components
  for (std::size_t i=0; i<U.nunk(); ++i) {
    out[2*nmat+1][i] = P(i, velocityDofIdx(nmat, 0, rdof, 0), offset);
    out[2*nmat+2][i] = P(i, velocityDofIdx(nmat, 1, rdof, 0), offset);
    out[2*nmat+3][i] = P(i, velocityDofIdx(nmat, 2, rdof, 0), offset);
  }

  // material pressures
  for (std::size_t k=0; k<nmat; ++k) {
    for (std::size_t i=0; i<U.nunk(); ++i) {
      out[2*nmat+4+k][i] = P(i, pressureDofIdx(nmat, k, rdof, 0), offset)
        /U(i, volfracDofIdx(nmat, k, rdof, 0), offset);
    }
  }

  // bulk pressure
  for (std::size_t i=0; i<U.nunk(); ++i) {
    for (std::size_t k=0; k<nmat; ++k)
      out[3*nmat+4][i] += P(i, pressureDofIdx(nmat, k, rdof, 0), offset);
  }

  // bulk total energy density
  for (std::size_t i=0; i<U.nunk(); ++i) {
    for (std::size_t k=0; k<nmat; ++k)
      out[3*nmat+5][i] += U(i, energyDofIdx(nmat, k, rdof, 0), offset );
  }

  return out;
}

} //inciter::
