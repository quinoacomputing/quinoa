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
  n.push_back( "density_numerical" );
  n.push_back( "x-velocity_numerical" );
  n.push_back( "y-velocity_numerical" );
  n.push_back( "z-velocity_numerical" );
  n.push_back( "pressure_numerical" );
  n.push_back( "total_energy_density_numerical" );

  return n;
}

std::vector< std::vector< tk::real > >
MultiMatFieldOutput(
  ncomp_t system,
  std::size_t nmat,
  ncomp_t offset,
  std::size_t rdof,
  tk::Fields& U )
// *****************************************************************************
//  Return field output going to file
//! \param[in] system Equation system index, i.e., which compressible
//!   flow equation system we operate on among the systems of PDEs
//! \param[in] nmat Number of materials in systen
//! \param[in] offset System offset specifying the position of the system of
//!   PDEs among other systems
//! \param[in] rdof Number of reconstructed degrees of freedom
//! \param[in] U Solution vector at recent time step
//! \return Vector of vectors to be output to file
// *****************************************************************************
{
  std::vector< std::vector< tk::real > > out;
  std::vector< std::vector< tk::real > > al, ar, ae;

  for (std::size_t k=0; k<nmat; ++k)
  {
    al.push_back( U.extract( volfracDofIdx(nmat, k, rdof, 0), offset ) );
    ar.push_back( U.extract( densityDofIdx(nmat, k, rdof, 0), offset ) );
    ae.push_back( U.extract( energyDofIdx(nmat, k, rdof, 0), offset ) );
  }
  const auto ru  = U.extract( momentumDofIdx(nmat, 0, rdof, 0), offset );
  const auto rv  = U.extract( momentumDofIdx(nmat, 1, rdof, 0), offset );
  const auto rw  = U.extract( momentumDofIdx(nmat, 2, rdof, 0), offset );

  //// mesh node coordinates
  //const auto& x = coord[0];
  //const auto& y = coord[1];
  //const auto& z = coord[2];

  // material volume-fractions
  for (std::size_t k=0; k<nmat; ++k)
    out.push_back( al[k] );

  // bulk density
  std::vector< tk::real > r( ru.size(), 0.0 );
  for (std::size_t i=0; i<r.size(); ++i) {
    for (std::size_t k=0; k<nmat; ++k)
      r[i] += ar[k][i];
  }
  out.push_back( r );

  // velocity components
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

  // bulk pressure
  std::vector< tk::real > P( r.size(), 0.0 );
  for (std::size_t i=0; i<P.size(); ++i) {
    for (std::size_t k=0; k<nmat; ++k)
      P[i] += eos_pressure< tag::multimat >( system, ar[k][i], u[i], v[i], w[i],
        ae[k][i], al[k][i], k );
  }
  out.push_back( P );

  // bulk total energy density
  std::vector< tk::real > E( r.size(), 0.0 );
  for (std::size_t i=0; i<E.size(); ++i) {
    for (std::size_t k=0; k<nmat; ++k)
      E[i] += ae[k][i];
  }
  out.push_back( E );

  return out;
}

} //inciter::
