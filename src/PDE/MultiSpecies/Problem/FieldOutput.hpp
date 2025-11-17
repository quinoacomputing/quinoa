// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/Problem/FieldOutput.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Field outputs for multi-species equation solver
  \details   This file defines functions for field quantites to be output to
    files for compressible multi-species equations.
*/
// *****************************************************************************
#ifndef FieldOutput_h
#define FieldOutput_h

#include "Types.hpp"
#include "Fields.hpp"
#include "EoS/EOS.hpp"
#include "FaceData.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "PDE/MultiSpecies/MultiSpeciesIndexing.hpp"
#include "ContainerUtil.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

using ncomp_t = tk::ncomp_t;

//! Return a map that associates user-specified strings to functions
std::map< std::string, tk::GetVarFn > MultiSpeciesOutVarFn();

//! Return multi-species field names to be output to file
std::vector< std::string >
MultiSpeciesFieldNames( std::size_t nspec );

//! Return surface field names to be output to file
std::vector< std::string > MultiSpeciesSurfNames();

//! Return element surface field output (on triangle faces) going to file
std::vector< std::vector< tk::real > >
MultiSpeciesSurfOutput(
  const std::size_t nspec,
  const std::size_t rdof,
  const FaceData& fd,
  const tk::Fields& U,
  const tk::Fields& P );

//! Return time history field names to be output to file
std::vector< std::string > MultiSpeciesHistNames();

//! Return diagnostic var names to be output to file
std::vector< std::string > MultiSpeciesDiagNames(std::size_t nspec);

/** @name Functions that compute physics variables from the numerical solution for MultiSpecies */
///@{

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-function"
#endif

namespace multispecies {

//! Compute mixture density for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Bulk density ready to be output to file
static tk::GetVarFn::result_type
mixDensityOutVar( const tk::Fields& U, std::size_t rdof )
{
  using tk::operator+=;
  auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();
  auto r = U.extract_comp( densityDofIdx(nspec,0,rdof,0) );
  for (std::size_t k=1; k<nspec; ++k)
    r += U.extract_comp( densityDofIdx(nspec,k,rdof,0) );
  return r;
}

////! Compute pressure for output to file
////! \note Must follow the signature in tk::GetVarFn
////! \param[in] U Numerical solution
////! \param[in] rdof Number of reconstructed solution DOFs
////! \return Pressure ready to be output to file
//static tk::GetVarFn::result_type
//pressureOutVar( const tk::Fields& U, std::size_t rdof )
//{
//  using tk::operator+=;
//  auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();
//  auto p = U.extract_comp( pressureDofIdx(nspec,0,rdof,0) );
//  for (std::size_t k=1; k<nspec; ++k)
//    p += U.extract_comp( pressureDofIdx(nspec,k,rdof,0) );
//  return p;
//}

//! Compute specific total energy (energy per unit volume) for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Specific total energy ready to be output to file
static tk::GetVarFn::result_type
specificTotalEnergyOutVar( const tk::Fields& U, std::size_t rdof )
{
  auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();
  return U.extract_comp( energyDofIdx(nspec,0,rdof,0) );
}

//! Compute velocity component for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \tparam dir Physical direction, encoded as 0:x, 1:y, 2:z
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Velocity component ready to be output to file
template< tk::ncomp_t dir >
tk::GetVarFn::result_type
velocityOutVar( const tk::Fields& U, std::size_t rdof )
{
  using tk::operator/=;
  using tk::operator+=;
  auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

  // mixture density
  auto r = U.extract_comp( densityDofIdx(nspec,0,rdof,0) );
  for (std::size_t k=1; k<nspec; ++k)
    r += U.extract_comp( densityDofIdx(nspec,k,rdof,0) );

  // momentum
  auto u = U.extract_comp( momentumDofIdx(nspec,dir,rdof,0) );

  // velocity
  u /= r;

  return u;
}

//! Compute mixture temperature for output to file
//! \note Must follow the signature in tk::GetVarFn
//! \param[in] U Numerical solution
//! \param[in] rdof Number of reconstructed solution DOFs
//! \return Mixture temperature ready to be output to file
static tk::GetVarFn::result_type
temperatureOutVar( const tk::Fields& U, std::size_t rdof )
{
  auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();
  auto r = U.extract_comp( temperatureDofIdx(nspec,0,rdof,0) );
  return r;
}

} // multispecies::

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//@}

} //inciter::

#endif // FieldOutput_h
