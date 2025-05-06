// *****************************************************************************
/*!
  \file      src/PDE/ConfigureMultiSpecies.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration for multi-species compressible
     flow PDE
*/
// *****************************************************************************
#ifndef ConfigureMultiSpecies_h
#define ConfigureMultiSpecies_h

#include <set>
#include <map>
#include <vector>

#include "PDEFactory.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Inciter/Options/PDE.hpp"
#include "PDE/MultiSpecies/MultiSpeciesIndexing.hpp"
#include "ContainerUtil.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! Register compressible flow PDEs into PDE factory
void
registerMultiSpecies( DGFactory& df, FVFactory& ff,
  std::set< ctr::PDEType >& fvt, std::set< ctr::PDEType >& dgt );

//! Return information on the multi-species compressible flow PDE
std::vector< std::pair< std::string, std::string > >
infoMultiSpecies( std::map< ctr::PDEType, tk::ncomp_t >& cnt );

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

} // inciter::

#endif // ConfigureMultiSpecies_h
