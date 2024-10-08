// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/MiscMultiSpeciesFns.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Misc multi-species system functions
  \details   This file defines functions that required for multi-species
    compressible fluid dynamics.
*/
// *****************************************************************************
#ifndef MiscMultiSpeciesFns_h
#define MiscMultiSpeciesFns_h

#include "Fields.hpp"
#include "UnsMesh.hpp"
#include "EoS/EOS.hpp"

namespace inciter {

void initializeSpeciesEoS( std::vector< EOS >& mat_blk );

//! Time step restriction for multi material cell-centered schemes
tk::real
timeStepSizeMultiSpecies(
  const std::vector< EOS >& mat_blk,
  const std::vector< int >& esuf,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const std::size_t nelem,
  std::size_t nspec,
  const tk::Fields& U,
  const tk::Fields& P );

} //inciter::

#endif // MiscMultiSpeciesFns_h
