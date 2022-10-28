// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/MiscMultiMatFns.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Misc multi-material system functions
  \details   This file defines functions that required for multi-material
    compressible hydrodynamics.
*/
// *****************************************************************************
#ifndef MiscMultiMatFns_h
#define MiscMultiMatFns_h

#include "Fields.hpp"
#include "UnsMesh.hpp"
#include "EoS/EoS_Base.hpp"

namespace inciter {

void initializeMaterialEoS( std::size_t system,
  std::vector< EoS_Base* >& mat_blk );

//! Clean up the state of trace materials for multi-material PDE system
bool
cleanTraceMultiMat(
  std::size_t nelem,
  const std::vector< EoS_Base* >& mat_blk,
  const tk::Fields& geoElem,
  std::size_t nmat,
  tk::Fields& U,
  tk::Fields& P );

//! Time step restriction for multi material cell-centered schemes
tk::real
timeStepSizeMultiMat(
  const std::vector< EoS_Base* >& mat_blk,
  const std::vector< int >& esuf,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const std::size_t nelem,
  std::size_t nmat,
  const tk::Fields& U,
  const tk::Fields& P );

//! Time step restriction for multi material cell-centered FV scheme
tk::real
timeStepSizeMultiMatFV(
  const std::vector< EoS_Base* >& mat_blk,
  const tk::Fields& geoElem,
  const std::size_t nelem,
  std::size_t system,
  std::size_t nmat,
  const int engSrcAd,
  const tk::Fields& U,
  const tk::Fields& P );

} //inciter::

#endif // MiscMultiMatFns_h
