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
#include "EoS/EOS.hpp"

namespace inciter {

void initializeMaterialEoS( std::vector< EOS >& mat_blk );

//! Clean up the state of trace materials for multi-material PDE system
bool
cleanTraceMultiMat(
  tk::real t,
  std::size_t nelem,
  const std::vector< EOS >& mat_blk,
  const tk::Fields& geoElem,
  std::size_t nmat,
  tk::Fields& U,
  tk::Fields& P );

//! Time step restriction for multi material cell-centered schemes
tk::real
timeStepSizeMultiMat(
  const std::vector< EOS >& mat_blk,
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
  const std::vector< EOS >& mat_blk,
  const tk::Fields& geoElem,
  std::size_t nelem,
  std::size_t nmat,
  const tk::Fields& U,
  const tk::Fields& P,
  std::vector< tk::real >& local_dte );

//! Reset the solid tensors
void
resetSolidTensors(
  std::size_t nmat,
  std::size_t k,
  std::size_t e,
  tk::Fields& U,
  tk::Fields& P );

//! Get the inverse deformation gradient tensor for a material at given location
std::array< std::array< tk::real, 3 >, 3 >
getDeformGrad(
  std::size_t nmat,
  std::size_t k,
  const std::vector< tk::real >& state );

//! Get the elastic Cauchy stress tensor for a material at given location
std::array< std::array< tk::real, 3 >, 3 >
getCauchyStress(
  std::size_t nmat,
  std::size_t k,
  std::size_t ncomp,
  const std::vector< tk::real >& state );

//! Check whether we have solid materials in our problem
bool
haveSolid(
  std::size_t nmat,
  const std::vector< std::size_t >& solidx );

//! Count total number of solid materials in the problem
std::size_t numSolids(
  std::size_t nmat,
  const std::vector< std::size_t >& solidx );

} //inciter::

#endif // MiscMultiMatFns_h
