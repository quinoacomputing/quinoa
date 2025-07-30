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
#include "MultiMatIndexing.hpp"
#include "EoS/EOS.hpp"
#include "Kokkos_Core.hpp"

using execution_space = Kokkos::DefaultExecutionSpace;
using memory_space = Kokkos::DefaultExecutionSpace::memory_space;

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

//! Compute the time step size restriction based on viscosity
tk::real
timeStepSizeViscousFV(
  const tk::Fields& geoElem,
  std::size_t nelem,
  std::size_t nmat,
  const tk::Fields& U );

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

KOKKOS_INLINE_FUNCTION
void getDeformGrad(
  std::size_t nmat,
  std::size_t k,
  Kokkos::View<const size_t*, memory_space> solidx,
  Kokkos::View<const tk::real*, memory_space> state, 
  Kokkos::View<tk::real***, memory_space> g)
// *****************************************************************************
//  Get the inverse deformation gradient tensor for a material at given location
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] k Material id whose deformation gradient is required
//! \param[in] state State vector at location
//! \return Inverse deformation gradient tensor (g_k) of material k
// *****************************************************************************
{

  if (solidx(k) > 0) {
    // deformation gradient for solids
    for (std::size_t i=0; i<3; ++i) {
      for (std::size_t j=0; j<3; ++j)
        g(k, i, j) = state(deformIdx(nmat, solidx(k),i,j));
    }
  }
}

//! Get the elastic Cauchy stress tensor for a material at given location
std::array< std::array< tk::real, 3 >, 3 >
getCauchyStress(
  std::size_t nmat,
  std::size_t k,
  std::size_t ncomp,
  const std::vector< tk::real >& state );

KOKKOS_INLINE_FUNCTION
void getCauchyStress(
  std::size_t nmat,
  std::size_t k,
  std::size_t ncomp,
  Kokkos::View<const size_t*, memory_space> solidx,
  Kokkos::View<const tk::real*, memory_space> state,
  Kokkos::View<tk::real***, memory_space> asigk)  
// *****************************************************************************
//  Get the elastic Cauchy stress tensor for a material at given location
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] k Material id whose deformation gradient is required
//! \param[in] ncomp Number of components in the PDE system
//! \param[in] state State vector at location
//! \return Elastic Cauchy stress tensor (alpha * \sigma_ij) of material k
// *****************************************************************************
{
  const Kokkos::Array<Kokkos::Array<size_t, 3>, 3> stressCmpKokkos = {{
  {{0, 3, 4}},
  {{3, 1, 5}},
  {{4, 5, 2}} }};

  // elastic Cauchy stress for solids
  if (solidx(k) > 0) {
    for (std::size_t i=0; i<3; ++i) {
      for (std::size_t j=0; j<3; ++j)
        asigk(k, i, j) = state[ncomp +
          stressIdx(nmat, solidx(k), stressCmpKokkos[i][j])];
    }
  }

}

//! Check whether we have solid materials in our problem
bool
haveSolid(
  std::size_t nmat,
  const std::vector< std::size_t >& solidx );

KOKKOS_INLINE_FUNCTION
bool haveSolid(
  std::size_t nmat,
  Kokkos::View<const size_t*, memory_space>  solidx )
// *****************************************************************************
//  Check whether we have solid materials in our problem
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] solidx Material index indicator
//! \return true if we have at least one solid, false otherwise.
// *****************************************************************************
{
  bool haveSolid = false;
  for (std::size_t k=0; k<nmat; ++k)
    if (solidx(k) > 0) haveSolid = true;

  return haveSolid;
}

//! Count total number of solid materials in the problem
std::size_t numSolids(
  std::size_t nmat,
  const std::vector< std::size_t >& solidx );

} //inciter::

#endif // MiscMultiMatFns_h
