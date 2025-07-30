// *****************************************************************************
/*!
  \file      src/PDE/Reconstruction.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Reconstruction for reconstructed Galerkin methods
  \details   This file contains functions that reconstruct an "n"th order
    polynomial to an "n+1"th order polynomial using a least-squares
    reconstruction procedure, used for reconstructed discontinuous Galerkin (DG)
    methods. It also contains functions used to compute reconstruction in 1D,
    used in edge-based continuous Galerkin methods.
*/
// *****************************************************************************
#ifndef Reconstruction_h
#define Reconstruction_h

#include "Types.hpp"
#include "Fields.hpp"
#include "FaceData.hpp"
#include "UnsMesh.hpp"
#include "FunctionPrototypes.hpp"
#include "Integrate/Basis.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "EoS/EOS.hpp"
#include "Kokkos_Core.hpp"

using execution_space = Kokkos::DefaultExecutionSpace;
using memory_space = Kokkos::DefaultExecutionSpace::memory_space;

namespace tk {

using ncomp_t = tk::ncomp_t;

//! \brief Reconstruct the second-order solution using least-squares approach
//!   from an extended stencil involving the node-neighbors
void
recoLeastSqExtStencil(
  std::size_t rdof,
  std::size_t e,
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const Fields& geoElem,
  Fields& W,
  const std::vector< std::size_t >& varList );

//! Transform the reconstructed P1-derivatives to the Dubiner dofs
void
transform_P0P1( std::size_t rdof,
                std::size_t e,
                const std::vector< std::size_t >& inpoel,
                const UnsMesh::Coords& coord,
                Fields& W,
                const std::vector< std::size_t >& varList );

//! Compute THINC reconstructions near material interfaces
void
THINCReco(
  std::size_t rdof,
  std::size_t nmat,
  std::size_t e,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const Fields& geoElem,
  const std::array< real, 3 >& xp,
  const Fields& U,
  const Fields& P,
  bool intInd,
  const std::vector< std::size_t >& matInt,
  const std::vector< real >& vfmin,
  const std::vector< real >& vfmax,
  std::vector< real >& state );

//! Compute THINC reconstructions for linear advection (transport)
void
THINCRecoTransport(
  std::size_t rdof,
  std::size_t,
  std::size_t e,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const Fields& geoElem,
  const std::array< real, 3 >& ref_xp,
  const Fields& U,
  const Fields&,
  [[maybe_unused]] const std::vector< real >& vfmin,
  [[maybe_unused]] const std::vector< real >& vfmax,
  std::vector< real >& state );

//! Kokkos Version of above
KOKKOS_INLINE_FUNCTION void
THINCReco( std::size_t rdof,
           std::size_t nmat,
           std::size_t e, 
           std::size_t ncomp,
           std::size_t m_nprop,
           std::size_t p_nprop,
           std::size_t geo_nprop,
            tk::real bparam,
           Kokkos::View<const size_t*, memory_space> inpoel,
           Kokkos::View<const real*, memory_space> cx,
           Kokkos::View<const real*, memory_space> cy, 
           Kokkos::View<const real*, memory_space> cz,
           Kokkos::View<const real*, memory_space> geoElem,
           const Kokkos::Array<real, 3>& ref_xp,
           Kokkos::View<const real*, memory_space> U,
           Kokkos::View<const real*, memory_space> P,
           bool intInd,
           Kokkos::View<const size_t*, memory_space> solidx,
           Kokkos::View<size_t*, memory_space> matInt,
           [[maybe_unused]] Kokkos::View<const real*, memory_space> vfmin,
           [[maybe_unused]] Kokkos::View<const real*, memory_space> vfmax,
           Kokkos::View<real*, memory_space> state, 
          Kokkos::View<real*, memory_space> alSol, 
        Kokkos::View<real*, memory_space> alReco,
        Kokkos::View<real**, memory_space> dBdx, 
        Kokkos::View<Kokkos::Array<real, 3>*, memory_space> ref_n);

//! Old THINC reconstruction function for volume fractions near interfaces
void
THINCFunction_old( std::size_t rdof,
  std::size_t nmat,
  std::size_t e,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const std::array< real, 3 >& ref_xp,
  real vol,
  real bparam,
  const std::vector< real >& alSol,
  bool intInd,
  const std::vector< std::size_t >& matInt,
  std::vector< real >& alReco );

KOKKOS_INLINE_FUNCTION
void THINCFunction( std::size_t rdof,
    std::size_t nmat,
    std::size_t e,
    Kokkos::View<const size_t*, memory_space> inpoel,
    Kokkos::View<const real*, memory_space> cx,
    Kokkos::View<const real*, memory_space> cy,
    Kokkos::View<const real*, memory_space> cz,
    const Kokkos::Array<real, 3>& ref_xp,
    real vol,
    real bparam,
    Kokkos::View<const real*, memory_space> alSol,
    bool intInd,
    Kokkos::View<const size_t*, memory_space> matInt,
    Kokkos::View<real*, memory_space> alReco,
    Kokkos::View<real**, memory_space> dBdx,
    Kokkos::View<Kokkos::Array<real, 3>*, memory_space> ref_n);

//! New THINC reconstruction function for volume fractions near interfaces
void
THINCFunction( std::size_t rdof,
  std::size_t nmat,
  std::size_t e,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const std::array< real, 3 >& ref_xp,
  real vol,
  real bparam,
  const std::vector< real >& alSol,
  bool intInd,
  const std::vector< std::size_t >& matInt,
  std::vector< real >& alReco );

//! Compute the temperatures based on FV conserved quantities
void
computeTemperaturesFV(
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t nmat,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const tk::Fields& geoElem,
  const tk::Fields& unk,
  const tk::Fields& prim,
  const std::vector< int >& srcFlag,
  tk::Fields& T );

//! Evaluate polynomial solution at quadrature point
std::vector< tk::real >
evalPolynomialSol(
  const std::vector< inciter::EOS >& mat_blk,
  int intsharp,
  std::size_t ncomp,
  std::size_t nprim,
  std::size_t rdof,
  std::size_t nmat, 
  std::size_t e,
  std::size_t dof_e,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const Fields& geoElem,
  const std::array< real, 3 >& ref_gp,
  const std::vector< real >& B,
  const Fields& U,
  const Fields& P);

  //! Kokkos evalPolySolution
KOKKOS_FUNCTION
void evalPolynomialSol( const std::vector< inciter::EOS >& mat_blk,
    int intsharp,
    std::size_t ncomp,
    std::size_t nprim,
    std::size_t rdof,
    std::size_t nmat,
    std::size_t e,
    std::size_t dof_e,
    std::size_t m_nprop,
    std::size_t p_nprop,
    std::size_t geo_nprop,
    tk::real bparam,
    Kokkos::View<const size_t*, memory_space> solidx,
    Kokkos::View<const size_t*, memory_space> inpoel,
    Kokkos::View<const real*, memory_space> cx,
    Kokkos::View<const real*, memory_space> cy,
    Kokkos::View<const real*, memory_space> cz,
    Kokkos::View<const real*, memory_space> geoElem,
    const Kokkos::Array<real, 3>& ref_gp,
     Kokkos::View<const tk::real*, memory_space> B,
    Kokkos::View<const real*, memory_space> U,
    Kokkos::View<const real*, memory_space> P,
  Kokkos::View<real*, memory_space> state, 
  Kokkos::View<size_t*, memory_space> matInt,
  Kokkos::View<real*, memory_space> alAvg, 
  Kokkos::View<real*, memory_space> vfmax, 
  Kokkos::View<real*, memory_space> vfmin,
  Kokkos::View<real*, memory_space> alSol, 
  Kokkos::View<real*, memory_space> alReco,
  Kokkos::View<real**, memory_space> dBdx, 
   Kokkos::View<Kokkos::Array<real, 3>*, memory_space> ref_n);

//! Evaluate second-order FV solution at quadrature point
std::vector< tk::real >
evalFVSol(
  const std::vector< inciter::EOS >& mat_blk,
  int intsharp,
  std::size_t ncomp,
  std::size_t nprim,
  std::size_t rdof,
  std::size_t nmat,
  std::size_t e,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const Fields& geoElem,
  const std::array< real, 3 >& ref_gp,
  const std::vector< real >& B,
  const Fields& U,
  const Fields& P,
  int srcFlag );

//! Enforce physical constraints on state at quadrature point
void
enforcePhysicalConstraints(
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t nmat,
  std::size_t ncomp,
  std::vector< tk::real >& state );

void
enforcePhysicalConstraints(
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t nmat,
  std::size_t ncomp,
  Kokkos::View<real*, memory_space> state);

//! Compute safe reconstructions near material interfaces
void
safeReco( std::size_t rdof,
          std::size_t nmat,
          std::size_t el,
          int er,
          const Fields& U,
          std::array< std::vector< real >, 2 >& state );

} // tk::

#endif // Reconstruction_h
