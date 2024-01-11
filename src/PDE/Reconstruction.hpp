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

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Compute lhs matrix for the least-squares reconstruction
void
lhsLeastSq_P0P1(
  const inciter::FaceData& fd,
  const Fields& geoElem,
  const Fields& geoFace,
  std::vector< std::array< std::array< real, 3 >, 3 > >& lhs_ls );

//! Compute internal surface contributions to the least-squares reconstruction
void
intLeastSq_P0P1( const std::size_t rdof,
                 const inciter::FaceData& fd,
                 const Fields& geoElem,
                 const Fields& W,
                 std::vector< std::vector< std::array< real, 3 > > >& rhs_ls,
                 const std::vector< std::size_t >& varList );

//! \brief Compute boundary surface contributions to rhs vector of the
//!   least-squares reconstruction of conserved quantities of the PDE system
void
bndLeastSqConservedVar_P0P1(
  ncomp_t ncomp,
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t rdof,
  const std::vector< bcconf_t >& bcconfig,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  real t,
  const StateFn& state,
  const Fields& P,
  const Fields& U,
  std::vector< std::vector< std::array< real, 3 > > >& rhs_ls,
  const std::vector< std::size_t >& varList,
  std::size_t nprim=0 );

//! Solve 3x3 system for least-squares reconstruction
void
solveLeastSq_P0P1(
  const std::size_t rdof,
  const std::vector< std::array< std::array< real, 3 >, 3 > >& lhs,
  const std::vector< std::vector< std::array< real, 3 > > >& rhs,
  Fields& W,
  const std::vector< std::size_t >& varList );

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
                std::size_t nelem,
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

//! Old THINC reconstruction function for volume fractions near interfaces
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

//! New THINC reconstruction function for volume fractions near interfaces
void
THINCFunction_new( std::size_t rdof,
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
