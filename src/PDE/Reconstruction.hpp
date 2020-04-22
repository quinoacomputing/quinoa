// *****************************************************************************
/*!
  \file      src/PDE/Reconstruction.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
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

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Compute lhs matrix for the least-squares reconstruction
void
lhsLeastSq_P0P1( const inciter::FaceData& fd,
  const Fields& geoElem,
  const Fields& geoFace,
  std::vector< std::array< std::array< real, 3 >, 3 > >& lhs_ls );

//! Compute internal surface contributions to the least-squares reconstruction
void
intLeastSq_P0P1( ncomp_t ncomp,
                 ncomp_t offset,
                 const std::size_t rdof,
                 const inciter::FaceData& fd,
                 const Fields& geoElem,
                 const Fields& W,
                 std::vector< std::vector< std::array< real, 3 > > >& rhs_ls );

//! \brief Compute boundary surface contributions to rhs vector of the
//!   least-squares reconstruction of conserved quantities of the PDE system
void
bndLeastSqConservedVar_P0P1( ncomp_t system,
  ncomp_t ncomp,
  ncomp_t offset,
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
  std::size_t nprim=0 );

//! \brief Compute boundary surface contributions to rhs vector of the
//!   least-squares reconstruction of primitive quantities of the PDE system
void
bndLeastSqPrimitiveVar_P0P1( ncomp_t system,
  ncomp_t nprim,
  ncomp_t offset,
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
  std::size_t ncomp );

//! Solve 3x3 system for least-squares reconstruction
void
solveLeastSq_P0P1(
  ncomp_t ncomp,
  ncomp_t offset,
  const std::size_t rdof,
  const std::vector< std::array< std::array< real, 3 >, 3 > >& lhs,
  const std::vector< std::vector< std::array< real, 3 > > >& rhs,
  Fields& W );

//! Transform the reconstructed P1-derivatives to the Dubiner dofs
void
transform_P0P1( ncomp_t ncomp,
                ncomp_t offset,
                std::size_t rdof,
                std::size_t nelem,
                const std::vector< std::size_t >& inpoel,
                const UnsMesh::Coords& coord,
                Fields& W );

//! \brief Compute MUSCL reconstruction in edge-end points using a MUSCL
//!   procedure with Van Leer limiting
#pragma omp declare simd
void
muscl( std::size_t p,
       std::size_t q,
       const UnsMesh::Coords& coord,
       const Fields& G,
       tk::real& rL, tk::real& uL, tk::real& vL, tk::real& wL, tk::real& eL,
       tk::real& rR, tk::real& uR, tk::real& vR, tk::real& wR, tk::real& eR,
       bool realizability = false );

} // tk::

#endif // Reconstruction_h
