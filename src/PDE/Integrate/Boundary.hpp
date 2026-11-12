// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Boundary.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing physical boundary surface integrals of a
     system of PDEs in DG methods
  \details   This file contains functionality for computing physical boundary
     surface integrals of a system of PDEs used in discontinuous Galerkin
     methods for various orders of numerical representation.
*/
// *****************************************************************************
#ifndef Boundary_h
#define Boundary_h

#include "Basis.hpp"
#include "Surface.hpp"
#include "Fields.hpp"
#include "FaceData.hpp"
#include "UnsMesh.hpp"
#include "FunctionPrototypes.hpp"
#include "EoS/EOS.hpp"

namespace tk {

using ncomp_t = tk::ncomp_t;

//! Compute boundary surface flux integrals for a given boundary type for DG
void
bndSurfInt( const bool pref,
            std::size_t nmat,
            const std::vector< inciter::EOS >& mat_blk,
            const std::size_t ndof,
            const std::size_t rdof,
            const std::vector< std::size_t >& bcconfig,
            const inciter::FaceData& fd,
            const Fields& geoFace,
            const Fields& geoElem,
            const std::vector< std::size_t >& inpoel,
            const UnsMesh::Coords& coord,
            real t,
            const RiemannFluxFn& flux,
            const VelFn& vel,
            const StateFn& state,
            const Fields& U,
            const Fields& P,
            const Fields& W,
            const std::vector< std::size_t >& ndofel,
            const std::vector< int >& srcFlag,
            Fields& R,
            std::vector< std::vector< tk::real > >& riemannDeriv,
            int intcompr=0 );

//! Update the rhs by adding the boundary surface integration term
void
update_rhs_bc ( ncomp_t ncomp,
                std::size_t nmat,
                const std::size_t ndof,
                const std::size_t ndof_l,
                const tk::real wt,
                const std::array< tk::real, 3 >& fn,
                const std::size_t el,
                const std::vector< tk::real >& fl,
                const std::vector< tk::real >& B_l,
                Fields& R,
                std::vector< std::vector< tk::real > >& riemannDeriv );

//! \brief Compute boundary surface flux integrals for a given boundary type for
//!   const-order DG (not p-adaptive)
void
bndSurfInt_constP(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t ndof,
  const std::size_t rdof,
  const std::vector< std::size_t >& bcconfig,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  real t,
  const RiemannFluxFn& flux,
  const VelFn& vel,
  const StateFn& state,
  const Fields& U,
  const Fields& P,
  const Fields& W,
  const std::vector< int >& srcFlag,
  Fields& R,
  std::vector< std::vector< tk::real > >& riemannDeriv,
  int intsharp=0 );

//! Compute boundary surface flux integrals for a given boundary type for FV
void
bndSurfIntFV(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t rdof,
  const std::vector< std::size_t >& bcconfig,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  real t,
  const RiemannFluxFn& flux,
  const VelFn& vel,
  const StateFn& state,
  const Fields& U,
  const Fields& P,
  const std::vector< int >& srcFlag,
  Fields& R,
  int intsharp );

//! Compute boundary surface flux integrals for a given boundary type for FV
void
bndSurfIntViscousFV(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t rdof,
  const std::vector< std::size_t >& bcconfig,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  real t,
  const StateFn& state,
  const StateFn& gradFn,
  const Fields& U,
  const Fields& P,
  const Fields& T,
  const std::vector< int >& srcFlag,
  Fields& R,
  int intsharp );

template< class ViscousTerms >
void
viscousBoundaryFaceIntDG(
            const std::vector< inciter::EOS >& mat_blk,
            const std::size_t ndof,
            const std::vector< std::size_t >& bcconfig,
            const std::vector< std::size_t >& inpoel,
            const UnsMesh::Coords& coord,
            const inciter::FaceData& fd,
            const Fields& geoFace,
            const Fields& geoElem,
            const Fields& U,
            const Fields& P,
            real t,
            const StateFn& state,
            const StateFn& gradFn,
            Fields& R );

// Compute boundary surface viscous flux integrals for multispecies flow
void
bndSurfIntViscousMultiSpecies(
  std::size_t nspec,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t ndof,
  const std::size_t rdof,
  const std::vector< std::size_t >& bcconfig,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  real t,
  const StateFn& state,
  const StateFn& gradFn,
  const Fields& U,
  const Fields& P,
  Fields& R );

} // tk::

#endif // Boundary_h
