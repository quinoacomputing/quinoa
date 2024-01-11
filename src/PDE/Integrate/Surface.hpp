// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Surface.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing internal surface integrals of a system
     of PDEs in DG methods
  \details   This file contains functionality for computing internal surface
     integrals of a system of PDEs used in discontinuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************
#ifndef Surface_h
#define Surface_h

#include "Basis.hpp"
#include "Types.hpp"
#include "Fields.hpp"
#include "FaceData.hpp"
#include "UnsMesh.hpp"
#include "MultiMatTerms.hpp"
#include "FunctionPrototypes.hpp"
#include "EoS/EOS.hpp"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;
using bcconf_t = kw::sideset::info::expect::type;

//! Compute internal surface flux integrals for DG
void
surfInt( std::size_t nmat,
         const std::vector< inciter::EOS >& mat_blk,
         real t,
         const std::size_t ndof,
         const std::size_t rdof,
         const std::vector< std::size_t >& inpoel,
         const std::vector< std::size_t >& solidx,
         const UnsMesh::Coords& coord,
         const inciter::FaceData& fd,
         const Fields& geoFace,
         const Fields& geoElem,
         const RiemannFluxFn& flux,
         const VelFn& vel,
         const Fields& U,
         const Fields& P,
         const std::vector< std::size_t >& ndofel,
         const tk::real dt,
         Fields& R,
         std::vector< std::vector< tk::real > >& vriem,
         std::vector< std::vector< tk::real > >& xcoord,
         std::vector< std::vector< tk::real > >& riemannDeriv,
         int intcompr=0 );

// Update the rhs by adding surface integration term
void
update_rhs_fa ( ncomp_t ncomp,
                std::size_t nmat,
                const std::size_t ndof,
                const std::size_t ndof_l,
                const std::size_t ndof_r,
                const tk::real wt,
                const std::array< tk::real, 3 >& fn,
                const std::size_t el,
                const std::size_t er,
                const std::vector< tk::real >& fl,
                const std::vector< tk::real >& B_l,
                const std::vector< tk::real >& B_r,
                Fields& R,
                std::vector< std::vector< tk::real > >& riemannDeriv );

// Compute internal surface flux integrals for second order FV
void
surfIntFV(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  real t,
  const std::size_t rdof,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const RiemannFluxFn& flux,
  const VelFn& vel,
  const Fields& U,
  const Fields& P,
  const std::vector< int >& srcFlag,
  Fields& R,
  int intsharp );

} // tk::

#endif // Surface_h
