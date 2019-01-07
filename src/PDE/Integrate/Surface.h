// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Surface.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for computing internal surface integrals of a system
     of PDEs in DG methods
  \details   This file contains functionality for computing internal surface
     integrals of a system of PDEs used in discontinuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************
#ifndef Surface_h
#define Surface_h

#include "Types.h"
#include "Fields.h"
#include "FaceData.h"
#include "UnsMesh.h"
#include "FunctionPrototypes.h"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;
using bcconf_t = kw::sideset::info::expect::type;

//! Compute internal surface flux integrals for DG(P0)
void
surfIntP0( ncomp_t system,
           ncomp_t ncomp,
           ncomp_t offset,
           const inciter::FaceData& fd,
           const Fields& geoFace,
           const RiemannFluxFn& flux,
           const VelFn& vel,
           const Fields& U,
           Fields& R );

//! Compute internal surface flux integrals for DG(P1)
void
surfIntP1( ncomp_t system,
           ncomp_t ncomp,
           ncomp_t offset,
           const std::vector< std::size_t >& inpoel,
           const UnsMesh::Coords& coord,
           const inciter::FaceData& fd,
           const Fields& geoFace,
           const RiemannFluxFn& flux,
           const VelFn& vel,
           const Fields& U,
           const Fields& limFunc,
           Fields& R );

//! Compute internal surface flux integrals for DG(P2)
void
surfIntP2( ncomp_t system,
           ncomp_t ncomp,
           ncomp_t offset,
           const std::vector< std::size_t >& inpoel,
           const UnsMesh::Coords& coord,
           const inciter::FaceData& fd,
           const Fields& geoFace,
           const RiemannFluxFn& flux,
           const VelFn& vel,
           const Fields& U,
           Fields& R );

} // tk::

#endif // Surface_h
