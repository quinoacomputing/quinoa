// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Surf.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for computing internal surface integrals of a system
     of PDEs in DG methods
  \details   This file contains functionality for computing internal surface
     integrals of a system of PDEs used in discontinuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************
#ifndef Surf_h
#define Surf_h

#include "Types.h"
#include "Fields.h"
#include "FaceData.h"
#include "UnsMesh.h"
#include "FunctionPrototypes.h"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Compute internal surface flux integrals for DG(P0)
void
surfIntP0( ncomp_t system,
           ncomp_t ncomp,
           ncomp_t offset,
           const inciter::FaceData& fd,
           const Fields& geoFace,
           const FluxFn& flux,
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
           const FluxFn& flux,
           const VelFn& vel,
           const Fields& U,
           const Fields& limFunc,
           Fields& R );

} // tk::

#endif // Surf_h
