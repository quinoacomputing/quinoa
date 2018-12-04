// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Boundary.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for computing boundary surface integrals of a system
     of PDEs in DG methods
  \details   This file contains functionality for computing boundary surface
     integrals of a system of PDEs used in discontinuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************
#ifndef Boundary_h
#define Boundary_h

#include "Types.h"
#include "Fields.h"
#include "FaceData.h"
#include "UnsMesh.h"
#include "FunctionPrototypes.h"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;
using bcconf_t = kw::sideset::info::expect::type;

//! Compute boundary surface integral for a number of faces for DG(P0)
void
bndSurfIntP0( ncomp_t system,
              ncomp_t ncomp,
              ncomp_t offset,
              const std::vector< std::size_t >& faces,
              const std::vector< int >& esuf,
              const tk::Fields& geoFace,
              tk::real t,
              const RiemannFluxFn& flux,
              const VelFn& vel,
              const StateFn& state,
              const tk::Fields& U,
              tk::Fields& R );

//! Compute boundary surface flux integrals for a given boundary type for DG(P0)
void
sidesetIntP0( ncomp_t system,
              ncomp_t ncomp,
              ncomp_t offset,
              const std::vector< bcconf_t >& bcconfig,
              const std::map< int, std::vector< std::size_t > >& bface,
              const std::vector< int >& esuf,
              const tk::Fields& geoFace,
              tk::real t,
              const RiemannFluxFn& flux,
              const VelFn& vel,
              const StateFn& state,
              const tk::Fields& U,
              tk::Fields& R );

//! Compute boundary surface integral for a number of faces for DG(P1)
void
bndSurfIntP1( ncomp_t system,
              ncomp_t ncomp,
              ncomp_t offset,
              const std::vector< std::size_t >& faces,
              const std::vector< int >& esuf,
              const tk::Fields& geoFace,
              const std::vector< std::size_t >& inpoel,
              const std::vector< std::size_t >& inpofa,
              const tk::UnsMesh::Coords& coord,
              tk::real t,
              const RiemannFluxFn& flux,
              const VelFn& vel,
              const StateFn& state,
              const tk::Fields& U,
              const tk::Fields& limFunc,
              tk::Fields& R );

//! Compute boundary surface flux integrals for a given boundary type for DG(P1)
void
sidesetIntP1( ncomp_t system,
              ncomp_t ncomp,
              ncomp_t offset,
              const std::vector< bcconf_t >& bcconfig,
              const std::map< int, std::vector< std::size_t > >& bface,
              const std::vector< int >& esuf,
              const tk::Fields& geoFace,
              const std::vector< std::size_t >& inpoel,
              const std::vector< std::size_t >& inpofa,
              const tk::UnsMesh::Coords& coord,
              tk::real t,
              const RiemannFluxFn& flux,
              const VelFn& vel,
              const StateFn& state,
              const tk::Fields& U,
              const tk::Fields& limFunc,
              tk::Fields& R );

} // tk::

#endif // Boundary_h
