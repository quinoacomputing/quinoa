// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Volume.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for computing volume integrals for a system of PDEs in DG
     methods
  \details   This file contains functionality for computing volume integrals for
     a system of PDEs used in discontinuous Galerkin methods for various orders
     of numerical representation.
*/
// *****************************************************************************
#ifndef Volume_h
#define Volume_h

#include "Types.h"
#include "Fields.h"
#include "UnsMesh.h"
#include "FunctionPrototypes.h"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Compute volume integrals for DG(P1)
void
volIntP1( ncomp_t system,
          ncomp_t ncomp,
          ncomp_t offset,
          const std::vector< std::size_t >& inpoel,
          const UnsMesh::Coords& coord,
          const Fields& geoElem,
          const FluxFn& flux,
          const VelFn& vel,
          const Fields& U,
          const Fields& limFunc,
          Fields& R );

//! Compute volume integrals for DG(P2)
void
volIntP2( ncomp_t system,
          ncomp_t ncomp,
          ncomp_t offset,
          const std::vector< std::size_t >& inpoel,
          const UnsMesh::Coords& coord,
          const Fields& geoElem,
          const FluxFn& flux,
          const VelFn& vel,
          const Fields& U,
          Fields& R );

} // tk::

#endif // Volume_h
