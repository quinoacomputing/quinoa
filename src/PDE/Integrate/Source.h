// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Source.h
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Functions for computing integrals of an arbitrary source term of a
     system of PDEs in DG methods
  \details   This file contains functionality for computing integrals of an
     arbitrary source term of a system of PDEs used in discontinuous Galerkin
     methods for various orders of numerical representation.
*/
// *****************************************************************************
#ifndef Source_h
#define Source_h

#include "Types.h"
#include "Fields.h"
#include "UnsMesh.h"
#include "FunctionPrototypes.h"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Compute source term integrals for DG(P0)
void
srcIntP0( ncomp_t system,
          ncomp_t ncomp,
          ncomp_t offset,
          real t,
          const Fields& geoElem,
          const SrcFn& src,
          Fields& R );

//! Compute source term integrals for DG(P1)
void
srcIntP1( ncomp_t system,
          ncomp_t ncomp,
          ncomp_t offset,
          real t,
          const std::vector< std::size_t >& inpoel,
          const UnsMesh::Coords& coord,
          const Fields& geoElem,
          const SrcFn& src,
          Fields& R );

//! Compute source term integrals for DG(P2)
void
srcIntP2( ncomp_t system,
          ncomp_t ncomp,
          ncomp_t offset,
          real t,
          const std::vector< std::size_t >& inpoel,
          const UnsMesh::Coords& coord,
          const Fields& geoElem,
          const SrcFn& src,
          Fields& R );

} // tk::

#endif // Source_h
