// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Source.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing integrals of an arbitrary source term of a
     for single-material compressible flow, CompFlow using DG methods
  \details   This file contains functionality for computing integrals of an
     arbitrary source term for single-material compressible flow, CompFlow with
     discontinuous Galerkin methods for various orders of numerical
     representation.
*/
// *****************************************************************************
#ifndef Source_h
#define Source_h

#include "Basis.hpp"
#include "Types.hpp"
#include "Fields.hpp"
#include "UnsMesh.hpp"
#include "FunctionPrototypes.hpp"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Compute source term integrals for DG
void
srcInt( ncomp_t system,
        ncomp_t offset,
        const std::vector< inciter::EoS_Base* >& mat_blk,
        real t,
        const std::size_t ndof,
        const std::size_t nelem,
        const std::vector< std::size_t >& inpoel,
        const UnsMesh::Coords& coord,
        const Fields& geoElem,
        const SrcFn& src,
        const std::vector< std::size_t >& ndofel,
        Fields& R,
        std::size_t nmat=1 );

//! Update the rhs by adding the source term integrals
void
update_rhs( ncomp_t offset,
            const std::size_t ndof,
            const std::size_t ndof_el,
            const tk::real wt,
            const std::size_t e,
            const std::vector< tk::real >& B,
            const std::vector< tk::real >& s,
            Fields& R );

} // tk::

#endif // Source_h
