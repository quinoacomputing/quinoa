// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Surface.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
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
#include "FunctionPrototypes.hpp"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;
using bcconf_t = kw::sideset::info::expect::type;

//! Compute internal surface flux integrals for DG
void
surfInt( ncomp_t system,
         std::size_t nmat,
         ncomp_t offset,
         const std::size_t ndof,
         const std::size_t rdof,
         const std::vector< std::size_t >& inpoel,
         const UnsMesh::Coords& coord,
         const inciter::FaceData& fd,
         const Fields& geoFace,
         const RiemannFluxFn& flux,
         const VelFn& vel,
         const Fields& U,
         const Fields& P,
         const std::vector< std::size_t >& ndofel,
         Fields& R,
         std::vector< std::vector< tk::real > >& riemannDeriv );

// Update the rhs by adding surface integration term
void
update_rhs_fa ( ncomp_t ncomp,
                std::size_t nmat,
                ncomp_t offset,
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

} // tk::

#endif // Surface_h
