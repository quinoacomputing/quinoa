// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Volume.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing volume integrals for a system of PDEs in DG
     methods
  \details   This file contains functionality for computing volume integrals for
     a system of PDEs used in discontinuous Galerkin methods for various orders
     of numerical representation.
*/
// *****************************************************************************
#ifndef Volume_h
#define Volume_h

#include "Basis.hpp"
#include "Fields.hpp"
#include "UnsMesh.hpp"
#include "FunctionPrototypes.hpp"
#include "EoS/EOS.hpp"

namespace tk {

using ncomp_t = tk::ncomp_t;

//! Compute volume integrals for DG
void
volInt( std::size_t nmat,
        real t,
        const std::vector< inciter::EOS >& mat_blk,
        const std::size_t ndof,
        const std::size_t rdof,
        const std::size_t nelem,
        const std::vector< std::size_t >& inpoel,
        const UnsMesh::Coords& coord,
        const Fields& geoElem,
        const FluxFn& flux,
        const VelFn& vel,
        const SrcFn& src,
        const Fields& U,
        const Fields& P,
        const std::vector< std::size_t >& ndofel,
        Fields& R,
        int intsharp=0 );

//! Update the rhs by adding the source term integrals
void
update_rhs( ncomp_t ncomp,
            const std::size_t ndof,
            const std::size_t ndof_el,
            const tk::real wt,
            const std::size_t e,
            const std::array< std::vector<tk::real>, 3 >& dBdx,
            const std::vector< std::array< tk::real, 3 > >& fl,
            Fields& R );

//! Update the rhs by adding the source term integrals
void
update_rhs_src( const std::size_t ndof,
                const std::size_t ndof_el,
                const tk::real wt,
                const std::size_t e,
                const std::vector< tk::real >& B,
                const std::vector< tk::real >& s,
                Fields& R );

//! Compute volume integrals for const-order DG (not p-adaptive)
void
volInt_constP(
  std::size_t nmat,
  real t,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t ndof,
  const std::size_t rdof,
  const std::size_t nelem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const Fields& geoElem,
  const FluxFn& flux,
  const VelFn& vel,
  const SrcFn& src,
  const Fields& U,
  const Fields& P,
  Fields& R,
  int intsharp=0 );

//! Compute source term integrals for FV
void
srcIntFV( const std::vector< inciter::EOS >& mat_blk,
          real t,
          const std::size_t nelem,
          const Fields& geoElem,
          const SrcFn& src,
          Fields& R,
          std::size_t nmat );

} // tk::

#endif // Volume_h
