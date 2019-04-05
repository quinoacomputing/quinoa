// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Volume.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
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

#include "Basis.h"
#include "Types.h"
#include "Fields.h"
#include "UnsMesh.h"
#include "FunctionPrototypes.h"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Compute volume integrals for DG
void
volInt( ncomp_t system,
        ncomp_t ncomp,
        ncomp_t offset,
        const std::vector< std::size_t >& inpoel,
        const UnsMesh::Coords& coord,
        const Fields& geoElem,
        const FluxFn& flux,
        const VelFn& vel,
        const Fields& U,
        const Fields& limFunc,
        const std::vector< std::size_t >& ndofel,
        Fields& R );

//! Update the rhs by adding the source term integrals
void
update_rhs( ncomp_t ncomp,
            ncomp_t offset,
            const std::size_t ndof,
            const std::size_t ndof_e,
            const tk::real wt,
            const std::size_t e,
            const std::array< std::vector<tk::real>, 3 >& dBdx,
            const std::vector< std::array< tk::real, 3 > >& fl,
            Fields& R );

} // tk::

#endif // Volume_h
