// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Mass.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing the mass matrix for a system of PDEs
  \details   This file contains functionality for computing the mass matrix for
     a system of PDEs used in discontinuous and continuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************
#ifndef Mass_h
#define Mass_h

#include "Types.hpp"
#include "Fields.hpp"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Compute the block-diagnoal mass matrix for DG
void
mass( ncomp_t ncomp, ncomp_t offset, ncomp_t ndof, const Fields& geoElem,
      Fields& l );

//! Compute lumped mass matrix for CG
tk::Fields
lump( ncomp_t ncomp,
      const std::array< std::vector< tk::real >, 3 >& coord,
      const std::vector< std::size_t >& inpoel );

//! Compute the diagonal mass matrix for DG with Dubiner basis functions
std::vector< tk::real >
massMatrixDubiner( std::size_t dof,
  tk::real vol );

} // tk::

#endif // Mass_h
