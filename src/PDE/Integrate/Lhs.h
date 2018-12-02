// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Lhs.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for computing the left hand side of system of PDEs in DG
     methods
  \details   This file contains functionality for computing the left hand side
     of PDEs used in discontinuous Galerkin methods for various orders of
     numerical representation.
*/
// *****************************************************************************
#ifndef Lhs_h
#define Lhs_h

#include "Types.h"
#include "Fields.h"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Compute the left hand side block-diagnoal mass matrix
void
lhs( ncomp_t ncomp, ncomp_t offset, const tk::Fields& geoElem, tk::Fields& l );

} // tk::

#endif // Lhs_h
