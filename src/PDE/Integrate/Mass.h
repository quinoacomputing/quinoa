// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Mass.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for computing the mass matrix for a system of PDEs in DG
     methods
  \details   This file contains functionality for computing the mass matrix for
     a system of PDEs used in discontinuous Galerkin methods for various orders
     of numerical representation.
*/
// *****************************************************************************
#ifndef Mass_h
#define Mass_h

#include "Types.h"
#include "Fields.h"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Compute the block-diagnoal mass matrix for DG
void
mass( ncomp_t ncomp, ncomp_t offset, const Fields& geoElem, Fields& l );

} // tk::

#endif // Mass_h
