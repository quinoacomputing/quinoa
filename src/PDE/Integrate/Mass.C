// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Mass.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for computing the mass matrix for a system of PDEs in DG
     methods
  \details   This file contains functionality for computing the mass matrix for
     a system of PDEs used in discontinuous Galerkin methods for various orders
     of numerical representation.
*/
// *****************************************************************************

#include "Mass.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

void
tk::mass( ncomp_t ncomp,
          ncomp_t offset,
          const Fields& geoElem,
          Fields& l )
// *****************************************************************************
//  Compute the block-diagonal mass matrix for DG
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] geoElem Element geometry array
//! \param[in,out] l Block diagonal mass matrix
// *****************************************************************************
{
  Assert( geoElem.nunk() == l.nunk(), "Size mismatch" );
  const auto nelem = geoElem.nunk();
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();

  // Compute LHS for DG(P0)
  for (std::size_t e=0; e<nelem; ++e)
    for (ncomp_t c=0; c<ncomp; ++c)
      l(e, c*ndof, offset) = geoElem(e,0,0);

  // Augment LHS for DG(P1)
  if (ndof > 1) {
    for (std::size_t e=0; e<nelem; ++e) {
      for (ncomp_t c=0; c<ncomp; ++c) {
        const auto mark = c * ndof;
        l(e, mark+1, offset) = geoElem(e,0,0) / 10.0;
        l(e, mark+2, offset) = geoElem(e,0,0) * 3.0/10.0;
        l(e, mark+3, offset) = geoElem(e,0,0) * 3.0/5.0;
      }
    }
  }
}
