// *****************************************************************************
/*!
  \file      src/PDE/Limiter.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Limiters for discontinous Galerkin methods
  \details   This file contains functions that provide limiter function
    calculations for maintaining monotonicity near solution discontinuities
    for the DG discretization.
*/
// *****************************************************************************
#ifndef Limiter_h
#define Limiter_h

#include "Types.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Fields.h"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

//! Weighted Essentially Non-Oscillatory (WENO) limiter for DGP1
void
WENO_P1( const std::vector< int >& esuel,
         inciter::ncomp_t offset,
         const tk::Fields& U,
         tk::Fields& limFunc );

} // inciter::

#endif // Limiter_h
