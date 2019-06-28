// *****************************************************************************
/*!
  \file      src/PDE/Reconstruction.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Reconstruction for reconstructed discontinuous Galerkin methods
  \details   This file contains functions that reconstruct an "n"th order
    polynomial to an "n+1"th order polynomial using a least-squares
    reconstruction procedure.
*/
// *****************************************************************************
#ifndef Reconstruction_h
#define Reconstruction_h

#include "Fields.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

//! Least-squares reconstruction for rDG(P0P1)
void
leastSquares_P0P1( const std::vector< int >& esuel,
                   inciter::ncomp_t offset,
                   const tk::Fields& geoElem,
                   tk::Fields& U );

} // inciter::

#endif // Reconstruction_h
