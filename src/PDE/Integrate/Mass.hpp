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

#include "Fields.hpp"

namespace tk {

using ncomp_t = tk::ncomp_t;

//! Compute lumped mass matrix for CG
tk::Fields
lump( ncomp_t ncomp,
      const std::array< std::vector< tk::real >, 3 >& coord,
      const std::vector< std::size_t >& inpoel );

} // tk::

#endif // Mass_h
