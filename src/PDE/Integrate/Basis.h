// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Basis.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for computing the Dubiner basis functions in DG methods
  \details   This file contains functionality for computing the basis functions
     and relating coordinates transformation functions used in discontinuous
     Galerkin methods for variaous orders of numerical representation. The basis
     functions chosen for the DG method are the Dubiner basis, which are Legendre
     polynomials modified for tetrahedra, which are defined only on the reference/master
     tetrahedron.
  \see [1] https://doi.org/10.1007/BF01060030
  \see [2] https://doi.org/10.1093/imamat/hxh111
*/
// *****************************************************************************

#ifndef Basis_h
#define Basis_h

#include "Types.h"
#include "Vector.h"
#include "Fields.h"
#include "FaceData.h"
#include "UnsMesh.h"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;
using bcconf_t = kw::sideset::info::expect::type;

//! Compute the determinant of Jacobian matrix for the tetrahedron element
tk::real
eval_det ( const std::size_t e,
           const std::vector< tk::real >& cx,
           const std::vector< tk::real >& cy,
           const std::vector< tk::real >& cz,
           const std::vector< std::size_t >& inpoel,
           std::array< std::array< tk::real, 3>, 4 >& coordel );

//! Compute the coordinates of quadrature points in physical space
std::array< tk::real, 3>
eval_gp ( const std::size_t igp,
          const std::array< std::array< tk::real, 3>, 3 >& coordfa,
          const std::array< std::vector< tk::real >, 2 >& coordgp );

//! Compute the Dubiner basis functions
std::vector< tk::real >
eval_basis( const std::size_t ndof,
            const std::array< std::array< tk::real, 3>, 4 >& coordel,
            const tk::real detT,
            const std::array < tk::real, 3 >& gp );

//! Compute the state variables for the tetrahedron element
std::vector< tk::real >
eval_state ( ncomp_t ncomp,
             ncomp_t offset,
             const std::size_t ndof,
             const std::size_t e,
             const Fields& U,
             const Fields& limFunc,
             const std::vector< tk::real >& B );

} // tk::

#endif // Basis_h
