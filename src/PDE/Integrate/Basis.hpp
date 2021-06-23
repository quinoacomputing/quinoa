// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Basis.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

#include "Types.hpp"
#include "Vector.hpp"
#include "Fields.hpp"
#include "FaceData.hpp"
#include "UnsMesh.hpp"
#include "Integrate/Quadrature.hpp"
#include "../MultiMat/MultiMatIndexing.hpp"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;
using bcconf_t = kw::sideset::info::expect::type;

//! Compute the coordinates of quadrature points for face integral
std::array< tk::real, 3 >
eval_gp ( const std::size_t igp,
          const std::array< std::array< tk::real, 3>, 3 >& coordfa,
          const std::array< std::vector< tk::real >, 2 >& coordgp );

//! Compute the coordinates of quadrature points for volume integral
std::array< tk::real, 3>
eval_gp ( const std::size_t igp,
          const std::array< std::array< tk::real, 3>, 4 >& coord,
          const std::array< std::vector< tk::real >, 3 >& coordgp );

//! Compute the derivatives of basis function for DG(P1)
std::array< std::vector<tk::real>, 3 >
eval_dBdx_p1( const std::size_t ndof,
              const std::array< std::array< tk::real, 3 >, 3 >& jacInv );

//! Compute the derivatives of basis function for DG(P2)
void
eval_dBdx_p2( const std::size_t igp,
              const std::array< std::vector< tk::real >, 3 >& coordgp,
              const std::array< std::array< tk::real, 3 >, 3 >& jacInv,
              std::array< std::vector<tk::real>, 3 >& dBdx );

//! Compute the Dubiner basis functions
std::vector< tk::real >
eval_basis( const std::size_t ndof,
            const tk::real xi,
            const tk::real eta,
            const tk::real zeta );

//! Compute the state variables for the tetrahedron element
std::vector< tk::real >
eval_state ( ncomp_t ncomp,
             ncomp_t offset,
             const std::size_t ndof,
             const std::size_t ndof_el,
             const std::size_t e,
             const Fields& U,
             const std::vector< tk::real >& B );

//! Compute the derivatives of basis function in physical domain
void
evaldBdx_p2(  const std::vector< tk::real >& coordgp,
              const std::array< std::array< tk::real, 3 >, 3 >& jacInv,
              std::array< std::vector<tk::real>, 3 >& dBdx );

//! Transform the solution with Dubiner basis to the solution with Taylor basis
void TransformBasis( ncomp_t ncomp,
                     ncomp_t offset,
                     const std::size_t e,
                     const std::size_t ndof,
                     const tk::Fields& U,
                     const std::vector< std::size_t >& inpoel,
                     const tk::UnsMesh::Coords& coord,
                     std::vector< std::vector< tk::real > >& unk );

//! Convert the solution with Taylor basis to the solution with Dubiner basis
void InverseBasis( ncomp_t ncomp,
                   ncomp_t offset,
                   std::size_t e,
                   std::size_t ndof,
                   const std::vector< std::size_t >& inpoel,
                   const tk::UnsMesh::Coords& coord,
                   const tk::Fields& geoElem,
                   tk::Fields& U,
                   std::vector< std::vector< tk::real > >& unk );

//! Evaluate the Taylor basis at points
std::vector< tk::real >
eval_TaylorBasis( const std::size_t ndof,
                  const std::array< tk::real, 3 >& x,
                  const std::array< tk::real, 3 >& x_c,
                  const std::array< std::array< tk::real, 3>, 4 >& coordel );
} // tk::

#endif // Basis_h
