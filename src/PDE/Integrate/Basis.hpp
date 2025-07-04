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

#include "Vector.hpp"
#include "Fields.hpp"
#include "FaceData.hpp"
#include "UnsMesh.hpp"
#include "Integrate/Quadrature.hpp"
#include "../MultiMat/MultiMatIndexing.hpp"
#include "Kokkos_Core.hpp"

using execution_space = Kokkos::Serial;
using memory_space = Kokkos::HostSpace;

namespace tk {

using ncomp_t = tk::ncomp_t;

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

//! Kokkos version of eval_gp volume integral
KOKKOS_INLINE_FUNCTION Kokkos::Array<real, 3>
tk::eval_gp ( const std::size_t igp,
              const Kokkos::Array<Kokkos::Array<real, 3>, 4>& coord,
              Kokkos::View<const double**, memory_space> coordgp );

//! Compute the derivatives of Dubiner basis wrt. reference coordinates
std::array< std::vector< tk::real >, 3 >
eval_dBdxi( const std::size_t ndof,
            const std::array< tk::real, 3 >& coordgp );

//! Compute the derivatives of basis function for DG(P1)
std::array< std::vector<tk::real>, 3 >
eval_dBdx_p1( const std::size_t ndof,
              const std::array< std::array< tk::real, 3 >, 3 >& jacInv );

//! Kokkos version of eval_dBdx_p1
KOKKOS_INLINE_FUNCTION auto
tk::eval_dBdx_p1( const std::size_t ndof,
                const Kokkos::Array<Kokkos::Array<real, 3>, 3>& jacInv, 
                Kokkos::View<real**, memory_space> dBdx);

//! Compute the derivatives of basis function for DG(P2)
void
eval_dBdx_p2( const std::size_t igp,
              const std::array< std::vector< tk::real >, 3 >& coordgp,
              const std::array< std::array< tk::real, 3 >, 3 >& jacInv,
              std::array< std::vector<tk::real>, 3 >& dBdx );

//! Kokkos version of eval_dBdx_p2
KOKKOS_INLINE_FUNCTION 
void tk::eval_dBdx_p2( const std::size_t igp,
              Kokkos::View<real**, memory_space> coordgp,
              const Kokkos::Array<Kokkos::Array<real, 3>, 3>& jacInv,
              Kokkos::View<real**, memory_space> dBdx);

//! Compute the Dubiner basis functions
std::vector< tk::real >
eval_basis( const std::size_t ndof,
            const tk::real xi,
            const tk::real eta,
            const tk::real zeta );
//! Kokks version of eval_basis

KOKKOS_INLINE_FUNCTION 
auto tk::eval_basis( const std::size_t ndof,
                const tk::real xi,
                const tk::real eta,
                const tk::real zeta, 
               const tk:real );

//! Compute the state variables for the tetrahedron element
std::vector< tk::real >
eval_state ( ncomp_t ncomp,
             const std::size_t ndof,
             const std::size_t ndof_el,
             const std::size_t e,
             const Fields& U,
             const std::vector< tk::real >& B );
  
//! Kokkos versioon of eval_state
template <typename BasisType>
KOKKOS_INLINE_FUNCTION 
void tk::eval_state ( ncomp_t ncomp,
                 const std::size_t ndof,
                 const std::size_t ndof_el,
                 const std::size_t e, size_t m_nprop,
                 Kokkos::View<const real*, memory_space> U,
                 BasisType B, 
                 Kokkos::View<real*, memory_space> state);

//! Transform the solution with Dubiner basis to the solution with Taylor basis
std::vector< std::vector< tk::real > >
DubinerToTaylor( ncomp_t ncomp,
                 const std::size_t e,
                 const std::size_t ndof,
                 const tk::Fields& U,
                 const std::vector< std::size_t >& inpoel,
                 const tk::UnsMesh::Coords& coord );

//! Convert the solution with Taylor basis to the solution with Dubiner basis
void
TaylorToDubiner( ncomp_t ncomp,
                 std::size_t e,
                 std::size_t ndof,
                 const std::vector< std::size_t >& inpoel,
                 const tk::UnsMesh::Coords& coord,
                 const tk::Fields& geoElem,
                 std::vector< std::vector< tk::real > >& unk );

//! Evaluate the Taylor basis at points
std::vector< tk::real >
eval_TaylorBasis( const std::size_t ndof,
                  const std::array< tk::real, 3 >& x,
                  const std::array< tk::real, 3 >& x_c,
                  const std::array< std::array< tk::real, 3>, 4 >& coordel );

// Reference element Taylor basis functions
// ----------------------------------------

//! Transform the solution from Dubiner basis to Taylor basis
std::vector< std::vector< tk::real > >
DubinerToTaylorRefEl( ncomp_t ncomp,
  const std::size_t e,
  const std::size_t ndof,
  const std::size_t ndof_el,
  const std::vector< std::vector< tk::real > >& mtInv,
  const tk::Fields& U );

//! Transform the solution from Taylor to Dubiner basis
void
TaylorToDubinerRefEl( ncomp_t ncomp,
  const std::size_t ndof,
  std::vector< std::vector< tk::real > >& unk );

//! Evaluate the Taylor basis at a point in the reference element
std::vector< tk::real >
eval_TaylorBasisRefEl( std::size_t ndof, tk::real x, tk::real y,
  tk::real z );

//! Obtain inverse mass matrix for Taylor basis in reference element
std::vector< std::vector< tk::real > >
invMassMatTaylorRefEl( std::size_t dof );

//! Obtain mass matrix for Taylor basis in reference element
std::vector< std::vector< tk::real > >
massMatrixTaylorRefEl(std::size_t dof);

} // tk::

#endif // Basis_h
