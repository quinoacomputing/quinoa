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

using execution_space = Kokkos::DefaultExecutionSpace;
using memory_space = Kokkos::DefaultExecutionSpace::memory_space;

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
KOKKOS_INLINE_FUNCTION Kokkos::Array<tk::real, 3>
eval_gp ( const std::size_t igp,
              const Kokkos::Array<Kokkos::Array<tk::real, 3>, 4>& coord,
              Kokkos::View<const tk::real**, memory_space> coordgp )
// *****************************************************************************
//  Compute the coordinates of quadrature points for volume integral in
//  physical space
//! \param[in] igp Index of quadrature points
//! \param[in] coord Array of nodal coordinates for tetrahedron element
//! \param[in] coordgp Array of coordinates for quadrature points in reference space
//! \return Array of coordinates for quadrature points in physical space
// *****************************************************************************
{
  // Barycentric coordinates for the tetradedron element
  auto shp1 = 1.0 - coordgp(0, igp) - coordgp(1, igp) - coordgp(2, igp);
  auto shp2 = coordgp(0, igp);
  auto shp3 = coordgp(1, igp);
  auto shp4 = coordgp(2, igp);

  // Transformation of the quadrature point from the reference/master
  // element to physical space, to obtain its physical (x,y,z) coordinates.
  return {{
   coord[0][0]*shp1 + coord[1][0]*shp2 + coord[2][0]*shp3 + coord[3][0]*shp4,
   coord[0][1]*shp1 + coord[1][1]*shp2 + coord[2][1]*shp3 + coord[3][1]*shp4,
   coord[0][2]*shp1 + coord[1][2]*shp2 + coord[2][2]*shp3 + coord[3][2]*shp4 }};
}

//! Compute the derivatives of Dubiner basis wrt. reference coordinates
std::array< std::vector< tk::real >, 3 >
eval_dBdxi( const std::size_t ndof,
            const std::array< tk::real, 3 >& coordgp );

//! Compute the derivatives of basis function for DG(P1)
std::array< std::vector<tk::real>, 3 >
eval_dBdx_p1( const std::size_t ndof,
              const std::array< std::array< tk::real, 3 >, 3 >& jacInv );

//! Kokkos version of eval_dBdx_p1
KOKKOS_INLINE_FUNCTION void
eval_dBdx_p1( const std::size_t ndof,
                  const Kokkos::Array<Kokkos::Array<real, 3>, 3>& jacInv, 
                  Kokkos::View<real**, memory_space> dBdx)
// *****************************************************************************
//  Compute the derivatives of basis functions for DG(P1)
//! \param[in] ndof Number of degrees of freedom
//! \param[in] jacInv Array of the inverse of Jacobian
//! \return Array of the derivatives of basis functions
// *****************************************************************************
{
  // The derivatives of the basis functions dB/dx are easily calculated
  // via a transformation to the reference space as,
  // dB/dx = dB/dxi . dxi/dx,
  // where, x = (x,y,z) are the physical coordinates, and
  //        xi = (xi, eta, zeta) are the reference coordinates.
  // The matrix dxi/dx is the inverse of the Jacobian of transformation
  // and the matrix vector product has to be calculated. This follows.


  auto db2dxi1 = 2.0;
  auto db2dxi2 = 1.0;
  auto db2dxi3 = 1.0;

  auto db3dxi1 = 0.0;
  auto db3dxi2 = 3.0;
  auto db3dxi3 = 1.0;

  auto db4dxi1 = 0.0;
  auto db4dxi2 = 0.0;
  auto db4dxi3 = 4.0;
  
  dBdx(0, 0) = 0.0;

  dBdx(0, 1) =  db2dxi1 * jacInv[0][0]
              + db2dxi2 * jacInv[1][0]
              + db2dxi3 * jacInv[2][0];

  dBdx(1, 1) =  db2dxi1 * jacInv[0][1]
              + db2dxi2 * jacInv[1][1]
              + db2dxi3 * jacInv[2][1];

  dBdx(2, 1) =  db2dxi1 * jacInv[0][2]
              + db2dxi2 * jacInv[1][2]
              + db2dxi3 * jacInv[2][2];

  dBdx(0, 2) =  db3dxi1 * jacInv[0][0]
              + db3dxi2 * jacInv[1][0]
              + db3dxi3 * jacInv[2][0];

  dBdx(1, 2) =  db3dxi1 * jacInv[0][1]
              + db3dxi2 * jacInv[1][1]
              + db3dxi3 * jacInv[2][1];

  dBdx(2, 2) =  db3dxi1 * jacInv[0][2]
              + db3dxi2 * jacInv[1][2]
              + db3dxi3 * jacInv[2][2];

  dBdx(0, 3) =  db4dxi1 * jacInv[0][0]
              + db4dxi2 * jacInv[1][0]
              + db4dxi3 * jacInv[2][0];

  dBdx(1, 3) =  db4dxi1 * jacInv[0][1]
              + db4dxi2 * jacInv[1][1]
              + db4dxi3 * jacInv[2][1];

  dBdx(2, 3) =  db4dxi1 * jacInv[0][2]
              + db4dxi2 * jacInv[1][2]
              + db4dxi3 * jacInv[2][2];
};

//! Compute the derivatives of basis function for DG(P2)
void
eval_dBdx_p2( const std::size_t igp,
              const std::array< std::vector< tk::real >, 3 >& coordgp,
              const std::array< std::array< tk::real, 3 >, 3 >& jacInv,
              std::array< std::vector<tk::real>, 3 >& dBdx );

//! Kokkos version of eval_dBdx_p2
KOKKOS_INLINE_FUNCTION 
void eval_dBdx_p2( const std::size_t igp,
                  Kokkos::View<real**, memory_space> coordgp,
                  const Kokkos::Array<Kokkos::Array<real, 3>, 3>& jacInv,
                  Kokkos::View<real**, memory_space> dBdx)
// *****************************************************************************
//  Compute the derivatives of Dubiner basis function for DG(P2)
//! \param[in] igp Index of quadrature points
//! \param[in] coordgp Gauss point coordinates for tetrahedron element
//! \param[in] jacInv Array of the inverse of Jacobian
//! \param[in,out] dBdx Array of the derivatives of basis function
// *****************************************************************************
{
  auto db5dxi1 = 12.0 * coordgp(0, igp) + 6.0 * coordgp(1, igp)
               +  6.0 * coordgp(2, igp) - 6.0;
  auto db5dxi2 =  6.0 * coordgp(0, igp) + 2.0 * coordgp(1, igp)
               +  2.0 * coordgp(2, igp) - 2.0;
  auto db5dxi3 =  6.0 * coordgp(0, igp) + 2.0 * coordgp(1, igp)
               +  2.0 * coordgp(2, igp) - 2.0;

  auto db6dxi1 = 10.0 * coordgp(1, igp) +  2.0 * coordgp(2, igp) - 2.0;
  auto db6dxi2 = 10.0 * coordgp(0, igp) + 10.0 * coordgp(1, igp)
               +  6.0 * coordgp(2, igp) - 6.0;
  auto db6dxi3 =  2.0 * coordgp(0, igp) +  6.0 * coordgp(1, igp)
               +  2.0 * coordgp(2, igp) - 2.0;

  auto db7dxi1 = 12.0 * coordgp(2, igp) - 2.0;
  auto db7dxi2 =  6.0 * coordgp(2, igp) - 1.0;
  auto db7dxi3 = 12.0 * coordgp(0, igp) + 6.0 * coordgp(1, igp)
               + 12.0 * coordgp(2, igp) - 7.0;

  auto db8dxi1 =  0;
  auto db8dxi2 = 20.0 * coordgp(1, igp) + 8.0 * coordgp(2, igp) - 8.0;
  auto db8dxi3 =  8.0 * coordgp(1, igp) + 2.0 * coordgp(2, igp) - 2.0;

  auto db9dxi1 =  0;
  auto db9dxi2 = 18.0 * coordgp(2, igp) -  3.0;
  auto db9dxi3 = 18.0 * coordgp(1, igp) + 12.0 * coordgp(2, igp) - 7.0;

  auto db10dxi1 =  0;
  auto db10dxi2 =  0;
  auto db10dxi3 = 30.0 * coordgp(2, igp) - 10.0;

  dBdx(0, 4) =  db5dxi1 * jacInv[0][0]
              + db5dxi2 * jacInv[1][0]
              + db5dxi3 * jacInv[2][0];

  dBdx(1, 4) =  db5dxi1 * jacInv[0][1]
              + db5dxi2 * jacInv[1][1]
              + db5dxi3 * jacInv[2][1];

  dBdx(2, 4) =  db5dxi1 * jacInv[0][2]
              + db5dxi2 * jacInv[1][2]
              + db5dxi3 * jacInv[2][2];

  dBdx(0, 5) =  db6dxi1 * jacInv[0][0]
              + db6dxi2 * jacInv[1][0]
              + db6dxi3 * jacInv[2][0];

  dBdx(1, 5) =  db6dxi1 * jacInv[0][1]
              + db6dxi2 * jacInv[1][1]
              + db6dxi3 * jacInv[2][1];

  dBdx(2, 5) =  db6dxi1 * jacInv[0][2]
              + db6dxi2 * jacInv[1][2]
              + db6dxi3 * jacInv[2][2];

  dBdx(0, 6) =  db7dxi1 * jacInv[0][0]
              + db7dxi2 * jacInv[1][0]
              + db7dxi3 * jacInv[2][0];

  dBdx(1, 6) =  db7dxi1 * jacInv[0][1]
              + db7dxi2 * jacInv[1][1]
              + db7dxi3 * jacInv[2][1];

  dBdx(2, 6) =  db7dxi1 * jacInv[0][2]
              + db7dxi2 * jacInv[1][2]
              + db7dxi3 * jacInv[2][2];

  dBdx(0, 7) =  db8dxi1 * jacInv[0][0]
              + db8dxi2 * jacInv[1][0]
              + db8dxi3 * jacInv[2][0];

  dBdx(1, 7) =  db8dxi1 * jacInv[0][1]
              + db8dxi2 * jacInv[1][1]
              + db8dxi3 * jacInv[2][1];

  dBdx(2, 7) =  db8dxi1 * jacInv[0][2]
              + db8dxi2 * jacInv[1][2]
              + db8dxi3 * jacInv[2][2];

  dBdx(0, 8) =  db9dxi1 * jacInv[0][0]
              + db9dxi2 * jacInv[1][0]
              + db9dxi3 * jacInv[2][0];

  dBdx(1, 8) =  db9dxi1 * jacInv[0][1]
              + db9dxi2 * jacInv[1][1]
              + db9dxi3 * jacInv[2][1];

  dBdx(2, 8) =  db9dxi1 * jacInv[0][2]
              + db9dxi2 * jacInv[1][2]
              + db9dxi3 * jacInv[2][2];

  dBdx(0, 9) =  db10dxi1 * jacInv[0][0]
              + db10dxi2 * jacInv[1][0]
              + db10dxi3 * jacInv[2][0];

  dBdx(1, 9) =  db10dxi1 * jacInv[0][1]
              + db10dxi2 * jacInv[1][1]
              + db10dxi3 * jacInv[2][1];

  dBdx(2, 9) =  db10dxi1 * jacInv[0][2]
              + db10dxi2 * jacInv[1][2]
              + db10dxi3 * jacInv[2][2];
}

//! Compute the Dubiner basis functions
std::vector< tk::real >
eval_basis( const std::size_t ndof,
            const tk::real xi,
            const tk::real eta,
            const tk::real zeta );

//! overloaded function for eval_basis for Kokkos 
KOKKOS_INLINE_FUNCTION 
void eval_basis( const std::size_t ndof,
                const tk::real xi,
                const tk::real eta,
                const tk::real zeta, 
                Kokkos::View<real*, memory_space> B)
// *****************************************************************************
//  Compute the Dubiner basis functions
//! \param[in] ndof Number of degrees of freedom
//! \param[in] xi,eta,zeta Coordinates for quadrature points in reference space
//! \return Vector of basis functions
// *****************************************************************************
{
  // Array of basis functions

  B(0) = 1.0;

  if ( ndof > 1 )           // DG(P1)
  {
    B(1) = 2.0 * xi + eta + zeta - 1.0;
    B(2) = 3.0 * eta + zeta - 1.0;
    B(3) = 4.0 * zeta - 1.0;

    if( ndof > 4 )         // DG(P2)
    {
      B(4) =  6.0 * xi * xi + eta * eta + zeta * zeta
            + 6.0 * xi * eta + 6.0 * xi * zeta + 2.0 * eta * zeta
            - 6.0 * xi - 2.0 * eta - 2.0 * zeta + 1.0;
      B(5) =  5.0 * eta * eta + zeta * zeta
            + 10.0 * xi * eta + 2.0 * xi * zeta + 6.0 * eta * zeta
            - 2.0 * xi - 6.0 * eta - 2.0 * zeta + 1.0;
      B(6) =  6.0 * zeta * zeta + 12.0 * xi * zeta + 6.0 * eta * zeta - 2.0 * xi
            - eta - 7.0 * zeta + 1.0;
      B(7) =  10.0 * eta * eta + zeta * zeta + 8.0 * eta * zeta
            - 8.0 * eta - 2.0 * zeta + 1.0;
      B(8) =  6.0 * zeta * zeta + 18.0 * eta * zeta - 3.0 * eta - 7.0 * zeta
            + 1.0;
      B(9) =  15.0 * zeta * zeta - 10.0 * zeta + 1.0;
    }
  }
}

//! Compute the state variables for the tetrahedron element
std::vector< tk::real >
eval_state ( ncomp_t ncomp,
             const std::size_t ndof,
             const std::size_t ndof_el,
             const std::size_t e,
             const Fields& U,
             const std::vector< tk::real >& B );
  
KOKKOS_INLINE_FUNCTION 
void eval_state ( ncomp_t ncomp,
                 const std::size_t ndof,
                 const std::size_t ndof_el,
                 const std::size_t e, 
                 size_t m_nprop,
                 Kokkos::View<const tk::real*, memory_space> U,
                 Kokkos::View<const tk::real*, memory_space> B, 
                 Kokkos::View<tk::real*, memory_space> state,
                 const size_t& idx)
// *****************************************************************************
//  Compute the state variables for the tetrahedron element
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_el Number of degrees of freedom for the local element
//! \param[in] e Index for the tetrahedron element
//! \param[in] U Solution vector at recent time step
//! \param[in] B Vector of basis functions
//! \return Vector of state variable for tetrahedron element
// *****************************************************************************
{
  // This is commented for now because that when p0/p1 adaptive with limiter
  // applied, the size of basis will be 10. However, ndof_el will be 4 which
  // leads to a size mismatch in limiter function.
  //Assert( B.size() == ndof_el, "Size mismatch" );

  // Array of state variable for tetrahedron element
  /*
   const tk::real&
    access( ncomp_t unknown, ncomp_t component, int2type< UnkEqComp > ) const
    {
      Assert( component < m_nprop, "Out-of-bounds access: "
              "component < number of properties" );
      Assert( unknown < m_nunk, "Out-of-bounds access: unknown < number of "
              "unknowns" );
      return m_vec[ unknown*m_nprop + component ];
      }
    */
    
  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    state(c + idx) = U(e * m_nprop + mark);
    // if (idx != 0)
    // {
    //   if (c == inciter::velocityIdx(2,0))
    //     printf("U = %e\n", U(e * m_nprop + mark));
    //   if (c == inciter::velocityIdx(2,1))
    //     printf("V = %e\n", U(e * m_nprop + mark));
    //   if (c == inciter::velocityIdx(2,2))
    //     printf("W = %e\n", U(e * m_nprop + mark));
    // }

    if(ndof_el > 1)        // Second order polynomial solution
    {
      state(c + idx) += U( e * m_nprop +mark+1 ) * B(1)
                + U( e * m_nprop +mark+2 ) * B(2)
                + U( e * m_nprop +mark+3 ) * B(3);
    }

    if(ndof_el > 4)        // Third order polynomial solution
    {
      state(c + idx) += U( e * m_nprop + mark+4 ) * B(4)
                + U( e * m_nprop + mark+5 ) * B(5)
                + U( e * m_nprop + mark+6 ) * B(6)
                + U( e * m_nprop + mark+7 ) * B(7)
                + U( e * m_nprop + mark+8 ) * B(8)
                + U( e * m_nprop + mark+9 ) * B(9);
    }
  }
}

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
