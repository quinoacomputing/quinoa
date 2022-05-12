// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Basis.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing the Dubiner basis functions in DG methods
  \details   This file contains functionality for computing the basis functions
     and relating coordinates transformation functions used in discontinuous
     Galerkin methods for variaous orders of numerical representation. The basis
     functions chosen for the DG method are the Dubiner basis, which are
     Legendre polynomials modified for tetrahedra, which are defined only on
     the reference/master tetrahedron.
  \see [1] https://doi.org/10.1007/BF01060030
  \see [2] https://doi.org/10.1093/imamat/hxh111
*/
// *****************************************************************************

#include <array>

#include "Basis.hpp"

std::array< tk::real, 3 >
tk::eval_gp ( const std::size_t igp,
              const std::array< std::array< tk::real, 3>, 3 >& coordfa,
              const std::array< std::vector< tk::real >, 2 >& coordgp )
// *****************************************************************************
//  Compute the coordinates of quadrature points for face integral in physical
//  space
//! \param[in] igp Index of quadrature points
//! \param[in] coordfa Array of nodal coordinates for face element
//! \param[in] coordgp Array of coordinates for quadrature points in reference
//!   space
//! \return Array of coordinates for quadrature points in physical space
// *****************************************************************************
{
  // Barycentric coordinates for the triangular face
  auto shp1 = 1.0 - coordgp[0][igp] - coordgp[1][igp];
  auto shp2 = coordgp[0][igp];
  auto shp3 = coordgp[1][igp];

  // Transformation of the quadrature point from the 2D reference/master
  // element to physical space, to obtain its physical (x,y,z) coordinates.
  return {{ coordfa[0][0]*shp1 + coordfa[1][0]*shp2 + coordfa[2][0]*shp3,
            coordfa[0][1]*shp1 + coordfa[1][1]*shp2 + coordfa[2][1]*shp3,
            coordfa[0][2]*shp1 + coordfa[1][2]*shp2 + coordfa[2][2]*shp3 }};
}

std::array< tk::real, 3 >
tk::eval_gp ( const std::size_t igp,
              const std::array< std::array< tk::real, 3>, 4 >& coord,
              const std::array< std::vector< tk::real >, 3 >& coordgp )
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
  auto shp1 = 1.0 - coordgp[0][igp] - coordgp[1][igp] - coordgp[2][igp];
  auto shp2 = coordgp[0][igp];
  auto shp3 = coordgp[1][igp];
  auto shp4 = coordgp[2][igp];

  // Transformation of the quadrature point from the reference/master
  // element to physical space, to obtain its physical (x,y,z) coordinates.
  return {{
   coord[0][0]*shp1 + coord[1][0]*shp2 + coord[2][0]*shp3 + coord[3][0]*shp4,
   coord[0][1]*shp1 + coord[1][1]*shp2 + coord[2][1]*shp3 + coord[3][1]*shp4,
   coord[0][2]*shp1 + coord[1][2]*shp2 + coord[2][2]*shp3 + coord[3][2]*shp4 }};
}

std::array< std::vector<tk::real>, 3 >
tk::eval_dBdx_p1( const std::size_t ndof,
                  const std::array< std::array< tk::real, 3 >, 3 >& jacInv )
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

  std::array< std::vector<tk::real>, 3 > dBdx;
  dBdx[0].resize( ndof, 0 );
  dBdx[1].resize( ndof, 0 );
  dBdx[2].resize( ndof, 0 );

  auto db2dxi1 = 2.0;
  auto db2dxi2 = 1.0;
  auto db2dxi3 = 1.0;

  auto db3dxi1 = 0.0;
  auto db3dxi2 = 3.0;
  auto db3dxi3 = 1.0;

  auto db4dxi1 = 0.0;
  auto db4dxi2 = 0.0;
  auto db4dxi3 = 4.0;

  dBdx[0][1] =  db2dxi1 * jacInv[0][0]
              + db2dxi2 * jacInv[1][0]
              + db2dxi3 * jacInv[2][0];

  dBdx[1][1] =  db2dxi1 * jacInv[0][1]
              + db2dxi2 * jacInv[1][1]
              + db2dxi3 * jacInv[2][1];

  dBdx[2][1] =  db2dxi1 * jacInv[0][2]
              + db2dxi2 * jacInv[1][2]
              + db2dxi3 * jacInv[2][2];

  dBdx[0][2] =  db3dxi1 * jacInv[0][0]
              + db3dxi2 * jacInv[1][0]
              + db3dxi3 * jacInv[2][0];

  dBdx[1][2] =  db3dxi1 * jacInv[0][1]
              + db3dxi2 * jacInv[1][1]
              + db3dxi3 * jacInv[2][1];

  dBdx[2][2] =  db3dxi1 * jacInv[0][2]
              + db3dxi2 * jacInv[1][2]
              + db3dxi3 * jacInv[2][2];

  dBdx[0][3] =  db4dxi1 * jacInv[0][0]
              + db4dxi2 * jacInv[1][0]
              + db4dxi3 * jacInv[2][0];

  dBdx[1][3] =  db4dxi1 * jacInv[0][1]
              + db4dxi2 * jacInv[1][1]
              + db4dxi3 * jacInv[2][1];

  dBdx[2][3] =  db4dxi1 * jacInv[0][2]
              + db4dxi2 * jacInv[1][2]
              + db4dxi3 * jacInv[2][2];

  return dBdx;
}

void
tk::eval_dBdx_p2( const std::size_t igp,
                  const std::array< std::vector< tk::real >, 3 >& coordgp,
                  const std::array< std::array< tk::real, 3 >, 3 >& jacInv,
                  std::array< std::vector<tk::real>, 3 >& dBdx )
// *****************************************************************************
//  Compute the derivatives of Dubiner basis function for DG(P2)
//! \param[in] igp Index of quadrature points
//! \param[in] coordgp Gauss point coordinates for tetrahedron element
//! \param[in] jacInv Array of the inverse of Jacobian
//! \param[in,out] dBdx Array of the derivatives of basis function
// *****************************************************************************
{
  auto db5dxi1 = 12.0 * coordgp[0][igp] + 6.0 * coordgp[1][igp]
               +  6.0 * coordgp[2][igp] - 6.0;
  auto db5dxi2 =  6.0 * coordgp[0][igp] + 2.0 * coordgp[1][igp]
               +  2.0 * coordgp[2][igp] - 2.0;
  auto db5dxi3 =  6.0 * coordgp[0][igp] + 2.0 * coordgp[1][igp]
               +  2.0 * coordgp[2][igp] - 2.0;

  auto db6dxi1 = 10.0 * coordgp[1][igp] +  2.0 * coordgp[2][igp] - 2.0;
  auto db6dxi2 = 10.0 * coordgp[0][igp] + 10.0 * coordgp[1][igp]
               +  6.0 * coordgp[2][igp] - 6.0;
  auto db6dxi3 =  2.0 * coordgp[0][igp] +  6.0 * coordgp[1][igp]
               +  2.0 * coordgp[2][igp] - 2.0;

  auto db7dxi1 = 12.0 * coordgp[2][igp] - 2.0;
  auto db7dxi2 =  6.0 * coordgp[2][igp] - 1.0;
  auto db7dxi3 = 12.0 * coordgp[0][igp] + 6.0 * coordgp[1][igp]
               + 12.0 * coordgp[2][igp] - 7.0;

  auto db8dxi1 =  0;
  auto db8dxi2 = 20.0 * coordgp[1][igp] + 8.0 * coordgp[2][igp] - 8.0;
  auto db8dxi3 =  8.0 * coordgp[1][igp] + 2.0 * coordgp[2][igp] - 2.0;

  auto db9dxi1 =  0;
  auto db9dxi2 = 18.0 * coordgp[2][igp] -  3.0;
  auto db9dxi3 = 18.0 * coordgp[1][igp] + 12.0 * coordgp[2][igp] - 7.0;

  auto db10dxi1 =  0;
  auto db10dxi2 =  0;
  auto db10dxi3 = 30.0 * coordgp[2][igp] - 10.0;

  dBdx[0][4] =  db5dxi1 * jacInv[0][0]
              + db5dxi2 * jacInv[1][0]
              + db5dxi3 * jacInv[2][0];

  dBdx[1][4] =  db5dxi1 * jacInv[0][1]
              + db5dxi2 * jacInv[1][1]
              + db5dxi3 * jacInv[2][1];

  dBdx[2][4] =  db5dxi1 * jacInv[0][2]
              + db5dxi2 * jacInv[1][2]
              + db5dxi3 * jacInv[2][2];

  dBdx[0][5] =  db6dxi1 * jacInv[0][0]
              + db6dxi2 * jacInv[1][0]
              + db6dxi3 * jacInv[2][0];

  dBdx[1][5] =  db6dxi1 * jacInv[0][1]
              + db6dxi2 * jacInv[1][1]
              + db6dxi3 * jacInv[2][1];

  dBdx[2][5] =  db6dxi1 * jacInv[0][2]
              + db6dxi2 * jacInv[1][2]
              + db6dxi3 * jacInv[2][2];

  dBdx[0][6] =  db7dxi1 * jacInv[0][0]
              + db7dxi2 * jacInv[1][0]
              + db7dxi3 * jacInv[2][0];

  dBdx[1][6] =  db7dxi1 * jacInv[0][1]
              + db7dxi2 * jacInv[1][1]
              + db7dxi3 * jacInv[2][1];

  dBdx[2][6] =  db7dxi1 * jacInv[0][2]
              + db7dxi2 * jacInv[1][2]
              + db7dxi3 * jacInv[2][2];

  dBdx[0][7] =  db8dxi1 * jacInv[0][0]
              + db8dxi2 * jacInv[1][0]
              + db8dxi3 * jacInv[2][0];

  dBdx[1][7] =  db8dxi1 * jacInv[0][1]
              + db8dxi2 * jacInv[1][1]
              + db8dxi3 * jacInv[2][1];

  dBdx[2][7] =  db8dxi1 * jacInv[0][2]
              + db8dxi2 * jacInv[1][2]
              + db8dxi3 * jacInv[2][2];

  dBdx[0][8] =  db9dxi1 * jacInv[0][0]
              + db9dxi2 * jacInv[1][0]
              + db9dxi3 * jacInv[2][0];

  dBdx[1][8] =  db9dxi1 * jacInv[0][1]
              + db9dxi2 * jacInv[1][1]
              + db9dxi3 * jacInv[2][1];

  dBdx[2][8] =  db9dxi1 * jacInv[0][2]
              + db9dxi2 * jacInv[1][2]
              + db9dxi3 * jacInv[2][2];

  dBdx[0][9] =  db10dxi1 * jacInv[0][0]
              + db10dxi2 * jacInv[1][0]
              + db10dxi3 * jacInv[2][0];

  dBdx[1][9] =  db10dxi1 * jacInv[0][1]
              + db10dxi2 * jacInv[1][1]
              + db10dxi3 * jacInv[2][1];

  dBdx[2][9] =  db10dxi1 * jacInv[0][2]
              + db10dxi2 * jacInv[1][2]
              + db10dxi3 * jacInv[2][2];
}

std::vector< tk::real >
tk::eval_basis( const std::size_t ndof,
                const tk::real xi,
                const tk::real eta,
                const tk::real zeta )
// *****************************************************************************
//  Compute the Dubiner basis functions
//! \param[in] ndof Number of degrees of freedom
//! \param[in] xi,eta,zeta Coordinates for quadrature points in reference space
//! \return Vector of basis functions
// *****************************************************************************
{
  // Array of basis functions
  std::vector< tk::real > B( ndof, 1.0 );

  if ( ndof > 1 )           // DG(P1)
  {
    B[1] = 2.0 * xi + eta + zeta - 1.0;
    B[2] = 3.0 * eta + zeta - 1.0;
    B[3] = 4.0 * zeta - 1.0;

    if( ndof > 4 )         // DG(P2)
    {
      B[4] =  6.0 * xi * xi + eta * eta + zeta * zeta
            + 6.0 * xi * eta + 6.0 * xi * zeta + 2.0 * eta * zeta
            - 6.0 * xi - 2.0 * eta - 2.0 * zeta + 1.0;
      B[5] =  5.0 * eta * eta + zeta * zeta
            + 10.0 * xi * eta + 2.0 * xi * zeta + 6.0 * eta * zeta
            - 2.0 * xi - 6.0 * eta - 2.0 * zeta + 1.0;
      B[6] =  6.0 * zeta * zeta + 12.0 * xi * zeta + 6.0 * eta * zeta - 2.0 * xi
            - eta - 7.0 * zeta + 1.0;
      B[7] =  10.0 * eta * eta + zeta * zeta + 8.0 * eta * zeta
            - 8.0 * eta - 2.0 * zeta + 1.0;
      B[8] =  6.0 * zeta * zeta + 18.0 * eta * zeta - 3.0 * eta - 7.0 * zeta
            + 1.0;
      B[9] =  15.0 * zeta * zeta - 10.0 * zeta + 1.0;
    }
  }

  return B;
}

std::vector< tk::real >
tk::eval_state ( ncomp_t ncomp,
                 ncomp_t offset,
                 const std::size_t ndof,
                 const std::size_t ndof_el,
                 const std::size_t e,
                 const Fields& U,
                 const std::vector< tk::real >& B,
                 const std::array< std::size_t, 2 >& VarRange )
// *****************************************************************************
//  Compute the state variables for the tetrahedron element
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_el Number of degrees of freedom for the local element
//! \param[in] e Index for the tetrahedron element
//! \param[in] U Solution vector at recent time step
//! \param[in] B Vector of basis functions
//! \param[in] VarRange Range of the variables to be evaluated
//! \return Vector of state variable for tetrahedron element
// *****************************************************************************
{
  // This is commented for now because that when p0/p1 adaptive with limiter
  // applied, the size of basis will be 10. However, ndof_el will be 4 which
  // leads to a size mismatch in limiter function.
  //Assert( B.size() == ndof_el, "Size mismatch" );

  if (U.empty()) return {};

  // Array of state variable for tetrahedron element
  std::vector< tk::real > state( ncomp, 0.0 );

  for (ncomp_t c=VarRange[0]; c<=VarRange[1]; ++c)
  {
    auto mark = c*ndof;
    state[c] = U( e, mark, offset );

    if(ndof_el > 1)        //DG(P1)
    {
      state[c] += U( e, mark+1, offset ) * B[1]
                + U( e, mark+2, offset ) * B[2]
                + U( e, mark+3, offset ) * B[3];
    }

    if(ndof_el > 4)        //DG(P2)
    {
      state[c] += U( e, mark+4, offset ) * B[4]
                + U( e, mark+5, offset ) * B[5]
                + U( e, mark+6, offset ) * B[6]
                + U( e, mark+7, offset ) * B[7]
                + U( e, mark+8, offset ) * B[8]
                + U( e, mark+9, offset ) * B[9];
    }
  }

  return state;
}

std::vector< std::vector< tk::real > >
tk::DubinerToTaylor( ncomp_t ncomp,
                     ncomp_t offset,
                     const std::size_t e,
                     const std::size_t ndof,
                     const tk::Fields& U,
                     const std::vector< std::size_t >& inpoel,
                     const tk::UnsMesh::Coords& coord )
// *****************************************************************************
//  Transform the solution with Dubiner basis to the solution with Taylor basis
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Index for equation systems
//! \param[in] e Id of element whose solution is to be limited
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] U High-order solution vector with Dubiner basis
//! \param[in] inpoel Element connectivity
//! \param[in] coord Array of nodal coordinates
//! \return High-order solution vector with Taylor basis
// *****************************************************************************
{
  std::vector< std::vector< tk::real > >
    unk(ncomp, std::vector<tk::real>(ndof, 0.0));

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  std::array< std::vector< tk::real >, 3 > center;
  center[0].resize(1, 0.25);
  center[1].resize(1, 0.25);
  center[2].resize(1, 0.25);

  // Evaluate the cell center solution
  for(ncomp_t icomp = 0; icomp < ncomp; icomp++)
  {
    auto mark = icomp * ndof;
    unk[icomp][0] = U(e, mark, offset);
  }

  // Evaluate the first order derivative
  std::array< std::array< tk::real, 3>, 4 > coordel {{
    {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
    {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
    {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
    {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }}
  }};

  auto jacInv =
              tk::inverseJacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

  // Compute the derivatives of basis function for DG(P1)
  auto dBdx = tk::eval_dBdx_p1( ndof, jacInv );

  if(ndof > 4) {
    tk::eval_dBdx_p2(0, center, jacInv, dBdx);
  }

  for(ncomp_t icomp = 0; icomp < ncomp; icomp++)
  {
    auto mark = icomp * ndof; 
    for(std::size_t idir = 0; idir < 3; idir++)
    {
      unk[icomp][idir+1] = 0;
      for(std::size_t idof = 1; idof < ndof; idof++)
        unk[icomp][idir+1] += U(e, mark+idof, offset) * dBdx[idir][idof];
    }
  }

  // Evaluate the second order derivative if DGP2 is applied
  // The basic idea of the computation follows
  //    d2Udx2 = /sum u_i * (d2B_i/dx2)
  // where d2B_i/dx2 = d( dB_i/dxi * dxi/dx ) / dxi * dxi/dx
  if(ndof > 4)
  {
    // Matrix to store the second order derivatives of basis functions in
    // reference domain
    tk::real d2Bdxi2[6][6] =
    { { 12.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
      {  2.0, 10.0,  0.0, 20.0,  0.0,  0.0 },
      {  2.0,  2.0, 12.0,  2.0, 12.0, 30.0 },
      {  6.0, 10.0,  0.0,  0.0,  0.0,  0.0 },
      {  6.0,  2.0, 12.0,  0.0,  0.0,  0.0 },
      {  2.0,  6.0,  6.0,  8.0, 18.0,  0.0 } };

    // Transform matrix to convert the second order derivatives of basis
    // function in reference domain to the one in physical domain
    tk::real d2xdxi2[6][6];

    d2xdxi2[0][0] = jacInv[0][0] * jacInv[0][0];
    d2xdxi2[0][1] = jacInv[1][0] * jacInv[1][0];
    d2xdxi2[0][2] = jacInv[2][0] * jacInv[2][0];
    d2xdxi2[0][3] = jacInv[0][0] * jacInv[1][0] * 2.0;
    d2xdxi2[0][4] = jacInv[0][0] * jacInv[2][0] * 2.0;
    d2xdxi2[0][5] = jacInv[1][0] * jacInv[2][0] * 2.0;

    d2xdxi2[1][0] = jacInv[0][1] * jacInv[0][1];
    d2xdxi2[1][1] = jacInv[1][1] * jacInv[1][1];
    d2xdxi2[1][2] = jacInv[2][1] * jacInv[2][1];
    d2xdxi2[1][3] = jacInv[0][1] * jacInv[1][1] * 2.0;
    d2xdxi2[1][4] = jacInv[0][1] * jacInv[2][1] * 2.0;
    d2xdxi2[1][5] = jacInv[1][1] * jacInv[2][1] * 2.0;

    d2xdxi2[2][0] = jacInv[0][2] * jacInv[0][2];
    d2xdxi2[2][1] = jacInv[1][2] * jacInv[1][2];
    d2xdxi2[2][2] = jacInv[2][2] * jacInv[2][2];
    d2xdxi2[2][3] = jacInv[0][2] * jacInv[1][2] * 2.0;
    d2xdxi2[2][4] = jacInv[0][2] * jacInv[2][2] * 2.0;
    d2xdxi2[2][5] = jacInv[1][2] * jacInv[2][2] * 2.0;

    d2xdxi2[3][0] = jacInv[0][0] * jacInv[0][1];
    d2xdxi2[3][1] = jacInv[1][0] * jacInv[1][1];
    d2xdxi2[3][2] = jacInv[2][0] * jacInv[2][1];
    d2xdxi2[3][3] = jacInv[0][0] * jacInv[1][1] + jacInv[1][0] * jacInv[0][1];
    d2xdxi2[3][4] = jacInv[0][0] * jacInv[2][1] + jacInv[2][0] * jacInv[0][1];
    d2xdxi2[3][5] = jacInv[1][0] * jacInv[2][1] + jacInv[2][0] * jacInv[1][1];

    d2xdxi2[4][0] = jacInv[0][0] * jacInv[0][2];
    d2xdxi2[4][1] = jacInv[1][0] * jacInv[1][2];
    d2xdxi2[4][2] = jacInv[2][0] * jacInv[2][2];
    d2xdxi2[4][3] = jacInv[0][0] * jacInv[1][2] + jacInv[1][0] * jacInv[0][2];
    d2xdxi2[4][4] = jacInv[0][0] * jacInv[2][2] + jacInv[2][0] * jacInv[0][2];
    d2xdxi2[4][5] = jacInv[1][0] * jacInv[2][2] + jacInv[2][0] * jacInv[1][2];

    d2xdxi2[5][0] = jacInv[0][1] * jacInv[0][2];
    d2xdxi2[5][1] = jacInv[1][1] * jacInv[1][2];
    d2xdxi2[5][2] = jacInv[2][1] * jacInv[2][2];
    d2xdxi2[5][3] = jacInv[0][1] * jacInv[1][2] + jacInv[1][1] * jacInv[0][2];
    d2xdxi2[5][4] = jacInv[0][1] * jacInv[2][2] + jacInv[2][1] * jacInv[0][2];
    d2xdxi2[5][5] = jacInv[1][1] * jacInv[2][2] + jacInv[2][1] * jacInv[1][2];

    // Matrix to store the second order derivatives of basis functions in
    // physical domain
    tk::real d2Bdx2[6][6];
    for(std::size_t ibasis = 0; ibasis < 6; ibasis++) {
      for(std::size_t idir = 0; idir < 6; idir++) {
        d2Bdx2[idir][ibasis] = 0;
        for(std::size_t k = 0; k < 6; k++)
          d2Bdx2[idir][ibasis] += d2xdxi2[idir][k] * d2Bdxi2[k][ibasis];
      }
    }

    for(ncomp_t icomp = 0; icomp < ncomp; icomp++)
    {
      auto mark = icomp * ndof;
      for(std::size_t idir = 0; idir < 6; idir++)
      {
        unk[icomp][idir+4] = 0;
        for(std::size_t ibasis = 0; ibasis < 6; ibasis++)
          unk[icomp][idir+4] += U(e, mark+4+ibasis, offset) * d2Bdx2[idir][ibasis];
      }
    }
  }
  return unk;
}

void
tk::TaylorToDubiner( ncomp_t ncomp,
                     std::size_t e,
                     std::size_t ndof,
                     const std::vector< std::size_t >& inpoel,
                     const tk::UnsMesh::Coords& coord,
                     const tk::Fields& geoElem,
                     std::vector< std::vector< tk::real > >& unk )
// *****************************************************************************
//  Convert the solution with Taylor basis to the solution with Dubiner basis by
//    projection method
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] e Id of element whose solution is to be limited
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] inpoel Element connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in, out] unk High-order solution vector with Taylor basis
// *****************************************************************************
{
  Assert( ncomp > 0, "Number of scalar components is incorrect" );
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  std::array< std::array< tk::real, 3>, 4 > coordel {{
    {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
    {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
    {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
    {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }} }};

  std::array< tk::real, 3 > dx;
  std::array< tk::real, 3 > dy;
  std::array< tk::real, 3 > dz;

  dx[0] = coordel[1][0] - coordel[0][0];
  dx[1] = coordel[2][0] - coordel[0][0];
  dx[2] = coordel[3][0] - coordel[0][0];

  dy[0] = coordel[1][1] - coordel[0][1];
  dy[1] = coordel[2][1] - coordel[0][1];
  dy[2] = coordel[3][1] - coordel[0][1];

  dz[0] = coordel[1][2] - coordel[0][2];
  dz[1] = coordel[2][2] - coordel[0][2];
  dz[2] = coordel[3][2] - coordel[0][2];

  std::array< tk::real, 3 > x0{coordel[0][0], coordel[0][1], coordel[0][2]};

  for (ncomp_t c=0; c<ncomp; ++c) {
    // Step 1: Rewritten the solution in (x,y,z) coordinate
    // Vector stores the coefficients for x^2 etc and the constant term is not
    // needed here.
    std::vector< tk::real > A(ndof-1, 0);
    A[0] = unk[c][4] / 2.0;     // Coefficient for x^2
    A[1] = unk[c][5] / 2.0;     // Coefficient for y^2
    A[2] = unk[c][6] / 2.0;     // Coefficient for z^2
    A[3] = unk[c][7];           // Coefficient for xy
    A[4] = unk[c][8];           // Coefficient for xz
    A[5] = unk[c][9];           // Coefficient for yz
    // Coefficient for x
    A[6] = unk[c][1] - unk[c][4]*geoElem(e, 1, 0) - unk[c][7]*geoElem(e, 2, 0)
         - unk[c][8]*geoElem(e, 3, 0);
    // Coefficient for y
    A[7] = unk[c][2] - unk[c][5]*geoElem(e, 2, 0) - unk[c][7]*geoElem(e, 1, 0)
         - unk[c][9]*geoElem(e, 3, 0);
    // Coefficient for z
    A[8] = unk[c][3] - unk[c][6]*geoElem(e, 3, 0) - unk[c][8]*geoElem(e, 1, 0)
         - unk[c][9]*geoElem(e, 2, 0);

    // Step 2: Transform the solution in (xi, eta, zeta) coordinate
    // Vector stores the coefficients for xi^2 etc and the constant term is not
    // needed here.
    std::vector< tk::real > T(ndof-1, 0);
    // Coefficient for xi^2
    T[0] = A[0]*dx[0]*dx[0] + A[1]*dy[0]*dy[0] + A[2]*dz[0]*dz[0]
         + A[3]*dx[0]*dy[0] + A[4]*dx[0]*dz[0] + A[5]*dy[0]*dz[0];
    // Coefficient for eta^2
    T[1] = A[0]*dx[1]*dx[1] + A[1]*dy[1]*dy[1] + A[2]*dz[1]*dz[1]
         + A[3]*dx[1]*dy[1] + A[4]*dx[1]*dz[1] + A[5]*dy[1]*dz[1];
    // Coefficient for zeta^2
    T[2] = A[0]*dx[2]*dx[2] + A[1]*dy[2]*dy[2] + A[2]*dz[2]*dz[2]
         + A[3]*dx[2]*dy[2] + A[4]*dx[2]*dz[2] + A[5]*dy[2]*dz[2];
    // Coefficient for xi*eta
    T[3] = 2.0 * (A[0]*dx[0]*dx[1] + A[1]*dy[0]*dy[1] + A[2]*dz[0]*dz[1])
         + A[3]*(dx[0]*dy[1] + dx[1]*dy[0]) + A[4]*(dx[0]*dz[1] + dx[1]*dz[0])
         + A[5]*(dy[0]*dz[1] + dy[1]*dz[0]);
    // Coefficient for xi*zeta
    T[4] = 2.0 * (A[0]*dx[0]*dx[2] + A[1]*dy[0]*dy[2] + A[2]*dz[0]*dz[2])
         + A[3]*(dx[0]*dy[2] + dx[2]*dy[0]) + A[4]*(dx[0]*dz[2] + dx[2]*dz[0])
         + A[5]*(dy[0]*dz[2] + dy[2]*dz[0]);
    // Coefficient for eta*zeta
    T[5] = 2.0 * (A[0]*dx[1]*dx[2] + A[1]*dy[1]*dy[2] + A[2]*dz[1]*dz[2])
         + A[3]*(dx[1]*dy[2] + dx[2]*dy[1]) + A[4]*(dx[1]*dz[2] + dx[2]*dz[1])
         + A[5]*(dy[1]*dz[2] + dy[2]*dz[1]);
    // Coefficient for xi
    T[6] = 2.0 * (A[0]*dx[0]*x0[0] + A[1]*dy[0]*x0[1] + A[2]*dz[0]*x0[2])
         + A[3]*(dx[0]*x0[1] + dy[0]*x0[0]) + A[4]*(dx[0]*x0[2] + dz[0]*x0[0])
         + A[5]*(dy[0]*x0[2] + dz[0]*x0[1]) + A[6]*dx[0] + A[7]*dy[0] + A[8]*dz[0];
    // Coefficient for eta
    T[7] = 2.0 * (A[0]*dx[1]*x0[0] + A[1]*dy[1]*x0[1] + A[2]*dz[1]*x0[2])
         + A[3]*(dx[1]*x0[1] + dy[1]*x0[0]) + A[4]*(dx[1]*x0[2] + dz[1]*x0[0])
         + A[5]*(dy[1]*x0[2] + dz[1]*x0[1]) + A[6]*dx[1] + A[7]*dy[1] + A[8]*dz[1];
    // Coefficient for zeta
    T[8] = 2.0 * (A[0]*dx[2]*x0[0] + A[1]*dy[2]*x0[1] + A[2]*dz[2]*x0[2])
         + A[3]*(dx[2]*x0[1] + dy[2]*x0[0]) + A[4]*(dx[2]*x0[2] + dz[2]*x0[0])
         + A[5]*(dy[2]*x0[2] + dz[2]*x0[1]) + A[6]*dx[2] + A[7]*dy[2] + A[8]*dz[2];

    // Step 3: Solve the linear system when coefficents are identical
    unk[c][4] = T[0] / 6.0;
    unk[c][5] = (T[3] - 6.0*unk[c][4]) / 10.0;
    unk[c][6] = (T[4] - 6.0*unk[c][4] - 2.0*unk[c][5]) / 12.0;
    unk[c][7] = (T[1] - unk[c][4] - 5.0*unk[c][5]) / 10.0;
    unk[c][8] = (T[5] - 2.0*unk[c][4] - 6.0*unk[c][5] - 6.0*unk[c][6]
              - 8.0*unk[c][7]) / 18.0;
    unk[c][9] = (T[2] - unk[c][4] - unk[c][5] - 6.0*unk[c][6] - unk[c][7]
              - 6.0*unk[c][8]) / 15.0;
    unk[c][1] = (T[6] + 6.0*unk[c][4] + 2.0*unk[c][5] + 2.0*unk[c][6]) / 2.0;
    unk[c][2] = (T[7] - unk[c][1] + 2.0*unk[c][4] + 6.0*unk[c][5] + unk[c][6]
              + 8.0*unk[c][7] + 3.0*unk[c][8]) / 3.0;
    unk[c][3] = (T[8] - unk[c][1] - unk[c][2] + 2.0*unk[c][4] + 2.0*unk[c][5]
              + 7.0*unk[c][6] + 2.0*unk[c][7] + 7.0*unk[c][8] + 10.0*unk[c][9])
              / 4.0;
  }
}
//{
//  Assert( ncomp > 0, "Number of scalar components is incorrect" );
//
//  // The diagonal of mass matrix
//  std::vector< tk::real > L(ndof, 0.0);
//
//  tk::real vol = 1.0 / 6.0;
//
//  L[0] = vol;
//
//  if(ndof > 1) {
//    Assert( (ndof == 4)||(ndof == 10),
//      "Mismatch in number of degrees of freedom" );
//    L[1] = vol / 10.0;
//    L[2] = vol * 3.0/10.0;
//    L[3] = vol * 3.0/5.0;
//  }
//
//  if(ndof > 4) {
//    Assert( ndof == 10, "Mismatch in number of degrees of freedom" );
//    L[4] = vol / 35.0;
//    L[5] = vol / 21.0;
//    L[6] = vol / 14.0;
//    L[7] = vol / 7.0;
//    L[8] = vol * 3.0/14.0;
//    L[9] = vol * 3.0/7.0;
//  }
//
//  // Coordinates of the centroid in physical domain
//  std::array< tk::real, 3 > x_c{geoElem(e,1,0), geoElem(e,2,0), geoElem(e,3,0)};
//
//  const auto& cx = coord[0];
//  const auto& cy = coord[1];
//  const auto& cz = coord[2];
//
//  std::array< std::array< tk::real, 3>, 4 > coordel {{
//    {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
//    {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
//    {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
//    {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }}
//  }};
//
//  // Number of quadrature points for volume integration
//  auto ng = tk::NGvol(ndof);
//
//  // arrays for quadrature points
//  std::array< std::vector< tk::real >, 3 > coordgp;
//  std::vector< tk::real > wgp;
//
//  coordgp[0].resize( ng );
//  coordgp[1].resize( ng );
//  coordgp[2].resize( ng );
//  wgp.resize( ng );
//
//  // get quadrature point weights and coordinates for triangle
//  tk::GaussQuadratureTet( ng, coordgp, wgp );
//
//  // right hand side vector
//  std::vector< tk::real > R( ncomp*ndof, 0.0 );
//
//  // Gaussian quadrature
//  for (std::size_t igp=0; igp<ng; ++igp)
//  {
//    auto wt = wgp[igp] * vol;
//
//    auto gp = tk::eval_gp( igp, coordel, coordgp );
//
//    auto B_taylor = eval_TaylorBasis( ndof, gp, x_c, coordel);
//
//    // Compute high order solution at gauss point
//    std::vector< tk::real > state( ncomp, 0.0 );
//    for (ncomp_t c=0; c<ncomp; ++c)
//    {
//      state[c] = unk[c][0];
//      state[c] += unk[c][1] * B_taylor[1]
//                + unk[c][2] * B_taylor[2]
//                + unk[c][3] * B_taylor[3];
//
//      if(ndof > 4)
//        state[c] += unk[c][4] * B_taylor[4] + unk[c][5] * B_taylor[5]
//                  + unk[c][6] * B_taylor[6] + unk[c][7] * B_taylor[7]
//                  + unk[c][8] * B_taylor[8] + unk[c][9] * B_taylor[9];
//    }
//
//    auto B = tk::eval_basis( ndof, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp] );
//
//    for (ncomp_t c=0; c<ncomp; ++c)
//    {
//      auto mark = c*ndof;
//      R[mark] += wt * state[c];
//
//      if(ndof > 1)
//      {
//        R[mark+1] += wt * state[c] * B[1];
//        R[mark+2] += wt * state[c] * B[2];
//        R[mark+3] += wt * state[c] * B[3];
//
//        if(ndof > 4)
//        {
//          R[mark+4] += wt * state[c] * B[4];
//          R[mark+5] += wt * state[c] * B[5];
//          R[mark+6] += wt * state[c] * B[6];
//          R[mark+7] += wt * state[c] * B[7];
//          R[mark+8] += wt * state[c] * B[8];
//          R[mark+9] += wt * state[c] * B[9];
//        }
//      }
//    }
//  }
//
//  for (ncomp_t c=0; c<ncomp; ++c)
//  {
//    auto mark = c*ndof;
//    for(std::size_t idof = 0; idof < ndof; idof++)
//      unk[c][idof] = R[mark+idof] / L[idof];
//  }
//}

std::vector< tk::real >
tk::eval_TaylorBasis( const std::size_t ndof,
                      const std::array< tk::real, 3 >& x,
                      const std::array< tk::real, 3 >& x_c,
                      const std::array< std::array< tk::real, 3>, 4 >& coordel )
// *****************************************************************************
//  Evaluate the Taylor basis at points
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] x Nodal coordinates
//! \param[in] x_c Coordinates of the centroid
//! \param[in] coordel Array of nodal coordinates for the tetrahedron
// *****************************************************************************
{
  std::vector< tk::real > avg( 6, 0.0 );
  if(ndof > 4)
  {
    Assert( ndof == 10, "Mismatch in number of degrees of freedom" );
    auto ng = tk::NGvol(ndof);

    std::array< std::vector< tk::real >, 3 > coordgp;
    std::vector< tk::real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    coordgp[2].resize( ng );
    wgp.resize( ng );

    tk::GaussQuadratureTet( ng, coordgp, wgp );

    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the coordinates of quadrature point at physical domain
      auto gp = tk::eval_gp( igp, coordel, coordgp );

      avg[0] += wgp[igp] * (gp[0] - x_c[0]) * (gp[0] - x_c[0]) * 0.5;
      avg[1] += wgp[igp] * (gp[1] - x_c[1]) * (gp[1] - x_c[1]) * 0.5;
      avg[2] += wgp[igp] * (gp[2] - x_c[2]) * (gp[2] - x_c[2]) * 0.5;
      avg[3] += wgp[igp] * (gp[0] - x_c[0]) * (gp[1] - x_c[1]);
      avg[4] += wgp[igp] * (gp[0] - x_c[0]) * (gp[2] - x_c[2]);
      avg[5] += wgp[igp] * (gp[1] - x_c[1]) * (gp[2] - x_c[2]);
    }
  }

  std::vector< tk::real > B( ndof, 1.0 );

  if(ndof > 1) {
    Assert( (ndof == 4)||(ndof == 10) ,
      "Mismatch in number of degrees of freedom" );
    B[1] = x[0] - x_c[0];
    B[2] = x[1] - x_c[1];
    B[3] = x[2] - x_c[2];
  }

  if(ndof > 4) {
    B[4] = B[1] * B[1] * 0.5 - avg[0];
    B[5] = B[2] * B[2] * 0.5 - avg[1];
    B[6] = B[3] * B[3] * 0.5 - avg[2];
    B[7] = B[1] * B[2] - avg[3];
    B[8] = B[1] * B[3] - avg[4];
    B[9] = B[2] * B[3] - avg[5];
  }

  return B;
}
