// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Basis.C
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

#include <array>

#include "Basis.h"

tk::real
tk::eval_det ( const std::size_t e,
               const std::vector< tk::real >& cx,
               const std::vector< tk::real >& cy,
               const std::vector< tk::real >& cz,
               const std::vector< std::size_t >& inpoel,
               std::array< std::array< tk::real, 3>, 4 >& coordel )
// *****************************************************************************
//  Compute the determinant of Jacobian matrix for the tetrahedron element
//! \param[in] e Index for the tetrahedron element
//! \param[in] cx Vector of x-coordinates of points
//! \param[in] cy Vector of y-coordinates of points
//! \param[in] cz Vector of z-coordinates of points
//! \param[in] inpoel Element-node connectivity
//! \param[in] coordel Array of nodal coordinates for tetrahedron element
//! \return The determinant of Jacobian matrix
// *****************************************************************************
{
  // nodal coordinates of the element
  coordel[0][0] = cx[ inpoel[4*e] ];
  coordel[0][1] = cy[ inpoel[4*e] ];
  coordel[0][2] = cz[ inpoel[4*e] ];

  coordel[1][0] = cx[ inpoel[4*e+1] ];
  coordel[1][1] = cy[ inpoel[4*e+1] ];
  coordel[1][2] = cz[ inpoel[4*e+1] ];

  coordel[2][0] = cx[ inpoel[4*e+2] ];
  coordel[2][1] = cy[ inpoel[4*e+2] ];
  coordel[2][2] = cz[ inpoel[4*e+2] ];

  coordel[3][0] = cx[ inpoel[4*e+3] ];
  coordel[3][1] = cy[ inpoel[4*e+3] ];
  coordel[3][2] = cz[ inpoel[4*e+3] ];

  return Jacobian( coordel[0], coordel[1], coordel[2], coordel[3] );
}

std::array< tk::real, 3 >
tk::eval_gp ( const std::size_t igp,
              const std::array< std::array< tk::real, 3>, 3 >& coordfa,
              const std::array< std::vector< tk::real >, 2 >& coordgp )
// *****************************************************************************
//  Compute the coordinates of quadrature points for face integral in physical space
//! \param[in] igp Index of quadrature points
//! \param[in] coordfa Array of nodal coordinates for face element
//! \param[in] coordgp Array of coordinates for quadrature points in reference space
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
//  Compute the coordinates of quadrature points for volume integral in physical space
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
  return {{ coord[0][0]*shp1 + coord[1][0]*shp2 + coord[2][0]*shp3 + coord[3][0]*shp4,
            coord[0][1]*shp1 + coord[1][1]*shp2 + coord[2][1]*shp3 + coord[3][1]*shp4,
            coord[0][2]*shp1 + coord[1][2]*shp2 + coord[2][2]*shp3 + coord[3][2]*shp4 }};
}

void
tk::eval_dBdx_p1( const std::size_t ndof,
                  const std::array< std::array< tk::real, 3 >, 3 >& jacInv, 
                  std::array< std::vector<tk::real>, 3 >& dBdx )
// *****************************************************************************
//  Compute the derivatives of basis function for DG(P1)
//! \param[in] ndof Number of degree of freedom
//! \param[in] jacInv Array of the inverse of Jacobian
//! \param[in,out] dBdx Array of the derivatives of basis function
// *****************************************************************************
{
  Assert( dBdx[0].size() == ndof, "Size mismatch" );
  Assert( dBdx[1].size() == ndof, "Size mismatch" );
  Assert( dBdx[2].size() == ndof, "Size mismatch" );

  // The derivatives of the basis functions dB/dx are easily calculated
  // via a transformation to the reference space as,
  // dB/dx = dB/dX . dx/dxi,
  // where, x = (x,y,z) are the physical coordinates, and
  //        xi = (xi, eta, zeta) are the reference coordinates.
  // The matrix dx/dxi is the inverse of the Jacobian of transformation
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

  dBdx[1][0] =  db2dxi1 * jacInv[0][0]
              + db2dxi2 * jacInv[1][0]
              + db2dxi3 * jacInv[2][0];

  dBdx[1][1] =  db2dxi1 * jacInv[0][1]
              + db2dxi2 * jacInv[1][1]
              + db2dxi3 * jacInv[2][1];

  dBdx[1][2] =  db2dxi1 * jacInv[0][2]
              + db2dxi2 * jacInv[1][2]
              + db2dxi3 * jacInv[2][2];

  dBdx[2][0] =  db3dxi1 * jacInv[0][0]
              + db3dxi2 * jacInv[1][0]
              + db3dxi3 * jacInv[2][0];

  dBdx[2][1] =  db3dxi1 * jacInv[0][1]
              + db3dxi2 * jacInv[1][1]
              + db3dxi3 * jacInv[2][1];

  dBdx[2][2] =  db3dxi1 * jacInv[0][2]
              + db3dxi2 * jacInv[1][2]
              + db3dxi3 * jacInv[2][2];

  dBdx[3][0] =  db4dxi1 * jacInv[0][0]
              + db4dxi2 * jacInv[1][0]
              + db4dxi3 * jacInv[2][0];

  dBdx[3][1] =  db4dxi1 * jacInv[0][1]
              + db4dxi2 * jacInv[1][1]
              + db4dxi3 * jacInv[2][1];

  dBdx[3][2] =  db4dxi1 * jacInv[0][2]
              + db4dxi2 * jacInv[1][2]
              + db4dxi3 * jacInv[2][2];
}

void
tk::eval_dBdx_p2( const std::size_t ndof,
                  const std::size_t igp,
                  const std::array< std::vector< tk::real >, 3 >& coordgp,
                  const std::array< std::array< tk::real, 3 >, 3 >& jacInv,
                  std::array< std::vector<tk::real>, 3 >& dBdx )
// *****************************************************************************
//  Compute the derivatives of basis function for DG(P2)
//! \param[in] ndof Number of degree of freedom
//! \param[in] igp Index of quadrature points
//! \param[in] coord Array of nodal coordinates for tetrahedron element
//! \param[in] jacInv Array of the inverse of Jacobian
//! \param[in,out] dBdx Array of the derivatives of basis function
// *****************************************************************************
{
  Assert( dBdx[0].size() == ndof, "Size mismatch" );
  Assert( dBdx[1].size() == ndof, "Size mismatch" );
  Assert( dBdx[2].size() == ndof, "Size mismatch" );

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

  dBdx[4][0] =  db5dxi1 * jacInv[0][0]
              + db5dxi2 * jacInv[1][0]
              + db5dxi3 * jacInv[2][0];

  dBdx[4][1] =  db5dxi1 * jacInv[0][1]
              + db5dxi2 * jacInv[1][1]
              + db5dxi3 * jacInv[2][1];

  dBdx[4][2] =  db5dxi1 * jacInv[0][2]
              + db5dxi2 * jacInv[1][2]
              + db5dxi3 * jacInv[2][2];

  dBdx[5][0] =  db6dxi1 * jacInv[0][0]
              + db6dxi2 * jacInv[1][0]
              + db6dxi3 * jacInv[2][0];

  dBdx[5][1] =  db6dxi1 * jacInv[0][1]
              + db6dxi2 * jacInv[1][1]
              + db6dxi3 * jacInv[2][1];

  dBdx[5][2] =  db6dxi1 * jacInv[0][2]
              + db6dxi2 * jacInv[1][2]
              + db6dxi3 * jacInv[2][2];

  dBdx[6][0] =  db7dxi1 * jacInv[0][0]
              + db7dxi2 * jacInv[1][0]
              + db7dxi3 * jacInv[2][0];

  dBdx[6][1] =  db7dxi1 * jacInv[0][1]
              + db7dxi2 * jacInv[1][1]
              + db7dxi3 * jacInv[2][1];

  dBdx[6][2] =  db7dxi1 * jacInv[0][2]
              + db7dxi2 * jacInv[1][2]
              + db7dxi3 * jacInv[2][2];

  dBdx[7][0] =  db8dxi1 * jacInv[0][0]
              + db8dxi2 * jacInv[1][0]
              + db8dxi3 * jacInv[2][0];

  dBdx[7][1] =  db8dxi1 * jacInv[0][1]
              + db8dxi2 * jacInv[1][1]
              + db8dxi3 * jacInv[2][1];

  dBdx[7][2] =  db8dxi1 * jacInv[0][2]
              + db8dxi2 * jacInv[1][2]
              + db8dxi3 * jacInv[2][2];

  dBdx[8][0] =  db9dxi1 * jacInv[0][0]
              + db9dxi2 * jacInv[1][0]
              + db9dxi3 * jacInv[2][0];

  dBdx[8][1] =  db9dxi1 * jacInv[0][1]
              + db9dxi2 * jacInv[1][1]
              + db9dxi3 * jacInv[2][1];

  dBdx[8][2] =  db9dxi1 * jacInv[0][2]
              + db9dxi2 * jacInv[1][2]
              + db9dxi3 * jacInv[2][2];

  dBdx[9][0] =  db10dxi1 * jacInv[0][0]
              + db10dxi2 * jacInv[1][0]
              + db10dxi3 * jacInv[2][0];

  dBdx[9][1] =  db10dxi1 * jacInv[0][1]
              + db10dxi2 * jacInv[1][1]
              + db10dxi3 * jacInv[2][1];

  dBdx[9][2] =  db10dxi1 * jacInv[0][2]
              + db10dxi2 * jacInv[1][2]
              + db10dxi3 * jacInv[2][2];
}

std::vector< tk::real >
tk::eval_basis( const std::size_t ndof,
                const std::array< std::array< tk::real, 3>, 4 >& coordel,
                const tk::real detT,
                const std::array < tk::real, 3 >& gp )
// *****************************************************************************
//  Compute the Dubiner basis functions for face integrals
//! \param[in] ndof Number of degree of freedom
//! \param[in] coordel Array of nodal coordinates for tetrahedron element
//! \param[in] detT Determination of Jacobian matrix for tetrahedron element
//! \param[in] gp Array of coordinates for quadrature points in physical space
//! \return Vector of basis functions
// *****************************************************************************
{
  // In order to determine the high-order solution from the left and right
  // elements at the surface quadrature points, the basis functions from
  // the left and right elements are needed. For this, a transformation to
  // the reference coordinates is necessary, since the basis functions are
  // defined on the reference tetrahedron only.

  tk::real detT_gp = 0.0;

  // Transformation of the physical coordinates of the quadrature point
  // to reference space for the element to be able to compute basis functions.
  detT_gp = Jacobian( coordel[0], gp, coordel[2], coordel[3] );
  auto xi = detT_gp / detT;
  detT_gp = Jacobian( coordel[0], coordel[1], gp, coordel[3] );
  auto eta = detT_gp / detT;
  detT_gp = Jacobian( coordel[0], coordel[1], coordel[2], gp );
  auto zeta = detT_gp / detT;

  // Array of basis functions
  std::vector< tk::real > B( ndof, 1.0 );

  if ( ndof > 1 )           // DG(P1)
  {
    // Basis functions (DGP1) at igp
    B[1] = 2.0 * xi + eta + zeta - 1.0;
    B[2] = 3.0 * eta + zeta - 1.0;
    B[3] = 4.0 * zeta - 1.0;

    if( ndof > 4 )         // DG(P2)
    {
      // Basis functions (DGP2) at igp
      auto xi_xi = xi * xi;
      auto xi_eta = xi * eta;
      auto xi_zeta = xi * zeta;
      auto eta_eta = eta * eta;
      auto eta_zeta = eta * zeta;
      auto zeta_zeta = zeta * zeta;

      B[4] = 6.0 * xi_xi + eta_eta + zeta_zeta
           + 6.0 * xi_eta + 6.0 * xi_zeta + 2.0 * eta_zeta
           - 6.0 * xi - 2.0 * eta - 2.0 * zeta + 1.0;
      B[5] = 5.0 * eta_eta + zeta_zeta
           + 10.0 * xi_eta + 2.0 * xi_zeta + 6.0 * eta_zeta
           - 2.0 * xi - 6.0 * eta - 2.0 * zeta + 1.0;
      B[6] = 6.0 * zeta_zeta + 12.0 * xi_zeta + 6.0 * eta_zeta
           - 2.0 * xi - eta - 7.0 * zeta + 1.0;
      B[7] = 10.0 * eta_eta + zeta_zeta + 8.0 * eta_zeta
           - 8.0 * eta - 2.0 * zeta + 1.0;
      B[8] = 6.0 * zeta_zeta + 18.0 * eta_zeta - 3.0 * eta
           - 7.0 * zeta + 1.0;
      B[9] = 15.0 * zeta_zeta - 10.0 * zeta + 1.0;
    }
  }

  return B;
}

std::vector< tk::real >
tk::eval_basis( const std::size_t ndof,
                const std::size_t igp,
                const std::array< std::vector< tk::real >, 3 >& coordgp )
// *****************************************************************************
//  Compute the Dubiner basis functions for volume integrals
//! \param[in] ndof Number of degree of freedom
//! \param[in] igp Index of gauss points
//! \param[in] coordgp Array of coordinates for quadrature points
//! \return Vector of basis functions
// *****************************************************************************
{
  // Array of basis functions
  std::vector< tk::real > B( ndof, 1.0 );

  if ( ndof > 1 )           // DG(P1)
  {
    auto xi   = coordgp[0][igp];
    auto eta  = coordgp[1][igp];
    auto zeta = coordgp[2][igp];

    // Basis functions (DGP1) at igp
    B[1] = 2.0 * xi + eta + zeta - 1.0;
    B[2] = 3.0 * eta + zeta - 1.0;
    B[3] = 4.0 * zeta - 1.0;

    if( ndof > 4 )         // DG(P2)
    {
      auto xi_xi   = coordgp[0][igp] * coordgp[0][igp];
      auto xi_eta  = coordgp[0][igp] * coordgp[1][igp];
      auto xi_zeta = coordgp[0][igp] * coordgp[2][igp];

      auto eta_eta  = coordgp[1][igp] * coordgp[1][igp];
      auto eta_zeta = coordgp[1][igp] * coordgp[2][igp];

      auto zeta_zeta = coordgp[2][igp] * coordgp[2][igp];

      // Basis functions (DGP2) at igp
      B[4] =  6.0 * xi_xi + eta_eta + zeta_zeta + 6.0 * xi_eta + 6.0 * xi_zeta + 2.0 * eta_zeta
            - 6.0 * xi - 2.0 * eta - 2.0 * zeta + 1.0;
      B[5] =  5.0 * eta_eta + zeta_zeta + 10.0 * xi_eta + 2.0 * xi_zeta + 6.0 * eta_zeta - 2.0 * xi
            - 6.0 * eta - 2.0 * zeta + 1.0;
      B[6] =  6.0 * zeta_zeta + 12.0 * xi_zeta + 6.0 * eta_zeta - 2.0 * xi
            - eta - 7.0 * zeta + 1.0;
      B[7] =  10.0 * eta_eta + zeta_zeta + 8.0 * eta_zeta - 8.0 * eta - 2.0 * zeta + 1.0;
      B[8] =  6.0 * zeta_zeta + 18.0 * eta_zeta - 3.0 * eta - 7.0 * zeta + 1.0;
      B[9] =  15.0 * zeta_zeta - 10.0 * zeta + 1.0;
    } 
  }

  return B;
}

std::vector< tk::real >
tk::eval_state ( ncomp_t ncomp,
                 ncomp_t offset,
                 const std::size_t ndof,
                 const std::size_t e,
                 const Fields& U,
                 const Fields& limFunc,
                 const std::vector< tk::real >& B )
// *****************************************************************************
//  Compute the state variables for the tetrahedron element
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] ndof Number of degree of freedom
//! \param[in] e Index for the tetrahedron element
//! \param[in] U Solution vector at recent time step
//! \param[in] limFunc Limiter function for higher-order solution dofs
//! \param[in] B Vector of basis functions
//! \return Vector of state variable for tetrahedron element
// *****************************************************************************
{
  Assert( B.size() == ndof, "Size mismatch" );

  // Array of state variable for tetrahedron element
  std::vector< tk::real > state( ncomp );

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    state[c] = U( e, mark, offset );

    if(ndof > 1)        //DG(P1)
    {
      auto lmark = c*(ndof-1);
      state[c] += limFunc( e, lmark  , 0 ) * U( e, mark+1, offset ) * B[1]
                + limFunc( e, lmark+1, 0 ) * U( e, mark+2, offset ) * B[2]
                + limFunc( e, lmark+2, 0 ) * U( e, mark+3, offset ) * B[3];
    }

    if(ndof > 4)        //DG(P2)
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
