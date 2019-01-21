// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Basis.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for computing the Dubiner basis functions in DG methods
  \details   This file contains functionality for computing the Dubiner basis
     functions and relating coordinates transformation functions used in 
     discontinuous Galerkin methods for variaous orders of numerical representation.
*/
// *****************************************************************************

#include <array>

#include "Basis.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

tk::real
tk::eval_det ( const std::size_t e,
               const std::vector< tk::real >& cx,
               const std::vector< tk::real >& cy,
               const std::vector< tk::real >& cz,
               const std::vector< std::size_t >& inpoel,
               std::array< std::array< tk::real, 3>, 4 >& coordel )
// *****************************************************************************
//  Compute the determination of Jacobian matrix for the tetrahedron element
//! \param[in] e Index for the tetrahedron element
//! \param[in] cx Vector of x-coordinates of points
//! \param[in] cy Vector of y-coordinates of points
//! \param[in] cz Vector of z-coordinates of points
//! \param[in] inpoel Element-node connectivity
//! \param[in] coordel Array of nodal coordinates for tetrahedron element
//! \return the determination of Jacobian matrix
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

void
tk::eval_gp ( const std::size_t igp,
              const std::array< std::array< tk::real, 3>, 3 >& coordfa,
              const std::array< std::vector< tk::real >, 2 >& coordgp,
              std::array < tk::real, 3 >& gp )
// *****************************************************************************
//  Compute the coordinates of quadrature points in physical domain
//! \param[in] igp Index of quadrature points
//! \param[in] coordfa Array of nodal coordinates for face element
//! \param[in] coordgp Array of coordinates for quadrature points in reference domain
//! \param[in,out] gp Array of coordinates for quadrature points in physical domian
// *****************************************************************************
{
  // Barycentric coordinates for the triangular face
  auto shp1 = 1.0 - coordgp[0][igp] - coordgp[1][igp];
  auto shp2 = coordgp[0][igp];
  auto shp3 = coordgp[1][igp];

  // transformation of the quadrature point from the 2D reference/master
  // element to physical domain, to obtain its physical (x,y,z)
  // coordinates.
  gp[0] = coordfa[0][0]*shp1 + coordfa[1][0]*shp2 + coordfa[2][0]*shp3;
  gp[1] = coordfa[0][1]*shp1 + coordfa[1][1]*shp2 + coordfa[2][1]*shp3;
  gp[2] = coordfa[0][2]*shp1 + coordfa[1][2]*shp2 + coordfa[2][2]*shp3;
}

void
tk::eval_xi ( const std::array< std::array< tk::real, 3>, 4 >& coordel,
              const tk::real detT,
              const std::array < tk::real, 3 >& gp,
              tk::real& xi,
              tk::real& eta,
              tk::real& zeta )
// *****************************************************************************
//  Compute the coordinates of quadrature points in reference domian
//! \param[in] coordel Array of nodal coordinates for tetrahedron element
//! \param[in] detT Determination of Jacobian matrix for tetrahedron element
//! \param[in] gp Array of coordinates for quadrature points in physical domain
//! \param[in] xi, eta, zeta Coordinates of quadrature points in reference domain
// *****************************************************************************
{
  // The basis functions chosen for the DG method are the Dubiner
  // basis, which are Legendre polynomials modified for tetrahedra,
  // which are defined only on the reference/master tetrahedron.
  // Thus, to determine the high-order solution from the left and right
  // elements at the surface quadrature points, the basis functions
  // from the left and right elements are needed. For this, a
  // transformation to the reference coordinates is necessary, since
  // the basis functions are defined on the reference tetrahedron only.
  // Ref: [1] https://doi.org/10.1007/BF01060030
  //      [2] https://doi.org/10.1093/imamat/hxh111

  tk::real detT_gp = 0.0;

  // transformation of the physical coordinates of the quadrature point
  // to reference space for the element to be able to compute basis functions.
  detT_gp = Jacobian( coordel[0], gp, coordel[2], coordel[3] );
  xi = detT_gp / detT;
  detT_gp = Jacobian( coordel[0], coordel[1], gp, coordel[3] );
  eta = detT_gp / detT;
  detT_gp = Jacobian( coordel[0], coordel[1], coordel[2], gp );
  zeta = detT_gp / detT;
}

void
tk::eval_basis( const tk::real xi, 
                const tk::real eta, 
                const tk::real zeta,
                std::array< tk::real, 10>& B )
// *****************************************************************************
//  Compute the Dubiner basis functions
//! \param[in] xi, eta, zeta Coordinates of quadrature points in reference domain
//! \param[in] B Array of basis functions
// *****************************************************************************
{
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();

  // basis functions (DGP1) at igp for the left element
  B[1] = 2.0 * xi + eta + zeta - 1.0;
  B[2] = 3.0 * eta + zeta - 1.0;
  B[3] = 4.0 * zeta - 1.0;

  if(ndof > 4)        //DG(P2)
  {
    // basis functions (DGP2) at igp for the left element
    auto xi_xi = xi * xi;
    auto xi_eta = xi * eta;
    auto xi_zeta = xi * zeta;
    auto eta_eta = eta * eta;
    auto eta_zeta = eta * zeta;
    auto zeta_zeta = zeta * zeta;

    B[4] =  6.0 * xi_xi + eta_eta + zeta_zeta
          + 6.0 * xi_eta + 6.0 * xi_zeta + 2.0 * eta_zeta
          - 6.0 * xi - 2.0 * eta - 2.0 * zeta + 1.0;
    B[5] =  5.0 * eta_eta + zeta_zeta
          + 10.0 * xi_eta + 2.0 * xi_zeta + 6.0 * eta_zeta
          - 2.0 * xi - 6.0 * eta - 2.0 * zeta + 1.0;
    B[6] =  6.0 * zeta_zeta + 12.0 * xi_zeta + 6.0 * eta_zeta
          - 2.0 * xi - eta - 7.0 * zeta + 1.0;
    B[7] =  10.0 * eta_eta + zeta_zeta + 8.0 * eta_zeta
          - 8.0 * eta - 2.0 * zeta + 1.0;
    B[8] =  6.0 * zeta_zeta + 18.0 * eta_zeta - 3.0 * eta
          - 7.0 * zeta + 1.0;
    B[9] =  15.0 * zeta_zeta - 10.0 * zeta + 1.0;
  }
}

void
tk::eval_state ( ncomp_t ncomp,
                 ncomp_t offset,
                 const std::size_t e,
                 const Fields& U,
                 const Fields& limFunc,
                 std::array< tk::real, 10>& B,
                 std::vector< tk::real >& state )
// *****************************************************************************
//  Compute the state variables for the tetrahedron element
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] e Index for the tetrahedron element
//! \param[in] U Solution vector at recent time step
//! \param[in] limFunc Limiter function for higher-order solution dofs
//! \param[in] B Array of basis functions
//! \param[in,out] state Array of state variable for tetrahedron element
// *****************************************************************************
{
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();

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
}
