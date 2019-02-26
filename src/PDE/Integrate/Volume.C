// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Volume.C
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Functions for computing volume integrals for a system of PDEs in DG
     methods
  \details   This file contains functionality for computing volume integrals for
     a system of PDEs used in discontinuous Galerkin methods for various orders
     of numerical representation.
*/
// *****************************************************************************

#include "Volume.h"
#include "Vector.h"
#include "Quadrature.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

void
tk::volIntP1( ncomp_t system,
              ncomp_t ncomp,
              ncomp_t offset,
              const std::vector< std::size_t >& inpoel,
              const UnsMesh::Coords& coord,
              const Fields& geoElem,
              const FluxFn& flux,
              const VelFn& vel,
              const Fields& U,
              const Fields& limFunc,
              Fields& R )
// *****************************************************************************
//  Compute volume integrals for DG(P1)
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] flux Flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] U Solution vector at recent time step
//! \param[in] limFunc Limiter function for higher-order solution dofs
//! \param[in,out] R Right-hand side vector added to
// *****************************************************************************
{
  using inciter::g_inputdeck;
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();

  // Number of integration points
  constexpr std::size_t NG = 5;

  // arrays for quadrature points
  std::array< std::array< real, NG >, 3 > coordgp;
  std::array< real, NG > wgp;

  // get quadrature point weights and coordinates for tetrahedron
  GaussQuadratureTet( coordgp, wgp );

  std::array< std::array< real, 3 >, 3 > jacInv;

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  // compute volume integrals
  for (std::size_t e=0; e<U.nunk(); ++e)
  {
    auto x1 = cx[ inpoel[4*e]   ];
    auto y1 = cy[ inpoel[4*e]   ];
    auto z1 = cz[ inpoel[4*e]   ];

    auto x2 = cx[ inpoel[4*e+1] ];
    auto y2 = cy[ inpoel[4*e+1] ];
    auto z2 = cz[ inpoel[4*e+1] ];

    auto x3 = cx[ inpoel[4*e+2] ];
    auto y3 = cy[ inpoel[4*e+2] ];
    auto z3 = cz[ inpoel[4*e+2] ];

    auto x4 = cx[ inpoel[4*e+3] ];
    auto y4 = cy[ inpoel[4*e+3] ];
    auto z4 = cz[ inpoel[4*e+3] ];

    jacInv = inverseJacobian( {{x1, y1, z1}},
                              {{x2, y2, z2}},
                              {{x3, y3, z3}},
                              {{x4, y4, z4}} );

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

    auto db2dx =  db2dxi1 * jacInv[0][0]
                + db2dxi2 * jacInv[1][0]
                + db2dxi3 * jacInv[2][0];

    auto db2dy =  db2dxi1 * jacInv[0][1]
                + db2dxi2 * jacInv[1][1]
                + db2dxi3 * jacInv[2][1];

    auto db2dz =  db2dxi1 * jacInv[0][2]
                + db2dxi2 * jacInv[1][2]
                + db2dxi3 * jacInv[2][2];

    auto db3dx =  db3dxi1 * jacInv[0][0]
                + db3dxi2 * jacInv[1][0]
                + db3dxi3 * jacInv[2][0];

    auto db3dy =  db3dxi1 * jacInv[0][1]
                + db3dxi2 * jacInv[1][1]
                + db3dxi3 * jacInv[2][1];

    auto db3dz =  db3dxi1 * jacInv[0][2]
                + db3dxi2 * jacInv[1][2]
                + db3dxi3 * jacInv[2][2];

    auto db4dx =  db4dxi1 * jacInv[0][0]
                + db4dxi2 * jacInv[1][0]
                + db4dxi3 * jacInv[2][0];

    auto db4dy =  db4dxi1 * jacInv[0][1]
                + db4dxi2 * jacInv[1][1]
                + db4dxi3 * jacInv[2][1];

    auto db4dz =  db4dxi1 * jacInv[0][2]
                + db4dxi2 * jacInv[1][2]
                + db4dxi3 * jacInv[2][2];

    // Gaussian quadrature
    for (std::size_t igp=0; igp<NG; ++igp)
    {
      auto B2 = 2.0 * coordgp[0][igp] + coordgp[1][igp] + coordgp[2][igp] - 1.0;
      auto B3 = 3.0 * coordgp[1][igp] + coordgp[2][igp] - 1.0;
      auto B4 = 4.0 * coordgp[2][igp] - 1.0;

      auto wt = wgp[igp] * geoElem(e, 0, 0);

      auto shp1 = 1.0 - coordgp[0][igp] - coordgp[1][igp] - coordgp[2][igp];
      auto shp2 = coordgp[0][igp];
      auto shp3 = coordgp[1][igp];
      auto shp4 = coordgp[2][igp];

      std::array< tk::real,3  > gp{{ x1*shp1 + x2*shp2 + x3*shp3 + x4*shp4,
                                     y1*shp1 + y2*shp2 + y3*shp3 + y4*shp4,
                                     z1*shp1 + z2*shp2 + z3*shp3 + z4*shp4 }};

      std::vector< real > ugp( ncomp );
      for (ncomp_t c=0; c<ncomp; ++c) {
        auto mark = c*ndof;
        auto lmark = c*(ndof-1);
        ugp[c] = U(e, mark, offset)
                 + limFunc(e, lmark+0, 0) * U(e, mark+1, offset) * B2
                 + limFunc(e, lmark+1, 0) * U(e, mark+2, offset) * B3
                 + limFunc(e, lmark+2, 0) * U(e, mark+3, offset) * B4;
      }

     // evaluate prescribed velocity (if any)
      auto v = vel( system, ncomp, gp[0], gp[1], gp[2] );
      // comput flux
      auto fl = flux( system, ncomp, ugp, v );

      for (ncomp_t c=0; c<ncomp; ++c)
      {
        auto mark = c*ndof;
        R(e, mark+1, offset) +=
          wt * (fl[c][0]*db2dx + fl[c][1]*db2dy + fl[c][2]*db2dz);
        R(e, mark+2, offset) +=
          wt * (fl[c][0]*db3dx + fl[c][1]*db3dy + fl[c][2]*db3dz);
        R(e, mark+3, offset) +=
          wt * (fl[c][0]*db4dx + fl[c][1]*db4dy + fl[c][2]*db4dz);
      }
    }
  }
}

void
tk::volIntP2( ncomp_t system,
              ncomp_t ncomp,
              ncomp_t offset,
              const std::vector< std::size_t >& inpoel,
              const UnsMesh::Coords& coord,
              const Fields& geoElem,
              const FluxFn& flux,
              const VelFn& vel,
              const Fields& U,
              Fields& R )
// *****************************************************************************
//  Compute volume integrals for DG(P2)
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] flux Flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] U Solution vector at recent time step
//! \param[in,out] R Right-hand side vector added to
// *****************************************************************************
{
  using inciter::g_inputdeck;
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();

  // Number of integration points
  constexpr std::size_t NG = 11;

  // arrays for quadrature points
  std::array< std::array< real, NG >, 3 > coordgp;
  std::array< real, NG > wgp;

  // get quadrature point weights and coordinates for tetrahedron
  GaussQuadratureTet( coordgp, wgp );

  std::array< std::array< real, 3 >, 3 > jacInv;

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  // compute volume integrals
  for (std::size_t e=0; e<U.nunk(); ++e)
  {
    auto x1 = cx[ inpoel[4*e]   ];
    auto y1 = cy[ inpoel[4*e]   ];
    auto z1 = cz[ inpoel[4*e]   ];

    auto x2 = cx[ inpoel[4*e+1] ];
    auto y2 = cy[ inpoel[4*e+1] ];
    auto z2 = cz[ inpoel[4*e+1] ];

    auto x3 = cx[ inpoel[4*e+2] ];
    auto y3 = cy[ inpoel[4*e+2] ];
    auto z3 = cz[ inpoel[4*e+2] ];

    auto x4 = cx[ inpoel[4*e+3] ];
    auto y4 = cy[ inpoel[4*e+3] ];
    auto z4 = cz[ inpoel[4*e+3] ];

    jacInv = inverseJacobian( {{x1, y1, z1}},
                              {{x2, y2, z2}},
                              {{x3, y3, z3}},
                              {{x4, y4, z4}} );

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

    auto db2dx =  db2dxi1 * jacInv[0][0]
                + db2dxi2 * jacInv[1][0]
                + db2dxi3 * jacInv[2][0];

    auto db2dy =  db2dxi1 * jacInv[0][1]
                + db2dxi2 * jacInv[1][1]
                + db2dxi3 * jacInv[2][1];

    auto db2dz =  db2dxi1 * jacInv[0][2]
                + db2dxi2 * jacInv[1][2]
                + db2dxi3 * jacInv[2][2];

    auto db3dx =  db3dxi1 * jacInv[0][0]
                + db3dxi2 * jacInv[1][0]
                + db3dxi3 * jacInv[2][0];

    auto db3dy =  db3dxi1 * jacInv[0][1]
                + db3dxi2 * jacInv[1][1]
                + db3dxi3 * jacInv[2][1];

    auto db3dz =  db3dxi1 * jacInv[0][2]
                + db3dxi2 * jacInv[1][2]
                + db3dxi3 * jacInv[2][2];

    auto db4dx =  db4dxi1 * jacInv[0][0]
                + db4dxi2 * jacInv[1][0]
                + db4dxi3 * jacInv[2][0];

    auto db4dy =  db4dxi1 * jacInv[0][1]
                + db4dxi2 * jacInv[1][1]
                + db4dxi3 * jacInv[2][1];

    auto db4dz =  db4dxi1 * jacInv[0][2]
                + db4dxi2 * jacInv[1][2]
                + db4dxi3 * jacInv[2][2];

    // Gaussian quadrature
    for (std::size_t igp=0; igp<NG; ++igp)
    {
      // Continue to calculate the derivatives of the basis functions dB/dx for DG(P2)
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

      auto db5dx =  db5dxi1 * jacInv[0][0]
                  + db5dxi2 * jacInv[1][0]
                  + db5dxi3 * jacInv[2][0];

      auto db5dy =  db5dxi1 * jacInv[0][1]
                  + db5dxi2 * jacInv[1][1]
                  + db5dxi3 * jacInv[2][1];

      auto db5dz =  db5dxi1 * jacInv[0][2]
                  + db5dxi2 * jacInv[1][2]
                  + db5dxi3 * jacInv[2][2];

      auto db6dx =  db6dxi1 * jacInv[0][0]
                  + db6dxi2 * jacInv[1][0]
                  + db6dxi3 * jacInv[2][0];

      auto db6dy =  db6dxi1 * jacInv[0][1]
                  + db6dxi2 * jacInv[1][1]
                  + db6dxi3 * jacInv[2][1];

      auto db6dz =  db6dxi1 * jacInv[0][2]
                  + db6dxi2 * jacInv[1][2]
                  + db6dxi3 * jacInv[2][2];

      auto db7dx =  db7dxi1 * jacInv[0][0]
                  + db7dxi2 * jacInv[1][0]
                  + db7dxi3 * jacInv[2][0];

      auto db7dy =  db7dxi1 * jacInv[0][1]
                  + db7dxi2 * jacInv[1][1]
                  + db7dxi3 * jacInv[2][1];

      auto db7dz =  db7dxi1 * jacInv[0][2]
                  + db7dxi2 * jacInv[1][2]
                  + db7dxi3 * jacInv[2][2];

      auto db8dx =  db8dxi1 * jacInv[0][0]
                  + db8dxi2 * jacInv[1][0]
                  + db8dxi3 * jacInv[2][0];

      auto db8dy =  db8dxi1 * jacInv[0][1]
                  + db8dxi2 * jacInv[1][1]
                  + db8dxi3 * jacInv[2][1];

      auto db8dz =  db8dxi1 * jacInv[0][2]
                  + db8dxi2 * jacInv[1][2]
                  + db8dxi3 * jacInv[2][2];

      auto db9dx =  db9dxi1 * jacInv[0][0]
                  + db9dxi2 * jacInv[1][0]
                  + db9dxi3 * jacInv[2][0];

      auto db9dy =  db9dxi1 * jacInv[0][1]
                  + db9dxi2 * jacInv[1][1]
                  + db9dxi3 * jacInv[2][1];

      auto db9dz =  db9dxi1 * jacInv[0][2]
                  + db9dxi2 * jacInv[1][2]
                  + db9dxi3 * jacInv[2][2];

      auto db10dx =  db10dxi1 * jacInv[0][0]
                   + db10dxi2 * jacInv[1][0]
                   + db10dxi3 * jacInv[2][0];

      auto db10dy =  db10dxi1 * jacInv[0][1]
                   + db10dxi2 * jacInv[1][1]
                   + db10dxi3 * jacInv[2][1];

      auto db10dz =  db10dxi1 * jacInv[0][2]
                   + db10dxi2 * jacInv[1][2]
                   + db10dxi3 * jacInv[2][2];

      // compute the basis function for DG(P2)
      auto xi_xi   = coordgp[0][igp] * coordgp[0][igp];
      auto xi_eta  = coordgp[0][igp] * coordgp[1][igp];
      auto xi_zeta = coordgp[0][igp] * coordgp[2][igp];

      auto eta_eta  = coordgp[1][igp] * coordgp[1][igp];
      auto eta_zeta = coordgp[1][igp] * coordgp[2][igp];

      auto zeta_zeta = coordgp[2][igp] * coordgp[2][igp];

      auto xi   = coordgp[0][igp];
      auto eta  = coordgp[1][igp];
      auto zeta = coordgp[2][igp];

      auto B2 = 2.0 * xi + eta + zeta - 1.0;
      auto B3 = 3.0 * eta + zeta - 1.0;
      auto B4 = 4.0 * zeta - 1.0;
      auto B5 = 6.0 * xi_xi + eta_eta + zeta_zeta + 6.0 * xi_eta + 6.0 * xi_zeta
              + 2.0 * eta_zeta - 6.0 * xi - 2.0 * eta - 2.0 * zeta + 1.0;
      auto B6 = 5.0 * eta_eta + zeta_zeta + 10.0 * xi_eta + 2.0 * xi_zeta
              + 6.0 * eta_zeta - 2.0 * xi - 6.0 * eta - 2.0 * zeta + 1.0;
      auto B7 = 6.0 * zeta_zeta + 12.0 * xi_zeta + 6.0 * eta_zeta
              - 2.0 * xi - eta - 7.0 * zeta + 1.0;
      auto B8 = 10.0 * eta_eta + zeta_zeta + 8.0 * eta_zeta
              - 8.0 * eta - 2.0 * zeta + 1.0;
      auto B9 = 6.0 * zeta_zeta + 18.0 * eta_zeta - 3.0 * eta
              - 7.0 * zeta + 1.0;
      auto B10 = 15.0 * zeta_zeta - 10.0 * zeta + 1.0;

      auto wt = wgp[igp] * geoElem(e, 0, 0);

      auto shp1 = 1.0 - xi - eta - zeta;
      auto shp2 = xi;
      auto shp3 = eta;
      auto shp4 = zeta;

      std::array< tk::real,3  > gp{{ x1*shp1 + x2*shp2 + x3*shp3 + x4*shp4,
                                     y1*shp1 + y2*shp2 + y3*shp3 + y4*shp4,
                                     z1*shp1 + z2*shp2 + z3*shp3 + z4*shp4 }};

      std::vector< real > ugp( ncomp );
      for (ncomp_t c=0; c<ncomp; ++c) {
        auto mark = c*ndof;
        ugp[c] =  U(e, mark, offset)
                + U(e, mark+1, offset) * B2
                + U(e, mark+2, offset) * B3
                + U(e, mark+3, offset) * B4
                + U(e, mark+4, offset) * B5
                + U(e, mark+5, offset) * B6
                + U(e, mark+6, offset) * B7
                + U(e, mark+7, offset) * B8
                + U(e, mark+8, offset) * B9
                + U(e, mark+9, offset) * B10;
      }

     // evaluate prescribed velocity (if any)
      auto v = vel( system, ncomp, gp[0], gp[1], gp[2] );
      // comput flux
      auto fl = flux( system, ncomp, ugp, v );

      for (ncomp_t c=0; c<ncomp; ++c)
      {
        auto mark = c*ndof;
        R(e, mark+1, offset) +=
          wt * (fl[c][0]*db2dx + fl[c][1]*db2dy + fl[c][2]*db2dz);
        R(e, mark+2, offset) +=
          wt * (fl[c][0]*db3dx + fl[c][1]*db3dy + fl[c][2]*db3dz);
        R(e, mark+3, offset) +=
          wt * (fl[c][0]*db4dx + fl[c][1]*db4dy + fl[c][2]*db4dz);
        R(e, mark+4, offset) += 
          wt * (fl[c][0]*db5dx + fl[c][1]*db5dy + fl[c][2]*db5dz);
        R(e, mark+5, offset) += 
          wt * (fl[c][0]*db6dx + fl[c][1]*db6dy + fl[c][2]*db6dz);
        R(e, mark+6, offset) += 
          wt * (fl[c][0]*db7dx + fl[c][1]*db7dy + fl[c][2]*db7dz);
        R(e, mark+7, offset) += 
          wt * (fl[c][0]*db8dx + fl[c][1]*db8dy + fl[c][2]*db8dz);
        R(e, mark+8, offset) += 
          wt * (fl[c][0]*db9dx + fl[c][1]*db9dy + fl[c][2]*db9dz);
        R(e, mark+9, offset) += 
          wt * (fl[c][0]*db10dx + fl[c][1]*db10dy + fl[c][2]*db10dz);
      }
    }
  }
}
