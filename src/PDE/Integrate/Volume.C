// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Volume.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
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
