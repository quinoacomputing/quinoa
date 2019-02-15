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
tk::volInt( ncomp_t system,
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
//  Compute volume integrals for DG
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

  // Number of quadrature points for volume integration
  auto ng = tk::NGvol(ndof);

  // arrays for quadrature points
  std::array< std::vector< real >, 3 > coordgp;
  std::vector< real > wgp;

  coordgp[0].resize( ng );
  coordgp[1].resize( ng );
  coordgp[2].resize( ng );
  wgp.resize( ng );

  // get quadrature point weights and coordinates for triangle
  GaussQuadratureTet( ng, coordgp, wgp );

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  // Nodal Coordinates of the tetrahedron element
  std::array< std::array< real, 3>, 4 > coordel;

  std::array< std::array< real, 3 >, 3 > jacInv;

  // compute volume integrals
  for (std::size_t e=0; e<U.nunk(); ++e)
  {
    coordel[0][0] = cx[ inpoel[4*e]   ];
    coordel[0][1] = cy[ inpoel[4*e]   ];
    coordel[0][2] = cz[ inpoel[4*e]   ];

    coordel[1][0] = cx[ inpoel[4*e+1] ];
    coordel[1][1] = cy[ inpoel[4*e+1] ];
    coordel[1][2] = cz[ inpoel[4*e+1] ];

    coordel[2][0] = cx[ inpoel[4*e+2] ];
    coordel[2][1] = cy[ inpoel[4*e+2] ];
    coordel[2][2] = cz[ inpoel[4*e+2] ];

    coordel[3][0] = cx[ inpoel[4*e+3] ];
    coordel[3][1] = cy[ inpoel[4*e+3] ];
    coordel[3][2] = cz[ inpoel[4*e+3] ];

    jacInv = inverseJacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

    std::array< std::vector<tk::real>, 3 > dBdx;
    dBdx[0].resize( ndof, 0 );
    dBdx[1].resize( ndof, 0 );
    dBdx[2].resize( ndof, 0 );

    // Compute the derivatives of basis function for DG(P1)
    eval_dBdx_p1( jacInv, dBdx );

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      if (ndof > 4)
        eval_dBdx_p2( igp, coordgp, jacInv, dBdx );

      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordel, coordgp );

      // Compute the basis function
      auto B =
        eval_basis( ndof, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp] );

      auto wt = wgp[igp] * geoElem(e, 0, 0);

      auto state = eval_state( ncomp, offset, ndof, e, U, limFunc, B );

      // evaluate prescribed velocity (if any)
      auto v = vel( system, ncomp, gp[0], gp[1], gp[2] );

      // comput flux
      auto fl = flux( system, ncomp, state, v );

      update_rhs( ncomp, offset, ndof, wt, e, dBdx, fl, R );
    }
  }
}

void
tk::update_rhs( ncomp_t ncomp,
                ncomp_t offset,
                const std::size_t ndof,
                const tk::real wt,
                const std::size_t e,
                const std::array< std::vector<tk::real>, 3 >& dBdx,
                const std::vector< std::array< tk::real, 3 > >& fl,
                Fields& R )
// *****************************************************************************
//  Update the rhs by adding the source term integrals
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] ndof Number of degree of freedom
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] e Element index
//! \param[in] B Vector of basis functions
//! \param[in] fl Vector of numerical flux
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  Assert( dBdx[0].size() == ndof, "Size mismatch for basis function derivatives" );
  Assert( dBdx[1].size() == ndof, "Size mismatch for basis function derivatives" );
  Assert( dBdx[2].size() == ndof, "Size mismatch for basis function derivatives" );
  Assert( fl.size() == ncomp, "Size mismatch for flux term" );

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    R(e, mark+1, offset) +=
      wt * (fl[c][0]*dBdx[0][1] + fl[c][1]*dBdx[1][1] + fl[c][2]*dBdx[2][1]);
    R(e, mark+2, offset) +=
      wt * (fl[c][0]*dBdx[0][2] + fl[c][1]*dBdx[1][2] + fl[c][2]*dBdx[2][2]);
    R(e, mark+3, offset) +=
      wt * (fl[c][0]*dBdx[0][3] + fl[c][1]*dBdx[1][3] + fl[c][2]*dBdx[2][3]);

    if( ndof > 4 )
    {
      R(e, mark+4, offset) +=
        wt * (fl[c][0]*dBdx[0][4] + fl[c][1]*dBdx[1][4] + fl[c][2]*dBdx[2][4]);
      R(e, mark+5, offset) +=
        wt * (fl[c][0]*dBdx[0][5] + fl[c][1]*dBdx[1][5] + fl[c][2]*dBdx[2][5]);
      R(e, mark+6, offset) +=
        wt * (fl[c][0]*dBdx[0][6] + fl[c][1]*dBdx[1][6] + fl[c][2]*dBdx[2][6]);
      R(e, mark+7, offset) +=
        wt * (fl[c][0]*dBdx[0][7] + fl[c][1]*dBdx[1][7] + fl[c][2]*dBdx[2][7]);
      R(e, mark+8, offset) +=
        wt * (fl[c][0]*dBdx[0][8] + fl[c][1]*dBdx[1][8] + fl[c][2]*dBdx[2][8]);
      R(e, mark+9, offset) +=
        wt * (fl[c][0]*dBdx[0][9] + fl[c][1]*dBdx[1][9] + fl[c][2]*dBdx[2][9]);
    }
  }
}
