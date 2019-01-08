// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Source.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for computing integrals of an arbitrary source term of a
     system of PDEs in DG methods
  \details   This file contains functionality for computing integrals of an
     arbitrary source term of a system of PDEs used in discontinuous Galerkin
     methods for various orders of numerical representation.
*/
// *****************************************************************************

#include <vector>

#include "Source.h"
#include "Quadrature.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

void
tk::srcIntP0( ncomp_t system,
              ncomp_t ncomp,
              ncomp_t offset,
              real t,
              const Fields& geoElem,
              const SrcFn& src,
              Fields& R )
// *****************************************************************************
//  Compute source term integrals for DG(P0)
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] t Physical time
//! \param[in] geoElem Element geometry array
//! \param[in] src Source function to use
//! \param[in,out] R Right-hand side vector added to
// *****************************************************************************
{
  for (std::size_t e=0; e<geoElem.nunk(); ++e) {
    auto vole = geoElem(e,0,0);
    auto xc = geoElem(e,1,0);
    auto yc = geoElem(e,2,0);
    auto zc = geoElem(e,3,0);
    auto s = src( system, ncomp, xc, yc, zc, t );
    for (ncomp_t c=0; c<ncomp; ++c)
      R(e, c, offset) += vole * s[c];
  }
}

void
tk::srcIntP1( ncomp_t system,
              ncomp_t ncomp,
              ncomp_t offset,
              real t,
              const std::vector< std::size_t >& inpoel,
              const UnsMesh::Coords& coord,
              const Fields& geoElem,
              const SrcFn& src,
              Fields& R )
// *****************************************************************************
//  Compute source term integrals for DG(P1)
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] t Physical time
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] src Source function to use
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();

  // Number of integration points
  constexpr std::size_t NG = 5;

  // arrays for quadrature points
  std::array< std::array< real, NG >, 3 > coordgp;
  std::array< real, NG > wgp;

  // get quadrature point weights and coordinates for tetrahedron
  GaussQuadratureTet( coordgp, wgp );

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  for (std::size_t e=0; e<geoElem.nunk(); ++e)
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

    for (std::size_t igp=0; igp<NG; ++igp)
    {
      auto shp1 = 1.0 - coordgp[0][igp] - coordgp[1][igp] - coordgp[2][igp];
      auto shp2 = coordgp[0][igp];
      auto shp3 = coordgp[1][igp];
      auto shp4 = coordgp[2][igp];

      auto wt = wgp[igp] * geoElem(e, 0, 0);

      auto xgp = x1*shp1 + x2*shp2 + x3*shp3 + x4*shp4;
      auto ygp = y1*shp1 + y2*shp2 + y3*shp3 + y4*shp4;
      auto zgp = z1*shp1 + z2*shp2 + z3*shp3 + z4*shp4;

      auto B2 = 2.0 * coordgp[0][igp] + coordgp[1][igp] + coordgp[2][igp] - 1.0;
      auto B3 = 3.0 * coordgp[1][igp] + coordgp[2][igp] - 1.0;
      auto B4 = 4.0 * coordgp[2][igp] - 1.0;

      auto s = src( system, ncomp, xgp, ygp, zgp, t );
      for (ncomp_t c=0; c<ncomp; ++c) {
        auto mark = c*ndof;
        R(e, mark, offset)   += wt * s[c];
        R(e, mark+1, offset) += wt * s[c] * B2;
        R(e, mark+2, offset) += wt * s[c] * B3;
        R(e, mark+3, offset) += wt * s[c] * B4;
      }
    }
  }
}

void
tk::srcIntP2( ncomp_t system,
              ncomp_t ncomp,
              ncomp_t offset,
              real t,
              const std::vector< std::size_t >& inpoel,
              const UnsMesh::Coords& coord,
              const Fields& geoElem,
              const SrcFn& src,
              Fields& R )
// *****************************************************************************
//  Compute source term integrals for DG(P2)
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] t Physical time
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] src Source function to use
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();

  // Number of integration points
  constexpr std::size_t NG = 11;

  // arrays for quadrature points
  std::array< std::array< real, NG >, 3 > coordgp;
  std::array< real, NG > wgp;

  // get quadrature point weights and coordinates for tetrahedron
  GaussQuadratureTet( coordgp, wgp );

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  for (std::size_t e=0; e<geoElem.nunk(); ++e)
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

    for (std::size_t igp=0; igp<NG; ++igp)
    {
      auto shp1 = 1.0 - coordgp[0][igp] - coordgp[1][igp] - coordgp[2][igp];
      auto shp2 = coordgp[0][igp];
      auto shp3 = coordgp[1][igp];
      auto shp4 = coordgp[2][igp];

      auto wt = wgp[igp] * geoElem(e, 0, 0);

      auto xgp = x1*shp1 + x2*shp2 + x3*shp3 + x4*shp4;
      auto ygp = y1*shp1 + y2*shp2 + y3*shp3 + y4*shp4;
      auto zgp = z1*shp1 + z2*shp2 + z3*shp3 + z4*shp4;

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
      auto B5 = 6.0 * xi_xi + eta_eta + zeta_zeta + 6.0 * xi_eta + 6.0 * xi_zeta + 2.0 * eta_zeta 
              - 6.0 * xi - 2.0 * eta - 2.0 * zeta + 1.0;
      auto B6 = 5.0 * eta_eta + zeta_zeta + 10.0 * xi_eta + 2.0 * xi_zeta + 6.0 * eta_zeta - 2.0 * xi 
              - 6.0 * eta - 2.0 * zeta + 1.0;
      auto B7 = 6.0 * zeta_zeta + 12.0 * xi_zeta + 6.0 * eta_zeta - 2.0 * xi 
              - eta - 7.0 * zeta + 1.0;
      auto B8 = 10.0 * eta_eta + zeta_zeta + 8.0 * eta_zeta - 8.0 * eta - 2.0 * zeta + 1.0;
      auto B9 = 6.0 * zeta_zeta + 18.0 * eta_zeta - 3.0 * eta - 7.0 * zeta + 1.0;
      auto B10 = 15.0 * zeta_zeta - 10.0 * zeta + 1.0;

      auto s = src( system, ncomp, xgp, ygp, zgp, t );
      for (ncomp_t c=0; c<ncomp; ++c) {
        auto mark = c*ndof;
        R(e, mark, offset)   += wt * s[c];
        R(e, mark+1, offset) += wt * s[c] * B2;
        R(e, mark+2, offset) += wt * s[c] * B3;
        R(e, mark+3, offset) += wt * s[c] * B4;
        R(e, mark+4, offset) += wt * s[c] * B5;
        R(e, mark+5, offset) += wt * s[c] * B6;
        R(e, mark+6, offset) += wt * s[c] * B7;
        R(e, mark+7, offset) += wt * s[c] * B8;
        R(e, mark+8, offset) += wt * s[c] * B9;
        R(e, mark+9, offset) += wt * s[c] * B10;
      }
    }
  }
}
