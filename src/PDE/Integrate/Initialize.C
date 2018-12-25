// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Initialize.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for initialization of system of PDEs in DG methods
  \details   This file contains functionality for setting initial conditions
     and evaluating known solutions used in discontinuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************

#include <array>

#include "Initialize.h"
#include "Quadrature.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

void
tk::initializeP0( ncomp_t system,
                  ncomp_t ncomp,
                  ncomp_t offset,
                  const std::vector< std::size_t >& inpoel,
                  const UnsMesh::Coords& coord,
                  const SolutionFn& solution,
                  Fields& unk,
                  real t )
// *****************************************************************************
//  Initalize a PDE system for DG(P0)
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] solution Function to call to evaluate known solution or initial
//!   conditions at x,y,z,t
//! \param[in,out] unk Array of unknowns
//! \param[in] t Physical time
// *****************************************************************************
{
  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t e=0; e<unk.nunk(); ++e) {    // for all tets
    // node ids
    const auto A = inpoel[e*4+0];
    const auto B = inpoel[e*4+1];
    const auto C = inpoel[e*4+2];
    const auto D = inpoel[e*4+3];
    // compute centroid
    auto xcc = (x[A]+x[B]+x[C]+x[D])/4.0;
    auto ycc = (y[A]+y[B]+y[C]+y[D])/4.0;
    auto zcc = (z[A]+z[B]+z[C]+z[D])/4.0;
    // evaluate solution at centroid
    const auto s = solution( system, ncomp, xcc, ycc, zcc, t );
    // initialize unknown vector with solution at centroids
    for (ncomp_t c=0; c<ncomp; ++c) unk(e, c, offset) = s[c];
  }
}

void
tk::initializeP1( ncomp_t system,
                  ncomp_t ncomp,
                  ncomp_t offset,
                  const Fields& L,
                  const std::vector< std::size_t >& inpoel,
                  const UnsMesh::Coords& coord,
                  const SolutionFn& solution,
                  Fields& unk,
                  real t )
// *****************************************************************************
//  Initalize a PDE system for DG(P1)
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] L Block diagonal mass matrix
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] solution Function to call to evaluate known solution or initial
//!   conditions at x,y,z,t
//! \param[in,out] unk Array of unknowns
//! \param[in] t Physical time
// *****************************************************************************
{
  Assert( L.nunk() == unk.nunk(), "Size mismatch" );

  // Number of integration points
  constexpr std::size_t NG = 5;

  // Number of solution degrees of freedom
  constexpr std::size_t ndof = 4;

  // arrays for quadrature points
  std::array< std::array< real, NG >, 3 > coordgp;
  std::array< real, NG > wgp;

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  // get quadrature point weights and coordinates for tetrahedron
  GaussQuadratureTet( coordgp, wgp );

  for (std::size_t e=0; e<unk.nunk(); ++e) {    // for all tets
    auto vole = L(e, 0, offset);

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

    // right hand side vector
    std::vector< real > R( unk.nprop(), 0.0 );

    // Gaussian quadrature
    for (std::size_t igp=0; igp<NG; ++igp)
    {
      auto B2 = 2.0 * coordgp[0][igp] + coordgp[1][igp] + coordgp[2][igp] - 1.0;
      auto B3 = 3.0 * coordgp[1][igp] + coordgp[2][igp] - 1.0;
      auto B4 = 4.0 * coordgp[2][igp] - 1.0;

      auto shp1 = 1.0 - coordgp[0][igp] - coordgp[1][igp] - coordgp[2][igp];
      auto shp2 = coordgp[0][igp];
      auto shp3 = coordgp[1][igp];
      auto shp4 = coordgp[2][igp];

      auto xgp = x1*shp1 + x2*shp2 + x3*shp3 + x4*shp4;
      auto ygp = y1*shp1 + y2*shp2 + y3*shp3 + y4*shp4;
      auto zgp = z1*shp1 + z2*shp2 + z3*shp3 + z4*shp4;

      auto wt = vole * wgp[igp];

      const auto s = solution( system, ncomp, xgp, ygp, zgp, t );

      for (ncomp_t c=0; c<ncomp; ++c) {
        auto mark = c*ndof;
        R[mark  ] += wt * s[c];
        R[mark+1] += wt * s[c]*B2;
        R[mark+2] += wt * s[c]*B3;
        R[mark+3] += wt * s[c]*B4;
      }
    }

    for (ncomp_t c=0; c<ncomp; ++c) {
      auto mark = c*ndof;
      unk(e, mark,   offset) = R[mark]   / L(e, mark,   offset);
      unk(e, mark+1, offset) = R[mark+1] / L(e, mark+1, offset);
      unk(e, mark+2, offset) = R[mark+2] / L(e, mark+2, offset);
      unk(e, mark+3, offset) = R[mark+3] / L(e, mark+3, offset);
    }
  }
}

void
tk::initializeP2( ncomp_t system,
                  ncomp_t ncomp,
                  ncomp_t offset,
                  const Fields& L,
                  const std::vector< std::size_t >& inpoel,
                  const UnsMesh::Coords& coord,
                  const SolutionFn& solution,
                  Fields& unk,
                  real t )
// *****************************************************************************
//  Initalize a PDE system for DG(P1)
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] L Block diagonal mass matrix
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] solution Function to call to evaluate known solution or initial
//!   conditions at x,y,z,t
//! \param[in,out] unk Array of unknowns
//! \param[in] t Physical time
// *****************************************************************************
{
  Assert( L.nunk() == unk.nunk(), "Size mismatch" );

  // Number of integration points
  constexpr std::size_t NG = 14;

  // Number of solution degrees of freedom
  constexpr std::size_t ndof = 10;

  // arrays for quadrature points
  std::array< std::array< real, NG >, 3 > coordgp;
  std::array< real, NG > wgp;

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  // get quadrature point weights and coordinates for tetrahedron
  GaussQuadratureTet( coordgp, wgp );

  for (std::size_t e=0; e<unk.nunk(); ++e) {    // for all tets
    auto vole = L(e, 0, offset);

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

    // right hand side vector
    std::vector< real > R( unk.nprop(), 0.0 );

    // Gaussian quadrature
    for (std::size_t igp=0; igp<NG; ++igp)
    {
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
      auto B9 = 6.0 * zeta_zeta + 18.0 * eta_zeta - 3.0 * eta - 7.0 * zeta
              + 1.0;
      auto B10 = 15.0 * zeta_zeta - 10.0 * zeta + 1.0;

      auto shp1 = 1.0 - coordgp[0][igp] - coordgp[1][igp] - coordgp[2][igp];
      auto shp2 = coordgp[0][igp];
      auto shp3 = coordgp[1][igp];
      auto shp4 = coordgp[2][igp];

      auto xgp = x1*shp1 + x2*shp2 + x3*shp3 + x4*shp4;
      auto ygp = y1*shp1 + y2*shp2 + y3*shp3 + y4*shp4;
      auto zgp = z1*shp1 + z2*shp2 + z3*shp3 + z4*shp4;

      auto wt = vole * wgp[igp];

      const auto s = solution( system, ncomp, xgp, ygp, zgp, t );

      for (ncomp_t c=0; c<ncomp; ++c) {
        auto mark = c*ndof;
        R[mark  ] += wt * s[c];
        R[mark+1] += wt * s[c]*B2;
        R[mark+2] += wt * s[c]*B3;
        R[mark+3] += wt * s[c]*B4;
        R[mark+4] += wt * s[c]*B5;
        R[mark+5] += wt * s[c]*B6;
        R[mark+6] += wt * s[c]*B7;
        R[mark+7] += wt * s[c]*B8;
        R[mark+8] += wt * s[c]*B9;
        R[mark+9] += wt * s[c]*B10;
      }
    }

    for (ncomp_t c=0; c<ncomp; ++c) {
      auto mark = c*ndof;
      unk(e, mark,   offset) = R[mark]   / L(e, mark,   offset);
      unk(e, mark+1, offset) = R[mark+1] / L(e, mark+1, offset);
      unk(e, mark+2, offset) = R[mark+2] / L(e, mark+2, offset);
      unk(e, mark+3, offset) = R[mark+3] / L(e, mark+3, offset);
      unk(e, mark+4, offset) = R[mark+4] / L(e, mark+4, offset);
      unk(e, mark+5, offset) = R[mark+5] / L(e, mark+5, offset);
      unk(e, mark+6, offset) = R[mark+6] / L(e, mark+6, offset);
      unk(e, mark+7, offset) = R[mark+7] / L(e, mark+7, offset);
      unk(e, mark+8, offset) = R[mark+8] / L(e, mark+8, offset);
      unk(e, mark+9, offset) = R[mark+9] / L(e, mark+9, offset);
    }
  }
}

void
tk::initialize( ncomp_t system,
                ncomp_t ncomp,
                ncomp_t offset,
                const Fields& L,
                const std::vector< std::size_t >& inpoel,
                const UnsMesh::Coords& coord,
                const SolutionFn& solution,
                Fields& unk,
                real t )
// *****************************************************************************
//! Initalize a system of DGPDEs
//! \details This is the public interface exposed to client code.
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] L Block diagonal mass matrix
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of node coordinates
//! \param[in] solution Function to call to evaluate known solution or initial
//!   conditions at x,y,z,t
//! \param[in,out] unk Array of unknowns
//! \param[in] t Physical time
// *****************************************************************************
{
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  switch(ndof) 
  {
    case 1:
      initializeP0( system, ncomp, offset, inpoel, coord, solution, unk, t );
      break;
    case 4:
      initializeP1( system, ncomp, offset, L, inpoel, coord, solution, unk, t );
      break;
    case 10:
      initializeP2( system, ncomp, offset, L, inpoel, coord, solution, unk, t );
      break;
    default:
      Throw( "initialize() not defined" );
  }
}
