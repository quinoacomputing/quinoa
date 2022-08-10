// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Initialize.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for initialization of system of PDEs in DG methods
  \details   This file contains functionality for setting initial conditions
     and evaluating known solutions used in discontinuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************

#include <array>
#include <vector>
#include <iostream>
#include "Data.hpp"
#include "Initialize.hpp"
#include "Quadrature.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "EoS/EoS_Base.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

void
tk::initialize( ncomp_t system,
                ncomp_t ncomp,
                ncomp_t offset,
                const std::vector< inciter::EoS_Base* >& mat_blk,
                const Fields& L,
                const std::vector< std::size_t >& inpoel,
                const UnsMesh::Coords& coord,
                const InitializeFn& solution,
                Fields& unk,
                real t,
                const std::size_t nielem )
// *****************************************************************************
//! Initalize a system of DGPDEs by projecting the exact solution in the DG
//! solution space
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
//! \param[in] nielem Number of internal elements
// *****************************************************************************
{
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();

  // Number of quadrature points for volume integration
  auto ng = tk::NGinit(ndof);

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

  for (std::size_t e=0; e<nielem; ++e) {    // for all tets
    // The volume of tetrahedron
    auto vole = L(e, 0, offset);

    // Extract the element coordinates
    std::array< std::array< real, 3>, 4 > coordel {{
      {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
      {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
      {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
      {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }} }};

    // right hand side vector
    std::vector< real > R( ncomp*ndof, 0.0 );

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordel, coordgp );

      // Compute the basis function
      auto B =
        eval_basis( ndof, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp] );

      const auto s = solution( system, ncomp, mat_blk, gp[0], gp[1], gp[2], t );

      auto wt = wgp[igp] * vole;

      update_rhs( ncomp, ndof, wt, B, s, R );
    }

    // Compute the initial conditions
    eval_init(ncomp, offset, ndof, rdof, e, R, L, unk);
  }
}

void
tk::update_rhs( ncomp_t ncomp,
                const std::size_t ndof,
                const tk::real wt,
                const std::vector< tk::real >& B,
                const std::vector< tk::real >& s,
                std::vector< tk::real >& R )
// *****************************************************************************
//  Update the rhs by adding the initial analytical solution term
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] ndof Number of degrees of freedom
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] B Vector of basis functions
//! \param[in] s Vector of analytical solution at quadrature point
//! \param[in,out] R Right-hand side vector
// *****************************************************************************
{
  Assert( B.size() == ndof, "Size mismatch for basis function" );
  Assert( s.size() >= ncomp, "Size mismatch for source term" );

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    // DG(P0)
    auto mark = c*ndof;
    R[mark] += wt * s[c];

    if(ndof > 1)         //DG(P1)
    {
      R[mark+1] += wt * s[c] * B[1];
      R[mark+2] += wt * s[c] * B[2];
      R[mark+3] += wt * s[c] * B[3];

      if(ndof > 4)      //DG(P2)
      {
        R[mark+4] += wt * s[c] * B[4];
        R[mark+5] += wt * s[c] * B[5];
        R[mark+6] += wt * s[c] * B[6];
        R[mark+7] += wt * s[c] * B[7];
        R[mark+8] += wt * s[c] * B[8];
        R[mark+9] += wt * s[c] * B[9];
      }
    }
  }
}

void
tk::eval_init( ncomp_t ncomp,
               ncomp_t offset,
               const std::size_t ndof,
               const std::size_t rdof,
               const std::size_t e,
               const std::vector< tk::real >& R,
               const Fields& L,
               Fields& unk )
// *****************************************************************************
//  Compute the initial conditions
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] ndof Number of degrees of freedom
//! \param[in] rdof Total number of reconstructed degrees of freedom
//! \param[in] e Element index
//! \param[in] R Right-hand side vector
//! \param[in] L Block diagonal mass matrix
//! \param[in,out] unk Array of unknowns
// *****************************************************************************
{
  for (ncomp_t c=0; c<ncomp; ++c)
  {
    // DG(P0)
    auto mark = c*ndof;
    auto rmark = c*rdof;
    unk(e, rmark, offset) = R[mark] / L(e, mark,   offset);

    // if P0P1, initialize higher dofs to 0
    if (rdof > ndof)
    {
      unk(e, rmark+1, offset) = 0.0;
      unk(e, rmark+2, offset) = 0.0;
      unk(e, rmark+3, offset) = 0.0;
    }

    if(ndof > 1)          // DG(P1)
    {
      unk(e, rmark+1, offset) = R[mark+1] / L(e, mark+1, offset);
      unk(e, rmark+2, offset) = R[mark+2] / L(e, mark+2, offset);
      unk(e, rmark+3, offset) = R[mark+3] / L(e, mark+3, offset);
 
      if(ndof > 4)        // DG(P2)
      {
        unk(e, rmark+4, offset) = R[mark+4] / L(e, mark+4, offset);
        unk(e, rmark+5, offset) = R[mark+5] / L(e, mark+5, offset);
        unk(e, rmark+6, offset) = R[mark+6] / L(e, mark+6, offset);
        unk(e, rmark+7, offset) = R[mark+7] / L(e, mark+7, offset);
        unk(e, rmark+8, offset) = R[mark+8] / L(e, mark+8, offset);
        unk(e, rmark+9, offset) = R[mark+9] / L(e, mark+9, offset);
      }
    }
  }
}
