// *****************************************************************************
/*!
  \file      src/PDE/Limiter.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Limiters for DG
  \details   This file contains functions that provide limiter function
    calculations for maintaining monotonicity near solution discontinuities
    for the DG discretization.
*/
// *****************************************************************************

#include <array>
#include <vector>

#include "Limiter.h"

namespace inciter { extern ctr::InputDeck g_inputdeck; }

void
WENO_P1( const std::vector< int >& esuel,
  inciter::ncomp_t offset,
  const tk::Fields& U,
  tk::Fields& limFunc )
// *****************************************************************************
//  Weighted Essentially Non-Oscillatory (WENO) limiter for DGP1
//! \param[in] esuel Elements surrounding elements
//! \param[in] offset
//! \param[in] U High-order solution vector
//! \param[in,out] limFunc Limiter function
// *****************************************************************************
{
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  const auto cweight = inciter::g_inputdeck.get< tag::discr, tag::cweight >();

  std::size_t ncomp = U.nprop()/ndof;

  tk::real wtotal;
  std::array< std::array< tk::real, 5 >, 3 > gradu;
  std::array< tk::real, 5 > wtStencil, osc, wtDof;

  for (inciter::ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;

    for (std::size_t e=0; e<U.nunk(); ++e)
    {
      // reset all stencil values to zero
      for (std::size_t is=0; is<5; ++is)
      {
        gradu[0][is] = 0.0;
        gradu[1][is] = 0.0;
        gradu[2][is] = 0.0;
        osc[is] = 0.0;
        wtDof[is] = 0.0;
        wtStencil[is] = 0.0;
      }

      // The WENO limiter uses solution data from the neighborhood in the form
      // of stencils to enforce non-oscillatory conditions. The immediate
      // (Von Neumann) neighborhood of a tetrahedral cell consists of the 4
      // cells that share faces with it. These are the 4 neighborhood-stencils
      // for the tetrahedron. The primary stencil is the tet itself. Weights are
      // assigned to these stencils, with the primary stencil usually assigned
      // the highest weight. The lower the primary/central weight, the more
      // dissipative the limiting effect. This central weight is usually problem
      // dependent. It is set higher for relatively weaker discontinuities, and
      // lower for stronger discontinuities.

      // primary stencil
      gradu[0][0] = U(e, mark+1, offset);
      gradu[1][0] = U(e, mark+2, offset);
      gradu[2][0] = U(e, mark+3, offset);
      wtStencil[0] = cweight;

      // stencils from the neighborhood
      for (std::size_t is=1; is<5; ++is)
      {
        auto nel = esuel[ 4*e+(is-1) ];

        // ignore physical domain ghosts
        if (nel == -1)
        {
          gradu[0][is] = 0.0;
          gradu[1][is] = 0.0;
          gradu[2][is] = 0.0;
          wtStencil[is] = 0.0;
          continue;
        }

        gradu[0][is] = U(nel, mark+1, offset);
        gradu[1][is] = U(nel, mark+2, offset);
        gradu[2][is] = U(nel, mark+3, offset);
        wtStencil[is] = 1.0;
      }

      // From these stencils, an oscillation indicator is calculated, which
      // determines the effective weights for the high-order solution DOFs.
      // These effective weights determine the contribution of each of the
      // stencils to the high-order solution DOFs of the current cell which are
      // being limited. If this indicator detects a large oscillation in the
      // solution of the current cell, it reduces the effective weight for the
      // central stencil contribution to its high-order DOFs. This results in
      // a more dissipative and well-behaved solution in the troubled cell.

      // oscillation indicators
      for (std::size_t is=0; is<5; ++is)
      {
        osc[is] = sqrt(  gradu[0][is] * gradu[0][is]
                       + gradu[1][is] * gradu[1][is]
                       + gradu[2][is] * gradu[2][is] );
      }

      wtotal = 0.0;

      // effective weights for dofs
      for (std::size_t is=0; is<5; ++is)
      {
        wtDof[is] = wtStencil[is] * pow( (1.0e-8 + osc[is]), -2 );
        wtotal += wtDof[is];
      }

      for (std::size_t is=0; is<5; ++is)
      {
        wtDof[is] = wtDof[is]/wtotal;
      }

      limFunc(e, mark+1, 0) = 0.0;
      limFunc(e, mark+2, 0) = 0.0;
      limFunc(e, mark+3, 0) = 0.0;

      // limiter function
      for (std::size_t is=0; is<5; ++is)
      {
        limFunc(e, mark+1, 0) += wtDof[is]*gradu[0][is] / (gradu[0][0]+1.0e-12);
        limFunc(e, mark+2, 0) += wtDof[is]*gradu[1][is] / (gradu[1][0]+1.0e-12);
        limFunc(e, mark+3, 0) += wtDof[is]*gradu[2][is] / (gradu[2][0]+1.0e-12);
      }
    }
  }
}
