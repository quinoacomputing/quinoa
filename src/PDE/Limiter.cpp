// *****************************************************************************
/*!
  \file      src/PDE/Limiter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Limiters for discontiunous Galerkin methods
  \details   This file contains functions that provide limiter function
    calculations for maintaining monotonicity near solution discontinuities
    for the DG discretization.
*/
// *****************************************************************************

#include <array>
#include <vector>

#include "Vector.hpp"
#include "Limiter.hpp"
#include "DerivedData.hpp"
#include "Integrate/Quadrature.hpp"
#include "Integrate/Basis.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

void
WENO_P1( const std::vector< int >& esuel,
         inciter::ncomp_t offset,
         tk::Fields& U )
// *****************************************************************************
//  Weighted Essentially Non-Oscillatory (WENO) limiter for DGP1
//! \param[in] esuel Elements surrounding elements
//! \param[in] offset Index for equation systems
//! \param[in,out] U High-order solution vector which gets limited
//! \note This limiter function is experimental and untested. Use with caution.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();
  const auto cweight = inciter::g_inputdeck.get< tag::discr, tag::cweight >();
  std::array< std::vector< tk::real >, 3 > limU;
  limU[0].resize( U.nunk() );
  limU[1].resize( U.nunk() );
  limU[2].resize( U.nunk() );

  std::size_t ncomp = U.nprop()/rdof;

  std::array< std::array< tk::real, 3 >, 5 > gradu;
  std::array< tk::real, 5 > wtStencil, osc, wtDof;

  for (inciter::ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*rdof;

    for (std::size_t e=0; e<esuel.size()/4; ++e)
    {
      // reset all stencil values to zero
      for (auto& g : gradu) g.fill(0.0);
      osc.fill(0);
      wtDof.fill(0);
      wtStencil.fill(0);

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
      gradu[0][1] = U(e, mark+2, offset);
      gradu[0][2] = U(e, mark+3, offset);
      wtStencil[0] = cweight;

      // stencils from the neighborhood
      for (std::size_t is=1; is<5; ++is)
      {
        auto nel = esuel[ 4*e+(is-1) ];

        // ignore physical domain ghosts
        if (nel == -1)
        {
          gradu[is].fill(0.0);
          wtStencil[is] = 0.0;
          continue;
        }

        std::size_t n = static_cast< std::size_t >( nel );
        gradu[is][0] = U(n, mark+1, offset);
        gradu[is][1] = U(n, mark+2, offset);
        gradu[is][2] = U(n, mark+3, offset);
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
        osc[is] = std::sqrt( tk::dot(gradu[is], gradu[is]) );

      tk::real wtotal = 0;

      // effective weights for dofs
      for (std::size_t is=0; is<5; ++is)
      {
        // A small number (1.0e-8) is needed here to avoid dividing by a zero in
        // the case of a constant solution, where osc would be zero. The number
        // is not set to machine zero because it is squared, and a number
        // between 1.0e-8 to 1.0e-6 is needed.
        wtDof[is] = wtStencil[is] * pow( (1.0e-8 + osc[is]), -2 );
        wtotal += wtDof[is];
      }

      for (std::size_t is=0; is<5; ++is)
      {
        wtDof[is] = wtDof[is]/wtotal;
      }

      limU[0][e] = 0.0;
      limU[1][e] = 0.0;
      limU[2][e] = 0.0;

      // limiter function
      for (std::size_t is=0; is<5; ++is)
      {
        limU[0][e] += wtDof[is]*gradu[is][0];
        limU[1][e] += wtDof[is]*gradu[is][1];
        limU[2][e] += wtDof[is]*gradu[is][2];
      }
    }

    for (std::size_t e=0; e<esuel.size()/4; ++e)
    {
      U(e, mark+1, offset) = limU[0][e];
      U(e, mark+2, offset) = limU[1][e];
      U(e, mark+3, offset) = limU[2][e];
    }
  }
}

void
Superbee_P1( const std::vector< int >& esuel,
             const std::vector< std::size_t >& inpoel,
             const std::vector< std::size_t >& ndofel,
             inciter::ncomp_t offset,
             const tk::UnsMesh::Coords& coord,
             tk::Fields& U,
             std::size_t nmat )
// *****************************************************************************
//  Superbee limiter for DGP1
//! \param[in] esuel Elements surrounding elements
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] offset Index for equation systems
//! \param[in] coord Array of nodal coordinates
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in] nmat Number of materials in this PDE system. Default is 1, so
//!   this argument can be left unspecified for single-material by the calling
//!   code.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  std::size_t ncomp = U.nprop()/rdof;

  auto beta_lim = 2.0;

  for (std::size_t e=0; e<esuel.size()/4; ++e)
  {
    // If an rDG method is set up (P0P1), then, currently we compute the P1
    // basis functions and solutions by default. This implies that P0P1 is
    // unsupported in the p-adaptive DG (PDG). This is a workaround until we
    // have rdofel, which is needed to distinguish between ndofs and rdofs per
    // element for pDG.
    std::size_t dof_el;
    if (rdof > ndof)
    {
      dof_el = rdof;
    }
    else
    {
      dof_el = ndofel[e];
    }

    if (dof_el > 1)
    {
      // Superbee is a TVD limiter, which uses min-max bounds that the
      // high-order solution should satisfy, to ensure TVD properties. For a
      // high-order method like DG, this involves 3 steps:
      // 1. Find min-max bounds in the immediate neighborhood of cell.
      // 2. Calculate the Superbee function for all the points where solution
      //    needs to be reconstructed to (all quadrature points). From these,
      //    use the minimum value of the limiter function.
      // 3. Limit the high-order terms using this function.

      std::vector< tk::real > uMin(ncomp, 0.0), uMax(ncomp, 0.0);

      for (inciter::ncomp_t c=0; c<ncomp; ++c)
      {
        auto mark = c*rdof;
        uMin[c] = U(e, mark, offset);
        uMax[c] = U(e, mark, offset);
      }

      // ----- Step-1: find min/max in the neighborhood
      for (std::size_t is=0; is<4; ++is)
      {
        auto nel = esuel[ 4*e+is ];

        // ignore physical domain ghosts
        if (nel == -1) continue;

        auto n = static_cast< std::size_t >( nel );
        for (inciter::ncomp_t c=0; c<ncomp; ++c)
        {
          auto mark = c*rdof;
          uMin[c] = std::min(uMin[c], U(n, mark, offset));
          uMax[c] = std::max(uMax[c], U(n, mark, offset));
        }
      }

      // ----- Step-2: loop over all quadrature points to get limiter function

      // to loop over all the quadrature points of all faces of element e,
      // coordinates of the quadrature points are needed.
      // Number of quadrature points for face integration
      auto ng = tk::NGfa(ndof);

      // arrays for quadrature points
      std::array< std::vector< tk::real >, 2 > coordgp;
      std::vector< tk::real > wgp;

      coordgp[0].resize( ng );
      coordgp[1].resize( ng );
      wgp.resize( ng );

      // get quadrature point weights and coordinates for triangle
      tk::GaussQuadratureTri( ng, coordgp, wgp );

      const auto& cx = coord[0];
      const auto& cy = coord[1];
      const auto& cz = coord[2];

      // Extract the element coordinates
      std::array< std::array< tk::real, 3>, 4 > coordel {{
        {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
        {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
        {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
        {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }} }};

      // Compute the determinant of Jacobian matrix
      auto detT =
        tk::Jacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

      // initialize limiter function
      std::vector< tk::real > phi(ncomp, 1.0);
      for (std::size_t lf=0; lf<4; ++lf)
      {
        // Extract the face coordinates
        std::array< std::size_t, 3 > inpofa_l {{ inpoel[4*e+tk::lpofa[lf][0]],
                                                 inpoel[4*e+tk::lpofa[lf][1]],
                                                 inpoel[4*e+tk::lpofa[lf][2]] }};

        std::array< std::array< tk::real, 3>, 3 > coordfa {{
          {{ cx[ inpofa_l[0] ], cy[ inpofa_l[0] ], cz[ inpofa_l[0] ] }},
          {{ cx[ inpofa_l[1] ], cy[ inpofa_l[1] ], cz[ inpofa_l[1] ] }},
          {{ cx[ inpofa_l[2] ], cy[ inpofa_l[2] ], cz[ inpofa_l[2] ] }} }};

        // Gaussian quadrature
        for (std::size_t igp=0; igp<ng; ++igp)
        {
          // Compute the coordinates of quadrature point at physical domain
          auto gp = tk::eval_gp( igp, coordfa, coordgp );

          //Compute the basis functions
          auto B_l = tk::eval_basis( rdof,
                tk::Jacobian( coordel[0], gp, coordel[2], coordel[3] ) / detT,
                tk::Jacobian( coordel[0], coordel[1], gp, coordel[3] ) / detT,
                tk::Jacobian( coordel[0], coordel[1], coordel[2], gp ) / detT );

          auto state = tk::eval_state( ncomp, offset, rdof, dof_el, e, U, B_l );

          Assert( state.size() == ncomp, "Size mismatch" );

          // compute the limiter function
          for (inciter::ncomp_t c=0; c<ncomp; ++c)
          {
            auto phi_gp = 1.0;
            auto mark = c*rdof;
            auto uNeg = state[c] - U(e, mark, offset);
            if (uNeg > 1.0e-14)
            {
              phi_gp = std::min( 1.0, (uMax[c]-U(e, mark, offset))/(2.0*uNeg) );
            }
            else if (uNeg < -1.0e-14)
            {
              phi_gp = std::min( 1.0, (uMin[c]-U(e, mark, offset))/(2.0*uNeg) );
            }
            else
            {
              phi_gp = 1.0;
            }
            phi_gp = std::max( 0.0,
                               std::max( std::min(beta_lim*phi_gp, 1.0),
                                         std::min(phi_gp, beta_lim) ) );
            phi[c] = std::min( phi[c], phi_gp );
          }
        }
      }

      if (nmat > 1)
        consistentMultiMatLimiting_P1(nmat, offset, rdof, e, U, phi);

      // ----- Step-3: apply limiter function
      for (inciter::ncomp_t c=0; c<ncomp; ++c)
      {
        auto mark = c*rdof;
        U(e, mark+1, offset) = phi[c] * U(e, mark+1, offset);
        U(e, mark+2, offset) = phi[c] * U(e, mark+2, offset);
        U(e, mark+3, offset) = phi[c] * U(e, mark+3, offset);
      }
    }
  }
}

void consistentMultiMatLimiting_P1( std::size_t nmat,
                                    ncomp_t offset,
                                    std::size_t rdof,
                                    std::size_t e,
                                    tk::Fields& U,
                                    std::vector< tk::real >& phi )
// *****************************************************************************
//  Consistent limiter modifications for P1 dofs
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] offset Index for equation system
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] e Element being checked for consistency
//! \param[in,out] U Second-order solution vector which gets modified near
//!   material interfaces for consistency
//! \param[in,out] phi Vector of limiter functions for the ncomp unknowns
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::densityIdx;
  using inciter::momentumIdx;
  using inciter::energyIdx;
  using inciter::volfracDofIdx;

  // find the limiter-function for volume-fractions
  auto phi_al(1.0), almax(0.0), dalmax(0.0);
  //std::size_t nmax(0);
  for (std::size_t k=0; k<nmat; ++k)
  {
    phi_al = std::min( phi_al, phi[volfracIdx(nmat, k)] );
    if (almax < U(e,volfracDofIdx(nmat, k, rdof, 0),offset))
    {
      //nmax = k;
      almax = U(e,volfracDofIdx(nmat, k, rdof, 0),offset);
    }
    auto dmax = std::max(
                  std::max(
                    std::abs(U(e,volfracDofIdx(nmat, k, rdof, 1),offset)),
                    std::abs(U(e,volfracDofIdx(nmat, k, rdof, 2),offset)) ),
                  std::abs(U(e,volfracDofIdx(nmat, k, rdof, 3),offset)) );
    dalmax = std::max( dalmax, dmax );
  }

  //phi_al = phi[nmax];

  // determine if cell is a material-interface cell based on ad-hoc tolerances.
  // if interface-cell, then modify high-order dofs of conserved unknowns
  // consistently and use same limiter for all equations
  if ( dalmax > 0.01 ||
       (almax > 0.0001 && almax < (1.0-0.0001)) )
  {
    // 1. consistent high-order dofs
    std::array< tk::real, 3 > drhob {{ 0.0, 0.0, 0.0 }};
    auto rhob(0.0), vel(0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      auto alk = std::max( 1.0e-14, U(e,volfracDofIdx(nmat, k, rdof, 0),offset) );
      auto rhok = U(e,densityDofIdx(nmat, k, rdof, 0),offset)/alk;
      auto rhoEk = U(e,energyDofIdx(nmat, k, rdof, 0),offset)/alk;
      for (std::size_t idir=1; idir<=3; ++idir)
      {
        U(e,densityDofIdx(nmat, k, rdof, idir),offset) = rhok *
          U(e,volfracDofIdx(nmat, k, rdof, idir),offset);
        U(e,energyDofIdx(nmat, k, rdof, idir),offset) = rhoEk *
          U(e,volfracDofIdx(nmat, k, rdof, idir),offset);
        drhob[idir-1] += U(e,densityDofIdx(nmat, k, rdof, idir),offset);
      }
      rhob += U(e,densityDofIdx(nmat, k, rdof, 0),offset);
    }
    for (std::size_t idir=1; idir<=3; ++idir)
    {
      for (std::size_t jdir=1; jdir<=3; ++jdir)
      {
        vel = U(e,momentumDofIdx(nmat, jdir-1, rdof, 0),offset)/rhob;
        U(e,momentumDofIdx(nmat, jdir-1, rdof, idir),offset) = vel * drhob[idir-1];
      }
    }

    // 2. same limiter for all equations
    for (auto& p : phi) p = phi_al;
  }
  else
  {
    // same limiter for all volume-fractions
    for (std::size_t k=0; k<nmat; ++k)
      phi[volfracIdx(nmat, k)] = phi_al;
  }
}

} // inciter::
