// *****************************************************************************
/*!
  \file      src/PDE/Integrate/SolidTerms.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing integrals for the RHS terms
             of the inverse deformation equations in the multi-material
             solid dynamics solver
  \details   This file contains routines to integrate right hand side terms
             to the inverse deformation equations for the multi-material
             solid dynamics solver.
*/
// *****************************************************************************

#include <array>
#include <vector>

#include "SolidTerms.hpp"
#include "Vector.hpp"
#include "Quadrature.hpp"
#include "Reconstruction.hpp"
#include "MultiMatTerms.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "MultiMat/MiscMultiMatFns.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {
extern ctr::InputDeck g_inputdeck;
}

namespace tk {

void
solidTermsVolInt(
               std::size_t nmat,
               const std::vector< inciter::EOS >& mat_blk,
               const std::size_t ndof,
               const std::size_t rdof,
               const std::size_t nelem,
               const std::vector< std::size_t >& inpoel,
               const UnsMesh::Coords& coord,
               const Fields& geoElem,
               const Fields& U,
               const Fields& P,
               const std::vector< std::size_t >& ndofel,
               const std::vector< tk::real >& rho0mat,
               const tk::real dt,
               Fields& R,
               int intsharp )
// *****************************************************************************
//  Compute all RHS volume terms in the inverse deformation equations to
//  satisfy the condition curl(g) = 0.
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] nelem Maximum number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] rho0mat Initial densities of all materials
//! \param[in] dt Delta time
//! \param[in,out] R Right-hand side vector computed
//! \param[in] intsharp Interface compression tag, an optional argument, with
//!   default 0, so that it is unused for single-material and transport.
// *****************************************************************************
{

  using inciter::velocityDofIdx;
  using inciter::volfracDofIdx;
  using inciter::deformDofIdx;

  const auto& solidx =
    inciter::g_inputdeck.get< tag::matidxmap, tag::solidx >();

  auto ncomp = R.nprop()/ndof;
  auto nprim = P.nprop()/rdof;
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  for (std::size_t e=0; e<nelem; ++e)
  {

    // Relaxation coefficient
    tk::real eta = 1.0/(6.0*dt);

    auto ng = tk::NGvol(ndofel[e]);

    // arrays for quadrature points
    std::array< std::vector< real >, 3 > coordgp;
    std::vector< real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    coordgp[2].resize( ng );
    wgp.resize( ng );

    GaussQuadratureTet( ng, coordgp, wgp );

    // Extract the element coordinates
    std::array< std::array< real, 3>, 4 > coordel {{
      {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
      {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
      {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
      {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }}}};

    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordel, coordgp );

      // Compute the basis function
      auto B =
        eval_basis(rdof,coordgp[0][igp],coordgp[1][igp],coordgp[2][igp]);

      // Get state
      std::vector< real > state;
      state = evalPolynomialSol(mat_blk, intsharp, ncomp, nprim, rdof,
        nmat, e, rdof, inpoel, coord, geoElem, gp, B, U, P);

      // Loop through materials
      for (std::size_t k=0; k<nmat; ++k)
      {
        if (solidx[k] > 0)
        {
          tk::real alpha = state[inciter::volfracIdx(nmat, k)];
          std::array< std::array< tk::real, 3 >, 3 > g;
          // Compute the source terms
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              g[i][j] = state[inciter::deformIdx(nmat,solidx[k],i,j)];

          // Compute rhs factor
          tk::real rho = state[inciter::densityIdx(nmat, k)]/alpha;
          tk::real rho0 = rho0mat[k];
          tk::real rfact = eta*(rho/(rho0*tk::determinant(g))-1.0);

          // Compute the source terms
          std::vector< real > s(9*ndof, 0.0);
          std::vector< real > deriv(6, 0.0);
          std::vector< std::size_t > defIdx(6, 0);
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              for (std::size_t idof=0; idof<ndof; ++idof)
                s[(i*3+j)*ndof+idof] = B[idof]*rfact*alpha*g[i][j];

          auto wt = wgp[igp] * geoElem(e, 0);

          update_rhs( nmat, ndof, wt, e, solidx[k], s, R );

        }
      }
     }
  }

}

void
update_rhs( std::size_t nmat,
            const std::size_t ndof,
            const tk::real wt,
            const std::size_t e,
            const std::size_t solidx,
            const std::vector< real >& s,
            Fields& R )
// *****************************************************************************
//  Update the rhs by adding the source term integrals
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] rdof Maximum number of degrees of freedom
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] e Element index
//! \param[in] solidx Material index indicator
//! \param[in] s Source term vector
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      for (std::size_t idof=0; idof<ndof; ++idof)
      {
        std::size_t dofId = inciter::deformDofIdx(nmat,solidx,i,j,ndof,idof);
        std::size_t srcId = (i*3+j)*ndof+idof;
        R(e, dofId) += wt * s[srcId];
      }
}

} // tk::
