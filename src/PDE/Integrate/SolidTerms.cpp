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
solidTermsVolInt( ncomp_t system,
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
               const tk::real dt,
               Fields& R,
               int intsharp )
// *****************************************************************************
//  Compute all RHS volume terms in the inverse deformation equations to
//  satisfy the condition curl(g) = 0.
//! \param[in] system Equation system index
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
//! \param[in] dt Delta time
//! \param[in,out] R Right-hand side vector computed
//! \param[in] intsharp Interface compression tag, an optional argument, with
//!   default 0, so that it is unused for single-material and transport.
// *****************************************************************************
{

  using inciter::velocityDofIdx;
  using inciter::volfracDofIdx;
  using inciter::deformDofIdx;

  const auto& solidx = inciter::g_inputdeck.get< tag::param, tag::multimat,
    tag::matidxmap >().template get< tag::solidx >();

  auto ncomp = R.nprop()/ndof;
  auto nprim = P.nprop()/rdof;
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  for (std::size_t e=0; e<nelem; ++e)
  {

    // Relaxation and diffusion coefficients
    tk::real dx = geoElem(0,4);
    tk::real D = dx*dx/(12.0*dt) * 1.0e-00;
    tk::real eta = 1.0/(6.0*dt) * 1.0e+00;

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

    // Compute the inverse of the Jacobian
    auto jacInv =
      tk::inverseJacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordel, coordgp );

      // Compute the basis function
      auto B =
        eval_basis(rdof,coordgp[0][igp],coordgp[1][igp],coordgp[2][igp]);

      // Compute derivatives
      std::array< std::vector<tk::real>, 3 > dBdx;
      dBdx[0].resize( ndof, 0 );
      dBdx[1].resize( ndof, 0 );
      dBdx[2].resize( ndof, 0 );
      if (rdof > 1)
      {
        dBdx = tk::eval_dBdx_p1( rdof, jacInv );
        if(ndof > 4) {
          tk::eval_dBdx_p2(igp, coordgp, jacInv, dBdx);
        }
      }

      // Get state
      std::vector< real > state;
      state = evalPolynomialSol(system, mat_blk, intsharp, ncomp, nprim, rdof,
        nmat, e, rdof, inpoel, coord, geoElem, gp, B, U, P);

      // Get velocity
      std::vector< real > v = {{state[ncomp+inciter::velocityIdx(nmat, 0)],
                                state[ncomp+inciter::velocityIdx(nmat, 1)],
                                state[ncomp+inciter::velocityIdx(nmat, 2)]}};

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
              g[i][j] = state[inciter::deformIdx(nmat,solidx[k],i,j)]/alpha;

          // ** HARDCODED **
          // Should be able to store rho_0 for every cell at every Gauss point.
          tk::real rho = state[inciter::densityIdx(nmat, k)];
          tk::real rho0;
          if (k==0) rho0 = 8900.0;
          else rho0 = 2700.0;
          tk::real rfact = 0.0; //eta*(rho/(rho0*tk::determinant(g))-1.0);

          // Compute the source terms
          std::vector< real > s(9*ndof, 0.0);
          std::vector< real > deriv(6, 0.0);
          std::vector< std::size_t > defIdx(6, 0);
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              for (std::size_t idof=0; idof<ndof; ++idof)
              {
                for (std::size_t l=0; l<deriv.size(); ++l)
		  deriv[l] = 0.0;
                for (std::size_t jdof=0; jdof<rdof; ++jdof)
                {
                  // Find indeces for all unknowns used
                  defIdx[0]  = deformDofIdx(nmat,solidx[k],i,(j+1)%3,rdof,jdof);
                  defIdx[1]  = deformDofIdx(nmat,solidx[k],i,      j,rdof,jdof);
                  defIdx[2]  = deformDofIdx(nmat,solidx[k],i,      j,rdof,jdof);
                  defIdx[3]  = deformDofIdx(nmat,solidx[k],i,(j+2)%3,rdof,jdof);
                  defIdx[4]  = volfracDofIdx(nmat,k,rdof,jdof);
                  defIdx[5]  = volfracDofIdx(nmat,k,rdof,jdof);
                  // Compute derivatives
                  deriv[0]  += U(e,defIdx[0]) *dBdx[      j][jdof];
                  deriv[1]  += U(e,defIdx[1]) *dBdx[(j+1)%3][jdof];
                  deriv[2]  += U(e,defIdx[2]) *dBdx[(j+2)%3][jdof];
                  deriv[3]  += U(e,defIdx[3]) *dBdx[      j][jdof];
                  deriv[4]  += U(e,defIdx[4]) *dBdx[(j+1)%3][jdof];
                  deriv[5]  += U(e,defIdx[5]) *dBdx[(j+2)%3][jdof];
                }
                s[(i*3+j)*ndof+idof] = B[idof]*rfact*alpha*g[i][j]
                  + alpha*B[idof]*(v[(j+1)%3]*(deriv[0]-deriv[1])
                                  -v[(j+2)%3]*(deriv[2]-deriv[3]))
                  + D*((alpha*dBdx[(j+1)%3][idof]+B[idof]*deriv[4])
                       *(deriv[0]-deriv[1])
		      -(alpha*dBdx[(j+2)%3][idof]+B[idof]*deriv[5])
                       *(deriv[2]-deriv[3]));
              }

        auto wt = wgp[igp] * geoElem(e, 0);

        update_rhs( nmat, ndof, wt, e, solidx[k], s, R );

        }
      }
     }
  }

}

void
solidTermsSurfInt( std::size_t nmat,
                   const std::size_t ndof,
                   const std::size_t rdof,
                   const std::array< tk::real, 3 >& fn,
                   const std::size_t el,
                   const std::size_t er,
                   const std::vector< std::size_t >& solidx,
                   const Fields& geoElem,
                   const Fields& U,
                   const std::array< std::array< tk::real, 3>, 4 > coordel_l,
                   const std::array< std::array< tk::real, 3>, 4 > coordel_r,
                   const std::size_t igp,
                   const std::array< std::vector< tk::real >, 2 >& coordgp,
                   const tk::real dt,
                   std::vector< tk::real >& fl )
// *****************************************************************************
//  Compute all RHS surface terms in the inverse deformation equations to
//  satisfy the condition curl(g) = 0.
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] fn Face/Surface normal
//! \param[in] el Left element index
//! \param[in] er Right element index
//! \param[in] B_l Basis function for the left element
//! \param[in] B_r Basis function for the right element
//! \param[in] solidx Solid material indicator
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] coordel_l Coordinates of left elements' nodes
//! \param[in] coordel_r Coordinates of right elements' nodes
//! \param[in] igp Index of quadrature points
//! \param[in] coordgp Gauss point coordinates for tetrahedron element
//! \param[in] dt Delta time
//! \param[in,out] fl Surface flux
// *****************************************************************************
{

  using inciter::deformIdx;
  using inciter::volfracDofIdx;
  using inciter::deformDofIdx;

  // Diffusion coefficient
  //tk::real dx = std::min(geoElem(el,4),geoElem(er,4));
  tk::real dx = geoElem(0,4);
  tk::real D = dx*dx/(12.0*dt) * 1.0e-00;
  // Compute the inverse of the Jacobian
  auto jacInv_l =
    tk::inverseJacobian(coordel_l[0],coordel_l[1],coordel_l[2],coordel_l[3]);
  auto jacInv_r =
    tk::inverseJacobian(coordel_r[0],coordel_r[1],coordel_r[2],coordel_r[3]);
  // Compute derivatives
  std::array< std::vector<tk::real>, 3 > dBdx_l, dBdx_r;
  dBdx_l[0].resize( ndof, 0 );
  dBdx_l[1].resize( ndof, 0 );
  dBdx_l[2].resize( ndof, 0 );
  dBdx_r[0].resize( ndof, 0 );
  dBdx_r[1].resize( ndof, 0 );
  dBdx_r[2].resize( ndof, 0 );
  if (rdof > 1)
  {
    dBdx_l = tk::eval_dBdx_p1( rdof, jacInv_l );
    dBdx_r = tk::eval_dBdx_p1( rdof, jacInv_r );
    if(ndof > 4) {
      // Since eval_dBdx_p2 takes a coordgp of dimension ng by 3
      // I need to add a row of zeros to my coordgp (UNTESTED)
      std::array< std::vector< tk::real >, 3 > coordgp_aux;
      coordgp_aux[0][igp] = coordgp[0][igp];
      coordgp_aux[1][igp] = coordgp[1][igp];
      coordgp_aux[2][igp] = 0.0;
      tk::eval_dBdx_p2(igp, coordgp_aux, jacInv_l, dBdx_l);
      tk::eval_dBdx_p2(igp, coordgp_aux, jacInv_r, dBdx_r);
    }
  }
  // Loop through materials
  for (std::size_t k=0; k<nmat; ++k)
  {
    if (solidx[k] > 0)
    {
      // Compute all derivatives
      std::array< std::array< std::array< tk::real, 3 >, 3 >, 3 >
        deriv_l, deriv_r;
      std::size_t defIdx;
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
        {
          deriv_l[0][i][j] = 0.0;
          deriv_l[1][i][j] = 0.0;
          deriv_l[2][i][j] = 0.0;
          deriv_r[0][i][j] = 0.0;
          deriv_r[1][i][j] = 0.0;
          deriv_r[2][i][j] = 0.0;
          for (std::size_t jdof=0; jdof<rdof; ++jdof)
          {
            defIdx = deformDofIdx(nmat,solidx[k],i,j,rdof,jdof);
            deriv_l[0][i][j] += U(el,defIdx)*dBdx_l[0][jdof];
            deriv_l[1][i][j] += U(el,defIdx)*dBdx_l[1][jdof];
            deriv_l[2][i][j] += U(el,defIdx)*dBdx_l[2][jdof];
            deriv_r[0][i][j] += U(er,defIdx)*dBdx_r[0][jdof];
            deriv_r[1][i][j] += U(er,defIdx)*dBdx_r[1][jdof];
            deriv_r[2][i][j] += U(er,defIdx)*dBdx_r[2][jdof];
          }
        }

      // Compute the source terms
      tk::real alpha_l = U(el,volfracDofIdx(nmat,k,rdof,0));
      tk::real alpha_r = U(er,volfracDofIdx(nmat,k,rdof,0));
      tk::real sl, sr;
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
        {
          // Compute source
          sl = alpha_l*D*(fn[(j+1)%3]*(deriv_l[j][i][(j+1)%3]
                                      -deriv_l[(j+1)%3][i][j])
                         -fn[(j+2)%3]*(deriv_l[(j+2)%3][i][j]
                                      -deriv_l[j][i][(j+2)%3]));
          sr = alpha_r*D*(fn[(j+1)%3]*(deriv_r[j][i][(j+1)%3]
                                      -deriv_r[(j+1)%3][i][j])
                         -fn[(j+2)%3]*(deriv_r[(j+2)%3][i][j]
                                      -deriv_r[j][i][(j+2)%3]));
          // Add to flux
          fl[deformIdx(nmat,solidx[k],i,j)] += 0.5*(sl+sr);
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
