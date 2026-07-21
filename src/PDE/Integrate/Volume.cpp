// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Volume.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing volume integrals for a system of PDEs in DG
     methods
  \details   This file contains functionality for computing volume integrals for
     a system of PDEs used in discontinuous Galerkin methods for various orders
     of numerical representation.
*/
// *****************************************************************************

#include "Volume.hpp"
#include "Vector.hpp"
#include "Quadrature.hpp"
#include "Reconstruction.hpp"
#include "ViscousTerms.hpp"

template< class ViscousTerms >
void
tk::volIntViscous(
            const ViscousTerms& viscousRhs,
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
            Fields& R)
// *****************************************************************************
//  Compute volume integrals for viscous DG
//! \tparam ViscousTerms Policy type that computes PDE-specific viscous RHS
//! \param[in] viscousRhs PDE-specific viscous residual policy
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Total number of degrees of freedom included reconstructed ones
//! \param[in] nelem Maximum number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  std::vector< tk::real > state(ncomp+nprim);

  Assert( ndof*ncomp == R.nprop(),
          "Mismatch in viscous RHS polynomial and component sizes" );

  std::vector<std::array<tk::real, 3>> visc_fl(ncomp);
  // compute volume integrals
  for (std::size_t e=0; e<nelem; ++e)
  {
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
      {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }}
    }};

    auto jacInv =
            inverseJacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

    auto dof_el = ndofel[e];

    // Compute the derivatives of basis function for second order terms
    std::array< std::vector<tk::real>, 3 > dBdx;
    for (std::size_t i=0; i<3; ++i) dBdx[i].resize( dof_el, 0 );
    eval_dBdx_p1( dof_el, jacInv, dBdx );

    std::vector< tk::real > B(dof_el);

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      if (dof_el > 4)
        eval_dBdx_p2( igp, coordgp, jacInv, dBdx );

      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordel, coordgp );

      // Compute the basis function
      eval_basis( dof_el, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp],
        B );

      auto wt = wgp[igp] * geoElem(e, 0);

      // volume fluxes
      if(dof_el > 1)
      {
        state = viscousRhs.stateAt(mat_blk, U, P, e, dof_el, B );
        // compute viscous flux
        viscousRhs.volumeFlux( mat_blk, ncomp, state, visc_fl ); //currently no-op

        update_rhs( ncomp, ndof, dof_el, wt, e, dBdx, visc_fl, R );
      }
    }
  }
}

void
tk::volInt( std::size_t nmat,
            real t,
            const std::vector< inciter::EOS >& mat_blk,
            const std::size_t ndof,
            const std::size_t rdof,
            const std::size_t nelem,
            const std::vector< std::size_t >& inpoel,
            const UnsMesh::Coords& coord,
            const Fields& geoElem,
            const FluxFn& flux,
            const VelFn& vel,
            const SrcFn& src,
            const Fields& U,
            const Fields& P,
            const std::vector< std::size_t >& ndofel,
            Fields& R,
            int intsharp )
// *****************************************************************************
//  Compute volume integrals for DG
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] t Physical time
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Total number of degrees of freedom included reconstructed ones
//! \param[in] nelem Maximum number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] flux Flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] src Source function to use
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in,out] R Right-hand side vector added to
//! \param[in] intsharp Interface compression tag, an optional argument, with
//!   default 0, so that it is unused for single-material and transport.
// *****************************************************************************
{
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  std::vector< tk::real > state(ncomp+nprim);

  // compute volume integrals
  for (std::size_t e=0; e<nelem; ++e)
  {
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
      {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }}
    }};

    auto jacInv =
            inverseJacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

    auto dof_el = ndofel[e];

    // Compute the derivatives of basis function for second order terms
    std::array< std::vector<tk::real>, 3 > dBdx;
    for (std::size_t i=0; i<3; ++i) dBdx[i].resize( dof_el, 0 );
    eval_dBdx_p1( dof_el, jacInv, dBdx );

    std::vector< tk::real > B(dof_el);

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      if (dof_el > 4)
        eval_dBdx_p2( igp, coordgp, jacInv, dBdx );

      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordel, coordgp );

      // Compute the basis function
      eval_basis( dof_el, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp],
        B );

      auto wt = wgp[igp] * geoElem(e, 0);

      // volume fluxes
      if(dof_el > 1)
      {
        evalPolynomialSol(mat_blk, intsharp, ncomp, nprim,
          rdof, nmat, e, ndofel[e], inpoel, coord, geoElem,
          {{coordgp[0][igp], coordgp[1][igp], coordgp[2][igp]}}, B, U, P, state);

        // evaluate prescribed velocity (if any)
        auto v = vel( ncomp, gp[0], gp[1], gp[2], t );

        // comput flux
        auto fl = flux( ncomp, mat_blk, state, v );

        update_rhs( ncomp, ndof, dof_el, wt, e, dBdx, fl, R );
      }

      // source terms
      std::vector< real > s(ncomp, 0.0);
      src( nmat, mat_blk, gp[0], gp[1], gp[2], t, s );

      update_rhs_src( ndof, ndofel[e], wt, e, B, s, R );
    }
  }
}

void
tk::update_rhs( ncomp_t ncomp,
                const std::size_t ndof,
                const std::size_t ndof_el,
                const tk::real wt,
                const std::size_t e,
                const std::array< std::vector<tk::real>, 3 >& dBdx,
                const std::vector< std::array< tk::real, 3 > >& fl,
                Fields& R )
// *****************************************************************************
//  Update the rhs by adding the source term integrals
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_el Number of degrees of freedom for local element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] e Element index
//! \param[in] dBdx Vector of basis function derivatives
//! \param[in] fl Vector of numerical flux
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  Assert( dBdx[0].size() == ndof_el,
    "Size mismatch for basis function derivatives" );
  Assert( dBdx[1].size() == ndof_el,
    "Size mismatch for basis function derivatives" );
  Assert( dBdx[2].size() == ndof_el,
    "Size mismatch for basis function derivatives" );
  Assert( fl.size() == ncomp, "Size mismatch for flux term" );

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    R(e, mark+1) +=
      wt * (fl[c][0]*dBdx[0][1] + fl[c][1]*dBdx[1][1] + fl[c][2]*dBdx[2][1]);
    R(e, mark+2) +=
      wt * (fl[c][0]*dBdx[0][2] + fl[c][1]*dBdx[1][2] + fl[c][2]*dBdx[2][2]);
    R(e, mark+3) +=
      wt * (fl[c][0]*dBdx[0][3] + fl[c][1]*dBdx[1][3] + fl[c][2]*dBdx[2][3]);

    if( ndof_el > 4 )
    {
      R(e, mark+4) +=
        wt * (fl[c][0]*dBdx[0][4] + fl[c][1]*dBdx[1][4] + fl[c][2]*dBdx[2][4]);
      R(e, mark+5) +=
        wt * (fl[c][0]*dBdx[0][5] + fl[c][1]*dBdx[1][5] + fl[c][2]*dBdx[2][5]);
      R(e, mark+6) +=
        wt * (fl[c][0]*dBdx[0][6] + fl[c][1]*dBdx[1][6] + fl[c][2]*dBdx[2][6]);
      R(e, mark+7) +=
        wt * (fl[c][0]*dBdx[0][7] + fl[c][1]*dBdx[1][7] + fl[c][2]*dBdx[2][7]);
      R(e, mark+8) +=
        wt * (fl[c][0]*dBdx[0][8] + fl[c][1]*dBdx[1][8] + fl[c][2]*dBdx[2][8]);
      R(e, mark+9) +=
        wt * (fl[c][0]*dBdx[0][9] + fl[c][1]*dBdx[1][9] + fl[c][2]*dBdx[2][9]);
    }
  }
}

void
tk::update_rhs_src(
  const std::size_t ndof,
  const std::size_t ndof_el,
  const tk::real wt,
  const std::size_t e,
  const std::vector< tk::real >& B,
  const std::vector< tk::real >& s,
  Fields& R )
// *****************************************************************************
//  Update the rhs by adding the source term integrals
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_el Number of degrees of freedom for local element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] e Element index
//! \param[in] B Vector of basis functions
//! \param[in] s Source term vector
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  //// Removed since this can be called for p0p1 where B.size() > ndof_el
  //Assert( B.size() == ndof_el, "Size mismatch for basis function" );

  for (ncomp_t c=0; c<s.size(); ++c)
  {
    auto mark = c*ndof;
    R(e, mark)   += wt * s[c];

    if ( ndof_el > 1 )
    {
      R(e, mark+1) += wt * s[c] * B[1];
      R(e, mark+2) += wt * s[c] * B[2];
      R(e, mark+3) += wt * s[c] * B[3];

      if( ndof_el > 4 )
      {
        R(e, mark+4) += wt * s[c] * B[4];
        R(e, mark+5) += wt * s[c] * B[5];
        R(e, mark+6) += wt * s[c] * B[6];
        R(e, mark+7) += wt * s[c] * B[7];
        R(e, mark+8) += wt * s[c] * B[8];
        R(e, mark+9) += wt * s[c] * B[9];
      }
    }
  }
}

void
tk::volInt_constP(
  std::size_t nmat,
  real t,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t ndof,
  const std::size_t rdof,
  const std::size_t nelem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const Fields& geoElem,
  const FluxFn& flux,
  const VelFn& vel,
  const SrcFn& src,
  const Fields& U,
  const Fields& P,
  Fields& R,
  int intsharp )
// *****************************************************************************
//  Compute volume integrals for const-order DG (not p-adaptive)
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] t Physical time
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Total number of degrees of freedom included reconstructed ones
//! \param[in] nelem Maximum number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] flux Flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] src Source function to use
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in,out] R Right-hand side vector added to
//! \param[in] intsharp Interface compression tag, an optional argument, with
//!   default 0, so that it is unused for single-material and transport.
// *****************************************************************************
{
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  // Quadrature pointd
  auto ng = tk::NGvol(ndof);

  // arrays for quadrature points
  std::array< std::vector< real >, 3 > coordgp;
  std::vector< real > wgp;

  coordgp[0].resize( ng );
  coordgp[1].resize( ng );
  coordgp[2].resize( ng );
  wgp.resize( ng );

  GaussQuadratureTet( ng, coordgp, wgp );

  // Allocate memory
  std::vector< tk::real > B(rdof), state(ncomp+nprim), s(ncomp, 0.0);
  std::array< std::vector<tk::real>, 3 > dBdx;
  for (std::size_t i=0; i<3; ++i) dBdx[i].resize( ndof, 0 );

  // compute volume integrals
  for (std::size_t e=0; e<nelem; ++e)
  {
    // Extract the element coordinates
    std::array< std::array< real, 3>, 4 > coordel {{
      {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
      {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
      {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
      {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }}
    }};

    auto jacInv =
            inverseJacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

    // Compute the derivatives of basis function for second order terms
    eval_dBdx_p1( ndof, jacInv, dBdx );

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      if (ndof > 4) eval_dBdx_p2( igp, coordgp, jacInv, dBdx );

      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordel, coordgp );

      // Compute the basis function
      eval_basis( rdof, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp], B );

      auto wt = wgp[igp] * geoElem(e, 0);

      // volume fluxes
      if (ndof > 1)
      {
        evalPolynomialSol(mat_blk, intsharp, ncomp, nprim,
          rdof, nmat, e, rdof, inpoel, coord, geoElem,
          {{coordgp[0][igp], coordgp[1][igp], coordgp[2][igp]}}, B, U, P, state);

        // evaluate prescribed velocity (if any)
        auto v = vel( ncomp, gp[0], gp[1], gp[2], t );

        // comput flux
        auto fl = flux( ncomp, mat_blk, state, v );

        update_rhs( ncomp, ndof, ndof, wt, e, dBdx, fl, R );
      }

      // source terms
      src( nmat, mat_blk, gp[0], gp[1], gp[2], t, s );

      update_rhs_src( ndof, ndof, wt, e, B, s, R );
    }
  }
}

void
tk::srcIntFV( const std::vector< inciter::EOS >& mat_blk,
              real t,
              const std::size_t nelem,
              const Fields& geoElem,
              const SrcFn& src,
              Fields& R,
              std::size_t nmat )
// *****************************************************************************
//  Compute source term integrals for DG
//! \param[in] mat_blk Material block EOS
//! \param[in] t Physical time
//! \param[in] nelem Maximum number of elements
//! \param[in] geoElem Element geometry array
//! \param[in] src Source function to use
//! \param[in,out] R Right-hand side vector computed
//! \param[in] nmat Number of materials. A default is set to 1, so that calling
//!   code for single material systems primitive quantities does not need to
//!   specify this argument.
// *****************************************************************************
{
  auto ncomp = R.nprop();

  for (std::size_t e=0; e<nelem; ++e)
  {
    // Compute the source term variable
    std::vector< real > s(ncomp, 0.0);
    src( nmat, mat_blk, geoElem(e,1), geoElem(e,2), geoElem(e,3), t, s );

    // Add the source term to the rhs
    for (ncomp_t c=0; c<ncomp; ++c)
    {
      R(e, c) += geoElem(e,0) * s[c];
    }
  }
}

void
tk::volIntViscousMultiSpecies(
  std::size_t nspec,
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
  Fields& R )
// *****************************************************************************
//  Compute volume integrals of viscous fluxes for multispecies flow
//! \param[in] nspec Number of species in this PDE system
//! \param[in] mat_blk Material EOS block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] nelem Maximum number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  if (ndof == 4) {
    MultiSpeciesViscousTermsDGP1 viscousRhs( nspec, rdof );
    volIntViscous( viscousRhs, mat_blk, ndof, rdof, nelem,
      inpoel, coord, geoElem, U, P, ndofel, R );
  }
  else
    Throw( "Viscous operators only implemented for scheme = 'dgp1'." );
}
