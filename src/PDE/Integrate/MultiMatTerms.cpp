// *****************************************************************************
/*!
  \file      src/PDE/Integrate/MultiMatTerms.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing volume integrals of multi-material terms
     using DG methods
  \details   This file contains functionality for computing volume integrals of
     non-conservative and pressure relaxation terms that appear in the
     multi-material hydrodynamic equations, using the discontinuous Galerkin
     method for various orders of numerical representation.
*/
// *****************************************************************************

#include "QuinoaConfig.hpp"

#include <lapacke.h>

#include "MultiMatTerms.hpp"
#include "Vector.hpp"
#include "Quadrature.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "Reconstruction.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "EoS/GetMatProp.hpp"

namespace inciter {
extern ctr::InputDeck g_inputdeck;
}

namespace tk {

void
nonConservativeInt( [[maybe_unused]] ncomp_t system,
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
                    const std::vector< std::vector< tk::real > >& riemannDeriv,
                    const std::vector< std::size_t >& ndofel,
                    Fields& R,
                    int intsharp )
// *****************************************************************************
//  Compute volume integrals for multi-material DG
//! \details This is called for multi-material DG, computing volume integrals of
//!   terms in the volume fraction and energy equations, which do not exist in
//!   the single-material flow formulation (for `CompFlow` DG). For further
//!   details see Pelanti, M., & Shyue, K. M. (2019). A numerical model for
//!   multiphase liquid–vapor–gas flows with interfaces and cavitation.
//!   International Journal of Multiphase Flow, 113, 208-230.
//! \param[in] system Equation system index
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] nelem Total number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitive quantities at recent time step
//! \param[in] riemannDeriv Derivatives of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms
//! \param[in] ndofel Vector of local number of degrees of freedome
//! \param[in,out] R Right-hand side vector added to
//! \param[in] intsharp Interface reconstruction indicator
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::densityIdx;
  using inciter::momentumIdx;
  using inciter::energyIdx;
  using inciter::velocityIdx;
  using inciter::deformIdx;

  const auto& solidx = inciter::g_inputdeck.get< tag::param, tag::multimat,
    tag::matidxmap >().template get< tag::solidx >();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

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

    // Compute the derivatives of basis function for DG(P1)
    std::array< std::vector<tk::real>, 3 > dBdx;
    if (ndofel[e] > 1)
      dBdx = eval_dBdx_p1( ndofel[e], jacInv );

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      if (ndofel[e] > 4)
        eval_dBdx_p2( igp, coordgp, jacInv, dBdx );

      // If an rDG method is set up (P0P1), then, currently we compute the P1
      // basis functions and solutions by default. This implies that P0P1 is
      // unsupported in the p-adaptive DG (PDG).
      std::size_t dof_el;
      if (rdof > ndof)
      {
        dof_el = rdof;
      }
      else
      {
        dof_el = ndofel[e];
      }

      // Compute the basis function
      auto B =
        eval_basis( dof_el, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp] );

      auto wt = wgp[igp] * geoElem(e, 0);

      auto state = evalPolynomialSol(system, mat_blk, intsharp, ncomp, nprim,
        rdof, nmat, e, dof_el, inpoel, coord, geoElem,
        {{coordgp[0][igp], coordgp[1][igp], coordgp[2][igp]}}, B, U, P);

      // get bulk properties
      tk::real rhob(0.0);
      for (std::size_t k=0; k<nmat; ++k)
          rhob += state[densityIdx(nmat, k)];

      // get the velocity vector
      std::array< tk::real, 3 > vel{{ state[ncomp+velocityIdx(nmat, 0)],
                                      state[ncomp+velocityIdx(nmat, 1)],
                                      state[ncomp+velocityIdx(nmat, 2)] }};

      std::vector< tk::real > ymat(nmat, 0.0);
      std::array< tk::real, 3 > dap{{0.0, 0.0, 0.0}};
      for (std::size_t k=0; k<nmat; ++k)
      {
        ymat[k] = state[densityIdx(nmat, k)]/rhob;

        for (std::size_t idir=0; idir<3; ++idir)
          dap[idir] += riemannDeriv[3*k+idir][e];
      }

      // compute non-conservative terms
      std::vector< std::vector< tk::real > > ncf
        (ncomp, std::vector<tk::real>(ndof,0.0));

      for (std::size_t idir=0; idir<3; ++idir)
        for(std::size_t idof=0; idof<ndof; ++idof)
          ncf[momentumIdx(nmat, idir)][idof] = 0.0;

      for (std::size_t k=0; k<nmat; ++k)
      {
        // evaluate non-conservative term for energy equation
        for(std::size_t idof=0; idof<ndof; ++idof)
        {
          ncf[densityIdx(nmat, k)][idof] = 0.0;
          for (std::size_t idir=0; idir<3; ++idir)
            ncf[energyIdx(nmat, k)][idof] -= vel[idir] * ( ymat[k]*dap[idir]
                                                  - riemannDeriv[3*k+idir][e] );
        }

        // Evaluate non-conservative term for volume fraction equation:
        // Here we make an assumption that the derivative of Riemann velocity
        // times the basis function is constant. Therefore, when P0P1/DGP1/DGP2
        // are used for constant velocity problems, the discretization is
        // consistent. However, for a general problem with varying velocity,
        // there will be errors since the said derivative is not constant.
        // A discretization that solves this issue has not been implemented yet.
        // Nevertheless, this does not affect high-order accuracy in
        // single material regions for problems with sharp interfaces. Since
        // volume fractions are nearly constant in such regions, using
        // high-order for volume fractions does not show any benefits over
        // THINC. Therefore, for such problems, we only use FV for the volume
        // fractions, and these non-conservative high-order terms do not need
        // to be computed.
        // In summary, high-order discretization for non-conservative terms in
        // volume fraction equations is avoided for sharp interface problems.
        if (ndof <= 4 || intsharp == 1) {
          for(std::size_t idof=0; idof<ndof; ++idof)
            ncf[volfracIdx(nmat, k)][idof] = state[volfracIdx(nmat, k)]
                                           * riemannDeriv[3*nmat+idof][e];
        } else if (intsharp == 0) {     // If DGP2 without THINC
          // DGP2 is discretized differently than DGP1/FV to guarantee 3rd order
          // convergence for the testcases with uniform and constant velocity.

          // P0 contributions for all equations
          for(std::size_t idof=0; idof<ndof; ++idof)
          ncf[volfracIdx(nmat, k)][idof] = state[volfracIdx(nmat, k)]
                                         * riemannDeriv[3*nmat][e] * B[idof];
          // High order contributions
          for(std::size_t idof=1; idof<ndof; ++idof)
            for(std::size_t idir=0; idir<3; ++idir)
            ncf[volfracIdx(nmat, k)][idof] += state[volfracIdx(nmat, k)]
                                            * vel[idir] * dBdx[idir][idof];
        }

        if (solidx[k] > 0)
        {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              for(std::size_t idof=0; idof<ndof; ++idof)
              {
                ncf[deformIdx(nmat,solidx[k],i,j)][idof] = 0.0;
                for (std::size_t l=0; l<3; ++l)
                {
                  ncf[deformIdx(nmat,solidx[k],i,j)][idof] +=
                    state[volfracIdx(nmat, k)]*(
                    riemannDeriv[3*nmat+ndof+3*3*9*k+3*3*(3*i+j)+3*l+l][e]
                   -riemannDeriv[3*nmat+ndof+3*3*9*k+3*3*(3*i+l)+3*l+j][e]);
                  }
              }
        }
      }

      updateRhsNonCons( ncomp, nmat, ndof, ndofel[e], wt, e, B, dBdx, ncf, R );
    }
  }
}

void
updateRhsNonCons(
  ncomp_t ncomp,
  const std::size_t nmat,
  const std::size_t ndof,
  [[maybe_unused]] const std::size_t ndof_el,
  const tk::real wt,
  const std::size_t e,
  const std::vector<tk::real>& B,
  [[maybe_unused]] const std::array< std::vector<tk::real>, 3 >& dBdx,
  const std::vector< std::vector< tk::real > >& ncf,
  Fields& R )
// *****************************************************************************
//  Update the rhs by adding the non-conservative term integrals
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] nmat Number of materials
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_el Number of degrees of freedom for local element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] e Element index
//! \param[in] B Basis function evaluated at local quadrature point
//! \param[in] dBdx Vector of basis function derivatives
//! \param[in] ncf Vector of non-conservative terms
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::energyIdx;
  using inciter::deformIdx;
  using inciter::volfracDofIdx;
  using inciter::energyDofIdx;
  using inciter::deformDofIdx;

  const auto& solidx = inciter::g_inputdeck.get< tag::param, tag::multimat,
    tag::matidxmap >().template get< tag::solidx >();

  //Assert( dBdx[0].size() == ndof_el,
  //        "Size mismatch for basis function derivatives" );
  //Assert( dBdx[1].size() == ndof_el,
  //        "Size mismatch for basis function derivatives" );
  //Assert( dBdx[2].size() == ndof_el,
  //        "Size mismatch for basis function derivatives" );
  //Assert( ncf.size() == ncomp,
  //        "Size mismatch for non-conservative term" );
  Assert( ncf.size() == ncomp, "Size mismatch for non-conservative term" );

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    R(e, mark) += wt * ncf[c][0];
  }

  if( ndof_el > 1)
  {
    // Update rhs with distributions from volume fraction and energy equations
    for (std::size_t k=0; k<nmat; ++k)
    {
      for(std::size_t idof = 1; idof < ndof; idof++)
      {
        R(e, volfracDofIdx(nmat,k,ndof,idof)) +=
          wt * ncf[volfracIdx(nmat,k)][idof];
        R(e, energyDofIdx(nmat,k,ndof,idof)) +=
          wt * ncf[energyIdx(nmat,k)][idof] * B[idof];
        for(std::size_t i=0; i<3; ++i)
          for(std::size_t j=0; j<3; ++j)
            R(e, deformDofIdx(nmat,solidx[k],i,j,ndof,idof)) +=
              wt * ncf[deformIdx(nmat,solidx[k],i,j)][idof] * B[idof];
      }
    }
  }
}

void
nonConservativeIntFV(
  ncomp_t system,
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t rdof,
  const std::size_t nelem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const Fields& geoElem,
  const Fields& U,
  const Fields& P,
  const std::vector< std::vector< tk::real > >& riemannDeriv,
  Fields& R )
// *****************************************************************************
//  Compute volume integrals of non-conservative terms for multi-material FV
//! \details This is called for multi-material FV, computing volume integrals of
//!   terms in the volume fraction and energy equations, which do not exist in
//!   the single-material flow formulation (for `CompFlow`). For further
//!   details see Pelanti, M., & Shyue, K. M. (2019). A numerical model for
//!   multiphase liquid–vapor–gas flows with interfaces and cavitation.
//!   International Journal of Multiphase Flow, 113, 208-230.
//! \param[in] system Equation system index
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] nelem Total number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitive quantities at recent time step
//! \param[in] riemannDeriv Derivatives of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms
//! \param[in,out] R Right-hand side vector added to
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::densityIdx;
  using inciter::momentumIdx;
  using inciter::energyIdx;
  using inciter::velocityIdx;

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  // compute volume integrals
  for (std::size_t e=0; e<nelem; ++e)
  {
    // Compute the basis function
    std::vector< tk::real > B(rdof, 0.0);
    B[0] = 1.0;

    auto state = evalPolynomialSol(system, mat_blk, 0, ncomp, nprim,
      rdof, nmat, e, rdof, inpoel, coord, geoElem,
      {{0.25, 0.25, 0.25}}, B, U, P);

    // get bulk properties
    tk::real rhob(0.0);
    for (std::size_t k=0; k<nmat; ++k)
        rhob += state[densityIdx(nmat, k)];

    // get the velocity vector
    std::array< tk::real, 3 > vel{{ state[ncomp+velocityIdx(nmat, 0)],
                                    state[ncomp+velocityIdx(nmat, 1)],
                                    state[ncomp+velocityIdx(nmat, 2)] }};

    std::vector< tk::real > ymat(nmat, 0.0);
    std::array< tk::real, 3 > dap{{0.0, 0.0, 0.0}};
    for (std::size_t k=0; k<nmat; ++k)
    {
      ymat[k] = state[densityIdx(nmat, k)]/rhob;

      for (std::size_t idir=0; idir<3; ++idir)
        dap[idir] += riemannDeriv[3*k+idir][e];
    }

    // compute non-conservative terms
    std::vector< tk::real > ncf(ncomp, 0.0);

    for (std::size_t idir=0; idir<3; ++idir)
      ncf[momentumIdx(nmat, idir)] = 0.0;

    for (std::size_t k=0; k<nmat; ++k)
    {
      // evaluate non-conservative term for energy equation
      ncf[densityIdx(nmat, k)] = 0.0;
      for (std::size_t idir=0; idir<3; ++idir)
        ncf[energyIdx(nmat, k)] -= vel[idir] * ( ymat[k]*dap[idir]
                                              - riemannDeriv[3*k+idir][e] );

      // evaluate non-conservative term for volume fraction equation
      ncf[volfracIdx(nmat, k)] = state[volfracIdx(nmat, k)]
        * riemannDeriv[3*nmat][e];
    }

    for (ncomp_t c=0; c<ncomp; ++c)
    {
      R(e, c) += geoElem(e,0) * ncf[c];
    }
  }
}

void
pressureRelaxationInt( ncomp_t system,
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
                       const tk::real ct,
                       Fields& R,
                       int intsharp )
// *****************************************************************************
//  Compute volume integrals of pressure relaxation terms in multi-material DG
//! \details This is called for multi-material DG to compute volume integrals of
//!   finite pressure relaxation terms in the volume fraction and energy
//!   equations, which do not exist in the single-material flow formulation (for
//!   `CompFlow` DG). For further details see Dobrev, V. A., Kolev, T. V.,
//!   Rieben, R. N., & Tomov, V. Z. (2016). Multi‐material closure model for
//!   high‐order finite element Lagrangian hydrodynamics. International Journal
//!   for Numerical Methods in Fluids, 82(10), 689-706.
//! \param[in] system Equation system index
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] nelem Total number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitive quantities at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedome
//! \param[in] ct Pressure relaxation time-scale for this system
//! \param[in,out] R Right-hand side vector added to
//! \param[in] intsharp Interface reconstruction indicator
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::densityIdx;
  using inciter::momentumIdx;
  using inciter::energyIdx;
  using inciter::pressureIdx;
  using inciter::velocityIdx;
  using inciter::deformIdx;

  const auto& solidx = inciter::g_inputdeck.get< tag::param, tag::multimat,
    tag::matidxmap >().template get< tag::solidx >();

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  // compute volume integrals
  for (std::size_t e=0; e<nelem; ++e)
  {
    auto dx = geoElem(e,4)/2.0;
    auto ng = NGvol(ndofel[e]);

    // arrays for quadrature points
    std::array< std::vector< real >, 3 > coordgp;
    std::vector< real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    coordgp[2].resize( ng );
    wgp.resize( ng );

    GaussQuadratureTet( ng, coordgp, wgp );

    // Compute the derivatives of basis function for DG(P1)
    std::array< std::vector<real>, 3 > dBdx;

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // If an rDG method is set up (P0P1), then, currently we compute the P1
      // basis functions and solutions by default. This implies that P0P1 is
      // unsupported in the p-adaptive DG (PDG).
      std::size_t dof_el;
      if (rdof > ndof)
      {
        dof_el = rdof;
      }
      else
      {
        dof_el = ndofel[e];
      }

      // Compute the basis function
      auto B =
        eval_basis( dof_el, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp] );

      auto wt = wgp[igp] * geoElem(e, 0);

      auto state = evalPolynomialSol(system, mat_blk, intsharp, ncomp, nprim,
        rdof, nmat, e, dof_el, inpoel, coord, geoElem,
        {{coordgp[0][igp], coordgp[1][igp], coordgp[2][igp]}}, B, U, P);

      // get bulk properties
      real rhob(0.0);
      for (std::size_t k=0; k<nmat; ++k)
        rhob += state[densityIdx(nmat, k)];

      // get pressures and bulk modulii
      real pb(0.0), nume(0.0), deno(0.0), trelax(0.0);
      std::vector< real > apmat(nmat, 0.0), kmat(nmat, 0.0);
      for (std::size_t k=0; k<nmat; ++k)
      {
        real arhomat = state[densityIdx(nmat, k)];
        real alphamat = state[volfracIdx(nmat, k)];
        apmat[k] = state[ncomp+pressureIdx(nmat, k)];
        real amat = 0.0;
        if (solidx[k] > 0)
        {
          std::array< std::array< tk::real, 3 >, 3 > ag;
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              ag[i][j] = state[deformIdx(nmat,solidx[k],i,j)];
          auto agrot = tk::rotateTensor(ag, {{1.0, 0.0, 0.0}});
          amat = mat_blk[k].compute< inciter::EOS::soundspeed >( arhomat,
            apmat[k], alphamat, k, 1.0, agrot);
          agrot = tk::rotateTensor(ag, {{0.0, 1.0, 0.0}});
          amat = std::max(amat, mat_blk[k].compute< inciter::EOS::soundspeed >(
            arhomat, apmat[k], alphamat, k, 1.0, agrot));
          agrot = tk::rotateTensor(ag, {{0.0, 0.0, 1.0}});
          amat = std::max(amat, mat_blk[k].compute< inciter::EOS::soundspeed >(
            arhomat, apmat[k], alphamat, k, 1.0, agrot));
        }
        else
        {
          amat = mat_blk[k].compute< inciter::EOS::soundspeed >( arhomat,
            apmat[k], alphamat, k );
        }
        kmat[k] = arhomat * amat * amat;
        pb += apmat[k];

        // relaxation parameters
        trelax = std::max(trelax, ct*dx/amat);
        nume += alphamat * apmat[k] / kmat[k];
        deno += alphamat * alphamat / kmat[k];
      }
      auto p_relax = nume/deno;

      // compute pressure relaxation terms
      std::vector< real > s_prelax(ncomp, 0.0);
      for (std::size_t k=0; k<nmat; ++k)
      {
        auto s_alpha = (apmat[k]-p_relax*state[volfracIdx(nmat, k)])
          * (state[volfracIdx(nmat, k)]/kmat[k]) / trelax;
        s_prelax[volfracIdx(nmat, k)] = s_alpha;
        s_prelax[energyIdx(nmat, k)] = - pb*s_alpha;
        if (solidx[k] > 0)
          for (size_t i=0; i<3; ++i)
            for (size_t j=0; j<3; ++j)
            {
              tk::real gij =
                state[deformIdx(nmat,solidx[k],i,j)]/state[volfracIdx(nmat, k)];
              s_prelax[deformIdx(nmat,solidx[k],i,j)] = gij * s_alpha;
            }
      }

      updateRhsPre( ncomp, ndof, dof_el, wt, e, B, s_prelax, R );
    }
  }
}

void
updateRhsPre(
  ncomp_t ncomp,
  const std::size_t ndof,
  [[maybe_unused]] const std::size_t ndof_el,
  const tk::real wt,
  const std::size_t e,
  const std::vector< tk::real >& B,
  std::vector< tk::real >& ncf,
  Fields& R )
// *****************************************************************************
//  Update the rhs by adding the pressure relaxation integrals
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_el Number of degrees of freedom for local element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] e Element index
//! \param[in] B Basis function evaluated at local quadrature point
//! \param[in] ncf Vector of non-conservative terms
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  //Assert( dBdx[0].size() == ndof_el,
  //        "Size mismatch for basis function derivatives" );
  //Assert( dBdx[1].size() == ndof_el,
  //        "Size mismatch for basis function derivatives" );
  //Assert( dBdx[2].size() == ndof_el,
  //        "Size mismatch for basis function derivatives" );
  //Assert( ncf.size() == ncomp,
  //        "Size mismatch for non-conservative term" );
  Assert( ncf.size() == ncomp, "Size mismatch for non-conservative term" );

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    for(std::size_t idof = 0; idof < ndof; idof++)
      R(e, mark+idof) += wt * ncf[c] * B[idof];
  }
}

void
pressureRelaxationIntFV(
  ncomp_t system,
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t rdof,
  const std::size_t nelem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const Fields& geoElem,
  const Fields& U,
  const Fields& P,
  const tk::real ct,
  Fields& R )
// *****************************************************************************
//  Compute volume integrals of pressure relaxation terms in multi-material FV
//! \details This is called for multi-material FV to compute volume integrals of
//!   finite pressure relaxation terms in the volume fraction and energy
//!   equations, which do not exist in the single-material flow formulation (for
//!   `CompFlow`). For further details see Dobrev, V. A., Kolev, T. V.,
//!   Rieben, R. N., & Tomov, V. Z. (2016). Multi‐material closure model for
//!   high‐order finite element Lagrangian hydrodynamics. International Journal
//!   for Numerical Methods in Fluids, 82(10), 689-706.
//! \param[in] system Equation system index
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] nelem Total number of elements
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitive quantities at recent time step
//! \param[in] ct Pressure relaxation time-scale for this system
//! \param[in,out] R Right-hand side vector added to
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::energyIdx;
  using inciter::pressureIdx;
  using inciter::velocityIdx;
  using inciter::densityIdx;

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  // compute volume integrals
  for (std::size_t e=0; e<nelem; ++e)
  {
    auto dx = geoElem(e,4)/2.0;

    // Compute the basis function
    std::vector< tk::real > B(rdof, 0.0);
    B[0] = 1.0;

    auto state = evalPolynomialSol(system, mat_blk, 0, ncomp, nprim,
      rdof, nmat, e, rdof, inpoel, coord, geoElem,
      {{0.25, 0.25, 0.25}}, B, U, P);

    // get bulk properties
    real rhob(0.0);
    for (std::size_t k=0; k<nmat; ++k)
      rhob += state[densityIdx(nmat, k)];

    // get pressures and bulk modulii
    real pb(0.0), nume(0.0), deno(0.0), trelax(0.0);
    std::vector< real > apmat(nmat, 0.0), kmat(nmat, 0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      real arhomat = state[densityIdx(nmat, k)];
      real alphamat = state[volfracIdx(nmat, k)];
      apmat[k] = state[ncomp+pressureIdx(nmat, k)];
      real amat = mat_blk[k].compute< inciter::EOS::soundspeed >( arhomat,
        apmat[k], alphamat, k );
      kmat[k] = arhomat * amat * amat;
      pb += apmat[k];

      // relaxation parameters
      trelax = std::max(trelax, ct*dx/amat);
      nume += alphamat * apmat[k] / kmat[k];
      deno += alphamat * alphamat / kmat[k];
    }
    auto p_relax = nume/deno;

    // compute pressure relaxation terms
    std::vector< real > s_prelax(ncomp, 0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      auto s_alpha = (apmat[k]-p_relax*state[volfracIdx(nmat, k)])
        * (state[volfracIdx(nmat, k)]/kmat[k]) / trelax;
      s_prelax[volfracIdx(nmat, k)] = s_alpha;
      s_prelax[energyIdx(nmat, k)] = - pb*s_alpha;
    }

    for (ncomp_t c=0; c<ncomp; ++c)
    {
      R(e, c) += geoElem(e,0) * s_prelax[c];
    }
  }
}

std::vector< std::vector< tk::real > >
solvevriem( std::size_t nelem,
            const std::vector< std::vector< tk::real > >& vriem,
            const std::vector< std::vector< tk::real > >& riemannLoc )
// *****************************************************************************
//  Solve the reconstruct velocity used for volume fraction equation by
//  Least square method
//! \param[in] nelem Numer of elements
//! \param[in,out] vriem Vector of the riemann velocity
//! \param[in,out] riemannLoc Vector of coordinates where Riemann velocity data
//!   is available
//! \return Vector of Riemann velocity polynomial solution
// *****************************************************************************
{
  std::vector< std::vector< tk::real > >
    vriempoly( nelem, std::vector<tk::real>(12,0.0) );

  for (std::size_t e=0; e<nelem; ++e)
  {
    // Use the normal method to construct the linear system A^T * A * x = u
    auto numgp = riemannLoc[e].size()/3;
    std::vector< std::vector< tk::real > > A(numgp,
                                             std::vector<tk::real>(4, 1.0));

    for(std::size_t k = 0; k < numgp; k++)
    {
      auto mark = k * 3;
      A[k][1] = riemannLoc[e][mark];
      A[k][2] = riemannLoc[e][mark+1];
      A[k][3] = riemannLoc[e][mark+2];
    }

    for(std::size_t idir = 0; idir < 3; idir++)
    {
      double AA_T[4*4], u[4];

      for(std::size_t i = 0; i < 4; i++)
        for(std::size_t j = 0; j < 4; j++)
        {
          auto id = 4 * i + j;
          AA_T[id] = 0;
          for(std::size_t k = 0; k < numgp; k++)
            AA_T[id] += A[k][i] * A[k][j];
        }

      std::vector<tk::real> vel(numgp, 1.0);
      for(std::size_t k = 0; k < numgp; k++)
      {
        auto mark = k * 3 + idir;
        vel[k] = vriem[e][mark];
      }
      for(std::size_t k = 0; k < 4; k++)
      {
        u[k] = 0;
        for(std::size_t i = 0; i < numgp; i++)
          u[k] += A[i][k] * vel[i];
      }
 
      lapack_int IPIV[4];
      #ifndef NDEBUG
      lapack_int info =
      #endif
        LAPACKE_dsysv( LAPACK_ROW_MAJOR, 'U', 4, 1, AA_T, 4, IPIV, u, 1 );
      Assert( info == 0, "Error in linear system solver" );

      auto idirmark = idir * 4;
      for(std::size_t k = 0; k < 4; k++)
        vriempoly[e][idirmark+k] = u[k];
    }
  }
  return vriempoly;
}

void evaluRiemann( ncomp_t ncomp,
                   const int e_left,
                   const int e_right,
                   const std::size_t nmat,
                   const std::vector< tk::real >& fl,
                   const std::array< tk::real, 3 >& fn,
                   const std::array< tk::real, 3 >& gp,
                   const std::array< std::vector< tk::real >, 2 >& state,
                   std::vector< std::vector< tk::real > >& vriem,
                   std::vector< std::vector< tk::real > >& riemannLoc )
// *****************************************************************************
//  Compute the riemann velocity at the interface
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] e_left Index for the left element
//! \param[in] e_right Index for the right element
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] fn Face/Surface normal
//! \param[in] gp Gauss points coordinates
//! \param[in] fl Surface flux
//! \param[in] state Vector of state variables for left and right side
//! \param[in,out] vriem Vector of the riemann velocity
//! \param[in,out] riemannLoc Vector of coordinates where Riemann velocity data
//!   is available
// *****************************************************************************
{
  using inciter::densityIdx;
  using inciter::momentumIdx;

  std::size_t el(0), er(0);
  el = static_cast< std::size_t >(e_left);
  if(e_right != -1)
    er = static_cast< std::size_t >(e_right);

  riemannLoc[el].push_back( gp[0] );
  riemannLoc[el].push_back( gp[1] );
  riemannLoc[el].push_back( gp[2] );

  if(e_right != -1)
  {
    riemannLoc[er].push_back( gp[0] );
    riemannLoc[er].push_back( gp[1] );
    riemannLoc[er].push_back( gp[2] );
  }

  tk::real rhobl(0.0), rhobr(0.0);
  for (std::size_t k=0; k<nmat; ++k)
  {
    rhobl += state[0][densityIdx(nmat, k)];
    rhobr += state[1][densityIdx(nmat, k)];
  }

  auto ul = state[0][momentumIdx(nmat, 0)] / rhobl;
  auto vl = state[0][momentumIdx(nmat, 1)] / rhobl;
  auto wl = state[0][momentumIdx(nmat, 2)] / rhobl;

  auto ur = state[1][momentumIdx(nmat, 0)] / rhobr;
  auto vr = state[1][momentumIdx(nmat, 1)] / rhobr;
  auto wr = state[1][momentumIdx(nmat, 2)] / rhobr;

  // Compute the normal velocities from left and right cells
  auto vnl = ul * fn[0] + vl * fn[1] + wl * fn[2];
  auto vnr = ur * fn[0] + vr * fn[1] + wr * fn[2];

  // The interface velocity is evaluated by adding the normal velocity which
  // is taken from the Riemann solver and the tangential velocity which is
  // evaluated as an average of the left and right cells
  auto urie = 0.5 * ((ul + ur) - fn[0] * (vnl + vnr)) + fl[ncomp+nmat] * fn[0];
  auto vrie = 0.5 * ((vl + vr) - fn[1] * (vnl + vnr)) + fl[ncomp+nmat] * fn[1];
  auto wrie = 0.5 * ((wl + wr) - fn[2] * (vnl + vnr)) + fl[ncomp+nmat] * fn[2];

  vriem[el].push_back(urie);
  vriem[el].push_back(vrie);
  vriem[el].push_back(wrie);

  if(e_right != -1)
  {
    vriem[er].push_back(urie);
    vriem[er].push_back(vrie);
    vriem[er].push_back(wrie);
  }
}

std::vector< std::array< tk::real, 3 > >
fluxTerms(
  std::size_t ncomp,
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const std::vector< tk::real >& ugp )
// *****************************************************************************
//  Compute the flux-function for the multimaterial PDEs
//! \param[in] ncomp Number of components in this PDE system
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] U Solution vector at recent time step
//! \return Flux vectors for all components in multi-material PDE system
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::densityIdx;
  using inciter::momentumIdx;
  using inciter::energyIdx;
  using inciter::velocityIdx;
  using inciter::pressureIdx;
  using inciter::deformIdx;

  const auto& solidx = inciter::g_inputdeck.get< tag::param, tag::multimat,
    tag::matidxmap >().template get< tag::solidx >();

  std::vector< std::array< tk::real, 3 > > fl( ncomp );

  if (inciter::haveSolid(nmat, solidx))
  {
    tk::real rho(0.0);
    for (std::size_t k=0; k<nmat; ++k)
      rho += ugp[densityIdx(nmat, k)];

    auto u = ugp[ncomp+velocityIdx(nmat,0)];
    auto v = ugp[ncomp+velocityIdx(nmat,1)];
    auto w = ugp[ncomp+velocityIdx(nmat,2)];

    std::vector< tk::real > al(nmat, 0.0);
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > ag, asig;
    std::array< std::array< tk::real, 3 >, 3 >
      sig {{ {{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}} }};
    for (std::size_t k=0; k<nmat; ++k)
    {
      al[k] = ugp[volfracIdx(nmat, k)];
      // inv deformation gradient and Cauchy stress tensors
      ag.push_back(inciter::getDeformGrad(nmat, k, ugp));
      asig.push_back(mat_blk[k].computeTensor< inciter::EOS::CauchyStress >(
        ugp[densityIdx(nmat, k)], u, v, w, ugp[energyIdx(nmat, k)],
        al[k], k, ag[k]));
      for (size_t i=0; i<3; ++i)
        for (size_t j=0; j<3; ++j)
          sig[i][j] += asig[k][i][j];
    }

    // conservative part of momentum flux
    fl[momentumIdx(nmat, 0)][0] = ugp[momentumIdx(nmat, 0)] * u - sig[0][0];
    fl[momentumIdx(nmat, 1)][0] = ugp[momentumIdx(nmat, 1)] * u - sig[0][1];
    fl[momentumIdx(nmat, 2)][0] = ugp[momentumIdx(nmat, 2)] * u - sig[0][2];

    fl[momentumIdx(nmat, 0)][1] = ugp[momentumIdx(nmat, 0)] * v - sig[1][0];
    fl[momentumIdx(nmat, 1)][1] = ugp[momentumIdx(nmat, 1)] * v - sig[1][1];
    fl[momentumIdx(nmat, 2)][1] = ugp[momentumIdx(nmat, 2)] * v - sig[1][2];

    fl[momentumIdx(nmat, 0)][2] = ugp[momentumIdx(nmat, 0)] * w - sig[2][0];
    fl[momentumIdx(nmat, 1)][2] = ugp[momentumIdx(nmat, 1)] * w - sig[2][1];
    fl[momentumIdx(nmat, 2)][2] = ugp[momentumIdx(nmat, 2)] * w - sig[2][2];

    for (std::size_t k=0; k<nmat; ++k)
    {
      // conservative part of volume-fraction flux
      fl[volfracIdx(nmat, k)][0] = 0.0;
      fl[volfracIdx(nmat, k)][1] = 0.0;
      fl[volfracIdx(nmat, k)][2] = 0.0;

      // conservative part of material continuity flux
      fl[densityIdx(nmat, k)][0] = u * ugp[densityIdx(nmat, k)];
      fl[densityIdx(nmat, k)][1] = v * ugp[densityIdx(nmat, k)];
      fl[densityIdx(nmat, k)][2] = w * ugp[densityIdx(nmat, k)];

      // conservative part of material total-energy flux
      fl[energyIdx(nmat, k)][0] = u * ugp[energyIdx(nmat, k)]
        - u * asig[k][0][0] - v * asig[k][1][0] - w * asig[k][2][0];
      fl[energyIdx(nmat, k)][1] = v * ugp[energyIdx(nmat, k)]
        - u * asig[k][0][1] - v * asig[k][1][1] - w * asig[k][2][1];
      fl[energyIdx(nmat, k)][2] = w * ugp[energyIdx(nmat, k)]
        - u * asig[k][0][2] - v * asig[k][1][2] - w * asig[k][2][2];

      // conservative part of material inverse deformation gradient
      if (solidx[k] > 0)
      {
        for (std::size_t i=0; i<3; ++i)
        {
          fl[deformIdx(nmat,solidx[k],i,0)][0] = u*ag[k][i][0];
          fl[deformIdx(nmat,solidx[k],i,0)][1] = v*ag[k][i][0];
          fl[deformIdx(nmat,solidx[k],i,0)][2] = w*ag[k][i][0];
          fl[deformIdx(nmat,solidx[k],i,1)][0] = u*ag[k][i][1];
          fl[deformIdx(nmat,solidx[k],i,1)][1] = v*ag[k][i][1];
          fl[deformIdx(nmat,solidx[k],i,1)][2] = w*ag[k][i][1];
          fl[deformIdx(nmat,solidx[k],i,2)][0] = u*ag[k][i][2];
          fl[deformIdx(nmat,solidx[k],i,2)][1] = v*ag[k][i][2];
          fl[deformIdx(nmat,solidx[k],i,2)][2] = w*ag[k][i][2];
        }
      }
    }
  }
  else
  {
    tk::real rho(0.0), p(0.0);
    for (std::size_t k=0; k<nmat; ++k)
      rho += ugp[densityIdx(nmat, k)];

    auto u = ugp[ncomp+velocityIdx(nmat,0)];
    auto v = ugp[ncomp+velocityIdx(nmat,1)];
    auto w = ugp[ncomp+velocityIdx(nmat,2)];

    std::vector< tk::real > apk( nmat, 0.0 );
    for (std::size_t k=0; k<nmat; ++k)
    {
      apk[k] = ugp[ncomp+pressureIdx(nmat,k)];
      p += apk[k];
    }

    // conservative part of momentum flux
    fl[momentumIdx(nmat, 0)][0] = ugp[momentumIdx(nmat, 0)] * u + p;
    fl[momentumIdx(nmat, 1)][0] = ugp[momentumIdx(nmat, 1)] * u;
    fl[momentumIdx(nmat, 2)][0] = ugp[momentumIdx(nmat, 2)] * u;

    fl[momentumIdx(nmat, 0)][1] = ugp[momentumIdx(nmat, 0)] * v;
    fl[momentumIdx(nmat, 1)][1] = ugp[momentumIdx(nmat, 1)] * v + p;
    fl[momentumIdx(nmat, 2)][1] = ugp[momentumIdx(nmat, 2)] * v;

    fl[momentumIdx(nmat, 0)][2] = ugp[momentumIdx(nmat, 0)] * w;
    fl[momentumIdx(nmat, 1)][2] = ugp[momentumIdx(nmat, 1)] * w;
    fl[momentumIdx(nmat, 2)][2] = ugp[momentumIdx(nmat, 2)] * w + p;

    for (std::size_t k=0; k<nmat; ++k)
    {
      // conservative part of volume-fraction flux
      fl[volfracIdx(nmat, k)][0] = 0.0;
      fl[volfracIdx(nmat, k)][1] = 0.0;
      fl[volfracIdx(nmat, k)][2] = 0.0;

      // conservative part of material continuity flux
      fl[densityIdx(nmat, k)][0] = u * ugp[densityIdx(nmat, k)];
      fl[densityIdx(nmat, k)][1] = v * ugp[densityIdx(nmat, k)];
      fl[densityIdx(nmat, k)][2] = w * ugp[densityIdx(nmat, k)];

      // conservative part of material total-energy flux
      auto hmat = ugp[energyIdx(nmat, k)] + apk[k];
      fl[energyIdx(nmat, k)][0] = u * hmat;
      fl[energyIdx(nmat, k)][1] = v * hmat;
      fl[energyIdx(nmat, k)][2] = w * hmat;
    }
  }

  return fl;
}

}// tk::
