// *****************************************************************************
/*!
  \file      src/PDE/Limiter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Limiters for discontiunous Galerkin methods
  \details   This file contains functions that provide limiter function
    calculations for maintaining monotonicity near solution discontinuities
    for the DG discretization.
*/
// *****************************************************************************

#include <array>
#include <vector>

#include "FaceData.hpp"
#include "Vector.hpp"
#include "Limiter.hpp"
#include "DerivedData.hpp"
#include "Integrate/Quadrature.hpp"
#include "Integrate/Basis.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "EoS/EoS.hpp"
#include "PrefIndicator.hpp"
#include "Reconstruction.hpp"
#include "Integrate/Mass.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

void
WENO_P1( const std::vector< int >& esuel,
         tk::Fields& U )
// *****************************************************************************
//  Weighted Essentially Non-Oscillatory (WENO) limiter for DGP1
//! \param[in] esuel Elements surrounding elements
//! \param[in,out] U High-order solution vector which gets limited
//! \details This WENO function should be called for transport and compflow
//! \note This limiter function is experimental and untested. Use with caution.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();
  const auto cweight = inciter::g_inputdeck.get< tag::discr, tag::cweight >();
  auto nelem = esuel.size()/4;
  std::array< std::vector< tk::real >, 3 >
    limU {{ std::vector< tk::real >(nelem),
            std::vector< tk::real >(nelem),
            std::vector< tk::real >(nelem) }};

  std::size_t ncomp = U.nprop()/rdof;

  for (inciter::ncomp_t c=0; c<ncomp; ++c)
  {
    for (std::size_t e=0; e<nelem; ++e)
    {
      WENOLimiting(U, esuel, e, c, rdof, cweight, limU);
    }

    auto mark = c*rdof;

    for (std::size_t e=0; e<nelem; ++e)
    {
      U(e, mark+1) = limU[0][e];
      U(e, mark+2) = limU[1][e];
      U(e, mark+3) = limU[2][e];
    }
  }
}

void
Superbee_P1( const std::vector< int >& esuel,
             const std::vector< std::size_t >& inpoel,
             const std::vector< std::size_t >& ndofel,
             const tk::UnsMesh::Coords& coord,
             tk::Fields& U )
// *****************************************************************************
//  Superbee limiter for DGP1
//! \param[in] esuel Elements surrounding elements
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] coord Array of nodal coordinates
//! \param[in,out] U High-order solution vector which gets limited
//! \details This Superbee function should be called for transport and compflow
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
      auto phi = SuperbeeLimiting(U, esuel, inpoel, coord, e, ndof, rdof,
                   dof_el, ncomp, beta_lim);

      // apply limiter function
      for (inciter::ncomp_t c=0; c<ncomp; ++c)
      {
        auto mark = c*rdof;
        U(e, mark+1) = phi[c] * U(e, mark+1);
        U(e, mark+2) = phi[c] * U(e, mark+2);
        U(e, mark+3) = phi[c] * U(e, mark+3);
      }
    }
  }
}

void
SuperbeeMultiMat_P1(
  const std::vector< int >& esuel,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t system,
  const tk::UnsMesh::Coords& coord,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat )
// *****************************************************************************
//  Superbee limiter for multi-material DGP1
//! \param[in] esuel Elements surrounding elements
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] system Index for equation systems
//! \param[in] coord Array of nodal coordinates
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in,out] P High-order vector of primitives which gets limited
//! \param[in] nmat Number of materials in this PDE system
//! \details This Superbee function should be called for multimat
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  const auto intsharp = inciter::g_inputdeck.get< tag::param, tag::multimat,
    tag::intsharp >()[system];
  std::size_t ncomp = U.nprop()/rdof;
  std::size_t nprim = P.nprop()/rdof;

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
      // limit conserved quantities
      auto phic = SuperbeeLimiting(U, esuel, inpoel, coord, e, ndof, rdof,
                    dof_el, ncomp, beta_lim);
      // limit primitive quantities
      auto phip = SuperbeeLimiting(P, esuel, inpoel, coord, e, ndof, rdof,
                    dof_el, nprim, beta_lim);

      std::vector< tk::real > phic_p2;
      std::vector< std::vector< tk::real > > unk, prim;
      if(ndof > 1)
        BoundPreservingLimiting(nmat, ndof, e, inpoel, coord, U, phic,
          phic_p2);

      // limits under which compression is to be performed
      std::vector< std::size_t > matInt(nmat, 0);
      std::vector< tk::real > alAvg(nmat, 0.0);
      for (std::size_t k=0; k<nmat; ++k)
        alAvg[k] = U(e, volfracDofIdx(nmat,k,rdof,0));
      auto intInd = interfaceIndicator(nmat, alAvg, matInt);
      if ((intsharp > 0) && intInd)
      {
        for (std::size_t k=0; k<nmat; ++k)
        {
          if (matInt[k])
            phic[volfracIdx(nmat,k)] = 1.0;
        }
      }
      else
      {
        if (!g_inputdeck.get< tag::discr, tag::accuracy_test >())
          consistentMultiMatLimiting_P1(nmat, rdof, e, U, P, phic,
            phic_p2);
      }

      // apply limiter function
      for (inciter::ncomp_t c=0; c<ncomp; ++c)
      {
        auto mark = c*rdof;
        U(e, mark+1) = phic[c] * U(e, mark+1);
        U(e, mark+2) = phic[c] * U(e, mark+2);
        U(e, mark+3) = phic[c] * U(e, mark+3);
      }
      for (inciter::ncomp_t c=0; c<nprim; ++c)
      {
        auto mark = c*rdof;
        P(e, mark+1) = phip[c] * P(e, mark+1);
        P(e, mark+2) = phip[c] * P(e, mark+2);
        P(e, mark+3) = phip[c] * P(e, mark+3);
      }
    }
  }
}

void
VertexBasedTransport_P1(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  std::size_t system,
  const tk::UnsMesh::Coords& coord,
  tk::Fields& U )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for transport DGP1
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] system Index for equation systems
//! \param[in] coord Array of nodal coordinates
//! \param[in,out] U High-order solution vector which gets limited
//! \details This vertex-based limiter function should be called for transport.
//!   For details see: Kuzmin, D. (2010). A vertex-based hierarchical slope
//!   limiter for p-adaptive discontinuous Galerkin methods. Journal of
//!   computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  const auto intsharp = inciter::g_inputdeck.get< tag::param, tag::transport,
    tag::intsharp >()[system];
  std::size_t ncomp = U.nprop()/rdof;

  for (std::size_t e=0; e<nelem; ++e)
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
      std::vector< tk::real > phi(ncomp, 1.0);
      // limit conserved quantities
      VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
        ncomp, phi, {0, ncomp-1});

      // limits under which compression is to be performed
      std::vector< std::size_t > matInt(ncomp, 0);
      std::vector< tk::real > alAvg(ncomp, 0.0);
      for (std::size_t k=0; k<ncomp; ++k)
        alAvg[k] = U(e,k*rdof);
      auto intInd = interfaceIndicator(ncomp, alAvg, matInt);
      if ((intsharp > 0) && intInd)
      {
        for (std::size_t k=0; k<ncomp; ++k)
        {
          if (matInt[k]) phi[k] = 1.0;
        }
      }

      // apply limiter function
      for (std::size_t c=0; c<ncomp; ++c)
      {
        auto mark = c*rdof;
        U(e, mark+1) = phi[c] * U(e, mark+1);
        U(e, mark+2) = phi[c] * U(e, mark+2);
        U(e, mark+3) = phi[c] * U(e, mark+3);
      }
    }
  }
}

void
VertexBasedCompflow_P1(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  std::size_t system,
  const std::vector< inciter::EoS_Base* >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const tk::FluxFn& flux,
  tk::Fields& U,
  std::vector< std::size_t >& shockmarker )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for single-material DGP1
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] system Index for equation systems
//! \param[in] mat_blk EOS material block
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
// //! \param[in] geoElem Element geometry array
//! \param[in] coord Array of nodal coordinates
//! \param[in] flux Riemann flux function to use
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in,out] shockmarker Shock detection marker array
//! \details This vertex-based limiter function should be called for compflow.
//!   For details see: Kuzmin, D. (2010). A vertex-based hierarchical slope
//!   limiter for p-adaptive discontinuous Galerkin methods. Journal of
//!   computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  std::size_t ncomp = U.nprop()/rdof;

  if (inciter::g_inputdeck.get< tag::discr, tag::shock_detector_coeff >()
    > 1e-6)
    MarkShockCells(nelem, 1, system, ndof, rdof, mat_blk, ndofel,
      inpoel, coord, fd, geoFace, geoElem, flux, U, U, shockmarker);

  for (std::size_t e=0; e<nelem; ++e)
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

    if (dof_el > 1 && shockmarker[e])
    {
      std::vector< tk::real > phi(ncomp, 1.0);
      // limit conserved quantities
      VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
        ncomp, phi, {0, ncomp-1});

      // apply limiter function
      for (std::size_t c=0; c<ncomp; ++c)
      {
        auto mark = c*rdof;
        U(e, mark+1) = phi[c] * U(e, mark+1);
        U(e, mark+2) = phi[c] * U(e, mark+2);
        U(e, mark+3) = phi[c] * U(e, mark+3);
      }
    }
  }
}

void
VertexBasedCompflow_P2(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  std::size_t system,
  const std::vector< inciter::EoS_Base* >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  [[maybe_unused]] const std::vector< std::size_t >& gid,
  [[maybe_unused]] const std::unordered_map< std::size_t, std::size_t >& bid,
  [[maybe_unused]] const std::vector< std::vector<tk::real> >& uNodalExtrm,
  [[maybe_unused]] const std::vector< std::vector<tk::real> >& mtInv,
  const tk::FluxFn& flux,
  tk::Fields& U,
  std::vector< std::size_t >& shockmarker )
// *****************************************************************************
//  Kuzmin's vertex-based limiter on reference element for single-material DGP2
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] system Index for equation systems
//! \param[in] mat_blk EOS material block
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
// //! \param[in] geoElem Element geometry array
//! \param[in] coord Array of nodal coordinates
//! \param[in] gid Local->global node id map
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!   global node ids (key)
//! \param[in] uNodalExtrm Chare-boundary nodal extrema for conservative
//!   variables
//! \param[in] mtInv Inverse of Taylor mass matrix
//! \param[in] flux Riemann flux function to use
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in,out] shockmarker Shock detection marker array
//! \details This vertex-based limiter function should be called for compflow.
//!   For details see: Kuzmin, D. (2010). A vertex-based hierarchical slope
//!   limiter for p-adaptive discontinuous Galerkin methods. Journal of
//!   computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  std::size_t ncomp = U.nprop()/rdof;

  if (inciter::g_inputdeck.get< tag::discr, tag::shock_detector_coeff >()
    > 1e-6)
    MarkShockCells(nelem, 1, system, ndof, rdof, mat_blk, ndofel,
      inpoel, coord, fd, geoFace, geoElem, flux, U, U, shockmarker);

  for (std::size_t e=0; e<nelem; ++e)
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

    if (dof_el > 1 && shockmarker[e])
    {
      // The vector of limiting coefficients for P1 and P2 coefficients
      std::vector< tk::real > phic_p1(ncomp, 1.0);

      // Removing 3rd order DOFs if discontinuity is detected, and applying
      // limiting to the 2nd order/P1 solution
      for (std::size_t c=0; c<ncomp; ++c) {
        for(std::size_t idof = 4; idof < rdof; idof++) {
          auto mark = c * rdof + idof;
          U(e, mark) = 0.0;
        }
      }

      // Obtain limiting coefficient for P1 coefficients
      VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
        ncomp, phic_p1, {0, ncomp-1});

      // apply limiter function to the solution with Taylor basis
      for (std::size_t c=0; c<ncomp; ++c) {
        auto mark = c * rdof;
        for(std::size_t idof=1; idof<4; idof++)
          U(e, mark+idof) = phic_p1[c] * U(e, mark+idof);
      }
    }
  }
}

void
VertexBasedMultiMat_P1(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  std::size_t system,
  const std::vector< inciter::EoS_Base* >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const tk::FluxFn& flux,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat,
  std::vector< std::size_t >& shockmarker )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for multi-material DGP1
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] system Index for equation systems
//! \param[in] mat_blk EOS material block
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
// //! \param[in] geoElem Element geometry array
//! \param[in] coord Array of nodal coordinates
//! \param[in] flux Riemann flux function to use
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in,out] P High-order vector of primitives which gets limited
//! \param[in] nmat Number of materials in this PDE system
//! \param[in,out] shockmarker Shock detection marker array
//! \details This vertex-based limiter function should be called for multimat.
//!   For details see: Kuzmin, D. (2010). A vertex-based hierarchical slope
//!   limiter for p-adaptive discontinuous Galerkin methods. Journal of
//!   computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  const auto intsharp = inciter::g_inputdeck.get< tag::param, tag::multimat,
    tag::intsharp >()[system];
  std::size_t ncomp = U.nprop()/rdof;
  std::size_t nprim = P.nprop()/rdof;

  // Evaluate the interface condition and mark the shock cells
  if (inciter::g_inputdeck.get< tag::discr, tag::shock_detector_coeff >()
    > 1e-6)
    MarkShockCells(nelem, nmat, system, ndof, rdof, mat_blk, ndofel,
      inpoel, coord, fd, geoFace, geoElem, flux, U, P, shockmarker);

  for (std::size_t e=0; e<nelem; ++e)
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
      std::vector< tk::real > phic(ncomp, 1.0);
      std::vector< tk::real > phip(nprim, 1.0);
      if(shockmarker[e]) {
        // When shockmarker is 1, there is discontinuity within the element.
        // Hence, the vertex-based limiter will be applied.

        // limit conserved quantities
        VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
          ncomp, phic, {0, ncomp-1});
        // limit primitive quantities
        VertexBasedLimiting(P, esup, inpoel, coord, e, rdof, dof_el,
          nprim, phip, {0, nprim-1});
      } else {
        // When shockmarker is 0, the volume fraction, density and energy
        // of minor material will still be limited to ensure a stable solution.
        VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
          ncomp, phic, {volfracIdx(nmat,0), volfracIdx(nmat,nmat-1)});

        for(std::size_t k=0; k<nmat; ++k) {
          if(U(e, volfracDofIdx(nmat,k,rdof,0)) < 1e-4) {
            // Vector to store the range of limited variables
            std::array< std::size_t, 2 > VarRange;

            // limit the density of minor materials
            VarRange[0] = densityIdx(nmat, k);
            VarRange[1] = VarRange[0];
            VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
              ncomp, phic, VarRange);

            // limit the energy of minor materials
            VarRange[0] = energyIdx(nmat, k);
            VarRange[1] = VarRange[0];
            VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
              ncomp, phic, VarRange);

            // limit the pressure of minor materials
            VarRange[0] = pressureIdx(nmat, k);
            VarRange[1] = VarRange[0];
            VertexBasedLimiting(P, esup, inpoel, coord, e, rdof, dof_el,
              nprim, phip, VarRange);
          }
        }
      }

      std::vector< tk::real > phic_p2, phip_p2;

      PositivityLimitingMultiMat(nmat, mat_blk, rdof, dof_el, e, inpoel,
        coord, U, P, phic, phic_p2, phip, phip_p2);

      // limits under which compression is to be performed
      std::vector< std::size_t > matInt(nmat, 0);
      std::vector< tk::real > alAvg(nmat, 0.0);
      for (std::size_t k=0; k<nmat; ++k)
        alAvg[k] = U(e, volfracDofIdx(nmat,k,rdof,0));
      auto intInd = interfaceIndicator(nmat, alAvg, matInt);
      if ((intsharp > 0) && intInd) {
        for (std::size_t k=0; k<nmat; ++k) {
          if (matInt[k]) {
            phic[volfracIdx(nmat,k)] = 1.0;
          }
        }
      }
      else {
        if(nmat > 1)
          BoundPreservingLimiting(nmat, ndof, e, inpoel, coord, U, phic,
            phic_p2);

        if (!g_inputdeck.get< tag::discr, tag::accuracy_test >())
          consistentMultiMatLimiting_P1(nmat, rdof, e, U, P, phic,
            phic_p2);
      }

      // apply limiter function
      for (std::size_t c=0; c<ncomp; ++c)
      {
        auto mark = c*rdof;
        U(e, mark+1) = phic[c] * U(e, mark+1);
        U(e, mark+2) = phic[c] * U(e, mark+2);
        U(e, mark+3) = phic[c] * U(e, mark+3);
      }
      for (std::size_t c=0; c<nprim; ++c)
      {
        auto mark = c*rdof;
        P(e, mark+1) = phip[c] * P(e, mark+1);
        P(e, mark+2) = phip[c] * P(e, mark+2);
        P(e, mark+3) = phip[c] * P(e, mark+3);
      }
    }
  }
}

void
VertexBasedMultiMat_P2(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  std::size_t system,
  const std::vector< inciter::EoS_Base* >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  [[maybe_unused]] const std::vector< std::size_t >& gid,
  [[maybe_unused]] const std::unordered_map< std::size_t, std::size_t >& bid,
  [[maybe_unused]] const std::vector< std::vector<tk::real> >& uNodalExtrm,
  [[maybe_unused]] const std::vector< std::vector<tk::real> >& pNodalExtrm,
  [[maybe_unused]] const std::vector< std::vector<tk::real> >& mtInv,
  const tk::FluxFn& flux,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat,
  std::vector< std::size_t >& shockmarker )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for multi-material DGP2
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] system Index for equation systems
//! \param[in] mat_blk EOS material block
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
// //! \param[in] geoElem Element geometry array
//! \param[in] coord Array of nodal coordinates
//! \param[in] gid Local->global node id map
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!   global node ids (key)
//! \param[in] uNodalExtrm Chare-boundary nodal extrema for conservative
//!   variables
//! \param[in] pNodalExtrm Chare-boundary nodal extrema for primitive
//!   variables
//! \param[in] mtInv Inverse of Taylor mass matrix
//! \param[in] flux Riemann flux function to use
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in,out] P High-order vector of primitives which gets limited
//! \param[in] nmat Number of materials in this PDE system
//! \param[in,out] shockmarker Shock detection marker array
//! \details This vertex-based limiter function should be called for multimat.
//!   For details see: Kuzmin, D. (2010). A vertex-based hierarchical slope
//!   limiter for p-adaptive discontinuous Galerkin methods. Journal of
//!   computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  const auto intsharp = inciter::g_inputdeck.get< tag::param, tag::multimat,
    tag::intsharp >()[system];
  std::size_t ncomp = U.nprop()/rdof;
  std::size_t nprim = P.nprop()/rdof;

  // Evaluate the interface condition and mark the shock cells
  if (inciter::g_inputdeck.get< tag::discr, tag::shock_detector_coeff >()
    > 1e-6)
    MarkShockCells(nelem, nmat, system, ndof, rdof, mat_blk, ndofel,
      inpoel, coord, fd, geoFace, geoElem, flux, U, P, shockmarker);

  for (std::size_t e=0; e<nelem; ++e)
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
      // The vector of limiting coefficients for P1
      std::vector< tk::real > phic_p1(ncomp, 1.0), phic_p2(ncomp, 1.0);
      std::vector< tk::real > phip_p1(nprim, 1.0), phip_p2(nprim, 1.0);

      if(shockmarker[e]) {
        // Removing 3rd order DOFs if discontinuity is detected, and applying
        // limiting to the 2nd order/P1 solution
        for (std::size_t c=0; c<ncomp; ++c) {
          auto mark = c * rdof;
          for(std::size_t idof = 4; idof < rdof; idof++)
            U(e, mark+idof) = 0.0;
        }
        for (std::size_t c=0; c<nprim; ++c) {
          auto mark = c * rdof;
          for(std::size_t idof = 4; idof < rdof; idof++)
            P(e, mark+idof) = 0.0;
        }

        // Obtain limiter coefficient for P1 conserved quantities
        VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
          ncomp, phic_p1, {0, ncomp-1});
        // Obtain limiter coefficient for P1 primitive quantities
        VertexBasedLimiting(P, esup, inpoel, coord, e, rdof, dof_el,
          nprim, phip_p1, {0, nprim-1});
      } else {
        // When shockmarker is 0, the volume fraction will still be limited to
        // ensure a stable solution. Since the limiting strategy for third order
        // solution will downgrade the accuracy to second order, the density,
        // energy and pressure of minor material will not be limited.
        VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
          ncomp, phic_p1, {volfracIdx(nmat,0), volfracIdx(nmat,nmat-1)});

        //for(std::size_t k=0; k<nmat; ++k) {
        //  if(U(e, volfracDofIdx(nmat,k,rdof,0)) < 1e-4) {
        //    // Vector to store the range of limited variables
        //    std::array< std::size_t, 2 > VarRange;

        //    // limit the density of minor materials
        //    VarRange[0] = densityIdx(nmat, k);
        //    VarRange[1] = VarRange[0];
        //    VertexBasedLimiting(unk, U, esup, inpoel, coord, e, rdof, dof_el,
        //      ncomp, phic_p1, VarRange);

        //    // limit the pressure of minor materials
        //    VarRange[0] = pressureIdx(nmat, k);
        //    VarRange[1] = VarRange[0];
        //    VertexBasedLimiting(prim, P, esup, inpoel, coord, e, rdof, dof_el,
        //      nprim, phip_p1, VarRange);
        //  }
        //}
      }

      PositivityLimitingMultiMat(nmat, mat_blk, ndof, dof_el, e, inpoel,
          coord, U, P, phic_p1, phic_p2, phip_p1, phic_p2);

      // limits under which compression is to be performed
      std::vector< std::size_t > matInt(nmat, 0);
      std::vector< tk::real > alAvg(nmat, 0.0);
      for (std::size_t k=0; k<nmat; ++k)
        alAvg[k] = U(e, volfracDofIdx(nmat,k,rdof,0));
      auto intInd = interfaceIndicator(nmat, alAvg, matInt);
      if ((intsharp > 0) && intInd) {
        for (std::size_t k=0; k<nmat; ++k) {
          if (matInt[k]) {
            phic_p1[volfracIdx(nmat,k)] = 1.0;
          }
        }
      }
      else {
        if(nmat > 1)
          BoundPreservingLimiting(nmat, ndof, e, inpoel, coord, U,
            phic_p1, phic_p2);

        if (!g_inputdeck.get< tag::discr, tag::accuracy_test >())
          consistentMultiMatLimiting_P1(nmat, rdof, e, U, P, phic_p1,
            phic_p2);
      }

      // apply limiing coefficient
      for (std::size_t c=0; c<ncomp; ++c)
      {
        auto mark = c * rdof;
        for(std::size_t idof=1; idof<4; idof++)
          U(e, mark+idof) = phic_p1[c] * U(e, mark+idof);
        for(std::size_t idof=4; idof<rdof; idof++)
          U(e, mark+idof) = phic_p2[c] * U(e, mark+idof);
      }
      for (std::size_t c=0; c<nprim; ++c)
      {
        auto mark = c * rdof;
        for(std::size_t idof=1; idof<4; idof++)
          P(e, mark+idof) = phip_p1[c] * P(e, mark+idof);
        for(std::size_t idof=4; idof<rdof; idof++)
          P(e, mark+idof) = phip_p2[c] * P(e, mark+idof);
      }
    }
  }
}

void
VertexBasedMultiMat_FV(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  std::size_t nelem,
  std::size_t system,
  const std::vector< inciter::EoS_Base* >& mat_blk,
  const tk::UnsMesh::Coords& coord,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for multi-material FV
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] nelem Number of elements
//! \param[in] system Index for equation systems
//! \param[in] coord Array of nodal coordinates
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in,out] P High-order vector of primitives which gets limited
//! \param[in] nmat Number of materials in this PDE system
//! \details This vertex-based limiter function should be called for multimat.
//!   For details see: Kuzmin, D. (2010). A vertex-based hierarchical slope
//!   limiter for p-adaptive discontinuous Galerkin methods. Journal of
//!   computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();
  const auto intsharp = inciter::g_inputdeck.get< tag::param, tag::multimat,
    tag::intsharp >()[system];
  std::size_t ncomp = U.nprop()/rdof;
  std::size_t nprim = P.nprop()/rdof;

  for (std::size_t e=0; e<nelem; ++e)
  {
    std::vector< tk::real > phic(ncomp, 1.0);
    std::vector< tk::real > phip(nprim, 1.0);
    // limit conserved quantities
    VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, rdof, ncomp,
      phic, {0, ncomp-1});
    // limit primitive quantities
    VertexBasedLimiting(P, esup, inpoel, coord, e, rdof, rdof, nprim,
      phip, {0, nprim-1});

    // limits under which compression is to be performed
    std::vector< std::size_t > matInt(nmat, 0);
    std::vector< tk::real > alAvg(nmat, 0.0);
    for (std::size_t k=0; k<nmat; ++k)
      alAvg[k] = U(e, volfracDofIdx(nmat,k,rdof,0));
    auto intInd = interfaceIndicator(nmat, alAvg, matInt);
    if ((intsharp > 0) && intInd)
    {
      for (std::size_t k=0; k<nmat; ++k)
      {
        if (matInt[k])
          phic[volfracIdx(nmat,k)] = 1.0;
      }
    }
    else
    {
      if (!g_inputdeck.get< tag::discr, tag::accuracy_test >())
        consistentMultiMatLimiting_P1(nmat, rdof, e, U, P, phic, phip);
    }

    // apply positivity preserving limiter
    std::vector< tk::real > phic_p2, phip_p2;
    PositivityLimitingMultiMat(nmat, mat_blk, rdof, rdof, e, inpoel,
      coord, U, P, phic, phic_p2, phip, phip_p2);

    // apply limiter function
    for (std::size_t c=0; c<ncomp; ++c)
    {
      auto mark = c*rdof;
      U(e, mark+1) = phic[c] * U(e, mark+1);
      U(e, mark+2) = phic[c] * U(e, mark+2);
      U(e, mark+3) = phic[c] * U(e, mark+3);
    }
    for (std::size_t c=0; c<nprim; ++c)
    {
      auto mark = c*rdof;
      P(e, mark+1) = phip[c] * P(e, mark+1);
      P(e, mark+2) = phip[c] * P(e, mark+2);
      P(e, mark+3) = phip[c] * P(e, mark+3);
    }
  }
}

void
WENOLimiting( const tk::Fields& U,
              const std::vector< int >& esuel,
              std::size_t e,
              inciter::ncomp_t c,
              std::size_t rdof,
              tk::real cweight,
              std::array< std::vector< tk::real >, 3 >& limU )
// *****************************************************************************
//  WENO limiter function calculation for P1 dofs
//! \param[in] U High-order solution vector which is to be limited
//! \param[in] esuel Elements surrounding elements
//! \param[in] e Id of element whose solution is to be limited
//! \param[in] c Index of component which is to be limited
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] cweight Weight of the central stencil
//! \param[in,out] limU Limited gradients of component c
// *****************************************************************************
{
  std::array< std::array< tk::real, 3 >, 5 > gradu;
  std::array< tk::real, 5 > wtStencil, osc, wtDof;

  auto mark = c*rdof;

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
  gradu[0][0] = U(e, mark+1);
  gradu[0][1] = U(e, mark+2);
  gradu[0][2] = U(e, mark+3);
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
    gradu[is][0] = U(n, mark+1);
    gradu[is][1] = U(n, mark+2);
    gradu[is][2] = U(n, mark+3);
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

std::vector< tk::real >
SuperbeeLimiting( const tk::Fields& U,
                  const std::vector< int >& esuel,
                  const std::vector< std::size_t >& inpoel,
                  const tk::UnsMesh::Coords& coord,
                  std::size_t e,
                  std::size_t ndof,
                  std::size_t rdof,
                  std::size_t dof_el,
                  inciter:: ncomp_t ncomp,
                  tk::real beta_lim )
// *****************************************************************************
//  Superbee limiter function calculation for P1 dofs
//! \param[in] U High-order solution vector which is to be limited
//! \param[in] esuel Elements surrounding elements
//! \param[in] inpoel Element connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] e Id of element whose solution is to be limited
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] dof_el Local number of degrees of freedom
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] beta_lim Parameter which is equal to 2 for Superbee and 1 for
//!   minmod limiter
//! \return phi Limiter function for solution in element e
// *****************************************************************************
{
  // Superbee is a TVD limiter, which uses min-max bounds that the
  // high-order solution should satisfy, to ensure TVD properties. For a
  // high-order method like DG, this involves the following steps:
  // 1. Find min-max bounds in the immediate neighborhood of cell.
  // 2. Calculate the Superbee function for all the points where solution
  //    needs to be reconstructed to (all quadrature points). From these,
  //    use the minimum value of the limiter function.

  std::vector< tk::real > uMin(ncomp, 0.0), uMax(ncomp, 0.0);

  for (inciter::ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*rdof;
    uMin[c] = U(e, mark);
    uMax[c] = U(e, mark);
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
      uMin[c] = std::min(uMin[c], U(n, mark));
      uMax[c] = std::max(uMax[c], U(n, mark));
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

      auto state =
        tk::eval_state(ncomp, rdof, dof_el, e, U, B_l, {0, ncomp-1});

      Assert( state.size() == ncomp, "Size mismatch" );

      // compute the limiter function
      for (inciter::ncomp_t c=0; c<ncomp; ++c)
      {
        auto phi_gp = 1.0;
        auto mark = c*rdof;
        auto uNeg = state[c] - U(e, mark);
        if (uNeg > 1.0e-14)
        {
          uNeg = std::max(uNeg, 1.0e-08);
          phi_gp = std::min( 1.0, (uMax[c]-U(e, mark))/(2.0*uNeg) );
        }
        else if (uNeg < -1.0e-14)
        {
          uNeg = std::min(uNeg, -1.0e-08);
          phi_gp = std::min( 1.0, (uMin[c]-U(e, mark))/(2.0*uNeg) );
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

  return phi;
}

void
VertexBasedLimiting(
  const tk::Fields& U,
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  std::size_t e,
  std::size_t rdof,
  std::size_t dof_el,
  std::size_t ncomp,
  std::vector< tk::real >& phi,
  const std::array< std::size_t, 2 >& VarRange )
// *****************************************************************************
//  Kuzmin's vertex-based limiter function calculation for P1 dofs
//! \param[in] U High-order solution vector which is to be limited
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] e Id of element whose solution is to be limited
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] dof_el Local number of degrees of freedom
//! \param[in] ncomp Number of scalar components in this PDE system
//! \return phi Limiter function for solution in element e
// *****************************************************************************
{
  // Kuzmin's vertex-based TVD limiter uses min-max bounds that the
  // high-order solution should satisfy, to ensure TVD properties. For a
  // high-order method like DG, this involves the following steps:
  // 1. Find min-max bounds in the nodal-neighborhood of cell.
  // 2. Calculate the limiter function (Superbee) for all the vertices of cell.
  //    From these, use the minimum value of the limiter function.

  // Prepare for calculating Basis functions
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

  std::vector< tk::real > uMin(VarRange[1]-VarRange[0]+1, 0.0),
                          uMax(VarRange[1]-VarRange[0]+1, 0.0);

  // loop over all nodes of the element e
  for (std::size_t lp=0; lp<4; ++lp)
  {
    // reset min/max
    for (std::size_t c=VarRange[0]; c<=VarRange[1]; ++c)
    {
      auto mark = c*rdof;
      auto cmark = c-VarRange[0];
      uMin[cmark] = U(e, mark);
      uMax[cmark] = U(e, mark);
    }
    auto p = inpoel[4*e+lp];
    const auto& pesup = tk::cref_find(esup, p);

    // ----- Step-1: find min/max in the neighborhood of node p
    // loop over all the internal elements surrounding this node p
    for (auto er : pesup)
    {
      for (std::size_t c=VarRange[0]; c<=VarRange[1]; ++c)
      {
        auto mark = c*rdof;
        auto cmark = c-VarRange[0];
        uMin[cmark] = std::min(uMin[cmark], U(er, mark));
        uMax[cmark] = std::max(uMax[cmark], U(er, mark));
      }
    }

    // ----- Step-2: compute the limiter function at this node
    // find high-order solution
    std::vector< tk::real > state;
    std::array< tk::real, 3 > gp{cx[p], cy[p], cz[p]};
    auto B_p = tk::eval_basis( rdof,
          tk::Jacobian( coordel[0], gp, coordel[2], coordel[3] ) / detT,
          tk::Jacobian( coordel[0], coordel[1], gp, coordel[3] ) / detT,
          tk::Jacobian( coordel[0], coordel[1], coordel[2], gp ) / detT );
    state = tk::eval_state(ncomp, rdof, dof_el, e, U, B_p, VarRange);

    Assert( state.size() == ncomp, "Size mismatch" );

    // compute the limiter function
    for (std::size_t c=VarRange[0]; c<=VarRange[1]; ++c)
    {
      auto phi_gp = 1.0;
      auto mark = c*rdof;
      auto uNeg = state[c] - U(e, mark);
      auto uref = std::max(std::fabs(U(e,mark)), 1e-14);
      auto cmark = c - VarRange[0];
      if (uNeg > 1.0e-06*uref)
      {
        phi_gp = std::min( 1.0, (uMax[cmark]-U(e, mark))/uNeg );
      }
      else if (uNeg < -1.0e-06*uref)
      {
        phi_gp = std::min( 1.0, (uMin[cmark]-U(e, mark))/uNeg );
      }
      else
      {
        phi_gp = 1.0;
      }

    // ----- Step-3: take the minimum of the nodal-limiter functions
      phi[c] = std::min( phi[c], phi_gp );
    }
  }
}

void
VertexBasedLimiting_P2( const std::vector< std::vector< tk::real > >& unk,
  const tk::Fields& U,
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  std::size_t e,
  std::size_t rdof,
  [[maybe_unused]] std::size_t dof_el,
  std::size_t ncomp,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& NodalExtrm,
  const std::array< std::size_t, 2 >& VarRange,
  std::vector< tk::real >& phi )
// *****************************************************************************
//  Kuzmin's vertex-based limiter function calculation for P2 dofs
//! \param[in] U High-order solution vector which is to be limited
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] e Id of element whose solution is to be limited
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] dof_el Local number of degrees of freedom
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] gid Local->global node id map
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!   global node ids (key)
//! \param[in] NodalExtrm Chare-boundary nodal extrema
//! \param[in] VarRange The range of limited variables
//! \param[out] phi Limiter function for solution in element e
//! \details This function limits the P2 dofs of P2 solution in a hierachical
//!   way to P1 dof limiting. Here we treat the first order derivatives the same
//!   way as cell average while second order derivatives represent the gradients
//!   to be limited in the P1 limiting procedure.
// *****************************************************************************
{
  const auto nelem = inpoel.size() / 4;

  std::vector< std::vector< tk::real > > uMin, uMax;
  uMin.resize( VarRange[1]-VarRange[0]+1, std::vector<tk::real>(3, 0.0) );
  uMax.resize( VarRange[1]-VarRange[0]+1, std::vector<tk::real>(3, 0.0) );

  // The coordinates of centroid in the reference domain
  std::array< std::vector< tk::real >, 3 > center;
  center[0].resize(1, 0.25);
  center[1].resize(1, 0.25);
  center[2].resize(1, 0.25);

  std::array< std::array< tk::real, 4 >, 3 > cnodes{{
    {{0, 1, 0, 0}},
    {{0, 0, 1, 0}},
    {{0, 0, 0, 1}} }};

  // loop over all nodes of the element e
  for (std::size_t lp=0; lp<4; ++lp)
  {
    // Find the max/min first-order derivatives for internal element
    for (std::size_t c=VarRange[0]; c<=VarRange[1]; ++c)
    {
      auto mark = c - VarRange[0];
      for (std::size_t idir=1; idir < 4; ++idir)
      {
        uMin[mark][idir-1] = unk[c][idir];
        uMax[mark][idir-1] = unk[c][idir];
      }
    }

    auto p = inpoel[4*e+lp];
    const auto& pesup = tk::cref_find(esup, p);

    // Step-1: find min/max first order derivative at the centroid in the
    // neighborhood of node p
    for (auto er : pesup)
    {
      if(er < nelem)      // If this is internal element
      {
        // Compute the derivatives of basis function in the reference domain
        auto dBdxi_er = tk::eval_dBdxi(rdof,
          {{center[0][0], center[1][0], center[2][0]}});

        for (std::size_t c=VarRange[0]; c<=VarRange[1]; ++c)
        {
          auto mark = c*rdof;
          auto cmark = c-VarRange[0];
          for (std::size_t idir = 0; idir < 3; ++idir)
          {
            // The first order derivative at the centroid of element er
            tk::real slope_er(0.0);
            for(std::size_t idof = 1; idof < rdof; idof++)
              slope_er += U(er, mark+idof) * dBdxi_er[idir][idof];

            uMin[cmark][idir] = std::min(uMin[cmark][idir], slope_er);
            uMax[cmark][idir] = std::max(uMax[cmark][idir], slope_er);

          }
        }
      }
    }
    // If node p is the chare-boundary node, find min/max by comparing with
    // the chare-boundary nodal extrema from vector NodalExtrm
    auto gip = bid.find( gid[p] );
    if(gip != end(bid))
    {
      auto ndof_NodalExtrm = NodalExtrm[0].size() / (ncomp * 2);
      for (std::size_t c=VarRange[0]; c<=VarRange[1]; ++c)
      {
        auto cmark = c-VarRange[0];
        for (std::size_t idir = 0; idir < 3; idir++)
        {
          auto max_mark = 2*c*ndof_NodalExtrm + 2*idir;
          auto min_mark = max_mark + 1;
          const auto& ex = NodalExtrm[gip->second];
          uMax[cmark][idir] = std::max(ex[max_mark], uMax[cmark][idir]);
          uMin[cmark][idir] = std::min(ex[min_mark], uMin[cmark][idir]);
        }
      }
    }

    //Step-2: compute the limiter function at this node
    std::array< tk::real, 3 > node{cnodes[0][lp], cnodes[1][lp], cnodes[2][lp]};

    // find high-order solution
    std::vector< std::array< tk::real, 3 > > state;
    state.resize(VarRange[1]-VarRange[0]+1);

    for (std::size_t c=VarRange[0]; c<=VarRange[1]; ++c)
    {
      auto cmark = c-VarRange[0];
      auto dx = node[0] - center[0][0];
      auto dy = node[1] - center[1][0];
      auto dz = node[2] - center[2][0];

      state[cmark][0] = unk[c][1] + unk[c][4]*dx + unk[c][7]*dy + unk[c][8]*dz;
      state[cmark][1] = unk[c][2] + unk[c][5]*dy + unk[c][7]*dx + unk[c][9]*dz;
      state[cmark][2] = unk[c][3] + unk[c][6]*dz + unk[c][8]*dx + unk[c][9]*dy;
    }

    // compute the limiter function
    for (std::size_t c=VarRange[0]; c<=VarRange[1]; ++c)
    {
      tk::real phi_dir(1.0);
      auto cmark = c-VarRange[0];
      for (std::size_t idir = 1; idir <= 3; ++idir)
      {
        phi_dir = 1.0;
        auto uNeg = state[cmark][idir-1] - unk[c][idir];
        auto uref = std::max(std::fabs(unk[c][idir]), 1e-14);
        if (uNeg > 1.0e-6*uref)
        {
          phi_dir =
            std::min( 1.0, ( uMax[cmark][idir-1] - unk[c][idir])/uNeg );
        }
        else if (uNeg < -1.0e-6*uref)
        {
          phi_dir =
            std::min( 1.0, ( uMin[cmark][idir-1] - unk[c][idir])/uNeg );
        }
        else
        {
          phi_dir = 1.0;
        }

        phi[c] = std::min( phi[c], phi_dir );
      }
    }
  }
}

void consistentMultiMatLimiting_P1(
  std::size_t nmat,
  std::size_t rdof,
  std::size_t e,
  tk::Fields& U,
  [[maybe_unused]] tk::Fields& P,
  std::vector< tk::real >& phic_p1,
  std::vector< tk::real >& phic_p2 )
// *****************************************************************************
//  Consistent limiter modifications for conservative variables
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] e Element being checked for consistency
//! \param[in] U Vector of conservative variables
//! \param[in] P Vector of primitive variables
//! \param[in,out] phic_p1 Vector of limiter functions for P1 dofs of the
//!   conserved quantities
//! \param[in,out] phip_p2 Vector of limiter functions for P2 dofs of the
//!   conserved quantities
// *****************************************************************************
{
  // find the limiter-function for volume-fractions
  auto phi_al_p1(1.0), phi_al_p2(1.0), almax(0.0), dalmax(0.0);
  //std::size_t nmax(0);
  for (std::size_t k=0; k<nmat; ++k)
  {
    phi_al_p1 = std::min( phi_al_p1, phic_p1[volfracIdx(nmat, k)] );
    if(rdof > 4)
      phi_al_p2 = std::min( phi_al_p2, phic_p2[volfracIdx(nmat, k)] );
    if (almax < U(e,volfracDofIdx(nmat, k, rdof, 0)))
    {
      //nmax = k;
      almax = U(e,volfracDofIdx(nmat, k, rdof, 0));
    }
    tk::real dmax(0.0);
    dmax = std::max(
             std::max(
               std::abs(U(e,volfracDofIdx(nmat, k, rdof, 1))),
               std::abs(U(e,volfracDofIdx(nmat, k, rdof, 2))) ),
               std::abs(U(e,volfracDofIdx(nmat, k, rdof, 3))) );
    dalmax = std::max( dalmax, dmax );
  }

  auto al_band = 1e-4;

  //phi_al = phic[nmax];

  // determine if cell is a material-interface cell based on ad-hoc tolerances.
  // if interface-cell, then modify high-order dofs of conserved unknowns
  // consistently and use same limiter for all equations.
  // Slopes of solution variables \alpha_k \rho_k and \alpha_k \rho_k E_k need
  // to be modified in interface cells, such that slopes in the \rho_k and
  // \rho_k E_k part are ignored and only slopes in \alpha_k are considered.
  // Ideally, we would like to not do this, but this is a necessity to avoid
  // limiter-limiter interactions in multiphase CFD (see "K.-M. Shyue, F. Xiao,
  // An Eulerian interface sharpening algorithm for compressible two-phase flow:
  // the algebraic THINC approach, Journal of Computational Physics 268, 2014,
  // 326–354. doi:10.1016/j.jcp.2014.03.010." and "A. Chiapolino, R. Saurel,
  // B. Nkonga, Sharpening diffuse interfaces with compressible fluids on
  // unstructured meshes, Journal of Computational Physics 340 (2017) 389–417.
  // doi:10.1016/j.jcp.2017.03.042."). This approximation should be applied in
  // as narrow a band of interface-cells as possible. The following if-test
  // defines this band of interface-cells. This tests checks the value of the
  // maximum volume-fraction in the cell (almax) and the maximum change in
  // volume-fraction in the cell (dalmax, calculated from second-order DOFs),
  // to determine the band of interface-cells where the aforementioned fix needs
  // to be applied. This if-test says that, the fix is applied when the change
  // in volume-fraction across a cell is greater than 0.1, *and* the
  // volume-fraction is between 0.1 and 0.9.
  if ( //dalmax > al_band &&
       (almax > al_band && almax < (1.0-al_band)) )
  {
    // 1. consistent high-order dofs
    for (std::size_t k=0; k<nmat; ++k)
    {
      auto alk =
        std::max( 1.0e-14, U(e,volfracDofIdx(nmat, k, rdof, 0)) );
      auto rhok = U(e,densityDofIdx(nmat, k, rdof, 0)) / alk;
      auto rhoE = U(e,energyDofIdx(nmat, k, rdof, 0)) / alk;
      for (std::size_t idof=1; idof<rdof; ++idof)
      {
          U(e,densityDofIdx(nmat, k, rdof, idof)) = rhok *
            U(e,volfracDofIdx(nmat, k, rdof, idof));
          U(e,energyDofIdx(nmat, k, rdof, idof)) = rhoE *
            U(e,volfracDofIdx(nmat, k, rdof, idof));
      }
    }

    // 2. same limiter for all volume-fractions and densities
    for (std::size_t k=0; k<nmat; ++k)
    {
      phic_p1[volfracIdx(nmat, k)] = phi_al_p1;
      phic_p1[densityIdx(nmat, k)] = phi_al_p1;
      phic_p1[energyIdx(nmat, k)] = phi_al_p1;
    }
    if(rdof > 4)
    {
      for (std::size_t k=0; k<nmat; ++k)
      {
        phic_p2[volfracIdx(nmat, k)] = phi_al_p2;
        phic_p2[densityIdx(nmat, k)] = phi_al_p2;
        phic_p2[energyIdx(nmat, k)] = phi_al_p2;
      }
    }
  }
  else
  {
    // same limiter for all volume-fractions
    for (std::size_t k=volfracIdx(nmat, 0); k<volfracIdx(nmat, nmat); ++k)
      phic_p1[volfracIdx(nmat, k)] = phi_al_p1;
    if(rdof > 4)
      for (std::size_t k=volfracIdx(nmat, 0); k<volfracIdx(nmat, nmat); ++k)
        phic_p2[volfracIdx(nmat, k)] = phi_al_p2;
  }
}

void BoundPreservingLimiting( std::size_t nmat,
                              std::size_t ndof,
                              std::size_t e,
                              const std::vector< std::size_t >& inpoel,
                              const tk::UnsMesh::Coords& coord,
                              const tk::Fields& U,
                              std::vector< tk::real >& phic_p1,
                              std::vector< tk::real >& phic_p2 )
// *****************************************************************************
//  Bound preserving limiter for volume fractions when MulMat scheme is selected
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] ndof Total number of reconstructed dofs
//! \param[in] e Element being checked for consistency
//! \param[in] inpoel Element connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in,out] U Second-order solution vector which gets modified near
//!   material interfaces for consistency
//! \param[in] unk Vector of conservative variables based on Taylor basis
//! \param[in,out] phic_p1 Vector of limiter functions for P1 dofs of the
//!   conserved quantities
//! \param[in,out] phic_p2 Vector of limiter functions for P2 dofs of the
//!   conserved quantities
//! \details This bound-preserving limiter is specifically meant to enforce
//!   bounds [0,1], but it does not suppress oscillations like the other 'TVD'
//!   limiters. TVD limiters on the other hand, do not preserve such bounds. A
//!   combination of oscillation-suppressing and bound-preserving limiters can
//!   obtain a non-oscillatory and bounded solution.
// *****************************************************************************
{
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

  std::vector< tk::real > phi_bound(nmat, 1.0);

  // Compute the upper and lower bound for volume fraction
  const tk::real min = 1e-14;
  const tk::real max = 1.0 - min;

  // loop over all faces of the element e
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

    auto ng = tk::NGfa(ndof);

    // arrays for quadrature points
    std::array< std::vector< tk::real >, 2 > coordgp;
    std::vector< tk::real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    wgp.resize( ng );

    // get quadrature point weights and coordinates for triangle
    tk::GaussQuadratureTri( ng, coordgp, wgp );

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the coordinates of quadrature point at physical domain
      auto gp = tk::eval_gp( igp, coordfa, coordgp );

      //Compute the basis functions
      auto B = tk::eval_basis( ndof,
            tk::Jacobian( coordel[0], gp, coordel[2], coordel[3] ) / detT,
            tk::Jacobian( coordel[0], coordel[1], gp, coordel[3] ) / detT,
            tk::Jacobian( coordel[0], coordel[1], coordel[2], gp ) / detT );

      auto state = eval_state( U.nprop()/ndof, ndof, ndof, e, U, B,
        {0, U.nprop()/ndof-1} );

      for(std::size_t imat = 0; imat < nmat; imat++)
      {
        auto phi = BoundPreservingLimitingFunction( min, max,
          state[volfracIdx(nmat, imat)],
          U(e,volfracDofIdx(nmat, imat, ndof, 0)) );
        phi_bound[imat] = std::min( phi_bound[imat], phi );
      }
    }
  }

  // If DG(P2), the bound-preserving limiter should also be applied to the gauss
  // point within the element
  if(ndof > 4)
  {
    auto ng = tk::NGvol(ndof);

    // arrays for quadrature points
    std::array< std::vector< tk::real >, 3 > coordgp;
    std::vector< tk::real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    coordgp[2].resize( ng );
    wgp.resize( ng );

    tk::GaussQuadratureTet( ng, coordgp, wgp );

    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the basis function
      auto B = tk::eval_basis( ndof, coordgp[0][igp], coordgp[1][igp],
        coordgp[2][igp] );

      auto state = tk::eval_state(U.nprop()/ndof, ndof, ndof, e, U, B,
        {0, U.nprop()/ndof-1} );

      for(std::size_t imat = 0; imat < nmat; imat++)
      {
        auto phi = BoundPreservingLimitingFunction(min, max,
          state[volfracIdx(nmat, imat)],
          U(e,volfracDofIdx(nmat, imat, ndof, 0)) );
        phi_bound[imat] = std::min( phi_bound[imat], phi );
      }
    }
  }

  for(std::size_t k=volfracIdx(nmat, 0); k<volfracIdx(nmat, nmat); k++)
    phic_p1[k] = std::min(phi_bound[k], phic_p1[k]);
  if(ndof > 4)
    for(std::size_t k=volfracIdx(nmat, 0); k<volfracIdx(nmat, nmat); k++)
      phic_p2[k] = std::min(phi_bound[k], phic_p2[k]);
}

tk::real
BoundPreservingLimitingFunction( const tk::real min,
                                 const tk::real max,
                                 const tk::real al_gp,
                                 const tk::real al_avg )
// *****************************************************************************
//  Bound-preserving limiter function for the volume fractions
//! \param[in] min Minimum bound for volume fraction
//! \param[in] max Maximum bound for volume fraction
//! \param[in] al_gp Volume fraction at the quadrature point
//! \param[in] al_avg Cell-average volume fraction
//! \return The limiting coefficient from the bound-preserving limiter function
// *****************************************************************************
{
  tk::real phi(1.0), al_diff(0.0);
  al_diff = al_gp - al_avg;
  if(al_gp > max && fabs(al_diff) > 1e-15)
    phi = std::fabs( (max - al_avg) / al_diff );
  else if(al_gp < min && fabs(al_diff) > 1e-15)
    phi = std::fabs( (min - al_avg) / al_diff );
  return phi;
}

void PositivityLimitingMultiMat( std::size_t nmat,
                               const std::vector< inciter::EoS_Base* >& mat_blk,
                                 std::size_t rdof,
                                 std::size_t ndof_el,
                                 std::size_t e,
                                 const std::vector< std::size_t >& inpoel,
                                 const tk::UnsMesh::Coords& coord,
                                 const tk::Fields& U,
                                 const tk::Fields& P,
                                 std::vector< tk::real >& phic_p1,
                                 std::vector< tk::real >& phic_p2,
                                 std::vector< tk::real >& phip_p1,
                                 std::vector< tk::real >& phip_p2 )
// *****************************************************************************
//  Positivity preserving limiter for multi-material solver
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] ndof_el Number of dofs for element e
//! \param[in] e Element being checked for consistency
//! \param[in] inpoel Element connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] U Vector of conservative variables
//! \param[in] P Vector of primitive variables
//! \param[in,out] phic_p1 Vector of limiter functions for P1 dofs of the
//!   conserved quantities
//! \param[in,out] phic_p2 Vector of limiter functions for P2 dofs of the
//!   conserved quantities
//! \param[in,out] phip_p1 Vector of limiter functions for P1 dofs of the
//!   primitive quantities
//! \param[in,out] phip_p2 Vector of limiter functions for P2 dofs of the
//!   primitive quantities
// *****************************************************************************
{
  const auto ncomp = U.nprop() / rdof;
  const auto nprim = P.nprop() / rdof;

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

  std::vector< tk::real > phic_bound(ncomp, 1.0);
  std::vector< tk::real > phip_bound(nprim, 1.0);

  const tk::real min = 1e-15;

  for (std::size_t lf=0; lf<4; ++lf)
  {
    std::array< std::size_t, 3 > inpofa_l {{ inpoel[4*e+tk::lpofa[lf][0]],
                                             inpoel[4*e+tk::lpofa[lf][1]],
                                             inpoel[4*e+tk::lpofa[lf][2]] }};

    std::array< std::array< tk::real, 3>, 3 > coordfa {{
      {{ cx[ inpofa_l[0] ], cy[ inpofa_l[0] ], cz[ inpofa_l[0] ] }},
      {{ cx[ inpofa_l[1] ], cy[ inpofa_l[1] ], cz[ inpofa_l[1] ] }},
      {{ cx[ inpofa_l[2] ], cy[ inpofa_l[2] ], cz[ inpofa_l[2] ] }} }};

    auto ng = tk::NGfa(ndof_el);

    std::array< std::vector< tk::real >, 2 > coordgp;
    std::vector< tk::real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    wgp.resize( ng );

    tk::GaussQuadratureTri( ng, coordgp, wgp );

    for (std::size_t igp=0; igp<ng; ++igp)
    {
      auto gp = tk::eval_gp( igp, coordfa, coordgp );
      auto B = tk::eval_basis( ndof_el,
            tk::Jacobian( coordel[0], gp, coordel[2], coordel[3] ) / detT,
            tk::Jacobian( coordel[0], coordel[1], gp, coordel[3] ) / detT,
            tk::Jacobian( coordel[0], coordel[1], coordel[2], gp ) / detT );

      auto state = eval_state(ncomp, rdof, ndof_el, e, U, B,
        {0, ncomp-1});
      auto sprim = eval_state(nprim, rdof, ndof_el, e, P, B,
        {0, nprim-1});

      for(std::size_t imat = 0; imat < nmat; imat++)
      {
        tk::real phi_rho(1.0), phi_rhoe(1.0), phi_pre(1.0);
        // Evaluate the limiting coefficient for material density
        auto rho = state[densityIdx(nmat, imat)];
        auto rho_avg = U(e, densityDofIdx(nmat, imat, rdof, 0));
        phi_rho = PositivityLimiting(min, rho, rho_avg);
        phic_bound[densityIdx(nmat, imat)] =
          std::min(phic_bound[densityIdx(nmat, imat)], phi_rho);
        // Evaluate the limiting coefficient for material energy
        auto rhoe = state[energyIdx(nmat, imat)];
        auto rhoe_avg = U(e, energyDofIdx(nmat, imat, rdof, 0));
        phi_rhoe = PositivityLimiting(min, rhoe, rhoe_avg);
        phic_bound[energyIdx(nmat, imat)] =
          std::min(phic_bound[energyIdx(nmat, imat)], phi_rhoe);
        // Evaluate the limiting coefficient for material pressure
        auto min_pre = mat_blk[imat]->min_eff_pressure(min);
        auto pre = sprim[pressureIdx(nmat, imat)];
        auto pre_avg = P(e, pressureDofIdx(nmat, imat, rdof, 0));
        phi_pre = PositivityLimiting(min_pre, pre, pre_avg);
        phip_bound[pressureIdx(nmat, imat)] =
          std::min(phip_bound[pressureIdx(nmat, imat)], phi_pre);
      }
    }
  }

  if(ndof_el > 4)
  {
    auto ng = tk::NGvol(ndof_el);
    std::array< std::vector< tk::real >, 3 > coordgp;
    std::vector< tk::real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    coordgp[2].resize( ng );
    wgp.resize( ng );

    tk::GaussQuadratureTet( ng, coordgp, wgp );

    for (std::size_t igp=0; igp<ng; ++igp)
    {
      auto B = tk::eval_basis( ndof_el, coordgp[0][igp], coordgp[1][igp],
        coordgp[2][igp] );

      auto state = eval_state(ncomp, rdof, ndof_el, e, U, B,
        {0, ncomp-1});
      auto sprim = eval_state(nprim, rdof, ndof_el, e, P, B,
        {0, nprim-1});

      for(std::size_t imat = 0; imat < nmat; imat++)
      {
        tk::real phi_rho(1.0), phi_rhoe(1.0), phi_pre(1.0);
        // Evaluate the limiting coefficient for material density
        auto rho = state[densityIdx(nmat, imat)];
        auto rho_avg = U(e, densityDofIdx(nmat, imat, rdof, 0));
        phi_rho = PositivityLimiting(min, rho, rho_avg);
        phic_bound[densityIdx(nmat, imat)] =
          std::min(phic_bound[densityIdx(nmat, imat)], phi_rho);
        // Evaluate the limiting coefficient for material energy
        auto rhoe = state[energyIdx(nmat, imat)];
        auto rhoe_avg = U(e, energyDofIdx(nmat, imat, rdof, 0));
        phi_rhoe = PositivityLimiting(min, rhoe, rhoe_avg);
        phic_bound[energyIdx(nmat, imat)] =
          std::min(phic_bound[energyIdx(nmat, imat)], phi_rhoe);
        // Evaluate the limiting coefficient for material pressure
        auto min_pre = mat_blk[imat]->min_eff_pressure(min);
        auto pre = sprim[pressureIdx(nmat, imat)];
        auto pre_avg = P(e, pressureDofIdx(nmat, imat, rdof, 0));
        phi_pre = PositivityLimiting(min_pre, pre, pre_avg);
        phip_bound[pressureIdx(nmat, imat)] =
          std::min(phip_bound[pressureIdx(nmat, imat)], phi_pre);
      }
    }
  }
  for(std::size_t icomp = volfracIdx(nmat, nmat); icomp < ncomp; icomp++)
    phic_p1[icomp] = std::min( phic_bound[icomp], phic_p1[icomp] );
  for(std::size_t icomp = pressureIdx(nmat, 0); icomp < pressureIdx(nmat, nmat);
      icomp++)
    phip_p1[icomp] = std::min( phip_bound[icomp], phip_p1[icomp] );
  if(ndof_el > 4) {
    for(std::size_t icomp = volfracIdx(nmat, nmat); icomp < ncomp; icomp++)
      phic_p2[icomp] = std::min( phic_bound[icomp], phic_p2[icomp] );
    for(std::size_t icomp = pressureIdx(nmat, 0); icomp < pressureIdx(nmat, nmat);
        icomp++)
      phip_p2[icomp] = std::min( phip_bound[icomp], phip_p2[icomp] );
  }
}

tk::real
PositivityLimiting( const tk::real min,
                    const tk::real u_gp,
                    const tk::real u_avg )
// *****************************************************************************
//  Positivity-preserving limiter function
//! \param[in] min Minimum bound for volume fraction
//! \param[in] u_gp Variable quantity at the quadrature point
//! \param[in] u_avg Cell-average variable quantitiy
//! \return The limiting coefficient from the positivity-preserving limiter
//!   function
// *****************************************************************************
{
  tk::real phi(1.0);
  tk::real diff = u_gp - u_avg;
  // Only when u_gp is less than minimum threshold and the high order
  // contribution is not zero, the limiting function will be applied
  if(u_gp < min)
    phi = std::fabs( (min - u_avg) / (diff+std::copysign(1e-15,diff)) );
  return phi;
}

bool
interfaceIndicator( std::size_t nmat,
  const std::vector< tk::real >& al,
  std::vector< std::size_t >& matInt )
// *****************************************************************************
//  Interface indicator function, which checks element for material interface
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] al Cell-averaged volume fractions
//! \param[in] matInt Array indicating which material has an interface
//! \return Boolean which indicates if the element contains a material interface
// *****************************************************************************
{
  bool intInd = false;

  // limits under which compression is to be performed
  auto al_eps = 1e-08;
  auto loLim = 2.0 * al_eps;
  auto hiLim = 1.0 - loLim;

  auto almax = 0.0;
  for (std::size_t k=0; k<nmat; ++k)
  {
    almax = std::max(almax, al[k]);
    matInt[k] = 0;
    if ((al[k] > loLim) && (al[k] < hiLim)) matInt[k] = 1;
  }

  if ((almax > loLim) && (almax < hiLim)) intInd = true;

  return intInd;
}

void MarkShockCells ( const std::size_t nelem,
                      const std::size_t nmat,
                      const std::size_t system,
                      const std::size_t ndof,
                      const std::size_t rdof,
                      const std::vector< inciter::EoS_Base* >& mat_blk,
                      const std::vector< std::size_t >& ndofel,
                      const std::vector< std::size_t >& inpoel,
                      const tk::UnsMesh::Coords& coord,
                      const inciter::FaceData& fd,
                      [[maybe_unused]] const tk::Fields& geoFace,
                      const tk::Fields& geoElem,
                      const tk::FluxFn& flux,
                      const tk::Fields& U,
                      const tk::Fields& P,
                      std::vector< std::size_t >& shockmarker )
// *****************************************************************************
//  Mark the cells that contain discontinuity according to the interface
//    condition
//! \param[in] nelem Number of elements
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] system Equation system index
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] mat_blk EOS material block
//! \param[in] ndofel Vector of local number of degrees of freedome
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] flux Flux function to use
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in, out] shockmarker Vector of the shock indicator
//! \details This function computes the discontinuity indicator based on
//!   interface conditon. It is based on the following paper:
//!   Hong L., Gianni A., Robert N. (2021) A moving discontinuous Galerkin
//!   finite element method with interface condition enforcement for
//!   compressible flows. Journal of Computational Physics,
//!   doi: https://doi.org/10.1016/j.jcp.2021.110618
// *****************************************************************************
{
  const auto coeff = g_inputdeck.get< tag::discr, tag::shock_detector_coeff >();

  std::vector< tk::real > IC(U.nunk(), 0.0);
  const auto& esuf = fd.Esuf();
  const auto& inpofa = fd.Inpofa();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  // The interface-conservation based indicator will only evaluate the flux jump
  // for the momentum equations
  std::array< std::size_t, 2 > VarRange;
  if(nmat > 1) {          // multi-material flow
    VarRange[0] = momentumIdx(nmat, 0);
    VarRange[1] = momentumIdx(nmat, 2);
  } else {                // single-material flow
    VarRange[0] = 1;
    VarRange[1] = 3;
  }

  // Loop over faces
  for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f) {
    Assert( esuf[2*f] > -1 && esuf[2*f+1] > -1, "Interior element detected "
            "as -1" );

    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    // When the number of gauss points for the left and right element are
    // different, choose the larger ng
    auto ng_l = tk::NGfa(ndofel[el]);
    auto ng_r = tk::NGfa(ndofel[er]);

    auto ng = std::max( ng_l, ng_r );

    std::array< std::vector< tk::real >, 2 > coordgp
      { std::vector<tk::real>(ng), std::vector<tk::real>(ng) };
    std::vector< tk::real > wgp( ng );

    tk::GaussQuadratureTri( ng, coordgp, wgp );

    // Extract the element coordinates
    std::array< std::array< tk::real, 3>, 4 > coordel_l {{
      {{ cx[ inpoel[4*el  ] ], cy[ inpoel[4*el  ] ], cz[ inpoel[4*el  ] ] }},
      {{ cx[ inpoel[4*el+1] ], cy[ inpoel[4*el+1] ], cz[ inpoel[4*el+1] ] }},
      {{ cx[ inpoel[4*el+2] ], cy[ inpoel[4*el+2] ], cz[ inpoel[4*el+2] ] }},
      {{ cx[ inpoel[4*el+3] ], cy[ inpoel[4*el+3] ], cz[ inpoel[4*el+3] ] }} }};

    std::array< std::array< tk::real, 3>, 4 > coordel_r {{
      {{ cx[ inpoel[4*er  ] ], cy[ inpoel[4*er  ] ], cz[ inpoel[4*er  ] ] }},
      {{ cx[ inpoel[4*er+1] ], cy[ inpoel[4*er+1] ], cz[ inpoel[4*er+1] ] }},
      {{ cx[ inpoel[4*er+2] ], cy[ inpoel[4*er+2] ], cz[ inpoel[4*er+2] ] }},
      {{ cx[ inpoel[4*er+3] ], cy[ inpoel[4*er+3] ], cz[ inpoel[4*er+3] ] }} }};

    // Compute the determinant of Jacobian matrix
    auto detT_l =
      tk::Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );
    auto detT_r =
      tk::Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], coordel_r[3] );

    std::array< std::array< tk::real, 3>, 3 > coordfa {{
      {{ cx[ inpofa[3*f  ] ], cy[ inpofa[3*f  ] ], cz[ inpofa[3*f  ] ] }},
      {{ cx[ inpofa[3*f+1] ], cy[ inpofa[3*f+1] ], cz[ inpofa[3*f+1] ] }},
      {{ cx[ inpofa[3*f+2] ], cy[ inpofa[3*f+2] ], cz[ inpofa[3*f+2] ] }} }};

    std::array< tk::real, 3 >
      fn{{ geoFace(f,1), geoFace(f,2), geoFace(f,3) }};

    // Numerator and denominator of the shock indicator
    tk::real numer(0.0), denom(0.0);
    std::vector< tk::real > fl_jump, fl_avg;
    fl_jump.resize(3, 0.0);
    fl_avg.resize(3, 0.0);

    for (std::size_t igp=0; igp<ng; ++igp) {
      auto gp = tk::eval_gp( igp, coordfa, coordgp );
      std::size_t dof_el, dof_er;
      if (rdof > ndof)
      {
        dof_el = rdof;
        dof_er = rdof;
      }
      else
      {
        dof_el = ndofel[el];
        dof_er = ndofel[er];
      }
      std::array< tk::real, 3> ref_gp_l{
        tk::Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
        tk::Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
        tk::Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l };
      std::array< tk::real, 3> ref_gp_r{
        tk::Jacobian( coordel_r[0], gp, coordel_r[2], coordel_r[3] ) / detT_r,
        tk::Jacobian( coordel_r[0], coordel_r[1], gp, coordel_r[3] ) / detT_r,
        tk::Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], gp ) / detT_r };
      auto B_l = tk::eval_basis( dof_el, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2] );
      auto B_r = tk::eval_basis( dof_er, ref_gp_r[0], ref_gp_r[1], ref_gp_r[2] );

      std::array< std::vector< tk::real >, 2 > state;

      // Evaluate the high order solution at the qudrature point
      state[0] = tk::evalPolynomialSol(system, mat_blk, 0, ncomp, nprim, rdof,
        nmat, el, dof_el, inpoel, coord, geoElem, ref_gp_l, B_l, U, P);
      state[1] = tk::evalPolynomialSol(system, mat_blk, 0, ncomp, nprim, rdof,
        nmat, er, dof_er, inpoel, coord, geoElem, ref_gp_r, B_r, U, P);

      Assert( state[0].size() == ncomp+nprim, "Incorrect size for "
              "appended boundary state vector" );
      Assert( state[1].size() == ncomp+nprim, "Incorrect size for "
              "appended boundary state vector" );

      // Evaluate the flux
      auto fl = flux( system, ncomp, mat_blk, state[0], {} );
      auto fr = flux( system, ncomp, mat_blk, state[1], {} );

      for(std::size_t icomp = VarRange[0]; icomp <= VarRange[1]; icomp++) {
        tk::real fn_l(0.0), fn_r(0.0);
        for(std::size_t idir = 0; idir < 3; idir++) {
          fn_l += fl[icomp][idir] * fn[idir];
          fn_r += fr[icomp][idir] * fn[idir];
        }
        auto mark = icomp - VarRange[0];
        fl_jump[mark] += wgp[igp] * (fn_l - fn_r) * (fn_l - fn_r);
        fl_avg[mark]  += wgp[igp] * (fn_l + fn_r) * (fn_l + fn_r) * 0.25;
      }
    }

    // Evaluate the numerator and denominator
    for(std::size_t idir = 0; idir < 3; idir++) {
      numer += std::sqrt(fl_jump[idir]);
      denom += std::sqrt(fl_avg[idir]);
    }

    tk::real Ind(0.0);
    if(denom > 1e-8)
      Ind = numer / denom;
    IC[el] = std::max(IC[el], Ind);
    IC[er] = std::max(IC[er], Ind);
  }

  tk::real power = 0.0;
  if(rdof == 10)  power = 1.5;
  else            power = 1.0;

  // Loop over element to mark shock cell
  for (std::size_t e=0; e<nelem; ++e) {
    // Evaluate the threshold
    auto thres = coeff * std::pow(geoElem(e, 4), power);
    if(IC[e] > thres)
      shockmarker[e] = 1;
    else
      shockmarker[e] = 0;
  }
}

void
correctLimConservMultiMat(
  std::size_t nelem,
  const std::vector< EoS_Base* >& mat_blk,
  std::size_t nmat,
  const tk::Fields& geoElem,
  const tk::Fields& prim,
  tk::Fields& unk )
// *****************************************************************************
//  Update the conservative quantities after limiting for multi-material systems
//! \param[in] nelem Number of internal elements
//! \param[in] mat_blk EOS material block
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] geoElem Element geometry array
//! \param[in] prim Array of primitive variables
//! \param[in,out] unk Array of conservative variables
//! \details This function computes the updated dofs for conservative
//!   quantities based on the limited primitive quantities
// *****************************************************************************
{
  const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
  std::size_t ncomp = unk.nprop()/rdof;
  std::size_t nprim = prim.nprop()/rdof;

  for (std::size_t e=0; e<nelem; ++e) {
    // Here we pre-compute the right-hand-side vector. The reason that the
    // lhs in DG.cpp is not used is that the size of this vector in this
    // projection procedure should be rdof instead of ndof.
    auto L = tk::massMatrixDubiner(rdof, geoElem(e,0));

    std::vector< tk::real > R((nmat+3)*rdof, 0.0);

    auto ng = tk::NGvol(rdof);

    // Arrays for quadrature points
    std::array< std::vector< tk::real >, 3 > coordgp;
    std::vector< tk::real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    coordgp[2].resize( ng );
    wgp.resize( ng );

    tk::GaussQuadratureTet( ng, coordgp, wgp );

    // Loop over quadrature points in element e
    for (std::size_t igp=0; igp<ng; ++igp) {
      // Compute the basis function
      auto B = tk::eval_basis( rdof, coordgp[0][igp], coordgp[1][igp],
                               coordgp[2][igp] );

      auto w = wgp[igp] * geoElem(e, 0);

      // Evaluate the solution at quadrature point
      auto U = tk::eval_state( ncomp, rdof, rdof, e, unk,  B,
                               {0, ncomp-1} );
      auto P = tk::eval_state( nprim, rdof, rdof, e, prim, B,
                               {0, nprim-1} );

      // Solution vector that stores the material energy and bulk momentum
      std::vector< tk::real > s(nmat+3, 0.0);

      // Bulk density at quadrature point
      tk::real rhob(0.0);
      for (std::size_t k=0; k<nmat; ++k)
        rhob += U[densityIdx(nmat, k)];

      // Velocity vector at quadrature point
      std::array< tk::real, 3 >
        vel{ P[velocityIdx(nmat, 0)],
             P[velocityIdx(nmat, 1)],
             P[velocityIdx(nmat, 2)] };

      // Compute and store the bulk momentum
      for(std::size_t idir = 0; idir < 3; idir++)
        s[nmat+idir] = rhob * vel[idir];

      // Compute and store material energy at quadrature point
      for(std::size_t imat = 0; imat < nmat; imat++) {
        auto alphamat = U[volfracIdx(nmat, imat)];
        auto rhomat = U[densityIdx(nmat, imat)]/alphamat;
        auto premat = P[pressureIdx(nmat, imat)]/alphamat;
        s[imat] = alphamat * mat_blk[imat]->eos_totalenergy( rhomat,
          vel[0], vel[1], vel[2], premat );
      }

      // Evaluate the righ-hand-side vector
      for(std::size_t k = 0; k < nmat+3; k++) {
        auto mark = k * rdof;
        for(std::size_t idof = 0; idof < rdof; idof++)
          R[mark+idof] += w * s[k] * B[idof];
      }
    }

    // Update the high order dofs of the material energy
    for(std::size_t imat = 0; imat < nmat; imat++) {
      auto mark = imat * rdof;
      for(std::size_t idof = 1; idof < rdof; idof++)
        unk(e, energyDofIdx(nmat, imat, rdof, idof)) =
          R[mark+idof] / L[idof];
    }

    // Update the high order dofs of the bulk momentum
    for(std::size_t idir = 0; idir < 3; idir++) {
      auto mark = (nmat + idir) * rdof;
      for(std::size_t idof = 1; idof < rdof; idof++)
        unk(e, momentumDofIdx(nmat, idir, rdof, idof)) =
          R[mark+idof] / L[idof];
    }
  }
}


//! Constrain material partial pressure (alpha_k * p_k)
//! \param[in] apr Material partial pressure (alpha_k * p_k)
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for the
//!   single-material system, this argument can be left unspecified by the
//!   calling code
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified by
//!   the calling code
//! \return Constrained material partial pressure (alpha_k * p_k)
tk::real constrain_pressure( const std::vector< EoS_Base* >& mat_blk,
  tk::real apr,
  tk::real alpha=1.0,
  std::size_t imat=0 )
{
  return std::max(apr, alpha*mat_blk[imat]->min_eff_pressure(1e-12));
}


} // inciter::
