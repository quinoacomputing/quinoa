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
#include "PrefIndicator.hpp"
#include "Reconstruction.hpp"
#include "Integrate/Mass.hpp"
#include "MultiMat/MiscMultiMatFns.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"
#include "MultiSpecies/Mixture/Mixture.hpp"

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
  const auto rdof = inciter::g_inputdeck.get< tag::rdof >();
  const auto cweight = inciter::g_inputdeck.get< tag::cweight >();
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
  const auto rdof = inciter::g_inputdeck.get< tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::ndof >();
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
  const tk::UnsMesh::Coords& coord,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat )
// *****************************************************************************
//  Superbee limiter for multi-material DGP1
//! \param[in] esuel Elements surrounding elements
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] coord Array of nodal coordinates
//! \param[in] solidx Solid material index indicator
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in,out] P High-order vector of primitives which gets limited
//! \param[in] nmat Number of materials in this PDE system
//! \details This Superbee function should be called for multimat
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::ndof >();
  const auto intsharp = inciter::g_inputdeck.get< tag::multimat,
    tag::intsharp >();
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
        if (!g_inputdeck.get< tag::accuracy_test >())
          consistentMultiMatLimiting_P1(nmat, rdof, e, solidx, U, P, phic,
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
  const tk::UnsMesh::Coords& coord,
  tk::Fields& U )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for transport DGP1
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] coord Array of nodal coordinates
//! \param[in,out] U High-order solution vector which gets limited
//! \details This vertex-based limiter function should be called for transport.
//!   For details see: Kuzmin, D. (2010). A vertex-based hierarchical slope
//!   limiter for p-adaptive discontinuous Galerkin methods. Journal of
//!   computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::ndof >();
  const auto intsharp = inciter::g_inputdeck.get< tag::transport,
    tag::intsharp >();
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
      std::vector< std::size_t > var;
      for (std::size_t c=0; c<ncomp; ++c) var.push_back(c);
      // limit conserved quantities
      VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
        ncomp, phi, var);

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
  const std::vector< inciter::EOS >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const tk::FluxFn& flux,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  std::vector< std::size_t >& shockmarker )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for single-material DGP1
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] mat_blk EOS material block
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
// //! \param[in] geoElem Element geometry array
//! \param[in] coord Array of nodal coordinates
//! \param[in] flux Riemann flux function to use
//! \param[in] solidx Solid material index indicator
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in,out] shockmarker Shock detection marker array
//! \details This vertex-based limiter function should be called for compflow.
//!   For details see: Kuzmin, D. (2010). A vertex-based hierarchical slope
//!   limiter for p-adaptive discontinuous Galerkin methods. Journal of
//!   computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::ndof >();
  std::size_t ncomp = U.nprop()/rdof;

  // Null field for MarkShockCells argument
  tk::Fields P;

  if (inciter::g_inputdeck.get< tag::shock_detector_coeff >() > 1e-6) {
    // Indicator based on momentum flux jump
    std::set< std::size_t > vars;
    for (std::size_t i=1; i<=3; ++i) vars.insert(i);
    MarkShockCells(false, nelem, 1, ndof, rdof, mat_blk, ndofel,
      inpoel, coord, fd, geoFace, geoElem, flux, solidx, U, P,
      vars, shockmarker);
  }

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
      std::vector< std::size_t > var;
      for (std::size_t c=0; c<ncomp; ++c) var.push_back(c);
      // limit conserved quantities
      VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
        ncomp, phi, var);

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
  const std::vector< inciter::EOS >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  [[maybe_unused]] const std::vector< std::size_t >& gid,
  [[maybe_unused]] const std::unordered_map< std::size_t, std::size_t >& bid,
  [[maybe_unused]] const std::vector< std::vector<tk::real> >& uNodalExtrm,
  [[maybe_unused]] const std::vector< std::vector<tk::real> >& mtInv,
  const tk::FluxFn& flux,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  std::vector< std::size_t >& shockmarker )
// *****************************************************************************
//  Kuzmin's vertex-based limiter on reference element for single-material DGP2
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
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
//! \param[in] solidx Solid material index indicator
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in,out] shockmarker Shock detection marker array
//! \details This vertex-based limiter function should be called for compflow.
//!   For details see: Kuzmin, D. (2010). A vertex-based hierarchical slope
//!   limiter for p-adaptive discontinuous Galerkin methods. Journal of
//!   computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::ndof >();
  std::size_t ncomp = U.nprop()/rdof;

  // Null field for MarkShockCells argument
  tk::Fields P;

  if (inciter::g_inputdeck.get< tag::shock_detector_coeff >() > 1e-6) {
    // Indicator based on momentum flux jump
    std::set< std::size_t > vars;
    for (std::size_t i=1; i<=3; ++i) vars.insert(i);
    MarkShockCells(false, nelem, 1, ndof, rdof, mat_blk, ndofel,
      inpoel, coord, fd, geoFace, geoElem, flux, solidx, U, P,
      vars, shockmarker);
  }

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
      std::vector< std::size_t > var;
      for (std::size_t c=0; c<ncomp; ++c) var.push_back(c);
      VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
        ncomp, phic_p1, var);

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
  const std::vector< inciter::EOS >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const tk::FluxFn& flux,
  const std::vector< std::size_t >& solidx,
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
//! \param[in] mat_blk EOS material block
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
// //! \param[in] geoElem Element geometry array
//! \param[in] coord Array of nodal coordinates
//! \param[in] flux Riemann flux function to use
//! \param[in] solidx Solid material index indicator
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
  const auto rdof = inciter::g_inputdeck.get< tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::ndof >();
  const auto intsharp = inciter::g_inputdeck.get< tag::multimat,
    tag::intsharp >();
  std::size_t ncomp = U.nprop()/rdof;
  std::size_t nprim = P.nprop()/rdof;

  // Evaluate the interface condition and mark the shock cells
  if (inciter::g_inputdeck.get< tag::shock_detector_coeff >()
    > 1e-6 && ndof > 1) {
    // Indicator based on momentum flux jump
    std::set< std::size_t > vars;
    for (std::size_t i=0; i<3; ++i) vars.insert(momentumIdx(nmat, i));
    MarkShockCells(false, nelem, nmat, ndof, rdof, mat_blk, ndofel,
      inpoel, coord, fd, geoFace, geoElem, flux, solidx, U, P,
      vars, shockmarker);
  }

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
        std::vector< std::size_t > varc;
        for (std::size_t c=0; c<ncomp; ++c) varc.push_back(c);
        VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
          ncomp, phic, varc);
        // limit primitive quantities
        std::vector< std::size_t > varp;
        for (std::size_t c=0; c<nprim; ++c) varp.push_back(c);
        VertexBasedLimiting(P, esup, inpoel, coord, e, rdof, dof_el,
          nprim, phip, varp);
      } else {
        // When shockmarker is 0, the volume fraction, density and energy
        // of minor material will still be limited to ensure a stable solution.
        std::vector< std::size_t > vars;
        for (std::size_t k=0; k<nmat; ++k) vars.push_back(volfracIdx(nmat,k));
        VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
          ncomp, phic, vars);

        for(std::size_t k=0; k<nmat; ++k) {
          if(U(e, volfracDofIdx(nmat,k,rdof,0)) < 1e-4) {
            // limit the density and energy of minor materials
            vars.clear();
            vars.push_back(densityIdx(nmat, k));
            vars.push_back(energyIdx(nmat, k));
            if (solidx[k] > 0) {
              for (std::size_t i=0; i<3; ++i)
                for (std::size_t j=0; j<3; ++j)
                  vars.push_back(deformIdx(nmat, solidx[k], i, j));
            }
            VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
              ncomp, phic, vars);

            // limit the pressure of minor materials
            VertexBasedLimiting(P, esup, inpoel, coord, e, rdof, dof_el,
              nprim, phip, std::vector< std::size_t >{pressureIdx(nmat, k)});
          }
        }
      }

      std::vector< tk::real > phic_p2, phip_p2;

      PositivityLimiting(nmat, 1, mat_blk, rdof, dof_el, ndofel, e,
        inpoel, coord, fd.Esuel(), U, P, phic, phic_p2, phip, phip_p2);

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

        if (!g_inputdeck.get< tag::accuracy_test >())
          consistentMultiMatLimiting_P1(nmat, rdof, e, solidx, U, P, phic,
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
  const bool pref,
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  const std::vector< inciter::EOS >& mat_blk,
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
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat,
  std::vector< std::size_t >& shockmarker )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for multi-material DGP2
//! \param[in] pref Indicator for p-adaptive algorithm
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
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
//! \param[in] solidx Solid material index indicator
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
  const auto rdof = inciter::g_inputdeck.get< tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::ndof >();
  const auto intsharp = inciter::g_inputdeck.get< tag::multimat,
    tag::intsharp >();
  std::size_t ncomp = U.nprop()/rdof;
  std::size_t nprim = P.nprop()/rdof;

  // Evaluate the interface condition and mark the shock cells
  if (inciter::g_inputdeck.get< tag::shock_detector_coeff >()
    > 1e-6) {
    // Indicator based on momentum flux jump
    std::set< std::size_t > vars;
    for (std::size_t i=0; i<3; ++i) vars.insert(momentumIdx(nmat, i));
    MarkShockCells(pref, nelem, nmat, ndof, rdof, mat_blk, ndofel,
      inpoel, coord, fd, geoFace, geoElem, flux, solidx, U, P,
      vars, shockmarker);
  }

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

    // For multi-material simulation, when dofel = 1, p0p1 is applied and ndof
    // for solution evaluation should be 4
    if(ncomp > 5 && dof_el == 1)
      dof_el = 4;

    if (dof_el > 1)
    {
      // The vector of limiting coefficients for P1
      std::vector< tk::real > phic_p1(ncomp, 1.0), phic_p2(ncomp, 1.0);
      std::vector< tk::real > phip_p1(nprim, 1.0), phip_p2(nprim, 1.0);

      // Only when the cell is marked with discontinuous solution or P0P1 scheme
      // is used, the vertex-based slope limiter will be applied.
      if(shockmarker[e] || dof_el == 4) {
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
        std::vector< std::size_t > varc;
        for (std::size_t c=0; c<ncomp; ++c) varc.push_back(c);
        VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
          ncomp, phic_p1, varc);
        // Obtain limiter coefficient for P1 primitive quantities
        std::vector< std::size_t > varp;
        for (std::size_t c=0; c<nprim; ++c) varp.push_back(c);
        VertexBasedLimiting(P, esup, inpoel, coord, e, rdof, dof_el,
          nprim, phip_p1, varp);
      } else {
        // When shockmarker is 0, the volume fraction will still be limited to
        // ensure a stable solution. Since the limiting strategy for third order
        // solution will downgrade the accuracy to second order, the density,
        // energy and pressure of minor material will not be limited.
        std::vector< std::size_t > vars;
        for (std::size_t k=0; k<nmat; ++k) vars.push_back(volfracIdx(nmat,k));
        VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
          ncomp, phic_p1, vars);

        //for(std::size_t k=0; k<nmat; ++k) {
        //  if(U(e, volfracDofIdx(nmat,k,rdof,0)) < 1e-4) {
        //    // limit the density of minor materials
        //    VertexBasedLimiting(unk, U, esup, inpoel, coord, e, rdof, dof_el,
        //      ncomp, phic_p1, std::vector< std::size_t >{densityIdx(nmat,k)});

        //    // limit the pressure of minor materials
        //    VertexBasedLimiting(prim, P, esup, inpoel, coord, e, rdof, dof_el,
        //      nprim, phip_p1, std::vector< std::size_t >{pressureIdx(nmat,k)});
        //  }
        //}
      }

      PositivityLimiting(nmat, 1, mat_blk, ndof, dof_el, ndofel, e,
        inpoel, coord, fd.Esuel(), U, P, phic_p1, phic_p2, phip_p1, phic_p2);

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

        if (!g_inputdeck.get< tag::accuracy_test >())
          consistentMultiMatLimiting_P1(nmat, rdof, e, solidx, U, P, phic_p1,
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
  const tk::UnsMesh::Coords& coord,
  const std::vector< int >& srcFlag,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for multi-material FV
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] nelem Number of elements
//! \param[in] coord Array of nodal coordinates
//! \param[in] srcFlag Whether the energy source was added
//! \param[in] solidx Solid material index indicator
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in,out] P High-order vector of primitives which gets limited
//! \param[in] nmat Number of materials in this PDE system
//! \details This vertex-based limiter function should be called for multimat.
//!   For details see: Kuzmin, D. (2010). A vertex-based hierarchical slope
//!   limiter for p-adaptive discontinuous Galerkin methods. Journal of
//!   computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::rdof >();
  const auto intsharp = inciter::g_inputdeck.get< tag::multimat,
    tag::intsharp >();
  std::size_t ncomp = U.nprop()/rdof;
  std::size_t nprim = P.nprop()/rdof;

  for (std::size_t e=0; e<nelem; ++e)
  {
    std::vector< tk::real > phic(ncomp, 1.0);
    std::vector< tk::real > phip(nprim, 1.0);
    // limit conserved quantities
    std::vector< std::size_t > var;
    for (std::size_t k=0; k<nmat; ++k) {
      var.push_back(volfracIdx(nmat,k));
      var.push_back(densityIdx(nmat,k));
    }
    VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, rdof, ncomp,
      phic, var);
    // limit primitive quantities
    var.clear();
    for (std::size_t c=0; c<nprim; ++c) var.push_back(c);
    VertexBasedLimiting(P, esup, inpoel, coord, e, rdof, rdof, nprim,
      phip, var);

    // limits under which compression is to be performed
    std::vector< std::size_t > matInt(nmat, 0);
    std::vector< tk::real > alAvg(nmat, 0.0);
    for (std::size_t k=0; k<nmat; ++k)
      alAvg[k] = U(e, volfracDofIdx(nmat,k,rdof,0));
    auto intInd = interfaceIndicator(nmat, alAvg, matInt);
    if ((intsharp > 0) && intInd && srcFlag[e] == 0)
    {
      for (std::size_t k=0; k<nmat; ++k)
      {
        if (matInt[k])
          phic[volfracIdx(nmat,k)] = 1.0;
      }
    }
    else
    {
      if (!g_inputdeck.get< tag::accuracy_test >()) {
        std::vector< tk::real > phic_p2(ncomp, 1.0);
        consistentMultiMatLimiting_P1(nmat, rdof, e, solidx, U, P, phic,
          phic_p2);
      }
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

void
VertexBasedMultiSpecies_P1(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  const std::vector< inciter::EOS >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const tk::FluxFn& flux,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nspec,
  std::vector< std::size_t >& shockmarker )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for multi-species DGP1
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] mat_blk EOS material block
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] coord Array of nodal coordinates
//! \param[in] flux Riemann flux function to use
//! \param[in] solidx Solid material index indicator
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in,out] P High-order primitive vector which gets limited
//! \param[in] nspec Number of species in this PDE system
//! \param[in,out] shockmarker Shock detection marker array
//! \details This vertex-based limiter function should be called for
//!   multispecies. For details see: Kuzmin, D. (2010). A vertex-based
//!   hierarchical slope limiter for p-adaptive discontinuous Galerkin methods.
//!   Journal of computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::ndof >();
  std::size_t ncomp = U.nprop()/rdof;
  std::size_t nprim = P.nprop()/rdof;

  // Evaluate the interface condition and mark the shock cells
  if (inciter::g_inputdeck.get< tag::shock_detector_coeff >()
    > 1e-6 && ndof > 1) {
    // Indicator based on momentum flux jump
    std::set< std::size_t > vars;
    for (std::size_t i=0; i<3; ++i)
      vars.insert(multispecies::momentumIdx(nspec, i));
    MarkShockCells(false, nelem, 1, ndof, rdof, mat_blk,
      ndofel, inpoel, coord, fd, geoFace, geoElem, flux, solidx, U, P,
      vars, shockmarker);
  }

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
      std::vector< std::size_t > varc, varp;
      if(shockmarker[e]) {
        // When shockmarker is 1, there is discontinuity within the element.
        // Hence, the vertex-based limiter will be applied.

        for (std::size_t c=0; c<ncomp; ++c) varc.push_back(c);
        for (std::size_t c=0; c<nprim; ++c) varp.push_back(c);
      } else {
        // When shockmarker is 0, the density of minor species will still be
        // limited to ensure a stable solution.

        tk::real rhob(0.0);
        for(std::size_t k=0; k<nspec; ++k)
          rhob += U(e, multispecies::densityDofIdx(nspec,k,rdof,0));
        for(std::size_t k=0; k<nspec; ++k) {
          if (U(e, multispecies::densityDofIdx(nspec,k,rdof,0))/rhob < 1e-4) {
            // limit the density of minor species
            varc.push_back(multispecies::densityIdx(nspec, k));
          }
        }
      }

      std::vector< tk::real > phic(ncomp, 1.0), phip(nprim, 1.0);
      VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
        ncomp, phic, varc);
      if (!varp.empty())
        VertexBasedLimiting(P, esup, inpoel, coord, e, rdof, dof_el,
          nprim, phip, varp);

      std::vector< tk::real > phic_p2, phip_p2;
      PositivityLimiting(1, nspec, mat_blk, rdof, dof_el, ndofel, e,
        inpoel, coord, fd.Esuel(), U, P, phic, phic_p2, phip, phip_p2);

      // TODO: Unit sum of mass fractions is maintained by using common limiter
      // for all species densities. Investigate better approaches.
      if (!g_inputdeck.get< tag::accuracy_test >()) {
        tk::real phi_rhos_p1(1.0);
        for (std::size_t k=0; k<nspec; ++k)
          phi_rhos_p1 = std::min( phi_rhos_p1,
            phic[multispecies::densityIdx(nspec, k)] );
        // same limiter for all densities
        for (std::size_t k=0; k<nspec; ++k)
          phic[multispecies::densityIdx(nspec, k)] = phi_rhos_p1;
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
VertexBasedMultiSpecies_P2(
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t nelem,
  const std::vector< inciter::EOS >& mat_blk,
  const inciter::FaceData& fd,
  const tk::Fields& geoFace,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const tk::FluxFn& flux,
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nspec,
  std::vector< std::size_t >& shockmarker )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for multi-species DGP2
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] mat_blk EOS material block
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] coord Array of nodal coordinates
//! \param[in] flux Riemann flux function to use
//! \param[in] solidx Solid material index indicator
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in,out] P High-order primitive vector which gets limited
//! \param[in] nspec Number of species in this PDE system
//! \param[in,out] shockmarker Shock detection marker array
//! \details This vertex-based limiter function should be called for
//!   multispecies. For details see: Kuzmin, D. (2010). A vertex-based
//!   hierarchical slope limiter for p-adaptive discontinuous Galerkin methods.
//!   Journal of computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::ndof >();
  std::size_t ncomp = U.nprop()/rdof;
  std::size_t nprim = P.nprop()/rdof;

  // Evaluate the interface condition and mark the shock cells
  if (inciter::g_inputdeck.get< tag::shock_detector_coeff >() > 1e-6) {
    std::set< std::size_t > vars;
    for (std::size_t i=0; i<3; ++i)
      vars.insert(multispecies::momentumIdx(nspec, i));
    MarkShockCells(false, nelem, 1, ndof, rdof, mat_blk,
      ndofel, inpoel, coord, fd, geoFace, geoElem, flux, solidx, U, P,
      vars, shockmarker);
  }

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
      std::vector< std::size_t > varc, varp;
      if(shockmarker[e]) {
        // When shockmarker is 1, there is discontinuity within the element.
        // Hence, the vertex-based limiter will be applied.

        for (std::size_t c=0; c<ncomp; ++c) varc.push_back(c);
        for (std::size_t c=0; c<nprim; ++c) varp.push_back(c);
      }
        // When shockmarker is 0, the density of minor species is not limited
        // since it will degrade the accuracy to second order.

      // Removing 3rd order DOFs if discontinuity is detected, and applying
      // limiting to the 2nd order/P1 solution
      for (std::size_t c=0; c<varc.size(); c++) {
        for(std::size_t idof = 4; idof < rdof; idof++) {
          auto mark = varc[c] * rdof + idof;
          U(e, mark) = 0.0;
        }
      }
      for (std::size_t c=0; c<varp.size(); c++) {
        for(std::size_t idof = 4; idof < rdof; idof++) {
          auto mark = varp[c] * rdof + idof;
          P(e, mark) = 0.0;
        }
      }

      std::vector< tk::real > phic_p1(ncomp, 1.0), phip_p1(nprim, 1.0);
      if (!varc.empty())
        VertexBasedLimiting(U, esup, inpoel, coord, e, rdof, dof_el,
          ncomp, phic_p1, varc);
      if (!varp.empty())
        VertexBasedLimiting(P, esup, inpoel, coord, e, rdof, dof_el,
          nprim, phip_p1, varp);

      std::vector< tk::real > phic_p2(ncomp, 1.0), phip_p2(nprim, 1.0);
      PositivityLimiting(1, nspec, mat_blk, ndof, dof_el, ndofel, e,
        inpoel, coord, fd.Esuel(), U, P, phic_p1, phic_p2, phip_p1, phip_p2);

      // TODO: Unit sum of mass fractions is maintained by using common limiter
      // for all species densities. Investigate better approaches.
      if (!g_inputdeck.get< tag::accuracy_test >()) {
        tk::real phi_rhos_p1(1.0);
        for (std::size_t k=0; k<nspec; ++k)
          phi_rhos_p1 = std::min( phi_rhos_p1,
            phic_p1[multispecies::densityIdx(nspec, k)] );
        // same limiter for all densities
        for (std::size_t k=0; k<nspec; ++k)
          phic_p1[multispecies::densityIdx(nspec, k)] = phi_rhos_p1;
      }

      // apply limiter function to the solution
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
        tk::eval_state(ncomp, rdof, dof_el, e, U, B_l);

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
  const std::vector< std::size_t >& VarList )
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
//! \param[in,out] phi Limiter function for solution in element e
//! \param[in] VarList List of variable indices to be limited
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

  std::vector< tk::real > uMin(VarList.size(), 0.0),
                          uMax(VarList.size(), 0.0);

  // loop over all nodes of the element e
  for (std::size_t lp=0; lp<4; ++lp)
  {
    // reset min/max
    for (std::size_t i=0; i<VarList.size(); ++i)
    {
      auto mark = VarList[i]*rdof;
      uMin[i] = U(e, mark);
      uMax[i] = U(e, mark);
    }
    auto p = inpoel[4*e+lp];
    const auto& pesup = tk::cref_find(esup, p);

    // ----- Step-1: find min/max in the neighborhood of node p
    // loop over all the internal elements surrounding this node p
    for (auto er : pesup)
    {
      for (std::size_t i=0; i<VarList.size(); ++i)
      {
        auto mark = VarList[i]*rdof;
        uMin[i] = std::min(uMin[i], U(er, mark));
        uMax[i] = std::max(uMax[i], U(er, mark));
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
    state = tk::eval_state(ncomp, rdof, dof_el, e, U, B_p);

    Assert( state.size() == ncomp, "Size mismatch" );

    // compute the limiter function
    for (std::size_t i=0; i<VarList.size(); ++i)
    {
      auto c = VarList[i];
      auto phi_gp = 1.0;
      auto mark = c*rdof;
      auto uNeg = state[c] - U(e, mark);
      auto uref = std::max(std::fabs(U(e,mark)), 1e-14);
      if (uNeg > 1.0e-06*uref)
      {
        phi_gp = std::min( 1.0, (uMax[i]-U(e, mark))/uNeg );
      }
      else if (uNeg < -1.0e-06*uref)
      {
        phi_gp = std::min( 1.0, (uMin[i]-U(e, mark))/uNeg );
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
  const std::vector< std::size_t >& VarList,
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
//! \param[in] VarList List of variable indices that need to be limited
//! \param[out] phi Limiter function for solution in element e
//! \details This function limits the P2 dofs of P2 solution in a hierachical
//!   way to P1 dof limiting. Here we treat the first order derivatives the same
//!   way as cell average while second order derivatives represent the gradients
//!   to be limited in the P1 limiting procedure.
// *****************************************************************************
{
  const auto nelem = inpoel.size() / 4;

  std::vector< std::vector< tk::real > > uMin, uMax;
  uMin.resize( VarList.size(), std::vector<tk::real>(3, 0.0) );
  uMax.resize( VarList.size(), std::vector<tk::real>(3, 0.0) );

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
    for (std::size_t i=0; i<VarList.size(); ++i)
    {
      for (std::size_t idir=1; idir < 4; ++idir)
      {
        uMin[i][idir-1] = unk[VarList[i]][idir];
        uMax[i][idir-1] = unk[VarList[i]][idir];
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

        for (std::size_t i=0; i<VarList.size(); ++i)
        {
          auto mark = VarList[i]*rdof;
          for (std::size_t idir = 0; idir < 3; ++idir)
          {
            // The first order derivative at the centroid of element er
            tk::real slope_er(0.0);
            for(std::size_t idof = 1; idof < rdof; idof++)
              slope_er += U(er, mark+idof) * dBdxi_er[idir][idof];

            uMin[i][idir] = std::min(uMin[i][idir], slope_er);
            uMax[i][idir] = std::max(uMax[i][idir], slope_er);

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
      for (std::size_t i=0; i<VarList.size(); ++i)
      {
        for (std::size_t idir = 0; idir < 3; idir++)
        {
          auto max_mark = 2*VarList[i]*ndof_NodalExtrm + 2*idir;
          auto min_mark = max_mark + 1;
          const auto& ex = NodalExtrm[gip->second];
          uMax[i][idir] = std::max(ex[max_mark], uMax[i][idir]);
          uMin[i][idir] = std::min(ex[min_mark], uMin[i][idir]);
        }
      }
    }

    //Step-2: compute the limiter function at this node
    std::array< tk::real, 3 > node{cnodes[0][lp], cnodes[1][lp], cnodes[2][lp]};

    // find high-order solution
    std::vector< std::array< tk::real, 3 > > state;
    state.resize(VarList.size());

    for (std::size_t i=0; i<VarList.size(); ++i)
    {
      auto dx = node[0] - center[0][0];
      auto dy = node[1] - center[1][0];
      auto dz = node[2] - center[2][0];

      auto c = VarList[i];

      state[i][0] = unk[c][1] + unk[c][4]*dx + unk[c][7]*dy + unk[c][8]*dz;
      state[i][1] = unk[c][2] + unk[c][5]*dy + unk[c][7]*dx + unk[c][9]*dz;
      state[i][2] = unk[c][3] + unk[c][6]*dz + unk[c][8]*dx + unk[c][9]*dy;
    }

    // compute the limiter function
    for (std::size_t i=0; i<VarList.size(); ++i)
    {
      auto c = VarList[i];
      tk::real phi_dir(1.0);
      for (std::size_t idir = 1; idir <= 3; ++idir)
      {
        phi_dir = 1.0;
        auto uNeg = state[i][idir-1] - unk[c][idir];
        auto uref = std::max(std::fabs(unk[c][idir]), 1e-14);
        if (uNeg > 1.0e-6*uref)
        {
          phi_dir =
            std::min( 1.0, ( uMax[i][idir-1] - unk[c][idir])/uNeg );
        }
        else if (uNeg < -1.0e-6*uref)
        {
          phi_dir =
            std::min( 1.0, ( uMin[i][idir-1] - unk[c][idir])/uNeg );
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
  const std::vector< std::size_t >& solidx,
  tk::Fields& U,
  [[maybe_unused]] tk::Fields& P,
  std::vector< tk::real >& phic_p1,
  std::vector< tk::real >& phic_p2 )
// *****************************************************************************
//  Consistent limiter modifications for conservative variables
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] e Element being checked for consistency
//! \param[in] solidx Solid material index indicator
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
      if (solidx[k] > 0)
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
          {
            for (std::size_t idof=1; idof<rdof; ++idof)
              U(e,deformDofIdx(nmat,solidx[k],i,j,rdof,idof)) = 0.0;
          }
    }

    // 2. same limiter for all volume-fractions and densities
    for (std::size_t k=0; k<nmat; ++k)
    {
      phic_p1[volfracIdx(nmat, k)] = phi_al_p1;
      phic_p1[densityIdx(nmat, k)] = phi_al_p1;
      phic_p1[energyIdx(nmat, k)] = phi_al_p1;
      if (solidx[k] > 0)
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            phic_p1[deformIdx(nmat,solidx[k],i,j)] = phi_al_p1;
    }
    if(rdof > 4)
    {
      for (std::size_t k=0; k<nmat; ++k)
      {
        phic_p2[volfracIdx(nmat, k)] = phi_al_p2;
        phic_p2[densityIdx(nmat, k)] = phi_al_p2;
        phic_p2[energyIdx(nmat, k)] = phi_al_p2;
        if (solidx[k] > 0)
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              phic_p2[deformIdx(nmat,solidx[k],i,j)] = phi_al_p2;
      }
    }
  }
  else
  {
    // same limiter for all volume-fractions
    for (std::size_t k=0; k<nmat; ++k)
      phic_p1[volfracIdx(nmat, k)] = phi_al_p1;
    if(rdof > 4)
      for (std::size_t k=0; k<nmat; ++k)
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
  const tk::real max = 1.0 - min * static_cast<tk::real>(nmat - 1);

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

      auto state = eval_state( U.nprop()/ndof, ndof, ndof, e, U, B );

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

      auto state = tk::eval_state(U.nprop()/ndof, ndof, ndof, e, U, B);

      for(std::size_t imat = 0; imat < nmat; imat++)
      {
        auto phi = BoundPreservingLimitingFunction(min, max,
          state[volfracIdx(nmat, imat)],
          U(e,volfracDofIdx(nmat, imat, ndof, 0)) );
        phi_bound[imat] = std::min( phi_bound[imat], phi );
      }
    }
  }

  for(std::size_t k=0; k<nmat; k++)
    phic_p1[volfracIdx(nmat, k)] = std::min(phi_bound[k],
      phic_p1[volfracIdx(nmat, k)]);

  if(ndof > 4)
    for(std::size_t k=0; k<nmat; k++)
      phic_p2[volfracIdx(nmat, k)] = std::min(phi_bound[k],
        phic_p2[volfracIdx(nmat, k)]);
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

void PositivityLimiting( std::size_t nmat,
                         std::size_t nspec,
                         const std::vector< inciter::EOS >& mat_blk,
                         std::size_t rdof,
                         std::size_t ndof_el,
                         const std::vector< std::size_t >& ndofel,
                         std::size_t e,
                         const std::vector< std::size_t >& inpoel,
                         const tk::UnsMesh::Coords& coord,
                         const std::vector< int >& esuel,
                         const tk::Fields& U,
                         const tk::Fields& P,
                         std::vector< tk::real >& phic_p1,
                         std::vector< tk::real >& phic_p2,
                         std::vector< tk::real >& phip_p1,
                         std::vector< tk::real >& phip_p2 )
// *****************************************************************************
//  Positivity preserving limiter for multi-material and multspecies solver
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] nspec Number of species in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] ndof_el Number of dofs for element e
//! \param[in] ndofel Vector of local number of degrees of freedome
//! \param[in] e Element being checked for consistency
//! \param[in] inpoel Element connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] esuel Elements surrounding elements
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

  for (std::size_t lf=0; lf<4; ++lf)
  {
    std::array< std::size_t, 3 > inpofa_l {{ inpoel[4*e+tk::lpofa[lf][0]],
                                             inpoel[4*e+tk::lpofa[lf][1]],
                                             inpoel[4*e+tk::lpofa[lf][2]] }};

    std::array< std::array< tk::real, 3>, 3 > coordfa {{
      {{ cx[ inpofa_l[0] ], cy[ inpofa_l[0] ], cz[ inpofa_l[0] ] }},
      {{ cx[ inpofa_l[1] ], cy[ inpofa_l[1] ], cz[ inpofa_l[1] ] }},
      {{ cx[ inpofa_l[2] ], cy[ inpofa_l[2] ], cz[ inpofa_l[2] ] }} }};

    std::size_t nel;
    if (esuel[4*e+lf] == -1) nel = e;
    else nel = static_cast< std::size_t >(esuel[4*e+lf]);

    auto ng = tk::NGfa(std::max(ndofel[e], ndofel[nel]));

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

      auto state = eval_state(ncomp, rdof, ndof_el, e, U, B);
      auto sprim = eval_state(nprim, rdof, ndof_el, e, P, B);

      if (nmat > 1) {
        // multi-material PDE bounds
        PositivityBoundsMultiMat(nmat, mat_blk, rdof, e, U, P, state, sprim,
          phic_bound, phip_bound);
      }
      else {
        // multispecies PDE bounds
        PositivityBoundsMultiSpecies(nspec, mat_blk, rdof, e, U, P, state,
          sprim, phic_bound, phip_bound);
      }
    }
  }

  if(ndofel[e] > 4)
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

      auto state = eval_state(ncomp, rdof, ndof_el, e, U, B);
      auto sprim = eval_state(nprim, rdof, ndof_el, e, P, B);

      if (nmat > 1) {
        // multi-material PDE bounds
        PositivityBoundsMultiMat(nmat, mat_blk, rdof, e, U, P, state, sprim,
          phic_bound, phip_bound);
      }
      else {
        // multispecies PDE bounds
        PositivityBoundsMultiSpecies(nspec, mat_blk, rdof, e, U, P, state,
          sprim, phic_bound, phip_bound);
      }
    }
  }

  // apply new bounds
  for (std::size_t c=0; c<ncomp; ++c) {
    phic_p1[c] = std::min( phic_bound[c], phic_p1[c] );
    if (ndof_el > 4) phic_p2[c] = std::min( phic_bound[c], phic_p2[c] );
  }
  for (std::size_t c=0; c<nprim; ++c) {
    phip_p1[c] = std::min( phip_bound[c], phip_p1[c] );
    if (ndof_el > 4) phip_p2[c] = std::min( phip_bound[c], phip_p2[c] );
  }
}

void PositivityBoundsMultiMat(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t rdof,
  std::size_t e,
  const tk::Fields& U,
  const tk::Fields& P,
  const std::vector< tk::real >& state,
  const std::vector< tk::real >& sprim,
  std::vector< tk::real >& phic_bound,
  std::vector< tk::real >& phip_bound )
// *****************************************************************************
//  Positivity bounds for multi-material PDE system
//! \param[in] nmat Number of materials in this system
//! \param[in] mat_blk EOS material block
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] e Element for which bounds are being calculated
//! \param[in] U Vector of conservative variables
//! \param[in] P Vector of primitive variables
//! \param[in] state Vector of state of conserved quantities at quadrature point
//! \param[in] sprim Vector of state of primitive quantities at quadrature point
//! \param[in,out] phic_bound Vector of bounding-limiter functions for dofs of
//!   the conserved quantities
//! \param[in,out] phip_bound Vector of bounding-limiter functions for dofs of
//!   the primitive quantities
// *****************************************************************************
{
  const tk::real min = 1e-15;

  for(std::size_t k = 0; k < nmat; k++)
  {
    tk::real phi_rho(1.0), phi_rhoe(1.0), phi_pre(1.0);
    // Evaluate the limiting coefficient for material density
    auto rho = state[densityIdx(nmat, k)];
    auto rho_avg = U(e, densityDofIdx(nmat, k, rdof, 0));
    phi_rho = PositivityFunction(min, rho, rho_avg);
    phic_bound[densityIdx(nmat, k)] =
      std::min(phic_bound[densityIdx(nmat, k)], phi_rho);
    // Evaluate the limiting coefficient for material energy
    auto rhoe = state[energyIdx(nmat, k)];
    auto rhoe_avg = U(e, energyDofIdx(nmat, k, rdof, 0));
    phi_rhoe = PositivityFunction(min, rhoe, rhoe_avg);
    phic_bound[energyIdx(nmat, k)] =
      std::min(phic_bound[energyIdx(nmat, k)], phi_rhoe);
    // Evaluate the limiting coefficient for material pressure
    auto min_pre = std::max(min, state[volfracIdx(nmat, k)] *
      mat_blk[k].compute< EOS::min_eff_pressure >(min, rho,
      state[volfracIdx(nmat, k)]));
    auto pre = sprim[pressureIdx(nmat, k)];
    auto pre_avg = P(e, pressureDofIdx(nmat, k, rdof, 0));
    phi_pre = PositivityFunction(min_pre, pre, pre_avg);
    phip_bound[pressureIdx(nmat, k)] =
      std::min(phip_bound[pressureIdx(nmat, k)], phi_pre);
  }
}

void PositivityBoundsMultiSpecies(
  std::size_t nspec,
  const std::vector< inciter::EOS >& /*mat_blk*/,
  std::size_t rdof,
  std::size_t e,
  const tk::Fields& /*U*/,
  const tk::Fields& P,
  const std::vector< tk::real >& /*state*/,
  const std::vector< tk::real >& sprim,
  std::vector< tk::real >& /*phic_bound*/,
  std::vector< tk::real >& phip_bound )
// *****************************************************************************
//  Positivity bounds for multispecies PDE system
//! \param[in] nmat Number of materials in this system
// //! \param[in] mat_blk EOS material block
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] e Element for which bounds are being calculated
// //! \param[in] U Vector of conservative variables
//! \param[in] P Vector of primitive variables
// //! \param[in] state Vector of state of conserved quantities at quadrature point
//! \param[in] sprim Vector of state of primitive quantities at quadrature point
// //! \param[in,out] phic_bound Vector of bounding-limiter functions for dofs of
// //!   the conserved quantities
//! \param[in,out] phip_bound Vector of bounding-limiter functions for dofs of
//!   the primitive quantities
// *****************************************************************************
{
  const tk::real min = 1e-8;

  tk::real phi_T(1.0);
  // Evaluate the limiting coefficient for mixture temperature
  auto T = sprim[multispecies::temperatureIdx(nspec, 0)];
  auto T_avg = P(e, multispecies::temperatureDofIdx(nspec, 0, rdof, 0));
  phi_T = PositivityFunction(min, T, T_avg);
  phip_bound[multispecies::temperatureIdx(nspec, 0)] =
    std::min(phip_bound[multispecies::temperatureIdx(nspec, 0)], phi_T);
}

void PositivityPreservingMultiMat_FV(
  const std::vector< std::size_t >& inpoel,
  std::size_t nelem,
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const tk::UnsMesh::Coords& coord,
  const tk::Fields& /*geoFace*/,
  tk::Fields& U,
  tk::Fields& P )
// *****************************************************************************
//  Positivity preserving limiter for the FV multi-material solver
//! \param[in] inpoel Element connectivity
//! \param[in] nelem Number of elements
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk Material EOS block
//! \param[in] coord Array of nodal coordinates
////! \param[in] geoFace Face geometry array
//! \param[in,out] U High-order solution vector which gets limited
//! \param[in,out] P High-order vector of primitives which gets limited
//! \details This positivity preserving limiter function should be called for
//!   FV multimat.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::rdof >();
  const auto ncomp = U.nprop() / rdof;
  const auto nprim = P.nprop() / rdof;

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  for (std::size_t e=0; e<nelem; ++e)
  {
    // Extract the element coordinates
    std::array< std::array< tk::real, 3>, 4 > coordel {{
      {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
      {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
      {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
      {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }} }};

    // Compute the determinant of Jacobian matrix
    auto detT =
      tk::Jacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

    std::vector< tk::real > phic(ncomp, 1.0);
    std::vector< tk::real > phip(nprim, 1.0);

    const tk::real min = 1e-15;

    // 1. Enforce positive density (total energy will be positive if pressure
    //    and density are positive)
    for (std::size_t lf=0; lf<4; ++lf)
    {
      std::array< std::size_t, 3 > inpofa_l {{ inpoel[4*e+tk::lpofa[lf][0]],
                                               inpoel[4*e+tk::lpofa[lf][1]],
                                               inpoel[4*e+tk::lpofa[lf][2]] }};

      // face coordinates
      std::array< std::array< tk::real, 3>, 3 > coordfa {{
        {{ cx[ inpofa_l[0] ], cy[ inpofa_l[0] ], cz[ inpofa_l[0] ] }},
        {{ cx[ inpofa_l[1] ], cy[ inpofa_l[1] ], cz[ inpofa_l[1] ] }},
        {{ cx[ inpofa_l[2] ], cy[ inpofa_l[2] ], cz[ inpofa_l[2] ] }} }};

      // face centroid
      std::array< tk::real, 3 > fc{{
        (coordfa[0][0]+coordfa[1][0]+coordfa[2][0])/3.0 ,
        (coordfa[0][1]+coordfa[1][1]+coordfa[2][1])/3.0 ,
        (coordfa[0][2]+coordfa[1][2]+coordfa[2][2])/3.0 }};

      auto B = tk::eval_basis( rdof,
            tk::Jacobian( coordel[0], fc, coordel[2], coordel[3] ) / detT,
            tk::Jacobian( coordel[0], coordel[1], fc, coordel[3] ) / detT,
            tk::Jacobian( coordel[0], coordel[1], coordel[2], fc ) / detT );
      auto state = eval_state(ncomp, rdof, rdof, e, U, B);

      for(std::size_t i=0; i<nmat; i++)
      {
        // Evaluate the limiting coefficient for material density
        auto rho = state[densityIdx(nmat, i)];
        auto rho_avg = U(e, densityDofIdx(nmat, i, rdof, 0));
        auto phi_rho = PositivityFunction(min, rho, rho_avg);
        phic[densityIdx(nmat, i)] =
          std::min(phic[densityIdx(nmat, i)], phi_rho);
      }
    }
    // apply limiter coefficient
    for(std::size_t i=0; i<nmat; i++)
    {
      U(e, densityDofIdx(nmat,i,rdof,1)) *= phic[densityIdx(nmat,i)];
      U(e, densityDofIdx(nmat,i,rdof,2)) *= phic[densityIdx(nmat,i)];
      U(e, densityDofIdx(nmat,i,rdof,3)) *= phic[densityIdx(nmat,i)];
    }

    // 2. Enforce positive pressure (assuming density is positive)
    for (std::size_t lf=0; lf<4; ++lf)
    {
      std::array< std::size_t, 3 > inpofa_l {{ inpoel[4*e+tk::lpofa[lf][0]],
                                               inpoel[4*e+tk::lpofa[lf][1]],
                                               inpoel[4*e+tk::lpofa[lf][2]] }};

      // face coordinates
      std::array< std::array< tk::real, 3>, 3 > coordfa {{
        {{ cx[ inpofa_l[0] ], cy[ inpofa_l[0] ], cz[ inpofa_l[0] ] }},
        {{ cx[ inpofa_l[1] ], cy[ inpofa_l[1] ], cz[ inpofa_l[1] ] }},
        {{ cx[ inpofa_l[2] ], cy[ inpofa_l[2] ], cz[ inpofa_l[2] ] }} }};

      // face centroid
      std::array< tk::real, 3 > fc{{
        (coordfa[0][0]+coordfa[1][0]+coordfa[2][0])/3.0 ,
        (coordfa[0][1]+coordfa[1][1]+coordfa[2][1])/3.0 ,
        (coordfa[0][2]+coordfa[1][2]+coordfa[2][2])/3.0 }};

      auto B = tk::eval_basis( rdof,
            tk::Jacobian( coordel[0], fc, coordel[2], coordel[3] ) / detT,
            tk::Jacobian( coordel[0], coordel[1], fc, coordel[3] ) / detT,
            tk::Jacobian( coordel[0], coordel[1], coordel[2], fc ) / detT );
      auto state = eval_state(ncomp, rdof, rdof, e, U, B);
      auto sprim = eval_state(nprim, rdof, rdof, e, P, B);

      for(std::size_t i=0; i<nmat; i++)
      {
        tk::real phi_pre(1.0);
        // Evaluate the limiting coefficient for material pressure
        auto rho = state[densityIdx(nmat, i)];
        auto min_pre = std::max(min, U(e,volfracDofIdx(nmat,i,rdof,0)) *
          mat_blk[i].compute< EOS::min_eff_pressure >(min, rho,
          U(e,volfracDofIdx(nmat,i,rdof,0))));
        auto pre = sprim[pressureIdx(nmat, i)];
        auto pre_avg = P(e, pressureDofIdx(nmat, i, rdof, 0));
        phi_pre = PositivityFunction(min_pre, pre, pre_avg);
        phip[pressureIdx(nmat, i)] =
          std::min(phip[pressureIdx(nmat, i)], phi_pre);
      }
    }
    // apply limiter coefficient
    for(std::size_t i=0; i<nmat; i++)
    {
      P(e, pressureDofIdx(nmat,i,rdof,1)) *= phip[pressureIdx(nmat,i)];
      P(e, pressureDofIdx(nmat,i,rdof,2)) *= phip[pressureIdx(nmat,i)];
      P(e, pressureDofIdx(nmat,i,rdof,3)) *= phip[pressureIdx(nmat,i)];
    }
  }
}

tk::real
PositivityFunction( const tk::real min,
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

void MarkShockCells ( const bool pref,
                      const std::size_t nelem,
                      const std::size_t nmat,
                      const std::size_t ndof,
                      const std::size_t rdof,
                      const std::vector< inciter::EOS >& mat_blk,
                      const std::vector< std::size_t >& ndofel,
                      const std::vector< std::size_t >& inpoel,
                      const tk::UnsMesh::Coords& coord,
                      const inciter::FaceData& fd,
                      [[maybe_unused]] const tk::Fields& geoFace,
                      const tk::Fields& geoElem,
                      const tk::FluxFn& flux,
                      const std::vector< std::size_t >& solidx,
                      const tk::Fields& U,
                      const tk::Fields& P,
                      const std::set< std::size_t >& vars,
                      std::vector< std::size_t >& shockmarker )
// *****************************************************************************
//  Mark the cells that contain discontinuity according to the interface
//    condition
//! \param[in] pref Indicator for p-adaptive algorithm
//! \param[in] nelem Number of elements
//! \param[in] nmat Number of materials in this PDE system
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
//! \param[in] solidx Solid material index indicator
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] vars Vector of variable indices to evaluate flux jump
//! \param[in, out] shockmarker Vector of the shock indicator
//! \details This function computes the discontinuity indicator based on
//!   interface conditon. It is based on the following paper:
//!   Hong L., Gianni A., Robert N. (2021) A moving discontinuous Galerkin
//!   finite element method with interface condition enforcement for
//!   compressible flows. Journal of Computational Physics,
//!   doi: https://doi.org/10.1016/j.jcp.2021.110618
// *****************************************************************************
{
  const auto coeff = g_inputdeck.get< tag::shock_detector_coeff >();

  std::vector< tk::real > IC(U.nunk(), 0.0);
  const auto& esuf = fd.Esuf();
  const auto& inpofa = fd.Inpofa();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

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

      // Evaluate the high order solution at the qudrature point
      std::array< std::vector< tk::real >, 2 > state;
      state[0] = tk::evalPolynomialSol(mat_blk, 0, ncomp, nprim, rdof,
        nmat, el, dof_el, inpoel, coord, geoElem, ref_gp_l, B_l, U, P);
      state[1] = tk::evalPolynomialSol(mat_blk, 0, ncomp, nprim, rdof,
        nmat, er, dof_er, inpoel, coord, geoElem, ref_gp_r, B_r, U, P);

      Assert( state[0].size() == ncomp+nprim, "Incorrect size for "
              "appended boundary state vector" );
      Assert( state[1].size() == ncomp+nprim, "Incorrect size for "
              "appended boundary state vector" );

      // Force deformation unknown to first order
      for (std::size_t k=0; k<nmat; ++k)
        if (solidx[k] > 0)
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
            {
              state[0][deformIdx(nmat, solidx[k], i, j)] = U(el,deformDofIdx(
                nmat, solidx[k], i, j, rdof, 0));
              state[1][deformIdx(nmat, solidx[k], i, j)] = U(er,deformDofIdx(
                nmat, solidx[k], i, j, rdof, 0));
            }

      // Evaluate the flux
      auto fl = flux( ncomp, mat_blk, state[0], {} );
      auto fr = flux( ncomp, mat_blk, state[1], {} );

      std::size_t i(0);
      for (const auto& c : vars) {
        tk::real fn_l(0.0), fn_r(0.0);
        for(std::size_t idir = 0; idir < 3; idir++) {
          fn_l += fl[c][idir] * fn[idir];
          fn_r += fr[c][idir] * fn[idir];
        }
        fl_jump[i] += wgp[igp] * (fn_l - fn_r) * (fn_l - fn_r);
        fl_avg[i]  += wgp[igp] * (fn_l + fn_r) * (fn_l + fn_r) * 0.25;
        ++i;
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

  // Loop over element to mark shock cell
  for (std::size_t e=0; e<nelem; ++e) {
    std::size_t dof_el = pref ? ndofel[e] : rdof;

    tk::real power = 0.0;
    if(dof_el == 10)  power = 1.5;
    else              power = 1.0;

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
  const std::vector< EOS >& mat_blk,
  std::size_t nmat,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const tk::Fields& geoElem,
  const tk::Fields& prim,
  tk::Fields& unk )
// *****************************************************************************
//  Update the conservative quantities after limiting for multi-material systems
//! \param[in] nelem Number of internal elements
//! \param[in] mat_blk EOS material block
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] prim Array of primitive variables
//! \param[in,out] unk Array of conservative variables
//! \details This function computes the updated dofs for conservative
//!   quantities based on the limited primitive quantities, to re-instate
//!   consistency between the limited primitive and evolved quantities. For
//!   further details, see Pandare et al. (2023). On the Design of Stable,
//!   Consistent, and Conservative High-Order Methods for Multi-Material
//!   Hydrodynamics. J Comp Phys, 112313.
// *****************************************************************************
{
  const auto rdof = g_inputdeck.get< tag::rdof >();
  std::size_t ncomp = unk.nprop()/rdof;
  std::size_t nprim = prim.nprop()/rdof;
  const auto intsharp = inciter::g_inputdeck.get< tag::multimat,
    tag::intsharp >();

  for (std::size_t e=0; e<nelem; ++e) {
    // Here we pre-compute the right-hand-side vector. The reason that the
    // lhs in DG.cpp is not used is that the size of this vector in this
    // projection procedure should be rdof instead of ndof.
    auto L = tk::massMatrixDubiner(rdof, geoElem(e,0));

    // The right-hand side vector is sized as nprim, i.e. the primitive quantity
    // vector. However, it stores the consistently obtained values of evolved
    // quantities, since nprim is the number of evolved quantities that need to
    // be evaluated consistently. For this reason, accessing R will require
    // the primitive quantity accessors. But this access is intended to give
    // the corresponding evolved quantites, as follows:
    // pressureIdx() - mat. total energy
    // velocityIdx() - bulk momentum components
    // stressIdx() - mat. inverse deformation gradient tensor components
    std::vector< tk::real > R(nprim*rdof, 0.0);

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
      auto state = evalPolynomialSol(mat_blk, intsharp, ncomp, nprim,
        rdof, nmat, e, rdof, inpoel, coord, geoElem,
        {{coordgp[0][igp], coordgp[1][igp], coordgp[2][igp]}}, B, unk, prim);

      // Solution vector that stores the material energy and bulk momentum
      std::vector< tk::real > s(nprim, 0.0);

      // Bulk density at quadrature point
      tk::real rhob(0.0);
      for (std::size_t k=0; k<nmat; ++k)
        rhob += state[densityIdx(nmat, k)];

      // Velocity vector at quadrature point
      std::array< tk::real, 3 >
        vel{ state[ncomp+velocityIdx(nmat, 0)],
             state[ncomp+velocityIdx(nmat, 1)],
             state[ncomp+velocityIdx(nmat, 2)] };

      // Compute and store the bulk momentum
      for(std::size_t idir = 0; idir < 3; idir++)
        s[velocityIdx(nmat, idir)] = rhob * vel[idir];

      // Compute and store material energy at quadrature point
      for(std::size_t imat = 0; imat < nmat; imat++) {
        auto alphamat = state[volfracIdx(nmat, imat)];
        auto rhomat = state[densityIdx(nmat, imat)]/alphamat;
        auto premat = state[ncomp+pressureIdx(nmat, imat)]/alphamat;
        auto gmat = getDeformGrad(nmat, imat, state);
        s[pressureIdx(nmat,imat)] = alphamat *
          mat_blk[imat].compute< EOS::totalenergy >( rhomat, vel[0], vel[1],
          vel[2], premat, gmat );
      }

      // Evaluate the righ-hand-side vector
      for(std::size_t k = 0; k < nprim; k++) {
        auto mark = k * rdof;
        for(std::size_t idof = 0; idof < rdof; idof++)
          R[mark+idof] += w * s[k] * B[idof];
      }
    }

    // Update the high order dofs of the material energy
    for(std::size_t imat = 0; imat < nmat; imat++) {
      for(std::size_t idof = 1; idof < rdof; idof++)
        unk(e, energyDofIdx(nmat, imat, rdof, idof)) =
          R[pressureDofIdx(nmat,imat,rdof,idof)] / L[idof];
    }

    // Update the high order dofs of the bulk momentum
    for(std::size_t idir = 0; idir < 3; idir++) {
      for(std::size_t idof = 1; idof < rdof; idof++)
        unk(e, momentumDofIdx(nmat, idir, rdof, idof)) =
          R[velocityDofIdx(nmat,idir,rdof,idof)] / L[idof];
    }
  }
}

void
correctLimConservMultiSpecies(
  std::size_t nelem,
  const std::vector< EOS >& mat_blk,
  std::size_t nspec,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const tk::Fields& geoElem,
  const tk::Fields& prim,
  tk::Fields& unk )
// *****************************************************************************
//  Update the conservative quantities after limiting for multispecies systems
//! \param[in] nelem Number of internal elements
//! \param[in] mat_blk EOS material block
//! \param[in] nspec Number of species in this PDE system
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] prim Array of primitive variables
//! \param[in,out] unk Array of conservative variables
//! \details This function computes the updated dofs for conservative
//!   quantities based on the limited primitive quantities, to re-instate
//!   consistency between the limited primitive and evolved quantities. For
//!   further details, see Pandare et al. (2023). On the Design of Stable,
//!   Consistent, and Conservative High-Order Methods for Multi-Material
//!   Hydrodynamics. J Comp Phys, 112313.
// *****************************************************************************
{
  const auto rdof = g_inputdeck.get< tag::rdof >();
  std::size_t ncomp = unk.nprop()/rdof;
  std::size_t nprim = prim.nprop()/rdof;

  for (std::size_t e=0; e<nelem; ++e) {
    // Here we pre-compute the right-hand-side vector. The reason that the
    // lhs in DG.cpp is not used is that the size of this vector in this
    // projection procedure should be rdof instead of ndof.
    auto L = tk::massMatrixDubiner(rdof, geoElem(e,0));

    // The right-hand side vector is sized as nprim, i.e. the primitive quantity
    // vector. However, it stores the consistently obtained values of evolved
    // quantities, since nprim is the number of evolved quantities that need to
    // be evaluated consistently. For this reason, accessing R will require
    // the primitive quantity accessors. But this access is intended to give
    // the corresponding evolved quantites, as follows:
    // multispecies::tempratureIdx() - total energy
    std::vector< tk::real > R(nprim*rdof, 0.0);

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
      auto state = evalPolynomialSol(mat_blk, 0, ncomp, nprim,
        rdof, 1, e, rdof, inpoel, coord, geoElem,
        {{coordgp[0][igp], coordgp[1][igp], coordgp[2][igp]}}, B, unk, prim);

      // Solution vector that stores the material energy and bulk momentum
      std::vector< tk::real > s(nprim, 0.0);

      // Mixture state at quadrature point
      Mixture mixgp(nspec, state, mat_blk);

      // Mixture density at quadrature point
      tk::real rhob = mixgp.get_mix_density();

      // velocity vector at quadrature point
      std::array< tk::real, 3 >
        vel{ state[multispecies::momentumIdx(nspec,0)]/rhob,
             state[multispecies::momentumIdx(nspec,1)]/rhob,
             state[multispecies::momentumIdx(nspec,2)]/rhob };

      // Compute and store total energy at quadrature point
      s[multispecies::temperatureIdx(nspec,0)] = mixgp.totalenergy(rhob,
        vel[0], vel[1], vel[2],
        state[ncomp+multispecies::temperatureIdx(nspec,0)], mat_blk);

      // Evaluate the right-hand-side vector
      for(std::size_t k = 0; k < nprim; k++) {
        auto mark = k * rdof;
        for(std::size_t idof = 0; idof < rdof; idof++)
          R[mark+idof] += w * s[k] * B[idof];
      }
    }

    // Update the high order dofs of the total energy
    for(std::size_t idof = 1; idof < rdof; idof++)
      unk(e, multispecies::energyDofIdx(nspec,0,rdof,idof)) =
        R[multispecies::temperatureDofIdx(nspec,0,rdof,idof)] / L[idof];
  }
}

tk::real
constrain_pressure( const std::vector< EOS >& mat_blk,
  tk::real apr,
  tk::real arho,
  tk::real alpha=1.0,
  std::size_t imat=0 )
// *****************************************************************************
//  Constrain material partial pressure (alpha_k * p_k)
//! \param[in] apr Material partial pressure (alpha_k * p_k)
//! \param[in] arho Material partial density (alpha_k * rho_k)
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for the
//!   single-material system, this argument can be left unspecified by the
//!   calling code
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified by
//!   the calling code
//! \return Constrained material partial pressure (alpha_k * p_k)
// *****************************************************************************
{
  return std::max(apr, alpha*mat_blk[imat].compute<
    EOS::min_eff_pressure >(1e-12, arho, alpha));
}


} // inciter::
