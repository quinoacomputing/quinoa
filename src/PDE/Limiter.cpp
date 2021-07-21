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

#include "Vector.hpp"
#include "Limiter.hpp"
#include "DerivedData.hpp"
#include "Integrate/Quadrature.hpp"
#include "Integrate/Basis.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Inciter/PrefIndicator.hpp"

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
      WENOLimiting(U, esuel, e, c, rdof, offset, cweight, limU);
    }

    auto mark = c*rdof;

    for (std::size_t e=0; e<nelem; ++e)
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
             tk::Fields& U )
// *****************************************************************************
//  Superbee limiter for DGP1
//! \param[in] esuel Elements surrounding elements
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] offset Index for equation systems
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
                   dof_el, offset, ncomp, beta_lim);

      // apply limiter function
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

void
SuperbeeMultiMat_P1(
  const std::vector< int >& esuel,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& ndofel,
  std::size_t system,
  inciter::ncomp_t offset,
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
//! \param[in] offset Offset this PDE system operates from
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
                    dof_el, offset, ncomp, beta_lim);
      // limit primitive quantities
      auto phip = SuperbeeLimiting(P, esuel, inpoel, coord, e, ndof, rdof,
                    dof_el, offset, nprim, beta_lim);

      std::vector< tk::real > phic_p2;
      std::vector< std::vector< tk::real > > unk, prim;
      if(ndof > 1)
        BoundPreservingLimiting(nmat, offset, ndof, e, inpoel, coord, U, phic,
          phic_p2);

      // limits under which compression is to be performed
      std::vector< std::size_t > matInt(nmat, 0);
      std::vector< tk::real > alAvg(nmat, 0.0);
      for (std::size_t k=0; k<nmat; ++k)
        alAvg[k] = U(e, volfracDofIdx(nmat,k,rdof,0), offset);
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
        consistentMultiMatLimiting_P1(nmat, offset, rdof, e, U, P, unk, prim,
          phic, phic_p2);
      }

      // apply limiter function
      for (inciter::ncomp_t c=0; c<ncomp; ++c)
      {
        auto mark = c*rdof;
        U(e, mark+1, offset) = phic[c] * U(e, mark+1, offset);
        U(e, mark+2, offset) = phic[c] * U(e, mark+2, offset);
        U(e, mark+3, offset) = phic[c] * U(e, mark+3, offset);
      }
      for (inciter::ncomp_t c=0; c<nprim; ++c)
      {
        auto mark = c*rdof;
        P(e, mark+1, offset) = phip[c] * P(e, mark+1, offset);
        P(e, mark+2, offset) = phip[c] * P(e, mark+2, offset);
        P(e, mark+3, offset) = phip[c] * P(e, mark+3, offset);
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
  std::size_t offset,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& uNodalExtrm,
  tk::Fields& U )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for transport DGP1
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] system Index for equation systems
//! \param[in] offset Index for equation systems
//! \param[in] geoElem Element geometry array
//! \param[in] coord Array of nodal coordinates
//! \param[in] gid Local->global node id map
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!   global node ids (key)
//! \param[in] uNodalExtrm Chare-boundary nodal extrema for conservative
//!   variables
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
      std::vector< std::vector< tk::real > > unk;
      // limit conserved quantities
      auto phi = VertexBasedLimiting(unk, U, esup, inpoel, coord, geoElem, e,
        rdof, dof_el, offset, ncomp, gid, bid, uNodalExtrm);

      // limits under which compression is to be performed
      std::vector< std::size_t > matInt(ncomp, 0);
      std::vector< tk::real > alAvg(ncomp, 0.0);
      for (std::size_t k=0; k<ncomp; ++k)
        alAvg[k] = U(e,k*rdof,offset);
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
        U(e, mark+1, offset) = phi[c] * U(e, mark+1, offset);
        U(e, mark+2, offset) = phi[c] * U(e, mark+2, offset);
        U(e, mark+3, offset) = phi[c] * U(e, mark+3, offset);
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
  std::size_t offset,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& uNodalExtrm,
  tk::Fields& U )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for single-material DGP1
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] offset Index for equation systems
//! \param[in] geoElem Element geometry array
//! \param[in] coord Array of nodal coordinates
//! \param[in] gid Local->global node id map
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!   global node ids (key)
//! \param[in] uNodalExtrm Chare-boundary nodal extrema for conservative
//!   variables
//! \param[in,out] U High-order solution vector which gets limited
//! \details This vertex-based limiter function should be called for compflow.
//!   For details see: Kuzmin, D. (2010). A vertex-based hierarchical slope
//!   limiter for p-adaptive discontinuous Galerkin methods. Journal of
//!   computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
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
      std::vector< std::vector< tk::real > > unk;
      // limit conserved quantities
      auto phi = VertexBasedLimiting(unk, U, esup, inpoel, coord, geoElem,
        e, rdof, dof_el, offset, ncomp, gid, bid, uNodalExtrm);

      // apply limiter function
      for (std::size_t c=0; c<ncomp; ++c)
      {
        auto mark = c*rdof;
        U(e, mark+1, offset) = phi[c] * U(e, mark+1, offset);
        U(e, mark+2, offset) = phi[c] * U(e, mark+2, offset);
        U(e, mark+3, offset) = phi[c] * U(e, mark+3, offset);
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
  std::size_t offset,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& uNodalExtrm,
  tk::Fields& U )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for single-material DGP2
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] offset Index for equation systems
//! \param[in] geoElem Element geometry array
//! \param[in] coord Array of nodal coordinates
//! \param[in] gid Local->global node id map
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!   global node ids (key)
//! \param[in] uNodalExtrm Chare-boundary nodal extrema for conservative
//!   variables
//! \param[in,out] U High-order solution vector which gets limited
//! \details This vertex-based limiter function should be called for compflow.
//!   For details see: Kuzmin, D. (2010). A vertex-based hierarchical slope
//!   limiter for p-adaptive discontinuous Galerkin methods. Journal of
//!   computational and applied mathematics, 233(12), 3077-3085.
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  std::size_t ncomp = U.nprop()/rdof;

  // Copy field data U to U_lim. U_lim will store the limited solution
  // temporarily, so that the limited solution is NOT used to find the
  // min/max bounds for the limiting function
  auto U_lim = U;

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

    bool shock_detec(false);

    // Evaluate the shock detection indicator
    auto Ind = evalDiscontinuityIndicator(e, ncomp, dof_el, ndofel[e], U);
    if(Ind > 1e-6)
      shock_detec = true;

    if (dof_el > 1 && shock_detec == true)
    {
      // Transform the solution with Dubiner basis to Taylor basis so that the
      // limiting function could be applied to physical derivatives in a
      // hierarchical manner
      auto unk =
        tk::DubinerToTaylor(ncomp, offset, e, dof_el, U, inpoel, coord);

      // The vector of limiting coefficients for P1 and P2 coefficients
      std::vector< tk::real > phic_p1(ncomp, 1.0);
      std::vector< tk::real > phic_p2(ncomp, 1.0);

      // If DGP2 is applied, apply the limiter function to the first derivative
      // to obtain the limiting coefficient for P2 coefficients
      if(dof_el > 4)
        phic_p2 = VertexBasedLimiting_P2(unk, U, esup, inpoel, coord, geoElem,
          e, rdof, dof_el, offset, ncomp, gid, bid, uNodalExtrm);

      // limit conserved quantities
      phic_p1 = VertexBasedLimiting(unk, U, esup, inpoel, coord, geoElem, e,
        rdof, dof_el, offset, ncomp, gid, bid, uNodalExtrm);

      if(dof_el > 4)
        for (std::size_t c=0; c<ncomp; ++c)
          phic_p1[c] = std::max(phic_p1[c], phic_p2[c]);

      // apply limiter function to the solution with Taylor basis
      for (std::size_t c=0; c<ncomp; ++c)
      {
        unk[c][1] = phic_p1[c] * unk[c][1];
        unk[c][2] = phic_p1[c] * unk[c][2];
        unk[c][3] = phic_p1[c] * unk[c][3];
      }
      if(dof_el > 4)
      {
        for (std::size_t c=0; c<ncomp; ++c)
        {
          unk[c][4] = phic_p2[c] * unk[c][4];
          unk[c][5] = phic_p2[c] * unk[c][5];
          unk[c][6] = phic_p2[c] * unk[c][6];
          unk[c][7] = phic_p2[c] * unk[c][7];
          unk[c][8] = phic_p2[c] * unk[c][8];
          unk[c][9] = phic_p2[c] * unk[c][9];
        }
      }

      // Convert the solution with Taylor basis to the solution with Dubiner
      // basis
      tk::TaylorToDubiner( ncomp, e, dof_el, inpoel, coord, geoElem, unk );

      // Store the limited solution in U_lim
      for(std::size_t c=0; c<ncomp; ++c)
      {
        auto mark = c*rdof;
        for(std::size_t idof = 1; idof < rdof; idof++)
          U_lim(e, mark+idof, offset) = unk[c][idof];
      }
    }
  }

  // Store the limited solution with Dubiner basis
  for (std::size_t e=0; e<nelem; ++e)
  {
    for (std::size_t c=0; c<ncomp; ++c)
    {
      auto mark = c*rdof;
      U(e, mark+1, offset) = U_lim(e, mark+1, offset);
      U(e, mark+2, offset) = U_lim(e, mark+2, offset);
      U(e, mark+3, offset) = U_lim(e, mark+3, offset);

      if(ndof > 4)
      {
        U(e, mark+4, offset) = U_lim(e, mark+4, offset);
        U(e, mark+5, offset) = U_lim(e, mark+5, offset);
        U(e, mark+6, offset) = U_lim(e, mark+6, offset);
        U(e, mark+7, offset) = U_lim(e, mark+7, offset);
        U(e, mark+8, offset) = U_lim(e, mark+8, offset);
        U(e, mark+9, offset) = U_lim(e, mark+9, offset);
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
  std::size_t offset,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& uNodalExtrm,
  const std::vector< std::vector<tk::real> >& pNodalExtrm,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for multi-material DGP1
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] system Index for equation systems
//! \param[in] offset Offset this PDE system operates from
//! \param[in] geoElem Element geometry array
//! \param[in] coord Array of nodal coordinates
//! \param[in] gid Local->global node id map
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!   global node ids (key)
//! \param[in] uNodalExtrm Chare-boundary nodal extrema for conservative
//!   variables
//! \param[in] pNodalExtrm Chare-boundary nodal extrema for primitive
//!   variables
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
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  const auto intsharp = inciter::g_inputdeck.get< tag::param, tag::multimat,
    tag::intsharp >()[system];
  std::size_t ncomp = U.nprop()/rdof;
  std::size_t nprim = P.nprop()/rdof;

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
      std::vector< std::vector< tk::real > > unk;
      std::vector< std::vector< tk::real > > prim;
      // limit conserved quantities
      auto phic = VertexBasedLimiting(unk, U, esup, inpoel, coord, geoElem, e,
        rdof, dof_el, offset, ncomp, gid, bid, uNodalExtrm);
      // limit primitive quantities
      auto phip = VertexBasedLimiting(unk, P, esup, inpoel, coord, geoElem, e,
        rdof, dof_el, offset, nprim, gid, bid, pNodalExtrm);

      std::vector< tk::real > phic_p2;

      if(ndof > 1 && intsharp == 0)
        BoundPreservingLimiting(nmat, offset, ndof, e, inpoel, coord, U, phic,
          phic_p2);

      // limits under which compression is to be performed
      std::vector< std::size_t > matInt(nmat, 0);
      std::vector< tk::real > alAvg(nmat, 0.0);
      for (std::size_t k=0; k<nmat; ++k)
        alAvg[k] = U(e, volfracDofIdx(nmat,k,rdof,0), offset);
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
        consistentMultiMatLimiting_P1(nmat, offset, rdof, e, U, P, unk, prim,
          phic, phic_p2);
      }

      // apply limiter function
      for (std::size_t c=0; c<ncomp; ++c)
      {
        auto mark = c*rdof;
        U(e, mark+1, offset) = phic[c] * U(e, mark+1, offset);
        U(e, mark+2, offset) = phic[c] * U(e, mark+2, offset);
        U(e, mark+3, offset) = phic[c] * U(e, mark+3, offset);
      }
      for (std::size_t c=0; c<nprim; ++c)
      {
        auto mark = c*rdof;
        P(e, mark+1, offset) = phip[c] * P(e, mark+1, offset);
        P(e, mark+2, offset) = phip[c] * P(e, mark+2, offset);
        P(e, mark+3, offset) = phip[c] * P(e, mark+3, offset);
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
  std::size_t offset,
  const tk::Fields& geoElem,
  const tk::UnsMesh::Coords& coord,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& uNodalExtrm,
  const std::vector< std::vector<tk::real> >& pNodalExtrm,
  tk::Fields& U,
  tk::Fields& P,
  std::size_t nmat )
// *****************************************************************************
//  Kuzmin's vertex-based limiter for multi-material DGP2
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in] nelem Number of elements
//! \param[in] system Index for equation systems
//! \param[in] offset Offset this PDE system operates from
//! \param[in] geoElem Element geometry array
//! \param[in] coord Array of nodal coordinates
//! \param[in] gid Local->global node id map
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!   global node ids (key)
//! \param[in] uNodalExtrm Chare-boundary nodal extrema for conservative
//!   variables
//! \param[in] pNodalExtrm Chare-boundary nodal extrema for primitive
//!   variables
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
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  const auto intsharp = inciter::g_inputdeck.get< tag::param, tag::multimat,
    tag::intsharp >()[system];
  std::size_t ncomp = U.nprop()/rdof;
  std::size_t nprim = P.nprop()/rdof;

  // Copy field data U to U_lim. U_lim will store the limited solution
  // temperally to avoid the limited solution will be used when finding
  // the min/max bounds for the limiting function
  auto U_lim = U;
  auto P_lim = P;

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
      // Transform the solution with Dubiner basis to Taylor basis so that the
      // limiting function could be applied to physical derivatives in a
      // hierarchical manner
      auto unk =
        tk::DubinerToTaylor(ncomp, offset, e, dof_el, U, inpoel, coord);
      auto prim =
        tk::DubinerToTaylor(nprim, offset, e, dof_el, P, inpoel, coord);

      // The vector of limiting coefficients for P1 and P2 coefficients
      std::vector< tk::real > phic_p1(ncomp, 1.0);
      std::vector< tk::real > phic_p2(ncomp, 1.0);
      std::vector< tk::real > phip_p1(nprim, 1.0);
      std::vector< tk::real > phip_p2(nprim, 1.0);

      // If DGP2 is applied, apply the limiter function to the first derivative
      // to obtain the limiting coefficient for P2 coefficients
      if(dof_el > 4)
      {
        phic_p2 = VertexBasedLimiting_P2(unk, U, esup, inpoel, coord, geoElem,
          e, rdof, dof_el, offset, ncomp, gid, bid, uNodalExtrm);
        phip_p2 = VertexBasedLimiting_P2(prim, P, esup, inpoel, coord, geoElem,
          e, rdof, dof_el, offset, nprim, gid, bid, pNodalExtrm);
      }

      // limit conserved quantities
      phic_p1 = VertexBasedLimiting(unk, U, esup, inpoel, coord, geoElem, e,
        rdof, dof_el, offset, ncomp, gid, bid, uNodalExtrm);
      // limit primitive quantities
      phip_p1 = VertexBasedLimiting(prim, P, esup, inpoel, coord, geoElem, e,
        rdof, dof_el, offset, nprim, gid, bid, pNodalExtrm);

      //for (std::size_t c=0; c<ncomp; ++c)
      //  phic_p1[c] = std::max(phic_p1[c], phic_p2[c]);
      //for (std::size_t c=0; c<nprim; ++c)
      //  phip_p1[c] = std::max(phip_p1[c], phip_p2[c]);

      if(ndof > 1 && intsharp == 0)
        BoundPreservingLimiting(nmat, offset, ndof, e, inpoel, coord, U, phic_p1,
          phic_p2);

      PositivityPreservingLimiting(nmat, offset, ndof, e, inpoel, coord, U, phic_p1,
          phic_p2);

      // limits under which compression is to be performed
      std::vector< std::size_t > matInt(nmat, 0);
      std::vector< tk::real > alAvg(nmat, 0.0);
      for (std::size_t k=0; k<nmat; ++k)
        alAvg[k] = U(e, volfracDofIdx(nmat,k,rdof,0), offset);
      auto intInd = interfaceIndicator(nmat, alAvg, matInt);
      if ((intsharp > 0) && intInd)
      {
        for (std::size_t k=0; k<nmat; ++k)
        {
          if (matInt[k])
          {
            phic_p1[volfracIdx(nmat,k)] = 1.0;
            phic_p2[volfracIdx(nmat,k)] = 1.0;
            for(std::size_t idof = 1; idof < rdof; idof++)
            {
              unk[densityIdx(nmat, k)][idof] = 0;
              unk[energyIdx(nmat, k)][idof] = 0;
            }
          }
        }
        for(std::size_t idir = 0; idir < 3; idir++)
          for(std::size_t idof = 1; idof < rdof; idof++)
            unk[momentumIdx(nmat, idir)][idof] = 0;
      }
      else
      {
        consistentMultiMatLimiting_P1(nmat, offset, rdof, e, U, P, unk, prim,
          phic_p1, phic_p2);
      }

      // apply limiter function
      for (std::size_t c=0; c<ncomp; ++c)
      {
        for(std::size_t idof=1; idof<4; idof++)
          unk[c][idof] = phic_p1[c] * unk[c][idof];
        for(std::size_t idof=4; idof<rdof; idof++)
          unk[c][idof] = phic_p2[c] * unk[c][idof];
      }
      for (std::size_t c=0; c<nprim; ++c)
      {
        for(std::size_t idof=1; idof<4; idof++)
          prim[c][idof] = phip_p1[c] * prim[c][idof];
        for(std::size_t idof=4; idof<rdof; idof++)
          prim[c][idof] = phip_p2[c] * prim[c][idof];
      }

      // Convert the solution with Taylor basis to the solution with Dubiner basis
      tk::TaylorToDubiner( ncomp, e, dof_el, inpoel, coord, geoElem, unk );
      tk::TaylorToDubiner( nprim, e, dof_el, inpoel, coord, geoElem, prim );

      // Store the limited solution in U_lim and P_lim
      for(std::size_t c=0; c<ncomp; ++c)
      {
        auto mark = c*rdof;
        for(std::size_t idof = 1; idof < rdof; idof++)
          U_lim(e, mark+idof, offset) = unk[c][idof];
      }
      for(std::size_t c=0; c<nprim; ++c)
      {
        auto mark = c*rdof;
        for(std::size_t idof = 1; idof < rdof; idof++)
          P_lim(e, mark+idof, offset) = prim[c][idof];
      }
    }
  }

  // Store the limited solution with Dubiner basis
  for (std::size_t e=0; e<nelem; ++e)
  {
    for (std::size_t c=0; c<ncomp; ++c)
    {
      auto mark = c*rdof;
      for(std::size_t idof=1; idof<rdof; idof++)
        U(e, mark+idof, offset) = U_lim(e, mark+idof, offset);
    }
    for (std::size_t c=0; c<nprim; ++c)
    {
      auto mark = c*rdof;
      for(std::size_t idof=1; idof<rdof; idof++)
        P(e, mark+idof, offset) = P_lim(e, mark+idof, offset);
    }
  }
}

void
WENOLimiting( const tk::Fields& U,
              const std::vector< int >& esuel,
              std::size_t e,
              inciter::ncomp_t c,
              std::size_t rdof,
              inciter::ncomp_t offset,
              tk::real cweight,
              std::array< std::vector< tk::real >, 3 >& limU )
// *****************************************************************************
//  WENO limiter function calculation for P1 dofs
//! \param[in] U High-order solution vector which is to be limited
//! \param[in] esuel Elements surrounding elements
//! \param[in] e Id of element whose solution is to be limited
//! \param[in] c Index of component which is to be limited
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] offset Index for equation systems
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

std::vector< tk::real >
SuperbeeLimiting( const tk::Fields& U,
                  const std::vector< int >& esuel,
                  const std::vector< std::size_t >& inpoel,
                  const tk::UnsMesh::Coords& coord,
                  std::size_t e,
                  std::size_t ndof,
                  std::size_t rdof,
                  std::size_t dof_el,
                  inciter::ncomp_t offset,
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
//! \param[in] offset Index for equation systems
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
          uNeg = std::max(uNeg, 1.0e-08);
          phi_gp = std::min( 1.0, (uMax[c]-U(e, mark, offset))/(2.0*uNeg) );
        }
        else if (uNeg < -1.0e-14)
        {
          uNeg = std::min(uNeg, -1.0e-08);
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

  return phi;
}

std::vector< tk::real >
VertexBasedLimiting( const std::vector< std::vector< tk::real > >& unk,
  const tk::Fields& U,
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const tk::Fields& geoElem,
  std::size_t e,
  std::size_t rdof,
  std::size_t dof_el,
  std::size_t offset,
  std::size_t ncomp,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& NodalExtrm )
// *****************************************************************************
//  Kuzmin's vertex-based limiter function calculation for P1 dofs
//! \param[in] U High-order solution vector which is to be limited
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] e Id of element whose solution is to be limited
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] dof_el Local number of degrees of freedom
//! \param[in] offset Index for equation systems
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] gid Local->global node id map
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!   global node ids (key)
//! \param[in] NodalExtrm Chare-boundary nodal extrema
//! \return phi Limiter function for solution in element e
// *****************************************************************************
{
  // Kuzmin's vertex-based TVD limiter uses min-max bounds that the
  // high-order solution should satisfy, to ensure TVD properties. For a
  // high-order method like DG, this involves the following steps:
  // 1. Find min-max bounds in the nodal-neighborhood of cell.
  // 2. Calculate the limiter function (Superbee) for all the vertices of cell.
  //    From these, use the minimum value of the limiter function.

  const auto nelem = inpoel.size() / 4;

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

  std::vector< tk::real > uMin(ncomp, 0.0), uMax(ncomp, 0.0), phi(ncomp, 1.0);

  // loop over all nodes of the element e
  for (std::size_t lp=0; lp<4; ++lp)
  {
    // reset min/max
    for (std::size_t c=0; c<ncomp; ++c)
    {
      auto mark = c*rdof;
      uMin[c] = U(e, mark, offset);
      uMax[c] = U(e, mark, offset);
    }
    auto p = inpoel[4*e+lp];
    const auto& pesup = tk::cref_find(esup, p);

    // ----- Step-1: find min/max in the neighborhood of node p
    // loop over all the internal elements surrounding this node p
    for (auto er : pesup)
    {
      if(er < nelem)
      {
        for (std::size_t c=0; c<ncomp; ++c)
        {
          auto mark = c*rdof;
          uMin[c] = std::min(uMin[c], U(er, mark, offset));
          uMax[c] = std::max(uMax[c], U(er, mark, offset));
        }
      }
    }

    // If node p is the chare-boundary node, find min/max by comparing with
    // the chare-boundary nodal extrema from vector NodalExtrm
    auto gip = bid.find( gid[p] );
    if(gip != end(bid))
    {
      auto ndof_NodalExtrm = NodalExtrm[0].size() / (ncomp * 2);
      for (std::size_t c=0; c<ncomp; ++c)
      {
        auto max_mark = 2*c*ndof_NodalExtrm;
        auto min_mark = max_mark + 1;
        uMax[c] = std::max(NodalExtrm[gip->second][max_mark], uMax[c]);
        uMin[c] = std::min(NodalExtrm[gip->second][min_mark], uMin[c]);
      }
    }

    // ----- Step-2: compute the limiter function at this node
    // find high-order solution
    std::vector< tk::real > state( ncomp, 0.0 );
    if(rdof == 4)
    {
      // If DG(P1), evaluate high order solution based on dubiner basis
      std::array< tk::real, 3 > gp{cx[p], cy[p], cz[p]};
      auto B_p = tk::eval_basis( rdof,
            tk::Jacobian( coordel[0], gp, coordel[2], coordel[3] ) / detT,
            tk::Jacobian( coordel[0], coordel[1], gp, coordel[3] ) / detT,
            tk::Jacobian( coordel[0], coordel[1], coordel[2], gp ) / detT );
      state = tk::eval_state( ncomp, offset, rdof, dof_el, e, U, B_p );
    }
    else {  // If DG(P2), evaluate high order solution based on Taylor basis
      // The nodal and central coordinates
      std::array< tk::real, 3 > node{cx[p], cy[p], cz[p]};
      std::array< tk::real, 3 > x_center
        { geoElem(e,1,0), geoElem(e,2,0), geoElem(e,3,0) };
      auto B_p = tk::eval_TaylorBasis( rdof, node, x_center, coordel );

      for (ncomp_t c=0; c<ncomp; ++c)
        for(std::size_t idof = 0; idof < 4; idof++)
          state[c] += unk[c][idof] * B_p[idof];
    }

    Assert( state.size() == ncomp, "Size mismatch" );

    // compute the limiter function
    for (std::size_t c=0; c<ncomp; ++c)
    {
      auto phi_gp = 1.0;
      auto mark = c*rdof;
      auto uNeg = state[c] - U(e, mark, offset);
      auto uref = std::max(std::fabs(U(e,mark,offset)), 1e-14);
      if (uNeg > 1.0e-06*uref)
      {
        phi_gp = std::min( 1.0, (uMax[c]-U(e, mark, offset))/uNeg );
      }
      else if (uNeg < -1.0e-06*uref)
      {
        phi_gp = std::min( 1.0, (uMin[c]-U(e, mark, offset))/uNeg );
      }
      else
      {
        phi_gp = 1.0;
      }

    // ----- Step-3: take the minimum of the nodal-limiter functions
      phi[c] = std::min( phi[c], phi_gp );
    }
  }

  return phi;
}

std::vector< tk::real >
VertexBasedLimiting_P2( const std::vector< std::vector< tk::real > >& unk,
  const tk::Fields& U,
  const std::map< std::size_t, std::vector< std::size_t > >& esup,
  const std::vector< std::size_t >& inpoel,
  const tk::UnsMesh::Coords& coord,
  const tk::Fields& geoElem,
  std::size_t e,
  std::size_t rdof,
  [[maybe_unused]] std::size_t dof_el,
  std::size_t offset,
  std::size_t ncomp,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& bid,
  const std::vector< std::vector<tk::real> >& NodalExtrm )
// *****************************************************************************
//  Kuzmin's vertex-based limiter function calculation for P2 dofs
//! \param[in] U High-order solution vector which is to be limited
//! \param[in] esup Elements surrounding points
//! \param[in] inpoel Element connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] e Id of element whose solution is to be limited
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] dof_el Local number of degrees of freedom
//! \param[in] offset Index for equation systems
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] gid Local->global node id map
//! \param[in] bid Local chare-boundary node ids (value) associated to
//!   global node ids (key)
//! \param[in] NodalExtrm Chare-boundary nodal extrema
//! \return phi Limiter function for solution in element e
//! \details This function limits the P2 dofs of P2 solution in a hierachical
//!   way to P1 dof limiting. Here we treat the first order derivatives the same
//!   way as cell average while second order derivatives represent the gradients
//!   to be limited in the P1 limiting procedure.
// *****************************************************************************
{
  const auto nelem = inpoel.size() / 4;

  // Prepare for calculating Basis functions
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  std::vector< tk::real > phi(ncomp, 1.0);
  std::vector< std::vector< tk::real > > uMin, uMax;
  uMin.resize( ncomp, std::vector<tk::real>(3, 0.0) );
  uMax.resize( ncomp, std::vector<tk::real>(3, 0.0) );

  // The coordinates of centroid in the reference domain
  std::array< std::vector< tk::real >, 3 > center;
  center[0].resize(1, 0.25);
  center[1].resize(1, 0.25);
  center[2].resize(1, 0.25);

  // loop over all nodes of the element e
  for (std::size_t lp=0; lp<4; ++lp)
  {
    // Find the max/min first-order derivatives for internal element
    for (std::size_t c=0; c<ncomp; ++c)
    {
      for (std::size_t idir=1; idir < 4; ++idir)
      {
        uMin[c][idir-1] = unk[c][idir];
        uMax[c][idir-1] = unk[c][idir];
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
        // Coordinates of the neighboring element
        std::array< std::array< tk::real, 3>, 4 > coorder {{
         {{ cx[ inpoel[4*er  ] ], cy[ inpoel[4*er  ] ], cz[ inpoel[4*er  ] ] }},
         {{ cx[ inpoel[4*er+1] ], cy[ inpoel[4*er+1] ], cz[ inpoel[4*er+1] ] }},
         {{ cx[ inpoel[4*er+2] ], cy[ inpoel[4*er+2] ], cz[ inpoel[4*er+2] ] }},
         {{ cx[ inpoel[4*er+3] ], cy[ inpoel[4*er+3] ], cz[ inpoel[4*er+3] ] }} }};

        auto jacInv_er = 
          tk::inverseJacobian( coorder[0], coorder[1], coorder[2], coorder[3] );

        // Compute the derivatives of basis function in the physical domain
        auto dBdx_er = tk::eval_dBdx_p1( rdof, jacInv_er );

        if(rdof > 4)
          tk::eval_dBdx_p2(0, center, jacInv_er, dBdx_er);

        for (std::size_t c=0; c<ncomp; ++c)
        {
          auto mark = c*rdof;
          for (std::size_t idir=0; idir < 3; ++idir)
          {
            // The first order derivative at the centroid of element er
            tk::real slope_er(0.0);
            for(std::size_t idof = 1; idof < rdof; idof++)
              slope_er += U(er, mark+idof, offset) * dBdx_er[idir][idof];

            uMin[c][idir] = std::min(uMin[c][idir], slope_er);
            uMax[c][idir] = std::max(uMax[c][idir], slope_er);

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
      for (std::size_t c=0; c<ncomp; ++c)
      {
        for (std::size_t idir = 1; idir < 4; idir++)
        {
          auto max_mark = 2*c*ndof_NodalExtrm + 2*idir;
          auto min_mark = max_mark + 1;
          uMax[c][idir-1] =
            std::max(NodalExtrm[gip->second][max_mark], uMax[c][idir-1]);
          uMin[c][idir-1] =
            std::min(NodalExtrm[gip->second][min_mark], uMin[c][idir-1]);
        }
      }
    }

    //Step-2: compute the limiter function at this node
    std::array< tk::real, 3 > node{cx[p], cy[p], cz[p]};
    std::array< tk::real, 3 >
      centroid_physical{geoElem(e,1,0), geoElem(e,2,0), geoElem(e,3,0)};

    // find high-order solution
    std::vector< std::vector< tk::real > > state;
    state.resize( ncomp, std::vector<tk::real>(3, 0.0) );

    for (ncomp_t c=0; c<ncomp; ++c)
    {
      auto dx = node[0] - centroid_physical[0];
      auto dy = node[1] - centroid_physical[1];
      auto dz = node[2] - centroid_physical[2];

      state[c][0] = unk[c][1] + unk[c][4] * dx + unk[c][7] * dy + unk[c][8] * dz;
      state[c][1] = unk[c][2] + unk[c][5] * dy + unk[c][7] * dx + unk[c][9] * dz;
      state[c][2] = unk[c][3] + unk[c][6] * dz + unk[c][8] * dx + unk[c][9] * dy;
    }

    // compute the limiter function
    for (std::size_t c=0; c<ncomp; ++c)
    {
      tk::real phi_dir(1.0);
      for (std::size_t idir=1; idir < 3; ++idir)
      {
        phi_dir = 1.0;
        auto uNeg = state[c][idir-1] - unk[c][idir];
        auto uref = std::max(std::fabs(unk[c][idir]), 1e-14);
        if (uNeg > 1.0e-6*uref)
        {
          phi_dir =
            std::min( 1.0, ( uMax[c][idir-1] - unk[c][idir])/uNeg );
        }
        else if (uNeg < -1.0e-6*uref)
        {
          phi_dir =
            std::min( 1.0, ( uMin[c][idir-1] - unk[c][idir])/uNeg );
        }
        else
        {
          phi_dir = 1.0;
        }

        phi[c] = std::min( phi[c], phi_dir );
      }
    }
  }

  return phi;
}


void consistentMultiMatLimiting_P1(
  const std::size_t nmat,
  const ncomp_t offset,
  const std::size_t rdof,
  const std::size_t e,
  tk::Fields& U,
  [[maybe_unused]] tk::Fields& P,
  std::vector< std::vector< tk::real > >& unk,
  [[maybe_unused]] std::vector< std::vector< tk::real > >& prim,
  std::vector< tk::real >& phic_p1, 
  std::vector< tk::real >& phic_p2 )
// *****************************************************************************
//  Consistent limiter modifications for conservative variables
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] offset Index for equation system
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] e Element being checked for consistency
//! \param[in] U Vector of conservative variables
//! \param[in] P Vector of primitive variables
//! \parame[in,out] unk Second-order solution vector based on Taylor basis which
//!   gets modified near material interfaces for consistency
//! \parame[in,out] prim Second-order solution vector of primitive variables
//!   based on Taylor basis which gets modified near material interfaces for
//!   consistency
//! \param[in,out] phic_p1 Vector of limiter functions for P1 dofs of the
//!   conserved quantities
//! \param[in,out] phip_p2 Vector of limiter functions for P2 dofs of the
//!   conserved quantities
// *****************************************************************************
{
  Assert(phic_p1.size()+phic_p2.size() == U.nprop()/rdof, "Number of unknowns "
    "in vector of conserved quantities incorrect");

  // find the limiter-function for volume-fractions
  auto phi_al_p1(1.0), phi_al_p2(1.0), almax(0.0), dalmax(0.0);
  //std::size_t nmax(0);
  for (std::size_t k=0; k<nmat; ++k)
  {
    phi_al_p1 = std::min( phi_al_p1, phic_p1[volfracIdx(nmat, k)] );
    if(rdof > 4)
      phi_al_p2 = std::min( phi_al_p2, phic_p2[volfracIdx(nmat, k)] );
    if (almax < U(e,volfracDofIdx(nmat, k, rdof, 0),offset))
    {
      //nmax = k;
      almax = U(e,volfracDofIdx(nmat, k, rdof, 0),offset);
    }
    tk::real dmax(0.0);
    if(rdof > 4)
      dmax = std::max(
               std::max( std::abs(unk[volfracIdx(nmat, k)][1]),
                         std::abs(unk[volfracIdx(nmat, k)][2]) ),
                         std::abs(unk[volfracIdx(nmat, k)][3]) );
    else
      dmax = std::max(
               std::max(
                 std::abs(U(e,volfracDofIdx(nmat, k, rdof, 1),offset)),
                 std::abs(U(e,volfracDofIdx(nmat, k, rdof, 2),offset)) ),
                 std::abs(U(e,volfracDofIdx(nmat, k, rdof, 3),offset)) );
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
  if ( dalmax > al_band &&
       (almax > al_band && almax < (1.0-al_band)) )
  {
    // 1. consistent high-order dofs
    for (std::size_t k=0; k<nmat; ++k)
    {
      auto alk =
        std::max( 1.0e-14, U(e,volfracDofIdx(nmat, k, rdof, 0),offset) );
      auto rhok = U(e,densityDofIdx(nmat, k, rdof, 0),offset)/alk;
      for (std::size_t idof=1; idof<rdof; ++idof)
      {
        // If DG(P2), modify the unk vector, otherwise the modification is
        // applied to solution vector U
        if(rdof > 4)
          unk[densityIdx(nmat, k)][idof] =
            rhok * unk[volfracIdx(nmat, k)][idof];
        else
          U(e,densityDofIdx(nmat, k, rdof, idof),offset) = rhok *
            U(e,volfracDofIdx(nmat, k, rdof, idof),offset);
      }
    }

    // 2. same limiter for all volume-fractions and densities
    for (std::size_t k=0; k<nmat; ++k)
    {
      phic_p1[volfracIdx(nmat, k)] = phi_al_p1;
      phic_p1[densityIdx(nmat, k)] = phi_al_p1;
    }
    if(rdof > 4)
    {
      for (std::size_t k=0; k<nmat; ++k)
      {
        phic_p2[volfracIdx(nmat, k)] = phi_al_p2;
        phic_p2[densityIdx(nmat, k)] = phi_al_p2;
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
                              ncomp_t offset,
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
//! \param[in] offset Index for equation system
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

      auto state = eval_state( U.nprop()/ndof, offset, ndof, ndof, e, U, B );

      for(std::size_t imat = 0; imat < nmat; imat++)
      {
        auto phi = BoundPreservingLimitingFunction( min, max,
          state[volfracIdx(nmat, imat)],
          U(e,volfracDofIdx(nmat, imat, ndof, 0),offset) );
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

      auto state = tk::eval_state(U.nprop()/ndof, offset, ndof, ndof, e, U, B);

      for(std::size_t imat = 0; imat < nmat; imat++)
      {
        auto phi = BoundPreservingLimitingFunction(min, max,
          state[volfracIdx(nmat, imat)],
          U(e,volfracDofIdx(nmat, imat, ndof, 0),offset) );
        phi_bound[imat] = std::min( phi_bound[imat], phi );
      }
    }
  }

  for(std::size_t imat = 0; imat < nmat; imat++)
    phic_p1[imat] = phi_bound[imat] * phic_p1[imat];
  if(ndof > 4)
    for(std::size_t imat = 0; imat < nmat; imat++)
      phic_p2[imat] = phi_bound[imat] * phic_p2[imat];
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
  tk::real phi(1.0);
  if(al_gp > max)
    phi = std::fabs( (max - al_avg) / (al_gp - al_avg) );
  else if(al_gp < min)
    phi = std::fabs( (min - al_avg) / (al_gp - al_avg) );
  return phi;
}

void PositivityPreservingLimiting( std::size_t nmat,
                                   ncomp_t offset,
                                   std::size_t ndof,
                                   std::size_t e,
                                   const std::vector< std::size_t >& inpoel,
                                   const tk::UnsMesh::Coords& coord,
                                   const tk::Fields& U,
                                   std::vector< tk::real >& phic_p1,
                                   std::vector< tk::real >& phic_p2 )
// *****************************************************************************
//  Positivity preserving limiter for density when MulMat DG(P2) scheme is
//    selected
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] offset Index for equation system
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

    auto ng = tk::NGfa(ndof);

    std::array< std::vector< tk::real >, 2 > coordgp;
    std::vector< tk::real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    wgp.resize( ng );

    tk::GaussQuadratureTri( ng, coordgp, wgp );

    for (std::size_t igp=0; igp<ng; ++igp)
    {
      auto gp = tk::eval_gp( igp, coordfa, coordgp );
      auto B = tk::eval_basis( ndof,
            tk::Jacobian( coordel[0], gp, coordel[2], coordel[3] ) / detT,
            tk::Jacobian( coordel[0], coordel[1], gp, coordel[3] ) / detT,
            tk::Jacobian( coordel[0], coordel[1], coordel[2], gp ) / detT );

      auto state = eval_state( U.nprop()/ndof, offset, ndof, ndof, e, U, B );

      for(std::size_t imat = 0; imat < nmat; imat++)
      {
        tk::real phi(1.0);
        auto al = state[densityIdx(nmat, imat)];
        auto al_avg = U(e, densityDofIdx(nmat, imat, ndof, 0), offset);
        if(al < min)
          phi = std::fabs( (min - al_avg) / (al  - al_avg) );
        phi_bound[imat] = std::min( phi_bound[imat], phi );
      }
    }
  }

  if(ndof > 4)
  {
    auto ng = tk::NGvol(ndof);
    std::array< std::vector< tk::real >, 3 > coordgp;
    std::vector< tk::real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    coordgp[2].resize( ng );
    wgp.resize( ng );

    tk::GaussQuadratureTet( ng, coordgp, wgp );

    for (std::size_t igp=0; igp<ng; ++igp)
    {
      auto B = tk::eval_basis( ndof, coordgp[0][igp], coordgp[1][igp],
        coordgp[2][igp] );

      auto state = tk::eval_state(U.nprop()/ndof, offset, ndof, ndof, e, U, B);

      for(std::size_t imat = 0; imat < nmat; imat++)
      {
        tk::real phi(1.0);
        auto al = state[densityIdx(nmat, imat)];
        auto al_avg = U(e, densityDofIdx(nmat, imat, ndof, 0), offset);
        if(al < min)
          phi = std::fabs( (min - al_avg) / (al  - al_avg) );
        phi_bound[imat] = std::min( phi_bound[imat], phi );
      }
    }
  }
  for(std::size_t imat = 0; imat < nmat; imat++)
    phic_p1[densityIdx(nmat,imat)] =
      phi_bound[imat] * phic_p1[densityIdx(nmat,imat)];
  if(ndof > 4)
    for(std::size_t imat = 0; imat < nmat; imat++)
      phic_p2[densityIdx(nmat,imat)] =
        phi_bound[imat] * phic_p2[densityIdx(nmat,imat)];
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

} // inciter::
