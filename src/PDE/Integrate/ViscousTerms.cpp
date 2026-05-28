// *****************************************************************************
/*!
  \file      src/PDE/Integrate/ViscousTerms.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing viscous surface terms of a system
     of PDEs in DG methods
  \details   This file contains functionality for computing viscous surface
     terms of a system of PDEs used in discontinuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************

#include "ViscousTerms.hpp"
#include "Basis.hpp"
#include "Reconstruction.hpp"
#include "EoS/GetMatProp.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"
#include "MultiSpecies/Mixture/Mixture.hpp"

namespace tk {

MultiSpeciesViscousTermsP0P1::MultiSpeciesViscousTermsP0P1(
  std::size_t nspec,
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t rdof,
  const Fields& U,
  const Fields& P ) :
  m_nspec( nspec ),
  m_mat_blk( mat_blk ),
  m_rdof( rdof ),
  m_U( U ),
  m_P( P )
// *****************************************************************************
//! Constructor
//! \param[in] nspec Number of species in this PDE system
//! \param[in] mat_blk Material EOS block
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
// *****************************************************************************
{}

std::vector< tk::real >
MultiSpeciesViscousTermsP0P1::stateAt( std::size_t e,
  std::size_t ndof,
  const std::vector< tk::real >& B ) const
// *****************************************************************************
//! \brief Reconstruct conserved variables and primitives at a face point
//! \param[in] e Element id
//! \param[in] ndof Number of local degrees of freedom
//! \param[in] B Basis functions evaluated at the face point
//! \return State vector containing conserved variables followed by
//!   multispecies primitives
// *****************************************************************************
{
  const auto ncomp = m_U.nprop()/m_rdof;
  const auto nprim = m_P.nprop()/m_rdof;
  std::vector< tk::real > state( ncomp+nprim, 0.0 );

  eval_state( ncomp, m_rdof, ndof, e, m_U, B, state.data() );
  eval_state( nprim, m_rdof, ndof, e, m_P, B, state.data()+ncomp );
  enforcePhysicalConstraints( m_mat_blk, ncomp, 1, state );

  return state;
}

std::vector< tk::real >
MultiSpeciesViscousTermsP0P1::interiorFlux(
  const std::array< std::size_t, 2 >& elem,
  const std::array< tk::real, 3 >& fn,
  const std::array< tk::real, 3 >&,
  const std::array< std::array< tk::real, 3 >, 2 >&,
  const std::array< std::array< tk::real, 3 >, 2 >& centroid,
  const std::array< std::vector< tk::real >, 2 >& B,
  const std::array< std::array< std::vector< tk::real >, 3 >, 2 >& dBdx ) const
// *****************************************************************************
//! \brief Compute the multispecies viscous flux at an interior face
//! \param[in] elem Left and right element ids
//! \param[in] fn Face normal
//! \param[in] centroid Left and right element centroids
//! \param[in] B Basis values at the face point
//! \param[in] dBdx Basis gradients in left and right elements
//! \return Numerical viscous flux using the Modified Gradient approach.
//! \details The average gradient is modified according to Weiss et al. to
//!   obtain a stable discretization (average results in unstable central
//!   central difference).
//!   Ref: Weiss, J. M., Maruszewski, J. P., & Smith, W. A. (1999). Implicit
//!   solution of preconditioned Navier-Stokes equations using algebraic
//!   multigrid. AIAA journal, 37(1), 29-36.
// *****************************************************************************
{
  using inciter::multispecies::momentumDofIdx;
  using inciter::multispecies::densityDofIdx;
  using inciter::multispecies::temperatureDofIdx;
  using inciter::multispecies::temperatureIdx;
  using inciter::multispecies::momentumIdx;
  using inciter::multispecies::energyIdx;

  const auto ncomp = m_U.nprop()/m_rdof;
  std::vector< tk::real > fl( ncomp, 0.0 );

  auto ul = stateAt( elem[0], m_rdof, B[0] );
  auto ur = stateAt( elem[1], m_rdof, B[1] );

  // 1. Compute average gradients
  // ---------------------------------------------------------------------------
  const auto& dBdx_l = dBdx[0];
  const auto& dBdx_r = dBdx[1];

  // d(\rho u)/dx
  std::array< std::array< real, 3 >, 3 > dudx_l, dudx_r;
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j) {
      dudx_l[i][j] =
        dBdx_l[j][1] * m_U(elem[0], momentumDofIdx(m_nspec,i,m_rdof,1))
      + dBdx_l[j][2] * m_U(elem[0], momentumDofIdx(m_nspec,i,m_rdof,2))
      + dBdx_l[j][3] * m_U(elem[0], momentumDofIdx(m_nspec,i,m_rdof,3));
      dudx_r[i][j] =
        dBdx_r[j][1] * m_U(elem[1], momentumDofIdx(m_nspec,i,m_rdof,1))
      + dBdx_r[j][2] * m_U(elem[1], momentumDofIdx(m_nspec,i,m_rdof,2))
      + dBdx_r[j][3] * m_U(elem[1], momentumDofIdx(m_nspec,i,m_rdof,3));
    }
  }

  // d(\rho)/dx
  std::array< real, 3 > drdx_l{{ 0.0, 0.0, 0.0 }}, drdx_r{{ 0.0, 0.0, 0.0 }};
  std::array< real, 2 > cellAvgRho{{ 0.0, 0.0 }};
  for (std::size_t k=0; k<m_nspec; ++k) {
    for (std::size_t j=0; j<3; ++j) {
      drdx_l[j] +=
        dBdx_l[j][1] * m_U(elem[0], densityDofIdx(m_nspec,k,m_rdof,1))
      + dBdx_l[j][2] * m_U(elem[0], densityDofIdx(m_nspec,k,m_rdof,2))
      + dBdx_l[j][3] * m_U(elem[0], densityDofIdx(m_nspec,k,m_rdof,3));

      drdx_r[j] +=
        dBdx_r[j][1] * m_U(elem[1], densityDofIdx(m_nspec,k,m_rdof,1))
      + dBdx_r[j][2] * m_U(elem[1], densityDofIdx(m_nspec,k,m_rdof,2))
      + dBdx_r[j][3] * m_U(elem[1], densityDofIdx(m_nspec,k,m_rdof,3));
    }
    cellAvgRho[0] += m_U(elem[0], densityDofIdx(m_nspec,k,m_rdof,0));
    cellAvgRho[1] += m_U(elem[1], densityDofIdx(m_nspec,k,m_rdof,0));
  }

  std::array< std::array< real, 3 >, 2 > cellAvgVel;
  for (std::size_t i=0; i<3; ++i) {
    cellAvgVel[0][i] = m_U(elem[0], momentumDofIdx(m_nspec,i,m_rdof,0))
      / cellAvgRho[0];
    cellAvgVel[1][i] = m_U(elem[1], momentumDofIdx(m_nspec,i,m_rdof,0))
      / cellAvgRho[1];
  }

  // du/dx
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j) {
      dudx_l[i][j] -= cellAvgVel[0][i] * drdx_l[j];
      dudx_l[i][j] /= cellAvgRho[0];

      dudx_r[i][j] -= cellAvgVel[1][i] * drdx_r[j];
      dudx_r[i][j] /= cellAvgRho[1];
    }
  }

  // 2. Modify the average gradients
  // ---------------------------------------------------------------------------
  std::array< real, 3 > r_f{{
    centroid[1][0] - centroid[0][0],
    centroid[1][1] - centroid[0][1],
    centroid[1][2] - centroid[0][2] }};
  real r_mag = std::sqrt(tk::dot(r_f, r_f));
  for (std::size_t i=0; i<3; ++i)
    r_f[i] /= r_mag;

  // average du_i/dx_j
  std::array< std::array< real, 3 >, 3 > dudx;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      dudx[i][j] = 0.5 * ( dudx_l[i][j] + dudx_r[i][j] );

  // modify du_i/dx_j
  auto dudx_m = dudx;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      dudx_m[i][j] = 0.5 * ( dudx_l[i][j] + dudx_r[i][j] )
        - ( tk::dot(dudx[i], r_f) -
        (cellAvgVel[1][i] - cellAvgVel[0][i])/r_mag ) * r_f[j];

  // average dT/dx_j
  std::array< real, 3 > dTdx{{0, 0, 0}};
  for (std::size_t j=0; j<3; ++j) {
    dTdx[j] = 0.5 * (
        dBdx_l[j][1] * m_P(elem[0], temperatureDofIdx(m_nspec,0,m_rdof,1))
      + dBdx_l[j][2] * m_P(elem[0], temperatureDofIdx(m_nspec,0,m_rdof,2))
      + dBdx_l[j][3] * m_P(elem[0], temperatureDofIdx(m_nspec,0,m_rdof,3))
      + dBdx_r[j][1] * m_P(elem[1], temperatureDofIdx(m_nspec,0,m_rdof,1))
      + dBdx_r[j][2] * m_P(elem[1], temperatureDofIdx(m_nspec,0,m_rdof,2))
      + dBdx_r[j][3] * m_P(elem[1], temperatureDofIdx(m_nspec,0,m_rdof,3)) );
  }

  // modified dT/dx_j
  auto dTdx_m = dTdx;
  for (std::size_t j=0; j<3; ++j)
    dTdx_m[j] -= ( tk::dot(dTdx, r_f) -
      (m_P(elem[1], temperatureDofIdx(m_nspec,0,m_rdof,0))
     - m_P(elem[0], temperatureDofIdx(m_nspec,0,m_rdof,0)))/r_mag ) * r_f[j];

  // 3. Compute viscous stress tensor
  // ---------------------------------------------------------------------------
  // establish mixture state
  inciter::Mixture mix_l(m_nspec, ul, m_mat_blk);
  inciter::Mixture mix_r(m_nspec, ur, m_mat_blk);

  // compute fluid properties (viscosity and conductivity)
  auto mu =
    0.5 * (mix_l.viscCoeff(m_mat_blk) + mix_r.viscCoeff(m_mat_blk));
  auto Cp = 0.5 * (mix_l.Cp(ul[ncomp+temperatureIdx(m_nspec,0)], m_mat_blk)
    + mix_r.Cp(ur[ncomp+temperatureIdx(m_nspec,0)], m_mat_blk));
  auto kTh = mu * Cp / 0.71; // TODO: make Prandtl number user-configurable

  // stress tensor
  std::array< real, 6 > tau;
  tau[0] = mu * ( 4.0 * dudx_m[0][0] - 2.0*(dudx_m[1][1] + dudx_m[2][2]) ) / 3.0;
  tau[1] = mu * ( 4.0 * dudx_m[1][1] - 2.0*(dudx_m[0][0] + dudx_m[2][2]) ) / 3.0;
  tau[2] = mu * ( 4.0 * dudx_m[2][2] - 2.0*(dudx_m[0][0] + dudx_m[1][1]) ) / 3.0;
  tau[3] = mu * ( dudx_m[0][1] + dudx_m[1][0] );
  tau[4] = mu * ( dudx_m[0][2] + dudx_m[2][0] );
  tau[5] = mu * ( dudx_m[1][2] + dudx_m[2][1] );

  // 3. Compute viscous flux terms
  // ---------------------------------------------------------------------------
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      fl[momentumIdx(m_nspec, i)] += tau[inciter::stressCmp[i][j]] * fn[j];

  auto rho_l = mix_l.get_mix_density();
  auto rho_r = mix_r.get_mix_density();
  std::array< real, 3 > energyFlux{{0.0, 0.0, 0.0}};
  std::array< real, 3 > uAvg{{
    0.5 * (ul[momentumIdx(m_nspec,0)]/rho_l + ur[momentumIdx(m_nspec,0)]/rho_r),
    0.5 * (ul[momentumIdx(m_nspec,1)]/rho_l + ur[momentumIdx(m_nspec,1)]/rho_r),
    0.5 * (ul[momentumIdx(m_nspec,2)]/rho_l + ur[momentumIdx(m_nspec,2)]/rho_r)
    }};

  for (std::size_t j=0; j<3; ++j) {
    energyFlux[j] = kTh * dTdx_m[j];
    for (std::size_t i=0; i<3; ++i)
      energyFlux[j] += uAvg[i] * tau[inciter::stressCmp[i][j]];
  }

  fl[energyIdx(m_nspec, 0)] = tk::dot(energyFlux, fn);

  return fl;
}

} // tk::
