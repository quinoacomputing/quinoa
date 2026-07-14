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
  std::size_t rdof ) :
  m_nspec( nspec ),
  m_rdof( rdof )
// *****************************************************************************
//! Constructor
//! \param[in] nspec Number of species in this PDE system
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
// *****************************************************************************
{}

std::vector< tk::real >
MultiSpeciesViscousTermsP0P1::stateAt(
  const std::vector< inciter::EOS >& mat_blk,
  const Fields& U,
  const Fields& P,
  std::size_t e,
  std::size_t ndof,
  const std::vector< tk::real >& B ) const
// *****************************************************************************
//! \brief Reconstruct conserved variables and primitives at a face point
//! \param[in] mat_blk Material EOS block
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] e Element id
//! \param[in] ndof Number of local degrees of freedom
//! \param[in] B Basis functions evaluated at the face point
//! \return State vector containing conserved variables followed by
//!   multispecies primitives
// *****************************************************************************
{
  const auto ncomp = U.nprop()/m_rdof;
  const auto nprim = P.nprop()/m_rdof;
  std::vector< tk::real > state( ncomp+nprim, 0.0 );

  eval_state( ncomp, m_rdof, ndof, e, U, B, state.data() );
  eval_state( nprim, m_rdof, ndof, e, P, B, state.data()+ncomp );
  enforcePhysicalConstraints( mat_blk, ncomp, 1, state );

  return state;
}

void
MultiSpeciesViscousTermsP0P1::gradientIntElem(
  const Fields& U,
  const Fields& P,
  std::size_t elem,
  const std::array< std::vector< tk::real >, 3 >& dBdx,
  std::array< std::array< tk::real, 3 >, 4>& grad ) const
// *****************************************************************************
//! Compute gradients of quantities for an interior element
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] elem element id
//! \param[in] dBdx Basis gradients in element
//! \param[in,out] grad Velocity and temperature gradients in element
// *****************************************************************************
{
  using inciter::multispecies::momentumDofIdx;
  using inciter::multispecies::densityDofIdx;
  using inciter::multispecies::temperatureDofIdx;

  // d(\rho u)/dx
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j) {
      grad[i][j] = //drudx[i][j]
        dBdx[j][1] * U(elem, momentumDofIdx(m_nspec,i,m_rdof,1))
      + dBdx[j][2] * U(elem, momentumDofIdx(m_nspec,i,m_rdof,2))
      + dBdx[j][3] * U(elem, momentumDofIdx(m_nspec,i,m_rdof,3));
    }
  }

  // d(\rho)/dx
  std::array< real, 3 > drdx{{ 0.0, 0.0, 0.0 }};
  real cellAvgRho(0.0);
  for (std::size_t k=0; k<m_nspec; ++k) {
    for (std::size_t j=0; j<3; ++j) {
      drdx[j] +=
        dBdx[j][1] * U(elem, densityDofIdx(m_nspec,k,m_rdof,1))
      + dBdx[j][2] * U(elem, densityDofIdx(m_nspec,k,m_rdof,2))
      + dBdx[j][3] * U(elem, densityDofIdx(m_nspec,k,m_rdof,3));
    }
    cellAvgRho += U(elem, densityDofIdx(m_nspec,k,m_rdof,0));
  }

  std::array< real, 3 > cellAvgVel;
  for (std::size_t i=0; i<3; ++i) {
    cellAvgVel[i] = U(elem, momentumDofIdx(m_nspec,i,m_rdof,0))
      / cellAvgRho;
  }

  // du/dx
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j) {
      // dudx[i][j]
      grad[i][j] -= cellAvgVel[i] * drdx[j];
      grad[i][j] /= cellAvgRho;
    }
  }

  // dT/dx_j
  for (std::size_t j=0; j<3; ++j) {
    // dTdx[j]
    grad[3][j] = dBdx[j][1] * P(elem, temperatureDofIdx(m_nspec,0,m_rdof,1))
               + dBdx[j][2] * P(elem, temperatureDofIdx(m_nspec,0,m_rdof,2))
               + dBdx[j][3] * P(elem, temperatureDofIdx(m_nspec,0,m_rdof,3));
  }
}

void
MultiSpeciesViscousTermsP0P1::interiorFlux(
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t ncomp,
  const std::array< std::vector< tk::real >, 2 >& state,
  const std::array< std::vector< tk::real >, 2 >& cellAvgState,
  const std::array< tk::real, 3 >& fn,
  const std::array< std::array< tk::real, 3 >, 2 >& centroid,
  const std::array< std::array< std::array< tk::real, 3 >, 4>, 2 >& grad,
  std::vector< tk::real >& fl ) const
// *****************************************************************************
//! \brief Compute the multispecies viscous flux at an interior face
//! \param[in] mat_blk Material EOS block
//! \param[in] ncomp Number of components in this system
//! \param[in] state Left and right high-order face states
//! \param[in] cellAvgState Left and right cell-average states
//! \param[in] fn Face normal
//! \param[in] centroid Left and right element centroids
//! \param[in] grad Velocity and temperature gradients in left and right elements
//! \param[in,out] fl Numerical viscous flux using the Modified Gradient approach
//! \details The average gradient is modified according to Weiss et al. to
//!   obtain a stable discretization (average results in unstable central
//!   central difference).
//!   Ref: Weiss, J. M., Maruszewski, J. P., & Smith, W. A. (1999). Implicit
//!   solution of preconditioned Navier-Stokes equations using algebraic
//!   multigrid. AIAA journal, 37(1), 29-36.
// *****************************************************************************
{
  using inciter::multispecies::densityIdx;
  using inciter::multispecies::temperatureIdx;
  using inciter::multispecies::momentumIdx;
  using inciter::multispecies::energyIdx;

  Assert( fl.size() == ncomp, "Incorrect viscous flux vector size" );
  for (auto& f : fl) f = 0.0;

  const auto& ul = state[0];
  const auto& ur = state[1];
  const auto& ucl = cellAvgState[0];
  const auto& ucr = cellAvgState[1];

  // 1. Compute average gradients
  // ---------------------------------------------------------------------------
  // du_i/dx_j
  std::array< std::array< real, 3 >, 3 > dudx;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      dudx[i][j] = 0.5 * ( grad[0][i][j] + grad[1][i][j] );

  // dT/dx_j
  std::array< real, 3 > dTdx{{0, 0, 0}};
  for (std::size_t j=0; j<3; ++j)
    dTdx[j] = 0.5 * ( grad[0][3][j] + grad[1][3][j] );

  // 2. Modify the average gradients
  // ---------------------------------------------------------------------------
  std::array< real, 3 > r_f{{
    centroid[1][0] - centroid[0][0],
    centroid[1][1] - centroid[0][1],
    centroid[1][2] - centroid[0][2] }};
  real r_mag = std::sqrt(tk::dot(r_f, r_f));
  for (std::size_t i=0; i<3; ++i)
    r_f[i] /= r_mag;

  // cell averages needed for modifications
  std::array< real, 2 > cellAvgRho{0, 0};
  for (std::size_t k=0; k<m_nspec; ++k) {
    cellAvgRho[0] += ucl[densityIdx(m_nspec,k)];
    cellAvgRho[1] += ucr[densityIdx(m_nspec,k)];
  }
  std::array< std::array< real, 3 >, 2 > cellAvgVel;
  for (std::size_t i=0; i<3; ++i) {
    cellAvgVel[0][i] = ucl[momentumIdx(m_nspec,i)] / cellAvgRho[0];
    cellAvgVel[1][i] = ucr[momentumIdx(m_nspec,i)] / cellAvgRho[1];
  }

  // modify du_i/dx_j
  auto dudx_m = dudx;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      dudx_m[i][j] -= ( tk::dot(dudx[i], r_f) -
        (cellAvgVel[1][i] - cellAvgVel[0][i])/r_mag ) * r_f[j];

  // modified dT/dx_j
  auto dTdx_m = dTdx;
  for (std::size_t j=0; j<3; ++j)
    dTdx_m[j] -= ( tk::dot(dTdx, r_f) -
      (ucr[ncomp+temperatureIdx(m_nspec,0)]
     - ucl[ncomp+temperatureIdx(m_nspec,0)])/r_mag ) * r_f[j];

  // 3. Compute viscous stress tensor
  // ---------------------------------------------------------------------------
  // establish mixture state
  inciter::Mixture mix_l(m_nspec, ul, mat_blk);
  inciter::Mixture mix_r(m_nspec, ur, mat_blk);

  auto mix_temp_l = ul[ncomp + temperatureIdx(m_nspec,0)];
  auto mix_temp_r = ur[ncomp + temperatureIdx(m_nspec,0)];

  // compute fluid properties (viscosity and conductivity)
  auto mu =
    0.5 * (mix_l.viscCoeff(mix_temp_l, mat_blk) + mix_r.viscCoeff(mix_temp_r, mat_blk));
  auto Cp = 0.5 * (mix_l.Cp(ul[ncomp+temperatureIdx(m_nspec,0)], mat_blk)
    + mix_r.Cp(ur[ncomp+temperatureIdx(m_nspec,0)], mat_blk));
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
}


MultiSpeciesViscousTermsDGP1::MultiSpeciesViscousTermsDGP1(
  std::size_t nspec,
  std::size_t rdof ) :
  m_nspec( nspec ),
  m_rdof( rdof )
// *****************************************************************************
//! Constructor
//! \param[in] nspec Number of species in this PDE system
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
// *****************************************************************************
{}

// DDG Functions will go here
std::vector< tk::real >
MultiSpeciesViscousTermsDGP1::stateAt(
  const std::vector< inciter::EOS >& mat_blk,
  const Fields& U,
  const Fields& P,
  std::size_t e,
  std::size_t ndof,
  const std::vector< tk::real >& B ) const
// *****************************************************************************
//! \brief Reconstruct conserved variables and primitives at a face point
//! \param[in] mat_blk Material EOS block
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] e Element id
//! \param[in] ndof Number of local degrees of freedom
//! \param[in] B Basis functions evaluated at the face point
//! \return State vector containing conserved variables followed by
//!   multispecies primitives
// *****************************************************************************
{
  const auto ncomp = U.nprop()/m_rdof;
  const auto nprim = P.nprop()/m_rdof;
  std::vector< tk::real > state( ncomp+nprim, 0.0 );

  eval_state( ncomp, m_rdof, ndof, e, U, B, state.data() );
  eval_state( nprim, m_rdof, ndof, e, P, B, state.data()+ncomp );
  enforcePhysicalConstraints( mat_blk, ncomp, 1, state );

  return state;
}

void
MultiSpeciesViscousTermsDGP1::gradientIntElem(
  const Fields& U,
  const Fields& P,
  std::size_t elem,
  const std::array< std::vector< tk::real >, 3 >& dBdx,
  std::array< std::array< tk::real, 3 >, 4>& grad ) const
// *****************************************************************************
//! Compute gradients of quantities for an interior element
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] elem element id
//! \param[in] dBdx Basis gradients in element
//! \param[in,out] grad Velocity and temperature gradients in element
// *****************************************************************************
{ 
  // no-op
}


void
MultiSpeciesViscousTermsDGP1::interiorFlux(
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t ncomp,
  const std::array< std::vector< tk::real >, 2 >& state,
  const std::array< std::vector< tk::real >, 2 >& cellAvgState,
  const std::array< tk::real, 3 >& fn,
  const std::array< std::array< tk::real, 3 >, 2 >& centroid,
  const std::array< std::array< std::array< tk::real, 3 >, 4>, 2 >& grad,
  std::vector< tk::real >& fl ) const
// *****************************************************************************
//! \brief Compute the multispecies viscous flux at an interior face
//! \param[in] mat_blk Material EOS block
//! \param[in] ncomp Number of components in this system
//! \param[in] state Left and right high-order face states
//! \param[in] cellAvgState Left and right cell-average states
//! \param[in] fn Face normal
//! \param[in] centroid Left and right element centroids
//! \param[in] grad Velocity and temperature gradients in left and right elements
//! \param[in,out] fl Numerical viscous flux
//! \details Will be implemented using DDG
// *****************************************************************************
{
  // no-op
}

} // tk::
