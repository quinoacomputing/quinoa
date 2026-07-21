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
  const std::array< std::vector< tk::real >, 6 >& d2Bdx2,
  std::array< std::array< tk::real, 3 >, 5>& grad,
  std::array< std::array< tk::real, 6 >, 5>& hess ) const
// *****************************************************************************
//! Compute gradients of conserved quantities for an interior element
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] elem element id
//! \param[in] dBdx Basis gradients in element
//! \param[in] d2Bdx2 Second derivatives of basis functions in element
//! \param[in,out] grad Density, momentum, energy gradients in element
//! \param[in,out] hess Hessian of conserved quantities in element 
//! \details Hessian matrix contains (d2u_i/dx_i^2, du_i/dx_jdx_i, etc. terms)
// *****************************************************************************
{
  using inciter::multispecies::momentumDofIdx;
  using inciter::multispecies::densityDofIdx;
  using inciter::multispecies::energyDofIdx;
  
  // Gradient 

  // d(\rho)/dx
  for (std::size_t k=0; k<m_nspec; ++k) {
    for (std::size_t j=0; j<3; ++j) {
      grad[0][j] +=
        dBdx[j][1] * U(elem, densityDofIdx(m_nspec,k,m_rdof,1))
      + dBdx[j][2] * U(elem, densityDofIdx(m_nspec,k,m_rdof,2))
      + dBdx[j][3] * U(elem, densityDofIdx(m_nspec,k,m_rdof,3));
    }
  }

  // d(\rho u)/dx
  for (std::size_t i=1; i<4; ++i) {
    for (std::size_t j=0; j<3; ++j) {
      grad[i][j] = //drudx[i][j]
        dBdx[j][1] * U(elem, momentumDofIdx(m_nspec,i,m_rdof,1))
      + dBdx[j][2] * U(elem, momentumDofIdx(m_nspec,i,m_rdof,2))
      + dBdx[j][3] * U(elem, momentumDofIdx(m_nspec,i,m_rdof,3));
    }
  }

  // dE/dx_j
  // double check that this indexing is correct. I think this needs to change to calc energy from temp
  // because quinoa stores temp not energy?
  for (std::size_t j=0; j<3; ++j) {
    // dEdx[j]
    grad[4][j] = dBdx[j][1] * U(elem, energyDofIdx(m_nspec,0,m_rdof,1))
               + dBdx[j][2] * U(elem, energyDofIdx(m_nspec,0,m_rdof,2))
               + dBdx[j][3] * U(elem, energyDofIdx(m_nspec,0,m_rdof,3));
  }

  // Hessian 

  // d2(rho)dx_jdx_k
  for (std::size_t k=0; k<m_nspec; ++k) {
    for (std::size_t j=0; j<6; ++j) {
      hess[0][j] +=
        d2Bdx2[j][1] * U(elem, densityDofIdx(m_nspec,k,m_rdof,1))
      + d2Bdx2[j][2] * U(elem, densityDofIdx(m_nspec,k,m_rdof,2))
      + d2Bdx2[j][3] * U(elem, densityDofIdx(m_nspec,k,m_rdof,3));
    }
  }

  // d2(rho*u_i)/dx_jdx_k
  for (std::size_t i=1; i<4; ++i) {
    for (std::size_t j=0; j<6; ++j) {
      hess[i][j] =
        d2Bdx2[j][1] * U(elem, momentumDofIdx(m_nspec,i,m_rdof,1))
      + d2Bdx2[j][2] * U(elem, momentumDofIdx(m_nspec,i,m_rdof,2))
      + d2Bdx2[j][3] * U(elem, momentumDofIdx(m_nspec,i,m_rdof,3));
    }
  }

  // d2E/dx_jdx_k
  for (std::size_t j=0; j<6; ++j) {
    hess[4][j] = d2Bdx2[j][1] * U(elem, energyDofIdx(m_nspec,0,m_rdof,1))
               + d2Bdx2[j][2] * U(elem, energyDofIdx(m_nspec,0,m_rdof,2))
               + d2Bdx2[j][3] * U(elem, energyDofIdx(m_nspec,0,m_rdof,3));
  }

}


std::vector<std::array<tk::real, 3>>
MultiSpeciesViscousTermsDGP1::interiorFlux(
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t ncomp,
  const std::array< std::vector< tk::real >, 2 >& state,
  const std::array< tk::real, 3 >& fn,
  const tk::real he,
  const std::array< std::array< std::array< tk::real, 3 >, 5 >, 2 >& grad,
  const std::array< std::array< std::array< tk::real, 6 >, 5 >, 2 >& hess,
  std::vector< tk::real >& fl ) const
// *****************************************************************************
//! \brief Compute the multispecies viscous flux at an interior face
//! \param[in] mat_blk Material EOS block
//! \param[in] ncomp Number of components in this system
//! \param[in] state Left and right high-order face states
//! \param[in] cellAvgState Left and right cell-average states
//! \param[in] fn Face normal
//! \param[in] he Average diameter of right and left cells
//! \param[in] centroid Left and right element centroids
//! \param[in] grad Gradients of conserved quantities in left and right elements
//! \param[in] hess Hessians of the conserved quantities in left and right elements
//! \param[in,out] fl Numerical viscous flux
//! \details The adjoint property is leveraged to circumvent the need for computing
//!   the highly nonlinear diffusion matrix inside the numerical flux. Instead, 
//!   new direction vectors are computed which contain information from the diffusion
//!   matrices. Averages of the face states are used to compute these new direction vectors.
//!   A standard numerical flux is then computed only in terms of conserved quantities.
//!   Ref: Danis, M. E., and Yan, M. (2022). A new direct discontinuous Galerkin method with 
//!   interface correction for two-dimensional compressible Navier Stokes equations. Journal
//!   of Computational Physics, 452.
// *****************************************************************************
{
  // Pseudocode for DDG with interface correction
  using inciter::multispecies::densityIdx;
  using inciter::multispecies::temperatureIdx;
  using inciter::multispecies::momentumIdx;
  using inciter::multispecies::energyIdx;

  // 1. Compute diffusion matrices for average of left and right states
  Assert ( fl.size() == ncomp, "incorrect viscous flux vector size" );

  // There are ncomp x ncomp diffusion matrices A^(lxm) for l,m = 1,...,ncomp
  // For 3D compressible Navier-Stokes, ncomp = 5, so there are 25 diffusion matrices
  // Each diffusion matrix is 3x3
  // Initialize container for these matrices
  // Initialize direction vectors for each diffusion matrix, which will be used to compute the numerical flux
  // constexpr std::size_t ncomp = 5;

  auto neq = ncomp;
  using Vector3 = std::array<tk::real, 3>;
  using Matrix3 = std::array<std::array<tk::real, 3>, 3>;
  std::vector<std::vector<Matrix3>> A(
  ncomp, std::vector<Matrix3>(ncomp)
  );

  std::vector<Vector3> dir(ncomp * ncomp, Vector3{});
  // All matrices A^(1m) are zero matrices
  
  // Compute averages of left and right conserved quantities
  // [0] index is left, [1] index is right
  const auto& ul = state[0];
  const auto& ur = state[1];

  // establish mixture state
  inciter::Mixture mix_l(m_nspec, ul, mat_blk);
  inciter::Mixture mix_r(m_nspec, ur, mat_blk);

  // Compute fluid properties at left and right faces
  auto mix_temp_l = ul[ncomp + temperatureIdx(m_nspec,0)];
  auto mix_temp_r = ur[ncomp + temperatureIdx(m_nspec,0)];
  auto rho_l = mix_l.get_mix_density();
  auto rho_r = mix_r.get_mix_density();
  auto R_l = mat_blk[0].compute<inciter::EOS::gas_constant>();
  auto R_r = mat_blk[0].compute<inciter::EOS::gas_constant>();
  auto R = 0.5 * (R_l + R_r);
  auto cp_l = mix_l.Cp(mix_temp_l, mat_blk);
  auto cp_r = mix_r.Cp(mix_temp_r, mat_blk);
  auto gamma_l = cp_l / (cp_l - R);
  auto gamma_r = cp_r / (cp_r - R);
  auto mu_l = mix_l.viscCoeff(mix_temp_l, mat_blk);
  auto mu_r = mix_r.viscCoeff(mix_temp_r, mat_blk);

  // Compute averages (means) of left and right conserved quantities
  auto e_m = 0.5 * (ul[energyIdx(m_nspec,0)]/rho_l + ur[energyIdx(m_nspec,0)]/rho_r);
  auto u_m = 0.5 * (ul[momentumIdx(m_nspec,0)]/rho_l + ur[momentumIdx(m_nspec,0)]/rho_r);
  auto v_m = 0.5 * (ul[momentumIdx(m_nspec,1)]/rho_l + ur[momentumIdx(m_nspec,1)]/rho_r);
  auto w_m = 0.5 * (ul[momentumIdx(m_nspec,2)]/rho_l + ur[momentumIdx(m_nspec,2)]/rho_r);
  auto mu_m = 0.5 * (mu_l + mu_r);
  auto rho_m = 0.5 * (rho_l + rho_r);
  auto nu_m = mu_m / rho_m;
  auto gamma_m = 0.5 * (gamma_l + gamma_r);
  auto Pr = 0.71; // TODO: make Prandtl number user-configurable

  // Build diffusion matrices A^(lm) for l,m = 1,...,ncomp

  // x-momentum (mathematically, l=2 but l = 1 in zero-based indexing)
  A[1][0] = { { 
    { -nu_m*(4.0/3.0)*u_m, -nu_m*(2.0/3.0)*v_m, -nu_m*(2.0/3.0)*w_m }, 
    { -nu_m*v_m,           -nu_m*u_m,                    0.0        }, 
    { -nu_m*w_m,               0.0,                  -nu_m*u_m      } } };

  A[1][1] = { { 
    { nu_m*(4.0/3.0), 0.0,  0.0 }, 
    { 0.0,            nu_m, 0.0 }, 
    { 0.0,            0.0,  nu_m} } };

  A[1][2] = { { 
    { 0.0,  -nu_m*(2.0/3.0), 0.0 }, 
    { nu_m,        0.0,      0.0 }, 
    { 0.0,         0.0,      0.0 } } };

  A[1][3] = { { 
    { 0.0,  0.0, -nu_m*(2.0/3.0) }, 
    { 0.0,  0.0,        0.0      }, 
    { nu_m, 0.0,        0.0      } } };

  A[1][4] = {};
  
  // y-momentum (mathematically, l=3, but l=2 in zero-based indexing)
  A[2][0] = { { 
    { -nu_m*v_m,             -nu_m*u_m,                0.0        }, 
    { nu_m*(2.0/3.0)*u_m, -nu_m*(4.0/3.0)*v_m, nu_m*(2.0/3.0)*w_m }, 
    { 0.0,                   -nu_m*w_m,            -nu_m*v_m      } } };

  A[2][1] = { { 
    { 0.0,             nu_m,  0.0 }, 
    { -nu_m*(2.0/3.0), 0.0,   0.0 }, 
    { 0.0,             0.0,  nu_m } } };

  A[2][2] = { { 
    { nu_m,      0.0,      0.0 }, 
    { 0.0, nu_m*(4.0/3.0), 0.0 }, 
    { 0.0,       0.0,     nu_m } } };

  A[2][3] = { { 
    { 0.0, 0.0,        0.0      }, 
    { 0.0, 0.0, -nu_m*(2.0/3.0) }, 
    { 0.0, nu_m,       0.0      } } };

  A[2][4] = {};

  // z-momentum (mathematically, l=4, but l=3 in zero-based indexing)
  A[3][0] = { { 
    { -nu_m*w_m,               0.0,                  -nu_m*u_m }, 
    { 0.0,                 -nu_m*w_m,                -nu_m*v_m }, 
    { nu_m*(2.0/3.0)*u_m, nu_m*(2.0/3.0)*v_m, -nu_m*(4.0/3.0)*w_m } } };

  A[3][1] = { { 
    { 0.0,             0.0, nu_m }, 
    { 0.0,             0.0, 0.0  }, 
    { -nu_m*(2.0/3.0), 0.0, 0.0  } } };

  A[3][2] = { { 
    { 0.0,        0.0,      0.0 }, 
    { 0.0,        0.0,     nu_m }, 
    { 0.0, -nu_m*(2.0/3.0), 0.0 } } };

  A[3][3] = { { 
    { nu_m, 0.0,      0.0      }, 
    { 0.0, nu_m,      0.0      }, 
    { 0.0, 0.0, nu_m*(4.0/3.0) } } };

  A[3][4] = {};

  // energy (mathematically, l=5, but l=4 in zero-based indexing)
  A[4][0] = { { 
    { nu_m*((gamma_m/Pr - 4.0/3.0)*u_m*u_m + (gamma_m/Pr - 1)*(v_m*v_m + w_m*w_m) + e_m*gamma_m/(rho_m*Pr)),
      -nu_m*(1.0/3.0)*u_m*v_m,
      -nu_m*(1.0/3.0)*u_m*w_m }, 
    { -nu_m*(1.0/3.0)*u_m*v_m,
      nu_m*((gamma_m/Pr - 4.0/3.0)*v_m*v_m + (gamma_m/Pr - 1)*(u_m*u_m + w_m*w_m) + e_m*gamma_m/(rho_m*Pr)),
      -nu_m*(1.0/3.0)*v_m*w_m }, 
    { -nu_m*(1.0/3.0)*u_m*w_m, 
      -nu_m*(1.0/3.0)*v_m*w_m,
      nu_m*((gamma_m/Pr - 4.0/3.0)*w_m*w_m
      + (gamma_m/Pr - 1)*(u_m*u_m + v_m*v_m)
      + e_m*gamma_m/(rho_m*Pr)) } } };

  A[4][2] = { { 
    { -nu_m*(gamma_m/Pr - 1)*v_m,       -nu_m*(2.0/3.0)*u_m,                    0.0           }, 
    { nu_m*u_m,                    nu_m*(4.0/3.0 - gamma_m/Pr)*v_m,          nu_m*w_m         }, 
    { 0.0,                              -nu_m*(2.0/3.0)*w_m,       -nu_m*(gamma_m/Pr - 1)*v_m } } };

  A[4][3] = { { 
    { -nu_m*(gamma_m/Pr - 1)*w_m,                0.0,                -nu_m*(2.0/3.0)*u_m      }, 
    { 0.0,                           -nu_m*(gamma_m/Pr - 1)*w_m,     -nu_m*(2.0/3.0)*v_m      }, 
    { nu_m*u_m,                               nu_m*v_m,       nu_m*(4.0/3.0 - gamma_m/Pr)*w_m } } };

  A[4][4] = { { 
    { nu_m*gamma_m/Pr,      0.0,              0.0      }, 
    { 0.0,             nu_m*gamma_m/Pr,       0.0      }, 
    { 0.0,                  0.0,       nu_m*gamma_m/Pr } } };  

  // 2. Compute new direction vectors using diffusion matrices and face normal
  // Take the transpose of the diffusion matrices to get A^T, then dot with face normal to get new direction vectors; /xi = A^T * n
  //auto lmIndex = [ncomp](std::size_t l, std::size_t m) { return l*ncomp + m; };
  for (std::size_t l = 0; l < ncomp; ++l) {
    for (std::size_t m = 0; m < ncomp; ++m) {
      for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
          const auto k = l*ncomp + m;
          dir[k][i] += A[l][m][j][i] * fn[j];
        }
      }
    }
  }
  
  // 3. Compute numerical flux from left and right states
  // double check: am I calculating averages of state gradients or of velocity gradients?
  std::vector<Vector3> numgradQ(ncomp, Vector3{});
  constexpr std::array<std::array<std::size_t, 3>, 3> hidx{{
  {{ 0, 3, 4 }},  // x row: xx, xy, xz
  {{ 3, 1, 5 }},  // y row: xy, yy, yz
  {{ 4, 5, 2 }}   // z row: xz, yz, zz
  }};
  auto beta0 = 4.0; // beta0 = (k+1)^2 for all k = polynomial degree. So k = 1 => beta0 = 4 for DGP1
  auto beta1 = 0.25; // beta1 = 1/(2k*(k+1)) for all k = polynomial degree. So k = 1 => beta1 = 0.25

  for (std::size_t m = 0; m<ncomp; ++m) {
    for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j=0; j<3; ++j) {
      const auto hij = hidx[i][j];
      numgradQ[m][i] += beta1*he*( hess[1][m][hij] - hess[0][m][hij])*fn[j];
    }
      numgradQ[m][i] += (beta0/he)*(state[1][m] - state[0][m])*fn[i] + 0.5*(grad[0][m][i] + grad[1][m][i]);
  }
}

  // 4. Dot product of numerical flux with new direction vectors to compute flux

  for (std::size_t l=1; l<ncomp; ++l) {
    fl[l] = 0.0;
    for (std::size_t m=0; m<ncomp; ++m) {
      const auto k = l*ncomp + m;
      for (std::size_t i=0; i<3; ++i)
        fl[l] += numgradQ[m][i] * dir[k][i];
    }
  }

  // 6. Return direction vectors
  return dir;
  
}

void
MultiSpeciesViscousTermsDGP1::interfaceCorrection(
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t ncomp,
  const std::array< std::vector< tk::real >, 2 >& state,
  const std::vector<std::array<tk::real, 3>>& dir,
  const std::array< std::array< std::array< tk::real, 3 >, 5>, 2 >& grad,
  const std::array< std::array< std::array< tk::real, 6 >, 5>, 2 >& hess,
  std::vector< tk::real >& ic ) const
// *****************************************************************************
//! \brief Compute the interface correction term for the multispecies viscous flux at an interior face
//! \param[in] mat_blk Material EOS block
//! \param[in] ncomp Number of components in this system
//! \param[in] state Left and right high-order face states
//! \param[in] dir direction vectors computed from interiorFlux
//! \param[in] grad Gradients of conserved quantities in left and right elements
//! \param[in] hess Hessians of the conserved quantities in left and right elements
//! \param[in,out] ic Interface correction term for the numerical viscous flux
// *****************************************************************************
{
  // 1. Calculate state jumps
  for (std::size_t l=1; l<ncomp; ++l) {
    for (std::size_t m=0; m<ncomp; ++m){
      ic[l] = 0.5 * (state[1][m] - state[0][m]) * dir[l][m];
    }
  }

}
void
MultiSpeciesViscousTermsDGP1::volumeFlux(
  const std::vector< inciter::EOS >& mat_blk,
  std::size_t ncomp,
  const std::vector<tk::real>& state,
  std::vector<std::array<tk::real, 3>>& visc_fl ) const
// *****************************************************************************
//! \brief Compute the multispecies viscous flux on a tet volume
//! \param[in] mat_blk Material EOS block
//! \param[in] ncomp Number of components in this system
//! \param[in] state vector of conserved quantities
//! \param[in,out] visc_fl viscous volume flux
// *****************************************************************************
{
 // no-op
}
} // tk::
