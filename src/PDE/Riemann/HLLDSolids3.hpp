// *****************************************************************************
/*!
  \file      src/PDE/Riemann/HLLDSolids.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2023 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     HLLD Riemann flux function for solids
  \details   This file implements the HLLD Riemann solver for solids.
             Ref. Barton, P. T. (2019). An interface-capturing Godunov method
             for the simulation of compressible solid-fluid problems. Journal of
             Computational Physics, 390, 25-50.
*/
// *****************************************************************************
#ifndef HLLDSolids_h
#define HLLDSolids_h

#include <vector>

#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EOS.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "MultiMat/MiscMultiMatFns.hpp"

namespace inciter {

//! HLLD approximate Riemann solver for solids
struct HLLDSolids {

  //! HLLD approximate Riemann solver flux function
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann flux solution according to HLLD, appended by
  //!   Riemann velocities and volume-fractions.
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::vector< EOS >& mat_blk,
        const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::vector< std::array< tk::real, 3 > >& )
  {
    auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
    const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

    auto nsld = numSolids(nmat, solidx);
    auto ncomp = u[0].size()-(3+nmat+nsld*6);
    std::vector< tk::real > flx(ncomp,0);
    std::array< std::array< tk::real, 3 >, 3 > auxflxg;
    // Outer states
    tk::real arho_l(0.0), arho_r(0.0);
    std::array< tk::real, 3 > vn_l, vn_r;
    std::vector< tk::real > al_l(nmat, 0.0), al_r(nmat, 0.0),
      pm_l(nmat, 0.0), pm_r(nmat, 0.0), am_l(nmat, 0.0), am_r(nmat, 0.0);
    tk::real amatl(0.0), amatr(0.0), asmatl(0.0), asmatr(0.0);
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > g_l, g_r,
      asig_l, asig_r;
    std::vector< std::array< tk::real, 3 > > asign_l, asign_r;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gn_l, gn_r,
      asignn_l, asignn_r;
    std::array< tk::real, 3 > sign_l {{0, 0, 0}}, sign_r {{0, 0, 0}};
    // Tilde states
    tk::real arho_t_l(0.0), arho_t_r(0.0);
    std::vector< tk::real > al_t_l(nmat, 0.0), al_t_r(nmat, 0.0),
      arhoe_t_l(nmat, 0.0), arhoe_t_r(nmat, 0.0);
    std::array< tk::real, 3 > vn_t_l, vn_t_r;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gn_t_l, gn_t_r;
    // Asterisk states
    tk::real arho_a_l(0.0), arho_a_r(0.0);
    std::vector< tk::real > al_a_l(nmat, 0.0), al_a_r(nmat, 0.0),
                            arhoe_a_l(nmat, 0.0), arhoe_a_r(nmat, 0.0);
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gn_a_l, gn_a_r;
    std::array< tk::real, 3 > vn_a_l, vn_a_r;

    auto c_plus(0.0), c_minus(0.0), p_plus(0.0), p_minus(0.0);

    // Independently limited velocities for advection
    auto ul = u[0][ncomp+velocityIdx(nmat, 0)];
    auto vl = u[0][ncomp+velocityIdx(nmat, 1)];
    auto wl = u[0][ncomp+velocityIdx(nmat, 2)];
    auto ur = u[1][ncomp+velocityIdx(nmat, 0)];
    auto vr = u[1][ncomp+velocityIdx(nmat, 1)];
    auto wr = u[1][ncomp+velocityIdx(nmat, 2)];

    // Rotated velocities from advective velocities
    vn_l = tk::rotateVector({ul, vl, wl}, fn);
    vn_r = tk::rotateVector({ur, vr, wr}, fn);

    // Compute mass fractions and bulk density
    tk::real brho_l(0.0), brho_r(0.0);
    std::vector< tk::real > mass_frac_l(nmat, 0.0), mass_frac_r(nmat, 0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      // Left state
      arho_l = u[0][densityIdx(nmat, k)];
      brho_l += arho_l;
      // Right state
      arho_r = u[1][densityIdx(nmat, k)];
      brho_r += arho_r;
    }
    for (std::size_t k=0; k<nmat; ++k)
    {
      // Left state
      arho_l = u[0][densityIdx(nmat, k)];
      mass_frac_l[k] = arho_l/brho_l;
      // Right state
      arho_r = u[1][densityIdx(nmat, k)];
      mass_frac_r[k] = arho_r/brho_r;
    }

    for (std::size_t k=0; k<nmat; ++k)
    {
      // Left outer state
      // -----------------------------------------------------------------------
      al_l[k] = u[0][volfracIdx(nmat, k)];
      pm_l[k] = u[0][ncomp+pressureIdx(nmat, k)];
      arho_l = u[0][densityIdx(nmat, k)];

      // inv deformation gradient and Cauchy stress tensors
      g_l.push_back(getDeformGrad(nmat, k, u[0]));
      asig_l.push_back(getCauchyStress(nmat, k, ncomp, u[0]));
      for (std::size_t i=0; i<3; ++i) asig_l[k][i][i] -= pm_l[k];

      // normal stress (traction) vector
      asign_l.push_back(tk::matvec(asig_l[k], fn));
      for (std::size_t i=0; i<3; ++i)
        sign_l[i] += asign_l[k][i];

      // rotate stress vector
      asignn_l.push_back(tk::rotateTensor(asig_l[k], fn));

      // rotate deformation gradient tensor
      gn_l.push_back(tk::rotateTensor(g_l[k], fn));
      amatl = mat_blk[k].compute< EOS::soundspeed >(
        u[0][densityIdx(nmat, k)], pm_l[k], al_l[k], k, gn_l[k] );

      am_l[k] = amatl;
      
      // Right outer state
      // -----------------------------------------------------------------------
      al_r[k] = u[1][volfracIdx(nmat, k)];
      pm_r[k] = u[1][ncomp+pressureIdx(nmat, k)];
      arho_r = u[1][densityIdx(nmat, k)];

      // inv deformation gradient and Cauchy stress tensors
      g_r.push_back(getDeformGrad(nmat, k, u[1]));
      asig_r.push_back(getCauchyStress(nmat, k, ncomp, u[1]));
      for (std::size_t i=0; i<3; ++i) asig_r[k][i][i] -= pm_r[k];

      // normal stress (traction) vector
      asign_r.push_back(tk::matvec(asig_r[k], fn));
      for (std::size_t i=0; i<3; ++i)
        sign_r[i] += asign_r[k][i];

      // rotate stress vector
      asignn_r.push_back(tk::rotateTensor(asig_r[k], fn));
      
      // rotate deformation gradient tensor
      gn_r.push_back(tk::rotateTensor(g_r[k], fn));
      amatr = mat_blk[k].compute< EOS::soundspeed >(
        u[1][densityIdx(nmat, k)], pm_r[k], al_r[k], k, gn_r[k] );

      am_r[k] = amatr;

      // Signal velocities
      auto Sl_minus = std::min((vn_l[0]-amatl), (vn_r[0]-amatr));
      auto Sl_plus = std::max((vn_l[0]+amatl), (vn_r[0]+amatr));
      auto Si = (arho_l*vn_l[0]*(Sl_minus-vn_l[0]) - arho_r*vn_r[0]*(Sl_plus-vn_r[0])
        + asignn_l[k][0][0]-asignn_r[k][0][0])
        / (arho_l*(Sl_minus-vn_l[0]) - arho_r*(Sl_plus-vn_r[0]));

      // Left tilde state
      // -----------------------------------------------------------------------
      al_t_l[k] = al_l[k]*(Sl_minus-vn_l[0])/(Sl_minus-Si);
      arho_t_l = arho_l*(Sl_minus-vn_l[0])/(Sl_minus-Si);
      arhoe_t_l[k] = u[0][energyIdx(nmat, k)]
        + (Si-vn_l[0])*(Si - asignn_l[k][0][0]/(arho_l*(Sl_minus-vn_l[0])));
      vn_t_l[0] = Si;
      vn_t_l[1] = vn_l[1];
      vn_t_l[2] = vn_l[2];

      v_t_l = tk::unrotateVector(vn_t_l, fn);

      // inv deformation gradient and Cauchy stress tensors
      gn_t_l.push_back(gn_l[k]);
      gn_t_l[k][0][0] *= (Sl_minus-vn_l[0])/(Sl_minus-Si);
      gn_t_l[k][1][0] *= (Sl_minus-vn_l[0])/(Sl_minus-Si);
      gn_t_l[k][2][0] *= (Sl_minus-vn_l[0])/(Sl_minus-Si);

      // rotate g back to original frame of reference
      g_t_l.push_back(tk::unrotateTensor(gn_t_l[k], fn));

      // compute pressure
      pm_t_l[k] = mat_blk[k].compute< EOS::pressure >(arho_t_l[k], v_t_l[0], v_t_l[1], v_t_l[2],
                                                      arhoe_t_l[k], al_t_l[k], k, g_t_l[k]);

      // compute sigma from equation of state
      asig_t_l = mat_blk[k].computeTensor< EOS::CauchyStress >(
        0.0, 0.0, 0.0, 0.0, 0.0, al_t_l[k], k, g_t_l[k] );
      asign_t_l.push_back(tk::matvec(asig_t_l, fn));

      // Compute shear speed of sounds with tilde variables
      asmatl = mat_blk[k].compute< EOS::shearspeed >(
        arho_t_l, al_t_l[k], k );

      // Right tilde state
      // -----------------------------------------------------------------------
      al_t_r[k] = al_r[k]*(Sl_plus-vn_r[0])/(Sl_plus-Si);
      arho_t_r = arho_r*(Sl_plus-vn_r[0])/(Sl_plus-Si);
      arhoe_t_r[k] = u[1][energyIdx(nmat, k)]
        + (Si-vn_r[0])*(Si-asignn_r[k][0][0]/(arho_r*(Sl_plus-vn_r[0])));
      vn_t_r[0] = Si;
      vn_t_r[1] = vn_r[1];
      vn_t_r[2] = vn_r[2];

      v_t_r = tk::unrotateVector(vn_t_l, fn);

      // inv deformation gradient and Cauchy stress tensors
      gn_t_r.push_back(gn_r[k]);
      gn_t_r[k][0][0] *= (Sl_plus-vn_r[0])/(Sl_plus-Si);
      gn_t_r[k][1][0] *= (Sl_plus-vn_r[0])/(Sl_plus-Si);
      gn_t_r[k][2][0] *= (Sl_plus-vn_r[0])/(Sl_plus-Si);

      // rotate g back to original frame of reference
      g_t_r.push_back(tk::unrotateTensor(gn_t_r[k], fn));

      // compute pressure
      pm_t_r[k] = mat_blk[k].compute< EOS::pressure >(arho_t_r[k], v_t_r[0], v_t_r[1], v_t_r[2],
                                                      arhoe_t_r[k], al_t_r[k], k, g_t_r[k]);
      // compute sigma from equation of state
      asig_t_r = mat_blk[k].computeTensor< EOS::CauchyStress >(
        0.0, 0.0, 0.0, 0.0, 0.0, al_t_r[k], k, g_t_r[k] );
      asign_t_r.push_back(tk::matvec(asig_t_r, fn));

      // Compute shear speed of sounds with tilde variables
      asmatr = mat_blk[k].compute< EOS::shearspeed >(
        arho_t_r, al_t_r[k], k );

      // Extra signal velocities
      auto Ss_minus = vn_l[0] - asmatl;
      auto Ss_plus = vn_r[0] + asmatr;

      // left asterisk state
      // -----------------------------------------------------------------------
      al_a_l[k] = al_t_l[k];
      arho_a_l = arho_t_l;
      arhoe_a_l[k] = arhoe_t_l[k]
        + (asignn_l[k][0][1]*vn_l[1]+asignn_l[k][0][2]*vn_l[2])/(arho_t_l*(Ss_minus-Si));
      vn_a_l[0] = Si;
      vn_a_l[1] = vn_l[1] + asignn_l[k][0][1]/(arho_t_l*(Ss_minus-Si));
      vn_a_l[2] = vn_l[2] + asignn_l[k][0][2]/(arho_t_l*(Ss_minus-Si));

      v_a_l = tk::unrotateVector(vn_a_l, fn);
      
      // inv deformation gradient and Cauchy stress tensors
      gn_a_l.push_back(gn_t_l[k]);
      gn_a_l[k][0][0] -=
        (gn_l[k][0][1]*(vn_l[1]-vn_a_l[1])-gn_l[k][0][2]*(vn_l[2]-vn_a_l[2]))/(Ss_minus-Si);
      gn_a_l[k][1][0] -=
        (gn_l[k][1][1]*(vn_l[1]-vn_a_l[1])-gn_l[k][1][2]*(vn_l[2]-vn_a_l[2]))/(Ss_minus-Si);
      gn_a_l[k][2][0] -=
        (gn_l[k][2][1]*(vn_l[1]-vn_a_l[1])-gn_l[k][2][2]*(vn_l[2]-vn_a_l[2]))/(Ss_minus-Si);

      // rotate g back to original frame of reference
      g_a_l.push_back(tk::unrotateTensor(gn_a_l[k], fn));

      // compute pressure
      pm_a_l[k] = mat_blk[k].compute< EOS::pressure >(arho_a_l[k], v_a_l[0], v_a_l[1], v_a_l[2],
                                                      arhoe_a_l[k], al_a_l[k], k, g_a_l[k]);

      // compute sigma from equation of state
      asig_a_l = mat_blk[k].computeTensor< EOS::CauchyStress >(
        0.0, 0.0, 0.0, 0.0, 0.0, al_a_l[k], k, g_a_l[k] );
      asign_a_l.push_back(tk::matvec(asig_a_l, fn));

      // right asterisk state
      // -----------------------------------------------------------------------
      al_a_r[k] = al_t_r[k];
      arho_a_r = arho_t_r;
      arhoe_a_r[k] = arhoe_t_r[k]
        + (asignn_r[k][0][1]*vn_r[1]+asignn_r[k][0][2]*vn_r[2])/(arho_t_r*(Ss_plus-Si));
      vn_a_r[0] = Si;
      vn_a_r[1] = vn_r[1] + asignn_r[k][0][1]/(arho_t_r*(Ss_plus-Si));
      vn_a_r[2] = vn_r[2] + asignn_r[k][0][2]/(arho_t_r*(Ss_plus-Si));

      v_a_r = tk::unrotateVector(vn_a_r, fn);

      // inv deformation gradient and Cauchy stress tensors
      gn_a_r.push_back(gn_t_r[k]);
      gn_a_r[k][0][0] -=
        (gn_r[k][0][1]*(vn_r[1]-vn_a_r[1])-gn_r[k][0][2]*(vn_r[2]-vn_a_r[2]))/(Ss_plus-Si);
      gn_a_r[k][1][0] -=
        (gn_r[k][1][1]*(vn_r[1]-vn_a_r[1])-gn_r[k][1][2]*(vn_r[2]-vn_a_r[2]))/(Ss_plus-Si);
      gn_a_r[k][2][0] -=
        (gn_r[k][2][1]*(vn_r[1]-vn_a_r[1])-gn_r[k][2][2]*(vn_r[2]-vn_a_r[2]))/(Ss_plus-Si);

      
      // rotate g back to original frame of reference
      g_a_r.push_back(tk::unrotateTensor(gn_a_r[k], fn));

      // compute pressure
      pm_a_r[k] = mat_blk[k].compute< EOS::pressure >(arho_a_r[k], v_a_r[0], v_a_r[1], v_a_r[2],
                                                      arhoe_a_r[k], al_a_r[k], k, g_a_r[k]);
      // compute sigma from equation of state
      asig_a_r = mat_blk[k].computeTensor< EOS::CauchyStress >(
        0.0, 0.0, 0.0, 0.0, 0.0, al_a_r[k], k, g_a_r[k] );
      asign_a_r.push_back(tk::matvec(asig_a_r, fn));
      
      // Numerical fluxes
      if (Sl_minus >= 0.0) {
        // Left fluxes
        flx[volfracIdx(nmat, k)] = vn_l[0] * al_l[k];
        flx[densityIdx(nmat, k)] = vn_l[0] * arho_l;
        flx[energyIdx(nmat, k)] = vn_l[0] * u[0][energyIdx(nmat, k)];
        for (std::size_t i=0; i<3; ++i) {
          flx[energyIdx(nmat, k)] -= u[0][ncomp+velocityIdx(nmat,i)] * asign_l[k][i];
        }
        // inv deformation gradient tensor
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              auxflxg[i][j] = (
                g_l[k][i][0] * ul +
                g_l[k][i][1] * vl +
                g_l[k][i][2] * wl ) * fn[j];
        }
        // Momentum flux contribution
        for (std::size_t idir=0; idir<3; ++idir)
        {
          flx[momentumIdx(nmat, idir)] += mass_frac_l[k] * (
            vn_l[0]*u[0][momentumIdx(nmat, idir)]
            - sign_l[idir]; );
        }
        
        // Quantities for non-conservative terms
        // Store Riemann-advected partial pressures
        for (std::size_t k=0; k<nmat; ++k) ????
          flx.push_back(pml[k]);
        // Store Riemann velocity
        flx.push_back( vn_l[0] );
        // Store Riemann asign_ij (3*nsld)
        for (std::size_t k=0; k<nmat; ++k) {
          if (solidx[k] > 0) {
            for (std::size_t i=0; i<3; ++i)
              flx.push_back(asign_l[k][i]);
          }
        }
      }
      // else if (Sl_minus <= 0.0 && Ss_minus > 0.0) {
      //   // Left fluxes
      //   flx[volfracIdx(nmat, k)] = vn_l[0] * al_l[k]
      //     + Sl_minus*(al_t_l[k] - al_l[k]);
      //   flx[densityIdx(nmat, k)] = vn_l[0] * arho_l
      //     + Sl_minus*(arho_t_l - arho_l);
      //   flx[energyIdx(nmat, k)] = vn_l[0] * u[0][energyIdx(nmat, k)]
      //     + Sl_minus*(arhoe_t_l[k] - u[0][energyIdx(nmat, k)]);
      //   for (std::size_t i=0; i<3; ++i) {
      //     flx[energyIdx(nmat, k)] -= vn_l[i] * asign_l[k][i];
      //   }
      //   // inv deformation gradient tensor
      //   if (solidx[k] > 0) {
      //     for (std::size_t i=0; i<3; ++i)
      //       for (std::size_t j=0; j<3; ++j)
      //         auxflxg[i][j] = (
      //           gn_l[k][i][0] * vn_l[0] +
      //           gn_l[k][i][1] * vn_l[1] +
      //           gn_l[k][i][2] * vn_l[2] ) * fn[j]
      //           + Sl_minus*(gn_t_l[k][i][j] - gn_l[k][i][j]);
      //   }
      //   // Momentum flux contribution
      //   for (std::size_t idir=0; idir<3; ++idir)
      //   {
      //     flx[momentumIdx(nmat, idir)] += mass_frac_l[k]*(
      //       vn_l[0]*u[0][momentumIdx(nmat, idir)] - sign_l[idir]
      //       + Sl_minus*(arho_t_l*vn_t_l[idir] - arho_l*vn_l[idir]));
      //   }
      // }
      // else if (Ss_minus <= 0.0 && Si > 0.0) {
      //   // Left fluxes
      //   flx[volfracIdx(nmat, k)] = vn_l[0] * al_l[k]
      //     + Sl_minus*(al_t_l[k] - al_l[k])
      //     + Ss_minus*(al_a_l[k] - al_t_l[k]);
      //   flx[densityIdx(nmat, k)] = vn_l[0] * arho_l
      //     + Sl_minus*(arho_t_l - arho_l) + Ss_minus*(arho_a_l - arho_t_l);
      //   flx[energyIdx(nmat, k)] = vn_l[0] * u[0][energyIdx(nmat, k)]
      //     + Sl_minus*(arhoe_t_l[k] - u[0][energyIdx(nmat, k)])
      //     + Ss_minus*(arhoe_a_l[k] - arhoe_t_l[k]);
      //   for (std::size_t i=0; i<3; ++i) {
      //     flx[energyIdx(nmat, k)] -= vn_l[i] * asign_l[k][i];
      //   }
      //   // inv deformation gradient tensor
      //   if (solidx[k] > 0) {
      //     for (std::size_t i=0; i<3; ++i)
      //       for (std::size_t j=0; j<3; ++j)
      //         auxflxg[i][j] = (
      //           gn_l[k][i][0] * vn_l[0] +
      //           gn_l[k][i][1] * vn_l[1] +
      //           gn_l[k][i][2] * vn_l[2] ) * fn[j]
      //           + Sl_minus*(gn_t_l[k][i][j] - gn_l[k][i][j])
      //           + Ss_minus*(gn_a_l[k][i][j] - gn_t_l[k][i][j]);
      //   }
      //   // Momentum flux contribution
      //   for (std::size_t idir=0; idir<3; ++idir)
      //   {
      //     flx[momentumIdx(nmat, idir)] += mass_frac_l[k]*(
      //       vn_l[0]*u[0][momentumIdx(nmat, idir)] - sign_l[idir]
      //       + Sl_minus*(arho_t_l*vn_t_l[idir] - arho_l*vn_l[idir])
      //       + Ss_minus*(arho_a_l*vn_a_l[idir] - arho_t_l*vn_t_l[idir]));
      //   }
      // }
      // else if (Si <= 0.0 && Ss_plus > 0.0) {
      //   // Right fluxes
      //   flx[volfracIdx(nmat, k)] = vn_r[0] * al_r[k]
      //     + Sl_plus*(al_t_r[k] - al_r[k])
      //     + Ss_plus*(al_a_r[k] - al_t_r[k]);
      //   flx[densityIdx(nmat, k)] = vn_r[0] * arho_r
      //     + Sl_plus*(arho_t_r - arho_r) + Ss_plus*(arho_a_r - arho_t_r);
      //   flx[energyIdx(nmat, k)] = vn_r[0] * u[1][energyIdx(nmat, k)]
      //     + Sl_plus*(arhoe_t_r[k] - u[1][energyIdx(nmat, k)])
      //     + Ss_plus*(arhoe_a_r[k] - arhoe_t_r[k]);
      //   for (std::size_t i=0; i<3; ++i) {
      //     flx[energyIdx(nmat, k)] -= vn_r[i] * asign_r[k][i];
      //   }
      //   // inv deformation gradient tensor
      //   if (solidx[k] > 0) {
      //     for (std::size_t i=0; i<3; ++i)
      //       for (std::size_t j=0; j<3; ++j)
      //         auxflxg[i][j] = (
      //           gn_r[k][i][0] * vn_r[0] +
      //           gn_r[k][i][1] * vn_r[1] +
      //           gn_r[k][i][2] * vn_r[2] ) * fn[j]
      //           + Sl_plus*(gn_t_r[k][i][j] - gn_r[k][i][j])
      //           + Ss_plus*(gn_a_r[k][i][j] - gn_t_r[k][i][j]);
      //   }
      //   // Momentum flux contribution
      //   for (std::size_t idir=0; idir<3; ++idir)
      //   {
      //     flx[momentumIdx(nmat, idir)] += mass_frac_r[k]*(
      //       vn_r[0]*u[1][momentumIdx(nmat, idir)] - sign_r[idir]
      //       + Sl_plus*(arho_t_r*vn_t_r[idir] - arho_r*vn_r[idir])
      //       + Ss_plus*(arho_a_r*vn_a_r[idir] - arho_t_r*vn_t_r[idir]));
      //   }
      // }
      // else if (Ss_plus <= 0.0 && Sl_minus > 0.0) {
      //   // Right fluxes
      //   flx[volfracIdx(nmat, k)] = vn_r[0] * al_r[k]
      //     + Sl_plus*(al_t_r[k] - al_r[k]);
      //   flx[densityIdx(nmat, k)] = vn_r[0] * arho_r
      //     + Sl_plus*(arho_t_r - arho_r);
      //   flx[energyIdx(nmat, k)] = vn_r[0] * u[1][energyIdx(nmat, k)]
      //     + Sl_plus*(arhoe_t_r[k] - u[1][energyIdx(nmat, k)]);
      //   for (std::size_t i=0; i<3; ++i) {
      //     flx[energyIdx(nmat, k)] -= vn_r[i] * asign_r[k][i];
      //   }
      //   // inv deformation gradient tensor
      //   if (solidx[k] > 0) {
      //     for (std::size_t i=0; i<3; ++i)
      //       for (std::size_t j=0; j<3; ++j)
      //         auxflxg[i][j] = (
      //           gn_r[k][i][0] * vn_r[0] +
      //           gn_r[k][i][1] * vn_r[1] +
      //           gn_r[k][i][2] * vn_r[2] ) * fn[j]
      //           + Sl_plus*(gn_t_r[k][i][j] - gn_r[k][i][j]);
      //   }
      //   // Momentum flux contribution
      //   for (std::size_t idir=0; idir<3; ++idir)
      //   {
      //     flx[momentumIdx(nmat, idir)] += mass_frac_r[k]*(
      //       vn_r[0]*u[1][momentumIdx(nmat, idir)] - sign_r[idir]
      //       + Sl_plus*(arho_t_r*vn_t_r[idir] - arho_r*vn_r[idir]));
      //   }
      // }
      else if (Sl_plus <= 0.0) {
        // Right fluxes
        flx[volfracIdx(nmat, k)] = vn_r[0] * al_r[k];
        flx[densityIdx(nmat, k)] = vn_r[0] * u[1][densityIdx(nmat, k)];
        flx[energyIdx(nmat, k)] = vn_r[0] * u[1][energyIdx(nmat, k)];
        for (std::size_t i=0; i<3; ++i) {
          flx[energyIdx(nmat, k)] -= u[1][ncomp+velocityIdx(nmat,i)] * asign_r[k][i];
        }
        // inv deformation gradient tensor
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              /*auxflxg[i][j] = (
                gn_r[k][i][0] * vn_r[0] +
                gn_r[k][i][1] * vn_r[1] +
                gn_r[k][i][2] * vn_r[2] ) * fn[j];*/
              auxflxg[i][j] = (
                g_r[k][i][0] * ur +
                g_r[k][i][1] * vr +
                g_r[k][i][2] * wr ) * fn[j];
        }
        // Momentum flux contribution
        for (std::size_t idir=0; idir<3; ++idir)
        {
          flx[momentumIdx(nmat, idir)] += mass_frac_r[k] * (
            vn_r[0]*vn_r[idir]*arho_r
            - sign_r[idir] );
        }
        // Other, temporary stuff
        c_minus = vn_r[0];
        p_minus = 1.0;
      }
      else // Intermediate state in order to recover HLL
      {
        flx[volfracIdx(nmat, k)] = Sl_plus * vn_l[0] * al_l[k]
          - Sl_minus * vn_r[0] * al_r[k]
          + Sl_minus * Sl_plus * (al_r[k] - al_l[k]) / (Sl_plus-Sl_minus);
        flx[densityIdx(nmat, k)] = Sl_plus * vn_l[0] * arho_l
          - Sl_minus * vn_r[0] * arho_r
          + Sl_minus * Sl_plus * (arho_r - arho_l) / (Sl_plus-Sl_minus);
        flx[energyIdx(nmat, k)] = Sl_plus * vn_l[0] * u[0][energyIdx(nmat, k)]
          - Sl_minus * vn_r[0] * u[1][energyIdx(nmat, k)]
          + Sl_minus * Sl_plus * (u[1][energyIdx(nmat, k)] - u[0][energyIdx(nmat, k)])
          / (Sl_plus-Sl_minus);
        for (std::size_t i=0; i<3; ++i) {
          flx[energyIdx(nmat, k)] -= Sl_plus * u[0][ncomp+velocityIdx(nmat,i)] * asign_l[k][i]
            - Sl_minus * u[1][ncomp+velocityIdx(nmat,i)] * asign_r[k][i];
        }
        // inv deformation gradient tensor
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              /*auxflxg[i][j] = Sl_plus * (
                gn_l[k][i][0] * vn_l[0] +
                gn_l[k][i][1] * vn_l[1] +
                gn_l[k][i][2] * vn_l[2] ) * fn[j]
                - Sl_minus * (
                gn_r[k][i][0] * vn_r[0] +
                gn_r[k][i][1] * vn_r[1] +
                gn_r[k][i][2] * vn_r[2] ) * fn[j]
                + Sl_minus * Sl_plus
                * (gn_r[k][i][j] - gn_l[k][i][j]) / (Sl_plus-Sl_minus);*/
              auxflxg[i][j] = Sl_plus * (
                g_l[k][i][0] * ul +
                g_l[k][i][1] * vl +
                g_l[k][i][2] * wl ) * fn[j]
                - Sl_minus * (
                g_r[k][i][0] * ur +
                g_r[k][i][1] * vr +
                g_r[k][i][2] * wr ) * fn[j]
                + Sl_minus * Sl_plus
                * (g_r[k][i][j] - g_l[k][i][j]) / (Sl_plus-Sl_minus);
        }
        // Momentum flux contribution
        for (std::size_t idir=0; idir<3; ++idir)
        {
          flx[momentumIdx(nmat, idir)] += mass_frac_l[k] * Sl_plus
            * ( vn_l[0]*vn_l[idir]*arho_l - sign_l[idir])
            - mass_frac_r[k] * Sl_minus
            * ( vn_r[0]*vn_r[idir]*arho_r - sign_r[idir])
            + Sl_minus * Sl_plus
            * (mass_frac_r[k]*arho_r*vn_r[idir]
               -mass_frac_l[k]*arho_l*vn_l[idir]) / (Sl_plus-Sl_minus);
        }
        // Other, temporary stuff
        c_plus  = (Sl_plus*vn_l[0] - Sl_plus*Sl_minus) / (Sl_plus-Sl_minus);
        c_minus = (Sl_plus*Sl_minus - Sl_minus*vn_r[0]) / (Sl_plus-Sl_minus);
        p_plus = Sl_plus / (Sl_plus-Sl_minus);
        p_minus = -Sl_minus / (Sl_plus-Sl_minus);
      }

      // Unrotate tensor fluxes to align with cartesian coordinates
      if (solidx[k] > 0) {
        //auxflxg = tk::unrotateTensor(auxflxg, fn);
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            flx[deformIdx(nmat,solidx[k],i,j)] = auxflxg[i][j];
      }
      
    }

    // Unrotate velocity flux to align with cartesian coordinates
    std::array< tk::real, 3 > auxflxv;
    auxflxv[0] = flx[momentumIdx(nmat, 0)];
    auxflxv[1] = flx[momentumIdx(nmat, 1)];
    auxflxv[2] = flx[momentumIdx(nmat, 2)];
    auxflxv = tk::unrotateVector(auxflxv, fn);
    flx[momentumIdx(nmat, 0)] = auxflxv[0];
    flx[momentumIdx(nmat, 1)] = auxflxv[1];
    flx[momentumIdx(nmat, 2)] = auxflxv[2];

    // Store Riemann-advected partial pressures
    for (std::size_t k=0; k<nmat; ++k)
      flx.push_back(p_plus*pm_l[k] + p_minus*pm_r[k]);

    // Store Riemann velocity
    flx.push_back( c_plus+c_minus );

    // Store Riemann asign_ij (3*nsld)
    for (std::size_t k=0; k<nmat; ++k) {
      if (solidx[k] > 0) {
        for (std::size_t i=0; i<3; ++i)
          flx.push_back(p_plus*asign_l[k][i] + p_minus*asign_r[k][i]);
      }
    }


    Assert( flx.size() == (ncomp+nmat+1+3*nsld), "Size of "
            "multi-material flux vector incorrect" );

    return flx;
  }

  ////! Flux type accessor
  ////! \return Flux type
  //static ctr::FluxType type() noexcept {
  //  return ctr::FluxType::HLLDSolids; }
};

} // inciter::

#endif // HLLDSolids_h
