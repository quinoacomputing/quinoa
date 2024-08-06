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
    std::vector< tk::real > flx(ncomp, 0), fl(ncomp, 0), fr(ncomp, 0),
      ftl(ncomp, 0), ftr(ncomp, 0), fal(ncomp, 0), far(ncomp, 0);

    // Primitive quantities
    tk::real rhol(0.0), rhor(0.0);
    tk::real amatl(0.0), amatr(0.0), ac_l(0.0), ac_r(0.0),
      acs_l(0.0), acs_r(0.0);
    std::vector< tk::real > al_l(nmat, 0.0), al_r(nmat, 0.0),
                            pml(nmat, 0.0), pmr(nmat, 0.0);
    std::array< tk::real, 3 > vn_l, vn_r;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > g_l, g_r,
      asig_l, asig_r;
    std::vector< std::array< tk::real, 3 > > asign_l, asign_r;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gn_l, gn_r,
      asignn_l, asignn_r;
    std::array< std::array< tk::real, 3 >, 3 > signn_l, signn_r, sig_l, sig_r;
    std::array< tk::real, 3 > sign_l {{0, 0, 0}}, sign_r {{0, 0, 0}};
    // Tilde states
    tk::real rho_t_l(0.0), rho_t_r(0.0);
    std::vector< tk::real > al_t_l(nmat, 0.0), al_t_r(nmat, 0.0),
      arho_t_l(nmat, 0.0), arho_t_r(nmat, 0.0),
      arhoe_t_l(nmat, 0.0), arhoe_t_r(nmat, 0.0),
      pm_t_l(nmat, 0.0), pm_t_r(nmat, 0.0);
    std::array< tk::real, 3 > vn_t_l, vn_t_r, v_t_l, v_t_r;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gn_t_l, gn_t_r,
      g_t_l, g_t_r;
    std::array< std::array< tk::real, 3 >, 3 > asig_t_l, asig_t_r;
    std::vector< std::array< tk::real, 3 > > asign_t_l, asign_t_r;
    // Asterisk states
    tk::real rho_a_l(0.0), rho_a_r(0.0);
    std::vector< tk::real > al_a_l(nmat, 0.0), al_a_r(nmat, 0.0),
      arho_a_l(nmat, 0.0), arho_a_r(nmat, 0.0),
      arhoe_a_l(nmat, 0.0), arhoe_a_r(nmat, 0.0),
      pm_a_l(nmat, 0.0), pm_a_r(nmat, 0.0);
    std::array< tk::real, 3 > vn_a_l, vn_a_r, v_a_l, v_a_r;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gn_a_l, gn_a_r,
      g_a_l, g_a_r;
    std::array< std::array< tk::real, 3 >, 3 > asig_a_l, asig_a_r;
    std::vector< std::array< tk::real, 3 > > asign_a_l, asign_a_r;

    // Initialize bulk tensors
    for (std::size_t i=0; i<3; ++i)
      for (std::size_t j=0; j<3; ++j)
      {
        sig_l[i][j] = 0.0;
        sig_r[i][j] = 0.0;
        signn_l[i][j] = 0.0;
        signn_r[i][j] = 0.0;
      }
    
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

    // Outer states
    for (std::size_t k=0; k<nmat; ++k)
    {
      // Left state
      // -----------------------------------------------------------------------
      al_l[k] = u[0][volfracIdx(nmat, k)];
      pml[k] = u[0][ncomp+pressureIdx(nmat, k)];
      rhol += u[0][densityIdx(nmat, k)];

      // inv deformation gradient and Cauchy stress tensors
      g_l.push_back(getDeformGrad(nmat, k, u[0]));
      asig_l.push_back(getCauchyStress(nmat, k, ncomp, u[0]));
      for (std::size_t i=0; i<3; ++i) asig_l[k][i][i] -= pml[k];

      // normal stress (traction) vector
      asign_l.push_back(tk::matvec(asig_l[k], fn));
      for (std::size_t i=0; i<3; ++i)
        sign_l[i] += asign_l[k][i];

      // rotate stress vector
      asignn_l.push_back(tk::rotateTensor(asig_l[k], fn));
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          signn_l[i][j] += asignn_l[k][i][j];

      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          sig_l[i][j] += asig_l[k][i][j];

      // rotate deformation gradient tensor for speed of sound in normal dir
      gn_l.push_back(tk::rotateTensor(g_l[k], fn));
      amatl = mat_blk[k].compute< EOS::soundspeed >(
        u[0][densityIdx(nmat, k)], pml[k], al_l[k], k, gn_l[k] );

      // Right state
      // -----------------------------------------------------------------------
      al_r[k] = u[1][volfracIdx(nmat, k)];
      pmr[k] = u[1][ncomp+pressureIdx(nmat, k)];
      rhor += u[1][densityIdx(nmat, k)];

      // inv deformation gradient and Cauchy stress tensors
      g_r.push_back(getDeformGrad(nmat, k, u[1]));
      asig_r.push_back(getCauchyStress(nmat, k, ncomp, u[1]));
      for (std::size_t i=0; i<3; ++i) asig_r[k][i][i] -= pmr[k];

      // normal stress (traction) vector
      asign_r.push_back(tk::matvec(asig_r[k], fn));
      for (std::size_t i=0; i<3; ++i)
        sign_r[i] += asign_r[k][i];

      // rotate stress vector
      asignn_r.push_back(tk::rotateTensor(asig_r[k], fn));
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          signn_r[i][j] += asignn_r[k][i][j];

      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          sig_r[i][j] += asig_r[k][i][j];

      // rotate deformation gradient tensor for speed of sound in normal dir
      gn_r.push_back(tk::rotateTensor(g_r[k], fn));
      amatr = mat_blk[k].compute< EOS::soundspeed >(
        u[1][densityIdx(nmat, k)], pmr[k], al_r[k], k, gn_r[k] );

      // Mixture speed of sound
      // -----------------------------------------------------------------------
      ac_l += u[0][densityIdx(nmat, k)] * amatl * amatl;
      ac_r += u[1][densityIdx(nmat, k)] * amatr * amatr;
    }

    ac_l = std::sqrt(ac_l/rhol);
    ac_r = std::sqrt(ac_r/rhor);
    
    // Signal velocities
    auto Sl = std::min((vn_l[0]-ac_l), (vn_r[0]-ac_r));
    auto Sr = std::max((vn_l[0]+ac_l), (vn_r[0]+ac_r));
    auto Si = (rhol*vn_l[0]*(Sl-vn_l[0])
               -rhor*vn_r[0]*(Sr-vn_r[0])
               +signn_l[0][0]-signn_r[0][0])
      / (rhol*(Sl-vn_l[0]) - rhor*(Sr-vn_r[0]));

    // Left tilde velocity
    vn_t_l[0] = Si;
    vn_t_l[1] = vn_l[1];
    vn_t_l[2] = vn_l[2];
    // Right tilde velocity
    vn_t_r[0] = Si;
    vn_t_r[1] = vn_r[1];
    vn_t_r[2] = vn_r[2];
    // Rotate velocity back
    v_t_l = tk::unrotateVector(vn_t_l, fn);
    v_t_r = tk::unrotateVector(vn_t_r, fn);

    // Tilde states
    for (std::size_t k=0; k<nmat; ++k)
    {
      // Left tilde state
      // -----------------------------------------------------------------------
      al_t_l[k] = al_l[k]*(Sl-vn_l[0])/(Sl-Si);
      arho_t_l[k] = u[0][densityIdx(nmat, k)]*(Sl-vn_l[0])/(Sl-Si);
      rho_t_l += arho_t_l[k];

      // inv deformation gradient
      gn_t_l.push_back(gn_l[k]);
      gn_t_l[k][0][0] *= (Sl-vn_l[0])/(Sl-Si);
      gn_t_l[k][1][0] *= (Sl-vn_l[0])/(Sl-Si);
      gn_t_l[k][2][0] *= (Sl-vn_l[0])/(Sl-Si);

      // rotate g back to original frame of reference
      g_t_l.push_back(tk::unrotateTensor(gn_t_l[k], fn));

      // // compute pressure
      // pm_t_l[k] = mat_blk[k].compute< EOS::pressure >(arho_t_l[k], v_t_l[0], v_t_l[1], v_t_l[2],
      //                                                 arhoe_t_l[k], al_t_l[k], k, g_t_l[k]);

      // compute pressure
      pm_t_l[k] = pml[k] + u[0][densityIdx(nmat, k)]*(Sl-vn_l[0])*(Si-vn_l[0]);

      // // compute sigma from equation of state
      // asig_t_l = mat_blk[k].computeTensor< EOS::CauchyStress >(
      //   0.0, 0.0, 0.0, 0.0, 0.0, al_t_l[k], k, g_t_l[k] );
      // for (std::size_t i=0; i<3; ++i) asig_t_l[i][i] -= pm_t_l[k];

      // compute sigma
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
        {
          if (i==j) {
            asig_t_l[i][j] = asig_l[k][i][j]
              + u[0][densityIdx(nmat, k)] * (Sl-vn_l[0]) * (Si-vn_l[0]);
          }
          else
          {
            asig_t_l[i][j] = asig_l[k][i][j];
          }
        }
      asign_t_l.push_back(tk::matvec(asig_t_l, fn));

      amatl = mat_blk[k].compute< EOS::shearspeed >(
        arho_t_l[k], al_t_l[k], k );

      arhoe_t_l[k] = u[0][energyIdx(nmat, k)]
        + (Si-vn_l[0])*(Si-asignn_l[k][0][0]/(u[0][densityIdx(nmat, k)]*(Sl-vn_l[0])));
      // arhoe_t_l[k] = u[0][energyIdx(nmat, k)] * (Sl-vn_l[0])/(Sl-Si)
      //   + (Si-vn_l[0])*(u[0][densityIdx(nmat, k)]*Si - asignn_l[k][0][0]/(Sl-vn_l[0]));

      // Right tilde state
      // -----------------------------------------------------------------------
      al_t_r[k] = al_r[k]*(Sr-vn_r[0])/(Sr-Si);
      arho_t_r[k] = u[1][densityIdx(nmat, k)]*(Sr-vn_r[0])/(Sr-Si);
      rho_t_r += arho_t_r[k];

      // inv deformation gradient
      gn_t_r.push_back(gn_r[k]);
      gn_t_r[k][0][0] *= (Sr-vn_r[0])/(Sr-Si);
      gn_t_r[k][1][0] *= (Sr-vn_r[0])/(Sr-Si);
      gn_t_r[k][2][0] *= (Sr-vn_r[0])/(Sr-Si);

      // rotate g back to original frame of reference
      g_t_r.push_back(tk::unrotateTensor(gn_t_r[k], fn));

      // // compute pressure
      // pm_t_r[k] = mat_blk[k].compute< EOS::pressure >(arho_t_r[k], v_t_r[0], v_t_r[1], v_t_r[2],
      //                                                 arhoe_t_r[k], al_t_r[k], k, g_t_r[k]);

      // compute pressure
      pm_t_r[k] = pmr[k] + u[1][densityIdx(nmat, k)]*(Sr-vn_r[0])*(Si-vn_r[0]);
      
      // // compute sigma from equation of state
      // asig_t_r = mat_blk[k].computeTensor< EOS::CauchyStress >(
      //   0.0, 0.0, 0.0, 0.0, 0.0, al_t_r[k], k, g_t_r[k] );
      // for (std::size_t i=0; i<3; ++i) asig_t_r[i][i] -= pm_t_r[k];
      
      // compute sigma
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
        {
          if (i==j) {
            asig_t_r[i][j] = asig_r[k][i][j]
              + u[1][densityIdx(nmat, k)] * (Sr-vn_r[0]) * (Si-vn_r[0]);
          }
          else
          {
            asig_t_r[i][j] = asig_r[k][i][j];
          }
        }
      
      asign_t_r.push_back(tk::matvec(asig_t_r, fn));

      amatr = mat_blk[k].compute< EOS::shearspeed >(
        arho_t_r[k], al_t_r[k], k );

      arhoe_t_r[k] = u[1][energyIdx(nmat, k)]
       + (Si-vn_r[0])*(Si-asignn_r[k][0][0]/(u[1][densityIdx(nmat, k)]*(Sr-vn_r[0])));
      // arhoe_t_r[k] = u[1][energyIdx(nmat, k)] * (Sr-vn_r[0])/(Sr-Si)
      //   + (Si-vn_r[0])*(u[1][densityIdx(nmat, k)]*Si - asignn_r[k][0][0]/(Sr-vn_r[0]));

      // Mixture speed of sound
      // -----------------------------------------------------------------------
      acs_l += arho_t_l[k] * amatl * amatl;
      acs_r += arho_t_r[k] * amatr * amatr;
    }

    acs_l = std::sqrt(acs_l/rho_t_l);
    acs_r = std::sqrt(acs_r/rho_t_r);

    // Extra signal velocities
    auto Ssl = vn_t_l[0] - acs_l;
    auto Ssr = vn_t_r[0] + acs_r;

    // Left asterisk velocity
    vn_a_l[0] = vn_t_l[0];
    vn_a_l[1] = vn_t_l[1] + sig_l[0][1]/(rho_t_l*(Ssl-Si));
    vn_a_l[2] = vn_t_l[2] + sig_l[0][2]/(rho_t_l*(Ssl-Si));
    // Right asterisk velocity
    vn_a_r[0] = vn_t_r[0];
    vn_a_r[1] = vn_t_r[1] + sig_r[0][1]/(rho_t_r*(Ssr-Si));
    vn_a_r[2] = vn_t_r[2] + sig_r[0][2]/(rho_t_r*(Ssr-Si));
    // Rotate velocity back
    v_a_l = tk::unrotateVector(vn_a_l, fn);
    v_a_r = tk::unrotateVector(vn_a_r, fn);

    // Asterisk states
    for (std::size_t k=0; k<nmat; ++k)
    {
      // left asterisk state
      al_a_l[k] = al_t_l[k];
      arho_a_l[k] = arho_t_l[k];
      rho_a_l += arho_a_l[k];

      // inv deformation gradient
      gn_a_l.push_back(gn_t_l[k]);
      gn_a_l[k][0][0] -=
        (gn_l[k][0][1]*(vn_l[1]-vn_a_l[1])-gn_l[k][0][2]*(vn_l[2]-vn_a_l[2]))/(Ssl-Si);
      gn_a_l[k][1][0] -=
        (gn_l[k][1][1]*(vn_l[1]-vn_a_l[1])-gn_l[k][1][2]*(vn_l[2]-vn_a_l[2]))/(Ssl-Si);
      gn_a_l[k][2][0] -=
        (gn_l[k][2][1]*(vn_l[1]-vn_a_l[1])-gn_l[k][2][2]*(vn_l[2]-vn_a_l[2]))/(Ssl-Si);

      // rotate g back to original frame of reference
      g_a_l.push_back(tk::unrotateTensor(gn_a_l[k], fn));

      // compute pressure
      pm_a_l[k] = mat_blk[k].compute< EOS::pressure >(arho_a_l[k], v_a_l[0], v_a_l[1], v_a_l[2],
                                                      arhoe_a_l[k], al_a_l[k], k, g_a_l[k]);
      // compute sigma from equation of state
      asig_a_l = mat_blk[k].computeTensor< EOS::CauchyStress >(
        0.0, 0.0, 0.0, 0.0, 0.0, al_a_l[k], k, g_a_l[k] );
      for (std::size_t i=0; i<3; ++i) asig_a_l[i][i] -= pm_a_l[k];
      asign_a_l.push_back(tk::matvec(asig_a_l, fn));

      arhoe_a_l[k] = arhoe_t_l[k]
        + vn_a_l[1]*asig_a_l[1][0] - vn_l[1]*asig_l[k][1][0]
        + vn_a_l[2]*asig_a_l[2][0] - vn_l[2]*asig_l[k][2][0];

      // right asterisk state
      al_a_r[k] = al_t_r[k];
      arho_a_r[k] = arho_t_r[k];
      rho_a_r += arho_a_r[k];

      // inv deformation gradient
      gn_a_r.push_back(gn_t_r[k]);
      gn_a_r[k][0][0] -=
        (gn_r[k][0][1]*(vn_r[1]-vn_a_r[1])-gn_r[k][0][2]*(vn_r[2]-vn_a_r[2]))/(Ssl-Si);
      gn_a_r[k][1][0] -=
        (gn_r[k][1][1]*(vn_r[1]-vn_a_r[1])-gn_r[k][1][2]*(vn_r[2]-vn_a_r[2]))/(Ssl-Si);
      gn_a_r[k][2][0] -=
        (gn_r[k][2][1]*(vn_r[1]-vn_a_r[1])-gn_r[k][2][2]*(vn_r[2]-vn_a_r[2]))/(Ssl-Si);

      // rotate g back to original frame of reference
      g_a_r.push_back(tk::unrotateTensor(gn_a_r[k], fn));

      // compute pressure
      pm_a_r[k] = mat_blk[k].compute< EOS::pressure >(arho_a_r[k], v_a_r[0], v_a_r[1], v_a_r[2],
                                                      arhoe_a_r[k], al_a_r[k], k, g_a_r[k]);
      // compute sigma from equation of state
      asig_a_r = mat_blk[k].computeTensor< EOS::CauchyStress >(
        0.0, 0.0, 0.0, 0.0, 0.0, al_a_r[k], k, g_a_r[k] );
      for (std::size_t i=0; i<3; ++i) asig_a_r[i][i] -= pm_a_r[k];
      asign_a_r.push_back(tk::matvec(asig_a_r, fn));

      arhoe_a_r[k] = arhoe_t_r[k]
        + vn_a_r[1]*asig_a_r[1][0] - vn_r[1]*asig_r[k][1][0]
        + vn_a_r[2]*asig_a_r[2][0] - vn_r[2]*asig_r[k][2][0];}

    // Conservative flux functions
    // -------------------------------------------------------------------------
    for (std::size_t k=0; k<nmat; ++k)
    {
      // Left fluxes
      fl[volfracIdx(nmat, k)] = vn_l[0] * al_l[k];
      fl[densityIdx(nmat, k)] = vn_l[0] * u[0][densityIdx(nmat, k)];
      fl[energyIdx(nmat, k)] = vn_l[0] * u[0][energyIdx(nmat, k)];
      for (std::size_t i=0; i<3; ++i) {
        fl[energyIdx(nmat, k)] -= u[0][ncomp+velocityIdx(nmat,i)] *
          asign_l[k][i];
      }

      // inv deformation gradient tensor
      if (solidx[k] > 0) {
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            fl[deformIdx(nmat,solidx[k],i,j)] = (
              g_l[k][i][0] * ul +
              g_l[k][i][1] * vl +
              g_l[k][i][2] * wl ) * fn[j];
      }

      // Right fluxes
      fr[volfracIdx(nmat, k)] = vn_r[0] * al_r[k];
      fr[densityIdx(nmat, k)] = vn_r[0] * u[1][densityIdx(nmat, k)];
      fr[energyIdx(nmat, k)] = vn_r[0] * u[1][energyIdx(nmat, k)];
      for (std::size_t i=0; i<3; ++i) {
        fr[energyIdx(nmat, k)] -= u[1][ncomp+velocityIdx(nmat,i)] *
          asign_r[k][i];
      }

      // inv deformation gradient tensor
      if (solidx[k] > 0) {
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            fr[deformIdx(nmat,solidx[k],i,j)] = (
              g_r[k][i][0] * ur +
              g_r[k][i][1] * vr +
              g_r[k][i][2] * wr ) * fn[j];
      }

      // Left tilde fluxes
      ftl = fl;
      ftl[volfracIdx(nmat, k)] += Sl * (al_t_l[k] - al_l[k]);
      ftl[densityIdx(nmat, k)] += Sl * (arho_t_l[k] - u[0][densityIdx(nmat, k)]);
      ftl[energyIdx(nmat, k)] += Sl * (arhoe_t_l[k] - u[0][energyIdx(nmat, k)]);

      // inv deformation gradient tensor
      if (solidx[k] > 0) {
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            ftl[deformIdx(nmat,solidx[k],i,j)] +=
              Sl * (g_t_l[k][i][j] - g_l[k][i][j]);
      }

      // Right tilde fluxes
      ftr = fr;
      ftr[volfracIdx(nmat, k)] += Sr * (al_t_r[k] - al_r[k]);
      ftr[densityIdx(nmat, k)] += Sr * (arho_t_r[k] - u[1][densityIdx(nmat, k)]);
      ftr[energyIdx(nmat, k)] += Sr * (arhoe_t_r[k] - u[1][energyIdx(nmat, k)]);

      // inv deformation gradient tensor
      if (solidx[k] > 0) {
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            ftr[deformIdx(nmat,solidx[k],i,j)] +=
              Sr * (g_t_r[k][i][j] - g_r[k][i][j]);
      }

      // Left asterisk fluxes
      fal = ftl;
      fal[volfracIdx(nmat, k)] += Ssl * (al_a_l[k] - al_t_l[k]);
      fal[densityIdx(nmat, k)] += Ssl * (arho_a_l[k] - arho_t_l[k]);
      fal[energyIdx(nmat, k)] += Ssl * (arhoe_a_l[k] - arhoe_t_l[k]);

      // inv deformation gradient tensor
      if (solidx[k] > 0) {
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            fal[deformIdx(nmat,solidx[k],i,j)] +=
              Ssl * (g_a_l[k][i][j] - g_t_l[k][i][j]);
      }
      
      // Right asterisk fluxes
      far = ftr;
      far[volfracIdx(nmat, k)] += Ssr * (al_a_r[k] - al_t_r[k]);
      far[densityIdx(nmat, k)] += Ssr * (arho_a_r[k] - arho_t_r[k]);
      far[energyIdx(nmat, k)] += Ssr * (arhoe_a_r[k] - arhoe_t_r[k]);

      // inv deformation gradient tensor
      if (solidx[k] > 0) {
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            far[deformIdx(nmat,solidx[k],i,j)] +=
              Ssr * (g_a_r[k][i][j] - g_t_r[k][i][j]);
      }
    }

    // bulk momentum
    for (std::size_t idir=0; idir<3; ++idir)
    {
      fl[momentumIdx(nmat, idir)] = vn_l[0]*u[0][momentumIdx(nmat, idir)]
        - sign_l[idir];
      fr[momentumIdx(nmat, idir)] = vn_r[0]*u[1][momentumIdx(nmat, idir)]
        - sign_r[idir];
      ftl[momentumIdx(nmat, idir)] = fl[momentumIdx(nmat, idir)]
        + Sl * (rho_t_l*v_t_l[idir] - u[0][momentumIdx(nmat, idir)]);
      ftr[momentumIdx(nmat, idir)] = fr[momentumIdx(nmat, idir)]
        + Sr * (rho_t_r*v_t_r[idir] - u[1][momentumIdx(nmat, idir)]);
      fal[momentumIdx(nmat, idir)] = ftl[momentumIdx(nmat, idir)]
        + Ssl * (rho_a_l*v_a_l[idir] - rho_t_l*v_t_l[idir]);
      far[momentumIdx(nmat, idir)] = ftr[momentumIdx(nmat, idir)]
        + Ssr * (rho_a_r*v_a_r[idir] - rho_t_r*v_t_r[idir]);
    }

    // Numerical fluxes
    // -------------------------------------------------------------------------

    // Quick debug parameters
    bool dbg = false;
    std::vector< std::size_t > hll_dofs = {{ 0, 1, 2, 3, 4, 5 }};
      //{{ 6, 7, 8, 9, 10, 11, 12, 13, 14 }};

    // Numerical flux functions and wave-speeds
    if (Sl >= 0.0)
    {
      flx = fl;
      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k)
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
    else if (Sl <= 0.0 && 0.0 <= Ssl)
    {
      flx = ftl;
      if (dbg)
      {
        for (std::size_t k=0; k<hll_dofs.size(); ++k)
        {
          int kdof = hll_dofs[k];
          flx[kdof] = (Sr*fl[kdof] - Sl*fr[kdof]
                       + Sl*Sr*(u[1][kdof]-u[0][kdof])) / (Sr-Sl);
        }
        // Quantities for non-conservative terms
        auto c_plus(0.0), c_minus(0.0), p_plus(0.0), p_minus(0.0);
        c_plus = (Sr*vn_l[0] - Sr*Sl) / (Sr-Sl);
        c_minus = (Sr*Sl - Sl*vn_r[0]) / (Sr-Sl);
        p_plus = Sr / (Sr-Sl);
        p_minus = -Sl / (Sr-Sl);
        auto vriem = c_plus+c_minus;
        // Store Riemann-advected partial pressures
        for (std::size_t k=0; k<nmat; ++k)
          flx.push_back(p_plus*pml[k] + p_minus*pmr[k]);
        // Store Riemann velocity
        flx.push_back( vriem );
        // Store Riemann asign_ij (3*nsld)
        for (std::size_t k=0; k<nmat; ++k) {
          if (solidx[k] > 0) {
            for (std::size_t i=0; i<3; ++i)
              flx.push_back(p_plus*asign_l[k][i] + p_minus*asign_r[k][i]);
          }
        }
      }
      else
      {
        // Store Riemann-advected partial pressures
        for (std::size_t k=0; k<nmat; ++k)
          flx.push_back( 0.5*(pm_t_l[k]+pm_t_r[k]) );
        // Store Riemann velocity
        flx.push_back( vn_t_l[0] );
        // Store Riemann asign_ij (3*nsld)
        for (std::size_t k=0; k<nmat; ++k) {
          if (solidx[k] > 0) {
            for (std::size_t i=0; i<3; ++i)
              flx.push_back( asign_t_l[k][i] );
          }
        }
      }
    }
    else if (Ssl <= 0.0 && 0.0 <= Si)
    {
      flx = fal;
      if (dbg)
      {
        for (std::size_t k=0; k<hll_dofs.size(); ++k)
        {
          int kdof = hll_dofs[k];
          flx[kdof] = (Sr*fl[kdof] - Sl*fr[kdof]
                       + Sl*Sr*(u[1][kdof]-u[0][kdof])) / (Sr-Sl);
        }
        // Quantities for non-conservative terms
        auto c_plus(0.0), c_minus(0.0), p_plus(0.0), p_minus(0.0);
        c_plus = (Sr*vn_l[0] - Sr*Sl) / (Sr-Sl);
        c_minus = (Sr*Sl - Sl*vn_r[0]) / (Sr-Sl);
        p_plus = Sr / (Sr-Sl);
        p_minus = -Sl / (Sr-Sl);
        auto vriem = c_plus+c_minus;
        // Store Riemann-advected partial pressures
        for (std::size_t k=0; k<nmat; ++k)
          flx.push_back(p_plus*pml[k] + p_minus*pmr[k]);
        // Store Riemann velocity
        flx.push_back( vriem );
        // Store Riemann asign_ij (3*nsld)
        for (std::size_t k=0; k<nmat; ++k) {
          if (solidx[k] > 0) {
            for (std::size_t i=0; i<3; ++i)
              flx.push_back(p_plus*asign_l[k][i] + p_minus*asign_r[k][i]);
          }
        }
      }
      else
      {
        // Store Riemann-advected partial pressures
        for (std::size_t k=0; k<nmat; ++k)
          flx.push_back( 0.5*(pm_a_l[k]+pm_a_r[k]) );
        // Store Riemann velocity
        flx.push_back( vn_a_l[0] );
        // Store Riemann asign_ij (3*nsld)
        for (std::size_t k=0; k<nmat; ++k) {
          if (solidx[k] > 0) {
            for (std::size_t i=0; i<3; ++i)
              flx.push_back( asign_a_l[k][i] );
          }
        }
      }
    }
    else if (Si <= 0.0 && 0.0 <= Ssr)
    {
      flx = far;
      if (dbg)
      {
        for (std::size_t k=0; k<hll_dofs.size(); ++k)
        {
          int kdof = hll_dofs[k];
          flx[kdof] = (Sr*fl[kdof] - Sl*fr[kdof]
                       + Sl*Sr*(u[1][kdof]-u[0][kdof])) / (Sr-Sl);
        }
        // Quantities for non-conservative terms
        auto c_plus(0.0), c_minus(0.0), p_plus(0.0), p_minus(0.0);
        c_plus = (Sr*vn_l[0] - Sr*Sl) / (Sr-Sl);
        c_minus = (Sr*Sl - Sl*vn_r[0]) / (Sr-Sl);
        p_plus = Sr / (Sr-Sl);
        p_minus = -Sl / (Sr-Sl);
        auto vriem = c_plus+c_minus;
        // Store Riemann-advected partial pressures
        for (std::size_t k=0; k<nmat; ++k)
          flx.push_back(p_plus*pml[k] + p_minus*pmr[k]);
        // Store Riemann velocity
        flx.push_back( vriem );
        // Store Riemann asign_ij (3*nsld)
        for (std::size_t k=0; k<nmat; ++k) {
          if (solidx[k] > 0) {
            for (std::size_t i=0; i<3; ++i)
              flx.push_back(p_plus*asign_l[k][i] + p_minus*asign_r[k][i]);
          }
        }
      }
      else
      {
        // Store Riemann-advected partial pressures
        for (std::size_t k=0; k<nmat; ++k)
          flx.push_back( 0.5*(pm_a_r[k]+pm_a_r[k]) );
        // Store Riemann velocity
        flx.push_back( vn_a_r[0] );
        // Store Riemann asign_ij (3*nsld)
        for (std::size_t k=0; k<nmat; ++k) {
          if (solidx[k] > 0) {
            for (std::size_t i=0; i<3; ++i)
              flx.push_back( asign_a_r[k][i] );
          }
        }
      }
    }
    else if (Ssr <= 0.0 && 0.0 <= Sr)
    {
      flx = ftr;
      if (dbg)
      {
        for (std::size_t k=0; k<hll_dofs.size(); ++k)
        {
          int kdof = hll_dofs[k];
          flx[kdof] = (Sr*fl[kdof] - Sl*fr[kdof]
                       + Sl*Sr*(u[1][kdof]-u[0][kdof])) / (Sr-Sl);
        }
        // Quantities for non-conservative terms
        auto c_plus(0.0), c_minus(0.0), p_plus(0.0), p_minus(0.0);
        c_plus = (Sr*vn_l[0] - Sr*Sl) / (Sr-Sl);
        c_minus = (Sr*Sl - Sl*vn_r[0]) / (Sr-Sl);
        p_plus = Sr / (Sr-Sl);
        p_minus = -Sl / (Sr-Sl);
        auto vriem = c_plus+c_minus;
        // Store Riemann-advected partial pressures
        for (std::size_t k=0; k<nmat; ++k)
          flx.push_back(p_plus*pml[k] + p_minus*pmr[k]);
        // Store Riemann velocity
        flx.push_back( vriem );
        // Store Riemann asign_ij (3*nsld)
        for (std::size_t k=0; k<nmat; ++k) {
          if (solidx[k] > 0) {
            for (std::size_t i=0; i<3; ++i)
              flx.push_back(p_plus*asign_l[k][i] + p_minus*asign_r[k][i]);
          }
        }
      }
      else
      {
        // Store Riemann-advected partial pressures
        for (std::size_t k=0; k<nmat; ++k)
          flx.push_back( 0.5*(pm_t_r[k]+pm_t_r[k]) );
        // Store Riemann velocity
        flx.push_back( vn_t_r[0] );
        // Store Riemann asign_ij (3*nsld)
        for (std::size_t k=0; k<nmat; ++k) {
          if (solidx[k] > 0) {
            for (std::size_t i=0; i<3; ++i)
              flx.push_back( asign_t_r[k][i] );
          }
        }
      }
    }
    else
    {
      flx = fr;
      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back(pmr[k]);
      // Store Riemann velocity
      flx.push_back( vn_r[0] );
      // Store Riemann asign_ij (3*nsld)
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(asign_r[k][i]);
        }
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
