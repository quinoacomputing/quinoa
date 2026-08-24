// *****************************************************************************
/*!
  \file      src/PDE/Riemann/HLLCMultiMat.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2023 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     HLLC Riemann flux function for solids
  \details   This file implements the HLLC Riemann solver for solids.
             Ref. Ndanou, S., Favrie, N., & Gavrilyuk, S. (2015). Multi-solid
             and multi-fluid diffuse interface model: Applications to dynamic
             fracture and fragmentation. Journal of Computational Physics, 295,
             523-555.
*/
// *****************************************************************************
#ifndef HLLCMultiMat_h
#define HLLCMultiMat_h

#include <vector>

#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EOS.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "MultiMat/MiscMultiMatFns.hpp"

namespace inciter {

//! HLLC approximate Riemann solver for solids
struct HLLCMultiMat {

//! HLLC approximate Riemann solver flux function
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \param[in] wn Mesh velocity normal to face
  //! \return Riemann solution according to Harten-Lax-van Leer-Contact
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::vector< EOS >& mat_blk,
        const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::vector< std::array< tk::real, 3 > >& = {},
        const tk::real wn = 0 )
  {
    auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
    const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

    auto nsld = numSolids(nmat, solidx);
    auto ncomp = u[0].size()-(3+nmat+nsld*6);
    std::vector< tk::real > flx(ncomp, 0);

    // Primitive variables
    // -------------------------------------------------------------------------
    tk::real rhol(0.0), rhor(0.0);
    for (size_t k=0; k<nmat; ++k) {
      rhol += u[0][densityIdx(nmat, k)];
      rhor += u[1][densityIdx(nmat, k)];
    }

    auto ul = u[0][ncomp+velocityIdx(nmat, 0)];
    auto vl = u[0][ncomp+velocityIdx(nmat, 1)];
    auto wl = u[0][ncomp+velocityIdx(nmat, 2)];
    auto ur = u[1][ncomp+velocityIdx(nmat, 0)];
    auto vr = u[1][ncomp+velocityIdx(nmat, 1)];
    auto wr = u[1][ncomp+velocityIdx(nmat, 2)];

    // Firewall: sanitize velocities and densities
    if (!std::isfinite(ul)) ul = 0.0;
    if (!std::isfinite(vl)) vl = 0.0;
    if (!std::isfinite(wl)) wl = 0.0;
    if (!std::isfinite(ur)) ur = 0.0;
    if (!std::isfinite(vr)) vr = 0.0;
    if (!std::isfinite(wr)) wr = 0.0;
    if (!std::isfinite(rhol) || rhol <= 0.0) rhol = 1.0e-12;
    if (!std::isfinite(rhor) || rhor <= 0.0) rhor = 1.0e-12;

    // Outer states
    // -------------------------------------------------------------------------
    [[maybe_unused]] tk::real pl(0.0), pr(0.0);
    tk::real acl(0.0), acr(0.0);
    std::vector< tk::real > apl(nmat, 0.0), apr(nmat, 0.0);
    std::array< tk::real, 3 > Tnl{{0, 0, 0}}, Tnr{{0, 0, 0}};
    std::vector< std::array< tk::real, 3 > > aTnl, aTnr;
    std::array< std::array< tk::real, 3 >, 3 > asigl, asigr;
    std::array< std::array< tk::real, 3 >, 3 >
      signnl{{{0,0,0},{0,0,0},{0,0,0}}};
    std::array< std::array< tk::real, 3 >, 3 >
      signnr{{{0,0,0},{0,0,0},{0,0,0}}};
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gl, gr;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gnl, gnr;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > asignnl, asignnr;

    for (std::size_t k=0; k<nmat; ++k) {
      // Left state
      apl[k] = u[0][ncomp+pressureIdx(nmat, k)];
      // Sanitize pressure
      if (!std::isfinite(apl[k])) apl[k] = 0.0;
      pl += apl[k];

      // inv deformation gradient and Cauchy stress tensors
      gl.push_back(getDeformGrad(nmat, k, u[0]));
      asigl = getCauchyStress(nmat, k, ncomp, u[0]);
      for (std::size_t i=0; i<3; ++i) asigl[i][i] -= apl[k];

      // normal stress (traction) vector
      aTnl.push_back(tk::matvec(asigl, fn));
      // Sanitize traction vector components
      for (std::size_t i=0; i<3; ++i) {
        if (!std::isfinite(aTnl[k][i])) aTnl[k][i] = 0.0;
        Tnl[i] += aTnl[k][i];
      }

      // rotate stress vector
      asignnl.push_back(tk::rotateTensor(asigl, fn));
      // Sanitize rotated stress tensor
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j) {
          if (!std::isfinite(asignnl[k][i][j])) asignnl[k][i][j] = 0.0;
          signnl[i][j] += asignnl[k][i][j];
        }

      // rotate deformation gradient tensor for speed of sound in normal dir
      gnl.push_back(tk::rotateTensor(gl[k], fn));
      // Safe damage calculation: damage/density with protection
      auto rhol_k = u[0][densityIdx(nmat, k)];
      auto damagel_k = u[0][damageIdx(nmat, nsld, solidx[k])];
      auto damage_ratiol = (std::fabs(rhol_k) > 1.0e-12 && std::isfinite(rhol_k) && std::isfinite(damagel_k))
                           ? std::max(0.0, std::min(1.0, damagel_k/rhol_k))
                           : 0.0;
      auto amatl = mat_blk[k].compute< EOS::soundspeed >(
        rhol_k, apl[k],
        u[0][volfracIdx(nmat, k)], k, gnl[k],
        damage_ratiol );

      // Right state
      apr[k] = u[1][ncomp+pressureIdx(nmat, k)];
      // Sanitize pressure
      if (!std::isfinite(apr[k])) apr[k] = 0.0;
      pr += apr[k];

      // inv deformation gradient and Cauchy stress tensors
      gr.push_back(getDeformGrad(nmat, k, u[1]));
      asigr = getCauchyStress(nmat, k, ncomp, u[1]);
      for (std::size_t i=0; i<3; ++i) asigr[i][i] -= apr[k];

      // normal stress (traction) vector
      aTnr.push_back(tk::matvec(asigr, fn));
      // Sanitize traction vector components
      for (std::size_t i=0; i<3; ++i) {
        if (!std::isfinite(aTnr[k][i])) aTnr[k][i] = 0.0;
        Tnr[i] += aTnr[k][i];
      }

      // rotate stress vector
      asignnr.push_back(tk::rotateTensor(asigr, fn));
      // Sanitize rotated stress tensor
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j) {
          if (!std::isfinite(asignnr[k][i][j])) asignnr[k][i][j] = 0.0;
          signnr[i][j] += asignnr[k][i][j];
        }

      // rotate deformation gradient tensor for speed of sound in normal dir
      gnr.push_back(tk::rotateTensor(gr[k], fn));
      // Safe damage calculation: damage/density with protection
      auto rhor_k = u[1][densityIdx(nmat, k)];
      auto damager_k = u[1][damageIdx(nmat, nsld, solidx[k])];
      auto damage_ratior = (std::fabs(rhor_k) > 1.0e-12 && std::isfinite(rhor_k) && std::isfinite(damager_k))
                           ? std::max(0.0, std::min(1.0, damager_k/rhor_k))
                           : 0.0;
      auto amatr = mat_blk[k].compute< EOS::soundspeed >(
        rhor_k, apr[k],
        u[1][volfracIdx(nmat, k)], k, gnr[k],
        damage_ratior );

      // Mixture speed of sound
      if (std::isfinite(amatl)) {
        acl += u[0][densityIdx(nmat, k)] * amatl * amatl;
      }
      if (std::isfinite(amatr)) {
        acr += u[1][densityIdx(nmat, k)] * amatr * amatr;
      }
    }
    // Safe mixture sound speed with division protection
    rhol = std::max(1.0e-12, rhol);
    rhor = std::max(1.0e-12, rhor);
    // Ensure acl and acr are finite before division
    acl = (std::isfinite(acl) && acl > 0.0) ? std::sqrt(acl/rhol) : 0.0;
    acr = (std::isfinite(acr) && acr > 0.0) ? std::sqrt(acr/rhor) : 0.0;

    // Rotated velocities from advective velocities
    auto vnl = tk::rotateVector({ul, vl, wl}, fn);
    auto vnr = tk::rotateVector({ur, vr, wr}, fn);

    // ALE mesh motion relative normal velocity
    vnl[0] -= wn;
    vnr[0] -= wn;

    // Signal velocities
    auto Sl = std::min((vnl[0]-acl), (vnr[0]-acr));
    auto Sr = std::max((vnl[0]+acl), (vnr[0]+acr));
    // Safe Sm calculation with division protection
    auto Sm_numer = rhor*vnr[0]*(Sr-vnr[0]) - rhol*vnl[0]*(Sl-vnl[0])
                    - signnl[0][0] + signnr[0][0];
    auto Sm_denom = rhor*(Sr-vnr[0]) - rhol*(Sl-vnl[0]);
    auto Sm = (std::fabs(Sm_denom) > 1.0e-12 && std::isfinite(Sm_numer) && std::isfinite(Sm_denom))
              ? Sm_numer / Sm_denom
              : 0.5 * (vnl[0] + vnr[0]);  // Fallback to average velocity

    // Middle-zone (star) variables
    // -------------------------------------------------------------------------
    // the stress in the star zone is theoretically equivalent when derived
    // from the left or right state. However, there might be differences
    // numerically due to truncation etc. Hence two separate evaluations
    // are used.
    std::vector< std::array< std::array< tk::real, 3 >, 3 > >
      asignnlStar, asignnrStar;
    std::array< std::array< tk::real, 3 >, 3 >
      signnlStar{{{0,0,0},{0,0,0},{0,0,0}}},
      signnrStar{{{0,0,0},{0,0,0},{0,0,0}}};
    std::array< std::array< tk::real, 3 >, 3 >
      siglStar{{{0,0,0},{0,0,0},{0,0,0}}}, sigrStar{{{0,0,0},{0,0,0},{0,0,0}}};
    asignnlStar.resize(nmat);
    asignnrStar.resize(nmat);
    std::array< tk::real, 3 > TnlStar{{0, 0, 0}}, TnrStar{{0, 0, 0}};
    std::vector< std::array< tk::real, 3 > > aTnlStar, aTnrStar;
    for (std::size_t k=0; k<nmat; ++k) {
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
        {
          asignnlStar[k][i][j] = asignnl[k][i][j];
          asignnrStar[k][i][j] = asignnr[k][i][j];
        }

      for (std::size_t i=0; i<1; ++i)
      {
        asignnlStar[k][i][i] -=
          u[0][densityIdx(nmat,k)]*(Sl-vnl[0])*(Sm-vnl[0]);
        asignnrStar[k][i][i] -=
          u[1][densityIdx(nmat,k)]*(Sr-vnr[0])*(Sm-vnr[0]);
      }

      for (std::size_t i=1; i<3; ++i)
      {
        // Safe star state shear stress calculation with division protection
        auto star_numer = (Sm-Sl)*u[0][densityIdx(nmat,k)]*asignnr[k][i][0]
                        - (Sm-Sr)*u[1][densityIdx(nmat,k)]*asignnl[k][i][0]
                        + (Sm-Sl)*u[0][densityIdx(nmat,k)]
                        * (Sm-Sr)*u[1][densityIdx(nmat,k)]
                        * (vnl[i]-vnr[i]);
        auto star_denom = (Sm-Sl)*u[0][densityIdx(nmat,k)]
                        - (Sm-Sr)*u[1][densityIdx(nmat,k)];

        if (std::fabs(star_denom) > 1.0e-12 && std::isfinite(star_numer) && std::isfinite(star_denom)) {
          asignnlStar[k][i][0] = star_numer / star_denom;
          asignnrStar[k][i][0] = star_numer / star_denom;
        } else {
          // Fallback to original values if division is unsafe
          asignnlStar[k][i][0] = asignnl[k][i][0];
          asignnrStar[k][i][0] = asignnr[k][i][0];
        }
      }
      // Symmetry
      asignnlStar[k][0][1] = asignnlStar[k][1][0];
      asignnlStar[k][0][2] = asignnlStar[k][2][0];
      asignnrStar[k][0][1] = asignnrStar[k][1][0];
      asignnrStar[k][0][2] = asignnrStar[k][2][0];

      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
        {
          signnlStar[i][j] += asignnlStar[k][i][j];
          signnrStar[i][j] += asignnrStar[k][i][j];
        }

      auto asiglStar = tk::unrotateTensor(asignnlStar[k], fn);
      auto asigrStar = tk::unrotateTensor(asignnrStar[k], fn);
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
        {
          siglStar[i][j] += asiglStar[i][j];
          sigrStar[i][j] += asigrStar[i][j];
        }
      aTnlStar.push_back(tk::matvec(asiglStar, fn));
      aTnrStar.push_back(tk::matvec(asigrStar, fn));
      for (std::size_t i=0; i<3; ++i)
      {
        TnlStar[i] += aTnlStar[k][i];
        TnrStar[i] += aTnrStar[k][i];
      }
    }

    // Safe weight calculations with division protection
    auto w_l_denom = Sm - Sl;
    auto w_r_denom = Sm - Sr;
    auto w_l = (std::fabs(w_l_denom) > 1.0e-12 && std::isfinite(w_l_denom))
               ? (vnl[0]-Sl)/w_l_denom : 1.0;
    auto w_r = (std::fabs(w_r_denom) > 1.0e-12 && std::isfinite(w_r_denom))
               ? (vnr[0]-Sr)/w_r_denom : 1.0;

    std::array< tk::real, 3 > vnlStar, vnrStar;
    // u*_L with safe divisions
    vnlStar[0] = Sm;
    auto denom_l = w_l*rhol*(vnl[0]-Sl);
    if (std::fabs(denom_l) > 1.0e-12 && std::isfinite(denom_l)) {
      vnlStar[1] = vnl[1] + (signnlStar[1][0] - signnl[1][0]) / denom_l;
      vnlStar[2] = vnl[2] + (signnlStar[2][0] - signnl[2][0]) / denom_l;
    } else {
      vnlStar[1] = vnl[1];
      vnlStar[2] = vnl[2];
    }
    // u*_R with safe divisions
    vnrStar[0] = Sm;
    auto denom_r = w_r*rhor*(vnr[0]-Sr);
    if (std::fabs(denom_r) > 1.0e-12 && std::isfinite(denom_r)) {
      vnrStar[1] = vnr[1] + (signnrStar[1][0] - signnr[1][0]) / denom_r;
      vnrStar[2] = vnr[2] + (signnrStar[2][0] - signnr[2][0]) / denom_r;
    } else {
      vnrStar[1] = vnr[1];
      vnrStar[2] = vnr[2];
    }

    auto vlStar = tk::unrotateVector(vnlStar, fn);
    auto vrStar = tk::unrotateVector(vnrStar, fn);

    const auto unl = vnl[0] + wn;
    const auto unr = vnr[0] + wn;
    const auto unStar = Sm + wn;

    auto uStar = u;

    tk::real rholStar(0.0), rhorStar(0.0);
    std::array< std::array< tk::real, 3 >, 3 > tempArray {{ {0,0,0}, {0,0,0}, {0,0,0} }};
    std::vector< std::array< std::array< tk::real, 3 >, 3 > >
      gnlStar(nmat, tempArray), gnrStar(nmat, tempArray),
      glStar(nmat, tempArray), grStar(nmat, tempArray);

    // Pre-compute safe inverse of wave speed differences
    auto inv_Sm_Sl = (std::fabs(Sm-Sl) > 1.0e-12 && std::isfinite(Sm-Sl))
                     ? 1.0/(Sm-Sl) : 0.0;
    auto inv_Sm_Sr = (std::fabs(Sm-Sr) > 1.0e-12 && std::isfinite(Sm-Sr))
                     ? 1.0/(Sm-Sr) : 0.0;

    for (std::size_t k=0; k<nmat; ++k) {
      // Left
      gnlStar[k] = gnl[k];
      if (solidx[k] > 0)
      {
        gnlStar[k][0][0] = w_l * gnl[k][0][0]
          + gnl[k][0][1]*(vnl[1]-vnlStar[1])*inv_Sm_Sl
          + gnl[k][0][2]*(vnl[2]-vnlStar[2])*inv_Sm_Sl;
        gnlStar[k][1][0] = w_l * gnl[k][1][0]
          + gnl[k][1][1]*(vnl[1]-vnlStar[1])*inv_Sm_Sl
          + gnl[k][1][2]*(vnl[2]-vnlStar[2])*inv_Sm_Sl;
        gnlStar[k][2][0] = w_l * gnl[k][2][0]
          + gnl[k][2][1]*(vnl[1]-vnlStar[1])*inv_Sm_Sl
          + gnl[k][2][2]*(vnl[2]-vnlStar[2])*inv_Sm_Sl;
        // rotate g back to original frame of reference
        glStar.push_back(tk::unrotateTensor(gnlStar[k], fn));
        // damage
        uStar[0][damageIdx(nmat, nsld, solidx[k])] = w_l * u[0][damageIdx(nmat, nsld, solidx[k])];
      }
      // rotate g back to original frame of reference
      glStar[k] = tk::unrotateTensor(gnlStar[k], fn);
      uStar[0][volfracIdx(nmat, k)] = u[0][volfracIdx(nmat, k)];
      uStar[0][densityIdx(nmat, k)] = w_l * u[0][densityIdx(nmat, k)];
      uStar[0][energyIdx(nmat, k)] = w_l * u[0][energyIdx(nmat, k)]
        + ( - asignnl[k][0][0]*unl
            - asignnl[k][1][0]*vnl[1]
            - asignnl[k][2][0]*vnl[2]
            + asignnlStar[k][0][0]*unStar
            + asignnlStar[k][1][0]*vnlStar[1]
            + asignnlStar[k][2][0]*vnlStar[2]
          ) * inv_Sm_Sl;
      rholStar += uStar[0][densityIdx(nmat, k)];

      // Right
      gnrStar[k] = gnr[k];
      if (solidx[k] > 0)
      {
        gnrStar[k][0][0] = w_r * gnr[k][0][0]
          + gnr[k][0][1]*(vnr[1]-vnrStar[1])*inv_Sm_Sr
          + gnr[k][0][2]*(vnr[2]-vnrStar[2])*inv_Sm_Sr;
        gnrStar[k][1][0] = w_r * gnr[k][1][0]
          + gnr[k][1][1]*(vnr[1]-vnrStar[1])*inv_Sm_Sr
          + gnr[k][1][2]*(vnr[2]-vnrStar[2])*inv_Sm_Sr;
        gnrStar[k][2][0] = w_r * gnr[k][2][0]
          + gnr[k][2][1]*(vnr[1]-vnrStar[1])*inv_Sm_Sr
          + gnr[k][2][2]*(vnr[2]-vnrStar[2])*inv_Sm_Sr;
        // rotate g back to original frame of reference
        grStar.push_back(tk::unrotateTensor(gnrStar[k], fn));
        // damage
        uStar[1][damageIdx(nmat, nsld, solidx[k])] = w_r * u[1][damageIdx(nmat, nsld, solidx[k])];
      }
      // rotate g back to original frame of reference
      grStar[k] = tk::unrotateTensor(gnrStar[k], fn);
      uStar[1][volfracIdx(nmat, k)] = u[1][volfracIdx(nmat, k)];
      uStar[1][densityIdx(nmat, k)] = w_r * u[1][densityIdx(nmat, k)];
      uStar[1][energyIdx(nmat, k)] = w_r * u[1][energyIdx(nmat, k)]
        + ( - asignnr[k][0][0]*unr
            - asignnr[k][1][0]*vnr[1]
            - asignnr[k][2][0]*vnr[2]
            + asignnrStar[k][0][0]*unStar
            + asignnrStar[k][1][0]*vnrStar[1]
            + asignnrStar[k][2][0]*vnrStar[2]
          ) * inv_Sm_Sr;
      rhorStar += uStar[1][densityIdx(nmat, k)];
    }
    // Safe momentum star state: reuse existing safe inversions
    // Note: 1/(Sl-Sm) = -1/(Sm-Sl) = -inv_Sm_Sl, same for Sr
    for (std::size_t idir=0; idir<3; ++idir) {
      uStar[0][momentumIdx(nmat, idir)] = w_l*u[0][momentumIdx(nmat, idir)]
        + (TnlStar[idir] - Tnl[idir])*inv_Sm_Sl;
      uStar[1][momentumIdx(nmat, idir)] = w_r*u[1][momentumIdx(nmat, idir)]
        + (TnrStar[idir] - Tnr[idir])*inv_Sm_Sr;
    }

    // Numerical fluxes
    // -------------------------------------------------------------------------
    if (0.0 <= Sl) {

      for (std::size_t idir=0; idir<3; ++idir)
        flx[momentumIdx(nmat, idir)] =
          u[0][momentumIdx(nmat, idir)] * vnl[0] - Tnl[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = u[0][volfracIdx(nmat, k)] * vnl[0];
        flx[densityIdx(nmat, k)] = u[0][densityIdx(nmat, k)] * vnl[0];
        flx[energyIdx(nmat, k)] = u[0][energyIdx(nmat, k)] * vnl[0]
          - ul * aTnl[k][0] - vl * aTnl[k][1] - wl * aTnl[k][2];
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx[deformIdx(nmat,solidx[k],i,j)] = (
                gl[k][i][0] * ul +
                gl[k][i][1] * vl +
                gl[k][i][2] * wl ) * fn[j]
                - wn * gl[k][i][j];
          flx[damageIdx(nmat, nsld, solidx[k])] = u[0][damageIdx(nmat, nsld, solidx[k])] * vnl[0];
        }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k) {
        auto mag_sq = aTnl[k][0]*aTnl[k][0] + aTnl[k][1]*aTnl[k][1] + aTnl[k][2]*aTnl[k][2];
        flx.push_back(std::isfinite(mag_sq) && mag_sq >= 0.0 ? std::sqrt(mag_sq) : 0.0);
      }
      // Store Riemann velocity
      flx.push_back((vnl[0]+wn));
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnl[k][i]);
        }
      }
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx.push_back(gl[k][i][j]);
        }
      }

    }

    else if (Sl < 0.0 && 0.0 <= Sm) {

      std::array< tk::real, 3 > ulStarPhysical{{
          vlStar[0] + wn*fn[0],
          vlStar[1] + wn*fn[1],
          vlStar[2] + wn*fn[2]
        }};

      for (std::size_t idir=0; idir<3; ++idir)
       flx[momentumIdx(nmat, idir)] =
         uStar[0][momentumIdx(nmat, idir)] * Sm - TnlStar[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = uStar[0][volfracIdx(nmat, k)] * Sm;
        flx[densityIdx(nmat, k)] = uStar[0][densityIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] = uStar[0][energyIdx(nmat, k)] * Sm
          - ulStarPhysical[0] * aTnlStar[k][0]
          - ulStarPhysical[1] * aTnlStar[k][1]
          - ulStarPhysical[2] * aTnlStar[k][2];
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx[deformIdx(nmat,solidx[k],i,j)] =
                ( glStar[k][i][0] * ulStarPhysical[0]
                + glStar[k][i][1] * ulStarPhysical[1]
                + glStar[k][i][2] * ulStarPhysical[2] ) * fn[j]
                - wn * glStar[k][i][j];
          flx[damageIdx(nmat, nsld, solidx[k])] = uStar[0][damageIdx(nmat, nsld, solidx[k])] * Sm;
        }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k) {
        auto mag_sq = aTnlStar[k][0]*aTnlStar[k][0] + aTnlStar[k][1]*aTnlStar[k][1] + aTnlStar[k][2]*aTnlStar[k][2];
        flx.push_back(std::isfinite(mag_sq) && mag_sq >= 0.0 ? std::sqrt(mag_sq) : 0.0);
      }
      // Store Riemann velocity
      flx.push_back(Sm+wn);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnlStar[k][i]);
        }
      }
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx.push_back(glStar[k][i][j]);
        }
      }

    }

    else if (Sm < 0.0 && 0.0 <= Sr) {

      std::array< tk::real, 3 > urStarPhysical{{
          vrStar[0] + wn*fn[0],
          vrStar[1] + wn*fn[1],
          vrStar[2] + wn*fn[2]
        }};

      for (std::size_t idir=0; idir<3; ++idir)
        flx[momentumIdx(nmat, idir)] =
          uStar[1][momentumIdx(nmat, idir)] * Sm - TnrStar[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = uStar[1][volfracIdx(nmat, k)] * Sm;
        flx[densityIdx(nmat, k)] = uStar[1][densityIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] = uStar[1][energyIdx(nmat, k)] * Sm
          - urStarPhysical[0] * aTnrStar[k][0]
          - urStarPhysical[1] * aTnrStar[k][1]
          - urStarPhysical[2] * aTnrStar[k][2];
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
              for (std::size_t j=0; j<3; ++j)
                flx[deformIdx(nmat,solidx[k],i,j)] =
                  ( grStar[k][i][0] * urStarPhysical[0]
                  + grStar[k][i][1] * urStarPhysical[1]
                  + grStar[k][i][2] * urStarPhysical[2] ) * fn[j]
                  - wn * grStar[k][i][j];
          flx[damageIdx(nmat, nsld, solidx[k])] = uStar[1][damageIdx(nmat, nsld, solidx[k])] * Sm;
        }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k) {
        auto mag_sq = aTnrStar[k][0]*aTnrStar[k][0] + aTnrStar[k][1]*aTnrStar[k][1] + aTnrStar[k][2]*aTnrStar[k][2];
        flx.push_back(std::isfinite(mag_sq) && mag_sq >= 0.0 ? std::sqrt(mag_sq) : 0.0);
      }
      // Store Riemann velocity
      flx.push_back(Sm+wn);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnrStar[k][i]);
        }
      }
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx.push_back(grStar[k][i][j]);
        }
      }

    }

    else {

      for (std::size_t idir=0; idir<3; ++idir)
        flx[momentumIdx(nmat, idir)] =
          u[1][momentumIdx(nmat, idir)] * vnr[0] - Tnr[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = u[1][volfracIdx(nmat, k)] * vnr[0];
        flx[densityIdx(nmat, k)] = u[1][densityIdx(nmat, k)] * vnr[0];
        flx[energyIdx(nmat, k)] = u[1][energyIdx(nmat, k)] * vnr[0]
          - ur * aTnr[k][0] - vr * aTnr[k][1] - wr * aTnr[k][2];
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
              for (std::size_t j=0; j<3; ++j)
                flx[deformIdx(nmat,solidx[k],i,j)] = (
                  gr[k][i][0] * ur +
                  gr[k][i][1] * vr +
                  gr[k][i][2] * wr ) * fn[j]
                  - wn * gr[k][i][j];
          flx[damageIdx(nmat, nsld, solidx[k])] = u[1][damageIdx(nmat, nsld, solidx[k])] * vnr[0];
        }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k) {
        auto mag_sq = aTnr[k][0]*aTnr[k][0] + aTnr[k][1]*aTnr[k][1] + aTnr[k][2]*aTnr[k][2];
        flx.push_back(std::isfinite(mag_sq) && mag_sq >= 0.0 ? std::sqrt(mag_sq) : 0.0);
      }
      // Store Riemann velocity
      flx.push_back(vnr[0]+wn);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnr[k][i]);
        }
      }
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx.push_back(gr[k][i][j]);
        }
      }

    }

    Assert( flx.size() == (ncomp+nmat+1+3*nsld+9*nsld), "Size of "
            "multi-material flux vector incorrect" );

    return flx;
  }

  ////! Flux type accessor
  ////! \return Flux type
  //static ctr::FluxType type() noexcept {
  //  return ctr::FluxType::HLLCMultiMat; }
};

} // inciter::

#endif // HLLCMultiMat_h
