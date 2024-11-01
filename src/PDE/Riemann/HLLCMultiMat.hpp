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
  //! \return Riemann solution according to Harten-Lax-van Leer-Contact
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::vector< EOS >& mat_blk,
        const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::vector< std::array< tk::real, 3 > >& = {} )
  {
    auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
    const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

    auto nsld = numSolids(nmat, solidx);
    auto ncomp = u[0].size()-(3+nmat+nsld*6);
    std::vector< tk::real > flx(ncomp, 0), fl(ncomp, 0), fr(ncomp, 0),
      ftl(ncomp, 0), ftr(ncomp, 0);

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

    // Outer states
    // -------------------------------------------------------------------------
    tk::real pl(0.0), pr(0.0);
    tk::real acl(0.0), acr(0.0);
    std::vector< tk::real > apl(nmat, 0.0), apr(nmat, 0.0);
    std::array< tk::real, 3 > Tnl{{0, 0, 0}}, Tnr{{0, 0, 0}};
    std::vector< std::array< tk::real, 3 > > aTnl, aTnr;
    std::array< std::array< tk::real, 3 >, 3 > asigl, asigr;
    std::array< std::array< tk::real, 3 >, 3 > signnl{{{0,0,0},{0,0,0},{0,0,0}}};
    std::array< std::array< tk::real, 3 >, 3 > signnr{{{0,0,0},{0,0,0},{0,0,0}}};
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gl, gr;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gnl, gnr;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > asignnl, asignnr;

    for (std::size_t k=0; k<nmat; ++k) {
      // Left state
      apl[k] = u[0][ncomp+pressureIdx(nmat, k)];
      pl += apl[k];

      // inv deformation gradient and Cauchy stress tensors
      gl.push_back(getDeformGrad(nmat, k, u[0]));
      asigl = getCauchyStress(nmat, k, ncomp, u[0]);
      // for (std::size_t i=0; i<3; ++i) asigl[i][i] -= apl[k];

      // normal stress (traction) vector
      aTnl.push_back(tk::matvec(asigl, fn));
      for (std::size_t i=0; i<3; ++i)
        Tnl[i] += aTnl[k][i];

      // rotate stress vector
      asignnl.push_back(tk::rotateTensor(asigl, fn));
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          signnl[i][j] += asignnl[k][i][j];

      // rotate deformation gradient tensor for speed of sound in normal dir
      gnl.push_back(tk::rotateTensor(gl[k], fn));
      auto amatl = mat_blk[k].compute< EOS::soundspeed >(
        u[0][densityIdx(nmat, k)], apl[k], u[0][volfracIdx(nmat, k)], k, gnl[k] );

      // Right state
      apr[k] = u[1][ncomp+pressureIdx(nmat, k)];
      pr += apr[k];

      // inv deformation gradient and Cauchy stress tensors
      gr.push_back(getDeformGrad(nmat, k, u[1]));
      asigr = getCauchyStress(nmat, k, ncomp, u[1]);
      // for (std::size_t i=0; i<3; ++i) asigr[i][i] -= apr[k];

      // normal stress (traction) vector
      aTnr.push_back(tk::matvec(asigr, fn));
      for (std::size_t i=0; i<3; ++i)
        Tnr[i] += aTnr[k][i];

      // rotate stress vector
      asignnr.push_back(tk::rotateTensor(asigr, fn));
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          signnr[i][j] += asignnr[k][i][j];

      // rotate deformation gradient tensor for speed of sound in normal dir
      gnr.push_back(tk::rotateTensor(gr[k], fn));
      auto amatr = mat_blk[k].compute< EOS::soundspeed >(
        u[1][densityIdx(nmat, k)], apr[k], u[1][volfracIdx(nmat, k)], k, gnr[k] );

      // Mixture speed of sound
      acl += u[0][densityIdx(nmat, k)] * amatl * amatl;
      acr += u[1][densityIdx(nmat, k)] * amatr * amatr;
    }
    acl = std::sqrt(acl/rhol);
    acr = std::sqrt(acr/rhor);

    // Rotated velocities from advective velocities
    auto vnl = tk::rotateVector({ul, vl, wl}, fn);
    auto vnr = tk::rotateVector({ur, vr, wr}, fn);

    // Signal velocities
    auto Sl = std::min((vnl[0]-acl), (vnr[0]-acr));
    auto Sr = std::max((vnl[0]+acl), (vnr[0]+acr));
    auto Sm = ( rhor*vnr[0]*(Sr-vnr[0]) - rhol*vnl[0]*(Sl-vnl[0])
                - signnl[0][0] + pl + signnr[0][0] - pr )
              /( rhor*(Sr-vnr[0]) - rhol*(Sl-vnl[0]) );

    // Middle-zone (star) variables
    // -------------------------------------------------------------------------
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > asignnStar;
    std::array< std::array< tk::real, 3 >, 3 >
      signnStar{{{0,0,0},{0,0,0},{0,0,0}}};
    std::array< std::array< tk::real, 3 >, 3 > sigStar{{{0,0,0},{0,0,0},{0,0,0}}};
    asignnStar.resize(nmat);
    std::array< tk::real, 3 > TnStar{{0, 0, 0}};
    std::vector< std::array< tk::real, 3 > > aTnStar;
    std::vector< tk::real > apStar(nmat, 0.0);
    tk::real pStar(0.0);
    for (std::size_t k=0; k<nmat; ++k) {
      // compute pressure
      apStar[k] = apl[k] + u[0][densityIdx(nmat, k)]*(vnl[0]-Sl)*(vnl[0]-Sm);
      pStar += apStar[k];
      // // compute asigma* (in normal direction
      // for (std::size_t i=0; i<3; ++i)
      //   for (std::size_t j=0; j<3; ++j) {
      //     asignnStar[k][i][j] = ( (vnr[0]-Sr)*u[1][densityIdx(nmat,k)]*asignnl[k][i][j]
      //                           - (vnl[0]-Sl)*u[0][densityIdx(nmat,k)]*asignnr[k][i][j]
      //                           + (vnl[0]-Sl)*u[0][densityIdx(nmat,k)]
      //                           * (vnr[0]-Sr)*u[1][densityIdx(nmat,k)]
      //                             * (vnr[j]-vnl[j]) ) /
      //                           ( (vnr[0]-Sr)*u[1][densityIdx(nmat,k)] -
      //                           (vnl[0]-Sl)*u[0][densityIdx(nmat,k)]);
      //     signnStar[i][j] += asignnStar[k][i][j];
      //   }
      // for (std::size_t i=0; i<3; ++i)
      //   for (std::size_t j=0; j<3; ++j)
      //     asignnStar[k][i][j] = 0.0;
      // for (std::size_t j=0; j<3; ++j)
      //   asignnStar[k][0][j] = ( (vnr[0]-Sr)*u[1][densityIdx(nmat,k)]*asignnl[k][0][j]
      //                           - (vnl[0]-Sl)*u[0][densityIdx(nmat,k)]*asignnr[k][0][j]
      //                           + (vnl[0]-Sl)*u[0][densityIdx(nmat,k)]
      //                           * (vnr[0]-Sr)*u[1][densityIdx(nmat,k)]
      //                           * (vnr[0]-vnl[0]) ) /
      //                           ( (vnr[0]-Sr)*u[1][densityIdx(nmat,k)] -
      //                             (vnl[0]-Sl)*u[0][densityIdx(nmat,k)]);
      // // asignnStar[k][1][0] = asignnStar[k][0][1];
      // // asignnStar[k][2][0] = asignnStar[k][0][2];
      // for (std::size_t i=0; i<3; ++i)
      //   for (std::size_t j=0; j<3; ++j)
      //     signnStar[i][j] += asignnStar[k][i][j];
      // Using the equation of state (On the WRONG state!)
      std::array< std::array< tk::real, 3 >, 3 >
        gavg{{{0,0,0},{0,0,0},{0,0,0}}};
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          gavg[i][j] = 0.5*(gl[k][i][j]+gr[k][i][j]);
      asignnStar.push_back(mat_blk[k].computeTensor< EOS::CauchyStress >(
        0.0, 0.0, 0.0, 0.0, 0.0, u[0][volfracIdx(nmat, k)], k, gavg));
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          signnStar[i][j] += asignnStar[k][i][j];

      auto asigStar = tk::unrotateTensor(asignnStar[k], fn);
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          sigStar[i][j] += asigStar[i][j];
      aTnStar.push_back(tk::matvec(asigStar, fn));
      for (std::size_t i=0; i<3; ++i)
        TnStar[i] += aTnStar[k][i];
    }

    // if (std::abs(pStar*fn[0]+TnStar[0]) > 1.0e-02) {
    //   printf("DEBUGGING\n");
    //   printf("vnl = [%e, %e, %e]\n", vnl[0], vnl[1], vnl[2]);
    //   printf("vnr = [%e, %e, %e]\n", vnr[0], vnr[1], vnr[2]);
    //   printf("vnr-vnl = [%e, %e, %e]\n", vnr[0]-vnl[0], vnr[1]-vnl[1], vnr[2]-vnl[2]);
    //   printf("asignnl[0][0][:] = [%e, %e, %e]\n", asignnl[0][0][0],
    //          asignnl[0][0][1], asignnl[0][0][2]);
    //   printf("asignnl[0][1][:] = [%e, %e, %e]\n", asignnl[0][1][0],
    //          asignnl[0][1][1], asignnl[0][1][2]);
    //   printf("asignnl[0][2][:] = [%e, %e, %e]\n", asignnl[0][2][0],
    //          asignnl[0][2][1], asignnl[0][2][2]);
    //   printf("asignnr[0][0][:] = [%e, %e, %e]\n", asignnr[0][0][0],
    //          asignnr[0][0][1], asignnr[0][0][2]);
    //   printf("asignnr[0][1][:] = [%e, %e, %e]\n", asignnr[0][1][0],
    //          asignnr[0][1][1], asignnr[0][1][2]);
    //   printf("asignnr[0][2][:] = [%e, %e, %e]\n", asignnr[0][2][0],
    //          asignnr[0][2][1], asignnr[0][2][2]);
    //   printf("asignnStar[0][0][:] = [%e, %e, %e]\n", asignnStar[0][0][0],
    //          asignnStar[0][0][1], asignnStar[0][0][2]);
    //   printf("asignnStar[0][1][:] = [%e, %e, %e]\n", asignnStar[0][1][0],
    //          asignnStar[0][1][1], asignnStar[0][1][2]);
    //   printf("asignnStar[0][2][:] = [%e, %e, %e]\n", asignnStar[0][2][0],
    //          asignnStar[0][2][1], asignnStar[0][2][2]);
    //   printf("term1, term2, term3a, term3b, term3c, term3d, term3e = %e, %e, %e, %e, %e, %e, %e\n",
    //          (vnr[0]-Sr)*u[1][densityIdx(nmat,0)]*asignnl[0][0][0],
    //          (vnl[0]-Sl)*u[0][densityIdx(nmat,0)]*asignnr[0][0][0],
    //          (vnl[0]-Sl),
    //          u[0][densityIdx(nmat,0)],
    //          (vnr[0]-Sr),
    //          u[1][densityIdx(nmat,0)],
    //          (vnr[0]-vnl[0]) );
    //   printf("TnStar = [%e, %e, %e]\n", TnStar[0], TnStar[1], TnStar[2]);
    //   printf("-pStar*n = [%e, %e, %e]\n", -pStar*fn[0], -pStar*fn[1], -pStar*fn[2]);
    // }

    std::array< tk::real, 3 > vnlStar, vnrStar;
    // u*_L
    vnlStar[0] = Sm;
    vnlStar[1] = vnl[1] + (signnStar[0][1] - signnl[0][1])
      / ((vnl[0]-Sl)*rhol);
    vnlStar[2] = vnl[2] + (signnStar[0][2] - signnl[0][2])
      / ((vnl[0]-Sl)*rhol);
    // u*_R
    vnrStar[0] = Sm;
    vnrStar[1] = vnr[1] + (signnStar[0][1] - signnr[0][1])
      / ((vnr[0]-Sr)*rhor);
    vnrStar[2] = vnr[2] + (signnStar[0][2] - signnr[0][2])
      / ((vnr[0]-Sr)*rhor);

    auto vlStar = tk::unrotateVector(vnlStar, fn);
    auto vrStar = tk::unrotateVector(vnrStar, fn);

    auto uStar = u;

    tk::real rholStar(0.0), rhorStar(0.0);
    std::array< std::array< tk::real, 3 >, 3 > gnlStar{{{0,0,0},{0,0,0},{0,0,0}}};
    std::array< std::array< tk::real, 3 >, 3 > gnrStar{{{0,0,0},{0,0,0},{0,0,0}}};
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > glStar, grStar;
    for (std::size_t k=0; k<nmat; ++k) {    
      // Left
      auto w_l = (Sl-vnl[0])/(Sl-Sm);
      if (solidx[k] > 0)
      {
        gnlStar = gnl[k];
        for (std::size_t i=0; i<3; ++i)
          gnlStar[i][0] = w_l*gnl[k][i][0] + (
            gnl[k][i][1]*(vnl[1]-vnlStar[1]) + gnl[k][i][2]*(vnl[2]-vnlStar[2])
            )/(Sm-Sl);
        // rotate g back to original frame of reference
        glStar.push_back(tk::unrotateTensor(gnlStar, fn));
        // for (std::size_t i=0; i<3; ++i)
        //   for (std::size_t j=0; j<3; ++j)
        //     glStar[k][i][j] = w_l*gl[k][i][j];
      }
      uStar[0][volfracIdx(nmat, k)] = u[0][volfracIdx(nmat, k)];
      uStar[0][densityIdx(nmat, k)] =
        (Sl-vnl[0]) * u[0][densityIdx(nmat, k)] / (Sl-Sm);
      // uStar[0][energyIdx(nmat, k)] =
      //   ((Sl-vnl[0]) * u[0][energyIdx(nmat, k)] - apl[k]*vnl[0] + apStar[k]*Sm) / (Sl-Sm);
      uStar[0][energyIdx(nmat, k)] = u[0][energyIdx(nmat, k)] * w_l + (
        - (asignnl[k][0][0]-apl[k])*vnl[0]
        - asignnl[k][0][1]*vnl[1]
        - asignnl[k][0][2]*vnl[2]
        + (asignnStar[k][0][0]-apStar[k])*vnlStar[0]
        + asignnStar[k][0][1]*vnlStar[1]
        + asignnStar[k][0][2]*vnlStar[2]
        ) / (Sm-Sl);
      rholStar += uStar[0][densityIdx(nmat, k)];

      // Right
      auto w_r = (Sr-vnr[0])/(Sr-Sm);
      if (solidx[k] > 0)
      {
        gnrStar = gnr[k];
        // gnrStar = gnlStar;
        for (std::size_t i=0; i<3; ++i)
          gnrStar[i][0] = w_r*gnr[k][i][0] + (
            gnr[k][i][1]*(vnr[1]-vnrStar[1]) + gnr[k][i][2]*(vnr[2]-vnrStar[2])
            )/(Sm-Sr);
        // rotate g back to original frame of reference
        grStar.push_back(tk::unrotateTensor(gnrStar, fn));
        // for (std::size_t i=0; i<3; ++i)
        //   for (std::size_t j=0; j<3; ++j)
        //     grStar[k][i][j] = w_r*gr[k][i][j];
      }
      uStar[1][volfracIdx(nmat, k)] = u[1][volfracIdx(nmat, k)];
      uStar[1][densityIdx(nmat, k)] =
        (Sr-vnr[0]) * u[1][densityIdx(nmat, k)] / (Sr-Sm);
      // uStar[1][energyIdx(nmat, k)] =
      //   ((Sr-vnr[0]) * u[1][energyIdx(nmat, k)] - apr[k]*vnr[0] + apStar[k]*Sm) / (Sr-Sm);
      uStar[1][energyIdx(nmat, k)] = u[1][energyIdx(nmat, k)] * w_r + (
        - (asignnr[k][0][0]-apr[k])*vnr[0]
        - asignnr[k][0][1]*vnr[1]
        - asignnr[k][0][2]*vnr[2]
        + (asignnStar[k][0][0]-apStar[k])*vnrStar[0]
        + asignnStar[k][0][1]*vnrStar[1]
        + asignnStar[k][0][2]*vnrStar[2]
        ) / (Sm-Sr);
      rhorStar += uStar[1][densityIdx(nmat, k)];

      // if (std::abs(glStar[k][0][1]-grStar[k][0][1]) > 1.0e-04) {
      //   printf("DBG\n");
      //   printf("glStar = %e, %e, %e\n", glStar[k][0][0], glStar[k][0][1], glStar[k][0][2]);
      //   printf("         %e, %e, %e\n", glStar[k][1][0], glStar[k][1][1], glStar[k][1][2]);
      //   printf("         %e, %e, %e\n", glStar[k][2][0], glStar[k][2][1], glStar[k][2][2]);
      //   printf("grStar = %e, %e, %e\n", grStar[k][0][0], grStar[k][0][1], grStar[k][0][2]);
      //   printf("         %e, %e, %e\n", grStar[k][1][0], grStar[k][1][1], grStar[k][1][2]);
      //   printf("         %e, %e, %e\n", grStar[k][2][0], grStar[k][2][1], grStar[k][2][2]);
      // }
    }

    // Numerical fluxes
    // -------------------------------------------------------------------------
    if (Sl > 0.0) {

      for (std::size_t idir=0; idir<3; ++idir)
        flx[momentumIdx(nmat, idir)] =
          u[0][momentumIdx(nmat, idir)] * vnl[0] + pl*fn[idir] - Tnl[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = u[0][volfracIdx(nmat, k)] * vnl[0];
        flx[densityIdx(nmat, k)] = u[0][densityIdx(nmat, k)] * vnl[0];
        flx[energyIdx(nmat, k)] = u[0][energyIdx(nmat, k)] * vnl[0];
        flx[energyIdx(nmat, k)] -= ul * (aTnl[k][0] - apl[k]*fn[0]);
        flx[energyIdx(nmat, k)] -= vl * (aTnl[k][1] - apl[k]*fn[1]);
        flx[energyIdx(nmat, k)] -= wl * (aTnl[k][2] - apl[k]*fn[2]);
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx[deformIdx(nmat,solidx[k],i,j)] = (
                gl[k][i][0] * ul +
                gl[k][i][1] * vl +
                gl[k][i][2] * wl ) * fn[j];
        }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back(apl[k]);
      // Store Riemann velocity
      flx.push_back(vnl[0]);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnl[k][i]-apl[k]*fn[i]);
        }
      }
    }

    else if (Sl <= 0.0 && Sm > 0.0) {

      for (std::size_t idir=0; idir<3; ++idir)
       flx[momentumIdx(nmat, idir)] =
         vlStar[idir] * rhorStar * Sm + pStar*fn[idir] - TnStar[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = uStar[0][volfracIdx(nmat, k)] * Sm;
        flx[densityIdx(nmat, k)] = uStar[0][densityIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] = uStar[0][energyIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] -= vlStar[0] * (aTnStar[k][0] - apStar[k]*fn[0]);
        flx[energyIdx(nmat, k)] -= vlStar[1] * (aTnStar[k][1] - apStar[k]*fn[1]);
        flx[energyIdx(nmat, k)] -= vlStar[2] * (aTnStar[k][2] - apStar[k]*fn[2]);
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx[deformIdx(nmat,solidx[k],i,j)] = (
                glStar[k][i][0] * vlStar[0] +
                glStar[k][i][1] * vlStar[1] +
                glStar[k][i][2] * vlStar[2] ) * fn[j];
        }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back(apStar[k]);
      // Store Riemann velocity
      flx.push_back(Sm);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnStar[k][i]-apStar[k]*fn[i]);
        }
      }
    }

    else if (Sm <= 0.0 && Sr >= 0.0) {

      for (std::size_t idir=0; idir<3; ++idir)
        flx[momentumIdx(nmat, idir)] =
          vrStar[idir] * rhorStar * Sm + pStar*fn[idir] - TnStar[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = uStar[1][volfracIdx(nmat, k)] * Sm;
        flx[densityIdx(nmat, k)] = uStar[1][densityIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] = uStar[1][energyIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] -= vrStar[0] * (aTnStar[k][0] - apStar[k]*fn[0]);
        flx[energyIdx(nmat, k)] -= vrStar[1] * (aTnStar[k][1] - apStar[k]*fn[1]);
        flx[energyIdx(nmat, k)] -= vrStar[2] * (aTnStar[k][2] - apStar[k]*fn[2]);
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx[deformIdx(nmat,solidx[k],i,j)] = (
                grStar[k][i][0] * vrStar[0] +
                grStar[k][i][1] * vrStar[1] +
                grStar[k][i][2] * vrStar[2] ) * fn[j];
        }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back(apStar[k]);
      // Store Riemann velocity
      flx.push_back(Sm);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnStar[k][i]-apStar[k]*fn[i]);
        }
      }
    }

    else {

      for (std::size_t idir=0; idir<3; ++idir)
        flx[momentumIdx(nmat, idir)] =
          u[1][momentumIdx(nmat, idir)] * vnr[0] + pr*fn[idir] - Tnr[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = u[1][volfracIdx(nmat, k)] * vnr[0];
        flx[densityIdx(nmat, k)] = u[1][densityIdx(nmat, k)] * vnr[0];
        flx[energyIdx(nmat, k)] = u[1][energyIdx(nmat, k)] * vnr[0];
        flx[energyIdx(nmat, k)] -= ur * (aTnr[k][0] - apr[k]*fn[0]);
        flx[energyIdx(nmat, k)] -= vr * (aTnr[k][1] - apr[k]*fn[1]);
        flx[energyIdx(nmat, k)] -= wr * (aTnr[k][2] - apr[k]*fn[2]);
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx[deformIdx(nmat,solidx[k],i,j)] = (
                gr[k][i][0] * ur +
                gr[k][i][1] * vr +
                gr[k][i][2] * wr ) * fn[j];
        }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back(apr[k]);
      // Store Riemann velocity
      flx.push_back(vnr[0]);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnr[k][i]-apr[k]*fn[i]);
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
  //  return ctr::FluxType::HLLCMultiMat; }
};

} // inciter::

#endif // HLLCMultiMat_h
