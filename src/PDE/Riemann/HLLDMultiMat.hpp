// *****************************************************************************
/*!
  \file      src/PDE/Riemann/HLLDMultiMat.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2023 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     HLLD Riemann flux function for solids
  \details   This file implements the HLLD Riemann solver for solids.
             Ref. Barton, Philip T. "An interface-capturing Godunov method for
             the simulation of compressible solid-fluid problems." Journal of
             Computational Physics 390 (2019): 25-50.
*/
// *****************************************************************************
#ifndef HLLDMultiMat_h
#define HLLDMultiMat_h

#include <vector>

#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EOS.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "MultiMat/MiscMultiMatFns.hpp"

namespace inciter {

//! HLLD approximate Riemann solver for solids
struct HLLDMultiMat {

//! HLLD approximate Riemann solver flux function
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
      for (std::size_t i=0; i<3; ++i) asigl[i][i] -= apl[k];

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
      for (std::size_t i=0; i<3; ++i) asigr[i][i] -= apr[k];

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
                - signnl[0][0] + signnr[0][0] )
              /( rhor*(Sr-vnr[0]) - rhol*(Sl-vnl[0]) );

    // Next zone inwards (star) variables
    // -------------------------------------------------------------------------
    tk::real acsl(0.0), acsr(0.0);
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > asignnStar;
    std::array< std::array< tk::real, 3 >, 3 >
      signnStar{{{0,0,0},{0,0,0},{0,0,0}}};
    std::array< std::array< tk::real, 3 >, 3 > sigStar{{{0,0,0},{0,0,0},{0,0,0}}};
    asignnStar.resize(nmat);
    std::array< tk::real, 3 > TnStar{{0, 0, 0}};
    std::vector< std::array< tk::real, 3 > > aTnStar;
    for (std::size_t k=0; k<nmat; ++k) {
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          asignnStar[k][i][j] = asignnl[k][i][j];

      for (std::size_t i=0; i<1; ++i)
        asignnStar[k][i][i] =
          ( (vnr[0]-Sr)*u[1][densityIdx(nmat,k)]*asignnl[k][i][i]
          - (vnl[0]-Sl)*u[0][densityIdx(nmat,k)]*asignnr[k][i][i]
          + (vnl[0]-Sl)*u[0][densityIdx(nmat,k)]
          * (vnr[0]-Sr)*u[1][densityIdx(nmat,k)]
          * (vnr[0]-vnl[0]) ) /
          ( (vnr[0]-Sr)*u[1][densityIdx(nmat,k)]
          - (vnl[0]-Sl)*u[0][densityIdx(nmat,k)]);
 
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

    std::array< tk::real, 3 > vnlStar, vnrStar;
    // u*_L
    vnlStar[0] = Sm;
    vnlStar[1] = vnl[1];
    vnlStar[2] = vnl[2];
    // u*_R
    vnrStar[0] = Sm;
    vnrStar[1] = vnr[1];
    vnrStar[2] = vnr[2];

    auto vlStar = tk::unrotateVector(vnlStar, fn);
    auto vrStar = tk::unrotateVector(vnrStar, fn);

    auto uStar = u;

    tk::real rholStar(0.0), rhorStar(0.0);
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gnlStar, gnrStar;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > glStar, grStar;
    for (std::size_t k=0; k<nmat; ++k) {
      // Left
      auto w_l = (vnl[0]-Sl)/(Sm-Sl);
      if (solidx[k] > 0)
      {
        gnlStar.push_back(gnl[k]);
        for (std::size_t i=1; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            gnlStar[k][i][j] *= w_l;
        // rotate g back to original frame of reference
        glStar.push_back(tk::unrotateTensor(gnlStar[k], fn));
      }
      uStar[0][volfracIdx(nmat, k)] = u[0][volfracIdx(nmat, k)];
      uStar[0][densityIdx(nmat, k)] = w_l * u[0][densityIdx(nmat, k)];
      uStar[0][energyIdx(nmat, k)] = w_l * u[0][energyIdx(nmat, k)]
        + (Sm-vnl[0]) * (u[0][densityIdx(nmat, k)]*Sm - asignnStar[k][0][0]/(Sl-vnl[0]));
      // uStar[0][energyIdx(nmat, k)] = u[0][energyIdx(nmat, k)]
      //   + (Sm-vnl[0])*(Sm-asignnl[k][0][0]/(u[0][densityIdx(nmat, k)]*(Sl-vnl[0])));
      rholStar += uStar[0][densityIdx(nmat, k)];

      auto amatl = mat_blk[k].compute< EOS::shearspeed >(
        uStar[0][densityIdx(nmat, k)], uStar[0][volfracIdx(nmat, k)], k );

      // Right
      auto w_r = (vnr[0]-Sr)/(Sm-Sr);
      if (solidx[k] > 0)
      {
        gnrStar.push_back(gnr[k]);
        for (std::size_t i=1; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            gnrStar[k][i][j] *= w_r;
        // rotate g back to original frame of reference
        grStar.push_back(tk::unrotateTensor(gnrStar[k], fn));
      }
      uStar[1][volfracIdx(nmat, k)] = u[1][volfracIdx(nmat, k)];
      uStar[1][densityIdx(nmat, k)] = w_r * u[1][densityIdx(nmat, k)];
      uStar[1][energyIdx(nmat, k)] = w_r * u[1][energyIdx(nmat, k)]
        + (Sm-vnr[0]) * (u[1][densityIdx(nmat, k)]*Sm - asignnStar[k][0][0]/(Sr-vnr[0]));
      // uStar[1][energyIdx(nmat, k)] = u[1][energyIdx(nmat, k)]
      //   + (Sm-vnr[0])*(Sm-asignnr[k][0][0]/(u[1][densityIdx(nmat, k)]*(Sr-vnr[0])));
      rhorStar += uStar[1][densityIdx(nmat, k)];

      auto amatr = mat_blk[k].compute< EOS::shearspeed >(
        uStar[1][densityIdx(nmat, k)], uStar[1][volfracIdx(nmat, k)], k );

      // Mixture speed of sound
      acsl += uStar[0][densityIdx(nmat, k)] * amatl * amatl;
      acsr += uStar[1][densityIdx(nmat, k)] * amatr * amatr;
    }
    acsl = std::sqrt(acsl/rholStar);
    acsr = std::sqrt(acsr/rhorStar);

    // Signal velocities
    auto Ssl = vnlStar[0]-acsl;
    auto Ssr = vnrStar[0]+acsr;

    // Middle-zone (StarStar) variables. Only relevant if solids are present.
    // -------------------------------------------------------------------------
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > asignnStarStar;
    std::array< std::array< tk::real, 3 >, 3 >
      signnStarStar{{{0,0,0},{0,0,0},{0,0,0}}};
    std::array< std::array< tk::real, 3 >, 3 >
      sigStarStar{{{0,0,0},{0,0,0},{0,0,0}}};
    asignnStarStar.resize(nmat);
    std::array< tk::real, 3 > TnStarStar{{0, 0, 0}};
    std::vector< std::array< tk::real, 3 > > aTnStarStar;
    std::array< tk::real, 3 > vnlStarStar, vnrStarStar;
    std::array< tk::real, 3 > vlStarStar, vrStarStar;
    tk::real rholStarStar(0.0), rhorStarStar(0.0);
    std::array< std::array< tk::real, 3 >, 3 > gnlStarStar{{{0,0,0},{0,0,0},{0,0,0}}};
    std::array< std::array< tk::real, 3 >, 3 > gnrStarStar{{{0,0,0},{0,0,0},{0,0,0}}};
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > glStarStar, grStarStar;
    auto uStarStar = uStar;
    if (nsld > 0) {
      for (std::size_t k=0; k<nmat; ++k) {
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            asignnStarStar[k][i][j] = asignnStar[k][i][j];

        if (solidx[k] > 0)
        {
          for (std::size_t i=1; i<3; ++i)
            asignnStarStar[k][i][0] =
              ( (Sm-Ssl)*uStar[0][densityIdx(nmat,k)]*asignnr[k][i][0]
              - (Sm-Ssr)*uStar[1][densityIdx(nmat,k)]*asignnl[k][i][0]
              + (Sm-Ssl)*uStar[0][densityIdx(nmat,k)]
              * (Sm-Ssr)*uStar[1][densityIdx(nmat,k)]
              * (vnl[i]-vnr[i]) ) /
              ( (Sm-Ssl)*uStar[0][densityIdx(nmat,k)]
              - (Sm-Ssr)*uStar[1][densityIdx(nmat,k)]);
          // Symmetry
          asignnStarStar[k][0][1] = asignnStarStar[k][1][0];
          asignnStarStar[k][0][2] = asignnStarStar[k][2][0];
        }

        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            signnStarStar[i][j] += asignnStarStar[k][i][j];

        auto asigStarStar = tk::unrotateTensor(asignnStarStar[k], fn);
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            sigStarStar[i][j] += asigStarStar[i][j];
        aTnStarStar.push_back(tk::matvec(asigStarStar, fn));
        for (std::size_t i=0; i<3; ++i)
          TnStarStar[i] += aTnStarStar[k][i];
      }

      // u*_L
      vnlStarStar[0] = Sm;
      vnlStarStar[1] = vnlStar[1]
       + (signnStarStar[1][0] - signnl[1][0]) / (rholStar*(Sm-Ssl));
      vnlStarStar[2] = vnlStar[2]
       + (signnStarStar[2][0] - signnl[2][0]) / (rholStar*(Sm-Ssl));
      // vnlStarStar[1] = (signnl[1][0] - signnr[1][0]
      //                  + rholStar*(Ssl-Sm)*vnlStar[1]
      //                  - rhorStar*(Ssr-Sm)*vnrStar[1])
      //                  / (rholStar*(Ssl-Sm) - rhorStar*(Ssr-Sm));
      // vnlStarStar[2] = (signnl[2][0] - signnr[2][0]
      //                  + rholStar*(Ssl-Sm)*vnlStar[2]
      //                  - rhorStar*(Ssr-Sm)*vnrStar[2])
      //                  / (rholStar*(Ssl-Sm) - rhorStar*(Ssr-Sm));
      // u*_R
      vnrStarStar[0] = Sm;
      vnrStarStar[1] = vnrStar[1]
       + (signnStarStar[1][0] - signnr[1][0]) / (rhorStar*(Sm-Ssr));
      vnrStarStar[2] = vnrStar[2]
       + (signnStarStar[2][0] - signnr[2][0]) / (rhorStar*(Sm-Ssr));
      // vnrStarStar[1] = (signnl[1][0] - signnr[1][0]
      //                  + rholStar*(Ssl-Sm)*vnlStar[1]
      //                  - rhorStar*(Ssr-Sm)*vnrStar[1])
      //                  / (rholStar*(Ssl-Sm) - rhorStar*(Ssr-Sm));
      // vnrStarStar[2] = (signnl[2][0] - signnr[2][0]
      //                  + rholStar*(Ssl-Sm)*vnlStar[2]
      //                  - rhorStar*(Ssr-Sm)*vnrStar[2])
      //                  / (rholStar*(Ssl-Sm) - rhorStar*(Ssr-Sm));

      vlStarStar = tk::unrotateVector(vnlStarStar, fn);
      vrStarStar = tk::unrotateVector(vnrStarStar, fn);

      for (std::size_t k=0; k<nmat; ++k) {
        // Left
        if (solidx[k] > 0)
        {
          gnlStarStar = gnlStar[k];
          for (std::size_t i=1; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              gnlStarStar[i][j] +=
                gnlStar[k][j][j] * (vnlStarStar[i]-vnl[i])/(Sm-Ssl);
          // rotate g back to original frame of reference
          glStarStar.push_back(tk::unrotateTensor(gnlStarStar, fn));
        }
        uStarStar[0][volfracIdx(nmat, k)] = uStar[0][volfracIdx(nmat, k)];
        uStarStar[0][densityIdx(nmat, k)] = uStar[0][densityIdx(nmat, k)];
        uStarStar[0][energyIdx(nmat, k)] = uStar[0][energyIdx(nmat, k)]
          + ( - asignnl[k][1][0]*vnl[1]
              - asignnl[k][2][0]*vnl[2]
              + asignnStarStar[k][1][0]*vnlStarStar[1]
              + asignnStarStar[k][2][0]*vnlStarStar[2]
            ) / (Sm-Ssl);
        rholStarStar += uStarStar[0][densityIdx(nmat, k)];

        // Right
        if (solidx[k] > 0)
        {
          gnrStarStar = gnrStar[k];
          for (std::size_t i=1; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              gnrStarStar[i][j] +=
                gnrStar[k][j][j] * (vnrStarStar[i]-vnr[i])/(Sm-Ssr);
          // rotate g back to original frame of reference
          grStarStar.push_back(tk::unrotateTensor(gnrStarStar, fn));
        }
        uStarStar[1][volfracIdx(nmat, k)] = uStar[1][volfracIdx(nmat, k)];
        uStarStar[1][densityIdx(nmat, k)] = uStar[1][densityIdx(nmat, k)];
        uStarStar[1][energyIdx(nmat, k)] = uStar[1][energyIdx(nmat, k)]
          + ( - asignnr[k][1][0]*vnr[1]
              - asignnr[k][2][0]*vnr[2]
              + asignnStarStar[k][1][0]*vnrStarStar[1]
              + asignnStarStar[k][2][0]*vnrStarStar[2]
            ) / (Sm-Ssr);
        rhorStarStar += uStarStar[1][densityIdx(nmat, k)];
      }
    }
    
    // Numerical fluxes
    // -------------------------------------------------------------------------
    if (Sl > 0.0) {

      for (std::size_t idir=0; idir<3; ++idir)
        flx[momentumIdx(nmat, idir)] =
          u[0][momentumIdx(nmat, idir)] * vnl[0] - Tnl[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = u[0][volfracIdx(nmat, k)] * vnl[0];
        flx[densityIdx(nmat, k)] = u[0][densityIdx(nmat, k)] * vnl[0];
        flx[energyIdx(nmat, k)] = u[0][energyIdx(nmat, k)] * vnl[0];
        flx[energyIdx(nmat, k)] -= ul * aTnl[k][0];
        flx[energyIdx(nmat, k)] -= vl * aTnl[k][1];
        flx[energyIdx(nmat, k)] -= wl * aTnl[k][2];
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
        flx.push_back(std::sqrt((aTnl[k][0]*aTnl[k][0]
                                +aTnl[k][1]*aTnl[k][1]
                                +aTnl[k][2]*aTnl[k][2])));
      // Store Riemann velocity
      flx.push_back(vnl[0]);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnl[k][i]);
        }
      }

    }

    else if (Sl <= 0.0 && Ssl > 0.0) {

      for (std::size_t idir=0; idir<3; ++idir)
       flx[momentumIdx(nmat, idir)] =
         vlStar[idir] * rholStar * Sm - TnStar[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = uStar[0][volfracIdx(nmat, k)] * Sm;
        flx[densityIdx(nmat, k)] = uStar[0][densityIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] = uStar[0][energyIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] -= vlStar[0] * aTnStar[k][0];
        flx[energyIdx(nmat, k)] -= vlStar[1] * aTnStar[k][1];
        flx[energyIdx(nmat, k)] -= vlStar[2] * aTnStar[k][2];
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
        flx.push_back(std::sqrt(aTnStar[k][0]*aTnStar[k][0]
                               +aTnStar[k][1]*aTnStar[k][1]
                               +aTnStar[k][2]*aTnStar[k][2]));
      // Store Riemann velocity
      flx.push_back(Sm); //(vnl[0] + Sl*((Sl-vnl[0])/(Sl-Sm)-1.0)); //(Sm);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnStar[k][i]);
        }
      }

    }

    else if (Ssl <= 0.0 && Sm > 0.0) {

      for (std::size_t idir=0; idir<3; ++idir)
       flx[momentumIdx(nmat, idir)] =
         vlStarStar[idir] * rholStarStar * Sm- TnStarStar[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = uStarStar[0][volfracIdx(nmat, k)] * Sm;
        flx[densityIdx(nmat, k)] = uStarStar[0][densityIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] = uStarStar[0][energyIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] -= vlStarStar[0] * aTnStarStar[k][0];
        flx[energyIdx(nmat, k)] -= vlStarStar[1] * aTnStarStar[k][1];
        flx[energyIdx(nmat, k)] -= vlStarStar[2] * aTnStarStar[k][2];
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
              for (std::size_t j=0; j<3; ++j)
                flx[deformIdx(nmat,solidx[k],i,j)] = (
                  glStarStar[k][i][0] * vlStarStar[0] +
                  glStarStar[k][i][1] * vlStarStar[1] +
                  glStarStar[k][i][2] * vlStarStar[2] ) * fn[j];
          }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back(std::sqrt(aTnStarStar[k][0]*aTnStarStar[k][0]
                               +aTnStarStar[k][1]*aTnStarStar[k][1]
                               +aTnStarStar[k][2]*aTnStarStar[k][2]));
      // Store Riemann velocity
      flx.push_back(Sm); //(vnl[0] + Sl*((Sl-vnl[0])/(Sl-Sm)-1.0)); //(Sm);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnStarStar[k][i]);
        }
      }

    }

    else if (Sm <= 0.0 && Ssr > 0.0) {

      for (std::size_t idir=0; idir<3; ++idir)
        flx[momentumIdx(nmat, idir)] =
          vrStarStar[idir] * rhorStarStar * Sm - TnStarStar[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = uStarStar[1][volfracIdx(nmat, k)] * Sm;
        flx[densityIdx(nmat, k)] = uStarStar[1][densityIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] = uStarStar[1][energyIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] -= vrStarStar[0] * aTnStarStar[k][0];
        flx[energyIdx(nmat, k)] -= vrStarStar[1] * aTnStarStar[k][1];
        flx[energyIdx(nmat, k)] -= vrStarStar[2] * aTnStarStar[k][2];
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
              for (std::size_t j=0; j<3; ++j)
                flx[deformIdx(nmat,solidx[k],i,j)] = (
                  grStarStar[k][i][0] * vrStarStar[0] +
                  grStarStar[k][i][1] * vrStarStar[1] +
                  grStarStar[k][i][2] * vrStarStar[2] ) * fn[j];
          }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back(std::sqrt(aTnStarStar[k][0]*aTnStarStar[k][0]
                               +aTnStarStar[k][1]*aTnStarStar[k][1]
                               +aTnStarStar[k][2]*aTnStarStar[k][2]));
      // Store Riemann velocity
      flx.push_back(Sm); //vnr[0] + Sr*((Sr-vnr[0])/(Sr-Sm)-1.0)); //(Sm);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnStarStar[k][i]);
        }
      }

    }

    else if (Ssr <= 0.0 && Sr > 0.0) {

      for (std::size_t idir=0; idir<3; ++idir)
        flx[momentumIdx(nmat, idir)] =
          vrStar[idir] * rhorStar * Sm - TnStar[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = uStar[1][volfracIdx(nmat, k)] * Sm;
        flx[densityIdx(nmat, k)] = uStar[1][densityIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] = uStar[1][energyIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] -= vrStar[0] * aTnStar[k][0];
        flx[energyIdx(nmat, k)] -= vrStar[1] * aTnStar[k][1];
        flx[energyIdx(nmat, k)] -= vrStar[2] * aTnStar[k][2];
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
        flx.push_back(std::sqrt(aTnStar[k][0]*aTnStar[k][0]
                               +aTnStar[k][1]*aTnStar[k][1]
                               +aTnStar[k][2]*aTnStar[k][2]));
      // Store Riemann velocity
      flx.push_back(Sm); //(vnr[0] + Sr*((Sr-vnr[0])/(Sr-Sm)-1.0)); //(Sm);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnStar[k][i]);
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
        flx[energyIdx(nmat, k)] = u[1][energyIdx(nmat, k)] * vnr[0];
        flx[energyIdx(nmat, k)] -= ur * aTnr[k][0];
        flx[energyIdx(nmat, k)] -= vr * aTnr[k][1];
        flx[energyIdx(nmat, k)] -= wr * aTnr[k][2];
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
        flx.push_back(std::sqrt(aTnr[k][0]*aTnr[k][0]
                               +aTnr[k][1]*aTnr[k][1]
                               +aTnr[k][2]*aTnr[k][2]));
      // Store Riemann velocity
      flx.push_back(vnr[0]);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTnr[k][i]);
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
  //  return ctr::FluxType::HLLDMultiMat; }
};

} // inciter::

#endif // HLLDMultiMat_h
