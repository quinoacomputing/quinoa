// *****************************************************************************
/*!
  \file      src/PDE/Riemann/HLLCMultiMatConstP.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2023 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     ConstP HLLC Riemann flux function for solids
  \details   This file implements the HLLC Riemann solver for constant P solids.
             Compared to the regular one the 3x3 tensor rotations use BLAS-free
             helpers (see TensorFast.hpp) rather than tk::(un)rotateTensor since
             those issue 2x cblas_dgemm calls for m=n=k=3 which is inefficient.
             Ref. Ndanou, S., Favrie, N., & Gavrilyuk, S. (2015). Multi-solid
             and multi-fluid diffuse interface model: Applications to dynamic
             fracture and fragmentation. Journal of Computational Physics, 295,
             523-555.
*/
// *****************************************************************************
#ifndef HLLCMultiMatConstP_h
#define HLLCMultiMatConstP_h

#include <vector>

#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EOS.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "MultiMat/MiscMultiMatFns.hpp"
#include "Riemann/TensorFast.hpp"
#include "Riemann/RiemannDevice.hpp"

namespace inciter {

//! ConstP HLLC approximate Riemann solver for solids
struct HLLCMultiMatConstP {

//! Caps for the fixed size scratch buffer matched to NMAT_MAX and NSTATE_MAX
//! in the KokkosDevice.hpp file. We expose at struct scope so callers can size fix
static constexpr std::size_t NMAT_MAX_FLUX=4;
static constexpr std::size_t NSTATE_MAX_FLUX=82;
static constexpr std::size_t NFLX_MAX=104;
// 51+4+1+12+36 = 104
// ncomp+nmat+1+3*nsld+9*nsld

//! HLLC approximate Riemann solver flux function
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \param[in] wn Mesh velocity normal to face
  //! \return Riemann solution according to Harten-Lax-van Leer-Contact
  //! \note The function signature must follow tk::RiemannFluxFn
  //! \note Templated on container types so one body serves both the host loop and the device kernel
  //! \note Because it's a template the KOKKOS_INLINE_FUNCTION will not force eager device codegen on hostside
  template< class MatBlkT, class FnT, class StateT, class SolidxT >
  KOKKOS_INLINE_FUNCTION static std::size_t
  flux( const MatBlkT& mat_blk,
        const FnT& fn,
        const StateT& u,
        const tk::real wn,
        std::size_t nmat,
        std::size_t nsld,
        std::size_t ncomp,
        std::size_t nstate,
        const SolidxT& solidx,
        Kokkos::Array< tk::real, NFLX_MAX >& flx )
  {
    using tk::Vec3Dev;
    using tk::Mat3Dev;

    // facenormal in device storage
    Vec3Dev fnd;
    for (std::size_t i=0; i<3; ++i)
      fnd[i] = fn[i];

    // flx is caller-owned so we have to zero consv block and append after it
    for (std::size_t i=0; i<ncomp; ++i) flx[i] = 0.0;
    std::size_t nflx = ncomp;

    // Primitive variables
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
    [[maybe_unused]] tk::real pl(0.0), pr(0.0);
    tk::real acl(0.0), acr(0.0);
    Kokkos::Array< tk::real, NMAT_MAX_FLUX > apl{}, apr{};
    Vec3Dev Tnl{}, Tnr{};
    Kokkos::Array< Vec3Dev, NMAT_MAX_FLUX > aTnl{}, aTnr{};
    Mat3Dev asigl, asigr;
    Mat3Dev signnl = tk::zeroMat3Dev();
    Mat3Dev signnr = tk::zeroMat3Dev();
    Kokkos::Array< Mat3Dev, NMAT_MAX_FLUX >
      gl{}, gr{}, gnl{}, gnr{}, asignnl{}, asignnr{};

    for (std::size_t k=0; k<nmat; ++k) {
      // Left state
      apl[k] = u[0][ncomp+pressureIdx(nmat, k)];
      pl += apl[k];

      // inv deformation gradient and Cauchy stress tensors
      gl[k] = getDeformGradDev(nmat, k, solidx, u[0]);
      asigl = getCauchyStressDev(nmat, k, ncomp, solidx, u[0]);
      for (std::size_t i=0; i<3; ++i) asigl[i][i] -= apl[k];

      // normal stress (traction) vector
      aTnl[k] = tk::matvecDev(asigl, fnd);
      for (std::size_t i=0; i<3; ++i)
        Tnl[i] += aTnl[k][i];

      // rotate stress vector
      asignnl[k] = tk::rotateTensorDev(asigl, fnd);
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          signnl[i][j] += asignnl[k][i][j];

      // rotate deformation gradient tensor for speed of sound in normal dir
      gnl[k] = tk::rotateTensorDev(gl[k], fnd);
      auto amatl = mat_blk[k].template compute< EOS::soundspeed >(
        u[0][densityIdx(nmat, k)], apl[k],
        u[0][volfracIdx(nmat, k)], k );

      // Right state
      apr[k] = u[1][ncomp+pressureIdx(nmat, k)];
      pr += apr[k];

      // inv deformation gradient and Cauchy stress tensors
      gr[k] = getDeformGradDev(nmat, k, solidx, u[1]);
      asigr = getCauchyStressDev(nmat, k, ncomp, solidx, u[1]);
      for (std::size_t i=0; i<3; ++i) asigr[i][i] -= apr[k];

      // normal stress (traction) vector
      aTnr[k] = tk::matvecDev(asigr, fnd);
      for (std::size_t i=0; i<3; ++i)
        Tnr[i] += aTnr[k][i];

      // rotate stress vector
      asignnr[k] = tk::rotateTensorDev(asigr, fnd);
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          signnr[i][j] += asignnr[k][i][j];

      // rotate deformation gradient tensor for speed of sound in normal dir
      gnr[k] = tk::rotateTensorDev(gr[k], fnd);
      auto amatr = mat_blk[k].template compute< EOS::soundspeed >(
        u[1][densityIdx(nmat, k)], apr[k],
        u[1][volfracIdx(nmat, k)], k );

      // Mixture speed of sound
      acl += u[0][densityIdx(nmat, k)] * amatl * amatl;
      acr += u[1][densityIdx(nmat, k)] * amatr * amatr;
    }
    acl = std::sqrt(acl/rhol);
    acr = std::sqrt(acr/rhor);

    // Rotated velocities from advective velocities
    Vec3Dev vvl, vvr;
    vvl[0] = ul; vvl[1] = vl; vvl[2] = wl;
    vvr[0] = ur; vvr[1] = vr; vvr[2] = wr;
    auto vnl = tk::rotateVectorDev(vvl, fnd);
    auto vnr = tk::rotateVectorDev(vvr, fnd);

    // ALE mesh motion relative normal velocity
    vnl[0] -= wn;
    vnr[0] -= wn;

    // Signal velocities
    auto Sl = fmin((vnl[0]-acl), (vnr[0]-acr));
    auto Sr = fmax((vnl[0]+acl), (vnr[0]+acr));
    auto Sm = ( rhor*vnr[0]*(Sr-vnr[0]) - rhol*vnl[0]*(Sl-vnl[0])
                - signnl[0][0] + signnr[0][0] )
              /( rhor*(Sr-vnr[0]) - rhol*(Sl-vnl[0]) );

    // Middle-zone (star) variables
    // the stress in the star zone is theoretically equivalent when derived
    // from the left or right state. However, there might be differences
    // numerically due to truncation etc. Hence two separate evaluations
    Kokkos::Array< Mat3Dev, NMAT_MAX_FLUX > asignnlStar{}, asignnrStar{};
    Mat3Dev signnlStar = tk::zeroMat3Dev();
    Mat3Dev signnrStar = tk::zeroMat3Dev();
    Mat3Dev siglStar = tk::zeroMat3Dev();
    Mat3Dev sigrStar = tk::zeroMat3Dev();
    Vec3Dev TnlStar{}, TnrStar{};
    Kokkos::Array< Vec3Dev, NMAT_MAX_FLUX > aTnlStar{}, aTnrStar{};
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
        asignnlStar[k][i][0] =
          ( (Sm-Sl)*u[0][densityIdx(nmat,k)]*asignnr[k][i][0]
          - (Sm-Sr)*u[1][densityIdx(nmat,k)]*asignnl[k][i][0]
          + (Sm-Sl)*u[0][densityIdx(nmat,k)]
          * (Sm-Sr)*u[1][densityIdx(nmat,k)]
          * (vnl[i]-vnr[i]) ) /
          ( (Sm-Sl)*u[0][densityIdx(nmat,k)]
          - (Sm-Sr)*u[1][densityIdx(nmat,k)]);
        asignnrStar[k][i][0] =
          ( (Sm-Sl)*u[0][densityIdx(nmat,k)]*asignnr[k][i][0]
          - (Sm-Sr)*u[1][densityIdx(nmat,k)]*asignnl[k][i][0]
          + (Sm-Sl)*u[0][densityIdx(nmat,k)]
          * (Sm-Sr)*u[1][densityIdx(nmat,k)]
          * (vnl[i]-vnr[i]) ) /
          ( (Sm-Sl)*u[0][densityIdx(nmat,k)]
          - (Sm-Sr)*u[1][densityIdx(nmat,k)]);
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

      auto asiglStar = tk::unrotateTensorDev(asignnlStar[k], fnd);
      auto asigrStar = tk::unrotateTensorDev(asignnrStar[k], fnd);
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
        {
          siglStar[i][j] += asiglStar[i][j];
          sigrStar[i][j] += asigrStar[i][j];
        }
      aTnlStar[k] = tk::matvecDev(asiglStar, fnd);
      aTnrStar[k] = tk::matvecDev(asigrStar, fnd);
      for (std::size_t i=0; i<3; ++i)
      {
        TnlStar[i] += aTnlStar[k][i];
        TnrStar[i] += aTnrStar[k][i];
      }
    }

    auto w_l = (vnl[0]-Sl)/(Sm-Sl);
    auto w_r = (vnr[0]-Sr)/(Sm-Sr);

    Vec3Dev vnlStar, vnrStar;
    // u*_L
    vnlStar[0] = Sm;
    vnlStar[1] = vnl[1]
      + (signnlStar[1][0] - signnl[1][0]) / (w_l*rhol*(vnl[0]-Sl));
    vnlStar[2] = vnl[2]
      + (signnlStar[2][0] - signnl[2][0]) / (w_l*rhol*(vnl[0]-Sl));
    // u*_R
    vnrStar[0] = Sm;
    vnrStar[1] = vnr[1]
      + (signnrStar[1][0] - signnr[1][0]) / (w_r*rhor*(vnr[0]-Sr));
    vnrStar[2] = vnr[2]
      + (signnrStar[2][0] - signnr[2][0]) / (w_r*rhor*(vnr[0]-Sr));

    auto vlStar = tk::unrotateVectorDev(vnlStar, fnd);
    auto vrStar = tk::unrotateVectorDev(vnrStar, fnd);

    const auto unl = vnl[0] + wn;
    const auto unr = vnr[0] + wn;
    const auto unStar = Sm + wn;

    Kokkos::Array< Kokkos::Array< tk::real, NSTATE_MAX_FLUX >, 2 > uStar;
    for (std::size_t s=0; s<2; ++s)
      for (std::size_t i=0; i<nstate; ++i) uStar[s][i] = u[s][i];

    tk::real rholStar(0.0), rhorStar(0.0);
    Kokkos::Array< Mat3Dev, NMAT_MAX_FLUX >
      gnlStar{}, gnrStar{}, glStar{}, grStar{};
    for (std::size_t k=0; k<nmat; ++k) {
      // Left
      gnlStar[k] = gnl[k];
      if (solidx[k] > 0)
      {
        gnlStar[k][0][0] = w_l * gnl[k][0][0]
          + gnl[k][0][1]*(vnl[1]-vnlStar[1])/(Sm-Sl)
          + gnl[k][0][2]*(vnl[2]-vnlStar[2])/(Sm-Sl);
        gnlStar[k][1][0] = w_l * gnl[k][1][0]
          + gnl[k][1][1]*(vnl[1]-vnlStar[1])/(Sm-Sl)
          + gnl[k][1][2]*(vnl[2]-vnlStar[2])/(Sm-Sl);
        gnlStar[k][2][0] = w_l * gnl[k][2][0]
          + gnl[k][2][1]*(vnl[1]-vnlStar[1])/(Sm-Sl)
          + gnl[k][2][2]*(vnl[2]-vnlStar[2])/(Sm-Sl);
      }
      // rotate g back to original frame of reference
      glStar[k] = tk::unrotateTensorDev(gnlStar[k], fnd);
      uStar[0][volfracIdx(nmat, k)] = u[0][volfracIdx(nmat, k)];
      uStar[0][densityIdx(nmat, k)] = w_l * u[0][densityIdx(nmat, k)];
      uStar[0][energyIdx(nmat, k)] = w_l * u[0][energyIdx(nmat, k)]
        + ( - asignnl[k][0][0]*unl
            - asignnl[k][1][0]*vnl[1]
            - asignnl[k][2][0]*vnl[2]
            + asignnlStar[k][0][0]*unStar
            + asignnlStar[k][1][0]*vnlStar[1]
            + asignnlStar[k][2][0]*vnlStar[2]
          ) / (Sm-Sl);
      rholStar += uStar[0][densityIdx(nmat, k)];

      // Right
      gnrStar[k] = gnr[k];
      if (solidx[k] > 0)
      {
        gnrStar[k][0][0] = w_r * gnr[k][0][0]
          + gnr[k][0][1]*(vnr[1]-vnrStar[1])/(Sm-Sr)
          + gnr[k][0][2]*(vnr[2]-vnrStar[2])/(Sm-Sr);
        gnrStar[k][1][0] = w_r * gnr[k][1][0]
          + gnr[k][1][1]*(vnr[1]-vnrStar[1])/(Sm-Sr)
          + gnr[k][1][2]*(vnr[2]-vnrStar[2])/(Sm-Sr);
        gnrStar[k][2][0] = w_r * gnr[k][2][0]
          + gnr[k][2][1]*(vnr[1]-vnrStar[1])/(Sm-Sr)
          + gnr[k][2][2]*(vnr[2]-vnrStar[2])/(Sm-Sr);
      }
      // rotate g back to original frame of reference
      grStar[k] = tk::unrotateTensorDev(gnrStar[k], fnd);
      uStar[1][volfracIdx(nmat, k)] = u[1][volfracIdx(nmat, k)];
      uStar[1][densityIdx(nmat, k)] = w_r * u[1][densityIdx(nmat, k)];
      uStar[1][energyIdx(nmat, k)] = w_r * u[1][energyIdx(nmat, k)]
        + ( - asignnr[k][0][0]*unr
            - asignnr[k][1][0]*vnr[1]
            - asignnr[k][2][0]*vnr[2]
            + asignnrStar[k][0][0]*unStar
            + asignnrStar[k][1][0]*vnrStar[1]
            + asignnrStar[k][2][0]*vnrStar[2]
          ) / (Sm-Sr);
      rhorStar += uStar[1][densityIdx(nmat, k)];
    }
    for (std::size_t idir=0; idir<3; ++idir) {
      uStar[0][momentumIdx(nmat, idir)] = w_l*u[0][momentumIdx(nmat, idir)]
        - (TnlStar[idir] - Tnl[idir])/(Sl-Sm);
      uStar[1][momentumIdx(nmat, idir)] = w_r*u[1][momentumIdx(nmat, idir)]
        - (TnrStar[idir] - Tnr[idir])/(Sr-Sm);
    }

    // Numerical fluxes
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
        }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k)
        flx[nflx++] = (std::sqrt((aTnl[k][0]*aTnl[k][0]
                                +aTnl[k][1]*aTnl[k][1]
                                +aTnl[k][2]*aTnl[k][2])));
      // Store Riemann velocity
      flx[nflx++] = ((vnl[0]+wn));
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx[nflx++] = (aTnl[k][i]);
        }
      }
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx[nflx++] = (gl[k][i][j]);
        }
      }

    }

    else if (Sl < 0.0 && 0.0 <= Sm) {
      
      Vec3Dev ulStarPhysical;
      for (std::size_t i=0; i<3; ++i)
        ulStarPhysical[i] = vlStar[i] + wn*fn[i];

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
        }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k)
        flx[nflx++] = (std::sqrt(aTnlStar[k][0]*aTnlStar[k][0]
                               +aTnlStar[k][1]*aTnlStar[k][1]
                               +aTnlStar[k][2]*aTnlStar[k][2]));
      // Store Riemann velocity
      flx[nflx++] = (Sm+wn);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx[nflx++] = (aTnlStar[k][i]);
        }
      }
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx[nflx++] = (glStar[k][i][j]);
        }
      }

    }

    else if (Sm < 0.0 && 0.0 <= Sr) {

      Vec3Dev urStarPhysical;
      for (std::size_t i=0; i<3; ++i)
        urStarPhysical[i] = vrStar[i] + wn*fn[i];

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
          }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k)
        flx[nflx++] = (std::sqrt(aTnrStar[k][0]*aTnrStar[k][0]
                               +aTnrStar[k][1]*aTnrStar[k][1]
                               +aTnrStar[k][2]*aTnrStar[k][2]));
      // Store Riemann velocity
      flx[nflx++] = (Sm+wn);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx[nflx++] = (aTnrStar[k][i]);
        }
      }
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx[nflx++] = (grStar[k][i][j]);
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
          }
      }

      // Quantities for non-conservative terms
      // Store Riemann-advected partial pressures
      for (std::size_t k=0; k<nmat; ++k)
        flx[nflx++] = (std::sqrt(aTnr[k][0]*aTnr[k][0]
                               +aTnr[k][1]*aTnr[k][1]
                               +aTnr[k][2]*aTnr[k][2]));
      // Store Riemann velocity
      flx[nflx++] = (vnr[0]+wn);
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx[nflx++] = (aTnr[k][i]);
        }
      }
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              flx[nflx++] = (gr[k][i][j]);
        }
      }

    }

    return nflx;
  }

  ////! Flux type accessor
  ////! \return Flux type
  //static ctr::FluxType type() noexcept {
  //  return ctr::FluxType::HLLCMultiMat; }
};

} // inciter::

#endif // HLLCMultiMatConstP_h
