// *****************************************************************************
/*!
  \file      src/PDE/Riemann/HLL.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Harten-Lax-vanLeer's (HLL) Riemann flux function
  \details   This file implements Harten-Lax-vanLeer's (HLL) Riemann solver.
*/
// *****************************************************************************
#ifndef HLL_h
#define HLL_h

#include <vector>

#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EOS.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

//! HLL approximate Riemann solver
struct HLL {

  //! HLL approximate Riemann solver flux function
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann flux solution according to HLL, appended by Riemann
  //!   velocities and volume-fractions.
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
    std::vector< tk::real > flx(ncomp, 0), fl(ncomp, 0), fr(ncomp, 0);

    // Primitive quantities
    tk::real rhol(0.0), rhor(0.0);
    tk::real amatl(0.0), amatr(0.0), ac_l(0.0), ac_r(0.0);
    std::vector< tk::real > al_l(nmat, 0.0), al_r(nmat, 0.0),
                            pml(nmat, 0.0), pmr(nmat, 0.0);
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > g_l, g_r,
      asig_l, asig_r;
    std::vector< std::array< tk::real, 3 > > asign_l, asign_r;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gn_l, gn_r;
    std::array< tk::real, 3 > sign_l {{0, 0, 0}}, sign_r {{0, 0, 0}};
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

    // Independently limited velocities for advection
    auto ul = u[0][ncomp+velocityIdx(nmat, 0)];
    auto vl = u[0][ncomp+velocityIdx(nmat, 1)];
    auto wl = u[0][ncomp+velocityIdx(nmat, 2)];
    auto ur = u[1][ncomp+velocityIdx(nmat, 0)];
    auto vr = u[1][ncomp+velocityIdx(nmat, 1)];
    auto wr = u[1][ncomp+velocityIdx(nmat, 2)];

    // Face-normal velocities from advective velocities
    auto vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];
    auto vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];

    // Conservative flux functions
    // -------------------------------------------------------------------------
    for (std::size_t k=0; k<nmat; ++k)
    {
      // Left fluxes
      fl[volfracIdx(nmat, k)] = vnl * al_l[k];
      fl[densityIdx(nmat, k)] = vnl * u[0][densityIdx(nmat, k)];
      fl[energyIdx(nmat, k)] = vnl * u[0][energyIdx(nmat, k)];
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
      fr[volfracIdx(nmat, k)] = vnr * al_r[k];
      fr[densityIdx(nmat, k)] = vnr * u[1][densityIdx(nmat, k)];
      fr[energyIdx(nmat, k)] = vnr * u[1][energyIdx(nmat, k)];
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
    }

    // bulk momentum
    for (std::size_t idir=0; idir<3; ++idir)
    {
      fl[momentumIdx(nmat, idir)] = vnl*u[0][momentumIdx(nmat, idir)]
        - sign_l[idir];
      fr[momentumIdx(nmat, idir)] = vnr*u[1][momentumIdx(nmat, idir)]
        - sign_r[idir];
    }

    // Numerical fluxes
    // -------------------------------------------------------------------------

    // Signal velocities
    auto Sl = std::min((vnl-ac_l), (vnr-ac_r));
    auto Sr = std::max((vnl+ac_l), (vnr+ac_r));
    //// Signal velocities by Einfeldt (HLLE)
    //auto Sl = std::min( 0.0, std::min((vnl-ac_l), 0.5*((vnl+vnr)-(ac_l+ac_r))) );
    //auto Sr = std::max( 0.0, std::max((vnr+ac_r), 0.5*((vnl+vnr)+(ac_l+ac_r))) );

    // Numerical flux functions and wave-speeds
    auto c_plus(0.0), c_minus(0.0), p_plus(0.0), p_minus(0.0);
    if (Sl >= 0.0)
    {
      flx = fl;
      c_plus = vnl;
      p_plus = 1.0;
    }
    else if (Sr <= 0.0)
    {
      flx = fr;
      c_minus = vnr;
      p_minus = 1.0;
    }
    else
    {
      for (std::size_t k=0; k<flx.size(); ++k)
        flx[k] = (Sr*fl[k] - Sl*fr[k] + Sl*Sr*(u[1][k]-u[0][k])) / (Sr-Sl);
      c_plus = (Sr*vnl - Sr*Sl) / (Sr-Sl);
      c_minus = (Sr*Sl - Sl*vnr) / (Sr-Sl);
      p_plus = Sr / (Sr-Sl);
      p_minus = -Sl / (Sr-Sl);
    }

    // Quantities for non-conservative terms
    // -------------------------------------------------------------------------

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

    // Store Riemann g_ij (9*nsld)
    for (std::size_t k=0; k<nmat; ++k) {
      if (solidx[k] > 0) {
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            flx.push_back( p_plus*g_l[k][i][j] + p_minus*g_r[k][i][j] );
      }
    }

    Assert( flx.size() == (ncomp+nmat+1+3*nsld+9*nsld), "Size of "
            "multi-material flux vector incorrect" );

    return flx;
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::HLL; }

};

} // inciter::

#endif // HLL_h
