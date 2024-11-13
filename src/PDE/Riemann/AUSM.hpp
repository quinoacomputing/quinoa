// *****************************************************************************
/*!
  \file      src/PDE/Riemann/AUSM.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Advection Upstream Splitting Method (AUSM+) Riemann flux function
  \details   This file implements the Advection Upstream Splitting Method (AUSM)
             Riemann solver, with the all-speed corrections.
             Ref. Liou, M. S. (2006). A sequel to AUSM, Part II: AUSM+-up for
             all speeds. Journal of computational physics, 214(1), 137-170.
*/
// *****************************************************************************
#ifndef AUSM_h
#define AUSM_h

#include <vector>

#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "SplitMachFns.hpp"
#include "EoS/EOS.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

//! AUSM+up approximate Riemann solver
struct AUSM {

  //! AUSM+up approximate Riemann solver flux function
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann flux solution according to AUSM+up, appended by Riemann
  //!   velocities and volume-fractions.
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::vector< EOS >& mat_blk,
        const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::vector< std::array< tk::real, 3 > >& = {} )
  {
    auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
    auto k_p = g_inputdeck.get< tag::lowspeed_kp >();
    const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

    auto nsld = numSolids(nmat, solidx);
    auto ncomp = u[0].size()-(3+nmat+nsld*6);
    std::vector< tk::real > flx( ncomp, 0 );

    // Primitive variables
    tk::real rhol(0.0), rhor(0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      rhol += u[0][densityIdx(nmat, k)];
      rhor += u[1][densityIdx(nmat, k)];
    }

    tk::real pl(0.0), pr(0.0), amatl(0.0), amatr(0.0), asmatl(0.0), asmatr(0.0);
    std::vector< tk::real > al_l(nmat, 0.0), al_r(nmat, 0.0),
                            pml(nmat, 0.0), pmr(nmat, 0.0),
                            arhom12(nmat, 0.0),
                            asecmat12(nmat, 0.0),
                            amat12(nmat, 0.0);
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > g_l, g_r,
      asig_l, asig_r, asigrot_l, asigrot_r;
    std::vector< std::array< tk::real, 3 > > aTn_l, aTn_r;
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > gn_l, gn_r;
    std::array< tk::real, 3 > Tn_l {{0, 0, 0}}, Tn_r {{0, 0, 0}}, Trot_l{{0, 0, 0}}, Trot_r{{0, 0, 0}};
    for (std::size_t k=0; k<nmat; ++k)
    {
      // Left state
      // -----------------------------------------------------------------------
      al_l[k] = u[0][volfracIdx(nmat, k)];
      pml[k] = u[0][ncomp+pressureIdx(nmat, k)];
      pl += pml[k];

      // inv deformation gradient and Cauchy stress tensors
      g_l.push_back(getDeformGrad(nmat, k, u[0]));
      asig_l.push_back(getCauchyStress(nmat, k, ncomp, u[0]));
      for (std::size_t i=0; i<3; ++i) asig_l[k][i][i] -= pml[k];

      // normal stress (traction) vector
      aTn_l.push_back(tk::matvec(asig_l[k], fn));
      for (std::size_t i=0; i<3; ++i)
        Tn_l[i] += aTn_l[k][i];

      // rotate stress vector
      asigrot_l.push_back(tk::rotateTensor(asig_l[k], fn));
      for (std::size_t i=0; i<3; ++i)
        Trot_l[i] += asigrot_l[k][0][i];

      // rotate deformation gradient tensor for speed of sound in normal dir
      gn_l.push_back(tk::rotateTensor(g_l[k], fn));
      amatl = mat_blk[k].compute< EOS::soundspeed >(
        u[0][densityIdx(nmat, k)], pml[k], al_l[k], k, gn_l[k] );
      asmatl = mat_blk[k].compute< EOS::shearspeed >(
        u[0][densityIdx(nmat, k)], al_l[k], k );

      // Right state
      // -----------------------------------------------------------------------
      al_r[k] = u[1][volfracIdx(nmat, k)];
      pmr[k] = u[1][ncomp+pressureIdx(nmat, k)];
      pr += pmr[k];

      // inv deformation gradient and Cauchy stress tensors
      g_r.push_back(getDeformGrad(nmat, k, u[1]));
      asig_r.push_back(getCauchyStress(nmat, k, ncomp, u[1]));
      for (std::size_t i=0; i<3; ++i) asig_r[k][i][i] -= pmr[k];

      // normal stress (traction) vector
      aTn_r.push_back(tk::matvec(asig_r[k], fn));
      for (std::size_t i=0; i<3; ++i)
        Tn_r[i] += aTn_r[k][i];

      // rotate stress vector
      asigrot_r.push_back(tk::rotateTensor(asig_r[k], fn));
      for (std::size_t i=0; i<3; ++i)
        Trot_r[i] += asigrot_r[k][0][i];

      // rotate deformation gradient tensor for speed of sound in normal dir
      gn_r.push_back(tk::rotateTensor(g_r[k], fn));
      amatr = mat_blk[k].compute< EOS::soundspeed >(
        u[1][densityIdx(nmat, k)], pmr[k], al_r[k], k, gn_r[k] );
      asmatr = mat_blk[k].compute< EOS::shearspeed >(
        u[1][densityIdx(nmat, k)], al_r[k], k );

      // Average states for mixture speed of sound
      arhom12[k] = 0.5*(u[0][densityIdx(nmat, k)] + u[1][densityIdx(nmat, k)]);
      amat12[k] = 0.5*(amatl+amatr);
      asecmat12[k] = 0.5*(asmatl+asmatr);
    }

    auto rho12 = 0.5*(rhol+rhor);

    // mixture speed of sound
    tk::real ac12(0.0), acsec12(0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      ac12 += (arhom12[k]*amat12[k]*amat12[k]);
      acsec12 += (arhom12[k]*asecmat12[k]*asecmat12[k]);
    }
    ac12 = std::sqrt( ac12/rho12 );
    acsec12 = std::max(1e-8, std::sqrt( acsec12/rho12 ));

    std::array< tk::real, 3 > acboth{{ac12, acsec12, acsec12}};

    // Independently limited velocities for advection
    auto ul = u[0][ncomp+velocityIdx(nmat, 0)];
    auto vl = u[0][ncomp+velocityIdx(nmat, 1)];
    auto wl = u[0][ncomp+velocityIdx(nmat, 2)];
    auto ur = u[1][ncomp+velocityIdx(nmat, 0)];
    auto vr = u[1][ncomp+velocityIdx(nmat, 1)];
    auto wr = u[1][ncomp+velocityIdx(nmat, 2)];

    // Coordinate transform velocities to r,s,t coordinate system where the
    // first coordinate (r) aligns with face-normal (fn)
    auto utl = tk::rotateVector({ul, vl, wl}, fn);
    auto utr = tk::rotateVector({ur, vr, wr}, fn);

    // Mach numbers
    std::array< tk::real, 3 > ml, mr;
    for (std::size_t i=0; i<3; ++i) {
      ml[i] = utl[i]/acboth[i];
      mr[i] = utr[i]/acboth[i];
    }

    // All-speed parameters
    // These parameters control the amount of all-speed diffusion necessary for
    // low-Mach flows. Setting k_u and k_p to zero does not add any all-speed
    // diffusion, whereas setting k_u and k_p to 1 adds maximum recommended
    // all-speed diffusion. See "Liou, M. S. (2006). A sequel to AUSM, Part II:
    // AUSM+-up for all speeds. Journal of computational physics, 214(1),
    // 137-170" for more mathematical explanation. k_u is the velocity diffusion
    // term and k_p is the pressure diffusion term. These two terms reduce
    // pressure-velocity decoupling (chequerboarding/odd-even oscillations).
    tk::real k_u(1.0), f_a(1.0);

    // Split Mach polynomials
    std::array< std::array< tk::real, 4 >, 3 > msl, msr;
    for (std::size_t i=0; i<3; ++i) {
      msl[i] = splitmach_ausm( f_a, ml[i] );
      msr[i] = splitmach_ausm( f_a, mr[i] );
    }

    // Riemann Mach number
    auto m0 = 1.0 - (0.5*(utl[0]*utl[0] + utr[0]*utr[0])/(ac12*ac12));
    auto mp = -k_p* std::max(m0, 0.0) * (pr-pl) / (f_a*rho12*ac12*ac12);
    auto m12 = msl[0][0] + msr[0][1] + mp;
    auto vriem = ac12 * m12;

    // Pressure correction
    auto pu = -k_u* msl[0][2] * msr[0][3] * f_a * rho12 * ac12 * (utr[0]-utl[0]);

    // Flux vector splitting
    auto l_plus = 0.5 * (vriem + std::fabs(vriem));
    auto l_minus = 0.5 * (vriem - std::fabs(vriem));

    // 3D flux vector splitting
    std::array< tk::real, 3 > mt12, lt_plus, lt_minus;
    for (std::size_t i=0; i<3; ++i) {
      mt12[i] = msl[i][0] + msr[i][1] + mp;
      auto vwave = acboth[i] * mt12[i];
      lt_plus[i] = 0.5 * (vwave + std::fabs(vwave));
      lt_minus[i] = 0.5 * (vwave - std::fabs(vwave));
    }

    // Conservative fluxes
    for (std::size_t k=0; k<nmat; ++k)
    {
      flx[volfracIdx(nmat, k)] = l_plus*al_l[k] + l_minus*al_r[k];
      flx[densityIdx(nmat, k)] = l_plus*u[0][densityIdx(nmat, k)]
                              + l_minus*u[1][densityIdx(nmat, k)];
      flx[energyIdx(nmat, k)] = l_plus*u[0][energyIdx(nmat, k)]
        + l_minus*u[1][energyIdx(nmat, k)];
      // stress-work term in energy fluxes
      //for (std::size_t i=0; i<3; ++i) {
      //  flx[energyIdx(nmat, k)] -= (
      //    lt_plus[i]*asigrot_l[k][i][0] + lt_minus[i]*asigrot_r[k][i][0] );
      //}
      flx[energyIdx(nmat, k)] -= (
        lt_plus[0]*asigrot_l[k][0][0] + lt_minus[0]*asigrot_r[k][0][0] +
        msl[1][2]*asigrot_l[k][1][0]*utl[1] + msr[1][3]*asigrot_r[k][1][0]*utr[1] +
        msl[2][2]*asigrot_l[k][2][0]*utl[2] + msr[2][3]*asigrot_r[k][2][0]*utr[2]
        );

      // inv deformation gradient tensor fluxes
      if (solidx[k] > 0) {

        //// 1. rotate entire flux vector (with the partial deriv)
        //// -------------------------------------------------------------------
        //std::array< std::array< tk::real, 3 >, 3 > flxn{{
        //  {{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}} }};

        ////// (a) Mach number splitting
        ////for (std::size_t i=0; i<3; ++i) {
        ////  // get fluxes for g_k in rotated frame
        ////  flxn[i][0] =
        ////      lt_plus[0] * gn_l[k][i][0] + lt_minus[0] * gn_r[k][i][0]
        ////    + lt_plus[1] * gn_l[k][i][1] + lt_minus[1] * gn_r[k][i][1]
        ////    + lt_plus[2] * gn_l[k][i][2] + lt_minus[2] * gn_r[k][i][2];
        ////  // flux second component (j) is always zero for j /= 0, since rotated
        ////  // frame is aligned to face-normal, i.e. fnt = (1, 0, 0)
        ////}

        //// (b) Mach number splitting in normal dir; pressure split in tangential
        //for (std::size_t i=0; i<3; ++i) {
        //  // get fluxes for g_k in rotated frame
        //  flxn[i][0] =
        //    lt_plus[0] * gn_l[k][i][0] + lt_minus[0] * gn_r[k][i][0] +
        //    msl[1][2] * gn_l[k][i][1] * utl[1] +
        //    msl[2][2] * gn_l[k][i][2] * utl[2] +
        //    msr[1][3] * gn_r[k][i][1] * utr[1] +
        //    msr[2][3] * gn_r[k][i][2] * utr[2] ;
        //  // flux second component (j) is always zero for j /= 0, since rotated
        //  // frame is aligned to face-normal, i.e. fnt = (1, 0, 0)
        //}

        ////// (c) Non-conservative form with the complete normal velocity:
        //////     Mach number splitting in normal dir
        ////for (std::size_t i=0; i<3; ++i) {
        ////  for (std::size_t j=0; j<3; ++j) {
        ////    // get fluxes for g_k in rotated frame
        ////    flxn[i][j] =
        ////      l_plus * gn_l[k][i][j] + l_minus * gn_r[k][i][j];
        ////  }
        ////}

        //// rotate fluxes back to Cartesian frame
        //auto flxx = tk::unrotateTensor(flxn, fn);

        //for (std::size_t i=0; i<3; ++i)
        //  for (std::size_t j=0; j<3; ++j)
        //    flx[deformIdx(nmat,solidx[k],i,j)] = flxx[i][j];
        //// -------------------------------------------------------------------

        //// 2. rotate only g_il.u_l
        //// -------------------------------------------------------------------
        //std::array< tk::real, 3 > gn_dot_un{{ 0, 0, 0 }};
        //for (std::size_t i=0; i<3; ++i) {
        //  gn_dot_un[i] =
        //    l_plus*gn_l[k][i][0] + l_minus*gn_r[k][i][0]
        //    + l_p * gn_l[k][i][1]*utl[1] + l_m * gn_r[k][i][1]*utr[1]
        //    + l_p * gn_l[k][i][2]*utl[2] + l_m * gn_r[k][i][2]*utr[2];
        //}

        //auto g_dot_u = tk::unrotateVector(gn_dot_un, fn);

        //for (std::size_t i=0; i<3; ++i)
        //  for (std::size_t j=0; j<3; ++j)
        //    flx[deformIdx(nmat,solidx[k],i,j)] = g_dot_u[i] * fn[j];
        //// -------------------------------------------------------------------

        //// 3. Pressure splitting
        //// -------------------------------------------------------------------
        //std::array< std::array< tk::real, 3 >, 3 > flxn{{
        //  {{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}} }};
        //for (std::size_t i=0; i<3; ++i) {
        //  std::size_t j=0;
        //  flxn[i][j] =
        //    msl[0][2] * (
        //    gn_l[k][i][0] * utl[0] +
        //    gn_l[k][i][1] * utl[1] +
        //    gn_l[k][i][2] * utl[2] ) +
        //    msr[0][3] * (
        //    gn_r[k][i][0] * utr[0] +
        //    gn_r[k][i][1] * utr[1] +
        //    gn_r[k][i][2] * utr[2] );
        //}

        //// rotate fluxes back to Cartesian frame
        //auto flxx = tk::unrotateTensor(flxn, fn);

        //for (std::size_t i=0; i<3; ++i)
        //  for (std::size_t j=0; j<3; ++j)
        //    flx[deformIdx(nmat,solidx[k],i,j)] = flxx[i][j];
        //// -------------------------------------------------------------------

        // 4. Lax Friedrichs
        // -------------------------------------------------------------------
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j) {
            auto fl = (
              g_l[k][i][0] * ul +
              g_l[k][i][1] * vl +
              g_l[k][i][2] * wl ) * fn[j];
            auto fr = (
              g_r[k][i][0] * ur +
              g_r[k][i][1] * vr +
              g_r[k][i][2] * wr ) * fn[j];

            auto ll = ac12 + std::max(std::abs(utl[0]), std::abs(utr[0]));

            flx[deformIdx(nmat,solidx[k],i,j)] =
              0.5 * (fl + fr - ll*(g_r[k][i][j] - g_l[k][i][j]));
          }
        // -------------------------------------------------------------------
      }
    }

    // Stress flux in momentum eq
    std::array< tk::real, 3 > Tn12{{ 0, 0, 0 }};
    for (std::size_t i=0; i<3; ++i) {
      Tn12[i] = msl[i][2]*Trot_l[i] + msr[i][3]*Trot_r[i];
    }

    auto T12 = tk::unrotateVector(Tn12, fn);

    for (std::size_t idir=0; idir<3; ++idir)
    {
      flx[momentumIdx(nmat, idir)] =
        l_plus* u[0][momentumIdx(nmat, idir)] +
        l_minus*u[1][momentumIdx(nmat, idir)]
        - (T12[idir] - pu*fn[idir]);
    }

    l_plus = l_plus/( std::fabs(vriem) + 1.0e-12 );
    l_minus = l_minus/( std::fabs(vriem) + 1.0e-12 );

    // Store Riemann-advected partial pressures (nmat)
    if (std::fabs(l_plus) > 1.0e-10)
    {
      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back( pml[k] );
    }
    else if (std::fabs(l_minus) > 1.0e-10)
    {
      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back( pmr[k] );
    }
    else
    {
      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back( 0.5*(pml[k] + pmr[k]) );
    }

    // Store Riemann velocity (1)
    flx.push_back( vriem );

    // Store Riemann aTn_ij (3*nsld)
    if (std::fabs(l_plus) > 1.0e-10)
    {
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTn_l[k][i]);
        }
      }
    }
    else if (std::fabs(l_minus) > 1.0e-10)
    {
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(aTn_r[k][i]);
        }
      }
    }
    else
    {
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back( 0.5 * (aTn_l[k][i] + aTn_r[k][i]) );
        }
      }
    }

    Assert( flx.size() == (ncomp+nmat+1+3*nsld), "Size of multi-material "
            "flux vector incorrect" );

    return flx;
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::AUSM; }
};

} // inciter::

#endif // AUSM_h
