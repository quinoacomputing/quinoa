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
      pl += pml[k];

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
      asign_r.push_back(tk::matvec(asig_r[k], fn));
      for (std::size_t i=0; i<3; ++i)
        sign_r[i] += asign_r[k][i];

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

    // If solids exist, coordinate transform velocities to r,s,t coordinate
    // system where the first coordinate (r) aligns with face-normal (fn)
    std::array< tk::real, 3 > utl, utr;
    if (nsld > 0) {
      utl = tk::rotateVectorToN({{ul, vl, wl}}, fn);
      utr = tk::rotateVectorToN({{ur, vr, wr}}, fn);
    }
    auto vtl = std::max(utl[1], utl[2]);
    auto vtr = std::max(utr[1], utr[2]);

    // Mach numbers
    auto ml = vnl/ac12;
    auto mr = vnr/ac12;
    auto mtl = vtl/acsec12;
    auto mtr = vtr/acsec12;

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
    auto msl = splitmach_ausm( f_a, ml );
    auto msr = splitmach_ausm( f_a, mr );
    auto mtsl = splitmach_ausm( f_a, mtl );  // shear-left
    auto mtsr = splitmach_ausm( f_a, mtr );  // shear-right

    // Riemann Mach number
    auto m0 = 1.0 - (0.5*(vnl*vnl + vnr*vnr)/(ac12*ac12));
    auto mp = -k_p* std::max(m0, 0.0) * (pr-pl) / (f_a*rho12*ac12*ac12);
    auto m12 = msl[0] + msr[1] + mp;
    auto vriem = ac12 * m12;
    auto mt12 = mtsl[0] + mtsr[1];  // shear
    auto vtriem = acsec12 * mt12;  // shear

    // Pressure correction
    auto pu = -k_u* msl[2] * msr[3] * f_a * rho12 * ac12 * (vnr-vnl);

    // Flux vector splitting
    auto l_plus = 0.5 * (vriem + std::fabs(vriem));
    auto l_minus = 0.5 * (vriem - std::fabs(vriem));
    auto lt_plus = 0.5 * (vtriem + std::fabs(vtriem));  // shear-left
    auto lt_minus = 0.5 * (vtriem - std::fabs(vtriem));  // shear-right

    // Normalized flux vectors
    auto l_p = l_plus/( vriem + std::copysign(1.0e-12, vriem) );
    auto l_m = l_minus/( vriem + std::copysign(1.0e-12, vriem) );
    auto lt_p = lt_plus/( vtriem + std::copysign(1.0e-12, vtriem) );  // Lshear
    auto lt_m = lt_minus/( vtriem + std::copysign(1.0e-12, vtriem) );  // Rshear

    // Conservative fluxes
    for (std::size_t k=0; k<nmat; ++k)
    {
      flx[volfracIdx(nmat, k)] = l_plus*al_l[k] + l_minus*al_r[k];
      flx[densityIdx(nmat, k)] = l_plus*u[0][densityIdx(nmat, k)]
                              + l_minus*u[1][densityIdx(nmat, k)];
      flx[energyIdx(nmat, k)] = l_plus*u[0][energyIdx(nmat, k)]
        + l_minus*u[1][energyIdx(nmat, k)];
      // stress-work term in energy fluxes
      for (std::size_t i=0; i<3; ++i) {
        flx[energyIdx(nmat, k)] -= (
          l_p * (u[0][ncomp+velocityIdx(nmat,i)]*asign_l[k][i]) +
          l_m * (u[1][ncomp+velocityIdx(nmat,i)]*asign_r[k][i]) );
      }

      // inv deformation gradient tensor fluxes
      if (solidx[k] > 0) {

        // 1. rotate entire flux vector (with the partial deriv)
        // -------------------------------------------------------------------
        std::array< std::array< tk::real, 3 >, 3 > flxn{{
          {{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}} }};
        for (std::size_t i=0; i<3; ++i) {
          // get fluxes for g_k in rotated frame
          flxn[i][0] = l_plus*gn_l[k][i][0] + l_minus*gn_r[k][i][0]
            + mtsl[2] * gn_l[k][i][1]*utl[1] + mtsr[3] * gn_r[k][i][1]*utr[1]
            + mtsl[2] * gn_l[k][i][2]*utl[2] + mtsr[3] * gn_r[k][i][2]*utr[2];
          // flux second component (j) is always zero for j /= 0, since rotated
          // frame is aligned to face-normal, i.e. fnt = (1, 0, 0)
        }

        // rotate fluxes back to Cartesian frame
        auto flxx = tk::rotateTensorBack(flxn, fn);

        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            flx[deformIdx(nmat,solidx[k],i,j)] = flxx[i][j];
        // -------------------------------------------------------------------

        //// 2. rotate only g_il.u_l
        //// -------------------------------------------------------------------
        //std::array< tk::real, 3 > gn_dot_un{{ 0, 0, 0 }};
        //for (std::size_t i=0; i<3; ++i) {
        //  gn_dot_un[i] =
        //    l_plus*gn_l[k][i][0] + l_minus*gn_r[k][i][0]
        //    + l_p * gn_l[k][i][1]*utl[1] + l_m * gn_r[k][i][1]*utr[1]
        //    + l_p * gn_l[k][i][2]*utl[2] + l_m * gn_r[k][i][2]*utr[2];
        //}

        //auto g_dot_u = tk::rotateVectorBackFromN(gn_dot_un, fn);

        //for (std::size_t i=0; i<3; ++i)
        //  for (std::size_t j=0; j<3; ++j)
        //    flx[deformIdx(nmat,solidx[k],i,j)] = g_dot_u[i] * fn[j];
        //// -------------------------------------------------------------------

        //// 3. Pressure splitting
        //// -------------------------------------------------------------------
        //for (std::size_t i=0; i<3; ++i)
        //  for (std::size_t j=0; j<3; ++j) {
        //    flx[deformIdx(nmat,solidx[k],i,j)] =
        //      msl[2] * (
        //      g_l[k][i][0] * ul +
        //      g_l[k][i][1] * vl +
        //      g_l[k][i][2] * wl ) * fn[j] +
        //      msr[3] * (
        //      g_r[k][i][0] * ur +
        //      g_r[k][i][1] * vr +
        //      g_r[k][i][2] * wr ) * fn[j];
        //  }
        //// -------------------------------------------------------------------

        //// . upwind everything
        //// -------------------------------------------------------------------
        //for (std::size_t i=0; i<3; ++i)
        //  for (std::size_t j=0; j<3; ++j) {
        //    flx[deformIdx(nmat,solidx[k],i,j)] =
        //      l_p * (
        //      g_l[k][i][0] * ul +
        //      g_l[k][i][1] * vl +
        //      g_l[k][i][2] * wl ) * fn[j] +
        //      l_m * (
        //      g_r[k][i][0] * ur +
        //      g_r[k][i][1] * vr +
        //      g_r[k][i][2] * wr ) * fn[j];
        //  }
        //// -------------------------------------------------------------------
      }
    }

    for (std::size_t idir=0; idir<3; ++idir)
    {
      flx[momentumIdx(nmat, idir)] =
        l_plus* u[0][momentumIdx(nmat, idir)] +
        l_minus*u[1][momentumIdx(nmat, idir)]
        - (msl[2]*sign_l[idir] + msr[3]*sign_r[idir] - pu*fn[idir]);
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

    // Store Riemann asign_ij (3*nsld)
    if (std::fabs(l_plus) > 1.0e-10)
    {
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(asign_l[k][i]);
        }
      }
    }
    else if (std::fabs(l_minus) > 1.0e-10)
    {
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back(asign_r[k][i]);
        }
      }
    }
    else
    {
      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          for (std::size_t i=0; i<3; ++i)
            flx.push_back( 0.5 * (asign_l[k][i] + asign_r[k][i]) );
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

  private:

  //! Split Mach polynomials for AUSM+ flux
  //! \param[in] fa All-speed parameter
  //! \param[in] mach Local Mach numner
  //! \return Values of the positive and negative split Mach and pressure
  //!   polynomials.
  //! \details This function returns a vector with positive and negative Mach
  //!   and pressure polynomials, as:
  //!   ms[0] = M_4(+),
  //!   ms[1] = M_4(-),
  //!   ms[2] = P_5(+), and
  //!   ms[3] = P_5(-).
  //!   For more details, ref. Liou, M. S. (2006). A sequel to AUSM, Part II:
  //!   AUSM+-up for all speeds. J. Comp. Phys., 214(1), 137-170.
  static std::array< tk::real, 4 > splitmach_ausm( tk::real fa,
                                                   tk::real mach )
  {
    std::array< tk::real, 4 > ms;

    std::array< tk::real, 3 > msplus, msminus;
    tk::real psplus, psminus;

    msplus[0] = 0.5*(mach + std::fabs(mach));
    msminus[0]= 0.5*(mach - std::fabs(mach));

    msplus[1] = +0.25*(mach + 1.0)*(mach + 1.0);
    msminus[1]= -0.25*(mach - 1.0)*(mach - 1.0);

    auto alph_fa = (3.0/16.0) * (-4.0 + 5.0*fa*fa);

    if (std::fabs(mach) >= 1.0)
    {
        msplus[2] = msplus[0];
        msminus[2]= msminus[0];
        psplus    = msplus[0]/mach;
        psminus   = msminus[0]/mach;
    }
    else
    {
        msplus[2] = msplus[1]* (1.0 - 2.0*msminus[1]);
        msminus[2]= msminus[1]* (1.0 + 2.0*msplus[1]);
        psplus    = msplus[1]*
                    ((+2.0 - mach) - (16.0 * alph_fa)*mach*msminus[1]);
        psminus   = msminus[1]*
                    ((-2.0 - mach) + (16.0 * alph_fa)*mach*msplus[1]);
    }

    ms[0] = msplus[2];
    ms[1] = msminus[2];
    ms[2] = psplus;
    ms[3] = psminus;

    return ms;
  }
};

} // inciter::

#endif // AUSM_h
