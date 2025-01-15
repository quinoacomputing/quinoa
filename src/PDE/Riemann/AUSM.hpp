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

    auto ncomp = u[0].size()-(3+nmat);
    std::vector< tk::real > flx( ncomp, 0 );

    // Primitive variables
    tk::real rhol(0.0), rhor(0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      rhol += u[0][densityIdx(nmat, k)];
      rhor += u[1][densityIdx(nmat, k)];
    }

    tk::real pl(0.0), pr(0.0), amatl(0.0), amatr(0.0);
    std::vector< tk::real > al_l(nmat, 0.0), al_r(nmat, 0.0),
                            hml(nmat, 0.0), hmr(nmat, 0.0),
                            pml(nmat, 0.0), pmr(nmat, 0.0),
                            arhom12(nmat, 0.0),
                            amat12(nmat, 0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      al_l[k] = u[0][volfracIdx(nmat, k)];
      pml[k] = u[0][ncomp+pressureIdx(nmat, k)];
      pl += pml[k];
      hml[k] = u[0][energyIdx(nmat, k)] + pml[k];
      amatl = mat_blk[k].compute< EOS::soundspeed >(
        u[0][densityIdx(nmat, k)], pml[k], al_l[k], k );

      al_r[k] = u[1][volfracIdx(nmat, k)];
      pmr[k] = u[1][ncomp+pressureIdx(nmat, k)];
      pr += pmr[k];
      hmr[k] = u[1][energyIdx(nmat, k)] + pmr[k];
      amatr = mat_blk[k].compute< EOS::soundspeed >(
        u[1][densityIdx(nmat, k)], pmr[k], al_r[k], k );

      // Average states for mixture speed of sound
      arhom12[k] = 0.5*(u[0][densityIdx(nmat, k)] + u[1][densityIdx(nmat, k)]);
      amat12[k] = 0.5*(amatl+amatr);
    }

    auto rho12 = 0.5*(rhol+rhor);

    // mixture speed of sound
    tk::real ac12(0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      ac12 += (arhom12[k]*amat12[k]*amat12[k]);
    }
    ac12 = std::sqrt( ac12/rho12 );

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

    // Mach numbers
    auto ml = vnl/ac12;
    auto mr = vnr/ac12;

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

    // Riemann Mach number
    auto m0 = 1.0 - (0.5*(vnl*vnl + vnr*vnr)/(ac12*ac12));
    auto mp = -k_p* std::max(m0, 0.0) * (pr-pl) / (f_a*rho12*ac12*ac12);
    auto m12 = msl[0] + msr[1] + mp;
    auto vriem = ac12 * m12;

    // Riemann pressure
    auto pu = -k_u* msl[2] * msr[3] * f_a * rho12 * ac12 * (vnr-vnl);
    auto p12 = msl[2]*pl + msr[3]*pr + pu;

    // Flux vector splitting
    auto l_plus = 0.5 * (vriem + std::fabs(vriem));
    auto l_minus = 0.5 * (vriem - std::fabs(vriem));

    // Conservative fluxes
    for (std::size_t k=0; k<nmat; ++k)
    {
      flx[volfracIdx(nmat, k)] = l_plus*al_l[k] + l_minus*al_r[k];
      flx[densityIdx(nmat, k)] = l_plus*u[0][densityIdx(nmat, k)]
                              + l_minus*u[1][densityIdx(nmat, k)];
      flx[energyIdx(nmat, k)] = l_plus*hml[k] + l_minus*hmr[k];
    }

    for (std::size_t idir=0; idir<3; ++idir)
    {
    flx[momentumIdx(nmat, idir)] = l_plus*u[0][momentumIdx(nmat, idir)]
                                 + l_minus*u[1][momentumIdx(nmat, idir)]
                                 + p12*fn[idir];
    }

    l_plus = l_plus/( vriem + std::copysign(1.0e-12,vriem) );
    l_minus = l_minus/( vriem + std::copysign(1.0e-12,vriem) );

    // Store Riemann-advected partial pressures
    for (std::size_t k=0; k<nmat; ++k)
      flx.push_back( l_plus*pml[k] + l_minus*pmr[k] );

    // Store Riemann velocity
    flx.push_back( vriem );

    Assert( flx.size() == (3*nmat+3+nmat+1), "Size of multi-material flux "
            "vector incorrect" );

    return flx;
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::AUSM; }
};

} // inciter::

#endif // AUSM_h
