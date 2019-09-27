// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Riemann/AUSM.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
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

#include "Types.hpp"
#include "Fields.hpp"
#include "Tags.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EoS.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

//! AUSM+up approximate Riemann solver
//! \details This class can be used polymorphically with inciter::RiemannSolver
struct AUSM {

  //! AUSM+up approximate Riemann solver flux function
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann flux solution according to AUSM+up, appended by Riemann
  //!   velocities and volume-fractions.
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::vector< std::array< tk::real, 3 > >& )
  {
    const auto nmat =
      g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[0];

    auto ncomp = u[0].size()-(3+nmat);
    std::vector< tk::real > flx( ncomp, 0 );

    // Primitive variables
    tk::real rhol(0.0), rhor(0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      rhol += u[0][densityIdx(nmat, k)];
      rhor += u[1][densityIdx(nmat, k)];
    }

    auto ul = u[0][momentumIdx(nmat, 0)]/rhol;
    auto vl = u[0][momentumIdx(nmat, 1)]/rhol;
    auto wl = u[0][momentumIdx(nmat, 2)]/rhol;

    auto ur = u[1][momentumIdx(nmat, 0)]/rhor;
    auto vr = u[1][momentumIdx(nmat, 1)]/rhor;
    auto wr = u[1][momentumIdx(nmat, 2)]/rhor;

    tk::real pl(0.0), pr(0.0), amatl(0.0), amatr(0.0);
    std::vector< tk::real > al_l(nmat, 0.0), al_r(nmat, 0.0),
                            hml(nmat, 0.0), hmr(nmat, 0.0),
                            pml(nmat, 0.0), pmr(nmat, 0.0),
                            al_12(nmat, 0.0), rhomat12(nmat, 0.0),
                            amat12(nmat, 0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      al_l[k] = u[0][volfracIdx(nmat, k)];
      pml[k] = u[0][ncomp+pressureIdx(nmat, k)];
      pl += al_l[k] * pml[k];
      hml[k] = u[0][energyIdx(nmat, k)] + al_l[k]*pml[k];
      amatl = eos_soundspeed< tag::multimat >( 0,
                                               u[0][densityIdx(nmat, k)]/al_l[k],
                                               pml[k], k );

      al_r[k] = u[1][volfracIdx(nmat, k)];
      pmr[k] = u[1][ncomp+pressureIdx(nmat, k)];
      pr += al_r[k] * pmr[k];
      hmr[k] = u[1][energyIdx(nmat, k)] + al_r[k]*pmr[k];
      amatr = eos_soundspeed< tag::multimat >( 0,
                                               u[1][densityIdx(nmat, k)]/al_r[k],
                                               pmr[k], k );

      // Average states for mixture speed of sound
      al_12[k] = 0.5*(al_l[k]+al_r[k]);
      rhomat12[k] = 0.5*(u[0][densityIdx(nmat, k)]/al_l[k]
                        + u[1][densityIdx(nmat, k)]/al_r[k]);
      amat12[k] = 0.5*(amatl+amatr);
    }

    auto rho12 = 0.5*(rhol+rhor);

    // mixture speed of sound
    tk::real ac12(0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      ac12 += (al_12[k]*rhomat12[k]*amat12[k]*amat12[k]);
    }
    ac12 = std::sqrt( ac12/rho12 );

    // Independently limited velocities for advection
    ul = u[0][ncomp+velocityIdx(nmat, 0)];
    vl = u[0][ncomp+velocityIdx(nmat, 1)];
    wl = u[0][ncomp+velocityIdx(nmat, 2)];
    ur = u[1][ncomp+velocityIdx(nmat, 0)];
    vr = u[1][ncomp+velocityIdx(nmat, 1)];
    wr = u[1][ncomp+velocityIdx(nmat, 2)];

    // Face-normal velocities from advective velocities
    auto vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];
    auto vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];

    // Mach numbers
    auto ml = vnl/ac12;
    auto mr = vnr/ac12;

    // All-speed parameters
    tk::real k_u(0.0), k_p(0.0), f_a(1.0);

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

    l_plus = l_plus/( std::fabs(vriem) + 1.0e-16 );
    l_minus = l_minus/( std::fabs(vriem) + 1.0e-16 );

    // Store Riemann-advected partial pressures
    if (std::fabs(l_plus) > 1.0e-10)
    {
      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back( al_l[k]*pml[k] );
    }
    else if (std::fabs(l_minus) > 1.0e-10)
    {
      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back( al_r[k]*pmr[k] );
    }
    else
    {
      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back( 0.5*(al_l[k]*pml[k] + al_r[k]*pmr[k]) );
    }

    // Store Riemann velocity
    flx.push_back( vriem );

    Assert( flx.size() == (3*nmat+3+nmat+1), "Size of multi-material flux "
            "vector incorrect" );

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
