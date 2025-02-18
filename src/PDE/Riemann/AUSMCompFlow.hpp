// *****************************************************************************
/*!
  \file      src/PDE/Riemann/AUSMCompFlow.hpp
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
#ifndef AUSMCompFlow_h
#define AUSMCompFlow_h

#include <vector>

#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "SplitMachFns.hpp"
#include "EoS/EOS.hpp"

namespace inciter {

//! AUSM+up approximate Riemann solver
struct AUSMCompFlow {

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
    auto k_p = g_inputdeck.get< tag::lowspeed_kp >();

    auto ncomp = u[0].size();
    std::vector< tk::real > flx( ncomp, 0 );

    // Primitive variables
    auto rhol = u[0][0];
    auto ul = u[0][1]/rhol;
    auto vl = u[0][2]/rhol;
    auto wl = u[0][3]/rhol;
    auto rhor = u[1][0];
    auto ur = u[1][1]/rhor;
    auto vr = u[1][2]/rhor;
    auto wr = u[1][3]/rhor;

    tk::real pl(0.0), pr(0.0), amatl(0.0), amatr(0.0);

    pl = mat_blk[0].compute< EOS::pressure >(rhol, ul, vl, wl, u[0][4]);
    amatl = mat_blk[0].compute< EOS::soundspeed >( rhol, pl );

    pr = mat_blk[0].compute< EOS::pressure >(rhor, ur, vr, wr, u[1][4]);
    amatr = mat_blk[0].compute< EOS::soundspeed >( rhor, pr );

    // Average states for mixture speed of sound
    auto ac12 = 0.5*(amatl+amatr);
    auto rho12 = 0.5*(rhol+rhor);

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
    auto msl = splitmach_ausm( ml, f_a );
    auto msr = splitmach_ausm( mr, f_a );

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
    flx[0] = l_plus*u[0][0] + l_minus*u[1][0];

    flx[1] = l_plus*u[0][1] + l_minus*u[1][1] + p12*fn[0];
    flx[2] = l_plus*u[0][2] + l_minus*u[1][2] + p12*fn[1];
    flx[3] = l_plus*u[0][3] + l_minus*u[1][3] + p12*fn[2];

    flx[4] = l_plus*(u[0][4] + pl) + l_minus*(u[1][4] + pr);

    return flx;
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::AUSM; }
};

} // inciter::

#endif // AUSMCompFlow_h
