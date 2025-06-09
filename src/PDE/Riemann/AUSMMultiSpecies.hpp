// *****************************************************************************
/*!
  \file      src/PDE/Riemann/AUSMMultiSpecies.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Advection Upstream Splitting Method (AUSM+) Riemann flux function
             for multi-species fluid dynamics.
  \details   This file implements the Advection Upstream Splitting Method (AUSM)
             Riemann solver, with the all-speed corrections.
             Ref. Liou, M. S. (2006). A sequel to AUSM, Part II: AUSM+-up for
             all speeds. Journal of computational physics, 214(1), 137-170.
*/
// *****************************************************************************
#ifndef AUSMMultiSpecies_h
#define AUSMMultiSpecies_h

#include <vector>

#include "Fields.hpp"
#include "SplitMachFns.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EOS.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"
#include "MultiSpecies/Mixture/Mixture.hpp"

namespace inciter {

//! AUSMMultiSpecies+up approximate Riemann solver
struct AUSMMultiSpecies {

  //! AUSM+up approximate Riemann solver flux function for multi-species flow
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
    auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

    // All-speed parameters
    // These parameters control the amount of all-speed diffusion necessary for
    // low-Mach flows. Setting k_u and k_p to zero does not add any all-speed
    // diffusion, whereas setting k_u and k_p to 1 adds maximum recommended
    // all-speed diffusion. See "Liou, M. S. (2006). A sequel to AUSM, Part II:
    // AUSM+-up for all speeds. Journal of computational physics, 214(1),
    // 137-170" for more mathematical explanation. k_u is the velocity diffusion
    // term and k_p is the pressure diffusion term. These two terms reduce
    // pressure-velocity decoupling (chequerboarding/odd-even oscillations).
    auto k_u = g_inputdeck.get< tag::lowspeed_ku >();
    auto k_p = g_inputdeck.get< tag::lowspeed_kp >();

    auto ncomp = u[0].size()-1;
    std::vector< tk::real > flx( ncomp, 0 );

    tk::real rhol(0.0), rhor(0.0), pl(0.0), pr(0.0), hl(0.0), hr(0.0),
      Tl(0.0), Tr(0.0), al(0.0), ar(0.0), a12(0.0), rho12(0.0);

    // initialize mixtures
    Mixture mixl(nspec, u[0], mat_blk);
    Mixture mixr(nspec, u[1], mat_blk);

    // Mixture densities
    rhol = mixl.get_mix_density();
    rhor = mixr.get_mix_density();
    Tl = u[0][ncomp+multispecies::temperatureIdx(nspec, 0)];
    Tr = u[1][ncomp+multispecies::temperatureIdx(nspec, 0)];

    // Velocities
    auto ul = u[0][multispecies::momentumIdx(nspec, 0)]/rhol;
    auto vl = u[0][multispecies::momentumIdx(nspec, 1)]/rhol;
    auto wl = u[0][multispecies::momentumIdx(nspec, 2)]/rhol;
    auto ur = u[1][multispecies::momentumIdx(nspec, 0)]/rhor;
    auto vr = u[1][multispecies::momentumIdx(nspec, 1)]/rhor;
    auto wr = u[1][multispecies::momentumIdx(nspec, 2)]/rhor;

    pl = mixl.pressure( rhol, Tl );
    hl = u[0][multispecies::energyIdx(nspec, 0)] + pl;
    al = mixl.frozen_soundspeed( rhol, Tl, mat_blk );

    pr = mixr.pressure( rhor, Tr );
    hr = u[1][multispecies::energyIdx(nspec, 0)] + pr;
    ar = mixr.frozen_soundspeed( rhor, Tr, mat_blk );

    // Average states for mixture speed of sound
    a12 = 0.5*(al+ar);
    rho12 = 0.5*(rhol+rhor);

    // Face-normal velocities from advective velocities
    auto vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];
    auto vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];

    // Mach numbers
    auto ml = vnl/a12;
    auto mr = vnr/a12;

    tk::real f_a(1.0);

    // Split Mach polynomials
    auto msl = splitmach_ausm( ml, f_a );
    auto msr = splitmach_ausm( mr, f_a );

    // Riemann Mach number
    auto m0 = 1.0 - (0.5*(vnl*vnl + vnr*vnr)/(a12*a12));
    auto mp = -k_p* std::max(m0, 0.0) * (pr-pl) / (f_a*rho12*a12*a12);
    auto m12 = msl[0] + msr[1] + mp;
    auto vriem = a12 * m12;

    // Riemann pressure
    auto pu = -k_u* msl[2] * msr[3] * f_a * rho12 * a12 * (vnr-vnl);
    auto p12 = msl[2]*pl + msr[3]*pr + pu;

    // Flux vector splitting
    auto l_plus = 0.5 * (vriem + std::fabs(vriem));
    auto l_minus = 0.5 * (vriem - std::fabs(vriem));

    // Conservative fluxes
    for (std::size_t k=0; k<nspec; ++k)
    {
      flx[multispecies::densityIdx(nspec, k)] =
        l_plus *u[0][multispecies::densityIdx(nspec, k)] +
        l_minus*u[1][multispecies::densityIdx(nspec, k)];
    }

    flx[multispecies::energyIdx(nspec, 0)] = l_plus*hl + l_minus*hr;

    for (std::size_t idir=0; idir<3; ++idir)
    {
      flx[multispecies::momentumIdx(nspec, idir)] =
        l_plus *u[0][multispecies::momentumIdx(nspec, idir)] +
        l_minus*u[1][multispecies::momentumIdx(nspec, idir)] + p12*fn[idir];
    }

    return flx;
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::AUSM; }
};

} // inciter::

#endif // AUSMMultiSpecies_h
