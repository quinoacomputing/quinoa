// *****************************************************************************
/*!
  \file      src/PDE/Riemann/LDFSS.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Low Diffusion Flux Splitting Scheme (LDFSS+) Riemann flux function
  \details   This file implements the modified Low Diffusion Flux Splitting
             Scheme (LDFSS) Riemann solver. See Edwards, J. (2001, June).
             Towards unified CFD simulations of real fluid flows. In 15th AIAA
             computational fluid dynamics conference (p. 2524).
*/
// *****************************************************************************
#ifndef LDFSS_h
#define LDFSS_h

#include <vector>

#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "SplitMachFns.hpp"
#include "EoS/EOS.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"
#include "MultiSpecies/Mixture/Mixture.hpp"

namespace inciter {

//! LDFSS approximate Riemann solver
struct LDFSS {

  //! LDFSS approximate Riemann solver flux function
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann flux solution according to LDFSS, appended by Riemann
  //!   velocities and volume-fractions.
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::vector< EOS >& mat_blk,
        const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::vector< std::array< tk::real, 3 > >& = {},
        const std::array< tk::real, 3>& = {} )
  {
    auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

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

    // Split Mach polynomials
    auto msl = splitmach_ausm( ml, 1.0, 0.0, 0.0 );
    auto msr = splitmach_ausm( mr, 1.0, 0.0, 0.0 );

    // Modified Mach polynomials
    // mtilde for general-order Mach polys (AUSM+)
    auto mtilde  = 0.5*(msl[0] - std::max(0.0,ml) - msr[1] + std::min(0.0,mr));
    //// Uncomment for LDFSS-2025M
    //// Modified mtilde for 2nd order Mach polys (vanLeer)
    //auto beta_l = -std::max(0.0, 1.0 - std::floor(std::abs(ml)));
    //auto beta_r = -std::max(0.0, 1.0 - std::floor(std::abs(mr)));
    //auto mtilde = 0.25*beta_l*beta_r*std::pow((std::pow(0.5*(ml*ml+mr*mr),0.25) - 1.0), 2.0);

    auto dp = pl - pr;
    auto delta = 4.0;

    // LDFSS-2001
    // Additional diffusion
    auto dpplus = std::abs(dp);
    //dpplus += 8.0*rhol*a12*std::sqrt(std::abs(vnl - vnr)*a12);  // 2025u mod
    dpplus += 8.0*rhol*a12*std::sqrt(std::abs(dp)/rho12);  // 2025p mod
    auto dpminus = std::abs(dp);
    //dpminus += 8.0*rhor*a12*std::sqrt(std::abs(vnl - vnr)*a12);  // 2025u mod
    dpminus += 8.0*rhor*a12*std::sqrt(std::abs(dp)/rho12);  // 2025p mod
    auto mtilde_plus = mtilde*
      std::max( 0.0, 1.0 - (dp+delta*dpplus)/(2.0*rhol*a12*a12) );
    auto mtilde_minus = mtilde*
      std::max( 0.0, 1.0 + (dp-delta*dpminus)/(2.0*rhor*a12*a12) );

    //// LDFSS(2) with the 2025p mod
    //auto dpplus = std::abs(dp)
    //  + 4.0*pl*std::sqrt(2.0*std::abs(dp)/(pl+pr));
    //auto dpminus = std::abs(dp)
    //  + 4.0*pr*std::sqrt(2.0*std::abs(dp)/(pl+pr));
    //auto mtilde_plus = mtilde*
    //  std::max( 0.0, 1.0 - dp/(pl+pr) - delta*dpplus/pl );
    //auto mtilde_minus = mtilde*
    //  std::max( 0.0, 1.0 + dp/(pl+pr) - delta*dpminus/pr );

    auto C_plus = msl[0] - mtilde_plus;
    auto C_minus = msr[1] + mtilde_minus;

    // Flux vector splitting
    auto l_plus = C_plus * a12;
    auto l_minus = C_minus * a12;

    // Riemann pressure
    auto p12 = 0.5*(pl+pr) + 0.5*(pl-pr)*(msl[2]-msr[3]) +
      rho12*a12*a12*(msl[2]+msr[3]-1.0);

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
  static ctr::FluxType type() noexcept { return ctr::FluxType::LDFSS; }
};

} // inciter::

#endif // LDFSS_h
