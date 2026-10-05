// *****************************************************************************
/*!
  \file      src/PDE/Riemann/HLLMultiSpecies.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Harten-Lax-van Leer (HLL) Riemann flux function for multi-species
             fluid dynamics
  \details   This file implements the Harten-Lax-van Leer (HLL) Riemann solver
             for multi-species fluid dynamics.
*/
// *****************************************************************************
#ifndef HLLMultiSpecies_h
#define HLLMultiSpecies_h

#include <vector>

#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EOS.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"
#include "MultiSpecies/Mixture/Mixture.hpp"

namespace inciter {

//! HLL approximate Riemann solver for multi-species flow
struct HLLMultiSpecies {

  //! HLL approximate Riemann solver flux function for multi-species flow
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann flux solution according to HLL
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::vector< EOS >& mat_blk,
        const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::vector< std::array< tk::real, 3 > >& = {},
        const tk::real = 0 )
  {
    auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();
    auto ncomp = u[0].size()-1;
    std::vector< tk::real > flx( ncomp, 0 ), fl( ncomp, 0 ), fr( ncomp, 0 );

    // Initialize mixtures and primitive variables
    Mixture mixl(nspec, u[0], mat_blk);
    Mixture mixr(nspec, u[1], mat_blk);

    auto rhol = mixl.get_mix_density();
    auto rhor = mixr.get_mix_density();
    auto Tl = u[0][ncomp+multispecies::temperatureIdx(nspec, 0)];
    auto Tr = u[1][ncomp+multispecies::temperatureIdx(nspec, 0)];

    auto ul = u[0][multispecies::momentumIdx(nspec, 0)] / rhol;
    auto vl = u[0][multispecies::momentumIdx(nspec, 1)] / rhol;
    auto wl = u[0][multispecies::momentumIdx(nspec, 2)] / rhol;
    auto ur = u[1][multispecies::momentumIdx(nspec, 0)] / rhor;
    auto vr = u[1][multispecies::momentumIdx(nspec, 1)] / rhor;
    auto wr = u[1][multispecies::momentumIdx(nspec, 2)] / rhor;

    auto pl = mixl.pressure( rhol, Tl );
    auto pr = mixr.pressure( rhor, Tr );
    auto al = mixl.frozen_soundspeed( rhol, Tl, mat_blk );
    auto ar = mixr.frozen_soundspeed( rhor, Tr, mat_blk );

    // Face-normal velocities
    auto vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];
    auto vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];

    // Conservative flux functions
    for (std::size_t k=0; k<nspec; ++k) {
      fl[multispecies::densityIdx(nspec, k)] =
        vnl * u[0][multispecies::densityIdx(nspec, k)];
      fr[multispecies::densityIdx(nspec, k)] =
        vnr * u[1][multispecies::densityIdx(nspec, k)];
    }

    for (std::size_t idir=0; idir<3; ++idir) {
      auto imom = multispecies::momentumIdx(nspec, idir);
      fl[imom] = vnl*u[0][imom] + pl*fn[idir];
      fr[imom] = vnr*u[1][imom] + pr*fn[idir];
    }

    auto iene = multispecies::energyIdx(nspec, 0);
    fl[iene] = vnl * (u[0][iene] + pl);
    fr[iene] = vnr * (u[1][iene] + pr);

    // Signal velocities
    auto Sl = std::min( vnl-al, vnr-ar );
    auto Sr = std::max( vnl+al, vnr+ar );

    // Numerical flux
    if (Sl >= 0.0) {
      flx = fl;
    }
    else if (Sr <= 0.0) {
      flx = fr;
    }
    else {
      for (std::size_t k=0; k<ncomp; ++k)
        flx[k] = (Sr*fl[k] - Sl*fr[k] + Sl*Sr*(u[1][k]-u[0][k])) / (Sr-Sl);
    }

    return flx;
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::HLL; }
};

} // inciter::

#endif // HLLMultiSpecies_h
