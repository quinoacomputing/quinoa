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

#include "Types.hpp"
#include "Fields.hpp"
#include "Tags.hpp"
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
    const auto nmat =
      g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[0];

    auto ncomp = u[0].size()-(3+nmat);
    std::vector< tk::real > flx(ncomp, 0), fl(ncomp, 0), fr(ncomp, 0);

    // Primitive variables
    tk::real rhol(0.0), rhor(0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      rhol += u[0][densityIdx(nmat, k)];
      rhor += u[1][densityIdx(nmat, k)];
    }

    tk::real pl(0.0), pr(0.0), amatl(0.0), amatr(0.0), ac_l(0.0), ac_r(0.0);
    std::vector< tk::real > al_l(nmat, 0.0), al_r(nmat, 0.0),
                            hml(nmat, 0.0), hmr(nmat, 0.0),
                            pml(nmat, 0.0), pmr(nmat, 0.0);
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

      // Mixture speed of sound
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

    // Flux functions
    for (std::size_t k=0; k<nmat; ++k)
    {
      fl[volfracIdx(nmat, k)] = vnl * al_l[k];
      fl[densityIdx(nmat, k)] = vnl * u[0][densityIdx(nmat, k)];
      fl[energyIdx(nmat, k)] = vnl * hml[k];

      fr[volfracIdx(nmat, k)] = vnr * al_r[k];
      fr[densityIdx(nmat, k)] = vnr * u[1][densityIdx(nmat, k)];
      fr[energyIdx(nmat, k)] = vnr * hmr[k];
    }

    for (std::size_t idir=0; idir<3; ++idir)
    {
      fl[momentumIdx(nmat, idir)] = vnl * u[0][momentumIdx(nmat, idir)]
        + pl*fn[idir];

      fr[momentumIdx(nmat, idir)] = vnr * u[1][momentumIdx(nmat, idir)]
        + pr*fn[idir];
    }

    // Signal velocities
    auto Sl = std::min((vnl-ac_l), (vnr-ac_r));
    auto Sr = std::min((vnl+ac_l), (vnr+ac_r));

    // Numerical flux functions and wave-speeds
    auto c_plus(0.0), c_minus(0.0), p_plus(0.0), p_minus(0.0);
    if (Sl >= 0.0)
    {
      for (std::size_t k=0; k<flx.size(); ++k)
        flx[k] = fl[k];
      c_plus = vnl;
      p_plus = 1.0;
    }
    else if (Sr <= 0.0)
    {
      for (std::size_t k=0; k<flx.size(); ++k)
        flx[k] = fr[k];
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

    auto vriem = c_plus+c_minus;

    c_plus = c_plus/( std::fabs(vriem) + 1.0e-16 );
    c_minus = c_minus/( std::fabs(vriem) + 1.0e-16 );

    // Store Riemann-advected partial pressures
    for (std::size_t k=0; k<nmat; ++k)
      flx.push_back(p_plus*pml[k] + p_minus*pmr[k]);

    // Store Riemann velocity
    flx.push_back( vriem );

    Assert( flx.size() == (3*nmat+3+nmat+1), "Size of multi-material flux "
            "vector incorrect" );

    return flx;
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::HLL; }

};

} // inciter::

#endif // HLL_h
