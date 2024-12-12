// *****************************************************************************
/*!
  \file      src/PDE/Riemann/SLAU.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Simple Low-dissipation AUSM (SLAU2) Riemann flux function
  \details   This file implements the (SLAU2) Riemann solver.
             Ref. Kitamura, K., & Shima, E. (2013). Towards shock-stable and
             accurate hypersonic heating computations: A new pressure flux for
             AUSM-family schemes. Journal of Computational Physics, 245, 62-83.
*/
// *****************************************************************************
#ifndef SLAU_h
#define SLAU_h

#include <vector>

#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "SplitMachFns.hpp"
#include "EoS/EOS.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

//! SLAU2 approximate Riemann solver
struct SLAU {

  //! SLAU2 approximate Riemann solver flux function
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann flux solution according to SLAU2, appended by Riemann
  //!   velocities and volume-fractions.
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::vector< EOS >& mat_blk,
        const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::vector< std::array< tk::real, 3 > >& = {} )
  {
    auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();

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

    // Split Mach polynomials
    auto msl = splitmach_ausm( 1.0, ml );
    auto msr = splitmach_ausm( 1.0, mr );

    auto mtilde = std::min( 1.0,
      std::sqrt(0.5*(ul*ul+vl*vl+wl*wl + ur*ur+vr*vr+wr*wr))/ac12 );
    auto kai = (1.0 - mtilde)*(1.0 - mtilde);

    // Riemann Mach number
    auto vnmod = (rhol*std::abs(vnl) + rhor*std::abs(vnr)) / (rhol+rhor);
    auto gg = -std::max(std::min(0.0,ml), -1.0) *
      std::min(std::max(0.0,mr), 1.0);
    auto vnplus = (1.0-gg)*vnmod + gg*std::abs(vnl);
    auto vnminus= (1.0-gg)*vnmod + gg*std::abs(vnr);
    auto massdot = 0.5*( rhol*(vnl+vnplus) + rhor*(vnr-vnminus)
      - kai*(pr-pl)/ac12 );

    // Riemann pressure
    auto p12 = 0.5*(pl+pr) + 0.5*(msl[2]-msr[3])*(pl-pr) +
      std::sqrt(0.5*(ul*ul+vl*vl+wl*wl + ur*ur+vr*vr+wr*wr))*(msl[2]+msr[3]-1.0)
      *rho12*ac12;

    // Flux vector splitting
    auto l_plus = 0.5 * (massdot + std::fabs(massdot));
    auto l_minus = 0.5 * (massdot - std::fabs(massdot));

    // Conservative fluxes
    for (std::size_t k=0; k<nmat; ++k)
    {
      //flx[volfracIdx(nmat, k)] = l_plus*al_l[k] + l_minus*al_r[k];
      flx[densityIdx(nmat, k)] = l_plus + l_minus;
      flx[energyIdx(nmat, k)] = l_plus*hml[k]/u[0][densityIdx(nmat, k)]
        + l_minus*hmr[k]/u[1][densityIdx(nmat, k)];
    }

    flx[momentumIdx(nmat, 0)] =
      l_plus*ul + l_minus*ur + p12*fn[0];
    flx[momentumIdx(nmat, 1)] =
      l_plus*vl + l_minus*vr + p12*fn[1];
    flx[momentumIdx(nmat, 2)] =
      l_plus*wl + l_minus*wr + p12*fn[2];

    l_plus = l_plus/( std::fabs(massdot) + 1.0e-12 );
    l_minus = l_minus/( std::fabs(massdot) + 1.0e-12 );

    // Store Riemann-advected partial pressures
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

    // Store Riemann velocity
    flx.push_back( 0.0 );

    Assert( flx.size() == (3*nmat+3+nmat+1), "Size of multi-material flux "
            "vector incorrect" );

    return flx;
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::SLAU; }
};

} // inciter::

#endif // SLAU_h
