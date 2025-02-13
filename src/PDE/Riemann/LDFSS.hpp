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
#include "MultiMat/MultiMatIndexing.hpp"

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
    auto msl = splitmach_ausm( ml, 1.0, 0.0, 0.0 );
    auto msr = splitmach_ausm( mr, 1.0, 0.0, 0.0 );

    // Modified Mach polynomials
    //// mtilde for 2nd order Mach polys (vanLeer)
    //auto beta_l = -std::max(0.0, 1.0 - std::floor(std::abs(ml)));
    //auto beta_r = -std::max(0.0, 1.0 - std::floor(std::abs(mr)));
    //auto mtilde = 0.25*beta_l*beta_r*std::pow((std::pow(0.5*(ml*ml+mr*mr),0.25) - 1.0), 2.0);
    // mtilde for general-order Mach polys (AUSM+)
    auto mtilde  = 0.5*(msl[0] - std::max(0.0,ml) - msr[1] + std::min(0.0,mr));
    auto dp = pl - pr;
    auto delta = 4.0;

    // LDFSS-2001
    // Additional diffusion
    auto dpplus = std::abs(dp)
      //+ 8.0*rhol*ac12*std::sqrt(std::abs(vnl - vnr)*ac12);
      + 8.0*rhol*ac12*std::sqrt(std::abs(dp)/rho12);
      //;
    auto dpminus = std::abs(dp)
      //+ 8.0*rhor*ac12*std::sqrt(std::abs(vnl - vnr)*ac12);
      + 8.0*rhor*ac12*std::sqrt(std::abs(dp)/rho12);
      //;
    auto mtilde_plus = mtilde*
      std::max( 0.0, 1.0 - (dp+delta*dpplus)/(2.0*rhol*ac12*ac12) );
    auto mtilde_minus = mtilde*
      std::max( 0.0, 1.0 + (dp-delta*dpminus)/(2.0*rhor*ac12*ac12) );

    //// LDFSS(2)
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
    auto l_plus = C_plus * ac12;
    auto l_minus = C_minus * ac12;
    auto vriem = l_plus + l_minus;

    // Riemann pressure
    auto p12 = 0.5*(pl+pr) + 0.5*(pl-pr)*(msl[2]-msr[3]) +
      rho12*ac12*ac12*(msl[2]+msr[3]-1.0);

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

    l_plus = l_plus/( vnl + std::copysign(1.0e-12, vnl) );
    l_minus = l_minus/( vnr + std::copysign(1.0e-12, vnr) );

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
  static ctr::FluxType type() noexcept { return ctr::FluxType::LDFSS; }
};

} // inciter::

#endif // LDFSS_h
