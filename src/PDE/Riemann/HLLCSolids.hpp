// *****************************************************************************
/*!
  \file      src/PDE/Riemann/HLLCSolids.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2023 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     HLLC Riemann flux function for solids
  \details   This file implements the HLLC Riemann solver for solids.
             Ref. Ndanou, S., Favrie, N., & Gavrilyuk, S. (2015). Multi-solid
             and multi-fluid diffuse interface model: Applications to dynamic
             fracture and fragmentation. Journal of Computational Physics, 295,
             523-555.
*/
// *****************************************************************************
#ifndef HLLCSolids_h
#define HLLCSolids_h

#include <vector>

#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EOS.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "MultiMat/MiscMultiMatFns.hpp"

namespace inciter {

//! HLLC approximate Riemann solver for solids
struct HLLCSolids {

//! HLLC approximate Riemann solver flux function
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann solution according to Harten-Lax-van Leer-Contact
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::vector< EOS >& mat_blk,
        const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::vector< std::array< tk::real, 3 > >& = {} )
  {
    auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
    const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

    auto nsld = numSolids(nmat, solidx);
    auto ncomp = u[0].size()-(3+nmat+nsld*6);
    std::vector< tk::real > flx(ncomp, 0), fl(ncomp, 0), fr(ncomp, 0),
      ftl(ncomp, 0), ftr(ncomp, 0);

    // Primitive variables
    tk::real rhol(0.0), rhor(0.0);
    for (size_t k=0; k<nmat; ++k) {
      rhol += u[0][densityIdx(nmat, k)];
      rhor += u[1][densityIdx(nmat, k)];
    }

    auto ul = u[0][ncomp+velocityIdx(nmat, 0)];
    auto vl = u[0][ncomp+velocityIdx(nmat, 1)];
    auto wl = u[0][ncomp+velocityIdx(nmat, 2)];
    auto ur = u[1][ncomp+velocityIdx(nmat, 0)];
    auto vr = u[1][ncomp+velocityIdx(nmat, 1)];
    auto wr = u[1][ncomp+velocityIdx(nmat, 2)];

    tk::real pl(0.0), pr(0.0);
    std::vector< tk::real > apl(nmat, 0.0), apr(nmat, 0.0);
    tk::real acl(0.0), acr(0.0);
    for (std::size_t k=0; k<nmat; ++k) {
      apl[k] = mat_blk[k].compute< EOS::pressure >( u[0][densityIdx(nmat, k)],
                                                    ul, vl, wl,
                                                    u[0][energyIdx(nmat, k)] );
      apr[k] = mat_blk[k].compute< EOS::pressure >( u[1][densityIdx(nmat, k)],
                                                    ur, vr, wr,
                                                    u[1][energyIdx(nmat, k)] );
      pl += apl[k];
      pr += apr[k];
      auto amatl = mat_blk[k].compute< EOS::soundspeed >(
        u[0][densityIdx(nmat, k)], apl[k], u[0][volfracIdx(nmat, k)], k );
      auto amatr = mat_blk[k].compute< EOS::soundspeed >(
        u[1][densityIdx(nmat, k)], apr[k], u[1][volfracIdx(nmat, k)], k );
      acl += u[0][densityIdx(nmat, k)] * amatl * amatl;
      acr += u[1][densityIdx(nmat, k)] * amatr * amatr;
    }
    acl = std::sqrt(acl/rhol);
    acr = std::sqrt(acr/rhor);

    // Face-normal velocities
    tk::real vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];
    tk::real vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];

    // Roe-averaged variables
    auto rlr = sqrt(rhor/rhol);
    auto rlr1 = 1.0 + rlr;

    auto vnroe = (vnr*rlr + vnl)/rlr1 ;
    auto aroe = (acr*rlr + acl)/rlr1 ;

    // Signal velocities
    auto Sl = fmin(vnl-acl, vnroe-aroe);
    auto Sr = fmax(vnr+acr, vnroe+aroe);
    auto Sm = ( rhor*vnr*(Sr-vnr) - rhol*vnl*(Sl-vnl) + pl-pr )
             /( rhor*(Sr-vnr) - rhol*(Sl-vnl) );

    // Middle-zone (star) variables
    tk::real pStar(0.0);
    std::vector< tk::real > apStar(nmat, 0.0);
    for (std::size_t k=0; k<nmat; ++k) {
      apStar[k] = u[0][densityIdx(nmat, k)]*(vnl-Sl)*(vnl-Sm) + apl[k];
      pStar += apStar[k];
    }
    auto uStar = u;

    uStar[0][ncomp+velocityIdx(nmat, 0)] =
      ((Sl-vnl)*u[0][ncomp+velocityIdx(nmat, 0)] + (pStar-pl)*fn[0])/(Sl-Sm);
    uStar[0][ncomp+velocityIdx(nmat, 1)] =
      ((Sl-vnl)*u[0][ncomp+velocityIdx(nmat, 1)] + (pStar-pl)*fn[1])/(Sl-Sm);
    uStar[0][ncomp+velocityIdx(nmat, 2)] =
      ((Sl-vnl)*u[0][ncomp+velocityIdx(nmat, 2)] + (pStar-pl)*fn[2])/(Sl-Sm);
    uStar[1][ncomp+velocityIdx(nmat, 0)] =
      ((Sr-vnr)*u[1][ncomp+velocityIdx(nmat, 0)] + (pStar-pr)*fn[0])/(Sr-Sm);
    uStar[1][ncomp+velocityIdx(nmat, 1)] =
      ((Sr-vnr)*u[1][ncomp+velocityIdx(nmat, 1)] + (pStar-pr)*fn[1])/(Sr-Sm);
    uStar[1][ncomp+velocityIdx(nmat, 2)] =
      ((Sr-vnr)*u[1][ncomp+velocityIdx(nmat, 2)] + (pStar-pr)*fn[2])/(Sr-Sm);

    for (std::size_t k=0; k<nmat; ++k) {
      uStar[0][volfracIdx(nmat, k)] = u[0][volfracIdx(nmat, k)];
      uStar[0][densityIdx(nmat, k)] =
        (Sl-vnl) * u[0][densityIdx(nmat, k)] / (Sl-Sm);
      uStar[0][energyIdx(nmat, k)] =
        ((Sl-vnl) * u[0][energyIdx(nmat, k)] - apl[k]*vnl + apStar[k]*Sm) / (Sl-Sm);

      uStar[1][volfracIdx(nmat, k)] = u[1][volfracIdx(nmat, k)];
      uStar[1][densityIdx(nmat, k)] =
        (Sr-vnr) * u[1][densityIdx(nmat, k)] / (Sr-Sm);
      uStar[1][energyIdx(nmat, k)] =
         ((Sr-vnr) * u[1][energyIdx(nmat, k)] - apr[k]*vnr + apStar[k]*Sm) / (Sr-Sm);
    }

    // Numerical fluxes
    if (Sl > 0.0) {

      for (std::size_t idir=0; idir<3; ++idir)
        flx[momentumIdx(nmat, idir)] =
          u[0][ncomp+velocityIdx(nmat, idir)] * vnl + pl*fn[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = u[0][volfracIdx(nmat, k)] * vnl;
        flx[densityIdx(nmat, k)] = u[0][densityIdx(nmat, k)] * vnl;
        flx[energyIdx(nmat, k)] = (u[0][energyIdx(nmat, k)] + apl[k]) * vnl;
      }

      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back(apl[k]);
      flx.push_back(vnl);
    }
    else if (Sl <= 0.0 && Sm > 0.0) {
      
      for (std::size_t idir=0; idir<3; ++idir)
        flx[momentumIdx(nmat, idir)] =
          uStar[0][ncomp+velocityIdx(nmat, idir)] * Sm + pStar*fn[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = uStar[0][volfracIdx(nmat, k)] * Sm;
        flx[densityIdx(nmat, k)] = uStar[0][densityIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] = (uStar[0][energyIdx(nmat, k)] + apStar[k]) * Sm;
      }

      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back(apStar[k]);
      flx.push_back(Sm);
    }
    else if (Sm <= 0.0 && Sr >= 0.0) {
      
      for (std::size_t idir=0; idir<3; ++idir)
        flx[momentumIdx(nmat, idir)] =
          uStar[1][ncomp+velocityIdx(nmat, idir)] * Sm + pStar*fn[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = uStar[1][volfracIdx(nmat, k)] * Sm;
        flx[densityIdx(nmat, k)] = uStar[1][densityIdx(nmat, k)] * Sm;
        flx[energyIdx(nmat, k)] = (uStar[1][energyIdx(nmat, k)] + apStar[k]) * Sm;
      }

      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back(apStar[k]);
      flx.push_back(Sm);
    }
    else {
      
      for (std::size_t idir=0; idir<3; ++idir)
        flx[momentumIdx(nmat, idir)] =
          u[1][ncomp+velocityIdx(nmat, idir)] * vnr + pr*fn[idir];

      for (std::size_t k=0; k<nmat; ++k) {
        flx[volfracIdx(nmat, k)] = u[1][volfracIdx(nmat, k)] * vnr;
        flx[densityIdx(nmat, k)] = u[1][densityIdx(nmat, k)] * vnr;
        flx[energyIdx(nmat, k)] = (u[1][energyIdx(nmat, k)] + apr[k]) * vnr;
      }

      for (std::size_t k=0; k<nmat; ++k)
        flx.push_back(apr[k]);
      flx.push_back(vnr);
    }

    return flx;
  }

  ////! Flux type accessor
  ////! \return Flux type
  //static ctr::FluxType type() noexcept {
  //  return ctr::FluxType::HLLCSolids; }
};

} // inciter::

#endif // HLLCSolids_h
