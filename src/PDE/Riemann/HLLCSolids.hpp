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

//! HLLC approximate Riemann solver for solids (SHOULD ONLY WORK FOR SINGLE
//! MATERIAL FLUIDS CURRENTLY)
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
    auto rhol = u[0][1]/u[0][0];
    auto rhor = u[1][1]/u[1][0];

    auto ul = u[0][2]/rhol;
    auto vl = u[0][3]/rhol;
    auto wl = u[0][4]/rhol;

    auto ur = u[1][2]/rhor;
    auto vr = u[1][3]/rhor;
    auto wr = u[1][4]/rhor;

    auto pl = mat_blk[0].compute< EOS::pressure >( rhol, ul, vl, wl,
      u[0][5] );
    auto pr = mat_blk[0].compute< EOS::pressure >( rhor, ur, vr, wr,
      u[1][5] );

    auto al = mat_blk[0].compute< EOS::soundspeed >( rhol, pl );
    auto ar = mat_blk[0].compute< EOS::soundspeed >( rhor, pr );

    // Face-normal velocities
    tk::real vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];
    tk::real vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];

    // Roe-averaged variables
    auto rlr = sqrt(rhor/rhol);
    auto rlr1 = 1.0 + rlr;

    auto vnroe = (vnr*rlr + vnl)/rlr1 ;
    auto aroe = (ar*rlr + al)/rlr1 ;

    // Signal velocities
    auto Sl = fmin(vnl-al, vnroe-aroe);
    auto Sr = fmax(vnr+ar, vnroe+aroe);
    auto Sm = ( rhor*vnr*(Sr-vnr) - rhol*vnl*(Sl-vnl) + pl-pr )
             /( rhor*(Sr-vnr) - rhol*(Sl-vnl) );

    // Middle-zone (star) variables
    auto pStar = rhol*(vnl-Sl)*(vnl-Sm) + pl;
    auto uStar = u;

    uStar[0][0] = u[0][0];
    uStar[0][1] = (Sl-vnl) * rhol / (Sl-Sm);
    uStar[0][2] = ((Sl-vnl) * u[0][2] + (pStar-pl)*fn[0]) / (Sl-Sm);
    uStar[0][3] = ((Sl-vnl) * u[0][3] + (pStar-pl)*fn[1]) / (Sl-Sm);
    uStar[0][4] = ((Sl-vnl) * u[0][4] + (pStar-pl)*fn[2]) / (Sl-Sm);
    uStar[0][5] = ((Sl-vnl) * u[0][5] - pl*vnl + pStar*Sm) / (Sl-Sm);

    uStar[1][0] = u[1][0];
    uStar[1][1] = (Sr-vnr) * rhor / (Sr-Sm);
    uStar[1][2] = ((Sr-vnr) * u[1][2] + (pStar-pr)*fn[0]) / (Sr-Sm);
    uStar[1][3] = ((Sr-vnr) * u[1][3] + (pStar-pr)*fn[1]) / (Sr-Sm);
    uStar[1][4] = ((Sr-vnr) * u[1][4] + (pStar-pr)*fn[2]) / (Sr-Sm);
    uStar[1][5] = ((Sr-vnr) * u[1][5] - pr*vnr + pStar*Sm) / (Sr-Sm);

    // Numerical fluxes
    if (Sl > 0.0) {
      flx[0] = u[0][0] * vnl;
      flx[1] = u[0][1] * vnl;
      flx[2] = u[0][2] * vnl + pl*fn[0];
      flx[3] = u[0][3] * vnl + pl*fn[1];
      flx[4] = u[0][4] * vnl + pl*fn[2];
      flx[5] = ( u[0][5] + pl ) * vnl;
      flx.push_back(pl);
      flx.push_back(vnl);
    }
    else if (Sl <= 0.0 && Sm > 0.0) {
      flx[0] = uStar[0][0] * Sm;
      flx[1] = uStar[0][1] * Sm;
      flx[2] = uStar[0][2] * Sm + pStar*fn[0];
      flx[3] = uStar[0][3] * Sm + pStar*fn[1];
      flx[4] = uStar[0][4] * Sm + pStar*fn[2];
      flx[5] = ( uStar[0][5] + pStar ) * Sm;
      flx.push_back(pStar);
      flx.push_back(Sm);
    }
    else if (Sm <= 0.0 && Sr >= 0.0) {
      flx[0] = uStar[1][0] * Sm;
      flx[1] = uStar[1][1] * Sm;
      flx[2] = uStar[1][2] * Sm + pStar*fn[0];
      flx[3] = uStar[1][3] * Sm + pStar*fn[1];
      flx[4] = uStar[1][4] * Sm + pStar*fn[2];
      flx[5] = ( uStar[1][5] + pStar ) * Sm;
      flx.push_back(pStar);
      flx.push_back(Sm);
    }
    else {
      flx[0] = u[1][0] * vnr;
      flx[1] = u[1][1] * vnr;
      flx[2] = u[1][2] * vnr + pr*fn[0];
      flx[3] = u[1][3] * vnr + pr*fn[1];
      flx[4] = u[1][4] * vnr + pr*fn[2];
      flx[5] = ( u[1][5] + pr ) * vnr;
      flx.push_back(pr);
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
