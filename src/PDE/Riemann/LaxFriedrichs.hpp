// *****************************************************************************
/*!
  \file      src/PDE/Riemann/LaxFriedrichs.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Lax-Friedrichs Riemann flux function
  \details   This file implements the Lax-Friedrichs Riemann solver.
*/
// *****************************************************************************
#ifndef LaxFriedrichs_h
#define LaxFriedrichs_h

#include <vector>

#include "Types.hpp"
#include "Fields.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EoS.hpp"
#include "EoS/EosVariant.hpp"

namespace inciter {

//! Lax-Friedrichs approximate Riemann solver
struct LaxFriedrichs {

  //! Lax-Friedrichs approximate Riemann solver flux function
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann solution according Lax and Friedrichs
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::vector< EOS >& mat_blk,
        const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::vector< std::array< tk::real, 3 > >& = {} )
  {
    std::vector< tk::real >  flx( u[0].size(), 0.0 ),
                            fluxl( u[0].size(), 0.0 ),
                            fluxr( u[0].size(), 0.0 );

    // Primitive variables
    auto rhol = u[0][0];
    auto rhor = u[1][0];

    auto ul = u[0][1]/rhol;
    auto vl = u[0][2]/rhol;
    auto wl = u[0][3]/rhol;

    auto ur = u[1][1]/rhor;
    auto vr = u[1][2]/rhor;
    auto wr = u[1][3]/rhor;

    auto pl = mat_blk[0].eosCall< EOS::eos_pressure >( rhol, ul, vl, wl,
      u[0][4] );
    auto pr = mat_blk[0].eosCall< EOS::eos_pressure >( rhor, ur, vr, wr,
      u[1][4] );

    auto al = mat_blk[0].eosCall< EOS::eos_soundspeed >( rhol, pl );
    auto ar = mat_blk[0].eosCall< EOS::eos_soundspeed >( rhor, pr );

    // Face-normal velocities
    auto vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];
    auto vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];

    // Flux functions
    fluxl[0] = u[0][0] * vnl;
    fluxl[1] = u[0][1] * vnl + pl*fn[0];
    fluxl[2] = u[0][2] * vnl + pl*fn[1];
    fluxl[3] = u[0][3] * vnl + pl*fn[2];
    fluxl[4] = ( u[0][4] + pl ) * vnl;

    fluxr[0] = u[1][0] * vnr;
    fluxr[1] = u[1][1] * vnr + pr*fn[0];
    fluxr[2] = u[1][2] * vnr + pr*fn[1];
    fluxr[3] = u[1][3] * vnr + pr*fn[2];
    fluxr[4] = ( u[1][4] + pr ) * vnr;

    auto lambda = std::max(al,ar) + std::max( std::abs(vnl), std::abs(vnr) );

    // Numerical flux function
    for(std::size_t c=0; c<5; ++c)
      flx[c] = 0.5 * (fluxl[c] + fluxr[c] - lambda*(u[1][c] - u[0][c]));

    return flx;
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::LaxFriedrichs; }
};

} // inciter::

#endif // LaxFriedrichs_h
