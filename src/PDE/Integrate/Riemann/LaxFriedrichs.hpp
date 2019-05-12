// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Riemann/LaxFriedrichs.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
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
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Inciter/Options/Flux.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! Lax-Friedrichs approximate Riemann solver
//! \details This class can be used polymorphically with inciter::RiemannSolver
struct LaxFriedrichs {

  //! Lax-Friedrichs approximate Riemann solver flux function
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann solution according Lax and Friedrichs
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::vector< std::array< tk::real, 3 > >& )
  {
    std::vector< tk::real >  flx( u[0].size(), 0.0 ),
                            fluxl( u[0].size(), 0.0 ),
                            fluxr( u[0].size(), 0.0 );

    // ratio of specific heats
    auto g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0];

    // Primitive variables
    auto rhol = u[0][0];
    auto rhor = u[1][0];

    auto pl = (g-1.0)*(u[0][4] - (u[0][1]*u[0][1] +
                                  u[0][2]*u[0][2] +
                                  u[0][3]*u[0][3]) / (2.0*rhol));

    auto pr = (g-1.0)*(u[1][4] - (u[1][1]*u[1][1] +
                                  u[1][2]*u[1][2] +
                                  u[1][3]*u[1][3]) / (2.0*rhor));

    auto al = sqrt(g * pl / rhol);
    auto ar = sqrt(g * pr / rhor);

    // Face-normal velocities
    auto ul = u[0][1]/rhol;
    auto vl = u[0][2]/rhol;
    auto wl = u[0][3]/rhol;

    tk::real vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];

    auto ur = u[1][1]/rhor;
    auto vr = u[1][2]/rhor;
    auto wr = u[1][3]/rhor;

    tk::real vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];

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

    auto lambda = fmax(al,ar) + fmax(fabs(vnl),fabs(vnr));

    // Numerical flux function
    for(std::size_t c=0; c<5; ++c)
    {
      flx[c] = 0.5 * ( fluxl[c] + fluxr[c]
                       - lambda * (u[1][c] - u[0][c]) );
    }

    return flx;
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::LaxFriedrichs; }
};

} // inciter::

#endif // LaxFriedrichs_h
