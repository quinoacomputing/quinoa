// *****************************************************************************
/*!
  \file      src/PDE/Riemann/Rusanov.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Rusanov Riemann flux function
  \details   This file implements the Rusanov Riemann solver.
*/
// *****************************************************************************
#ifndef Rusanov_h
#define Rusanov_h

#include <vector>

#include "Types.hpp"
#include "Fields.hpp"
#include "Tags.hpp"
#include "FunctionPrototypes.hpp"
#include "Inciter/Options/Flux.hpp"
#include "EoS/EoS.hpp"

namespace inciter {

//! Rusanov approximate Riemann solver
//! \details This class can be used polymorphically with inciter::RiemannSolver
struct Rusanov {

  //! Rusanov approximate Riemann solver flux function
  //! \param[in] fn Face/Surface normal
  //! \param[in] u Left and right unknown/state vector
  //! \return Riemann solution according to Rusanov
  //! \note The function signature must follow tk::RiemannFluxFn
  static tk::RiemannFluxFn::result_type
  flux( const std::array< tk::real, 3 >& fn,
        const std::array< std::vector< tk::real >, 2 >& u,
        const std::array< std::vector< tk::real >, 2 >& = {} )
  {
    std::vector< tk::real > flx( u[0].size(), 0 );

    const auto & UL = u[0];
    const auto & UR = u[1];

    // Primitive variables
    auto rhol = UL[0];
    auto rhor = UR[0];

    auto ul = UL[1]/rhol;
    auto vl = UL[2]/rhol;
    auto wl = UL[3]/rhol;

    auto ur = UR[1]/rhor;
    auto vr = UR[2]/rhor;
    auto wr = UR[3]/rhor;

    auto pl = eos_pressure< tag::compflow >( 0, rhol, ul, vl, wl, UL[4] );
    auto pr = eos_pressure< tag::compflow >( 0, rhor, ur, vr, wr, UR[4] );

    auto al = eos_soundspeed< tag::compflow >( 0, rhol, pl );
    auto ar = eos_soundspeed< tag::compflow >( 0, rhor, pr );

    // Face-normal velocities
    tk::real vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];
    tk::real vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];


    // Numerical fluxes
    flx[0]  = 0.5 * ( UL[0] * vnl );
    flx[1]  = 0.5 * ( UL[1] * vnl + pl*fn[0] );
    flx[2]  = 0.5 * ( UL[2] * vnl + pl*fn[1] );
    flx[3]  = 0.5 * ( UL[3] * vnl + pl*fn[2] );
    flx[4]  = 0.5 * ( ( UL[4] + pl ) * vnl );
    
    flx[0] += 0.5 * ( UR[0] * vnr );
    flx[1] += 0.5 * ( UR[1] * vnr + pr*fn[0] );
    flx[2] += 0.5 * ( UR[2] * vnr + pr*fn[1] );
    flx[3] += 0.5 * ( UR[3] * vnr + pr*fn[2] );
    flx[4] += 0.5 * ( ( UR[4] + pr ) * vnr );
    
    auto sl = std::abs(vnl) + al;
    auto sr = std::abs(vnr) + ar;
    auto smax = std::max( sl, sr );

    flx[0] -= 0.5 * smax * ( UR[0] - UL[0] ); 
    flx[1] -= 0.5 * smax * ( UR[1] - UL[1] ); 
    flx[2] -= 0.5 * smax * ( UR[2] - UL[2] ); 
    flx[3] -= 0.5 * smax * ( UR[3] - UL[3] ); 
    flx[4] -= 0.5 * smax * ( UR[4] - UL[4] ); 

    return flx;
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::Rusanov; }
};

} // inciter::

#endif // Rusanov_h
