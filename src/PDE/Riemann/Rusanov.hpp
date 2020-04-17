// *****************************************************************************
/*!
  \file      src/PDE/Riemann/Rusanov.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Rusanov Riemann flux function
  \details   This file implements the Rusanov Riemann solver, specific to ALECG.
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
  //! \param[in] uL Left unknown/state vector
  //! \param[in] uR right unknown/state vector
  //! \param[in] aux Auxiliary vector, used here to pass in normal vectors
  //!    weighted by the number of contributions to the edge
  //! \return Riemann solution according to Rusanov
  static std::array< tk::real, 5 >
  flux( const std::array< tk::real, 3 >& fn,
        const std::vector< tk::real >& uL,
        const std::vector< tk::real >& uR,
        const std::array< tk::real, 3 > & aux )
  {
    Assert( uL.size() == 5 && uR.size() == 5, "Size mismatch" );

    std::array< tk::real, 5 > flx;

    // Primitive variables
    auto rhol = uL[0];
    auto rhor = uR[0];

    auto ul = uL[1]/rhol;
    auto vl = uL[2]/rhol;
    auto wl = uL[3]/rhol;

    auto ur = uR[1]/rhor;
    auto vr = uR[2]/rhor;
    auto wr = uR[3]/rhor;

    auto pl = eos_pressure< tag::compflow >( 0, rhol, ul, vl, wl, uL[4] );
    auto pr = eos_pressure< tag::compflow >( 0, rhor, ur, vr, wr, uR[4] );

    auto al = eos_soundspeed< tag::compflow >( 0, rhol, pl );
    auto ar = eos_soundspeed< tag::compflow >( 0, rhor, pr );

    // Face-normal velocities
    tk::real vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];
    tk::real vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];

    // Numerical fluxes
    flx[0]  = 0.5 * ( uL[0] * vnl );
    flx[1]  = 0.5 * ( uL[1] * vnl + pl*fn[0] );
    flx[2]  = 0.5 * ( uL[2] * vnl + pl*fn[1] );
    flx[3]  = 0.5 * ( uL[3] * vnl + pl*fn[2] );
    flx[4]  = 0.5 * ( ( uL[4] + pl ) * vnl );
    
    flx[0] += 0.5 * ( uR[0] * vnr );
    flx[1] += 0.5 * ( uR[1] * vnr + pr*fn[0] );
    flx[2] += 0.5 * ( uR[2] * vnr + pr*fn[1] );
    flx[3] += 0.5 * ( uR[3] * vnr + pr*fn[2] );
    flx[4] += 0.5 * ( ( uR[4] + pr ) * vnr );
    
    // dissipation term
    const auto& n2 = aux;
    auto len = tk::length(n2);
    vnl = ( ul*n2[0] + vl*n2[1] + wl*n2[2] );
    vnr = ( ur*n2[0] + vr*n2[1] + wr*n2[2] );
    auto sl = std::abs(vnl) + al*len;
    auto sr = std::abs(vnr) + ar*len;
    auto smax = std::max( sl, sr );

    flx[0] -= 0.5 * smax * ( uR[0] - uL[0] );
    flx[1] -= 0.5 * smax * ( uR[1] - uL[1] );
    flx[2] -= 0.5 * smax * ( uR[2] - uL[2] );
    flx[3] -= 0.5 * smax * ( uR[3] - uL[3] );
    flx[4] -= 0.5 * smax * ( uR[4] - uL[4] );

    return flx;
  }

  //! Flux type accessor
  //! \return Flux type
  static ctr::FluxType type() noexcept { return ctr::FluxType::Rusanov; }
};

} // inciter::

#endif // Rusanov_h
