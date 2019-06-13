// *****************************************************************************
/*!
  \file      src/PDE/EoS/EoS.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Equation of state class
  \details   This file defines functions for equations of state for the
    compressible flow equations.
*/
// *****************************************************************************
#ifndef EoS_h
#define EoS_h

#include "Data.hpp"

namespace inciter {

using ncomp_t = kw::ncomp::info::expect::type;

  //! \brief Calculate pressure from the material density, momentum and total energy
  //!   using the stiffened-gas equation of state
  tk::real eos_pressure( ncomp_t system,
                         tk::real rho,
                         tk::real rhou,
                         tk::real rhov,
                         tk::real rhow,
                         tk::real rhoE,
                         int imat=0 );

  //! Calculate speed of sound from the material density and material pressure
  tk::real eos_soundspeed( ncomp_t system, tk::real rho, tk::real pr,
                           int imat=0 );


  //! \brief Calculate material specific total energy from the material density,
  //!   momentum and material pressure
  tk::real eos_totalenergy( ncomp_t system,
                            tk::real rho,
                            tk::real rhou,
                            tk::real rhov,
                            tk::real rhow,
                            tk::real pr,
                            int imat=0 );

} //inciter::

#endif // EoS_h
