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
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

using ncomp_t = kw::ncomp::info::expect::type;

  //! \brief Calculate pressure from the material density, momentum and total
  //!   energy using the stiffened-gas equation of state
  template< class T > tk::real eos_pressure( ncomp_t system,
                         tk::real rho,
                         tk::real u,
                         tk::real v,
                         tk::real w,
                         tk::real rhoE,
                         std::size_t imat=0 )
  // ***************************************************************************
  //! \brief Calculate pressure from the material density, momentum and total
  //!   energy using the stiffened-gas equation of state
  //! \param[in] system Equation system index
  //! \param[in] rho Material density
  //! \param[in] u X-velocity
  //! \param[in] v Y-velocity
  //! \param[in] w Z-velocity
  //! \param[in] rhoE Material total energy
  //! \param[in] imat Material-id who's EoS is required. Default is 0, so that
  //!   for the single-material system, this argument can be left unspecified by
  //!   the calling code
  //! \return Material pressure calculated using the stiffened-gas EoS
  // ***************************************************************************
  {
    // query input deck to get gamma, p_c
    auto g =
      g_inputdeck.get< tag::param, T, tag::gamma >()[ system ][imat];
    auto p_c =
      g_inputdeck.get< tag::param, T, tag::pstiff >()[ system ][imat];

    tk::real pressure = (rhoE - 0.5 * rho * (u*u + v*v + w*w) - p_c)
                        * (g-1.0) - p_c;
    return pressure;
  }

  //! Calculate speed of sound from the material density and material pressure
  template< class T > tk::real eos_soundspeed( ncomp_t system,
                                               tk::real rho, tk::real pr,
                                               std::size_t imat=0 )
  // ***************************************************************************
  //  Calculate speed of sound from the material density and material pressure
  //! \param[in] system Equation system index
  //! \param[in] rho Material density
  //! \param[in] pr Material pressure
  //! \param[in] imat Material-id who's EoS is required. Default is 0, so that
  //!   for the single-material system, this argument can be left unspecified by
  //!   the calling code
  //! \return Material speed of sound using the stiffened-gas EoS
  // ***************************************************************************
  {
    // query input deck to get gamma, p_c
    auto g =
      g_inputdeck.get< tag::param, T, tag::gamma >()[ system ][imat];
    auto p_c =
      g_inputdeck.get< tag::param, T, tag::pstiff >()[ system ][imat];

    tk::real a = std::sqrt( g * (pr+p_c) / rho );
    return a;
  }

  //! \brief Calculate material specific total energy from the material density,
  //!   momentum and material pressure
  template< class T > tk::real eos_totalenergy( ncomp_t system,
                            tk::real rho,
                            tk::real u,
                            tk::real v,
                            tk::real w,
                            tk::real pr,
                            std::size_t imat=0 )
  // ***************************************************************************
  //  \brief Calculate material specific total energy from the material density,
  //    momentum and material pressure using the stiffened-gas equation of state
  //! \param[in] system Equation system index
  //! \param[in] rho Material density
  //! \param[in] u X-velocity
  //! \param[in] v Y-velocity
  //! \param[in] w Z-velocity
  //! \param[in] pr Material pressure
  //! \param[in] imat Material-id who's EoS is required. Default is 0, so that
  //!   for the single-material system, this argument can be left unspecified by
  //!   the calling code
  //! \return Material specific total energy using the stiffened-gas EoS
  // ***************************************************************************
  {
    // query input deck to get gamma, p_c
    auto g =
      g_inputdeck.get< tag::param, T, tag::gamma >()[ system ][imat];
    auto p_c =
      g_inputdeck.get< tag::param, T, tag::pstiff >()[ system ][imat];

    tk::real rhoE = (pr + p_c) / (g-1.0) + 0.5 * rho * (u*u + v*v + w*w) + p_c;
    return rhoE;
  }

} //inciter::

#endif // EoS_h
