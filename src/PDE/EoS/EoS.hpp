// *****************************************************************************
/*!
  \file      src/PDE/EoS/EoS.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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

//! Get the ratio of specific heats (gamma) for a material
//! \tparam Eq Equation type to operate on, e.g., tag::compflow, tag::multimat
//! \param[in] system Equation system index
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified by
//!   the calling code
//! \return Material ratio of specific heats (gamma)
template< class Eq >
tk::real gamma( ncomp_t system,
  std::size_t imat=0 )
{
  const auto& matprop =
    g_inputdeck.get< tag::param, Eq, tag::material >()[system];
  const auto& meos =
    g_inputdeck.get< tag::param, Eq, tag::matidxmap >().template get<
    tag::eosidx >()[imat];
  const auto& midx =
    g_inputdeck.get< tag::param, Eq, tag::matidxmap >().template get<
    tag::matidx >()[imat];

  return matprop[meos].template get< tag::gamma >()[midx];
}

//! Get the specific heat at constant volume (cv) for a material
//! \tparam Eq Equation type to operate on, e.g., tag::compflow, tag::multimat
//! \param[in] system Equation system index
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified by
//!   the calling code
//! \return Material specific heat at constant volume (cv)
template< class Eq >
tk::real cv( ncomp_t system,
  std::size_t imat=0 )
{
  const auto& matprop =
    g_inputdeck.get< tag::param, Eq, tag::material >()[system];
  const auto& meos =
    g_inputdeck.get< tag::param, Eq, tag::matidxmap >().template get<
    tag::eosidx >()[imat];
  const auto& midx =
    g_inputdeck.get< tag::param, Eq, tag::matidxmap >().template get<
    tag::matidx >()[imat];

  return matprop[meos].template get< tag::cv >()[midx];
}

//! Get the stiffness parameter (pstiff) for a material
//! \tparam Eq Equation type to operate on, e.g., tag::compflow, tag::multimat
//! \param[in] system Equation system index
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified by
//!   the calling code
//! \return Material stiffness parameter (pstiff)
template< class Eq >
tk::real pstiff( ncomp_t system,
  std::size_t imat=0 )
{
  const auto& matprop =
    g_inputdeck.get< tag::param, Eq, tag::material >()[system];
  const auto& meos =
    g_inputdeck.get< tag::param, Eq, tag::matidxmap >().template get<
    tag::eosidx >()[imat];
  const auto& midx =
    g_inputdeck.get< tag::param, Eq, tag::matidxmap >().template get<
    tag::matidx >()[imat];

  return matprop[meos].template get< tag::pstiff >()[midx];
}

//! Get the thermal conductivity (k) for a material
//! \tparam Eq Equation type to operate on, e.g., tag::compflow, tag::multimat
//! \param[in] system Equation system index
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified by
//!   the calling code
//! \return Material thermal conductivity (k)
template< class Eq >
tk::real k( ncomp_t system,
  std::size_t imat=0 )
{
  const auto& matprop =
    g_inputdeck.get< tag::param, Eq, tag::material >()[system];
  const auto& meos =
    g_inputdeck.get< tag::param, Eq, tag::matidxmap >().template get<
    tag::eosidx >()[imat];
  const auto& midx =
    g_inputdeck.get< tag::param, Eq, tag::matidxmap >().template get<
    tag::matidx >()[imat];

  return matprop[meos].template get< tag::k >()[midx];
}

//! Get the dynamic viscosity (mu) for a material
//! \tparam Eq Equation type to operate on, e.g., tag::compflow, tag::multimat
//! \param[in] system Equation system index
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified by
//!   the calling code
//! \return Material dynamic viscosity (mu)
template< class Eq >
tk::real mu( ncomp_t system,
  std::size_t imat=0 )
{
  const auto& matprop =
    g_inputdeck.get< tag::param, Eq, tag::material >()[system];
  const auto& meos =
    g_inputdeck.get< tag::param, Eq, tag::matidxmap >().template get<
    tag::eosidx >()[imat];
  const auto& midx =
    g_inputdeck.get< tag::param, Eq, tag::matidxmap >().template get<
    tag::matidx >()[imat];

  return matprop[meos].template get< tag::mu >()[midx];
}

//! \brief Calculate density from the material pressure and temperature using
//!   the stiffened-gas equation of state
//! \tparam Eq Equation type to operate on, e.g., tag::compflow, tag::multimat
//! \param[in] system Equation system index
//! \param[in] pr Material pressure
//! \param[in] temp Material temperature
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified by
//!   the calling code
//! \return Material density calculated using the stiffened-gas EoS
template< class Eq >
tk::real eos_density( ncomp_t system,
                      tk::real pr,
                      tk::real temp,
                      std::size_t imat=0 )
{
  // query input deck to get gamma, p_c, cv
  auto g = gamma< Eq >(system, imat);
  auto p_c = pstiff< Eq >(system, imat);
  auto c_v = cv< Eq >(system, imat);

  tk::real rho = (pr + p_c) / ((g-1.0) * c_v * temp);
  return rho;
}

//! \brief Calculate pressure from the material density, momentum and total
//!   energy using the stiffened-gas equation of state
//! \tparam Eq Equation type to operate on, e.g., tag::compflow, tag::multimat
//! \param[in] system Equation system index
//! \param[in] arho Material partial density (alpha_k * rho_k)
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] arhoE Material total energy (alpha_k * rho_k * E_k)
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for the
//!   single-material system, this argument can be left unspecified by the
//!   calling code
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified by
//!   the calling code
//! \return Material partial pressure (alpha_k * p_k) calculated using the
//!   stiffened-gas EoS
#pragma omp declare simd
template< class Eq >
tk::real eos_pressure( ncomp_t system,
                       tk::real arho,
                       tk::real u,
                       tk::real v,
                       tk::real w,
                       tk::real arhoE,
                       tk::real alpha=1.0,
                       std::size_t imat=0 )
{
  // query input deck to get gamma, p_c
  auto g = gamma< Eq >(system, imat);
  auto p_c = pstiff< Eq >(system, imat);

  tk::real partpressure = (arhoE - 0.5 * arho * (u*u + v*v + w*w) - alpha*p_c)
                          * (g-1.0) - alpha*p_c;
  return partpressure;
}

//! Calculate speed of sound from the material density and material pressure
//! \tparam Eq Equation type to operate on, e.g., tag::compflow, tag::multimat
//! \param[in] system Equation system index
//! \param[in] arho Material partial density (alpha_k * rho_k)
//! \param[in] apr Material partial pressure (alpha_k * p_k)
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for the
//!   single-material system, this argument can be left unspecified by the
//!   calling code
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified by
//!   the calling code
//! \return Material speed of sound using the stiffened-gas EoS
template< class Eq >
tk::real eos_soundspeed( ncomp_t system,
                         tk::real arho, tk::real apr,
                         tk::real alpha=1.0, std::size_t imat=0 )
{
  // query input deck to get gamma, p_c
  auto g = gamma< Eq >(system, imat);
  auto p_c = pstiff< Eq >(system, imat);

  auto p_eff = std::max( 1.0e-15, apr+(alpha*p_c) );

  tk::real a = std::sqrt( g * p_eff / arho );
  return a;
}

//! \brief Calculate material specific total energy from the material density,
//!   momentum and material pressure
//! \tparam Eq Equation type to operate on, e.g., tag::compflow, tag::multimat
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
template< class Eq >
tk::real eos_totalenergy( ncomp_t system,
                          tk::real rho,
                          tk::real u,
                          tk::real v,
                          tk::real w,
                          tk::real pr,
                          std::size_t imat=0 )
{
  // query input deck to get gamma, p_c
  auto g = gamma< Eq >(system, imat);
  auto p_c = pstiff< Eq >(system, imat);

  tk::real rhoE = (pr + p_c) / (g-1.0) + 0.5 * rho * (u*u + v*v + w*w) + p_c;
  return rhoE;
}

//! \brief Calculate material temperature from the material density, and
//!   material specific total energy
//! \tparam Eq Equation type to operate on, e.g., tag::compflow, tag::multimat
//! \param[in] system Equation system index
//! \param[in] arho Material partial density (alpha_k * rho_k)
//! \param[in] u X-velocity
//! \param[in] v Y-velocity
//! \param[in] w Z-velocity
//! \param[in] arhoE Material total energy (alpha_k * rho_k * E_k)
//! \param[in] alpha Material volume fraction. Default is 1.0, so that for the
//!   single-material system, this argument can be left unspecified by the
//!   calling code
//! \param[in] imat Material-id who's EoS is required. Default is 0, so that
//!   for the single-material system, this argument can be left unspecified by
//!   the calling code
//! \return Material temperature using the stiffened-gas EoS
template< class Eq >
tk::real eos_temperature( ncomp_t system,
                          tk::real arho,
                          tk::real u,
                          tk::real v,
                          tk::real w,
                          tk::real arhoE,
                          tk::real alpha=1.0,
                          std::size_t imat=0 )
{
  // query input deck to get p_c, cv
  auto c_v = cv< Eq >(system, imat);
  auto p_c = pstiff< Eq >(system, imat);

  tk::real t = (arhoE - 0.5 * arho * (u*u + v*v + w*w) - alpha*p_c) / (arho*c_v);
  return t;
}

} //inciter::

#endif // EoS_h
