// *****************************************************************************
/*!
  \file      src/PDE/EoS/GodunovRomenski.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Godunov-Romenski equation of state for solids
  \details   This file declares functions for the Godunov-Romenski equation of
             state for solids and a hydro EoS for aluminum. These functions were
             taken from Example 1 of Barton, Philip T. "An interface-capturing
             Godunov method for the simulation of compressible solid-fluid
             problems." Journal of Computational Physics 390 (2019): 25-50.
*/
// *****************************************************************************
#ifndef GodunovRomenski_h
#define GodunovRomenski_h

#include "Data.hpp"

namespace inciter {

class GodunovRomenski {

  private:
    tk::real m_gamma, m_mu, m_rho0, m_alpha, m_K0;

    //! \brief Calculate elastic contribution to material energy from the
    //!   material density, and deformation gradient tensor
    tk::real elasticEnergy(
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
      std::array< std::array< tk::real, 3 >, 3 >& devH ) const;

    //! \brief Calculate cold-compression contribution to material energy from
    //!   the material density
    tk::real coldcomprEnergy( tk::real rho ) const;

    //! Calculate the derivative of the cold compression pressure wrt. density
    tk::real DpccDrho( tk::real rho ) const;

  public:
    //! Default constructor
    GodunovRomenski() = default;

    //! Constructor
    GodunovRomenski(
      tk::real gamma,
      tk::real mu,
      tk::real rho0,
      tk::real alpha,
      tk::real K0 );

    //! Set rho0 EOS parameter; no-op since provided by user
    void setRho0(tk::real) {}

    //! Calculate density from the material pressure and temperature
    tk::real density( tk::real pr,
                      tk::real temp ) const;

    //! Calculate pressure from the material density, momentum and total energy
    tk::real pressure(
      tk::real arho,
      tk::real u,
      tk::real v,
      tk::real w,
      tk::real arhoE,
      tk::real alpha=1.0,
      std::size_t imat=0,
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}} ) const;

    //! \brief Calculate cold-compression contribution to material pressure from
    //!   the material density
    tk::real pressure_coldcompr(
      tk::real arho,
      tk::real alpha=1.0 ) const;

    //! \brief Calculate the elastic Cauchy stress tensor from the material
    //!   density, momentum, total energy, and inverse deformation gradient
    //!   tensor using the GodunovRomenski equation of state
    std::array< std::array< tk::real, 3 >, 3 >
    CauchyStress(
      tk::real arho,
      tk::real,
      tk::real,
      tk::real,
      tk::real,
      tk::real alpha,
      std::size_t /*imat*/,
      const std::array< std::array< tk::real, 3 >, 3 >& adefgrad ) const;

    //! Calculate speed of sound from the material density and material pressure
    tk::real soundspeed(
      tk::real arho,
      tk::real apr,
      tk::real alpha=1.0,
      std::size_t imat=0,
      const std::array< std::array< tk::real, 3 >, 3 >& adefgrad={{}} ) const;

    //! Calculate speed of shear waves
    tk::real shearspeed(
      tk::real arho,
      tk::real alpha=1.0,
      std::size_t imat=0 ) const;

    //! \brief Calculate material specific total energy from the material
    //!   density, momentum and material pressure
    tk::real totalenergy(
      tk::real rho,
      tk::real u,
      tk::real v,
      tk::real w,
      tk::real pr,
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}} ) const;

    //! \brief Calculate material temperature from the material density, and
    //!   material specific total energy
    tk::real temperature(
      tk::real arho,
      tk::real u,
      tk::real v,
      tk::real w,
      tk::real arhoE,
      tk::real alpha=1.0,
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}} ) const;

    //! Compute the minimum allowed pressure
    tk::real min_eff_pressure(
      tk::real min,
      tk::real arho,
      tk::real alpha=1.0 ) const;

    //! Compute the reference density
    tk::real refDensity() const { return density(refPressure(), 300.0); }

    //! Compute the reference pressure
    tk::real refPressure() const { return 1.0e5; }

    //! Return initial density
    tk::real rho0() const { return m_rho0; }

    //! Return gas constant (no-op)
    tk::real gas_constant() { return 0.0; }

    //! Return internal energy (no-op)
    tk::real calc_e( [[maybe_unused]] tk::real temp) const { return 0.0; }

    //! Return specific heat (no-op)
    tk::real calc_cv( [[maybe_unused]] tk::real temp) const { return 0.0; }

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) /*override*/ {
      p | m_gamma;
      p | m_mu;
      p | m_rho0;
      p | m_alpha;
      p | m_K0;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i GodunovRomenski object reference
    friend void operator|( PUP::er& p, GodunovRomenski& i ) { i.pup(p); }
    //@}
};

} //inciter::

#endif // GodunovRomenski_h
