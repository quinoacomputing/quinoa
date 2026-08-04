// *****************************************************************************
/*!
  \file      src/PDE/EoS/LinearMieGruneisen.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Linear Mie-Gruneisen equation of state (a.k.a. shock-wave EOS)
  \details   This file declares functions for the LinearMieGruneisen equation of
             state for solids. The elastic contribution follows SmallShearSolid.
             The hydrodynamic contribution uses a linear Us-Up Hugoniot and a
             density-dependent Gruneisen coefficient. Ref. Shyue, K. M. (2001).
             A fluid-mixture type algorithm for compressible multicomponent flow
             with Mie–Grüneisen equation of state. J Comp Phys, 171(2), 678-707.
*/
// *****************************************************************************
#ifndef LinearMieGruneisen_h
#define LinearMieGruneisen_h

#include "Data.hpp"

namespace inciter {

class LinearMieGruneisen {

  private:
    tk::real m_gamma0, m_rho0, m_alpha, m_c0, m_s1, m_cv, m_mu;

    //! \brief Calculate elastic contribution to material energy from the
    //!   material density, and deformation gradient tensor
    tk::real elasticEnergy(
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
      tk::real& eps2 ) const;

    //! Density-dependent Gruneisen coefficient
    tk::real gruneisen( tk::real rho ) const;

    //! Derivative of density-dependent Gruneisen coefficient wrt density
    tk::real dGruneisenDrho( tk::real rho ) const;

    //! Compression measure based on reference density
    tk::real eta( tk::real rho ) const;

    //! Derivative of compression measure wrt density
    tk::real dEtaDrho( tk::real rho ) const;

    //! Linear Us-Up Hugoniot pressure
    tk::real hugoniotPressure( tk::real rho ) const;

    //! Derivative of Hugoniot pressure wrt density
    tk::real dHugoniotPressureDrho( tk::real rho ) const;

    //! Hugoniot internal energy
    tk::real hugoniotEnergy( tk::real rho ) const;

    //! Derivative of Hugoniot internal energy wrt density
    tk::real dHugoniotEnergyDrho( tk::real rho ) const;

  public:
    //! Default constructor
    LinearMieGruneisen() = default;

    //! Constructor
    LinearMieGruneisen(
      tk::real gamma0,
      tk::real rho0,
      tk::real alpha,
      tk::real c0,
      tk::real s1,
      tk::real cv,
      tk::real mu );

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

    //! Calculate cold-compression component of pressure
    tk::real pressure_coldcompr(
      tk::real arho,
      tk::real alpha=1.0 ) const;

    //! \brief Calculate the elastic Cauchy stress tensor from the material
    //!   inverse deformation gradient tensor using the LinearMieGruneisen EOS
    std::array< std::array< tk::real, 3 >, 3 >
    CauchyStress(
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
      tk::real arho,
      tk::real u,
      tk::real v,
      tk::real w,
      tk::real apr,
      tk::real alpha=1.0,
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
      tk::real,
      tk::real ) const;

    //! Compute the reference density
    tk::real refDensity() const { return density(refPressure(), 300.0); }

    //! Compute the reference pressure
    tk::real refPressure() const { return 1.0e5; }

    //! Return initial density
    tk::real rho0() const { return m_rho0; }

    //! Return gas constant (no-op)
    tk::real gas_constant() const { return 0.0; }

    //! Return internal energy (no-op)
    tk::real internalenergy(tk::real temp) const { return m_cv * temp; }

    //! Return specific heat (no-op)
    tk::real cv( [[maybe_unused]] tk::real temp) const { return m_cv; }

    //! Return specific heat at constant pressure
    tk::real cp( tk::real ) const { return 0.0; }

    //! Return dynamic viscosity coefficient
    tk::real viscCoeff( tk::real ) const { return 0.0; }

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) /*override*/ {
      p | m_gamma0;
      p | m_rho0;
      p | m_alpha;
      p | m_c0;
      p | m_s1;
      p | m_cv;
      p | m_mu;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i LinearMieGruneisen object reference
    friend void operator|( PUP::er& p, LinearMieGruneisen& i ) { i.pup(p); }
    //@}
};

} //inciter::

#endif // LinearMieGruneisen_h
