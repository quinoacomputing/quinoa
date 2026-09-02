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
#include <cmath>
#include <iostream>
#include "EoS/EOSDeviceFn.hpp"
#include "EoS/TensorEOSDev.hpp"

namespace inciter {

class GodunovRomenski {

  private:
    tk::real m_gamma, m_mu, m_rho0, m_alpha, m_K0;

    //! \brief Calculate elastic contribution to material energy from the
    //!   material density, and deformation gradient tensor
    tk::real elasticEnergy(
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
      std::array< std::array< tk::real, 3 >, 3 >& devH ) const;

    //! \brief Device version of elasticEnergy
    //! \details Mirror of GodunovRomenski::elasticEnergy in
    //!   GodunovRomenski.cpp, with tk::getDevHencky -> tk::devHenckyEOS.
    //!   Distinct name rather than an overload, for the reason given in
    //!   SmallShearSolid.hpp.
    //! \warning Not bit-identical to the host version: devHenckyEOS inverts
    //!   g.g^T by adjugate rather than by LAPACK LU. See the warning in
    //!   EoS/TensorEOSDev.hpp for measured deviations.
    EOS_FN tk::real elasticEnergyDev(
      const tk::real defgrad[3][3],
      tk::real devH[3][3] ) const
    {
      // Compute deviatoric part of Hencky tensor
      tk::devHenckyEOS(defgrad, devH);

      // Compute elastic energy
      tk::real rhoEe = 0.0;
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          rhoEe += m_mu*devH[i][j]*devH[i][j];

      return rhoEe;
    }

    //! \brief Calculate cold-compression contribution to material energy from
    //!   the material density
    //! \details Moved inline (from a .cpp definition) so the device overload
    //!   of totalenergy() below can call it from device code.
    EOS_FN tk::real coldcomprEnergy( tk::real rho ) const
    {
      auto rrho0a = std::pow(rho/m_rho0, m_alpha);
      return ( rho * m_K0/(2.0*m_rho0*m_alpha*m_alpha) * (rrho0a-1.0)*(rrho0a-1.0) );
    }

    //! Calculate the derivative of the cold compression pressure wrt. density
    EOS_FN tk::real DpccDrho( tk::real rho ) const
    {
      auto rrho0a = std::pow(rho/m_rho0, m_alpha);
      return m_K0/(m_rho0*m_alpha) * ((2.0*m_alpha+1.0)*(rrho0a*rrho0a) - (m_alpha+1.0)*rrho0a);
    }

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
    EOS_FN tk::real pressure_coldcompr(
      tk::real arho,
      tk::real alpha=1.0 ) const//;
    {
      auto rho = arho/alpha;
      auto rrho0a = std::pow(rho/m_rho0, m_alpha);
      return alpha * (m_K0/m_alpha * (rrho0a*rho/m_rho0) * (rrho0a-1.0));
    }

    //! \brief Calculate the elastic Cauchy stress tensor from the material
    //!   inverse deformation gradient tensor using the GodunovRomenski EOS
    std::array< std::array< tk::real, 3 >, 3 >
    CauchyStress(
      tk::real alpha,
      std::size_t /*imat*/,
      const std::array< std::array< tk::real, 3 >, 3 >& adefgrad ) const;

    //! Calculate speed of sound from the material density and material pressure
    EOS_FN tk::real soundspeed(
      tk::real arho,
      tk::real apr,
      tk::real alpha=1.0,
      std::size_t imat=0,
      const std::array< std::array< tk::real, 3 >, 3 >& adefgrad={{}} ) const//;
    {
      tk::real a = 0.0;
      auto al_eff = fmax( 1e-14, alpha );
      tk::real rho = arho/al_eff;
      auto p_cc = pressure_coldcompr(arho, al_eff);
      a += fmax( 1.0e-15, DpccDrho(rho) + (m_gamma+1.0) * (apr - p_cc)/arho );
      a += (4.0/3.0) * al_eff * m_mu / arho;
      a = std::sqrt(a);
#if !defined(__CUDA_ARCH__)
      if (!std::isfinite(a)) {
        std::cout << "Material-id: " << imat << std::endl;
        std::cout << "Volume fraction: " << alpha << std::endl;
        std::cout << "Partial density: " << arho << std::endl;
        std::cout << "Partial pressure: " << apr << std::endl;
        Throw("Material " + std::to_string(imat) + " has nan/inf sound speed.");
      }
#else
      (void)imat;
#endif
      return a;
    }

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

    //! \brief Device overload of totalenergy
    //! \details Takes defgrad as a raw C array; see the note on the
    //!   StiffenedGas equivalent. Arithmetic and its ordering are identical to
    //!   GodunovRomenski::totalenergy in GodunovRomenski.cpp.
    //! \warning Not bit-identical to the host version, via elasticEnergyDev.
    EOS_FN tk::real totalenergy(
      tk::real arho,
      tk::real u,
      tk::real v,
      tk::real w,
      tk::real apr,
      tk::real alpha,
      const tk::real defgrad[3][3] ) const
    {
      // obtain thermal and kinetic energy
      auto apt = apr - pressure_coldcompr(arho, alpha);
      // in the above expression, shear pressure is not included in pr in the
      // first place, so should not subtract it
      auto arhoEh = apt/m_gamma + 0.5*arho*(u*u + v*v + w*w);
      // obtain elastic contribution to energy
      tk::real devH[3][3];
      auto arhoEe = alpha*elasticEnergyDev(defgrad, devH);
      auto arhoEc = alpha*coldcomprEnergy(arho/alpha);

      return (arhoEh + arhoEe + arhoEc);
    }

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
      tk::real /*arho*/,
      tk::real /*alpha=1.0*/ ) const;

    //! Compute the reference density
    tk::real refDensity() const { return density(refPressure(), 300.0); }

    //! Compute the reference pressure
    tk::real refPressure() const { return 1.0e5; }

    //! Return initial density
    tk::real rho0() const { return m_rho0; }

    //! Return gas constant (no-op)
    tk::real gas_constant() const { return 0.0; }

    //! Return internal energy (no-op)
    tk::real internalenergy( [[maybe_unused]] tk::real temp) const { return 0.0; }

    //! Return specific heat (no-op)
    tk::real cv( [[maybe_unused]] tk::real temp) const { return 0.0; }

    //! Return specific heat at constant pressure
    tk::real cp( tk::real ) const { return 0.0; }

    //! Return dynamic viscosity coefficient
    tk::real viscCoeff( tk::real ) const { return 0.0; }

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
