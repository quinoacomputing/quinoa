// *****************************************************************************
/*!
  \file      src/PDE/EoS/WilkinsAluminum.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Wilkins equation of state for aluminum
  \details   This file declares functions for the Wilkins equation of
             state for solids and a hydro EoS for aluminum. These functions were
             taken from Example 4 of Barton, Philip T. "An interface-capturing
             Godunov method for the simulation of compressible solid-fluid
             problems." Journal of Computational Physics 390 (2019): 25-50.
*/
// *****************************************************************************
#ifndef WilkinsAluminum_h
#define WilkinsAluminum_h

#include "Data.hpp"
#include <cmath>
#include <iostream>
#include "EoS/EOSDeviceFn.hpp"
#include "EoS/TensorEOSDev.hpp"

namespace inciter {

class WilkinsAluminum {

  private:
    tk::real m_gamma, m_cv, m_mu, m_rho0;

    //! Wilkins-aluminum fit constants (Barton 2019, Example 4)
    //! \details Moved here from file-scope statics in WilkinsAluminum.cpp so
    //!   they're visible to the device compile of soundspeed() below. This is
    //!   now their one definition; WilkinsAluminum.cpp's uses of e1..e5
    //!   resolve to these members unqualified, same as before.
    static constexpr tk::real e1 = -13.0e+09;
    static constexpr tk::real e2 = 20.0e+09;
    static constexpr tk::real e3 = 52.0e+09;
    static constexpr tk::real e4 = -59.0e+09;
    static constexpr tk::real e5 = 151.0e+09;

    //! \brief Calculate elastic contribution to material energy from the
    //!   material density, and deformation gradient tensor
    tk::real elasticEnergy(
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
      std::array< std::array< tk::real, 3 >, 3 >& devH ) const;

    //! \brief Device version of elasticEnergy
    //! \details Mirror of WilkinsAluminum::elasticEnergy in
    //!   WilkinsAluminum.cpp, with tk::getDevHencky -> tk::devHenckyEOS.
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

  public:
    //! Default constructor
    WilkinsAluminum() = default;

    //! Constructor
    WilkinsAluminum(tk::real gamma, tk::real cv, tk::real mu );

    //! Set rho0 EOS parameter; i.e. the initial density
    void setRho0(tk::real rho0);

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

    //! Calculate cold-compression component of pressure (no-op)
    tk::real pressure_coldcompr(
      tk::real,
      tk::real ) const
    { return 0.0; }

    //! \brief Calculate the elastic Cauchy stress tensor from the material
    //!   inverse deformation gradient tensor using the WilkinsAluminum EOS
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
      const std::array< std::array< tk::real, 3 >, 3 >& /*adefgrad*/={{}} ) const
    {
      tk::real a = 0.0;

      // Hydro contribution
      auto al_eff = fmax( 1.0e-14, alpha );
      tk::real rho0 = m_rho0;
      tk::real rho = arho/al_eff;
      a += fmax( 1.0e-15, 6*e2*std::pow(rho/rho0,2.0)/rho0
                     + 2*e3*rho/(rho0*rho0) - e5/rho0 );

      // Shear contribution
      a += (4.0/3.0) * m_mu / (arho/al_eff);

      // Compute square root
      a = std::sqrt(a);

#if !defined(__CUDA_ARCH__)
      // check sound speed divergence
      if (!std::isfinite(a)) {
        std::cout << "Material-id:      " << imat << std::endl;
        std::cout << "Volume-fraction:  " << alpha << std::endl;
        std::cout << "Partial density:  " << arho << std::endl;
        std::cout << "Partial pressure: " << apr << std::endl;
        Throw("Material-" + std::to_string(imat) + " has nan/inf sound speed: "
          + std::to_string(a) + ", material volume fraction: " +
          std::to_string(alpha));
      }
#else
      (void)imat; (void)apr;
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
    //!   StiffenedGas equivalent. apr is unused, as in the host version
    //!   (where the parameter is unnamed). Arithmetic and its ordering are
    //!   identical to WilkinsAluminum::totalenergy in WilkinsAluminum.cpp.
    //! \warning Not bit-identical to the host version, via elasticEnergyDev.
    EOS_FN tk::real totalenergy(
      tk::real arho,
      tk::real u,
      tk::real v,
      tk::real w,
      tk::real /*apr*/,
      tk::real alpha,
      const tk::real defgrad[3][3] ) const
    {
      // obtain hydro contribution to energy
      tk::real rho0 = m_rho0;
      tk::real rho = arho/alpha;
      tk::real rhoEh = (e1+e2*std::pow(rho/rho0,2.0)+e3*(rho/rho0)
                        +e4*std::pow(rho/rho0,-1.0)-e5*std::log(rho/rho0))/rho0
                       + 0.5*rho*(u*u + v*v + w*w);
      // obtain elastic contribution to energy
      tk::real devH[3][3];
      tk::real rhoEe = elasticEnergyDev(defgrad, devH);

      return alpha*(rhoEh + rhoEe);
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
    tk::real cp( tk::real ) const { return m_gamma*m_cv; }

    //! Return dynamic viscosity coefficient
    tk::real viscCoeff( tk::real ) const { return 0.0; }

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) /*override*/ {
      p | m_gamma;
      p | m_cv;
      p | m_mu;
      p | m_rho0;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i WilkinsAluminum object reference
    friend void operator|( PUP::er& p, WilkinsAluminum& i ) { i.pup(p); }
    //@}
};

} //inciter::

#endif // WilkinsAluminum_h
