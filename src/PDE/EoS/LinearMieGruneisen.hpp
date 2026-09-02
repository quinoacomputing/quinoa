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
#include <cmath>
#include <iostream>
#include "EoS/EOSDeviceFn.hpp"
#include "EoS/TensorEOSDev.hpp"

namespace inciter {

class LinearMieGruneisen {

  private:
    tk::real m_gamma0, m_rho0, m_alpha, m_c0, m_s1, m_cv, m_mu;

    //! \brief Calculate elastic contribution to material energy from the
    //!   material density, and deformation gradient tensor
    tk::real elasticEnergy(
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
      tk::real& eps2 ) const;

    //! \brief Device version of elasticEnergy
    //! \details Mirror of LinearMieGruneisen::elasticEnergy in
    //!   LinearMieGruneisen.cpp, with tk::getIsochorRightCauchyGreen ->
    //!   tk::isochorRightCauchyGreenEOS. Distinct name rather than an
    //!   overload, for the reason given in SmallShearSolid.hpp.
    //! \warning Not bit-identical to the host version: the replacement inverts
    //!   g.g^T by adjugate rather than by LAPACK LU. See the warning in
    //!   EoS/TensorEOSDev.hpp for measured deviations.
    EOS_FN tk::real elasticEnergyDev(
      const tk::real defgrad[3][3],
      tk::real& eps2 ) const
    {
      // compute volume-preserving part of Right Cauchy-Green strain tensor
      tk::real Ct[3][3];
      tk::isochorRightCauchyGreenEOS(defgrad, Ct);

      // compute elastic shear distortion
      eps2 = 0.5 * (Ct[0][0]+Ct[1][1]+Ct[2][2] - 3.0);

      // compute elastic energy
      auto rhoEe = m_mu * eps2;

      return rhoEe;
    }

    //! Density-dependent Gruneisen coefficient
    //! \details Moved inline (from a .cpp definition), along with the six
    //!   helpers below, so soundspeed() can call them from device code.
    EOS_FN tk::real gruneisen( tk::real rho ) const
    {
      return m_gamma0 * std::pow(m_rho0/rho, m_alpha);
    }

    //! Derivative of density-dependent Gruneisen coefficient wrt density
    EOS_FN tk::real dGruneisenDrho( tk::real rho ) const
    {
      return (- m_alpha * gruneisen(rho) / rho);
    }

    //! Compression measure based on reference density
    EOS_FN tk::real eta( tk::real rho ) const
    {
      return 1.0 - m_rho0/rho;
    }

    //! Derivative of compression measure wrt density
    EOS_FN tk::real dEtaDrho( tk::real rho ) const
    {
      return m_rho0/(rho*rho);
    }

    //! Linear Us-Up Hugoniot pressure
    EOS_FN tk::real hugoniotPressure( tk::real rho ) const
    {
      const auto et = eta(rho);
      const auto den = 1.0 - m_s1*et;
      return m_rho0*m_c0*m_c0*et/(den*den);
    }

    //! Derivative of Hugoniot pressure wrt density
    EOS_FN tk::real dHugoniotPressureDrho( tk::real rho ) const
    {
      const auto et = eta(rho);
      const auto den = 1.0 - m_s1*et;
      return m_rho0*m_c0*m_c0*(1.0 + m_s1*et)/(den*den*den)*
        dEtaDrho(rho);
    }

    //! Hugoniot internal energy
    //! \details Moved inline now that the device overload of totalenergy()
    //!   needs it. (An earlier revision left this out-of-line on the grounds
    //!   that only host code called it; that is no longer true.)
    EOS_FN tk::real hugoniotEnergy( tk::real rho ) const
    {
      return 0.5*hugoniotPressure(rho)*eta(rho)/m_rho0;
    }

    //! Derivative of Hugoniot internal energy wrt density
    EOS_FN tk::real dHugoniotEnergyDrho( tk::real rho ) const
    {
      return 0.5/m_rho0*( dHugoniotPressureDrho(rho)*eta(rho) +
        hugoniotPressure(rho)*dEtaDrho(rho) );
    }

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
    //! \details Moved inline as EOS_FN; see the note on the StiffenedGas
    //!   equivalent. The 50-iteration Newton loop is preserved verbatim apart
    //!   from std::abs/std::max -> fabs/fmax (both are constexpr in libstdc++
    //!   and nvcc rejects them without --expt-relaxed-constexpr). All the
    //!   Hugoniot/Gruneisen helpers it calls are already device-callable.
    EOS_FN tk::real density( tk::real pr,
                             tk::real temp ) const
    {
      auto rho = m_rho0;
      const std::size_t maxiter = 50;
      const tk::real tol = 1.0e-10;

      for (std::size_t iter=0; iter<maxiter; ++iter) {
        const auto p = hugoniotPressure(rho) + gruneisen(rho)*rho*m_cv*temp;
        const auto dpdrho = dHugoniotPressureDrho(rho) +
          (gruneisen(rho) + rho*dGruneisenDrho(rho))*m_cv*temp;
        const auto rhoold = rho;
        const auto delta = (p - pr)/dpdrho;
        rho -= delta;
        if (rho <= 0.0) rho = 0.5*rhoold;
        if (fabs(delta) <= tol*fmax(1.0, fabs(rho))) break;
      }

      return rho;
    }

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
      const std::array< std::array< tk::real, 3 >, 3 >& adefgrad={{}} ) const;

    //! Calculate speed of sound from the material density and material pressure
    EOS_FN tk::real soundspeed(
      tk::real arho,
      tk::real apr,
      tk::real alpha=1.0,
      std::size_t imat=0,
      const std::array< std::array< tk::real, 3 >, 3 >& /*adefgrad*/={{}} ) const
    {
      // Approximated elastic contribution, from Barton, P. T. (2019).
      // An interface-capturing Godunov method for the simulation of
      // compressible solid-fluid problems. Journal of Computational Physics,
      // 390, 25-50
      auto al_eff = fmax( 1.0e-14, alpha );
      tk::real a = (4.0/3.0) * m_mu * al_eff / arho;

      // Hydrodynamic contribution from the Mie-Gruneisen derivatives.
      const auto rho = arho/al_eff;
      const auto gamma = gruneisen(rho);
      const auto e_thermal =
        (apr - al_eff*hugoniotPressure(rho))/(gamma*arho);
      a += dHugoniotPressureDrho(rho)
         + (gamma + rho*dGruneisenDrho(rho))*e_thermal
         - gamma*rho*dHugoniotEnergyDrho(rho)
         + gamma*apr/arho;

      // Compute square root
      a = std::sqrt(fmax(1.0e-15, a));

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
    //!   LinearMieGruneisen::totalenergy in LinearMieGruneisen.cpp.
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
      const auto rho = arho/alpha;

      // obtain hydrodynamic and kinetic energy
      tk::real arhoEh = arho*hugoniotEnergy(rho) +
        (apr - alpha*hugoniotPressure(rho))/gruneisen(rho) +
        0.5 * arho * (u*u + v*v + w*w);
      // obtain elastic contribution to energy
      tk::real eps2;
      tk::real arhoEe = alpha*elasticEnergyDev(defgrad, eps2);

      return (arhoEh + arhoEe);
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

    //! \brief Device overload of temperature
    //! \details Takes defgrad as a raw C array; see the note on the device
    //!   totalenergy overload. Arithmetic and its ordering are identical to
    //!   LinearMieGruneisen::temperature in LinearMieGruneisen.cpp.
    //! \warning Not bit-identical to the host version, via elasticEnergyDev.
    EOS_FN tk::real temperature(
      tk::real arho,
      tk::real u,
      tk::real v,
      tk::real w,
      tk::real arhoE,
      tk::real alpha,
      const tk::real defgrad[3][3] ) const
    {
      // obtain elastic contribution to energy
      tk::real eps2;
      auto arhoEe = alpha*elasticEnergyDev(defgrad, eps2);
      // obtain hydrodynamic internal energy
      auto arhoEi = arhoE - arhoEe - 0.5*arho*(u*u + v*v + w*w);

      const auto rho = arho/alpha;
      const auto t = (arhoEi/arho - hugoniotEnergy(rho))/m_cv;
      return t;
    }

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
