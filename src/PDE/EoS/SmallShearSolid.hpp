// *****************************************************************************
/*!
  \file      src/PDE/EoS/SmallShearSolid.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Small shear strain equation of state for solids
  \details   This file declares functions for the SmallShearSolid equation of
             state for the compressible flow equations. These functions are
             taken from Plohr, J. N., & Plohr, B. J. (2005). Linearized analysis
             of Richtmyer–Meshkov flow for elastic materials. Journal of Fluid
             Mechanics, 537, 55-89. The SmallShearSolid EOS uses a small-shear
             approximation for the elastic contribution, and a stiffened gas EOS
             for the hydrodynamic contribution of the internal energy.
*/
// *****************************************************************************
#ifndef SmallShearSolid_h
#define SmallShearSolid_h

#include "Data.hpp"
#include <cmath>
#include <iostream>
#include "EoS/EOSDeviceFn.hpp"
#include "EoS/TensorEOSDev.hpp"

namespace inciter {

class SmallShearSolid {

  private:
    tk::real m_gamma, m_pstiff, m_cv, m_mu, m_rho0;

    //! \brief Calculate elastic contribution to material energy from the
    //!   material density, and deformation gradient tensor
    tk::real elasticEnergy(
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
      tk::real& eps2 ) const;

    //! \brief Device version of elasticEnergy
    //! \details Deliberately a distinct name rather than an overload: the
    //!   host/device overload pairs elsewhere in this codebase have repeatedly
    //!   caused a count=1 string replace to hit the wrong call site. Mirror of
    //!   SmallShearSolid::elasticEnergy in SmallShearSolid.cpp, with
    //!   tk::getIsochorRightCauchyGreen -> tk::isochorRightCauchyGreenEOS.
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

  public:
    //! Default constructor
    SmallShearSolid() = default;

    //! Constructor
    SmallShearSolid(tk::real gamma, tk::real pstiff, tk::real cv, tk::real mu );

    //! Set rho0 EOS parameter; i.e. the initial density
    void setRho0(tk::real rho0);

    //! Calculate density from the material pressure and temperature
    //! \details Moved inline as EOS_FN; see the note on the StiffenedGas
    //!   equivalent. Takes only scalars, so no device overload is needed.
    EOS_FN tk::real density( tk::real pr,
                             tk::real temp ) const
    {
      tk::real g = m_gamma;
      tk::real p_c = m_pstiff;
      tk::real c_v = m_cv;

      tk::real rho = (pr + p_c) / ((g-1.0) * c_v * temp);
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

    //! Calculate cold-compression component of pressure (no-op)
    tk::real pressure_coldcompr(
      tk::real,
      tk::real ) const
    { return 0.0; }

    //! \brief Calculate the elastic Cauchy stress tensor from the material
    //!   inverse deformation gradient tensor using the SmallShearSolid EOS
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
      auto al_eff = fmax( 1.0e-14, alpha );
      tk::real a = (4.0/3.0) * m_mu * al_eff / arho;
      auto p_eff = fmax( 1.0e-15, apr+(al_eff*m_pstiff) );
      a += m_gamma * p_eff / arho;
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
    //!   SmallShearSolid::totalenergy in SmallShearSolid.cpp.
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
      // obtain hydro contribution to energy
      tk::real arhoEh = (apr + alpha*m_gamma*m_pstiff) / (m_gamma-1.0) + 0.5 * arho *
        (u*u + v*v + w*w);
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
    //!   SmallShearSolid::temperature in SmallShearSolid.cpp.
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
      // obtain hydro contribution to energy
      auto arhoEh = arhoE - arhoEe;

      tk::real t = (arhoEh - 0.5 * arho * (u*u + v*v + w*w) - alpha*m_pstiff)
                   / (arho*m_cv);
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
    tk::real gas_constant() const { return (m_gamma-1.0)*m_cv; }

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
      p | m_pstiff;
      p | m_cv;
      p | m_mu;
      p | m_rho0;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i SmallShearSolid object reference
    friend void operator|( PUP::er& p, SmallShearSolid& i ) { i.pup(p); }
    //@}
};

} //inciter::

#endif // SmallShearSolid_h
