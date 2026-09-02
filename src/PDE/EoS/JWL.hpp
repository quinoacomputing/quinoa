// *****************************************************************************
/*!
  \file      src/PDE/EoS/JWL.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Jones, Wilkins, and Lee (JWL) equation of state
  \details   This file declares functions for the JWL equation of
             state for the compressible flow equations. These functions are
             taken from 'JWL Equation of State', Menikoff, LA-UR-15-29536.
*/
// *****************************************************************************
#ifndef JWL_h
#define JWL_h

#include "Data.hpp"
#include <cmath>
#include <iostream>
#include "EoS/EOSDeviceFn.hpp"

namespace inciter {

class JWL {

  private:
    tk::real m_w, m_cv, m_rho0, m_de, m_rhor, m_tr, m_pr, m_a, m_b, m_r1, m_r2;

    //! Calculate specific internal energy
    //! \details Moved inline (from a .cpp definition) so the device overload
    //!   of totalenergy() below can call it from device code.
    EOS_FN tk::real intEnergy( tk::real rho, tk::real pr ) const
    {
      tk::real e = - m_de + 1.0/m_w/rho*( pr
                    - m_a*(1.0 - m_w*rho/m_r1/m_rho0)*exp(-m_r1*m_rho0/rho)
                    - m_b*(1.0 - m_w*rho/m_r2/m_rho0)*exp(-m_r2*m_rho0/rho) );

      return e;
    }

    //! \brief Calculate density from known pressure and temperature using
    //!   bisection root finding method
    tk::real bisection( tk::real a, tk::real b, tk::real p_known,
      tk::real t_known ) const;

    //! Calculate pressure from density and temperature
    tk::real PfromRT( tk::real rho, tk::real T) const;

  public:
    //! Default constructor
    JWL() = default;

    //! Constructor
    JWL( tk::real w, tk::real cv, tk::real rho0, tk::real de, tk::real rhor,
         tk::real tr, tk::real pr, tk::real A, tk::real B, tk::real R1,
         tk::real R2 );

    //! Set rho0 EOS parameter. No-op since rho0 is set in JWL ctor
    void setRho0(tk::real) {}

    //! Calculate density from the material pressure and temperature
    //! \details Deliberately NOT ported to the device. The host
    //!   implementation goes through bisection(), which runs up to 1000
    //!   iterations of PfromRT, is preceded by two unbounded bound-expansion
    //!   while loops, and throws on non-convergence. None of that belongs in
    //!   a boundary kernel: the iteration count would diverge badly across a
    //!   warp and the Throw cannot exist on the device.
    //!
    //!   The device branch therefore returns a NaN, mirroring the stub
    //!   pattern the other materials' soundspeed used before it was ported.
    //!   This is unreachable in practice: checkDeviceEOSSupport() permits JWL,
    //!   so a JWL case CAN take the device path, and any device call to
    //!   compute< EOS::density > on a JWL material will silently produce NaN.
    //!   Porting this, or gating JWL out of device boundary conditions, is
    //!   outstanding work.
    EOS_FN tk::real density( tk::real pr,
                             tk::real temp ) const
    {
#if defined(__CUDA_ARCH__)
      (void)pr; (void)temp;
      tk::real z=0.0;
      return z/z; //NaN
#else
      return densityHost( pr, temp );
#endif
    }

    //! Host implementation of density (outofline)
    tk::real densityHost( tk::real pr,
                          tk::real temp ) const;

    //! Calculate pressure from the material density, momentum and total energy
    tk::real pressure( tk::real arho,
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

    //! \brief Calculate the Cauchy stress tensor from the material
    //!   inverse deformation gradient tensor
    std::array< std::array< tk::real, 3 >, 3 >
    CauchyStress(
      tk::real,
      std::size_t,
      const std::array< std::array< tk::real, 3 >, 3 >& adefgrad={{}} ) const;

    //! Calculate speed of sound from the material density and material pressure
    EOS_FN tk::real soundspeed( tk::real arho,
                         tk::real apr,
                         tk::real alpha=1.0,
                         std::size_t imat=0,
      const std::array< std::array< tk::real, 3 >, 3 >& /*adefgrad*/={{}} ) const
    {
      alpha = fmax(1e-14,alpha);
      // limiting pressure to near-zero
      auto apr_eff = fmax(alpha*
        min_eff_pressure(1e-4*fabs(apr/alpha), arho, alpha), apr);

      auto co1 = m_rho0*alpha*alpha/(arho*arho);
      auto co2 = alpha*(1.0+m_w)/arho;

      tk::real ss = m_a*(m_r1*co1 - co2) * exp(-m_r1*alpha*m_rho0/arho)
                  + m_b*(m_r2*co1 - co2) * exp(-m_r2*alpha*m_rho0/arho)
                  + (1.0+m_w)*apr_eff/arho;

      auto ss2 = ss;
      ss = std::sqrt(ss);

#if !defined(__CUDA_ARCH__)
      // check sound speed divergence
      if (!std::isfinite(ss)) {
        std::cout << "Material-id:      " << imat << std::endl;
        std::cout << "Volume-fraction:  " << alpha << std::endl;
        std::cout << "Partial density:  " << arho << std::endl;
        std::cout << "Partial pressure: " << apr << std::endl;
        std::cout << "Min allowed pres: " << alpha*min_eff_pressure(0.0, arho,
          alpha) << std::endl;
        Throw("Material-" + std::to_string(imat) + " has nan/inf sound speed. "
          "ss^2: " + std::to_string(ss2) + ", material volume fraction: " +
          std::to_string(alpha));
      }
#else
      (void)imat; (void)ss2;
#endif

      return ss;
    }

    //! Calculate speed of shear waves
    tk::real shearspeed(
      tk::real,
      tk::real,
      std::size_t ) const { return 0.0; }

    //! \brief Calculate material specific total energy from the material
    //!   density, momentum and material pressure
    tk::real totalenergy( tk::real arho,
                          tk::real u,
                          tk::real v,
                          tk::real w,
                          tk::real apr,
                          tk::real alpha=1.0,
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}} ) const;

    //! \brief Device overload of totalenergy
    //! \details Takes defgrad as a raw C array; see the note on the
    //!   StiffenedGas equivalent. defgrad is unused here, as in the host
    //!   version. Arithmetic and its ordering are identical to
    //!   JWL::totalenergy in JWL.cpp.
    EOS_FN tk::real totalenergy( tk::real arho,
                                 tk::real u,
                                 tk::real v,
                                 tk::real w,
                                 tk::real apr,
                                 tk::real alpha,
      const tk::real /*defgrad*/[3][3] ) const
    {
      tk::real arhoE = arho*intEnergy( arho/alpha, apr/alpha )
                    + 0.5*arho*(u*u + v*v + w*w);

      return arhoE;
    }

    //! \brief Calculate material temperature from the material density, and
    //!   material specific total energy
    tk::real temperature( tk::real arho,
                          tk::real u,
                          tk::real v,
                          tk::real w,
                          tk::real arhoE,
                          tk::real alpha=1.0,
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}} ) const;

    //! \brief Device overload of temperature
    //! \details Takes defgrad as a raw C array; see the note on the device
    //!   totalenergy overload. defgrad is unused here, as in the host version.
    //!   std::max -> fmax. Arithmetic otherwise identical to JWL::temperature
    //!   in JWL.cpp.
    EOS_FN tk::real temperature( tk::real arho,
                                 tk::real u,
                                 tk::real v,
                                 tk::real w,
                                 tk::real arhoE,
                                 tk::real alpha,
      const tk::real /*defgrad*/[3][3] ) const
    {
      alpha = fmax(1e-14,alpha);
      tk::real rho = arho/alpha;

      tk::real t = ((arhoE - 0.5*arho*(u*u + v*v + w*w))/arho + m_de -
        1.0/m_rho0*( m_a/m_r1*exp(-m_r1*m_rho0/rho)
                   + m_b/m_r2*exp(-m_r2*m_rho0/rho) ))/m_cv;

      return t;
    }

    //! Compute the minimum allowed pressure
    //! \details Moved inline (from a .cpp definition) so soundspeed() below
    //!   can call it from device code.
    EOS_FN tk::real min_eff_pressure(
      tk::real min,
      tk::real arho,
      tk::real alpha ) const
    {
      alpha = fmax(1e-14,alpha);
      auto co1 = m_rho0*alpha*alpha/(arho*arho);
      auto co2 = alpha*(1.0+m_w)/arho;

      // minimum pressure is constrained by zero soundspeed.
      auto apr = -(m_a*(m_r1*co1 - co2) * exp(-m_r1*alpha*m_rho0/arho)
                 + m_b*(m_r2*co1 - co2) * exp(-m_r2*alpha*m_rho0/arho))
        * arho/(1.0+m_w);

      return apr/alpha + min;
    }

    //! Compute the reference density
    //! \details Returns the reference density
    tk::real refDensity() const { return m_rhor; }

    //! Compute the reference pressure
    //! \details Returns the reference pressure
    tk::real refPressure() const { return m_pr; }

    //! Return initial density
    tk::real rho0() const { return m_rho0; }

    //! Return gas constant (no-op)
    tk::real gas_constant() const { return 0.0; }

    //! Return internal energy (no-op)
    tk::real internalenergy(tk::real temp) const { return m_cv * temp; }

    //! Return specific heat (no-op)
    tk::real cv( [[maybe_unused]] tk::real temp) const { return m_cv; }

    //! Return specific heat at constant pressure
    tk::real cp( tk::real ) const { return m_cv; }

    //! Return dynamic viscosity coefficient
    tk::real viscCoeff( tk::real ) const { return 0.0; }

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) /*override*/ {
      p | m_w;
      p | m_cv;
      p | m_rho0;
      p | m_de;
      p | m_rhor;
      p | m_pr;
      p | m_a;
      p | m_b;
      p | m_r1;
      p | m_r2;
      p | m_tr;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i JWL object reference
    friend void operator|( PUP::er& p, JWL& i ) { i.pup(p); }
    //@}
};

} //inciter::

#endif // JWL_h
