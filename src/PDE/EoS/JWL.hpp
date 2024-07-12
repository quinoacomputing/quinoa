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

namespace inciter {

class JWL {

  private:
    tk::real m_w, m_cv, m_rho0, m_de, m_rhor, m_tr, m_pr, m_a, m_b, m_r1, m_r2;

    //! Calculate specific internal energy
    tk::real intEnergy( tk::real rho, tk::real pr ) const;

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

    //! Calculate density from the material pressure and temperature
    tk::real density( tk::real pr,
                      tk::real temp,
                      tk::real rho0=1.0 ) const;

    //! Calculate pressure from the material density, momentum and total energy
    tk::real pressure( tk::real arho,
                       tk::real u,
                       tk::real v,
                       tk::real w,
                       tk::real arhoE,
                       tk::real alpha=1.0,
                       std::size_t imat=0,
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}},
                       tk::real rho0=1.0 ) const;

    //! \brief Calculate the Cauchy stress tensor from the material density,
    //!   momentum, and total energy
    std::array< std::array< tk::real, 3 >, 3 >
    CauchyStress(
      tk::real,
      tk::real,
      tk::real,
      tk::real,
      tk::real,
      tk::real,
      std::size_t,
      const std::array< std::array< tk::real, 3 >, 3 >& adefgrad={{}} ) const;

    //! Calculate speed of sound from the material density and material pressure
    tk::real soundspeed( tk::real arho,
                         tk::real apr,
                         tk::real alpha=1.0,
                         std::size_t imat=0,
      const std::array< std::array< tk::real, 3 >, 3 >& adefgrad={{}},
      const std::array< tk::real, 3 >& adefgradn={{}},
      const std::array< tk::real, 3 >& asigman={{}},
                         tk::real rho0=1.0 ) const;

    //! Calculate speed of shear waves
    tk::real shearspeed(
      tk::real,
      tk::real,
      std::size_t ) const { return 0.0; }

    //! \brief Calculate material specific total energy from the material
    //!   density, momentum and material pressure
    tk::real totalenergy( tk::real rho,
                          tk::real u,
                          tk::real v,
                          tk::real w,
                          tk::real pr,
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}},
                          tk::real rho0=1.0 ) const;

    //! \brief Calculate material temperature from the material density, and
    //!   material specific total energy
    tk::real temperature( tk::real arho,
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
      tk::real alpha ) const;

    //! Compute the reference density
    //! \details Returns the reference density
    tk::real refDensity() const { return m_rhor; }

    //! Compute the reference pressure
    //! \details Returns the reference pressure
    tk::real refPressure() const { return m_pr; }

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
