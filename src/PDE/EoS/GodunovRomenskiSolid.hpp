// *****************************************************************************
/*!
  \file      src/PDE/EoS/GodunovRomenskiSolid.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Godunov-Romenski equation of state for solids
  \details   This file declares functions for the Godunov-Romenski equation of
             state for solids. These function were mostly taken from Barton,
             Philip T. "An interface-capturing Godunov method for the simulation
             of compressible solid-fluid problems." Journal of Computational
             Physics 390 (2019): 25-50. The elastic energy and stress is
             obtained from a the deviatoric part of the Hencky strain, while the
             hydrodynamics contributions are obtained from a stiffened gas EOS.
*/
// *****************************************************************************
#ifndef GodunovRomenskiSolid_h
#define GodunovRomenskiSolid_h

#include "Data.hpp"

namespace inciter {

class GodunovRomenskiSolid {

  private:
    tk::real m_gamma, m_pstiff, m_cv, m_mu;

    //! \brief Calculate elastic contribution to material energy from the
    //!   material density, and deformation gradient tensor
    tk::real elasticEnergy(
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad,
      std::array< std::array< tk::real, 3 >, 3 >& devH ) const;

  public:
    //! Default constructor
    GodunovRomenskiSolid() = default;

    //! Constructor
    GodunovRomenskiSolid(tk::real gamma, tk::real pstiff, tk::real cv, tk::real mu );


    //! Calculate density from the material pressure and temperature
    tk::real density( tk::real pr,
                      tk::real temp,
                      tk::real rho0=1.0 ) const;

    //! Calculate pressure from the material density, momentum and total energy
    tk::real pressure(
      tk::real arho,
      tk::real u,
      tk::real v,
      tk::real w,
      tk::real arhoE,
      tk::real alpha=1.0,
      std::size_t imat=0,
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}},
      tk::real rho0=1.0 ) const;

    //! \brief Calculate the elastic Cauchy stress tensor from the material
    //!   density, momentum, total energy, and inverse deformation gradient
    //!   tensor using the GodunovRomenskiSolid equation of state
    std::array< std::array< tk::real, 3 >, 3 >
    CauchyStress(
      tk::real,
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
      const std::array< std::array< tk::real, 3 >, 3 >& adefgrad={{}},
      const std::array< tk::real, 3 >& adefgradn={{}},
      const std::array< tk::real, 3 >& asigman={{}},
      tk::real rho0=1.0 ) const;

    //! \brief Calculate material specific total energy from the material
    //!   density, momentum and material pressure
    tk::real totalenergy(
      tk::real rho,
      tk::real u,
      tk::real v,
      tk::real w,
      tk::real pr,
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}},
      tk::real rho0=1.0 ) const;

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

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) /*override*/ {
      p | m_gamma;
      p | m_pstiff;
      p | m_cv;
      p | m_mu;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i GodunovRomenskiSolid object reference
    friend void operator|( PUP::er& p, GodunovRomenskiSolid& i ) { i.pup(p); }
    //@}
};

} //inciter::

#endif // GodunovRomenskiSolid_h
