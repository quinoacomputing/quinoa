// *****************************************************************************
/*!
  \file      src/PDE/EoS/TPG.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Thermally perfect gas equation of state
  \details   This file declares functions for the thermally perfect gas equation
             of state for the compressible flow equations.
*/
// *****************************************************************************
#ifndef TPG_h
#define TPG_h

#include "Data.hpp"

namespace inciter {

class TPG {

  private:
    tk::real m_gamma;
    tk::real m_R;
    std::vector< tk::real > m_cp_TPG{std::vector< tk::real >(8)};

  public:
    //! Default constructor
    TPG() = default;

    //! Constructor
    TPG(tk::real gamma, tk::real R, std::vector< tk::real > cp_TPG);

    //! Set rho0 EOS parameter. No-op.
    void setRho0(tk::real) {}

    //! Calculate density from the material pressure and temperature
    tk::real density( tk::real pr,
                      tk::real temp ) const;

    //! Calculate pressure from the material density, momentum and total energy
    tk::real pressure( tk::real rho,
                       tk::real u,
                       tk::real v,
                       tk::real w,
                       tk::real rhoE,
                       tk::real,
                       std::size_t,
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}}) const;

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
    tk::real soundspeed( tk::real rho,
                         tk::real pr,
                         tk::real,
                         std::size_t,
      const std::array< std::array< tk::real, 3 >, 3 >& adefgrad={{}},
      const std::array< tk::real, 3 >& asigman={{}} ) const;

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
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}} ) const;

    //! \brief Calculate material temperature from the material density, and
    //!   material specific total energy
    tk::real temperature( tk::real rho,
                          tk::real u,
                          tk::real v,
                          tk::real w,
                          tk::real rhoE,
                          tk::real,
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
    tk::real rho0() const { return density(1.0e5, 300.0); }

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) /*override*/ {
      p | m_gamma;
      p | m_R;
      p | m_cp_TPG;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i TPG object reference
    friend void operator|( PUP::er& p, TPG& i ) { i.pup(p); }
    //@}
};

} //inciter::

#endif // TPG_h
