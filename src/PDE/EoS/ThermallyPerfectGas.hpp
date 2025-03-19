// *****************************************************************************
/*!
  \file      src/PDE/EoS/ThermallyPerfectGas.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Thermally perfect gas equation of state
  \details   This file declares functions for the thermally perfect gas equation
             of state for the compressible flow equations.
*/
// *****************************************************************************
#ifndef ThermallyPerfectGas_h
#define ThermallyPerfectGas_h

#include "Data.hpp"

namespace inciter {

class ThermallyPerfectGas {

  private:
    tk::real m_gamma;
    tk::real m_R;
    std::vector< std::vector< tk::real > > m_cp_coeff{3, std::vector< tk::real >(8)};
    std::vector< tk::real > m_t_range{std::vector< tk::real >(4)};
    tk::real m_dH_ref;

    std::size_t get_t_range(tk::real temp) const
    // *************************************************************************
    //! Check what temperature range, if any, the given temperature is in
    //! \param[in] temp Given temperature to be checked for range
    //! \return Index of temperature range the given temperature is in
    // *************************************************************************
    {
      tk::real fdg = 0.1; // Fudge factor to accomodate numerical overshoot
      if (temp < m_t_range[0] * (1 - fdg) || temp > m_t_range.back() * (1 + fdg)) {
        Throw("ThermallyPerfectGas totalenergy temperature outside t_range bounds: "
        + std::to_string(temp));
      }

      std::size_t t_rng_idx(0);
      for (std::size_t k = 0; k < m_t_range.size() - 1; k++) {
        // Apply fudge factor to max/min bounds
        tk::real fdgl = 1., fdgu = 1.;
        if (k == 0) {
            fdgl = 1 - fdg;
        } else if (k == m_t_range.size() - 2) {
            fdgu = 1 + fdg;
        }
        if (temp >= m_t_range[k] * fdgl && temp < m_t_range[k+1] * fdgu) {
          t_rng_idx = k;
          break;
        }
      }
      return t_rng_idx;
    }


  public:
    //! Default constructor
    ThermallyPerfectGas() = default;

    //! Constructor
    ThermallyPerfectGas(
      tk::real gamma,
      tk::real R,
      std::vector< std::vector< tk::real > > cp_coeff,
      std::vector< tk::real > t_range,
      tk::real dH_ref);

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
                       tk::real alpha=1.0,
                       std::size_t imat=0,
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
                         tk::real alpha=1.0,
                         std::size_t imat=0,
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
                          tk::real alpha=1.0,
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}} ) const;

    //! Compute the minimum allowed pressure
    tk::real min_eff_pressure(
      tk::real min,
      tk::real,
      tk::real ) const
    { return min; }

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
      p | m_cp_coeff;
      p | m_t_range;
      p | m_dH_ref;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i ThermallyPerfectGas object reference
    friend void operator|( PUP::er& p, ThermallyPerfectGas& i ) { i.pup(p); }
    //@}
};

} //inciter::

#endif // ThermallyPerfectGas_h
