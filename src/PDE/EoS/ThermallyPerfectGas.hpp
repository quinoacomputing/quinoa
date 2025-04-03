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
    //! Check what temperature range the given temperature is in.
    //! A check should be performed before this to make sure the bounds are not
    //! violated.
    //! \param[in] temp Given temperature to be checked for range
    //! \return Index of temperature range the given temperature is in
    // *************************************************************************
    {
      std::size_t t_rng_idx(0);
      for (std::size_t k = 0; k < m_t_range.size() - 1; k++) {
        if (temp >= m_t_range[k] && temp < m_t_range[k+1]) {
          t_rng_idx = k;
          break;
        }
      }
      return t_rng_idx;
    }

    tk::real calc_h(tk::real temp) const
    // *************************************************************************
    //! Calculate dimensionless enthalpy according to the NASA-9 polynomial
    //! \param[in] temp temperature at which to calculate enthalpy
    //! \return dimensionless enthalpy, h / (R * T)
    // *************************************************************************
    {
      // Identify what temperature range this falls in. If it falls outside the
      // bounds, some corrections must be applied.
      std::size_t t_rng_idx; // Reference the correct polynomial.
      tk::real temp_poly = temp; // temp to use in polynomial expression
      bool correction_needed = false; // If bounds are violated
      if (temp < m_t_range[0]) {
        t_rng_idx = 0;
        temp_poly = m_t_range[0];
        correction_needed = true;
      } else if (temp > m_t_range.back()) {
        t_rng_idx = m_t_range.size() - 1;
        temp_poly = m_t_range.back();
        correction_needed = true;
      } else {
        t_rng_idx = get_t_range(temp);
      }

      // h = h_poly(T) + h_ref
      tk::real h = -m_cp_coeff[t_rng_idx][0] * std::pow(temp_poly, -2) +
          m_cp_coeff[t_rng_idx][1] * std::log(temp_poly) / temp_poly +
          m_cp_coeff[t_rng_idx][2] +
          m_cp_coeff[t_rng_idx][3] * temp_poly / 2. +
          m_cp_coeff[t_rng_idx][4] * std::pow(temp_poly, 2) / 3. +
          m_cp_coeff[t_rng_idx][5] * std::pow(temp_poly, 3) / 4. +
          m_cp_coeff[t_rng_idx][6] * std::pow(temp_poly, 4) / 5. +
          m_cp_coeff[t_rng_idx][7] / temp_poly;

      // If bounds exceeded, apply correction to enthalpy
      if (correction_needed) {
        tk::real cp_star = calc_cp(temp_poly);
        h = h * temp_poly / temp + (temp - temp_poly) / temp * cp_star;
      }

      return h;
    }

    tk::real calc_cp(tk::real temp) const
    // *************************************************************************
    //! Calculate dimensionless specific heat according to the NASA-9 polynomial
    //! \param[in] temp temperature at which to calculate specific heat
    //! \return dimensionless enthalpy, c_p / R
    // *************************************************************************
    {
      // Identify what temperature range this falls in. If it falls outside the
      // bounds, some corrections must be applied.
      std::size_t t_rng_idx; // Reference the correct polynomial.
      tk::real temp_poly = temp; // temp to use in polynomial expression
      if (temp < m_t_range[0]) {
        t_rng_idx = 0;
        temp_poly = m_t_range[0];
      } else if (temp > m_t_range.back()) {
        t_rng_idx = m_t_range.size() - 1;
        temp_poly = m_t_range.back();
      } else {
        t_rng_idx = get_t_range(temp);
      }

      tk::real cp = m_cp_coeff[t_rng_idx][0] * std::pow(temp_poly, -2) +
          m_cp_coeff[t_rng_idx][1] / temp_poly +
          m_cp_coeff[t_rng_idx][2] +
          m_cp_coeff[t_rng_idx][3] * temp_poly +
          m_cp_coeff[t_rng_idx][4] * std::pow(temp_poly, 2) +
          m_cp_coeff[t_rng_idx][5] * std::pow(temp_poly, 3) +
          m_cp_coeff[t_rng_idx][6] * std::pow(temp_poly, 4);

      return cp;
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
