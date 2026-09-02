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
#include "EoS/EOSDeviceFn.hpp"

namespace inciter {

class ThermallyPerfectGas {

  private:
    // Fixed NASA-9 extents; hard-coded in LuaParser.cpp
    static constexpr std::size_t s_nintervals = 3;
    static constexpr std::size_t s_ncoeff = 8;
    static constexpr std::size_t s_ntrange = 4;

    tk::real m_R;
    // Used raw C arrays here to avoid constexpr issue
    tk::real m_cp_coeff[s_nintervals][s_ncoeff]{};
    tk::real m_t_range[s_ntrange]{};

    //std::vector< std::vector< tk::real > > m_cp_coeff{3, std::vector< tk::real >(8)};
    //std::vector< tk::real > m_t_range{std::vector< tk::real >(4)};
    tk::real m_dH_ref;
    tk::real m_mu;
    tk::real m_temp_ref;
    tk::real m_mu_ref;
    tk::real m_C;
    bool m_Sutherland;

    void get_t_range( tk::real &temp_poly,
                             std::size_t &t_rng_idx ) const
    // *************************************************************************
    //! \brief Check what temperature range the given temperature is in. If it
    //!   exceeds the bounds, reset the temp to the bounds.
    //! \param[in] temp Given temperature to be checked for range
    //! \return Index of temperature range the given temperature is in
    // *************************************************************************
    {
      // First, bounds check
      if (temp_poly < m_t_range[0]) {
        t_rng_idx = 0;
        temp_poly = m_t_range[0];
      } else if (temp_poly > m_t_range[s_ntrange-1]) {
        t_rng_idx = s_ntrange - 2;  //m_t_range.size() - 2;
        temp_poly = m_t_range[s_ntrange-1]; //.back();
      } else {
      // Valid bounds
        for (std::size_t k = 0; k < s_ntrange - 1; k++) {
          if (temp_poly >= m_t_range[k] && temp_poly <= m_t_range[k+1]) {
            t_rng_idx = k;
            break;
          }
        }
      }
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
      std::size_t t_rng_idx(0); // Reference the correct polynomial.
      tk::real temp_poly = temp; // temp to use in polynomial expression
      get_t_range(temp_poly, t_rng_idx); // Bounds check performed inside

      // h = h_poly(T) + h_ref
      tk::real h = -m_cp_coeff[t_rng_idx][0] * std::pow(temp_poly, -2) +
          m_cp_coeff[t_rng_idx][1] * std::log(temp_poly) / temp_poly +
          m_cp_coeff[t_rng_idx][2] +
          m_cp_coeff[t_rng_idx][3] * temp_poly / 2. +
          m_cp_coeff[t_rng_idx][4] * std::pow(temp_poly, 2) / 3. +
          m_cp_coeff[t_rng_idx][5] * std::pow(temp_poly, 3) / 4. +
          m_cp_coeff[t_rng_idx][6] * std::pow(temp_poly, 4) / 5. +
          m_cp_coeff[t_rng_idx][7] / temp_poly;

      // If bounds exceeded, temp_poly will be different than temp. Apply correction.
      if (std::abs(temp_poly - temp) > std::numeric_limits< tk::real >::epsilon()) {
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
      std::size_t t_rng_idx(0); // Reference the correct polynomial.
      tk::real temp_poly = temp; // temp to use in polynomial expression
      get_t_range(temp_poly, t_rng_idx); // Bounds check performed inside

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
      tk::real R,
      std::vector< std::vector< tk::real > > cp_coeff,
      std::vector< tk::real > t_range,
      tk::real dH_ref,
      tk::real mu,
      tk::real temp_ref,
      tk::real mu_ref,
      tk::real C,
      bool Sutherland);

    //! Set rho0 EOS parameter. No-op.
    void setRho0(tk::real) {}

    //! Calculate density from the material pressure and temperature
    [[noreturn]] tk::real density( tk::real pr,
                      tk::real temp ) const;

    //! Calculate pressure from the material density, momentum and total energy
    [[noreturn]] tk::real pressure( tk::real rho,
                       tk::real u,
                       tk::real v,
                       tk::real w,
                       tk::real rhoE,
                       tk::real alpha=1.0,
                       std::size_t imat=0,
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}}) const;

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
    EOS_FN tk::real soundspeed( tk::real rho,
                         tk::real pr,
                         tk::real alpha=1.0,
                         std::size_t imat=0,
      const std::array< std::array< tk::real, 3 >, 3 >& adefgrad={{}},
      const std::array< tk::real, 3 >& asigman={{}} ) const//;
    {
#if defined(__CUDA_ARCH__)
      // Multimat never constructs a TPG so the host overload will throw
      (void)rho; (void)pr; (void)alpha; (void)imat;
      (void)adefgrad; (void)asigman;
      tk::real z=0.0;
      return z/z; //loud NaN will yell at us if reached
#else
      soundspeedHost( rho, pr, alpha, imat, adefgrad, asigman );
#endif
    }

    [[noreturn]] tk::real soundspeedHost( tk::real rho,
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
    [[noreturn]] tk::real totalenergy( tk::real arho,
                          tk::real u,
                          tk::real v,
                          tk::real w,
                          tk::real apr,
                          tk::real alpha=1.0,
      const std::array< std::array< tk::real, 3 >, 3 >& defgrad={{}} ) const;

    //! \brief Device overload of totalenergy: not implemented
    //! \details Exists only so that EOS::compute< EOS::totalenergy > can be
    //!   instantiated for device code. compute() selects the material with a
    //!   runtime if, so every branch must compile even though multimat never
    //!   constructs a ThermallyPerfectGas. There is no TPG totalenergy to port
    //!   -- the host version above is [[noreturn]] and throws.
    //!   checkDeviceEOSSupport() rejects this material on the host before any
    //!   device kernel runs, so the NaN below should be unreachable.
    EOS_FN tk::real totalenergy( tk::real arho,
                                 tk::real u,
                                 tk::real v,
                                 tk::real w,
                                 tk::real apr,
                                 tk::real alpha,
      const tk::real /*defgrad*/[3][3] ) const
    {
#if defined(__CUDA_ARCH__)
      (void)arho; (void)u; (void)v; (void)w; (void)apr; (void)alpha;
      tk::real z=0.0;
      return z/z; //loud NaN will yell at us if reached
#else
      Throw("Direct call to TPG totalenergy should not occur. Use Mixture class.");
#endif
    }

    //! \brief Calculate material temperature from the material density, and
    //!   material specific total energy
    [[noreturn]] tk::real temperature( tk::real rho,
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

    //! Return gas constant (species specific)
    tk::real gas_constant() const { return m_R; }

    //! Return species internal energy
    tk::real internalenergy(tk::real temp) const;

    //! Return species specific heat (constant volume)
    tk::real cv(tk::real temp) const;

    //! Return species specific heat at constant pressure
    tk::real cp(tk::real temp) const;

    //! Return dynamic viscosity coefficient
    tk::real viscCoeff(tk::real temp) const;

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) /*override*/ {
      p | m_R;
      //p | m_cp_coeff;
      for (std::size_t i=0; i<s_nintervals; ++i){
        for (std::size_t j=0; j<s_ncoeff; ++j)
          p | m_cp_coeff[i][j];
      }
      for (std::size_t i=0; i<s_ntrange; ++i)
        p | m_t_range[i];
      //p | m_t_range;
      p | m_dH_ref;
      p | m_mu;
      p | m_temp_ref;
      p | m_mu_ref;
      p | m_C;
      p | m_Sutherland;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i ThermallyPerfectGas object reference
    friend void operator|( PUP::er& p, ThermallyPerfectGas& i ) { i.pup(p); }
    //@}
};

static_assert( std::is_trivially_copyable_v< ThermallyPerfectGas >, "ThermallyPerfectGas must stay POD for the EOS union." );

} //inciter::

#endif // ThermallyPerfectGas_h
