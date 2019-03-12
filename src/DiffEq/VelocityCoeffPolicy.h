// *****************************************************************************
/*!
  \file      src/DiffEq/VelocityCoeffPolicy.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Velocity equation coefficients policies
  \details   This file defines coefficients policy classes for the velocity
     equation for the fluctuating velocity in variable-density turbulence,
     defined in DiffEq/Velocity.h.

    General requirements on velocity equation coefficients policy classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficient, C0. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          kw::sde_c0::info::expect::type C0_,
          kw::sde_c0::info::expect::type& C0,
          std::array< tk::real, 9 >& dU )
      \endcode
      where
      - C0_ denotes a real value used to initialize the velocity system.
      - The reference C0 is to be initialized based on C0_.
      - _dU_ is an optionally prescribed mean velocity gradient.

    - Must define the function _update()_, called from Velocity::advance(),
      updating the model coefficients.
      Required signature:
      \code{.cpp}
        void update( char depvar,
                     char dissipation_depvar,
                     const std::map< tk::ctr::Product, tk::real >& moments,
                     const tk::Table& hts,
                     ctr::DepvarType solve,
                     ctr::VelocityVariantType variant,
                     kw::sde_c0::info::expect::type C0,
                     tk::real t,
                     tk::real& eps,
                     std::array< tk::real, 9 >& G ) const
      \endcode
      where _depvar_ is the dependent variable of the velocity equation,
      _dissipation_depvar_ is the dependent variable of the coupled dissipation
      equation, _moments_ if the map of computed statistical moments, _hts_ is
      a ctr::Table containing the inverse hydrodynamic timescale, _variant_ is
      the velocity model variant, _solve_ is the the lable of the dependent
      variable to solve for (full variable or fluctuation), _C0_ is the Langevin
      eq constat to use, _t_ is the physical time, _eps_ is the dissipation rate
      of turbulent kinetic energy to update, and _G_ is the G_{ij} tensor in the
      Langevin model to update.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::CoeffPolicyType type() noexcept {
          return ctr::CoeffPolicyType::HYDROTIMESCALE;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef VelocityCoeffPolicy_h
#define VelocityCoeffPolicy_h

#include <array>

#include <brigand/sequences/list.hpp>

#include "Types.h"
#include "SystemComponents.h"
#include "Walker/Options/CoeffPolicy.h"
#include "Langevin.h"

namespace walker {

//! Velocity equation coefficients policy with prescribed mean shear
//! \details C0 is user-defined and we prescibe a hard-coded mean shear in the x
//!   direction
//! \see kw::const_shear_info
class VelocityCoeffConstShear {

  public:
    //! Constructor: initialize coefficients
    //! \param[in] C0_ Value of C0 parameter in the Langevin model
    //! \param[in,out] C0 Value of to set the C0 parameter in the Langevin model
    //! \param[in,out] dU Prescribed mean velocity gradient1
    VelocityCoeffConstShear( kw::sde_c0::info::expect::type C0_,
                             kw::sde_c0::info::expect::type& C0,
                             std::array< tk::real, 9 >& dU ) :
      m_dU( {{ 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0 }} )
    {
      C0 = C0_;
      dU = m_dU;
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONST_SHEAR; }

    //! Update the model coefficients (prescribing shear)
    //! \details Update the dissipation rate (eps) and G_{ij} based on the
    //!   turbulent kinetic energy (k) for a prescribed honmogeneous shear flow.
    void update( char depvar,
                 char dissipation_depvar,
                 const std::map< tk::ctr::Product, tk::real >& moments,
                 const tk::Table&,
                 ctr::DepvarType solve,
                 ctr::VelocityVariantType variant,
                 kw::sde_c0::info::expect::type C0,
                 tk::real,
                 tk::real& eps,
                 std::array< tk::real, 9 >& G ) const
    {
      using tk::ctr::lookup;
      using tk::ctr::mean;

      // Compute turbulent kinetic energy
      auto rs = reynoldsStress( depvar, solve, moments );

      // Compute turbulent kinetic energy
      auto k = (rs[0] + rs[1] + rs[2]) / 2.0;

      // Access mean turbulence frequency
      tk::real O = lookup( mean(dissipation_depvar,0), moments );

      // compute turbulent kinetic energy dissipation rate
      eps = O*k;

      // update drift tensor based on the Langevin model variant configured
      if (variant == ctr::VelocityVariantType::SLM)     // simplified
        G = slm( O, C0 );
      else if (variant == ctr::VelocityVariantType::GLM)// generalized
        G = glm( O, C0, rs, m_dU );
      else Throw( "Velocity variant type not implemented" );
    }

  private:
    //! Mean velocity gradient prescribed for simpled 1D homogeneous shear
    std::array< tk::real, 9 > m_dU;
};

//! \brief Velocity equation coefficients policy yielding a statistically
//!   stationary state
//! \see kw::stationary
class VelocityCoeffStationary {

  public:
    //! Constructor: initialize coefficients
    //! \param[in] C0_ Value of C0 parameter in the Langevin model
    //! \param[in,out] C0 Value of to set the C0 parameter in the Langevin model
    //! \param[in,out] dU Prescribed mean velocity gradient1
    VelocityCoeffStationary( kw::sde_c0::info::expect::type C0_,
                               kw::sde_c0::info::expect::type& C0,
                               std::array< tk::real, 9 >& dU ) :
      m_dU( {{ 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0 }} )
    {
      C0 = C0_;
      dU = m_dU;
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::STATIONARY; }

    //! Update the model coefficients (prescribing shear)
    //! \details Update the dissipation rate (eps) and G_{ij} based on the
    //!   turbulent kinetic energy (k) for a prescribed honmogeneous shear flow.
    void update( char depvar,
                 char,
                 const std::map< tk::ctr::Product, tk::real >& moments,
                 const tk::Table&,
                 ctr::DepvarType solve,
                 ctr::VelocityVariantType variant,
                 kw::sde_c0::info::expect::type C0,
                 tk::real,
                 tk::real& eps,
                 std::array< tk::real, 9 >& G ) const
    {
      // Compute turbulent kinetic energy
      auto rs = reynoldsStress( depvar, solve, moments );

      // Override turbulent kinetic energy to keep PDF stationary
      tk::real k = 1.0;

      // Override mean turbulence frequency to keep PDF stationary
      tk::real O = 1.0;

      // Compute turbulent kinetic energy dissipation rate
      eps = O*k;

      // update drift tensor based on the Langevin model variant configured
      if (variant == ctr::VelocityVariantType::SLM)     // simplified
        G = slm( O, C0 );
      else if (variant == ctr::VelocityVariantType::GLM)// generalized
        G = glm( O, C0, rs, m_dU );
      else Throw( "Velocity variant type not implemented" );
    }

  private:
    //! Mean velocity gradient prescribed for simpled 1D homogeneous shear
    std::array< tk::real, 9 > m_dU;
};

//! Velocity equation coefficients policy with DNS hydrodynamics time scale
//! \details C0 is user-defined and we pull in a hydrodynamic timescale from an
//!   external function (from DNS).
//! \see kw::hydrotimescale_info
class VelocityCoeffHydroTimeScale {

  public:
    //! Constructor: initialize coefficients
    //! \param[in] C0_ Value of C0 parameter in the Langevin model
    //! \param[in,out] C0 Value of to set the C0 parameter in the Langevin model
    VelocityCoeffHydroTimeScale( kw::sde_c0::info::expect::type C0_,
                                 kw::sde_c0::info::expect::type& C0,
                                 std::array< tk::real, 9 >& )
    {
      C0 = C0_;
    }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::HYDROTIMESCALE; }

    //! \brief Update the model coefficients
    //! \details Update the dissipation rate (eps) based on eps/k (from DNS) and
    //!    the turbulent kinetic energy (k) (from the SDE)
    void update( char depvar,
                 char,
                 const std::map< tk::ctr::Product, tk::real >& moments,
                 const tk::Table& hts,
                 ctr::DepvarType solve,
                 ctr::VelocityVariantType,
                 kw::sde_c0::info::expect::type C0,
                 tk::real t,
                 tk::real& eps,
                 std::array< tk::real, 9 >& G ) const
    {
      // Compute turbulent kinetic energy
      auto k = tke( depvar, solve, moments );

      // Sample hydrodynamics timescale and prod/diss at time t
      auto ts = hydrotimescale( t, hts );  // eps/k

      // compute turbulent kinetic energy dissipation rate
      eps = ts * k;

      // update drift tensor based on the simplified Langevin model
      G.fill( 0.0 );
      G[0] = G[4] = G[8] = -(0.5+0.75*C0) * ts;
    }

    //! Sample the inverse hydrodynamics time scale at time t
    //! \param[in] t Time at which to sample inverse hydrodynamics time scale
    //! \param[in] ts Hydro time scale table to sample
    //! \return Sampled value from discrete table of inverse hydro time scale
    tk::real hydrotimescale( tk::real t, const tk::Table& ts ) const
    { return tk::sample( t, ts ); }
};

//! List of all Velocity's coefficients policies
using VelocityCoeffPolicies = brigand::list< VelocityCoeffConstShear
                                           , VelocityCoeffStationary
                                           , VelocityCoeffHydroTimeScale
                                           >;

} // walker::

#endif // VelocityCoeffPolicy_h
