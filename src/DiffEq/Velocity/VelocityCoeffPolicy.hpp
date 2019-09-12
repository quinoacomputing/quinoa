// *****************************************************************************
/*!
  \file      src/DiffEq/Velocity/VelocityCoeffPolicy.hpp
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

#include "Types.hpp"
#include "Table.hpp"
#include "SystemComponents.hpp"
#include "Walker/Options/CoeffPolicy.hpp"
#include "Walker/Options/VelocityVariant.hpp"
#include "Langevin.hpp"

namespace walker {

//! Velocity equation coefficients policy with prescribed mean shear
//! \see kw::const_shear_info
class VelocityCoeffConstShear {

  public:
    //! Constructor: initialize coefficients
    VelocityCoeffConstShear( kw::sde_c0::info::expect::type C0_,
                             kw::sde_c0::info::expect::type& C0,
                             std::array< tk::real, 9 >& dU );
      
    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONST_SHEAR; }

    //! Update the model coefficients prescribing shear
    void update( char depvar,
                 char dissipation_depvar,
                 const std::map< tk::ctr::Product, tk::real >& moments,
                 const tk::Table&,
                 ctr::DepvarType solve,
                 ctr::VelocityVariantType variant,
                 kw::sde_c0::info::expect::type C0,
                 tk::real,
                 tk::real& eps,
                 std::array< tk::real, 9 >& G ) const;

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
    VelocityCoeffStationary( kw::sde_c0::info::expect::type C0_,
                             kw::sde_c0::info::expect::type& C0,
                             std::array< tk::real, 9 >& );

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::STATIONARY; }

    //! Update the model coefficients forcing a statistically stationary PDF
    void update( char /*depvar*/,
                 char,
                 const std::map< tk::ctr::Product, tk::real >&,
                 const tk::Table&,
                 ctr::DepvarType,
                 ctr::VelocityVariantType,
                 kw::sde_c0::info::expect::type C0,
                 tk::real,
                 tk::real& eps,
                 std::array< tk::real, 9 >& G ) const;
};

//! Velocity equation coefficients policy with DNS hydrodynamics time scale
//! \see kw::hydrotimescale_info
class VelocityCoeffHydroTimeScale {

  public:
    //! Constructor: initialize coefficients
    VelocityCoeffHydroTimeScale( kw::sde_c0::info::expect::type C0_,
                                 kw::sde_c0::info::expect::type& C0,
                                 std::array< tk::real, 9 >& );

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::HYDROTIMESCALE; }

    //! \brief Update the model coefficients sampling the hydrodynamics time
    //!   scale from a prescribed function table
    void update( char depvar,
                 char,
                 const std::map< tk::ctr::Product, tk::real >& moments,
                 const tk::Table& hts,
                 ctr::DepvarType solve,
                 ctr::VelocityVariantType,
                 kw::sde_c0::info::expect::type C0,
                 tk::real t,
                 tk::real& eps,
                 std::array< tk::real, 9 >& G ) const;
};

//! List of all Velocity's coefficients policies
using VelocityCoeffPolicies = brigand::list< VelocityCoeffConstShear
                                           , VelocityCoeffStationary
                                           , VelocityCoeffHydroTimeScale
                                           >;

} // walker::

#endif // VelocityCoeffPolicy_h
