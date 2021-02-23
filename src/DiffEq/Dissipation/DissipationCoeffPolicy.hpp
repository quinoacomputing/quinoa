// *****************************************************************************
/*!
  \file      src/DiffEq/Dissipation/DissipationCoeffPolicy.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Particle dissipation equation coefficients policies
  \details   This file defines coefficients policy classes for the Lagrangian
    particle dissipation equation defined in DiffEq/Dissipation.h.

    General requirements on dissipation equation coefficients policy classes:

    - Must define a _constructor_, which is used to initialize the SDE
      coefficient, C0. Required signature:
      \code{.cpp}
        CoeffPolicyName(
          kw::sde_c3::info::expect::type c3_,
          kw::sde_c4::info::expect::type c4_,
          kw::sde_com1::info::expect::type com1_,
          kw::sde_com2::info::expect::type com2_,
          kw::sde_c3::info::expect::type& c3,
          kw::sde_c4::info::expect::type& c4,
          kw::sde_com1::info::expect::type& com1,
          kw::sde_com2::info::expect::type& com2 )
      \endcode

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::CoeffPolicyType type() noexcept {
          return ctr::CoeffPolicyType::CONST_COEFF;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.

    - Must define the static function _src()_, with the signature
      \code{.cpp}
        void src( tk::real& Som () {}
      \endcode
      which sets its argument _Som_ to the level of the source of the
      dissipation equation. This member function must be static as it is called
      without an object instance.
*/
// *****************************************************************************
#ifndef DissipationCoeffPolicy_h
#define DissipationCoeffPolicy_h

#include <array>

#include <brigand/sequences/list.hpp>

#include "Types.hpp"
#include "SystemComponents.hpp"
#include "Walker/Options/CoeffPolicy.hpp"

namespace walker {

//! Dissipation equation coefficients policy keeping the coefficients constant
class DissipationCoeffConst {

  public:
    //! Constructor: initialize coefficients
    DissipationCoeffConst(
      kw::sde_c3::info::expect::type c3_,
      kw::sde_c4::info::expect::type c4_,
      kw::sde_com1::info::expect::type com1_,
      kw::sde_com2::info::expect::type com2_,
      kw::sde_c3::info::expect::type& c3,
      kw::sde_c4::info::expect::type& c4,
      kw::sde_com1::info::expect::type& com1,
      kw::sde_com2::info::expect::type& com2 );

    //! Update turbulence frequency source (no-op for const-coeff policy)
    static void src( tk::real& ) {}

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONST_COEFF; }
};

//! \brief Dissipation equation coefficients policy keeping the dissipation
//!   rate in a constant statistically stationary state
class DissipationCoeffStationary {

  public:
    //! Constructor: initialize coefficients
    DissipationCoeffStationary(
      kw::sde_c3::info::expect::type c3_,
      kw::sde_c4::info::expect::type c4_,
      kw::sde_com1::info::expect::type com1_,
      kw::sde_com2::info::expect::type com2_,
      kw::sde_c3::info::expect::type& c3,
      kw::sde_c4::info::expect::type& c4,
      kw::sde_com1::info::expect::type& com1,
      kw::sde_com2::info::expect::type& com2 );

    //! Update turbulence frequency source (zero for stationary coeff policy)
    static void src( tk::real& Som ) { Som = 0.0; }

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::STATIONARY; }
};

//! List of all dissipation eq coefficients policies
using DissipationCoeffPolicies =
  brigand::list< DissipationCoeffConst
               , DissipationCoeffStationary >;

} // walker::

#endif // DissipationCoeffPolicy_h
