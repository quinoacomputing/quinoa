// *****************************************************************************
/*!
  \file      src/DiffEq/DissipationCoeffPolicy.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
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
*/
// *****************************************************************************
#ifndef DissipationCoeffPolicy_h
#define DissipationCoeffPolicy_h

#include <array>

#include <brigand/sequences/list.hpp>

#include "Types.h"
#include "SystemComponents.h"
#include "Walker/Options/CoeffPolicy.h"

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
      kw::sde_com2::info::expect::type& com2 )
    {
      c3 = c3_;
      c4 = c4_;
      com1 = com1_;
      com2 = com2_;
    }

    //! Update turbulence frequency source (no-op for const-coeff policy)
    static void src( tk::real& ) {}

    //! Coefficients policy type accessor
    static ctr::CoeffPolicyType type() noexcept
    { return ctr::CoeffPolicyType::CONST_COEFF; }
};

//! List of all dissipation eq coefficients policies
using DissipationCoeffPolicies = brigand::list< DissipationCoeffConst >;

} // walker::

#endif // DissipationCoeffPolicy_h
