//******************************************************************************
/*!
  \file      src/DiffEq/InitPolicy.h
  \author    J. Bakosi
  \date      Mon 26 Jan 2015 01:23:51 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Initialization policies
  \details   This file defines initialization policy classes. As opposed to
    coefficients policies, see, e.g., DiffEq/BetaCoeffPolicy.h, initialization
    policies are not SDE-specific -- at least at this time.

    General requirements on initialization policy classes:

    - Must define a _constructor_, which is used to do the initialization.
      Required signature:
      \code{.cpp}
        InitPolicyName( ParProps& particles )
      \endcode
      where particles denotes the particle properties array to be initialized.

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::InitPolicyType type() noexcept {
          return ctr::InitPolicyType::RAW;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for initialization policies.
*/
//******************************************************************************
#ifndef InitPolicy_h
#define InitPolicy_h

#include <cstring>

#include <boost/mpl/vector.hpp>

#include <Macro.h>
#include <Types.h>
#include <ParticleProperties.h>
#include <Options/InitPolicy.h>

namespace tk {

//! Raw initialization policy: leave memory uninitialized
struct InitRaw {

  //! Constructor: initialize particle properties (raw: no-op)
  InitRaw( ParProps& particles ) {}

  static ctr::InitPolicyType type() noexcept
  { return ctr::InitPolicyType::RAW; }
};

//! Zero initialization policy: zero particle properties
struct InitZero {

  //! Constructor: initialize particle properties (zero)
  InitZero( ParProps& particles )
  { memset( particles.ptr(), 0, particles.size()*sizeof(real) ); }

  static ctr::InitPolicyType type() noexcept
  { return ctr::InitPolicyType::ZERO; }
};

//! List of all initialization policies
using InitPolicies = boost::mpl::vector< InitRaw, InitZero >;

} // tk::

#endif // InitPolicy_h
