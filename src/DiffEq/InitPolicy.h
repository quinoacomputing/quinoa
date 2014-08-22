//******************************************************************************
/*!
  \file      src/DiffEq/InitPolicy.h
  \author    J. Bakosi
  \date      Fri 15 Aug 2014 09:52:57 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Initialization policies
  \details   Initialization policies
*/
//******************************************************************************
#ifndef InitPolicy_h
#define InitPolicy_h

#include <cstring>

#include <boost/mpl/vector.hpp>

#include <Macro.h>
#include <Types.h>
#include <ParticleProperties.h>
#include <Quinoa/Options/InitPolicy.h>

namespace quinoa {

//! Raw initialization policy: leave memory uninitialized
struct InitRaw {

  //! Constructor: default for accessing policy name, type, etc.
  InitRaw() = default;
  //! Constructor: initialize particle properties (raw: no-op)
  InitRaw( ParProps& particles ) {}

  std::string policy() const noexcept
  { return ctr::InitPolicy().name( ctr::InitPolicyType::RAW ); }

  ctr::InitPolicyType type() const noexcept
  { return ctr::InitPolicyType::RAW; }
};

//! Zero initialization policy: zero particle properties
struct InitZero {

  //! Constructor: default for accessing policy name, type, etc.
  InitZero() = default;
  //! Constructor: initialize particle properties (zero)
  InitZero( ParProps& particles )
  { memset( particles.ptr(), 0, particles.size()*sizeof(tk::real) ); }

  std::string policy() const noexcept
  { return ctr::InitPolicy().name( ctr::InitPolicyType::ZERO ); }

  ctr::InitPolicyType type() const noexcept
  { return ctr::InitPolicyType::ZERO; }
};

//! List of all initialization policies
using InitPolicies = boost::mpl::vector< InitRaw, InitZero >;

} // quinoa::

#endif // InitPolicy_h
