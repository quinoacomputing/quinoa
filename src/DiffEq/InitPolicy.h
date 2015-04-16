//******************************************************************************
/*!
  \file      src/DiffEq/InitPolicy.h
  \author    J. Bakosi
  \date      Thu 12 Mar 2015 09:38:54 PM MDT
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
#include <algorithm>

#include <boost/mpl/vector.hpp>

#include <Macro.h>
#include <Types.h>
#include <ParticleProperties.h>
#include <Walker/Options/InitPolicy.h>

namespace walker {

//! Raw initialization policy: leave memory uninitialized
struct InitRaw {

  //! Initialize particle properties (raw: no-op)
  template< class eq, class InputDeck >
  static void init( const InputDeck& deck,
                    tk::ParProps& particles,
                    tk::ctr::ncomp_type e,
                    tk::ctr::ncomp_type ncomp,
                    tk::ctr::ncomp_type offset ) {}

  static ctr::InitPolicyType type() noexcept
  { return ctr::InitPolicyType::RAW; }
};

//! Zero initialization policy: zero particle properties
struct InitZero {

  //! Initialize particle properties (zero)
  template< class eq, class InputDeck >
  static void init( const InputDeck& deck,
                    tk::ParProps& particles,
                    tk::ctr::ncomp_type e,
                    tk::ctr::ncomp_type ncomp,
                    tk::ctr::ncomp_type offset )
  {
    std::memset( particles.ptr(), 0, particles.size()*sizeof(tk::real) );
  }

  static ctr::InitPolicyType type() noexcept
  { return ctr::InitPolicyType::ZERO; }
};

//! Delta initialization policy: put in delta-spikes as the joint PDF
struct InitDelta {

  //! Initialize particle properties (zero)
  template< class eq, class InputDeck >
  static void init( const InputDeck& deck,
                    tk::ParProps& particles,
                    tk::ctr::ncomp_type e,
                    tk::ctr::ncomp_type ncomp,
                    tk::ctr::ncomp_type offset )
  {
    using ncomp_t = kw::ncomp::info::expect::type;

    const auto& spike = deck.template get< tag::param, eq, tag::spike >().at(e);

    // use only the first ncomp spikes if there are more than the equation is
    // configured for
    const ncomp_t size = std::min( ncomp, spike.size() );

    for (ncomp_t c=0; c<size; ++c) {
      const auto& sc = spike[c];        // vector of spikes for component c

      ncomp_t i = 0;
      for (ncomp_t s=0; s<sc.size(); s+=2) {
        // compute number of samples to be set at relative probability height
        const auto npar =
          static_cast< ncomp_t >(
            static_cast< tk::real >( particles.npar() ) * sc[s+1] );
        // assign sample values
        for (ncomp_t p=0; p<npar; ++p) particles( i+p, c, offset ) = sc[s];
        i += npar;
      }
    }
  }

  static ctr::InitPolicyType type() noexcept
  { return ctr::InitPolicyType::DELTA; }
};

//! List of all initialization policies
using InitPolicies = boost::mpl::vector< InitRaw, InitZero, InitDelta >;

} // walker::

#endif // InitPolicy_h
