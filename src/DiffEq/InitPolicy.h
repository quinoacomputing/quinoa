// *****************************************************************************
/*!
  \file      src/DiffEq/InitPolicy.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Initialization policies
  \details   This file defines initialization policy classes. As opposed to
    coefficients policies, see, e.g., DiffEq/BetaCoeffPolicy.h, initialization
    policies are not SDE-specific -- at least at this time.

    General requirements on initialization policy classes:

    - Must define the member function _init_, which is used to do the
      initialization. Required signature:
      \code{.cpp}
        template< class eq >
        static void init( const ctr::InputDeck& deck,
                          const tk::RNG& rng,
                          int stream,
                          tk::Particles& particles,
                          tk::ctr::ncomp_type e,
                          tk::ctr::ncomp_type ncomp,
                          tk::ctr::ncomp_type offset );
      \endcode
      where _deck_ is the input deck from which configuration is read, _rng_ is
      a reference to a random number generator to use, _stream_ is the thread
      (or stream) id, _particles_ denotes the particle properties array to be
      initialized, _e_ is the component index selecting which equation is to be
      initialized in the system, _ncomp_ is the total number of equations in the
      system, and _offset_ is the offset in the particle array at which
      initialization should be done.

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
// *****************************************************************************
#ifndef InitPolicy_h
#define InitPolicy_h

#include <algorithm>

#include <boost/mpl/vector.hpp>

#include "Macro.h"
#include "Types.h"
#include "Particles.h"
#include "Walker/Options/InitPolicy.h"
#include "SystemComponents.h"
#include "RNG.h"

namespace walker {

//! Raw initialization policy: leave memory uninitialized
struct InitRaw {

  //! Initialize particle properties
  template< class eq >
  static void init( const ctr::InputDeck& deck,
                    const tk::RNG& rng,
                    int stream,
                    tk::Particles& particles,
                    tk::ctr::ncomp_type e,
                    tk::ctr::ncomp_type ncomp,
                    tk::ctr::ncomp_type offset )
  {
    IGNORE( deck );
    IGNORE( rng );
    IGNORE( stream );
    IGNORE( particles );
    IGNORE( e );
    IGNORE( ncomp );
    IGNORE( offset );
  }

  static ctr::InitPolicyType type() noexcept
  { return ctr::InitPolicyType::RAW; }
};

//! Zero initialization policy: zero particle properties
struct InitZero {

  //! Initialize particle properties
  template< class eq >
  static void init( const ctr::InputDeck& deck,
                    const tk::RNG& rng,
                    int stream,
                    tk::Particles& particles,
                    tk::ctr::ncomp_type e,
                    tk::ctr::ncomp_type ncomp,
                    tk::ctr::ncomp_type offset )
  {
    IGNORE( deck );
    IGNORE( rng );
    IGNORE( stream );
    IGNORE( e );
    IGNORE( ncomp );
    IGNORE( offset );
    particles.fill( 0.0 );
  }

  static ctr::InitPolicyType type() noexcept
  { return ctr::InitPolicyType::ZERO; }
};

//! Delta initialization policy: put in delta-spikes as the joint PDF
struct InitDelta {

  //! Initialize particle properties
  template< class eq >
  static void init( const ctr::InputDeck& deck,
                    const tk::RNG& rng,
                    int stream,
                    tk::Particles& particles,
                    tk::ctr::ncomp_type e,
                    tk::ctr::ncomp_type ncomp,
                    tk::ctr::ncomp_type offset )
  {
    IGNORE( rng );
    IGNORE( stream );
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
            static_cast< tk::real >( particles.nunk() ) * sc[s+1] );
        // assign sample values
        for (ncomp_t p=0; p<npar; ++p) particles( i+p, c, offset ) = sc[s];
        i += npar;
      }
    }
  }

  static ctr::InitPolicyType type() noexcept
  { return ctr::InitPolicyType::JOINTDELTA; }
};

//! Beta initialization policy: generate samples from a joint beta PDF
struct InitBeta {

  //! Initialize particle properties (zero)
  template< class eq >
  static void init( const ctr::InputDeck& deck,
                    const tk::RNG& rng,
                    int stream,
                    tk::Particles& particles,
                    tk::ctr::ncomp_type e,
                    tk::ctr::ncomp_type ncomp,
                    tk::ctr::ncomp_type offset )
  {
    using ncomp_t = kw::ncomp::info::expect::type;

    const auto& betapdf =
      deck.template get< tag::param, eq, tag::betapdf >().at(e);

    // use only the first ncomp betapdfs if there are more than the equation is
    // configured for
    const ncomp_t size = std::min( ncomp, betapdf.size() );

    for (ncomp_t c=0; c<size; ++c) {
      // get vector of betapdf parameters for component c
      const auto& bc = betapdf[c];

      for (ncomp_t s=0; s<bc.size(); s+=4) {
        // generate beta random numbers for all particles using parameters in bc
        for (ncomp_t p=0; p<particles.nunk(); ++p)
          rng.beta( stream, 1, bc[s], bc[s+1], bc[s+2], bc[s+3],
                    &particles( p, c, offset ) );
      }
    }

  }

  static ctr::InitPolicyType type() noexcept
  { return ctr::InitPolicyType::JOINTBETA; }
};

//! List of all initialization policies
using InitPolicies = boost::mpl::vector< InitRaw
                                       , InitZero
                                       , InitDelta
                                       , InitBeta >;

} // walker::

#endif // InitPolicy_h
