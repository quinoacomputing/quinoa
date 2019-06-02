// *****************************************************************************
/*!
  \file      tests/unit/RNG/TestRNG.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for RNG/RNG.hpp
  \details   Unit tests for RNG/RNG.hpp
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "QuinoaConfig.hpp"

#include "NoWarning/threefry.hpp"
#include "NoWarning/philox.hpp"

#ifdef HAS_MKL
  #include <mkl_vsl_types.h>
  #include "MKLRNG.hpp"
#endif

#ifdef HAS_RNGSSE2
  #include <gm19.h>
  #include <gm29.h>
  #include <gm31.h>
  #include <gm55.h>
  #include <gm61.h>
  #include <gq58x1.h>
  #include <gq58x3.h>
  #include <gq58x4.h>
  #include <mt19937.h>
  #include <lfsr113.h>
  #include <mrg32k3a.h>
  #include "RNGSSE.hpp"
#endif

#include "Random123.hpp"
#include "TestRNG.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! Test group shortcuts
using RNG_group = test_group< RNG_common, MAX_TESTS_IN_GROUP >;
using RNG_object = RNG_group::object;

//! Define test group
static RNG_group RNG( "RNG/RNG" );

//! Test definitions for group

//! Test constructor taking an object modeling Concept in tk::RNG
template<> template<>
void RNG_object::test< 1 >() {
  set_test_name( "ctor( rng() ) & nthreads()" );
  for (const auto& r : rngs)
    ensure_equals( "nthreads() via polymorphic tk::RNG call incorrect",
                   r.nthreads(), 4 );
}

//! \brief Test constructor taking a function pointer to a constructor of an
//!   object modeling Concept in tk::RNG
template<> template<>
void RNG_object::test< 2 >() {
  set_test_name( "ctor( std::function<rng>(rng) )" );

  std::vector< tk::RNG > v;

  #ifdef HAS_MKL
  add< tk::MKLRNG >( v, 4, VSL_BRNG_MCG31 );
  #endif

  #ifdef HAS_RNGSSE2
  add< tk::RNGSSE< gm19_state, unsigned, gm19_generate_ > >
     ( v, 4, gm19_init_sequence_ );
  add< tk::RNGSSE< gm29_state, unsigned, gm29_generate_ > >
     ( v, 4, gm29_init_short_sequence_ );
  add< tk::RNGSSE< gm31_state, unsigned, gm31_generate_ > >
     ( v, 4, gm31_init_short_sequence_ );
  add< tk::RNGSSE< gm55_state, unsigned long long, gm55_generate_ > >
     ( v, 4, gm55_init_short_sequence_ );
  add< tk::RNGSSE< gm61_state, unsigned long long, gm61_generate_ > >
     ( v, 4, gm61_init_sequence_ );
  add< tk::RNGSSE< gq58x1_state, unsigned, gq58x1_generate_ > >
     ( v, 4, gq58x1_init_short_sequence_ );
  add< tk::RNGSSE< gq58x3_state, unsigned, gq58x3_generate_ > >
     ( v, 4, gq58x3_init_short_sequence_ );
  add< tk::RNGSSE< gq58x4_state, unsigned, gq58x4_generate_ > >
     ( v, 4, gq58x4_init_short_sequence_ );
  add< tk::RNGSSE< mt19937_state, unsigned long long, mt19937_generate_ > >
     ( v, 4, mt19937_init_sequence_ );
  add< tk::RNGSSE< lfsr113_state, unsigned long long, lfsr113_generate_ > >
     ( v, 4, lfsr113_init_long_sequence_ );
  add< tk::RNGSSE< mrg32k3a_state, unsigned long long, mrg32k3a_generate_ > >
     ( v, 4, mrg32k3a_init_sequence_ );
  #endif

  add< tk::Random123< r123::Threefry2x64 > >( v, 4 );
  add< tk::Random123< r123::Philox2x64 > >( v, 4 );

  for (const auto& r : v)
    ensure_equals( "nthreads() via polymorphic tk::RNG call incorrect",
                   r.nthreads(), 4 );
}

//! Test Gaussian generator statistics via polymorphic call in tk::RNG
template<> template<>
void RNG_object::test< 3 >() {
  set_test_name( "Gaussian from 4 emulated streams" );
  for (const auto& r : rngs) test_gaussian( r );
}

//! Test beta generator statistics via polymorphic call in tk::RNG
template<> template<>
void RNG_object::test< 4 >() {
  set_test_name( "beta from 4 emulated streams" );
  for (const auto& r : rngs) test_beta( r );
}

//! Test uniform generator statistics via polymorphic call in tk::RNG
template<> template<>
void RNG_object::test< 5 >() {
  set_test_name( "uniform from 4 emulated streams" );
  for (const auto& r : rngs) test_uniform( r );
}

//! Test copy constructor
template<> template<>
void RNG_object::test< 6 >() {
  set_test_name( "copy constructor" );
  for (const auto& r : rngs) test_copy_ctor( r );
}

//! Test move constructor
template<> template<>
void RNG_object::test< 7 >() {
  set_test_name( "move constructor" );
  for (const auto& r : rngs) test_move_ctor( r );
}

//! Test copy assignment
template<> template<>
void RNG_object::test< 8 >() {
  set_test_name( "copy assignment" );
  for (const auto& r : rngs) test_copy_assignment( r );
}

//! Test move assignment
template<> template<>
void RNG_object::test< 9 >() {
  set_test_name( "move assignment" );
  for (const auto& r : rngs) test_move_assignment( r );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
