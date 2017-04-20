// *****************************************************************************
/*!
  \file      src/UnitTest/tests/RNG/TestRNGSSE.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for RNG/RNGSSE.h
  \details   Unit tests for RNG/RNGSSE.h
*/
// *****************************************************************************
#ifndef test_RNGSSE_h
#define test_RNGSSE_h

#include "NoWarning/tut.h"

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

#include "RNGSSE.h"

namespace tut {

//! All tests in group inherited from this base
struct RNGSSE_common {};

//! Test group shortcuts
using RNGSSE_group = test_group< RNGSSE_common, MAX_TESTS_IN_GROUP >;
using RNGSSE_object = RNGSSE_group::object;

//! Define test group
static RNGSSE_group RNGSSE( "RNG/RNGSSE" );

//! Test definitions for group

//! Test that constructor throws with zero number of threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 1 >() {
  set_test_name( "constructor throws with zero threads" );

  // Test all possible specializations

  // Attempt to instantiate RNGSSE RNG with zero number of threads
  try {
    tk::RNGSSE< gm19_state, unsigned, gm19_generate_ >
      r( 0, gm19_init_sequence_ );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Need at least one thread" ) !=
              std::string::npos );
  }

  // Attempt to instantiate RNGSSE RNG with zero number of threads
  try {
    tk::RNGSSE< gm29_state, unsigned, gm29_generate_ >
      r( 0, gm29_init_short_sequence_ );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Need at least one thread" ) !=
              std::string::npos );
  }

  // Attempt to instantiate RNGSSE RNG with zero number of threads
  try {
    tk::RNGSSE< gm31_state, unsigned, gm31_generate_ >
      r( 0, gm31_init_short_sequence_ );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Need at least one thread" ) !=
              std::string::npos );
  }

  // Attempt to instantiate RNGSSE RNG with zero number of threads
  try {
    tk::RNGSSE< gm55_state, unsigned long long, gm55_generate_ >
      r( 0, gm55_init_short_sequence_ );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Need at least one thread" ) !=
              std::string::npos );
  }

  // Attempt to instantiate RNGSSE RNG with zero number of threads
  try {
    tk::RNGSSE< gm61_state, unsigned long long, gm61_generate_ >
      r( 0, gm61_init_sequence_ );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Need at least one thread" ) !=
              std::string::npos );
  }

  // Attempt to instantiate RNGSSE RNG with zero number of threads
  try {
    tk::RNGSSE< gq58x1_state, unsigned, gq58x1_generate_ >
      r( 0, gq58x1_init_short_sequence_ );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Need at least one thread" ) !=
              std::string::npos );
  }

  // Attempt to instantiate RNGSSE RNG with zero number of threads
  try {
    tk::RNGSSE< gq58x3_state, unsigned, gq58x3_generate_ >
      r( 0, gq58x3_init_short_sequence_ );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Need at least one thread" ) !=
              std::string::npos );
  }

  // Attempt to instantiate RNGSSE RNG with zero number of threads
  try {
    tk::RNGSSE< gq58x4_state, unsigned, gq58x4_generate_ >
      r( 0, gq58x4_init_short_sequence_ );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Need at least one thread" ) !=
              std::string::npos );
  }

  // Attempt to instantiate RNGSSE RNG with zero number of threads
  try {
    tk::RNGSSE< mt19937_state, unsigned long long, mt19937_generate_ >
      r( 0, mt19937_init_sequence_ );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Need at least one thread" ) !=
              std::string::npos );
  }

  // Attempt to instantiate RNGSSE RNG with zero number of threads
  try {
    tk::RNGSSE< lfsr113_state, unsigned long long, lfsr113_generate_ >
      r( 0, lfsr113_init_long_sequence_ );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Need at least one thread" ) !=
              std::string::npos );
  }

  // Attempt to instantiate RNGSSE RNG with zero number of threads
  try {
    tk::RNGSSE< mrg32k3a_state, unsigned long long, mrg32k3a_generate_ >
      r( 0, mrg32k3a_init_sequence_ );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Need at least one thread" ) !=
              std::string::npos );
  }
}

//! Test Gaussian generator statistics from gm19 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 2 >() {
  set_test_name( "Gaussian gm19 from 4 emulated streams" );

  tk::RNGSSE< gm19_state, unsigned, gm19_generate_ >
    r( 4, gm19_init_sequence_ );
  RNG_common::test_gaussian( r );
}

//! Test Gaussian generator statistics from gm29 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 3 >() {
  set_test_name( "Gaussian gm29 from 4 emulated streams" );

  tk::RNGSSE< gm29_state, unsigned, gm29_generate_ >
    r( 4, gm29_init_short_sequence_ );
  RNG_common::test_gaussian( r );
}

//! Test Gaussian generator statistics from gm31 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 4 >() {
  set_test_name( "Gaussian gm29 from 4 emulated streams" );

  tk::RNGSSE< gm31_state, unsigned, gm31_generate_ >
    r( 4, gm31_init_short_sequence_ );
  RNG_common::test_gaussian( r );
}

//! Test Gaussian generator statistics from gm55 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 5 >() {
  set_test_name( "Gaussian gm55 from 4 emulated streams" );

  tk::RNGSSE< gm55_state, unsigned long long, gm55_generate_ >
    r( 4, gm55_init_short_sequence_ );
  RNG_common::test_gaussian( r );
}

//! Test Gaussian generator statistics from gm61 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 6 >() {
  set_test_name( "Gaussian gm61 from 4 emulated streams" );

  tk::RNGSSE< gm61_state, unsigned long long, gm61_generate_ >
    r( 4, gm61_init_sequence_ );
  RNG_common::test_gaussian( r );
}

//! Test Gaussian generator statistics from gq58.1 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 7 >() {
  set_test_name( "Gaussian gq58.1 from 4 emulated streams" );

  tk::RNGSSE< gq58x1_state, unsigned, gq58x1_generate_ >
    r( 4, gq58x1_init_short_sequence_ );
  RNG_common::test_gaussian( r );
}

//! Test Gaussian generator statistics from gq58.3 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 8 >() {
  set_test_name( "Gaussian gq58.3 from 4 emulated streams" );

  tk::RNGSSE< gq58x3_state, unsigned, gq58x3_generate_ >
    r( 4, gq58x3_init_short_sequence_ );
  RNG_common::test_gaussian( r );
}

//! Test Gaussian generator statistics from gq58.4 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 9 >() {
  set_test_name( "Gaussian gq58.4 from 4 emulated streams" );

  tk::RNGSSE< gq58x4_state, unsigned, gq58x4_generate_ >
    r( 4, gq58x4_init_short_sequence_ );
  RNG_common::test_gaussian( r );
}

//! Test Gaussian generator statistics from mt19937 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 10 >() {
  set_test_name( "Gaussian mt19937 from 4 emulated streams" );

  tk::RNGSSE< mt19937_state, unsigned long long, mt19937_generate_ >
    r( 4, mt19937_init_sequence_ );
  RNG_common::test_gaussian( r );
}

//! Test Gaussian generator statistics from lfsr113 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 11 >() {
  set_test_name( "Gaussian lfsr113 from 4 emulated streams" );

  tk::RNGSSE< lfsr113_state, unsigned long long, lfsr113_generate_ >
    r( 4, lfsr113_init_long_sequence_ );
  RNG_common::test_gaussian( r );
}

//! Test Gaussian generator statistics from mrg32k3a using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 12 >() {
  set_test_name( "Gaussian mrg32k3a from 4 emulated streams" );

  tk::RNGSSE< mrg32k3a_state, unsigned long long, mrg32k3a_generate_ >
    r( 4, mrg32k3a_init_sequence_ );
  RNG_common::test_gaussian( r );
}

//! Test beta generator statistics from gm19 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 13 >() {
  set_test_name( "beta gm19 from 4 emulated streams" );

  tk::RNGSSE< gm19_state, unsigned, gm19_generate_ >
    r( 4, gm19_init_sequence_ );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from gm29 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 14 >() {
  set_test_name( "beta gm29 from 4 emulated streams" );

  tk::RNGSSE< gm29_state, unsigned, gm29_generate_ >
    r( 4, gm29_init_short_sequence_ );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from gm31 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 15 >() {
  set_test_name( "beta gm29 from 4 emulated streams" );

  tk::RNGSSE< gm31_state, unsigned, gm31_generate_ >
    r( 4, gm31_init_short_sequence_ );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from gm55 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 16 >() {
  set_test_name( "beta gm55 from 4 emulated streams" );

  tk::RNGSSE< gm55_state, unsigned long long, gm55_generate_ >
    r( 4, gm55_init_short_sequence_ );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from gm61 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 17 >() {
  set_test_name( "beta gm61 from 4 emulated streams" );

  tk::RNGSSE< gm61_state, unsigned long long, gm61_generate_ >
    r( 4, gm61_init_sequence_ );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from gq58.1 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 18 >() {
  set_test_name( "beta gq58.1 from 4 emulated streams" );

  tk::RNGSSE< gq58x1_state, unsigned, gq58x1_generate_ >
    r( 4, gq58x1_init_short_sequence_ );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from gq58.3 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 19 >() {
  set_test_name( "beta gq58.3 from 4 emulated streams" );

  tk::RNGSSE< gq58x3_state, unsigned, gq58x3_generate_ >
    r( 4, gq58x3_init_short_sequence_ );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from gq58.4 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 20 >() {
  set_test_name( "beta gq58.4 from 4 emulated streams" );

  tk::RNGSSE< gq58x4_state, unsigned, gq58x4_generate_ >
    r( 4, gq58x4_init_short_sequence_ );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from mt19937 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 21 >() {
  set_test_name( "beta mt19937 from 4 emulated streams" );

  tk::RNGSSE< mt19937_state, unsigned long long, mt19937_generate_ >
    r( 4, mt19937_init_sequence_ );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from lfsr113 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 22 >() {
  set_test_name( "beta lfsr113 from 4 emulated streams" );

  tk::RNGSSE< lfsr113_state, unsigned long long, lfsr113_generate_ >
    r( 4, lfsr113_init_long_sequence_ );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from mrg32k3a using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 23 >() {
  set_test_name( "beta mrg32k3a from 4 emulated streams" );

  tk::RNGSSE< mrg32k3a_state, unsigned long long, mrg32k3a_generate_ >
    r( 4, mrg32k3a_init_sequence_ );
  RNG_common::test_beta( r );
}

//! Test uniform generator statistics from gm19 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 24 >() {
  set_test_name( "uniform gm19 from 4 emulated streams" );

  tk::RNGSSE< gm19_state, unsigned, gm19_generate_ >
    r( 4, gm19_init_sequence_ );
  RNG_common::test_uniform( r );
}

//! Test uniform generator statistics from gm29 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 25 >() {
  set_test_name( "uniform gm29 from 4 emulated streams" );

  tk::RNGSSE< gm29_state, unsigned, gm29_generate_ >
    r( 4, gm29_init_short_sequence_ );
  RNG_common::test_uniform( r );
}

//! Test uniform generator statistics from gm31 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 26 >() {
  set_test_name( "uniform gm29 from 4 emulated streams" );

  tk::RNGSSE< gm31_state, unsigned, gm31_generate_ >
    r( 4, gm31_init_short_sequence_ );
  RNG_common::test_uniform( r );
}

//! Test uniform generator statistics from gm55 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 27 >() {
  set_test_name( "uniform gm55 from 4 emulated streams" );

  tk::RNGSSE< gm55_state, unsigned long long, gm55_generate_ >
    r( 4, gm55_init_short_sequence_ );
  RNG_common::test_uniform( r );
}

//! Test uniform generator statistics from gm61 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 28 >() {
  set_test_name( "uniform gm61 from 4 emulated streams" );

  tk::RNGSSE< gm61_state, unsigned long long, gm61_generate_ >
    r( 4, gm61_init_sequence_ );
  RNG_common::test_uniform( r );
}

//! Test uniform generator statistics from gq58.1 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 29 >() {
  set_test_name( "uniform gq58.1 from 4 emulated streams" );

  tk::RNGSSE< gq58x1_state, unsigned, gq58x1_generate_ >
    r( 4, gq58x1_init_short_sequence_ );
  RNG_common::test_uniform( r );
}

//! Test uniform generator statistics from gq58.3 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 30 >() {
  set_test_name( "uniform gq58.3 from 4 emulated streams" );

  tk::RNGSSE< gq58x3_state, unsigned, gq58x3_generate_ >
    r( 4, gq58x3_init_short_sequence_ );
  RNG_common::test_uniform( r );
}

//! Test uniform generator statistics from gq58.4 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 31 >() {
  set_test_name( "uniform gq58.4 from 4 emulated streams" );

  tk::RNGSSE< gq58x4_state, unsigned, gq58x4_generate_ >
    r( 4, gq58x4_init_short_sequence_ );
  RNG_common::test_uniform( r );
}

//! Test uniform generator statistics from mt19937 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 32 >() {
  set_test_name( "uniform mt19937 from 4 emulated streams" );

  tk::RNGSSE< mt19937_state, unsigned long long, mt19937_generate_ >
    r( 4, mt19937_init_sequence_ );
  RNG_common::test_uniform( r );
}

//! Test uniform generator statistics from lfsr113 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 33 >() {
  set_test_name( "uniform lfsr113 from 4 emulated streams" );

  tk::RNGSSE< lfsr113_state, unsigned long long, lfsr113_generate_ >
    r( 4, lfsr113_init_long_sequence_ );
  RNG_common::test_uniform( r );
}

//! Test uniform generator statistics from mrg32k3a using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 34 >() {
  set_test_name( "uniform mrg32k3a from 4 emulated streams" );

  tk::RNGSSE< mrg32k3a_state, unsigned long long, mrg32k3a_generate_ >
    r( 4, mrg32k3a_init_sequence_ );
  RNG_common::test_uniform( r );
}

//! Test copy constructor for gm19
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 35 >() {
  set_test_name( "copy constructor with gm19" );

  tk::RNGSSE< gm19_state, unsigned, gm19_generate_ >
    r( 4, gm19_init_sequence_ );
  RNG_common::test_copy_ctor( r );
}

//! Test copy constructor for gm29 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 36 >() {
  set_test_name( "copy constructor with gm29" );

  tk::RNGSSE< gm29_state, unsigned, gm29_generate_ >
    r( 4, gm29_init_short_sequence_ );
  RNG_common::test_copy_ctor( r );
}

//! Test copy constructor for gm31 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 37 >() {
  set_test_name( "copy constructor with gm29" );

  tk::RNGSSE< gm31_state, unsigned, gm31_generate_ >
    r( 4, gm31_init_short_sequence_ );
  RNG_common::test_copy_ctor( r );
}

//! Test copy constructor for gm55 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 38 >() {
  set_test_name( "copy constructor with gm55" );

  tk::RNGSSE< gm55_state, unsigned long long, gm55_generate_ >
    r( 4, gm55_init_short_sequence_ );
  RNG_common::test_copy_ctor( r );
}

//! Test copy constructor for gm61 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 39 >() {
  set_test_name( "copy constructor with gm61" );

  tk::RNGSSE< gm61_state, unsigned long long, gm61_generate_ >
    r( 4, gm61_init_sequence_ );
  RNG_common::test_copy_ctor( r );
}

//! Test copy constructor for gq58.1 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 40 >() {
  set_test_name( "copy constructor with gq58.1" );

  tk::RNGSSE< gq58x1_state, unsigned, gq58x1_generate_ >
    r( 4, gq58x1_init_short_sequence_ );
  RNG_common::test_copy_ctor( r );
}

//! Test copy constructor for gq58.3 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 41 >() {
  set_test_name( "copy constructor with gq58.3" );

  tk::RNGSSE< gq58x3_state, unsigned, gq58x3_generate_ >
    r( 4, gq58x3_init_short_sequence_ );
  RNG_common::test_copy_ctor( r );
}

//! Test copy constructor for gq58.4 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 42 >() {
  set_test_name( "copy constructor with gq58.4" );

  tk::RNGSSE< gq58x4_state, unsigned, gq58x4_generate_ >
    r( 4, gq58x4_init_short_sequence_ );
  RNG_common::test_copy_ctor( r );
}

//! Test copy constructor for mt19937 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 43 >() {
  set_test_name( "copy constructor with mt19937" );

  tk::RNGSSE< mt19937_state, unsigned long long, mt19937_generate_ >
    r( 4, mt19937_init_sequence_ );
  RNG_common::test_copy_ctor( r );
}

//! Test copy constructor for lfsr113 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 44 >() {
  set_test_name( "copy constructor with lfsr113" );

  tk::RNGSSE< lfsr113_state, unsigned long long, lfsr113_generate_ >
    r( 4, lfsr113_init_long_sequence_ );
  RNG_common::test_copy_ctor( r );
}

//! Test copy constructor for mrg32k3a using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 45 >() {
  set_test_name( "copy constructor with mrg32k3a" );

  tk::RNGSSE< mrg32k3a_state, unsigned long long, mrg32k3a_generate_ >
    r( 4, mrg32k3a_init_sequence_ );
  RNG_common::test_copy_ctor( r );
}

//! Test move constructor for gm19
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 46 >() {
  set_test_name( "move constructor with gm19" );

  tk::RNGSSE< gm19_state, unsigned, gm19_generate_ >
    r( 4, gm19_init_sequence_ );
  RNG_common::test_move_ctor( r );
}

//! Test move constructor for gm29 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 47 >() {
  set_test_name( "move constructor with gm29" );

  tk::RNGSSE< gm29_state, unsigned, gm29_generate_ >
    r( 4, gm29_init_short_sequence_ );
  RNG_common::test_move_ctor( r );
}

//! Test move constructor for gm31 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 48 >() {
  set_test_name( "move constructor with gm29" );

  tk::RNGSSE< gm31_state, unsigned, gm31_generate_ >
    r( 4, gm31_init_short_sequence_ );
  RNG_common::test_move_ctor( r );
}

//! Test move constructor for gm55 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 49 >() {
  set_test_name( "move constructor with gm55" );

  tk::RNGSSE< gm55_state, unsigned long long, gm55_generate_ >
    r( 4, gm55_init_short_sequence_ );
  RNG_common::test_move_ctor( r );
}

//! Test move constructor for gm61 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 50 >() {
  set_test_name( "move constructor with gm61" );

  tk::RNGSSE< gm61_state, unsigned long long, gm61_generate_ >
    r( 4, gm61_init_sequence_ );
  RNG_common::test_move_ctor( r );
}

//! Test move constructor for gq58.1 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 51 >() {
  set_test_name( "move constructor with gq58.1" );

  tk::RNGSSE< gq58x1_state, unsigned, gq58x1_generate_ >
    r( 4, gq58x1_init_short_sequence_ );
  RNG_common::test_move_ctor( r );
}

//! Test move constructor for gq58.3 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 52 >() {
  set_test_name( "move constructor with gq58.3" );

  tk::RNGSSE< gq58x3_state, unsigned, gq58x3_generate_ >
    r( 4, gq58x3_init_short_sequence_ );
  RNG_common::test_move_ctor( r );
}

//! Test move constructor for gq58.4 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 53 >() {
  set_test_name( "move constructor with gq58.4" );

  tk::RNGSSE< gq58x4_state, unsigned, gq58x4_generate_ >
    r( 4, gq58x4_init_short_sequence_ );
  RNG_common::test_move_ctor( r );
}

//! Test move constructor for mt19937 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 54 >() {
  set_test_name( "move constructor with mt19937" );

  tk::RNGSSE< mt19937_state, unsigned long long, mt19937_generate_ >
    r( 4, mt19937_init_sequence_ );
  RNG_common::test_move_ctor( r );
}

//! Test move constructor for lfsr113 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 55 >() {
  set_test_name( "move constructor with lfsr113" );

  tk::RNGSSE< lfsr113_state, unsigned long long, lfsr113_generate_ >
    r( 4, lfsr113_init_long_sequence_ );
  RNG_common::test_move_ctor( r );
}

//! Test move constructor for mrg32k3a using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 56 >() {
  set_test_name( "move constructor with mrg32k3a" );

  tk::RNGSSE< mrg32k3a_state, unsigned long long, mrg32k3a_generate_ >
    r( 4, mrg32k3a_init_sequence_ );
  RNG_common::test_move_ctor( r );
}

//! Test copy assignment for gm19
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 57 >() {
  set_test_name( "copy assignment with gm19" );

  tk::RNGSSE< gm19_state, unsigned, gm19_generate_ >
    r( 4, gm19_init_sequence_ );
  RNG_common::test_copy_assignment( r );
}

//! Test copy assignment for gm29 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 58 >() {
  set_test_name( "copy assignment with gm29" );

  tk::RNGSSE< gm29_state, unsigned, gm29_generate_ >
    r( 4, gm29_init_short_sequence_ );
  RNG_common::test_copy_assignment( r );
}

//! Test copy assignment for gm31 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 59 >() {
  set_test_name( "copy assignment with gm29" );

  tk::RNGSSE< gm31_state, unsigned, gm31_generate_ >
    r( 4, gm31_init_short_sequence_ );
  RNG_common::test_copy_assignment( r );
}

//! Test copy assignment for gm55 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 60 >() {
  set_test_name( "copy assignment with gm55" );

  tk::RNGSSE< gm55_state, unsigned long long, gm55_generate_ >
    r( 4, gm55_init_short_sequence_ );
  RNG_common::test_copy_assignment( r );
}

//! Test copy assignment for gm61 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 61 >() {
  set_test_name( "copy assignment with gm61" );

  tk::RNGSSE< gm61_state, unsigned long long, gm61_generate_ >
    r( 4, gm61_init_sequence_ );
  RNG_common::test_copy_assignment( r );
}

//! Test copy assignment for gq58.1 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 62 >() {
  set_test_name( "copy assignment with gq58.1" );

  tk::RNGSSE< gq58x1_state, unsigned, gq58x1_generate_ >
    r( 4, gq58x1_init_short_sequence_ );
  RNG_common::test_copy_assignment( r );
}

//! Test copy assignment for gq58.3 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 63 >() {
  set_test_name( "copy assignment with gq58.3" );

  tk::RNGSSE< gq58x3_state, unsigned, gq58x3_generate_ >
    r( 4, gq58x3_init_short_sequence_ );
  RNG_common::test_copy_assignment( r );
}

//! Test copy assignment for gq58.4 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 64 >() {
  set_test_name( "copy assignment with gq58.4" );

  tk::RNGSSE< gq58x4_state, unsigned, gq58x4_generate_ >
    r( 4, gq58x4_init_short_sequence_ );
  RNG_common::test_copy_assignment( r );
}

//! Test copy assignment for mt19937 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 65 >() {
  set_test_name( "copy assignment with mt19937" );

  tk::RNGSSE< mt19937_state, unsigned long long, mt19937_generate_ >
    r( 4, mt19937_init_sequence_ );
  RNG_common::test_copy_assignment( r );
}

//! Test copy assignment for lfsr113 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 66 >() {
  set_test_name( "copy assignment with lfsr113" );

  tk::RNGSSE< lfsr113_state, unsigned long long, lfsr113_generate_ >
    r( 4, lfsr113_init_long_sequence_ );
  RNG_common::test_copy_assignment( r );
}

//! Test copy assignment for mrg32k3a using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 67 >() {
  set_test_name( "copy assignment with mrg32k3a" );

  tk::RNGSSE< mrg32k3a_state, unsigned long long, mrg32k3a_generate_ >
    r( 4, mrg32k3a_init_sequence_ );
  RNG_common::test_copy_assignment( r );
}

//! Test move assignment for gm19
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 68 >() {
  set_test_name( "move assignment with gm19" );

  tk::RNGSSE< gm19_state, unsigned, gm19_generate_ >
    r( 4, gm19_init_sequence_ );
  RNG_common::test_move_assignment( r );
}

//! Test move assignment for gm29 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 69 >() {
  set_test_name( "move assignment with gm29" );

  tk::RNGSSE< gm29_state, unsigned, gm29_generate_ >
    r( 4, gm29_init_short_sequence_ );
  RNG_common::test_move_assignment( r );
}

//! Test move assignment for gm31 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 70 >() {
  set_test_name( "move assignment with gm29" );

  tk::RNGSSE< gm31_state, unsigned, gm31_generate_ >
    r( 4, gm31_init_short_sequence_ );
  RNG_common::test_move_assignment( r );
}

//! Test move assignment for gm55 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 71 >() {
  set_test_name( "move assignment with gm55" );

  tk::RNGSSE< gm55_state, unsigned long long, gm55_generate_ >
    r( 4, gm55_init_short_sequence_ );
  RNG_common::test_move_assignment( r );
}

//! Test move assignment for gm61 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 72 >() {
  set_test_name( "move assignment with gm61" );

  tk::RNGSSE< gm61_state, unsigned long long, gm61_generate_ >
    r( 4, gm61_init_sequence_ );
  RNG_common::test_move_assignment( r );
}

//! Test move assignment for gq58.1 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 73 >() {
  set_test_name( "move assignment with gq58.1" );

  tk::RNGSSE< gq58x1_state, unsigned, gq58x1_generate_ >
    r( 4, gq58x1_init_short_sequence_ );
  RNG_common::test_move_assignment( r );
}

//! Test move assignment for gq58.3 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 74 >() {
  set_test_name( "move assignment with gq58.3" );

  tk::RNGSSE< gq58x3_state, unsigned, gq58x3_generate_ >
    r( 4, gq58x3_init_short_sequence_ );
  RNG_common::test_move_assignment( r );
}

//! Test move assignment for gq58.4 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 75 >() {
  set_test_name( "move assignment with gq58.4" );

  tk::RNGSSE< gq58x4_state, unsigned, gq58x4_generate_ >
    r( 4, gq58x4_init_short_sequence_ );
  RNG_common::test_move_assignment( r );
}

//! Test move assignment for mt19937 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 76 >() {
  set_test_name( "move assignment with mt19937" );

  tk::RNGSSE< mt19937_state, unsigned long long, mt19937_generate_ >
    r( 4, mt19937_init_sequence_ );
  RNG_common::test_move_assignment( r );
}

//! Test move assignment for lfsr113 using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 77 >() {
  set_test_name( "move assignment with lfsr113" );

  tk::RNGSSE< lfsr113_state, unsigned long long, lfsr113_generate_ >
    r( 4, lfsr113_init_long_sequence_ );
  RNG_common::test_move_assignment( r );
}

//! Test move assignment for mrg32k3a using multiple threads
//! \author J. Bakosi
template<> template<>
void RNGSSE_object::test< 78 >() {
  set_test_name( "move assignment with mrg32k3a" );

  tk::RNGSSE< mrg32k3a_state, unsigned long long, mrg32k3a_generate_ >
    r( 4, mrg32k3a_init_sequence_ );
  RNG_common::test_move_assignment( r );
}

} // tut::

#endif // test_RNGSSE_h
