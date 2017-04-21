// *****************************************************************************
/*!
  \file      src/UnitTest/tests/RNG/TestRandom123.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for RNG/Random123.h
  \details   Unit tests for RNG/Random123.h
*/
// *****************************************************************************
#ifndef test_Random123_h
#define test_Random123_h

#include "NoWarning/tut.h"

#include "NoWarning/threefry.h"
#include "NoWarning/philox.h"

#include "Random123.h"

namespace tut {

//! All tests in group inherited from this base
struct Random123_common {};

//! Test group shortcuts
using Random123_group = test_group< Random123_common, MAX_TESTS_IN_GROUP >;
using Random123_object = Random123_group::object;

//! Define test group
static Random123_group Random123( "RNG/Random123" );

//! Test definitions for group

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif

//! Test that constructor throws with zero number of threads
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 1 >() {
  set_test_name( "constructor throws with zero threads" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  // Attempt to instantiate Random123 with zero number of threads
  try {
    tk::Random123< r123::Threefry2x64 > r( 0 );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Need at least one thread" ) !=
              std::string::npos );
  }
  #endif
}

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

//! Test uniform generator statistics from threefry using a single thread
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 2 >() {
  set_test_name( "uniform threefry from a single stream" );

  tk::Random123< r123::Threefry2x64 > r( 1 );
  RNG_common::test_uniform( r );
}

//! Test uniform generator statistics from threefry using multiple threads
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 3 >() {
  set_test_name( "uniform threefry from 4 emulated streams" );

  tk::Random123< r123::Threefry2x64 > r( 4 );
  RNG_common::test_uniform( r );
}

//! Test Gaussian generator statistics from threefry using a single thread
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 4 >() {
  set_test_name( "Gaussian threefry from a single stream" );

  tk::Random123< r123::Threefry2x64 > r( 1 );
  RNG_common::test_gaussian( r );
}

//! Test Gaussian generator statistics from threefry using multiple threads
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 5 >() {
  set_test_name( "Gaussian threefry from 4 emulated streams" );

  tk::Random123< r123::Threefry2x64 > r( 4 );
  RNG_common::test_gaussian( r );
}

//! Test beta generator statistics from threefry using a single thread
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 6 >() {
  set_test_name( "beta threefry from a single stream" );

  tk::Random123< r123::Threefry2x64 > r( 1 );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from threefry using multiple threads
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 7 >() {
  set_test_name( "beta threefry from 4 emulated streams" );

  tk::Random123< r123::Threefry2x64 > r( 4 );
  RNG_common::test_beta( r );
}

//! Test copy constructor for threefry
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 8 >() {
  set_test_name( "copy constructor with threefry" );

  tk::Random123< r123::Threefry2x64 > p( 1 );   // one thread
  RNG_common::test_copy_ctor( p );
  tk::Random123< r123::Threefry2x64 > r( 4 );   // 4 threads
  RNG_common::test_copy_ctor( r );
}

//! Test move constructor for threefry
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 9 >() {
  set_test_name( "move constructor with threefry" );

  tk::Random123< r123::Threefry2x64 > r( 4 );
  RNG_common::test_move_ctor( r );
}

//! Test copy assignment for threefry
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 10 >() {
  set_test_name( "copy assignment with threefry" );

  tk::Random123< r123::Threefry2x64 > p( 1 );   // one thread
  RNG_common::test_copy_assignment( p );
  tk::Random123< r123::Threefry2x64 > r( 4 );   // 4 threads
  RNG_common::test_copy_assignment( r );
}

//! Test move assignment for threefry
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 11 >() {
  set_test_name( "move assignment with threefry" );

  tk::Random123< r123::Threefry2x64 > p( 1 );   // one thread
  RNG_common::test_move_assignment( p );
  tk::Random123< r123::Threefry2x64 > r( 4 );   // 4 threads
  RNG_common::test_move_assignment( r );
}

//! Test uniform generator statistics from philox using a single thread
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 12 >() {
  set_test_name( "uniform philox from a single stream" );

  tk::Random123< r123::Philox2x64 > r( 1 );
  RNG_common::test_uniform( r );
}

//! Test uniform generator statistics from philox using multiple threads
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 13 >() {
  set_test_name( "uniform philox from 4 emulated streams" );

  tk::Random123< r123::Philox2x64 > r( 4 );
  RNG_common::test_uniform( r );
}

//! Test Gaussian generator statistics from philox using a single thread
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 14 >() {
  set_test_name( "Gaussian philox from a single stream" );

  tk::Random123< r123::Philox2x64 > r( 1 );
  RNG_common::test_gaussian( r );
}

//! Test Gaussian generator statistics from philox using multiple threads
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 15 >() {
  set_test_name( "Gaussian philox from 4 emulated streams" );

  tk::Random123< r123::Philox2x64 > r( 4 );
  RNG_common::test_gaussian( r );
}

//! Test beta generator statistics from philox using a single thread
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 16 >() {
  set_test_name( "beta philox from a single stream" );

  tk::Random123< r123::Philox2x64 > r( 1 );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from philox using multiple threads
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 17 >() {
  set_test_name( "beta philox from 4 emulated streams" );

  tk::Random123< r123::Philox2x64 > r( 4 );
  RNG_common::test_beta( r );
}

//! Test copy constructor for philox
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 18 >() {
  set_test_name( "copy constructor with philox" );

  tk::Random123< r123::Philox2x64 > p( 1 );   // one thread
  RNG_common::test_copy_ctor( p );
  tk::Random123< r123::Philox2x64 > r( 4 );   // 4 threads
  RNG_common::test_copy_ctor( r );
}

//! Test move constructor for philox
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 19 >() {
  set_test_name( "move constructor with philox" );

  tk::Random123< r123::Philox2x64 > r( 4 );
  RNG_common::test_move_ctor( r );
}

//! Test copy assignment for philox
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 20 >() {
  set_test_name( "copy assignment with philox" );

  tk::Random123< r123::Philox2x64 > p( 1 );   // one thread
  RNG_common::test_copy_assignment( p );
  tk::Random123< r123::Philox2x64 > r( 4 );   // 4 threads
  RNG_common::test_copy_assignment( r );
}

//! Test move assignment for philox
//! \author J. Bakosi
template<> template<>
void Random123_object::test< 21 >() {
  set_test_name( "move assignment with philox" );

  tk::Random123< r123::Philox2x64 > p( 1 );   // one thread
  RNG_common::test_move_assignment( p );
  tk::Random123< r123::Philox2x64 > r( 4 );   // 4 threads
  RNG_common::test_move_assignment( r );
}

} // tut::

#endif // test_Random123_h
