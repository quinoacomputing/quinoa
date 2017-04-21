// *****************************************************************************
/*!
  \file      src/UnitTest/tests/RNG/TestMKLRNG.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for RNG/MKLRNG.h
  \details   Unit tests for RNG/MKLRNG.h
*/
// *****************************************************************************
#ifndef test_MKLRNG_h
#define test_MKLRNG_h

#include "NoWarning/tut.h"

#include "MKLRNG.h"

namespace tut {

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif

//! All tests in group inherited from this base
struct MKLRNG_common {};

//! Test group shortcuts
using MKLRNG_group = test_group< MKLRNG_common, MAX_TESTS_IN_GROUP >;
using MKLRNG_object = MKLRNG_group::object;

//! Define test group
static MKLRNG_group MKLRNG( "RNG/MKLRNG" );

//! Test definitions for group

//! Test that constructor throws with zero number of threads
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 1 >() {
  set_test_name( "constructor throws with zero threads" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  // Attempt to instantiate MKLRNG with zero number of threads
  try {
    tk::MKLRNG r( 0 );
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

//! Test that constructor throws with crap data for basic RNG
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 2 >() {
  set_test_name( "constructor throws with crap brng" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  // Attempt to instantiate MKLRNG with crap brng
  try {
    tk::MKLRNG r( 1, -1 );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Asserts skipped in RELEASE mode
    // since we feed crap to the 2nd argument above, we expect an exception,
    // test if the correct exception is thrown
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Basic RNG MKL parameter" ) !=
              std::string::npos ||
            std::string( e.what() ).find( "MKL VSL Error Code" ) !=
              std::string::npos );
  }
  #endif
}

//! \brief Test that constructor throws with multiple threads if basic RNG does
//!   not support leapfrogging
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 3 >() {
  set_test_name( "ctor throws w/ multiple threads if leapfrog unsupported" );

  // Attempt to instantiate MKLRNG with multiple threads and a basic RNG that
  // does not support leapfrogging, i.e., multiple streams
  try {
    tk::MKLRNG r( 2, VSL_BRNG_MRG32K3A );
    fail( "should throw exception in both DEBUG and RELEASE mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "MKL VSL Error Code" ) !=
              std::string::npos );
  }
}

//! Test Gaussian generator statistics from mcg31 using multiple threads
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 4 >() {
  set_test_name( "Gaussian mcg31 from 4 emulated streams" );

  tk::MKLRNG r( 4, VSL_BRNG_MCG31 );
  RNG_common::test_gaussian( r );
}

//! Test Gaussian generator statistics from mcg59 using multiple threads
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 5 >() {
  set_test_name( "Gaussian mcg59 from 4 emulated streams" );

  tk::MKLRNG r( 4, VSL_BRNG_MCG59 );
  RNG_common::test_gaussian( r );
}

// //! Test Gaussian generator statistics from wh using multiple threads
// //! \note For some reason this generator consistently fails to generate
// //!   reasonable Gaussian random numbers
// //! \author J. Bakosi
// template<> template<>
// void MKLRNG_object::test< 6 >() {
//   set_test_name( "Gaussian wh from 4 emulated streams" );
//
//   tk::MKLRNG r( 4, VSL_BRNG_WH );
//   RNG_common::test_gaussian( r );
// }

//! Test beta generator statistics from mcg31 using multiple threads
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 7 >() {
  set_test_name( "beta mcg31 from 4 emulated streams" );

  tk::MKLRNG r( 4, VSL_BRNG_MCG31 );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from mcg59 using multiple threads
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 8 >() {
  set_test_name( "beta mcg59 from 4 emulated streams" );

  tk::MKLRNG r( 4, VSL_BRNG_MCG59 );
  RNG_common::test_beta( r );
}

//! Test beta generator statistics from wh using multiple threads
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 9 >() {
  set_test_name( "beta wh from 4 emulated streams" );

  tk::MKLRNG r( 4, VSL_BRNG_WH );
  RNG_common::test_beta( r );
}

//! Test uniform generator statistics from mcg31 using multiple threads
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 10 >() {
  set_test_name( "uniform mcg31 from 4 emulated streams" );

  tk::MKLRNG r( 4, VSL_BRNG_MCG31 );
  RNG_common::test_uniform( r );
}

//! Test uniform generator statistics from mcg59 using multiple threads
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 11 >() {
  set_test_name( "uniform mcg59 from 4 emulated streams" );

  tk::MKLRNG r( 4, VSL_BRNG_MCG59 );
  RNG_common::test_uniform( r );
}

// //! Test uniform generator statistics from wh using multiple threads
// //! \note For some reason this generator consistently fails to generate
// //!   reasonable uniform random numbers: it generates samples outside the
// //!   bounds (-eps, 1.0+eps)
// //! \author J. Bakosi
// template<> template<>
// void MKLRNG_object::test< 12 >() {
//   set_test_name( "uniform wh from 4 emulated streams" );
//
//   tk::MKLRNG r( 4, VSL_BRNG_WH );
//   RNG_common::test_uniform( r );
// }

//! Test copy constructor for mcg31
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 13 >() {
  set_test_name( "copy constructor with mcg31" );

  tk::MKLRNG p( 1 );                    // one thread
  RNG_common::test_copy_ctor( p );
  tk::MKLRNG r( 4, VSL_BRNG_MCG31 );    // 4 threads
  RNG_common::test_copy_ctor( r );
}

//! Test copy constructor for mcg59
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 14 >() {
  set_test_name( "copy constructor with mcg59" );

  tk::MKLRNG r( 4, VSL_BRNG_MCG59 );
  RNG_common::test_copy_ctor( r );
}

// //! Test copy constructor for wh
// //! \note For some reason this generator consistently fails to generate
// //!   reasonable Gaussian random numbers
// //! \author J. Bakosi
// template<> template<>
// void MKLRNG_object::test< 15 >() {
//   set_test_name( "copy constructor with wh" );
//
//   tk::MKLRNG r( 4, VSL_BRNG_WH );
//   RNG_common::test_copy_ctor( r );
// }

//! Test move constructor for mcg31
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 16 >() {
  set_test_name( "move constructor with mcg31" );

  tk::MKLRNG r( 4, VSL_BRNG_MCG31 );
  RNG_common::test_move_ctor( r );
}

//! Test move constructor for mcg59
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 17 >() {
  set_test_name( "move constructor with mcg59" );

  tk::MKLRNG r( 4, VSL_BRNG_MCG59 );
  RNG_common::test_move_ctor( r );
}

// //! Test move constructor for wh
// //! \note For some reason this generator consistently fails to generate
// //!   reasonable Gaussian random numbers
// //! \author J. Bakosi
// template<> template<>
// void MKLRNG_object::test< 18 >() {
//   set_test_name( "move constructor with wh" );
//
//   tk::MKLRNG r( 4, VSL_BRNG_WH );
//   RNG_common::test_move_ctor( r );
// }

//! Test copy assignment for mcg31
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 19 >() {
  set_test_name( "copy assignment with mcg31" );

  tk::MKLRNG p( 1 );                    // one thread
  RNG_common::test_copy_assignment( p );
  tk::MKLRNG r( 4, VSL_BRNG_MCG31 );    // 4 threads
  RNG_common::test_copy_assignment( r );
}

//! Test copy assignment for mcg59
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 20 >() {
  set_test_name( "copy assignment with mcg59" );

  tk::MKLRNG r( 4, VSL_BRNG_MCG59 );
  RNG_common::test_copy_assignment( r );
}

// //! Test copy assignment for wh
// //! \note For some reason this generator consistently fails to generate
// //!   reasonable Gaussian random numbers
// //! \author J. Bakosi
// template<> template<>
// void MKLRNG_object::test< 21 >() {
//   set_test_name( "copy assignment with wh" );
//
//   tk::MKLRNG r( 4, VSL_BRNG_WH );
//   RNG_common::test_copy_assignment( r );
// }

//! Test move assignment for mcg31
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 22 >() {
  set_test_name( "move assignment with mcg31" );

  tk::MKLRNG p( 1 );                    // one thread
  RNG_common::test_move_assignment( p );
  tk::MKLRNG r( 4, VSL_BRNG_MCG31 );    // 4 threads
  RNG_common::test_move_assignment( r );
}

//! Test move assignment for mcg59
//! \author J. Bakosi
template<> template<>
void MKLRNG_object::test< 23 >() {
  set_test_name( "move assignment with mcg59" );

  tk::MKLRNG r( 4, VSL_BRNG_MCG59 );
  RNG_common::test_move_assignment( r );
}

// //! Test move assignment for wh
// //! \note For some reason this generator consistently fails to generate
// //!   reasonable Gaussian random numbers
// //! \author J. Bakosi
// template<> template<>
// void MKLRNG_object::test< 24 >() {
//   set_test_name( "move assignment with wh" );
//
//   tk::MKLRNG r( 4, VSL_BRNG_WH );
//   RNG_common::test_move_assignment( r );
// }

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

} // tut::

#endif // test_MKLRNG_h
