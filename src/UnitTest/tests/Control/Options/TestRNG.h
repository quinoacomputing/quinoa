// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Control/Options/TestRNG.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Control/Options/RNG
  \details   Unit tests for Control/Options/RNG
*/
// *****************************************************************************
#ifndef test_RNGOptions_h
#define test_RNGOptions_h

#include "NoWarning/tut.h"

#include "Options/RNG.h"
#include "RNGParam.h"
#include "QuinoaConfig.h"

namespace tut {

//! All tests in group inherited from this base
struct RNGOptions_common {
  RNGOptions_common() : m() {}
  const tk::ctr::RNG m;
};

//! Test group shortcuts
using RNGOptions_group = test_group< RNGOptions_common, MAX_TESTS_IN_GROUP >;
using RNGOptions_object = RNGOptions_group::object;

//! Define test group
static RNGOptions_group RNGOptions( "Control/Options/RNGOptions" );

//! Test definitions for group

//! Test that member function param() finds RNG parameter for method type
//! \author J. Bakosi
template<> template<>
void RNGOptions_object::test< 1 >() {
  set_test_name( "param() finds RNG parameter" );
  ensure( "cannot find parameter",
          m.param( tk::ctr::RNGType::R123_THREEFRY ) == 0 );
}

//! Test that member function param() throws in DEBUG mode if can't find param
//! \author J. Bakosi
template<> template<>
void RNGOptions_object::test< 2 >() {
  set_test_name( "param() throws if can't find" );

  try {

    m.param( static_cast< tk::ctr::RNGType >( 250 ) );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif

  } catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Cannot find parameter" ) !=
              std::string::npos );
  }
}

//! Test copy constructor
//! \author J. Bakosi
template<> template<>
void RNGOptions_object::test< 3 >() {
  set_test_name( "copy constructor" );

  tk::ctr::RNG c( m );
  std::vector< tk::ctr::RNG > v;
  v.push_back( c );
  ensure( "copy constructor used to push_back a RNG object to a std::vector",
          v[0].param( tk::ctr::RNGType::R123_PHILOX ) == 1 );
}

//! Test move constructor
//! \author J. Bakosi
template<> template<>
void RNGOptions_object::test< 4 >() {
  set_test_name( "move constructor" );

  tk::ctr::RNG c( m );
  std::vector< tk::ctr::RNG > v;
  v.emplace_back( std::move(c) );
  ensure( "move constructor used to emplace_back a RNG object to a std::vector",
           v[0].param( tk::ctr::RNGType::R123_PHILOX ) == 1 );
}

//! Test copy assignment
//! \author J. Bakosi
template<> template<>
void RNGOptions_object::test< 5 >() {
  set_test_name( "copy assignment" );

  tk::ctr::RNG c;
  c = m;
  ensure( "find param of copy-assigned RNG",
          c.param( tk::ctr::RNGType::R123_PHILOX ) == 1 );
}

//! Test move assignment
//! \author J. Bakosi
template<> template<>
void RNGOptions_object::test< 6 >() {
  set_test_name( "move assignment" );

  tk::ctr::RNG c;
  c = std::move( m );
  ensure( "find param of move-assigned RNG",
          c.param( tk::ctr::RNGType::R123_PHILOX ) == 1 );
}

#ifdef HAS_MKL
//! Test that member function lib() finds MKL library type for a MKL rng
//! \author J. Bakosi
template<> template<>
void RNGOptions_object::test< 7 >() {
  set_test_name( "lib() finds MKL library type" );
  ensure( "cannot find library type",
          m.lib( tk::ctr::RNGType::MKL_R250 ) == tk::ctr::RNGLibType::MKL );
}
#endif

#ifdef HAS_RNGSSE2
//! Test that member function lib() finds RNGSSE library type for a RNGSSE rng
//! \author J. Bakosi
template<> template<>
void RNGOptions_object::test< 8 >() {
  set_test_name( "lib() finds RNGSSE library type" );
  ensure( "cannot find library type",
          m.lib( tk::ctr::RNGType::RNGSSE_GM29 ) ==
            tk::ctr::RNGLibType::RNGSSE );
}
#endif

//! \brief Test that member function lib() finds Random123 library type for a
//!   Random123 rng
//! \author J. Bakosi
template<> template<>
void RNGOptions_object::test< 9 >() {
  set_test_name( "lib() finds Random123 library type" );
  ensure( "cannot find library type",
          m.lib( tk::ctr::RNGType::R123_THREEFRY ) ==
            tk::ctr::RNGLibType::R123 );
}

#ifdef HAS_RNGSSE2
//! Test that member function supportSeq() returns true for an RNGSSE rng
//! \author J. Bakosi
template<> template<>
void RNGOptions_object::test< 10 >() {
  set_test_name( "supportsSeq() true for RNGSSE" );
  ensure_equals( "cannot find RNGSSE rng in support map",
                 m.supportsSeq( tk::ctr::RNGType::RNGSSE_GM29 ), true );
}
#endif

//! \brief Test that member function supportSeq() returns false for an
//!    non-RNGSSE rng
//! \author J. Bakosi
template<> template<>
void RNGOptions_object::test< 11 >() {
  set_test_name( "supportsSeq() false for non-RNGSSE" );
  ensure_equals( "cannot find non-RNGSSE rng in support map",
                 m.supportsSeq( tk::ctr::RNGType::R123_PHILOX ), false );
}

#ifdef HAS_RNGSSE2
//! \brief Test that member function param<>() returns default for non-specified
//!   parameter
//! \author J. Bakosi
template<> template<>
void RNGOptions_object::test< 12 >() {
  set_test_name( "param() correctly returns default" );

  // empty bundle: no parameter specified
  std::map< tk::ctr::RNGType, tk::ctr::RNGSSEParam > b;
  ensure_equals( "does not return default seed for no parameters",
                 m.param< tag::seed >
                        ( tk::ctr::RNGType::RNGSSE_GM31, 0U, b ), 0 );
}
#endif

#ifdef HAS_RNGSSE2
//! \brief Test that member function param<>() returns parameter for specified
//!   parameter
//! \author J. Bakosi
template<> template<>
void RNGOptions_object::test< 13 >() {
  set_test_name( "param() returns specified param" );

  // specify sequence length parameter for RNGSSE rng
  std::map< tk::ctr::RNGType, tk::ctr::RNGSSEParam > b {
    { tk::ctr::RNGType::RNGSSE_GQ581, { 12, tk::ctr::RNGSSESeqLenType::LONG } }
  };
  ensure( "does not return specified sequence length for RNGSSE rng",
          m.param< tag::seqlen >                    // query this field
                 ( tk::ctr::RNGType::RNGSSE_GQ581,      // query this rng
                   tk::ctr::RNGSSESeqLenType::SHORT,    // default if not spec'd
                   b ) ==                               // query this bundle
            tk::ctr::RNGSSESeqLenType::LONG );
}
#endif

} // tut::

#endif // test_RNGOptions_h
