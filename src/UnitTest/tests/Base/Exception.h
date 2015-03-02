//******************************************************************************
/*!
  \file      src/UnitTest/tests/Base/Exception.h
  \author    J. Bakosi
  \date      Sat 28 Feb 2015 10:17:08 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unit tests for Base/Exception.h
  \details   Unit tests for Base/Exception.h
*/
//******************************************************************************
#ifndef test_Exception_h
#define test_Exception_h

#include <tut/tut.hpp>
#include <Exception.h>

namespace tut {

//! All tests in group inherited from this base
struct Exception_common {};

//! Test group shortcuts
using Exception_group = test_group< Exception_common >;
using Exception_object = Exception_group::object;

//! Define test group
Exception_group Exception( "Base/Exception" );

//! Test definitions for group

//! Test constructor with line number info
template<> template<>
void Exception_object::test< 1 >() {
  set_test_name( "constructor message w/ line number info" );

  tk::Exception e( "msg", "file", "func", 12 );
  ensure_equals( "get exception message",
                 std::string( e.what() ),
                 std::string( "msg\n>>> Exception in file:12: func" ) );
}

//! Test constructor without line number info
template<> template<>
void Exception_object::test< 2 >() {
  set_test_name( "constructor message w/o line number info" );

  tk::Exception e( "msg", "file", "func", 0 );
  ensure_equals( "get exception message",
                 std::string( e.what() ),
                 std::string( "msg\n>>> No file:line:func information from "
                              "exception" ) );
}

//! Test Throw macro
template<> template<>
void Exception_object::test< 3 >() {
  set_test_name( "Throw macro" );

  try {
    Throw( "msg" );
    fail( "should throw exception" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown, test ok
  }
}

//! Test that Assert macro throws if condition is false
template<> template<>
void Exception_object::test< 4 >() {
  set_test_name( "Assert macro throws if condition is false" );

  try {
    Assert( 0 == 1, "msg" );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }
}

//! Test that Assert macro does not throw if condition is true
template<> template<>
void Exception_object::test< 5 >() {
  set_test_name( "Assert macro doesn't throw if cond is true" );

  try {
    Assert( 1 == 1, "msg" );
  }
  catch ( tk::Exception& e ) {
    fail( "should not throw exception" );
  }
}

//! Test that ErrChk macro throws if condition is false
template<> template<>
void Exception_object::test< 6 >() {
  set_test_name( "ErrChk macro throws if condition is false" );

  try {
    ErrChk( 0 == 1, "msg" );
    fail( "should throw excecption" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown, test ok
  }
}

//! Test that ErrChk macro does not throw if condition is true
template<> template<>
void Exception_object::test< 7 >() {
  set_test_name( "ErrChk macro doesn't throw if cond is true" );

  try {
    ErrChk( 0 != 1, "msg" );
  }
  catch ( tk::Exception& e ) {
    fail( "should not throw excecption" );
  }
}

} // tut::

#endif // test_Exception_h
