// *****************************************************************************
/*!
  \file      tests/unit/Base/TestException.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Base/Exception.hpp
  \details   Unit tests for Base/Exception.hpp
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "Exception.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct Exception_common {};

//! Test group shortcuts
using Exception_group = test_group< Exception_common, MAX_TESTS_IN_GROUP >;
using Exception_object = Exception_group::object;

//! Define test group
static Exception_group Exception( "Base/Exception" );

//! Test definitions for group

//! Test constructor with line number info
template<> template<>
void Exception_object::test< 1 >() {
  set_test_name( "constructor message w/ line number info" );

  tk::Exception e( "msg", "file", "func", 12 );
  ensure_equals( "get exception message",
                 std::string( e.what() ),
                 std::string( "msg\n>>> Exception at file:12: func" ) );
}

//! Test constructor without line number info
template<> template<>
void Exception_object::test< 2 >() {
  set_test_name( "constructor message w/o line number info" );

  tk::Exception e( "msg", "file", "func", 0 );
  ensure_equals( "get exception message",
                 std::string( e.what() ),
                 std::string( "msg" ) );
}

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunreachable-code"
  #pragma clang diagnostic ignored "-Wcomma"
#endif

//! Test Throw macro
template<> template<>
void Exception_object::test< 3 >() {
  set_test_name( "Throw macro" );

  try {
    // Use the comma operator to make cppcheck ignore unreachable code error
    Throw( "msg" ), fail( "should throw exception" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
}

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

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
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }
}

//! Test that Assert macro does not throw if condition is true
template<> template<>
void Exception_object::test< 5 >() {
  set_test_name( "Assert macro doesn't throw if cond is true" );

  try {
    // cppcheck-suppress duplicateExpression
    Assert( 1 == 1, "msg" );
  }
  catch ( tk::Exception& ) {
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
  catch ( tk::Exception& ) {
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
  catch ( tk::Exception& ) {
    fail( "should not throw excecption" );
  }
}

// //! Test tk::Exception::handleException()
// //! \note Disabled until we find a way to redirect/divert printf's output to
// //!   from stdout to a file (buffer, etc.) and be able to safely restore the
// //!   output to stdout in a way that correctly interoperates with Charm++.
// //!   Currently, the redirection to file works using the code below, but
// //!   restoration happens in an unpredictable way: sometimes works sometimes
// //!   does not, swallowing the rest of the output from the unit test harness.
// template<> template<>
// void Exception_object::test< 8 >() {
//   set_test_name( "handleException" );
//
//   try {
//     ErrChk( 0 == 1, "msg" );    // will throw tk::Exception
//   }
//   catch ( tk::Exception& e ) {
//
//     // redirect printf's output to file (requires POSIX)
//     // see http://stackoverflow.com/a/11110451
//     int stdout_fd = dup( STDOUT_FILENO );
//     freopen( "handle_output", "w", stdout );
//
//     e.handleException();
//
//     // restore stdout to its original state (requires POSIX)
//     fclose( stdout );
//     dup2( stdout_fd, STDOUT_FILENO );
//     stdout = fdopen( STDOUT_FILENO, "w" );
//     close( stdout_fd );
//     // remove handle_output from disk
//     tk::rm( "handle_output" );
//
//   }
//
// }

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
