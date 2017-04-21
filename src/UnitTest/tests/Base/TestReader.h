// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestReader.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Base/Reader
  \details   Unit tests for Base/Reader
*/
// *****************************************************************************
#ifndef test_Reader_h
#define test_Reader_h

#include "NoWarning/tut.h"

#include "Reader.h"

namespace unittest {

extern std::string g_executable;

} // unittest::

namespace tut {

//! All tests in group inherited from this base
struct Reader_common {};

//! Test group shortcuts
using Reader_group = test_group< Reader_common, MAX_TESTS_IN_GROUP >;
using Reader_object = Reader_group::object;

//! Define test group
static Reader_group Reader( "Base/Reader" );

//! Test definitions for group

//! Test if constructor finds and can open an existing file (the executable)
//! \author J. Bakosi
template<> template<>
void Reader_object::test< 1 >() {
  set_test_name( "ctor finds and opens executable" );

  // Throws exception if something goes wrong, which Template Unit Test catches
  tk::Reader r( unittest::g_executable );
}

//! Test if constructor throws an exception if empty filename is given
//! \author J. Bakosi
template<> template<>
void Reader_object::test< 2 >() {
  set_test_name( "ctor throws if filename empty" );

  // Correctly throws exception in DEBUG mode if empty filename string is given
  try {

    tk::Reader r( "" );
    fail( "should throw exception" );

  } catch( tk::Exception& e ) {
    // exception thrown, test ok
    // if any other type of exception is thrown, test fails with except
    // find out if exception was thrown due to the correct reason: testing on
    // whether the filename is empty is an Assert and compiled away in RELEASE
    // mode, in which case the ErrChk throws when it fails to open the file
    #ifdef NDEBUG
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Failed to open file" ) !=
              std::string::npos );
    #else
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "No filename specified" ) !=
              std::string::npos );
    #endif
  }
}

//! Test if constructor throws exception if file does not exist
//! \author J. Bakosi
template<> template<>
void Reader_object::test< 3 >() {
  set_test_name( "ctor throws if file doesn't exist" );

  // Correctly throws exception if file does not exist
  try {

    tk::Reader p( "very_little_chance_that_a_file_with_this_name_exists" );
    fail( "should throw exception" );

  } catch( tk::Exception& e ) {
    // exception thrown, test ok
    // if any other type of exception is thrown, test fails with except
    // find  out if exception was thrown due to the correct reason
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Failed to open file" ) !=
              std::string::npos ||
            std::string( e.what() ).find( "Failed to read from file" ) !=
              std::string::npos );
  }
}

//! Test if constructor throws exception if cannot read from file
//! \author J. Bakosi
template<> template<>
void Reader_object::test< 4 >() {
  set_test_name( "ctor throws if cannot read from file" );

  // Correctly throws exception if cannot read from file
  try {

    tk::Reader p( "." );   // Nasty: trying open an existing directory as a file

  } catch( tk::Exception& e ) {
    // exception thrown, test ok
    // if any other type of exception is thrown, test fails with except
    // find  out if exception was thrown due to the correct reason
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Failed to read from file" ) !=
              std::string::npos );
  }
}

//! Test if function firstline() can read a line
//! \author J. Bakosi
template<> template<>
void Reader_object::test< 5 >() {
  set_test_name( "firstline() can read a line" );

  tk::Reader r( unittest::g_executable );
  r.firstline();        // if throws, TUT catches it, throw away return value
}

} // tut::

#endif // test_Reader_h
