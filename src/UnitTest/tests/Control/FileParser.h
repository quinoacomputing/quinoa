//******************************************************************************
/*!
  \file      src/UnitTest/tests/Control/FileParser.h
  \author    J. Bakosi
  \date      Tue 05 Aug 2014 04:02:51 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Unit tests for Control/FileParser
  \details   Unit tests for Control/FileParser
*/
//******************************************************************************
#ifndef test_FileParser_h
#define test_FileParser_h

#include <tut/tut.hpp>
#include <FileParser.h>

namespace unittest {

extern std::string g_executable;

} // unittest::

namespace tut {

//! All tests in group inherited from this base
struct FileParser_common {
  // tk::FileParser only has a protected constructor: designed to be used as a
  // base class
  struct parser : tk::FileParser {
    parser( const std::string& f ) : FileParser( f ) {}
  };
};

//! Test group shortcuts
using FileParser_group = test_group< FileParser_common >;
using FileParser_object = FileParser_group::object;

//! Define test group
FileParser_group FileParser( "Control/FileParser" );

//! Test definitions for group

//! Test if constructor finds and can open an existing file (the executable)
template<> template<>
void FileParser_object::test< 1 >() {
  set_test_name( "ctor finds and opens executable" );

  // Throws exception if something goes wrong, which Template Unit Test catches
  parser p( unittest::g_executable );
}

//! Test if constructor throws an exception if empty filename is given
template<> template<>
void FileParser_object::test< 2 >() {
  set_test_name( "ctor throws if filename empty" );

  // Correctly throws exception in DEBUG mode if empty filename string is given
  try {

    parser p( "" );

  } catch( tk::Exception& e ) {
    // exception thrown, test ok
    // if any other type of exception is thrown, test fails with except
    // find  out if exception was thrown due to the correct reason: testing on
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
template<> template<>
void FileParser_object::test< 3 >() {
  set_test_name( "ctor throws if file doesn't exist" );

  // Correctly throws exception if file does not exist
  try {

    parser p( "very_little_chance_that_a_file_with_this_name_exists" );

  } catch( tk::Exception& e ) {
    // exception thrown, test ok
    // if any other type of exception is thrown, test fails with except
    // find  out if exception was thrown due to the correct reason
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Failed to open file" ) !=
              std::string::npos );
  }
}

//! Test if constructor throws exception if it cannot read from file
template<> template<>
void FileParser_object::test< 4 >() {
  set_test_name( "ctor throws if it cannot read from file" );

  // Correctly throws exception if cannot read from file
  try {

    parser p( "." );    // Nasty: trying open an existing directory as a file

  } catch( tk::Exception& e ) {
    // exception thrown, test ok
    // if any other type of exception is thrown, test fails with except
    // find  out if exception was thrown due to the correct reason
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Failed to read from file" ) !=
              std::string::npos );
  }
}

} // tut::

#endif // test_FileParser_h
