//******************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestWriter.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 03:11:58 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unit tests for Base/Writer
  \details   Unit tests for Base/Writer
*/
//******************************************************************************
#ifndef test_Writer_h
#define test_Writer_h

#include <tut/tut.hpp>

#include "Writer.h"

namespace unittest {

extern std::string g_executable;

} // unittest::

namespace tut {

//! All tests in group inherited from this base
struct Writer_common {
  struct W : public tk::Writer {
     W( const std::string& filename ) : tk::Writer( filename ) {}
  };
};

//! Test group shortcuts
using Writer_group = test_group< Writer_common, MAX_TESTS_IN_GROUP >;
using Writer_object = Writer_group::object;

//! Define test group
Writer_group Writer( "Base/Writer" );

//! Test definitions for group

//! Test if constructor can open a file
//! \author J. Bakosi
template<> template<>
void Writer_object::test< 1 >() {
  set_test_name( "ctor can open file" );

  // throws exception if something goes wrong, which Template Unit Test catches
  W w( "very_little_chance_that_a_file_with_this_name_exists" );
  // clean up
  std::remove( "very_little_chance_that_a_file_with_this_name_exists" );
}

//! Test if constructor does not throw if empty filename is given
//! \author J. Bakosi
template<> template<>
void Writer_object::test< 2 >() {
  set_test_name( "ctor does not throw if filename empty" );

  // throws exception if something goes wrong, which Template Unit Test catches
  W w( "" );
}

} // tut::

#endif // test_Writer_h
