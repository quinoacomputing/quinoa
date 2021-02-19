// *****************************************************************************
/*!
  \file      tests/unit/Base/TestWriter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Base/Writer.hpp
  \details   Unit tests for Base/Writer.hpp
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "Writer.hpp"

namespace unittest {

extern std::string g_executable;

} // unittest::

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct Writer_common {
  struct W : public tk::Writer {
    explicit W( const std::string& filename ) : tk::Writer( filename ) {}
  };
};

//! Test group shortcuts
using Writer_group = test_group< Writer_common, MAX_TESTS_IN_GROUP >;
using Writer_object = Writer_group::object;

//! Define test group
static Writer_group Writer( "Base/Writer" );

//! Test definitions for group

//! Test if constructor can open a file
template<> template<>
void Writer_object::test< 1 >() {
  set_test_name( "ctor can open file" );

  // throws exception if something goes wrong, which Template Unit Test catches
  W w( "very_little_chance_that_a_file_with_this_name_exists" );
  // clean up
  std::remove( "very_little_chance_that_a_file_with_this_name_exists" );
}

//! Test if constructor does not throw if empty filename is given
template<> template<>
void Writer_object::test< 2 >() {
  set_test_name( "ctor does not throw if filename empty" );

  // throws exception if something goes wrong, which Template Unit Test catches
  W w( "" );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
