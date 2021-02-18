// *****************************************************************************
/*!
  \file      tests/unit/Control/TestStringParser.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Control/StringParser
  \details   Unit tests for Control/StringParser
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "StringParser.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct StringParser_common {
  // tk::StringParser only has a protected constructor: designed to be used as a
  // base class
  struct parser : tk::StringParser {
    explicit parser( const std::string& f ) : StringParser( f ) {}
    explicit parser( int argc, char** argv ) : StringParser( argc, argv ) {}
    const std::string& string() const { return m_string; }
  };
};

//! Test group shortcuts
using StringParser_group =
  test_group< StringParser_common, MAX_TESTS_IN_GROUP >;
using StringParser_object = StringParser_group::object;

//! Define test group
static StringParser_group StringParser( "Control/StringParser" );

//! Test definitions for group

//! Test if constructor stores string passed in as std::string
template<> template<>
void StringParser_object::test< 1 >() {
  set_test_name( "ctor stores string passed in" );

  parser p( "blah" );
  ensure_equals( "string passed in doesn't match", p.string(), "blah" );
}

//! Test if constructor stores string passed in as argc, argv
template<> template<>
void StringParser_object::test< 2 >() {
  set_test_name( "ctor stores argc,argv passed in" );

  int argc = 3;
  char s1[] = "executable";
  char s2[] = "banana";
  char s3[] = "apple";
  char* argv[] = { s1, s2, s3 };
  parser p( argc, argv );
  // Constructor ignores argv[0]
  ensure_equals( "string passed in doesn't match", p.string(), "banana apple " );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
