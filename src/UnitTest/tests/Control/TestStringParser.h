// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Control/TestStringParser.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Control/StringParser
  \details   Unit tests for Control/StringParser
*/
// *****************************************************************************
#ifndef test_StringParser_h
#define test_StringParser_h

#include "NoWarning/tut.h"

#include "StringParser.h"

namespace tut {

//! All tests in group inherited from this base
struct StringParser_common {
  // tk::StringParser only has a protected constructor: designed to be used as a
  // base class
  struct parser : tk::StringParser {
    parser( const std::string& f ) : StringParser( f ) {}
    parser( int argc, char** argv ) : StringParser( argc, argv ) {}
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
//! \author J. Bakosi
template<> template<>
void StringParser_object::test< 1 >() {
  set_test_name( "ctor stores string passed in" );

  parser p( "blah" );
  ensure_equals( "string passed in doesn't match", p.string(), "blah" );
}

//! Test if constructor stores string passed in as argc, argv
//! \author J. Bakosi
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

#endif // test_StringParser_h
