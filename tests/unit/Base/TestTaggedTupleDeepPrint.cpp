// *****************************************************************************
/*!
  \file      tests/unit/Base/TestTaggedTupleDeepPrint.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for tk::TaggedTupleDeepPrint
  \details   Unit tests for tk::TaggedTupleDeepPrint
*/
// *****************************************************************************

#include <sstream>

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "TaggedTupleDeepPrint.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct TaggedTupleDeepPrint_common {
  // Tags
  struct nam { static std::string name() { return "name"; } };
  struct age { static std::string name() { return "age"; } };
  struct email { static std::string name() {return "email"; } };
  struct tag1 { static std::string name() { return "tag1"; } };
  struct tag2 { static std::string name() { return "tag2"; } };
  struct tag3 { static std::string name() { return "tag3"; } };

  // Define a tk::TaggedTuple by inheriting from TaggedTuple (optionally)
  // defining type ignore
  struct Cmd : public tk::TaggedTuple< brigand::list<
                        nam,  std::string,
                        age,   int,
                        email, std::string,
                        tag1,  tk::TaggedTuple< brigand::list <
                                  tag2, std::string,
                                  tag3, std::string > > > >
  {
    // ignore these tags when printing to a stream
    using ignore = brigand::set< email >;
  };

  // Constructor
  TaggedTupleDeepPrint_common() {
    cmd.get< nam >() = "Bob";
    cmd.get< age >() = 32;
    cmd.get< email >() = "bob@google.com";
    auto& t1 = cmd.get< tag1 >();
    t1.get< tag2 >() = "string2";
    t1.get< tag3 >() = "string3";
  }

  Cmd cmd;
};

//! Test group shortcuts
using TaggedTupleDeepPrint_group =
  test_group< TaggedTupleDeepPrint_common, MAX_TESTS_IN_GROUP >;
using TaggedTupleDeepPrint_object = TaggedTupleDeepPrint_group::object;

//! Define test group
static TaggedTupleDeepPrint_group
  TaggedTupleDeepPrint( "Base/TaggedTupleDeepPrint" );

//! Test definitions for group

//! Test print() of TaggedTuple with depth/indentation
template<> template<>
void TaggedTupleDeepPrint_object::test< 1 >() {
  set_test_name( "print()" );

  std::stringstream s;
  tk::print( s, "cmd", cmd );
  ensure_equals( "print()", s.str(),
R"(# vim: filetype=sh:
#
# Contents of a tagged tuple.
#
# A string in single quotes denotes the name/tag of a (nested)
# tagged tuple. The contents of tuples are enclosed within braces.
# Vectors are enclosed within square brackets. Keys of associative
# containers are in paretheses.

'cmd' {
  name            : Bob
  age             : 32
  'tag1' {
    tag2            : string2
    tag3            : string3
  }
})" );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
