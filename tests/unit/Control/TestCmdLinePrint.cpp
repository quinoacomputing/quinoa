// *****************************************************************************
/*!
  \file      tests/unit/Control/TestCmdLinePrint.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for tk::CmdLinePrint
  \details   Unit tests for tk::CmdLinePrint
*/
// *****************************************************************************

#include <sstream>

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "CmdLinePrint.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct CmdLinePrint_common {
  // Tags
  struct nam { static std::string name() { return "name"; } };
  struct age { static std::string name() { return "age"; } };
  struct email { static std::string name() {return "email"; } };
  struct tag1 { static std::string name() { return "tag1"; } };
  struct tag2 { static std::string name() { return "tag2"; } };
  struct tag3 { static std::string name() { return "tag3"; } };

  // Define a CmdLine: must inherit from TaggedTuple and must define ignore
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
  CmdLinePrint_common() {
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
using CmdLinePrint_group =
  test_group< CmdLinePrint_common, MAX_TESTS_IN_GROUP >;
using CmdLinePrint_object = CmdLinePrint_group::object;

//! Define test group
static CmdLinePrint_group CmdLinePrint( "Control/CmdLinePrint" );

//! Test definitions for group

//! Test print() of CmdLine
template<> template<>
void CmdLinePrint_object::test< 1 >() {
  set_test_name( "print()" );

  std::stringstream s;
  tk::print( s, cmd );
  ensure_equals( "print()", s.str(),
R"(# vim: filetype=sh:
#
# Contents of a tagged tuple.
#
# A string in single quotes denotes the name/tag of a (nested)
# tagged tuple. The contents of tuples are enclosed within braces,
# indented, and aligned compared to the parent tuple.

'cmdline' {
  name       : Bob
  age        : 32
  'tag1' {
    tag2       : string2
    tag3       : string3 } })" );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
