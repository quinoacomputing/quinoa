// *****************************************************************************
/*!
  \file      tests/unit/Base/TestTaggedTuplePrint.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for tk::TaggedTuplePrint
  \details   Unit tests for tk::TaggedTuplePrint
*/
// *****************************************************************************

#include <sstream>

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "TaggedTuplePrint.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct TaggedTuplePrint_common {
  // Tags
  struct nam { static std::string name() { return "name"; } };
  struct age { static std::string name() { return "age"; } };
  struct email { static std::string name() {return "email"; } };
  struct tag1 { static std::string name() { return "tag1"; } };
  struct tag2 { static std::string name() { return "tag2"; } };
  struct tag3 { static std::string name() { return "tag3"; } };

  using MemberList = brigand::list<
    nam,  std::string,
    age,   int,
    email, std::string,
    tag1,  tk::TaggedTuple< brigand::list <
              tag2, std::string,
              tag3, std::string > > >;

  // Define a tagged tuple: odd template arguments are tags, even ones are types
  using record = tk::TaggedTuple< MemberList >;

  // Constructor
  TaggedTuplePrint_common() {
    tup.get< nam >() = "Bob";
    tup.get< age >() = 32;
    tup.get< email >() = "bob@google.com";
    auto& t1 = tup.get< tag1 >();
    t1.get< tag2 >() = "string2";
    t1.get< tag3 >() = "string3";
  }

  record tup;
};

//! Test group shortcuts
using TaggedTuplePrint_group =
  test_group< TaggedTuplePrint_common, MAX_TESTS_IN_GROUP >;
using TaggedTuplePrint_object = TaggedTuplePrint_group::object;

//! Define test group
static TaggedTuplePrint_group TaggedTuplePrint( "Base/TaggedTuplePrint" );

//! Test definitions for group

//! Test operator<< of TaggedTuple
template<> template<>
void TaggedTuplePrint_object::test< 1 >() {
  set_test_name( "operator<<" );

  std::stringstream s;
  s << tup;
  ensure_equals( "operator<<(TaggedTuple)", s.str(), "name: 'Bob' age: '32' "
    "email: 'bob@google.com' tag1: { tag2: 'string2' tag3: 'string3' } " );
}

//! Test print() with empty ignore list
template<> template<>
void TaggedTuplePrint_object::test< 2 >() {
  set_test_name( "print() with empty ignore list" );

  std::stringstream s;
  tk::print( s, tup );
  ensure_equals( "print()", s.str(), "name: 'Bob' age: '32' "
    "email: 'bob@google.com' tag1: { tag2: 'string2' tag3: 'string3' } " );
}

//! Test print() with non-empty ignore list
template<> template<>
void TaggedTuplePrint_object::test< 3 >() {
  set_test_name( "print() with non-empty ignore list" );

  std::stringstream s;
  tk::print< MemberList, brigand::set< email > >( s, tup );
  ensure_equals( "print<ignore>()", s.str(), "name: 'Bob' age: '32' "
    "tag1: { tag2: 'string2' tag3: 'string3' } " );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
