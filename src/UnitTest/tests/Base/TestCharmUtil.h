//******************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestCharmUtil.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 03:12:09 PM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Unit tests for Base/CharmUtil.h
  \details   Unit tests for Base/CharmUtil.h
*/
//******************************************************************************
#ifndef test_CharmUtil_h
#define test_CharmUtil_h

#include <unistd.h>

#include <tut/tut.hpp>

#include "CharmUtil.h"

namespace tut {

//! All tests in group inherited from this base
struct CharmUtil_common {
  struct noProxy {};
  struct yesProxy { using Proxy = int; };
};

//! Test group shortcuts
using CharmUtil_group = test_group< CharmUtil_common, MAX_TESTS_IN_GROUP >;
using CharmUtil_object = CharmUtil_group::object;

//! Define test group
CharmUtil_group CharmUtil( "Base/CharmUtil" );

//! Test definitions for group

//! Test if is_enum_class correctly detects a strongly-typed enum
//! \author J. Bakosi
template<> template<>
void CharmUtil_object::test< 1 >() {
  set_test_name( "is_enum_class detects a strongly-typed enum" );

  enum class A { FIELD1, FIELD2 };
  struct yes { bool value = tk::is_enum_class< A >::value; };
  ensure_equals( "enum is strongly typed", yes().value, true );
}

//! Test if is_enum_class correctly detects a C-style enum
//! \author J. Bakosi
template<> template<>
void CharmUtil_object::test< 2 >() {
  set_test_name( "is_enum_class detects a C-style enum" );

  enum A { FIELD1, FIELD2 };
  struct no { bool value = tk::is_enum_class< A >::value; };
  ensure_equals( "enum is C-style", no().value, false );
}

} // tut::

#endif // test_CharmUtil_h
