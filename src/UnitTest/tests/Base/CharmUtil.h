//******************************************************************************
/*!
  \file      src/UnitTest/tests/Base/CharmUtil.h
  \author    J. Bakosi
  \date      Mon 28 Jul 2014 02:11:25 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Unit tests for Base/CharmUtil.h
  \details   Unit tests for Base/CharmUtil.h
*/
//******************************************************************************
#ifndef test_CharmUtil_h
#define test_CharmUtil_h

#include <unistd.h>
#include <tut/tut.hpp>
#include <CharmUtil.h>

namespace tut {

//! All tests in group inherited from this base
struct CharmUtil_common {
  struct noProxy {};
  struct yesProxy { using Proxy = int; };
};

//! Test group shortcuts
using CharmUtil_group = test_group< CharmUtil_common >;
using CharmUtil_object = CharmUtil_group::object;

//! Define test group
CharmUtil_group CharmUtil( "Base/CharmUtil" );

//! Test definitions for group

//! Test if tk::HasProxy correctly detects the presence of typedef Proxy
template<> template<>
void CharmUtil_object::test< 1 >() {
  set_test_name( "HasProxy detects absence of typedef Proxy" );

  struct no { bool value = tk::HasProxy< noProxy >::value; };
  ensure_equals( "struct has no Proxy", no().value, false );
}

//! Test if tk::HasProxy correctly detects the absence of typedef Proxy
template<> template<>
void CharmUtil_object::test< 2 >() {
  set_test_name( "HasProxy detects presence of typedef Proxy" );

  struct yes { bool value = tk::HasProxy< yesProxy >::value; };
  ensure_equals( "struct has Proxy", yes().value, true );
}

//! Test if is_enum_class correctly detects a strongly-typed enum
template<> template<>
void CharmUtil_object::test< 3 >() {
  set_test_name( "is_enum_class detects a strongly-typed enum" );

  enum class A { FIELD1, FIELD2 };
  struct yes { bool value = tk::is_enum_class< A >::value; };
  ensure_equals( "enum is strongly typed", yes().value, true );
}

//! Test if is_enum_class correctly detects a C-style enum
template<> template<>
void CharmUtil_object::test< 4 >() {
  set_test_name( "is_enum_class detects a C-style enum" );

  enum A { FIELD1, FIELD2 };
  struct no { bool value = tk::is_enum_class< A >::value; };
  ensure_equals( "enum is C-style", no().value, false );
}

} // tut::

#endif // test_CharmUtil_h
