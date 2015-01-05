//******************************************************************************
/*!
  \file      src/UnitTest/tests/Base/Has.h
  \author    J. Bakosi
  \date      Sat 17 Jan 2015 07:05:25 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Unit tests for Base/Has.h
  \details   Unit tests for Base/Has.h
*/
//******************************************************************************
#ifndef test_Has_h
#define test_Has_h

#include <unistd.h>
#include <tut/tut.hpp>
#include <Has.h>

namespace tut {

//! All tests in group inherited from this base
struct Has_common {
  struct noProxy {};
  struct yesProxy { using Proxy = int; };
};

//! Test group shortcuts
using Has_group = test_group< Has_common >;
using Has_object = Has_group::object;

//! Define test group
Has_group Has( "Base/Has" );

//! Test definitions for group

//! Test if tk::HasTypedefProxy correctly detects the presence of typedef Proxy
template<> template<>
void Has_object::test< 1 >() {
  set_test_name( "HasTypedefProxy detects absence of typedef Proxy" );

  struct no { bool value = tk::HasTypedefProxy< noProxy >::value; };
  ensure_equals( "struct has no Proxy", no().value, false );
}

//! Test if tk::HasTypedefProxy correctly detects the absence of typedef Proxy
template<> template<>
void Has_object::test< 2 >() {
  set_test_name( "HasTypedefProxy detects presence of typedef Proxy" );

  struct yes { bool value = tk::HasTypedefProxy< yesProxy >::value; };
  ensure_equals( "struct has Proxy", yes().value, true );
}

} // tut::

#endif // test_Has_h
