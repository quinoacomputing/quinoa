// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestHas.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Base/Has.h
  \details   Unit tests for Base/Has.h
*/
// *****************************************************************************
#ifndef test_Has_h
#define test_Has_h

#include <unistd.h>

#include "NoWarning/tut.h"

#include "Has.h"

namespace tut {

//! All tests in group inherited from this base
struct Has_common {
  struct noProxy {};
  struct yesProxy { using Proxy = int; };
  struct noAlias {};
  struct yesAlias { using alias = char; };
  struct noCode {};
  struct yesCode { using code = char; };
  struct noExpectType {};
  struct yesExpectType { struct expect{ using type = char; }; };
  struct noExpectDescription {};
  struct yesExpectDescription { struct expect{ void description(){} }; };
  struct noExpectChoices {};
  struct yesExpectChoices { struct expect{ void choices(){} }; };
};

//! Test group shortcuts
using Has_group = test_group< Has_common, MAX_TESTS_IN_GROUP >;
using Has_object = Has_group::object;

//! Define test group
static Has_group Has( "Base/Has" );

//! Test definitions for group

//! Test if tk::HasTypedefProxy correctly detects the absence of typedef Proxy
//! \author J. Bakosi
template<> template<>
void Has_object::test< 1 >() {
  set_test_name( "HasTypedefProxy detects absence" );

  struct no { bool value = tk::HasTypedefProxy< noProxy >::value; };
  ensure_equals( "struct has no Proxy", no().value, false );
}

//! Test if tk::HasTypedefProxy correctly detects the presence of typedef Proxy
//! \author J. Bakosi
template<> template<>
void Has_object::test< 2 >() {
  set_test_name( "HasTypedefProxy detects presence" );

  struct yes { bool value = tk::HasTypedefProxy< yesProxy >::value; };
  ensure_equals( "struct has Proxy", yes().value, true );
}

//! Test if tk::HasTypedefAlias correctly detects the absence of typedef alias
//! \author J. Bakosi
template<> template<>
void Has_object::test< 3 >() {
  set_test_name( "HasTypedefAlias detects absence" );

  struct no { bool value = tk::HasTypedefAlias< noAlias >::value; };
  ensure_equals( "struct has no alias", no().value, false );
}

//! Test if tk::HasTypedefAlias correctly detects the presence of typedef alias
//! \author J. Bakosi
template<> template<>
void Has_object::test< 4 >() {
  set_test_name( "HasTypedefAlias detects presence" );

  struct yes { bool value = tk::HasTypedefAlias< yesAlias >::value; };
  ensure_equals( "struct has alias", yes().value, true );
}

//! Test if tk::HasTypedefCode correctly detects the absence of typedef code
//! \author J. Bakosi
template<> template<>
void Has_object::test< 5 >() {
  set_test_name( "HasTypedefCode detects absence" );

  struct no { bool value = tk::HasTypedefCode< noCode >::value; };
  ensure_equals( "struct has no code", no().value, false );
}

//! Test if tk::HasTypedefCode correctly detects the presence of typedef code
//! \author J. Bakosi
template<> template<>
void Has_object::test< 6 >() {
  set_test_name( "HasTypedefCode detects presence" );

  struct yes { bool value = tk::HasTypedefCode< yesCode >::value; };
  ensure_equals( "struct has code", yes().value, true );
}

//! \brief Test if tk::HasTypedefExpectType correctly detects the absence of
//!   typedef expect::type
//! \author J. Bakosi
template<> template<>
void Has_object::test< 7 >() {
  set_test_name( "HasTypedefExpectType: absence" );

  struct no { bool value = tk::HasTypedefExpectType< noExpectType >::value; };
  ensure_equals( "struct has no expect::type", no().value, false );
}

//! \brief Test if tk::HasTypedefExpectType correctly detects the presence of
//!   typedef expect::type
//! \author J. Bakosi
template<> template<>
void Has_object::test< 8 >() {
  set_test_name( "HasTypedefExpectType: presence" );

  struct yes { bool value = tk::HasTypedefExpectType< yesExpectType >::value; };
  ensure_equals( "struct has expect::type", yes().value, true );
}

//! \brief Test if tk::HasFunctionExpectDescription correctly detects the
//!   absence of function expect::description()
//! \author J. Bakosi
template<> template<>
void Has_object::test< 9 >() {
  set_test_name( "HasFunctionExpectDescription: absence" );

  struct no {
    bool value = tk::HasFunctionExpectDescription< noExpectDescription >::value;
  };
  ensure_equals( "struct has no expect::description", no().value, false );
}

//! \brief Test if tk::HasFunctionExpectDescription correctly detects the
//!   presence of function expect::description()
//! \author J. Bakosi
template<> template<>
void Has_object::test< 10 >() {
  set_test_name( "HasFunctionExpectDescription: presence" );

  struct yes {
    bool value = tk::HasFunctionExpectDescription< yesExpectDescription >::value;
  };
  ensure_equals( "struct has expect::description", yes().value, true );
}

//! \brief Test if tk::HasFunctionExpectChoices correctly detects the absence of
//!   function expect::choices()
//! \author J. Bakosi
template<> template<>
void Has_object::test< 11 >() {
  set_test_name( "HasFunctionExpectChoices: absence" );

  struct no {
    bool value = tk::HasFunctionExpectChoices< noExpectChoices >::value;
  };
  ensure_equals( "struct has no expect::choices", no().value, false );
}

//! \brief Test if tk::HasFunctionExpectChoices correctly detects the presence
//!   of function expect::choices()
//! \author J. Bakosi
template<> template<>
void Has_object::test< 12 >() {
  set_test_name( "HasFunctionExpectChoices: presence" );

  struct yes {
    bool value = tk::HasFunctionExpectChoices< yesExpectChoices >::value;
  };
  ensure_equals( "struct has expect::choices", yes().value, true );
}

} // tut::

#endif // test_Has_h
