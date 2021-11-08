// *****************************************************************************
/*!
  \file      tests/unit/Base/TestHas.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Base/Has.hpp
  \details   Unit tests for Base/Has.hpp
*/
// *****************************************************************************

#include <unistd.h>

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "Has.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! All tests in group inherited from this base
struct Has_common {
  struct noAlias {};
  struct yesAlias { using alias = char; };
  struct noCode {};
  struct yesCode { using code = char; };
  struct noExpectDescription {};
  struct yesExpectDescription { struct expect{ void description(){} }; };
  struct noExpectLower {};
  struct yesExpectLower { struct expect{ int lower = 1; }; };
  struct noExpectUpper {};
  struct yesExpectUpper { struct expect{ int upper = 1; }; };
  struct noExpectChoices {};
  struct yesExpectChoices { struct expect{ void choices(){} }; };
};

//! Test group shortcuts
using Has_group = test_group< Has_common, MAX_TESTS_IN_GROUP >;
using Has_object = Has_group::object;

//! Define test group
static Has_group Has( "Base/Has" );

//! Test definitions for group

//! Test if tk::HasTypedef_alias correctly detects the absence of typedef alias
template<> template<>
void Has_object::test< 1 >() {
  set_test_name( "HasTypedef_alias detects absence" );

  struct no { bool value = tk::HasTypedef_alias_v< noAlias >; };
  ensure_equals( "struct has no alias", no().value, false );
}

//! Test if tk::HasTypedef_alias correctly detects the presence of typedef alias
template<> template<>
void Has_object::test< 2 >() {
  set_test_name( "HasTypedef_alias detects presence" );

  struct yes { bool value = tk::HasTypedef_alias_v< yesAlias >; };
  ensure_equals( "struct has alias", yes().value, true );
}

//! \brief Test if tk::HasFunction_expect_description correctly detects the
//!   absence of function expect::description()
template<> template<>
void Has_object::test< 3 >() {
  set_test_name( "HasFunction_expect_description: absence" );

  struct no {
    bool value = tk::HasFunction_expect_description_v< noExpectDescription >;
  };
  ensure_equals( "struct has no expect::description", no().value, false );
}

//! \brief Test if tk::HasFunction_expect_description correctly detects the
//!   presence of function expect::description()
template<> template<>
void Has_object::test< 4 >() {
  set_test_name( "HasFunction_expect_description: presence" );

  struct yes {
    bool value = tk::HasFunction_expect_description_v< yesExpectDescription >;
  };
  ensure_equals( "struct has expect::description", yes().value, true );
}

//! \brief Test if tk::HasVar_expect_lower correctly detects the
//!   absence of function expect::lower()
template<> template<>
void Has_object::test< 5 >() {
  set_test_name( "HasVar_expect_lower: absence" );

  struct no {
    bool value = tk::HasVar_expect_lower_v< noExpectLower >;
  };
  ensure_equals( "struct has no expect::lower", no().value, false );
}

//! \brief Test if tk::HasVar_expect_lower correctly detects the
//!   presence of variable expect::lower()
template<> template<>
void Has_object::test< 6 >() {
  set_test_name( "HasVar_expect_lower: presence" );

  struct yes {
    bool value = tk::HasVar_expect_lower_v< yesExpectLower >;
  };
  ensure_equals( "struct has expect::lower", yes().value, true );
}

//! \brief Test if tk::HasVar_expect_upper correctly detects the
//!   absence of function expect::upper()
template<> template<>
void Has_object::test< 7 >() {
  set_test_name( "HasVar_expect_upper: absence" );

  struct no {
    bool value = tk::HasVar_expect_upper_v< noExpectUpper >;
  };
  ensure_equals( "struct has no expect::upper", no().value, false );
}

//! \brief Test if tk::HasVar_expect_upper correctly detects the
//!   presence of variable expect::upper()
template<> template<>
void Has_object::test< 8 >() {
  set_test_name( "HasVar_expect_upper: presence" );

  struct yes {
    bool value = tk::HasVar_expect_upper_v< yesExpectUpper >;
  };
  ensure_equals( "struct has expect::upper", yes().value, true );
}

//! \brief Test if tk::HasFunction_expect_choices correctly detects the absence
//!   of function expect::choices()
template<> template<>
void Has_object::test< 9 >() {
  set_test_name( "HasFunction_expect_choices: absence" );

  struct no {
    bool value = tk::HasFunction_expect_choices_v< noExpectChoices >;
  };
  ensure_equals( "struct has no expect::choices", no().value, false );
}

//! \brief Test if tk::HasFunction_expect_choices correctly detects the presence
//!   of function expect::choices()
template<> template<>
void Has_object::test< 10 >() {
  set_test_name( "HasFunction_expect_choices: presence" );

  struct yes {
    bool value = tk::HasFunction_expect_choices_v< yesExpectChoices >;
  };
  ensure_equals( "struct has expect::choices", yes().value, true );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
