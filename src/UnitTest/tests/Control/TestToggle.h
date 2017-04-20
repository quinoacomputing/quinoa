// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Control/TestToggle.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Control/Toggle
  \details   Unit tests for Control/Toggle
*/
// *****************************************************************************
#ifndef test_Toggle_h
#define test_Toggle_h

#include "NoWarning/tut.h"

#include "Toggle.h"

namespace tut {

//! All tests in group inherited from this base
struct Toggle_common {
  enum class Enum { F1, F2, F3, F4 };
  //! Toggle only has a protected constructor, designed to be base class only
  struct toggle : tk::Toggle< Enum > {
     // The child constructor must initialize all state data Toggle holds: a
     // std::string, a std::map< Enum, std::string >, and its inverse, a
     // std::map< std::string, Enum >. The easiest is to use initializer lists
     // and initialize in place. The size of the two maps must be equal.
     toggle( std::string&& g,
             std::map< Enum, std::string >&& n,
             std::map< std::string, Enum >&& v ) :
      tk::Toggle< Enum >( std::move(g), std::move(n), std::move(v) ) {}
  };
  //! Example toggle
  const toggle p = toggle( "switch name",
                           //! map Enums to names
                           { { Enum::F1, "1st option name" },
                             { Enum::F2, "2nd option name" },
                             { Enum::F3, "3rd option name" } },
                           //! map Enums to names
                           { { "1st option keyword", Enum::F1 },
                             { "2nd option keyword", Enum::F2 },
                             { "3rd option keyword", Enum::F3 } } );
};

//! Test group shortcuts
using Toggle_group = test_group< Toggle_common, MAX_TESTS_IN_GROUP >;
using Toggle_object = Toggle_group::object;

//! Define test group
static Toggle_group Toggle( "Control/Toggle" );

//! Test definitions for group

//! Test that constructor throws if maps aren't the same size
//! \author J. Bakosi
template<> template<>
void Toggle_object::test< 1 >() {
  set_test_name( "ctor throws if maps aren't the same size" );

  try {

    toggle t( "some switch name describing the option",
              //! map Enums to names
              { { Enum::F1, "1st option name" },
                { Enum::F2, "2nd option name" },
                { Enum::F3, "3rd option name" } },
              //! map Enums to names
              { { "1st option keyword", Enum::F1 },
                { "2nd option keyword", Enum::F2 } } );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif

  } catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "map sizes differ in Toggle" ) !=
              std::string::npos );
  }
}

//! Test that member function group() returns the correct groupname
//! \author J. Bakosi
template<> template<>
void Toggle_object::test< 2 >() {
  set_test_name( "group() returns groupname" );
  ensure_equals( "groupname incorrect", p.group(), "switch name" );
}

//! Test that member function value() finds value for keyword
//! \author J. Bakosi
template<> template<>
void Toggle_object::test< 3 >() {
  set_test_name( "value() finds value for keyword" );
  ensure( "cannot find keyword", p.value("2nd option keyword") == Enum::F2 );
}

//! Test that member function value() throws in DEBUG mode if can't find keyword
//! \author J. Bakosi
template<> template<>
void Toggle_object::test< 4 >() {
  set_test_name( "value() throws if can't find keyword" );

  try {

    p.value("bogus keyword");
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif

  } catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Cannot find value for keyword" ) !=
              std::string::npos );

  }
}

//! Test that member function name() finds name for value
//! \author J. Bakosi
template<> template<>
void Toggle_object::test< 5 >() {
  set_test_name( "name() finds name for value" );
  ensure( "cannot find value", p.name( Enum::F2 ) == "2nd option name" );
}

//! Test that member function name() throws in DEBUG mode if can't find value
//! \author J. Bakosi
template<> template<>
void Toggle_object::test< 6 >() {
  set_test_name( "name() throws if can't find value" );

  try {

    p.name( Enum::F4 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif

  } catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok, Assert skipped in RELEASE mode
    // if any other type of exception is thrown, test fails with except
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "Cannot find name for value" ) !=
              std::string::npos );

  }
}

//! Test that member function exist() finds existing keyword
//! \author J. Bakosi
template<> template<>
void Toggle_object::test< 7 >() {
  set_test_name( "exist() finds existing keyword" );
  ensure_equals( "cannot find existing keyword",
                 p.exist( "1st option keyword" ), true );
}

//! Test that member function exists() does not find non-existing keyword
//! \author J. Bakosi
template<> template<>
void Toggle_object::test< 8 >() {
  set_test_name( "exist() doesn't find non-existing keyword" );
  ensure_equals( "finds non-existing keyword",
                 p.exist( "bogus keyword" ), false );
}

//! Test copy constructor
//! \author J. Bakosi
template<> template<>
void Toggle_object::test< 9 >() {
  set_test_name( "copy constructor" );

  toggle c( p );
  std::vector< toggle > v;
  v.push_back( c );
  ensure_equals( "copy constructor used to push_back a child of Toggle to a "
                 "std::vector", v[0].exist( "2nd option keyword" ), true );
}

//! Test move constructor
//! \author J. Bakosi
template<> template<>
void Toggle_object::test< 10 >() {
  set_test_name( "move constructor" );

  toggle c( p );
  std::vector< toggle > v;
  v.emplace_back( std::move(c) );
  ensure_equals( "move constructor used to emplace_back a child of Toggle to a "
                 "std::vector", v[0].group(), "switch name" );
}

//! Test copy assignment
//! \author J. Bakosi
template<> template<>
void Toggle_object::test< 11 >() {
  set_test_name( "copy assignment" );

  // Only map first two fields of enum
  toggle c = toggle( "second switch name",
                     //! map Enums to names
                     { { Enum::F1, "1st option name" },
                       { Enum::F2, "2nd option name" } },
                     //! map Enums to names
                     { { "1st option keyword", Enum::F1 },
                       { "2nd option keyword", Enum::F2 } } );
  // Copy p to c
  c = p;
  // c now should have a different group name and should have 3 options mapped
  ensure_equals( "new group name of copy-assigned Toggle",
                 c.group(), "switch name" );
  ensure_equals( "find new key of copy-assigned Toggle",
                 c.exist( "3rd option keyword" ), true );
}

//! Test move assignment
//! \author J. Bakosi
template<> template<>
void Toggle_object::test< 12 >() {
  set_test_name( "move assignment" );

  // Only map first two fields of enum
  toggle c = toggle( "second switch name",
                     //! map Enums to names
                     { { Enum::F1, "1st option name" },
                       { Enum::F2, "2nd option name" } },
                     //! map Enums to names
                     { { "1st option keyword", Enum::F1 },
                       { "2nd option keyword", Enum::F2 } } );
  // Move p to c
  c = std::move( p );
  // c now should have a different group name and should have 3 options mapped
  ensure_equals( "new group name of move-assigned Toggle",
                 c.group(), "switch name" );
  ensure_equals( "find new key of moved-assigned Toggle",
                 c.exist( "3rd option keyword" ), true );
}

} // tut::

#endif // test_Toggle_h
