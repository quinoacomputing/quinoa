// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestTaggedTuple.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Base/TaggedTuple.h
  \details   Unit tests for Base/TaggedTuple.h
*/
// *****************************************************************************
#ifndef test_TaggedTuple_h
#define test_TaggedTuple_h

#include "NoWarning/tut.h"

#include "TaggedTuple.h"

namespace tut {

//! All tests in group inherited from this base
struct TaggedTuple_common {
  // Tags
  struct name {};
  struct age {};
  struct email {};

  // Define a tagged tuple: odd template arguments are tags, even ones are types
  using record = tk::tuple::tagged_tuple< name,  std::string,
                                          age,   int,
                                          email, std::string >;

  // Declare a record, initialize with an std::initializer_list
  record tup{ "Bob", 32, "bob@bob.bob" };
};

//! Test group shortcuts
using TaggedTuple_group = test_group< TaggedTuple_common, MAX_TESTS_IN_GROUP >;
using TaggedTuple_object = TaggedTuple_group::object;

//! Define test group
static TaggedTuple_group TaggedTuple( "Base/TaggedTuple" );

//! Test definitions for group

//! Test const-ref accessor of tagged_tuple
//! \author J. Bakosi
template<> template<>
void TaggedTuple_object::test< 1 >() {
  set_test_name( "const-ref accessors" );

  // Defina a struct holding a record with const-ref accessors
  struct A {
    A( const record& r ) : m_rec( r ) {}
    const std::string& getName() const { return m_rec.get< name >(); }
    const int& getAge() const { return m_rec.get< age >(); }
    const std::string& getEmail() const { return m_rec.get< email >(); }
    const record m_rec;
  };

  // Use const-ref accessors to read out the initialized data
  A a( tup );
  ensure_equals( "const-ref accessor to std::string", a.getName(), "Bob" );
  ensure_equals( "const-ref accessor to int", a.getAge(), 32 );
  ensure_equals( "const-ref accessor to std::string",
                 a.getEmail(), "bob@bob.bob" );
}

//! Test non-const-ref accessors of tagged_tuple
//! \author J. Bakosi
template<> template<>
void TaggedTuple_object::test< 2 >() {
  set_test_name( "non-const-ref accessors" );

  // Use (non-const) ref accessors to read out the initialized data
  ensure_equals( "non-const-ref accessor to std::string",
                 tup.get< name >(), "Bob" );
  ensure_equals( "non-const-ref accessor to int",
                 tup.get< age >(), 32 );
  ensure_equals( "non-const-ref accessor to std::string",
                 tup.get< email >(), "bob@bob.bob" );
}

//! Test set(const char*) with const rvalue argument of tagged_tuple
//! \author J. Bakosi
template<> template<>
void TaggedTuple_object::test< 3 >() {
  set_test_name( "set(const char*) - const rvalue arg" );

  record t{ "Bob", 32, "bob@bob.bob" };
  t.set< name >( "Alice" );     // const rvalue
  ensure_equals( "get() after set(const char*)", t.get< name >(), "Alice" );
}

//! Test set(std::string) with non-const lvalue argument of tagged_tuple
//! \author J. Bakosi
template<> template<>
void TaggedTuple_object::test< 4 >() {
  set_test_name( "set(std::string) - non-const lvalue arg" );

  record t{ "Bob", 32, "bob@bob.bob" };
  std::string n( "Alice" );
  t.set< name >( n );           // non-const lvalue
  ensure_equals( "get() after set(std::string)", t.get< name >(), "Alice" );
  ensure_equals( "original source kept intact", n, "Alice" );
}

//! Test set(std::string) with const lvalue argument of tagged_tuple
//! \author J. Bakosi
template<> template<>
void TaggedTuple_object::test< 5 >() {
  set_test_name( "set(const std::string) - const lvalue" );

  record t{ "Bob", 32, "bob@bob.bob" };
  const std::string n( "Alice" );
  t.set< name >( n );           // const lvalue
  ensure_equals( "get() after set(std::string)", t.get< name >(), "Alice" );
  ensure_equals( "original source kept intact", n, "Alice" );
}

//! Test set(std::string) with rvalue argument of tagged_tuple
//! \author J. Bakosi
template<> template<>
void TaggedTuple_object::test< 6 >() {
  set_test_name( "set(std::string&&) - rvalue ref arg" );

  record t{ "Bob", 32, "bob@bob.bob" };
  std::string n( "Alice" );
  t.set< name >( std::move(n) );           // rvalue reference
  ensure_equals( "get() after set(std::string)", t.get< name >(), "Alice" );
  // n here should still be in a valid but unspecified state, so we try to use
  // it, if that's not the case, the constructor or the copy assignment throws
  std::string k( "abc" );
  k = n;
  ensure_equals( "original source in valid but unspecified state", k, n );
}

//! Test size of tagged_tuple
//! \author J. Bakosi
template<> template<>
void TaggedTuple_object::test< 7 >() {
  set_test_name( "tagged_tuple_size()" );
  ensure_equals( "tagged_tuple_size = std::tuple_size/2",
                 tk::tuple::tagged_tuple_size< record >::value, 3UL );
}

} // tut::

#endif // test_TaggedTuple_h
