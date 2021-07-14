// *****************************************************************************
/*!
  \file      tests/unit/Base/TestContainerUtil.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Base/ContainerUtil.hpp
  \details   Unit tests for Base/ContainerUtil.hpp
*/
// *****************************************************************************

#include <map>
#include <vector>
#include <string>
#include <unordered_map>

#include <unistd.h>

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "ContainerUtil.hpp"
#include "StatCtr.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif

//! All tests in group inherited from this base
struct ContainerUtil_common {
  // cppcheck-suppress unusedStructMember
  double precision = 1.0e-15;    // required floating-point precision
};

//! Test group shortcuts
using ContainerUtil_group =
  test_group< ContainerUtil_common, MAX_TESTS_IN_GROUP >;
using ContainerUtil_object = ContainerUtil_group::object;

//! Define test group
static ContainerUtil_group ContainerUtil( "Base/ContainerUtil" );

//! Test definitions for group

//! Test uniquecopy making the elements of a container unique (on copy)
template<> template<>
void ContainerUtil_object::test< 1 >() {
  set_test_name( "uniquecopy" );

  // std::vector
  std::vector< int > v{{ 1, 1, 2, 6, -2, 3, 5, 5 }},
                     correct{{ -2, 1, 2, 3, 5, 6 }};
  auto v1 = tk::uniquecopy( v );
  ensure( "make vector (c)unique incorrect", v1 == correct );

  // std::string
  std::string s( "blahblah" );
  auto s1 = tk::uniquecopy( s );
  ensure( "make string (c)unique incorrect", s1 == "abhl" );

  // std::vector< std::vector< tk::ctr::Term > >
  std::vector< std::vector< tk::ctr::Term > > p{ tk::ctr::mean('c',1),
                                                 tk::ctr::variance('y',2),
                                                 tk::ctr::variance('y',2),
                                                 tk::ctr::variance('y',2) };
  std::vector< std::vector< tk::ctr::Term > > c{ tk::ctr::mean('c',1),
                                                 tk::ctr::variance('y',2) };
  auto p1 = tk::uniquecopy( p );
  ensure( "make vector<vector<tk::ctr::Term>> (c)unique incorrect", p1 == c );
}

//! Test unique making the elements of a container unique
template<> template<>
void ContainerUtil_object::test< 2 >() {
  set_test_name( "unique" );

  // std::vector
  std::vector< int > v{{ 1, 1, 2, 6, -2, 3, 5, 5 }},
                     correct{{ -2, 1, 2, 3, 5, 6 }};
  tk::unique( v );
  ensure( "make vector unique incorrect", v == correct );

  // std::string
  std::string s( "blahblah" );
  tk::unique( s );
  ensure( "make string unique incorrect", s == "abhl" );

  // std::vector< std::vector< tk::ctr::Term > >
  std::vector< std::vector< tk::ctr::Term > > p{ tk::ctr::mean('c',1),
                                                 tk::ctr::variance('y',2),
                                                 tk::ctr::variance('y',2),
                                                 tk::ctr::variance('y',2) };
  std::vector< std::vector< tk::ctr::Term > > c{ tk::ctr::mean('c',1),
                                                 tk::ctr::variance('y',2) };
  tk::unique( p );
  ensure( "make vector<vector<tk::ctr::Term>> unique incorrect", c == p );
}

//! Test cref_find returning a const-ref to value found for key in map
template<> template<>
void ContainerUtil_object::test< 3 >() {
  set_test_name( "[c]ref_find" );

  std::map< int, std::string > m{ {1,"one"}, {2,"two"} };
  std::unordered_map< int, std::string > u{ {1,"one"}, {2,"two"} };

  const auto one = tk::cref_find(m,1);
  ensure_equals( "cref_find incorrect", one, "one" );
  auto two = tk::ref_find(u,2);
  ensure_equals( "ref_find incorrect", two, "two" );

  // test if cref_find throws in DEBUG when cannot find key
  // skipped in RELEASE mode, would yield segmentation fault
  #ifdef NDEBUG
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    auto three = tk::ref_find(u,3);
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test vector extents
template<> template<>
void ContainerUtil_object::test< 4 >() {
  set_test_name( "vector extents" );

  std::vector< int > v{{ 1, 1, 2, 6, -2, 3, 5, 5 }};
  auto e = tk::extents( v );
  ensure( "vector extents incorrect", e == std::array< int, 2 >{{ -2, 6 }} );

  std::vector< std::size_t > w{{ 1, 1, 2, 6, 2, 3, 5, 5 }};
  auto f = tk::extents( w );
  ensure( "vector extents incorrect",
          f == std::array< std::size_t, 2 >{{ 1, 6 }} );
}

//! Test associative container (value) extents
template<> template<>
void ContainerUtil_object::test< 5 >() {
  set_test_name( "[hash]map value extents" );

  std::map< int, tk::real > m{ {1,34.2}, {2,12.3}, {4,-3.0} };
  auto me = tk::extents( m );
  ensure_equals( "map min incorrect", me[0], -3.0, precision );
  ensure_equals( "map max incorrect", me[1], 34.2, precision );

  std::unordered_map< int, tk::real > u{ {1,34.2}, {2,12.3}, {4,-3.0} };
  auto mu = tk::extents( u );
  ensure_equals( "hash map min incorrect", mu[0], -3.0, precision );
  ensure_equals( "hash map max incorrect", mu[1], 34.2, precision );
}

//! Test operator += adding values of a std::vector to another one
template<> template<>
void ContainerUtil_object::test< 6 >() {
  set_test_name( "operator+= std::vector" );

  using tk::operator+=;

  // add non-empty vector to empty one: copy src to dst (main intended use-case)
  std::vector< tk::real > v1, v2{{ 4.0, 9.0, 2.0 }};
  v1 += v2;
  ensure_equals( "add non-empty vector to empty one, dst size incorrect",
                 v1.size(), 3 );
  ensure_equals( "add non-empty vector to empty one, src size incorrect",
                 v2.size(), 3 );
  ensure_equals( "add non-empty vector to empty one, src[0]==dst[0], incorrect",
                 v1[0], v2[0], precision );
  ensure_equals( "add non-empty vector to empty one, src[1]==dst[1], incorrect",
                 v1[1], v2[1], precision );
  ensure_equals( "add non-empty vector to empty one, src[2]==dst[2], incorrect",
                 v1[2], v2[2], precision );

  // add empty vector to non-empty one: throw in DEBUG to warn on no-op
  // skipped in RELEASE mode, would yield segmentation fault
  #ifndef NDEBUG        // exception only thrown in DEBUG mode
  try {
    std::vector< tk::real > r1{{ 4.0, 9.0, 2.0 }}, r2;
    // cppcheck-suppress unreadVariable
    r1 += r2;
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif

  // add empty vector to empty one: throw in DEBUG to warn on no-op
  // skipped in RELEASE mode, would yield segmentation fault
  #ifndef NDEBUG        // exception only thrown in DEBUG mode
  try {
    std::vector< tk::real > q1, q2;
    // cppcheck-suppress unreadVariable
    q1 += q2;
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif

  // add non-empty vector to non-empty one with src.size() == dst.size():
  // dst += src for all components, leave src unchanged
  std::vector< tk::real > s1{{ 4.0, 9.0, 2.0 }}, s2{{ 3.0, -3.0, 1.0 }};
  s1 += s2;
  ensure_equals( "add non-empty vector to non-empty one, dst size incorrect",
                 s1.size(), 3 );
  ensure_equals( "add non-empty vector to non-empty one, src size incorrect",
                 s2.size(), 3 );
  ensure_equals( "add non-empty vector to non-empty one, dst[0], incorrect",
                 s1[0], 7.0, precision );
  ensure_equals( "add non-empty vector to non-empty one, dst[1], incorrect",
                 s1[1], 6.0, precision );
  ensure_equals( "add non-empty vector to non-empty one, dst[2], incorrect",
                 s1[2], 3.0, precision );
  ensure_equals( "add non-empty vector to non-empty one, src[0], incorrect",
                 s2[0], 3.0, precision );
  ensure_equals( "add non-empty vector to non-empty one, src[1], incorrect",
                 s2[1], -3.0, precision );
  ensure_equals( "add non-empty vector to non-empty one, src[2], incorrect",
                 s2[2], 1.0, precision );

  // add non-empty vector to non-empty one with src.size() > dst.size(): grow
  // dst to that of src.size() (padding with zeros), dst += src for all
  // components, leave src unchanged
  std::vector< tk::real > w1{{ 4.0, 9.0 }}, w2{{ 3.0, -3.0, 1.0 }};
  w1 += w2;
  ensure_equals( "add non-empty vector to non-empty one, dst size incorrect",
                 w1.size(), 3 );
  ensure_equals( "add non-empty vector to non-empty one, src size incorrect",
                 w2.size(), 3 );
  ensure_equals( "add non-empty vector to non-empty one, dst[0], incorrect",
                 w1[0], 7.0, precision );
  ensure_equals( "add non-empty vector to non-empty one, dst[1], incorrect",
                 w1[1], 6.0, precision );
  ensure_equals( "add non-empty vector to non-empty one, dst[2], incorrect",
                 w1[2], 1.0, precision );
  ensure_equals( "add non-empty vector to non-empty one, src[0], incorrect",
                 w2[0], 3.0, precision );
  ensure_equals( "add non-empty vector to non-empty one, src[1], incorrect",
                 w2[1], -3.0, precision );
  ensure_equals( "add non-empty vector to non-empty one, src[2], incorrect",
                 w2[2], 1.0, precision );

  // add non-empty vector to non-empty one with src.size() < dst.size(): thrown
  // in DEBUG to warn on loosing data
  // skipped in RELEASE mode, would yield segmentation fault
  #ifndef NDEBUG        // exception only thrown in DEBUG mode
  try {
    std::vector< tk::real > n1{{ 4.0, 9.0, 2.0 }}, n2{{ 3.0, -3.0 }};
    // cppcheck-suppress unreadVariable
    n1 += n2;
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test keyEqual()
template<> template<>
void ContainerUtil_object::test< 7 >() {
  set_test_name( "keyEqual" );

  // Test if keys are equal
  std::map< int, tk::real > t1{ {1,4.0}, {2,2.0} }, t2{ {1,4.0}, {2,3.0} };
  ensure_equals( "keys are not equal", tk::keyEqual(t1,t2), true );

  // Test if keys are unequal
  std::map< int, tk::real > q1{ {3,4.0}, {2,2.0} }, q2{ {1,4.0}, {2,3.0} };
  ensure_equals( "keys are equal", tk::keyEqual(q1,q2), false );

  // test if throws in DEBUG to warn on unequal-size containers
  // skipped in RELEASE mode, would yield segmentation fault
  #ifndef NDEBUG        // exception only thrown in DEBUG mode
  try {
    std::map< int, tk::real > r1{ {1,4.0}, {2,2.0} }, r2{ {1,4.0} };
    tk::keyEqual( r1, r2 );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test sumsize()
template<> template<>
void ContainerUtil_object::test< 8 >() {
  set_test_name( "sumsize" );

  // Test sum of the sizes of a vector of vectors
  std::vector< std::vector< tk::real > > c{ {3,4.0}, {2,2.0} };
  ensure_equals( "sum of sizes incorrect", tk::sumsize(c), 4 );
}

//! Test destroy()
template<> template<>
void ContainerUtil_object::test< 9 >() {
  set_test_name( "destroy" );

  // Test destroying a vector of vectors
  std::vector< std::vector< tk::real > > c{ {3.2,4.0}, {2.1,2.0} };
  tk::destroy( c );
  ensure_equals( "destroy yields empty vec of vec", c.empty(), true );

  // Test destroying a map of vectors
  std::map< int, std::vector< tk::real > > m{ {3,{4.0,5.0}}, {2,{1.0,2.0}} };
  tk::destroy( m );
  ensure_equals( "destroy yields empty map of vec", m.empty(), true );
}

//! Test numunique()
template<> template<>
void ContainerUtil_object::test< 10 >() {
  set_test_name( "numunique" );

  // Test a vector of vector of ints
  std::vector< std::vector< int > > v{ {3,-4}, {2,3,6} };
  ensure_equals(
    "number of unique values in container of containers incorrect",
    tk::numunique( v ), 4 );

  // Test a set of vector of std::size_ts
  std::set< std::vector< std::size_t > > s{ {1,2}, {3,4,5} };
  ensure_equals(
    "number of unique values in container of containers incorrect",
    tk::numunique( s ), 5 );

  // Test a vector of set of std::size_ts
  std::vector< std::set< std::size_t > > q{ {1,2}, {1,4,5} };
  ensure_equals(
    "number of unique values in container of containers incorrect",
    tk::numunique( q ), 4 );
}

//! Test erase_if()
template<> template<>
void ContainerUtil_object::test< 11 >() {
  set_test_name( "erase_if" );

  std::vector< int > v{ 3,-4,2,3,6 };
  tk::erase_if( v, []( int i ){ return i%2; } );
  ensure( "erase_if on vector incorrect", v == std::vector<int>{-4,2,6} );

  std::vector< int > w{ 3,-4,2,3,6 };
  tk::erase_if( w, []( int i ){ return !(i%2); } );
  ensure( "erase_if on vector incorrect", w == std::vector<int>{3,3} );

  std::map< int, std::vector< std::size_t > > b{ {3,{4,5,6}}, {-2,{3,2,3}} };
  std::map< int, std::vector< std::size_t > > correct_result{ {3,{4,5,6}} };
  tk::erase_if( b, []( decltype(b)::value_type& p ){ return p.first<0; } );
  ensure( "erase_if on map incorrect", b == correct_result );
}

//! Test operator += adding values of a std::vector to another one
template<> template<>
void ContainerUtil_object::test< 12 >() {
  set_test_name( "operator+= std::array" );

  using tk::operator+=;

  // add array to other one, dst += src for all components, leave src unchanged
  std::array< tk::real, 3 > s1{{ 4.0, 9.0, 2.0 }}, s2{{ 3.0, -3.0, 1.0 }};
  s1 += s2;
  ensure_equals( "add array to other one, dst[0], incorrect",
                 s1[0], 7.0, precision );
  ensure_equals( "add array to other one, dst[1], incorrect",
                 s1[1], 6.0, precision );
  ensure_equals( "add array to other one, dst[2], incorrect",
                 s1[2], 3.0, precision );
  ensure_equals( "add array to other one, src[0], incorrect",
                 s2[0], 3.0, precision );
  ensure_equals( "add array to other one, src[1], incorrect",
                 s2[1], -3.0, precision );
  ensure_equals( "add array to other one, src[2], incorrect",
                 s2[2], 1.0, precision );
}

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
