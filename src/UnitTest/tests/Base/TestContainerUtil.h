// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestContainerUtil.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Base/ContainerUtil.h
  \details   Unit tests for Base/ContainerUtil.h
*/
// *****************************************************************************
#ifndef test_ContainerUtil_h
#define test_ContainerUtil_h

#include <map>
#include <vector>
#include <string>
#include <unordered_map>

#include <unistd.h>

#include "NoWarning/tut.h"

#include "ContainerUtil.h"
#include "StatCtr.h"

namespace tut {

//! All tests in group inherited from this base
struct ContainerUtil_common {
  double precision = 1.0e-15;    // required floating-point precision
};

//! Test group shortcuts
using ContainerUtil_group =
  test_group< ContainerUtil_common, MAX_TESTS_IN_GROUP >;
using ContainerUtil_object = ContainerUtil_group::object;

//! Define test group
static ContainerUtil_group ContainerUtil( "Base/ContainerUtil" );

//! Test definitions for group

//! Test unique making the elements of a container unique
//! \author J. Bakosi
template<> template<>
void ContainerUtil_object::test< 1 >() {
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
//! \author J. Bakosi
template<> template<>
void ContainerUtil_object::test< 2 >() {
  set_test_name( "[c]ref_find" );

  std::map< int, std::string > m{ {1,"one"}, {2,"two"} };
  std::unordered_map< int, std::string > u{ {1,"one"}, {2,"two"} };

  const auto one = tk::cref_find(m,1);
  ensure_equals( "cref_find incorrect", one, "one" );
  auto two = tk::ref_find(u,2);
  ensure_equals( "ref_find incorrect", two, "two" );

  // test if cref_find throws in DEBUG when cannot find key
  try {
    auto three = tk::ref_find(u,3);
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }
}

//! Test vector extents
//! \author J. Bakosi
template<> template<>
void ContainerUtil_object::test< 3 >() {
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
//! \author J. Bakosi
template<> template<>
void ContainerUtil_object::test< 4 >() {
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
//! \author J. Bakosi
template<> template<>
void ContainerUtil_object::test< 5 >() {
  set_test_name( "operator+=" );

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
  try {
    std::vector< tk::real > r1{{ 4.0, 9.0, 2.0 }}, r2;
    r1 += r2;
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  // add empty vector to empty one: throw in DEBUG to warn on no-op
  try {
    std::vector< tk::real > q1, q2;
    q1 += q2;
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

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
  try {
    std::vector< tk::real > n1{{ 4.0, 9.0, 2.0 }}, n2{{ 3.0, -3.0 }};
    n1 += n2;
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }
}

//! Test keyEqual()
//! \author J. Bakosi
template<> template<>
void ContainerUtil_object::test< 6 >() {
  set_test_name( "keyEqual" );

  // test if throws in DEBUG to warn on unequal-size containers
  try {
    std::map< int, tk::real > r1{ {1,4.0}, {2,2.0} }, r2{ {1,4.0} };
    tk::keyEqual( r1, r2 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, test ok
  }

  // Test if keys are equal
  std::map< int, tk::real > r1{ {1,4.0}, {2,2.0} }, r2{ {1,4.0}, {2,3.0} };
  ensure_equals( "keys are not equal", tk::keyEqual(r1,r2), true );
  
  // Test if keys are unequal
  std::map< int, tk::real > q1{ {3,4.0}, {2,2.0} }, q2{ {1,4.0}, {2,3.0} };
  ensure_equals( "keys are equal", tk::keyEqual(q1,q2), false );
}

//! Test sumsize()
//! \author J. Bakosi
template<> template<>
void ContainerUtil_object::test< 7 >() {
  set_test_name( "sumsize" );

  // Test sum of the sizes of a vector of vectors
  std::vector< std::vector< tk::real > > c{ {3,4.0}, {2,2.0} };
  ensure_equals( "sum of sizes incorrect", tk::sumsize(c), 4 );
}

//! Test destroy()
//! \author J. Bakosi
template<> template<>
void ContainerUtil_object::test< 8 >() {
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

} // tut::

#endif // test_ContainerUtil_h
