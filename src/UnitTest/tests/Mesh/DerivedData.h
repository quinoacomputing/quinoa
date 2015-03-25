//******************************************************************************
/*!
  \file      src/UnitTest/tests/Mesh/DerivedData.h
  \author    J. Bakosi
  \date      Wed 25 Mar 2015 12:43:50 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unit tests for Mesh/DerivedData
  \details   Unit tests for Mesh/DerivedData
*/
//******************************************************************************
#ifndef test_DerivedData_h
#define test_DerivedData_h

#include <tut/tut.hpp>
#include <DerivedData.h>

namespace tut {

//! All tests in group inherited from this base
struct DerivedData_common {};

//! Test group shortcuts
using DerivedData_group = test_group< DerivedData_common >;
using DerivedData_object = DerivedData_group::object;

//! Define test group
DerivedData_group DerivedData( "Mesh/DerivedData" );

//! Test definitions for group

//! Attempt to shift empty container using shiftToZero
template<> template<>
void DerivedData_object::test< 1 >() {
  set_test_name( "shiftToZero graceful on empty container" );

  // Attempt to shift node IDs on emptycontainer. If some error happens or an
  // exception is throw that will go to the screen; no further tests are
  // necessary.
  std::vector< int > empty;
  tk::shiftToZero( empty );
}

//! Shift node ids to zero in line mesh
template<> template<>
void DerivedData_object::test< 2 >() {
  set_test_name( "shiftToZero for lines" );

  // Mesh connectivity for simple line-only mesh
  std::vector< int > inpoel { 1, 2,
                              2, 3,
                              3, 4,
                              4, 1,
                              5, 6,
                              6, 7,
                              7, 8,
                              8, 5,
                              1, 5,
                              2, 6,
                              3, 7,
                              4, 8 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Test new extents of node IDs in element connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  ensure_equals( "node ids should start from zero", *minmax.first, 0 );
}

//! Shift node ids to zero in triangle mesh
template<> template<>
void DerivedData_object::test< 3 >() {
  set_test_name( "shiftToZero for triangles" );

  // Mesh connectivity for simple triangle-only mesh
  std::vector< int > inpoel { 2,   9,  3,
                              1,   9,  2,
                              3,   9,  4,
                              1,   4,  9,
                              7,  10,  6,
                              6,  10,  5,
                              8,  10,  7,
                              10,  8,  5,
                              6,  11,  2,
                              2,  11,  1,
                              11,  6,  5,
                              11,  5,  1,
                              3,  12,  2,
                              12,  6,  2,
                              7,  12,  3,
                              12,  7,  6,
                              13,  7,  3,
                              4,  13,  3,
                              13,  8,  7,
                              8,  13,  4,
                              5,   8, 14,
                              1,   5, 14,
                              4,  14,  8,
                              1,  14,  4 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Test new extents of node IDs in element connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  ensure_equals( "node ids should start from zero", *minmax.first, 0 );
}

//! Shift node ids to zero in tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 4 >() {
  set_test_name( "shiftToZero for tetrahedra" );

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< int > inpoel { 12, 14,  9, 11,
                              10, 14, 13, 12,
                              14, 13, 12,  9,
                              10, 14, 12, 11,
                              1,  14,  5, 11,
                              7,   6, 10, 12,
                              14,  8,  5, 10,
                              8,   7, 10, 13,
                              7,  13,  3, 12,
                              1,   4, 14,  9,
                              13,  4,  3,  9,
                              3,   2, 12,  9,
                              4,   8, 14, 13,
                              6,   5, 10, 11,
                              1,   2,  9, 11,
                              2,   6, 12, 11,
                              6,  10, 12, 11,
                              2,  12,  9, 11,
                              5,  14, 10, 11,
                              14,  8, 10, 13,
                              13,  3, 12,  9,
                              7,  10, 13, 12,
                              14,  4, 13,  9,
                              14,  1,  9, 11 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Test new extents of node IDs in element connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  ensure_equals( "node ids should start from zero", *minmax.first, 0 );
}

//! Attempt to generate elements surrounding points on empty connectivity
template<> template<>
void DerivedData_object::test< 5 >() {
  set_test_name( "genEsup throws on empty container" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  // Attempt to generate elements surrounding points on empty container
  try {
    std::vector< int > empty;
    tk::genEsup( empty, 4 );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genEsup if it throws on non-zero based element connectivity
template<> template<>
void DerivedData_object::test< 6 >() {
  set_test_name( "genEsup throws on non-positive nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    // Partial mesh non-zero based mesh connectivity for tetrahedron-mesh
    std::vector< int > inpoel { 12, 14,  9, 11,
                                14,  4, 13,  9 };

    // Attempt to generate elements surrounding points passing partial inpoel
    auto esup = tk::genEsup( inpoel, 0 );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Generate and test elements surrounding points for tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 7 >() {
  set_test_name( "genEsup for tetrahedra" );

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< int > inpoel { 12, 14,  9, 11,
                              10, 14, 13, 12,
                              14, 13, 12,  9,
                              10, 14, 12, 11,
                              1,  14,  5, 11,
                              7,   6, 10, 12,
                              14,  8,  5, 10,
                              8,   7, 10, 13,
                              7,  13,  3, 12,
                              1,   4, 14,  9,
                              13,  4,  3,  9,
                              3,   2, 12,  9,
                              4,   8, 14, 13,
                              6,   5, 10, 11,
                              1,   2,  9, 11,
                              2,   6, 12, 11,
                              6,  10, 12, 11,
                              2,  12,  9, 11,
                              5,  14, 10, 11,
                              14,  8, 10, 13,
                              13,  3, 12,  9,
                              7,  10, 13, 12,
                              14,  4, 13,  9,
                              14,  1,  9, 11 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Generate elements surrounding points
  auto esup = tk::genEsup( inpoel, 4 );

  // Generate correct solution for elements surrounding points
  std::map< std::size_t, std::vector< std::size_t > > correct_esup {
    { { 0 }, { 4, 9, 14, 23 } },
    { { 1 }, { 11, 14, 15, 17 } },
    { { 2 }, { 8, 10, 11, 20 } },
    { { 3 }, { 9, 10, 12, 22 } },
    { { 4 }, { 4, 6, 13, 18 } },
    { { 5 }, { 5, 13, 15, 16 } },
    { { 6 }, { 5, 7, 8, 21 } },
    { { 7 }, { 6, 7, 12, 19 } },
    { { 8 }, { 0, 2, 9, 10, 11, 14, 17, 20, 22, 23 } },
    { { 9 }, { 1, 3, 5, 6, 7, 13, 16, 18, 19, 21 } },
    { { 10 }, { 0, 3, 4, 13, 14, 15, 16, 17, 18, 23 } },
    { { 11 }, { 0, 1, 2, 3, 5, 8, 11, 15, 16, 17, 20, 21 } },
    { { 12 }, { 1, 2, 7, 8, 10, 12, 19, 20, 21, 22 } },
    { { 13 }, { 0, 1, 2, 3, 4, 6, 9, 12, 18, 19, 22, 23 } }
  };

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = static_cast< std::size_t >( *minmax.second + 1 );

  ensure_equals( "number of points in esup incorrect",
                 npoin, correct_esup.size() );

  // test generated derived data structure, elements surrounding points
  for (std::size_t p=0; p<npoin; ++p) {
    // extract element ids from generated elements surrounding point p
    std::vector< std::size_t > points;
    for (auto i=esup.second[p]+1; i<=esup.second[p+1]; ++i)
      points.push_back( esup.first[i] );
    // find correct element ids surrounding point p
    auto it = correct_esup.find( p );
    // test if element ids exist surrounding point p
    ensure( "node id '" + std::to_string(p) + "' generated into esup but not "
            "in correct esup",
            it != correct_esup.end() );
    // test if element ids surrounding point p are correct
    if (it != correct_esup.end())
      ensure( "element ids surrounding point '" + std::to_string(p) +
              "' incorrect", points == it->second );
  }
}

// genEsup should also be tested for triangles

//! Attempt to generate points surrounding points on empty connectivity
template<> template<>
void DerivedData_object::test< 8 >() {
  set_test_name( "genPsup throws on empty container" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  // Attempt to generate points surrounding points on empty container
  try {
    std::vector< int > inpoel { 0, 1, 2, 3 };
    std::vector< int > empty;
    tk::genPsup( empty, 4, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genPsup if it throws on non-zero based element connectivity
template<> template<>
void DerivedData_object::test< 9 >() {
  set_test_name( "genPsup throws on non-positive nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< int > inpoel { 0, 1, 2, 3 };
    // Attempt to generate elements surrounding points passing partial inpoel
    tk::genPsup( inpoel, 0, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genPsup if it throws on empty element surrounding points
template<> template<>
void DerivedData_object::test< 10 >() {
  set_test_name( "genPsup throws on empty esup" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< int > inpoel { 0, 1, 2, 3 };
    // Attempt to generate elements surrounding points passing partial inpoel
    tk::genPsup( inpoel, 4, {} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown, test ok
  }
  #endif
}

//! Generate and test points surrounding points for tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 11 >() {
  set_test_name( "genPsup for tetrahedra" );

  // mesh connectivity for simple tetrahedron-only mesh
  std::vector< int > inpoel { 12, 14,  9, 11,
                              10, 14, 13, 12,
                              14, 13, 12,  9,
                              10, 14, 12, 11,
                              1,  14,  5, 11,
                              7,   6, 10, 12,
                              14,  8,  5, 10,
                              8,   7, 10, 13,
                              7,  13,  3, 12,
                              1,   4, 14,  9,
                              13,  4,  3,  9,
                              3,   2, 12,  9,
                              4,   8, 14, 13,
                              6,   5, 10, 11,
                              1,   2,  9, 11,
                              2,   6, 12, 11,
                              6,  10, 12, 11,
                              2,  12,  9, 11,
                              5,  14, 10, 11,
                              14,  8, 10, 13,
                              13,  3, 12,  9,
                              7,  10, 13, 12,
                              14,  4, 13,  9,
                              14,  1,  9, 11 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Generate elements surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );

  // Generate correct solution for elements surrounding points
  std::map< std::size_t, std::vector< std::size_t > > correct_psup {
    { { 0 }, { 1, 3, 4, 8, 10, 13 } },
    { { 1 }, { 0, 2, 5, 8, 10, 11 } },
    { { 2 }, { 1, 3, 6, 8, 11, 12 } },
    { { 3 }, { 0, 2, 7, 8, 12, 13 } },
    { { 4 }, { 0, 5, 7, 9, 10, 13 } },
    { { 5 }, { 1, 4, 6, 9, 10, 11 } },
    { { 6 }, { 2, 5, 7, 9, 11, 12 } },
    { { 7 }, { 3, 4, 6, 9, 12, 13 } },
    { { 8 }, { 0, 1, 2, 3, 10, 11, 12, 13 } },
    { { 9 }, { 4, 5, 6, 7, 10, 11, 12, 13 } },
    { { 10 }, { 0, 1, 4, 5, 8, 9, 11, 13 } },
    { { 11 }, { 1, 2, 5, 6, 8, 9, 10, 12, 13 } },
    { { 12 }, { 2, 3, 6, 7, 8, 9, 11, 13 } },
    { { 13 }, { 0, 3, 4, 7, 8, 9, 10, 11, 12 } }
  };

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = static_cast< std::size_t >( *minmax.second + 1 );

  ensure_equals( "number of points in psup incorrect",
                 npoin, correct_psup.size() );

  // test generated derived data structure, elements surrounding points
  for (std::size_t p=0; p<npoin; ++p) {
    // extract element ids from generated elements surrounding point p
    std::vector< std::size_t > points;
    for (auto i=psup.second[p]+1; i<=psup.second[p+1]; ++i)
      points.push_back( psup.first[i] );
    // find correct element ids surrounding point p
    auto it = correct_psup.find( p );
    // test if element ids exist surrounding point p
    ensure( "node id '" + std::to_string(p) + "' generated into psup but not "
            "in correct psup",
            it != correct_psup.end() );
    // test if element ids surrounding point p are correct
    if (it != correct_psup.end())
      ensure( "point ids surrounding point '" + std::to_string(p) +
              "' incorrect", points == it->second );
  }
}

// genPsup should also be tested for triangles

} // tut::

#endif // test_DerivedData_h
