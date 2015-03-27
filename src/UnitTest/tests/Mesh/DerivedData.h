//******************************************************************************
/*!
  \file      src/UnitTest/tests/Mesh/DerivedData.h
  \author    J. Bakosi
  \date      Fri 27 Mar 2015 10:48:36 AM MDT
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
  set_test_name( "shiftToZero graceful with empty inpoel" );

  // Attempt to shift node IDs with empty connectivity. If some error happens or an
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

//! Attempt to generate elements surrounding points with empty connectivity
template<> template<>
void DerivedData_object::test< 5 >() {
  set_test_name( "genEsup throws with empty inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
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

//! Test genEsup if it throws on non-positive nodes per elements
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

  // this is more of a test on this test
  ensure_equals( "number of points in 'correct' esup incorrect",
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

//! Attempt to generate points surrounding points with empty connectivity
template<> template<>
void DerivedData_object::test< 8 >() {
  set_test_name( "genPsup throws with empty inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
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

//! Test genPsup if it throws on non-positive number of nodes per elements
template<> template<>
void DerivedData_object::test< 9 >() {
  set_test_name( "genPsup throws on non-positive nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< int > inpoel { 0, 1, 2, 3 };
    tk::genPsup( inpoel, 0, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genPsup if it throws with empty element surrounding points
template<> template<>
void DerivedData_object::test< 10 >() {
  set_test_name( "genPsup throws with empty esup" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< int > inpoel { 0, 1, 2, 3 };
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

  // Generate correct solution for points surrounding points
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

  // this is more of a test on this test
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

//! \brief Attempt to generate elements surrounding points of elements with
//!   empty connectivity
template<> template<>
void DerivedData_object::test< 12 >() {
  set_test_name( "genEsupel throws with empty inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< int > inpoel { 0, 1, 2, 3 };
    std::vector< int > empty;
    tk::genEsupel( empty, 4, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genEsupel if it throws on non-positive nodes per elements
template<> template<>
void DerivedData_object::test< 13 >() {
  set_test_name( "genEsupel throws on non-positive nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< int > inpoel { 0, 1, 2, 3 };
    tk::genEsupel( inpoel, 0, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genEsupel if it throws with empty element surrounding points
template<> template<>
void DerivedData_object::test< 14 >() {
  set_test_name( "genEsupel throws with empty esup" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< int > inpoel { 0, 1, 2, 3 };
    tk::genEsupel( inpoel, 4, {} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown, test ok
  }
  #endif
}

//! \brief Generate and test elements surrounding points of elements for
//!   tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 15 >() {
  set_test_name( "genEsupel for tetrahedra" );

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
  auto esupel = tk::genEsupel( inpoel, 4, tk::genEsup(inpoel,4) );

  // Generate correct solution for elements surrounding points of elements
  std::map< std::size_t, std::vector< std::size_t > > correct_esupel {
    { { 0 }, { 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
               20, 21, 22, 23 } },
    { { 1 }, { 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19,
               20, 21, 22, 23 } },
    { { 2 }, { 0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19,
               20, 21, 22, 23 } },
    { { 3 }, { 0, 1, 2, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19,
               20, 21, 22, 23 } },
    { { 4 }, { 0, 1, 2, 3, 6, 9, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23 } },
    { { 5 }, { 0, 1, 2, 3, 6, 7, 8, 11, 13, 15, 16, 17, 18, 19, 20, 21 } },
    { { 6 }, { 0, 1, 2, 3, 4, 5, 7, 9, 12, 13, 16, 18, 19, 21, 22, 23 } },
    { { 7 }, { 1, 2, 3, 5, 6, 8, 10, 12, 13, 16, 18, 19, 20, 21, 22 } },
    { { 8 }, { 0, 1, 2, 3, 5, 7, 10, 11, 12, 15, 16, 17, 19, 20, 21, 22 } },
    { { 9 }, { 0, 1, 2, 3, 4, 6, 10, 11, 12, 14, 17, 18, 19, 20, 22, 23 } },
    { { 10 }, { 0, 1, 2, 7, 8, 9, 11, 12, 14, 17, 19, 20, 21, 22, 23 } },
    { { 11 }, { 0, 1, 2, 3, 5, 8, 9, 10, 14, 15, 16, 17, 20, 21, 22, 23 } },
    { { 12 }, { 0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 18, 19, 20, 21, 22, 23 } },
    { { 13 }, { 0, 1, 3, 4, 5, 6, 7, 14, 15, 16, 17, 18, 19, 21, 23 } },
    { { 14 }, { 0, 2, 3, 4, 9, 10, 11, 13, 15, 16, 17, 18, 20, 22, 23 } },
    { { 15 }, { 0, 1, 2, 3, 4, 5, 8, 11, 13, 14, 16, 17, 18, 20, 21, 23 } },
    { { 16 }, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 14, 15, 17, 18, 19, 20, 21,
                23 } },
    { { 17 }, { 0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 13, 14, 15, 16, 18, 20, 21, 22,
                23 } },
    { { 18 }, { 0, 1, 2, 3, 4, 5, 6, 7, 9, 12, 13, 14, 15, 16, 17, 19, 21, 22,
                23 } },
    { { 19 }, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 16, 18, 20, 21, 22,
                23 } },
    { { 20 }, { 0, 1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 19, 21, 22,
                23 } },
    { { 21 }, { 0, 1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20,
                22 } },
    { { 22 }, { 0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 14, 17, 18, 19, 20, 21,
                23 } },
    { { 23 }, { 0, 1, 2, 3, 4, 6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                22 } }
  };

  // find out number of elements from mesh connectivity
  auto nelem = inpoel.size()/4;

  // this is more of a test on this test
  ensure_equals( "number of elements in esupel incorrect",
                 nelem, correct_esupel.size() );

  // test generated derived data structure, elements surrounding points of
  // elements
  for (std::size_t e=0; e<nelem; ++e) {
    // extract element ids from generated elements surrounding points of
    // elements
    std::vector< std::size_t > elements;
    for (auto i=esupel.second[e]+1; i<=esupel.second[e+1]; ++i)
      elements.push_back( esupel.first[i] );
    // find correct element ids surrounding points of elements e
    auto it = correct_esupel.find( e );
    // test if element ids exist surrounding points of element e
    ensure( "elem id '" + std::to_string(e) + "' generated into esupel but not "
            "in correct esupel",
            it != correct_esupel.end() );
    // test if element ids surrounding points of element e are correct
    if (it != correct_esupel.end())
      ensure( "elem ids surrounding points of element '" + std::to_string(e) +
              "' incorrect", elements == it->second );
  }
}

// genEsupel should also be tested for triangles

//! Attempt to generate elements surrounding elements with empty connectivity
template<> template<>
void DerivedData_object::test< 16 >() {
  set_test_name( "genEsuel throws with empty inpoel" );

  try {
    std::vector< int > inpoel { 0, 1, 2, 3 };
    std::vector< int > empty;
    tk::genEsuel( empty, tk::genEsupel( inpoel, 4, tk::genEsup(inpoel,4)) );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, genEsuel still graceful, test ok
  }
}

//! Attempt to generate elements surrounding elements with non-tet connectivity
template<> template<>
void DerivedData_object::test< 17 >() {
  set_test_name( "genEsuel throws with non-tet inpoel" );

  try {
    std::vector< int > inpoel { 0, 1, 2 };  // size non-divisible by 4 = non-tet
    tk::genEsuel( inpoel, tk::genEsupel( inpoel, 4, tk::genEsup(inpoel,4)) );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, genEsuel still graceful, test ok
  }
}

//! Test genEsuel if it throws with empty element surrounding points of elements
template<> template<>
void DerivedData_object::test< 18 >() {
  set_test_name( "genEsuel throws with empty esupel" );

 #ifdef NDEBUG        // exception only thrown in DEBUG mode
   skip( "in RELEASE mode, would yield segmentation fault" );
 #else
  try {
    std::vector< int > inpoel { 0, 1, 2, 3 };
    tk::genEsuel( inpoel, {} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown, test ok
  }
  #endif
}

//! Generate and test elements surrounding elements for tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 19 >() {
  set_test_name( "genEsuel for tetrahedra" );

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
  auto esuel = tk::genEsuel( inpoel,
                             tk::genEsupel( inpoel, 4, tk::genEsup(inpoel,4)) );

  // Generate correct solution for elements surrounding points of elements
  std::map< std::size_t, std::vector< long int > > correct_esuel {
    { { 0 }, { 2, 3, 17, 23 } },
    { { 1 }, { 2, 3, 19, 21 } },
    { { 2 }, { 0, 1, 20, 22 } },
    { { 3 }, { 0, 1, 16, 18 } },
    { { 4 }, { 18, 23, -1, -1 } },
    { { 5 }, { 16, 21, -1, -1 } },
    { { 6 }, { 18, 19, -1, -1 } },
    { { 7 }, { 19, 21, -1, -1 } },
    { { 8 }, { 20, 21, -1, -1 } },
    { { 9 }, { 22, 23, -1, -1 } },
    { { 10 }, { 20, 22, -1, -1 } },
    { { 11 }, { 17, 20, -1, -1 } },
    { { 12 }, { 19, 22, -1, -1 } },
    { { 13 }, { 16, 18, -1, -1 } },
    { { 14 }, { 17, 23, -1, -1 } },
    { { 15 }, { 16, 17, -1, -1 } },
    { { 16 }, { 3, 5, 13, 15 } },
    { { 17 }, { 0, 11, 14, 15 } },
    { { 18 }, { 3, 4, 6, 13 } },
    { { 19 }, { 1, 6, 7, 12 } },
    { { 20 }, { 2, 8, 10, 11 } },
    { { 21 }, { 1, 5, 7, 8 } },
    { { 22 }, { 2, 9, 10, 12 } },
    { { 23 }, { 0, 4, 9, 14 } }
  };

  // find out number of elements from mesh connectivity
  auto nelem = inpoel.size()/4;

  // this is more of a test on this test
  ensure_equals( "number of elements in esuel incorrect",
                 nelem, correct_esuel.size() );

  // test generated derived data structure, elements surrounding elements
  for (std::size_t e=0; e<nelem; ++e) {
    // extract element ids from generated elements surrounding elements
    std::vector< long int > elements;
    for (std::size_t i=0; i<4; ++i) elements.push_back( esuel[e*4+i] );
    // find correct element ids surrounding elements e
    auto it = correct_esuel.find( e );
    // test if element ids exist surrounding element e
    ensure( "element id '" + std::to_string(e) + "' generated into esuel but not "
            "in correct esuel",
            it != correct_esuel.end() );
    // test if element ids surrounding element e are correct
    if (it != correct_esuel.end()) {
      ensure_equals( "number of elements surrounding element " +
                     std::to_string(e) + " in generated esuel incorrect",
                     it->second.size(), 4 );
      ensure( "element ids surrounding element '" + std::to_string(e) +
              "' incorrect", elements == it->second );
    }
  }
}

//! Attempt to generate edge connectivity with empty connectivity
template<> template<>
void DerivedData_object::test< 20 >() {
  set_test_name( "genInpoed throws with empty inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< int > inpoel { 0, 1, 2, 3 };
    std::vector< int > empty;
    tk::genInpoed( empty, 4, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genInpoed if it throws on non-positive number of nodes per elements
template<> template<>
void DerivedData_object::test< 21 >() {
  set_test_name( "genInpoed throws on non-positive nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< int > inpoel { 0, 1, 2, 3 };
    tk::genInpoed( inpoel, 0, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genInpoed if it throws with empty element surrounding points
template<> template<>
void DerivedData_object::test< 22 >() {
  set_test_name( "genInpoed throws with empty esup" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< int > inpoel { 0, 1, 2, 3 };
    tk::genPsup( inpoel, 4, {} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& e ) {
    // exception thrown, test ok
  }
  #endif
}

//! Generate and test edge connectivity for tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 23 >() {
  set_test_name( "genInpoed for tetrahedra" );

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

  // Generate edge connectivity
  auto inpoed = tk::genInpoed( inpoel, 4, tk::genEsup(inpoel,4) );

  // Generate correct solution for edge connectivity
  std::vector< std::size_t > correct_inpoed {
    0, 1, 0, 3, 0, 4, 0, 8, 0, 10, 0, 13,
    1, 2, 1, 5, 1, 8, 1, 10, 1, 11,
    2, 3, 2, 6, 2, 8, 2, 11, 2, 12,
    3, 7, 3, 8, 3, 12, 3, 13,
    4, 5, 4, 7, 4, 9, 4, 10, 4, 13,
    5, 6, 5, 9, 5, 10, 5, 11,
    6, 7, 6, 9, 6, 11, 6, 12,
    7, 9, 7, 12, 7, 13,
    8, 10, 8, 11, 8, 12, 8, 13,
    9, 10, 9, 11, 9, 12, 9, 13,
    10, 11, 10, 13,
    11, 12, 11, 13,
    12, 13 };

  // this is more of a test on this test
  ensure_equals( "number of edges in 'correct' edge connectivity non-divisble "
                 "by 2",
                 correct_inpoed.size() % 2, 0 );

  // test if edge connectivity is the correct size
  ensure_equals( "number of edges in generated edge connectivity non-divisble "
                 "by 2",
                 inpoed.size() % 2, 0 );
  ensure_equals( "number of edges in edge connectivity incorrect",
                 inpoed.size(), correct_inpoed.size() );

  // this if edge connectivity correct
  ensure( "generatd edge connectivity incorrect", inpoed == correct_inpoed );
}

// genInpoed should also be tested for triangles

} // tut::

#endif // test_DerivedData_h
