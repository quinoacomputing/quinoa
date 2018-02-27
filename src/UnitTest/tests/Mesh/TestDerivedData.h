// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Mesh/TestDerivedData.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Unit tests for Mesh/DerivedData
  \details   Unit tests for Mesh/DerivedData. All unit tests start from simple
     mesh connectivities defined in the code. The tetrahedron mesh in Gmsh ASCII
     format is as follows. Note that ids start from zero in the code, but from
     one in Gmsh.
     \code{.sh}
       $MeshFormat
       2.2 0 8
       $EndMeshFormat
       $Nodes
       14
       1 0 0 0
       2 1 0 0
       3 1 1 0
       4 0 1 0
       5 0 0 1
       6 1 0 1
       7 1 1 1
       8 0 1 1
       9 0.5 0.5 0
       10 0.5 0.5 1
       11 0.5 0 0.5
       12 1 0.5 0.5
       13 0.5 1 0.5
       14 0 0.5 0.5
       $EndNodes
       $Elements
       24
       1 4 1 0 12 14 9 11
       2 4 1 0 10 14 13 12
       3 4 1 0 14 13 12 9
       4 4 1 0 10 14 12 11
       5 4 1 0 1 14 5 11
       6 4 1 0 7 6 10 12
       7 4 1 0 14 8 5 10
       8 4 1 0 8 7 10 13
       9 4 1 0 7 13 3 12
       10 4 1 0 1 4 14 9
       11 4 1 0 13 4 3 9
       12 4 1 0 3 2 12 9
       13 4 1 0 4 8 14 13
       14 4 1 0 6 5 10 11
       15 4 1 0 1 2 9 11
       16 4 1 0 2 6 12 11
       17 4 1 0 6 10 12 11
       18 4 1 0 2 12 9 11
       19 4 1 0 5 14 10 11
       20 4 1 0 14 8 10 13
       21 4 1 0 13 3 12 9
       22 4 1 0 7 10 13 12
       23 4 1 0 14 4 13 9
       24 4 1 0 14 1 9 11
       $EndElements
     \endcode
     Here is the simple triangle mesh used below by the unit tests in Gmsh ASCII
     format. Note that ids start from zero in the code, but from one in Gmsh.
     \code{.sh}
       $MeshFormat
       2.2 0 8
       $EndMeshFormat
       $Nodes
       14
       1 0 0 0
       2 1 0 0
       3 1 1 0
       4 0 1 0
       5 0 0 1
       6 1 0 1
       7 1 1 1
       8 0 1 1
       9 0.5 0.5 0
       10 0.5 0.5 1
       11 0.5 0 0.5
       12 1 0.5 0.5
       13 0.5 1 0.5
       14 0 0.5 0.5
       $EndNodes
       $Elements
       24
       1 2 2 0 1 1 9 2
       2 2 2 0 1 1 4 9
       3 2 2 0 1 2 9 3
       4 2 2 0 1 3 9 4
       5 2 2 0 2 5 6 10
       6 2 2 0 2 5 10 8
       7 2 2 0 2 6 7 10
       8 2 2 0 2 7 8 10
       9 2 2 0 3 1 2 11
       10 2 2 0 3 1 11 5
       11 2 2 0 3 2 6 11
       12 2 2 0 3 5 11 6
       13 2 2 0 4 2 3 12
       14 2 2 0 4 2 12 6
       15 2 2 0 4 3 7 12
       16 2 2 0 4 6 12 7
       17 2 2 0 5 3 4 13
       18 2 2 0 5 3 13 7
       19 2 2 0 5 4 8 13
       20 2 2 0 5 7 13 8
       21 2 2 0 6 1 14 4
       22 2 2 0 6 1 5 14
       23 2 2 0 6 4 14 8
       24 2 2 0 6 5 8 14
       $EndElements
     \endcode
*/
// *****************************************************************************
#ifndef test_DerivedData_h
#define test_DerivedData_h

#include "NoWarning/tut.h"

#include "DerivedData.h"

namespace tut {

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif

//! All tests in group inherited from this base
struct DerivedData_common {};

// Test group shortcuts
// The 2nd template argument is the max number of tests in this group. If
// omitted, the default is 50, specified in tut/tut.hpp.
using DerivedData_group = test_group< DerivedData_common, MAX_TESTS_IN_GROUP >;
using DerivedData_object = DerivedData_group::object;

//! Define test group
static DerivedData_group DerivedData( "Mesh/DerivedData" );

//! Test definitions for group

//! Attempt to generate elements surrounding points with empty connectivity
template<> template<>
void DerivedData_object::test< 1 >() {
  set_test_name( "genEsup throws with empty inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > empty;
    tk::genEsup( empty, 4 );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genEsup if it throws on non-positive nodes per elements
template<> template<>
void DerivedData_object::test< 2 >() {
  set_test_name( "genEsup throws on non-positive nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    // Partial mesh non-zero based mesh connectivity for tetrahedron-mesh
    std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                        14,  4, 13,  9 };
    auto esup = tk::genEsup( inpoel, 0 );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! \brief Test genEsup if it throws on inpoel non-divisible by the number of
//!   nodes per elements
template<> template<>
void DerivedData_object::test< 3 >() {
  set_test_name( "genEsup throws on inpoel non-div nnpe" );

  try {
    // Partial mesh mesh connectivity
    std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                        14,  4, 13 };
    tk::genEsup( inpoel, 4 );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, genEsup still graceful, test ok
  }
}

//! Generate and test elements surrounding points for tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 4 >() {
  set_test_name( "genEsup for tetrahedra" );

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12, 14,  9, 11,
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

  // find out number of points from mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

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

//! Generate and test elements surrounding points for triangle-only mesh
template<> template<>
void DerivedData_object::test< 5 >() {
  set_test_name( "genEsup for triangles" );

  // Mesh connectivity for simple triangle-only mesh
  std::vector< std::size_t > inpoel { 1,  9,  2,
                                      1,  4,  9,
                                      2,  9,  3,
                                      3,  9,  4,
                                      5,  6, 10,
                                      5, 10,  8,
                                      6,  7, 10,
                                      7,  8, 10,
                                      1,  2, 11,
                                      1, 11,  5,
                                      2,  6, 11,
                                      5, 11,  6,
                                      2,  3, 12,
                                      2, 12,  6,
                                      3,  7, 12,
                                      6, 12,  7,
                                      3,  4, 13,
                                      3, 13,  7,
                                      4,  8, 13,
                                      7, 13,  8,
                                      1, 14,  4,
                                      1,  5, 14,
                                      4, 14,  8,
                                      5,  8, 14 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Generate elements surrounding points
  auto esup = tk::genEsup( inpoel, 3 );

  // Generate correct solution for elements surrounding points
  std::map< std::size_t, std::vector< std::size_t > > correct_esup {
    { { 0 }, { 0, 1, 8, 9, 20, 21 } },
    { { 1 }, { 0, 2, 8, 10, 12, 13 } },
    { { 2 }, { 2, 3, 12, 14, 16, 17 } },
    { { 3 }, { 1, 3, 16, 18, 20, 22 } },
    { { 4 }, { 4, 5, 9, 11, 21, 23 } },
    { { 5 }, { 4, 6, 10, 11, 13, 15 } },
    { { 6 }, { 6, 7, 14, 15, 17, 19 } },
    { { 7 }, { 5, 7, 18, 19, 22, 23 } },
    { { 8 }, { 0, 1, 2, 3 } },
    { { 9 }, { 4, 5, 6, 7 } },
    { { 10 }, { 8, 9, 10, 11 } },
    { { 11 }, { 12, 13, 14, 15 } },
    { { 12 }, { 16, 17, 18, 19 } },
    { { 13 }, { 20, 21, 22, 23 } }
  };

  // find out number of points from mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

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

//! Attempt to generate points surrounding points with empty connectivity
template<> template<>
void DerivedData_object::test< 6 >() {
  set_test_name( "genPsup throws with empty inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    std::vector< std::size_t > empty;
    tk::genPsup( empty, 4, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genPsup if it throws on non-positive number of nodes per elements
template<> template<>
void DerivedData_object::test< 7 >() {
  set_test_name( "genPsup throws on non-positive nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genPsup( inpoel, 0, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genPsup if it throws with empty element surrounding points
template<> template<>
void DerivedData_object::test< 8 >() {
  set_test_name( "genPsup throws with empty esup.first" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genPsup( inpoel, 4, {} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! Test genPsup if it throws with empty element surrounding points
template<> template<>
void DerivedData_object::test< 9 >() {
  set_test_name( "genPsup throws with empty esup.second" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genPsup( inpoel, 4, {{1},{}} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! \brief Test genPsup if it throws on inpoel non-divisible by the number of
//!   nodes per elements
template<> template<>
void DerivedData_object::test< 10 >() {
  set_test_name( "genPsup throws on inpoel non-div nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield invalid read" );
  #else
  try {
    // Partial mesh mesh connectivity
    std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                        14,  4, 13 };
    tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Generate and test points surrounding points for tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 11 >() {
  set_test_name( "genPsup for tetrahedra" );

  // mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12, 14,  9, 11,
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

  // find out number of points from mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

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

//! Generate and test points surrounding points for triangle-only mesh
template<> template<>
void DerivedData_object::test< 12 >() {
  set_test_name( "genPsup for triangles" );

  // Mesh connectivity for simple triangle-only mesh
  std::vector< std::size_t > inpoel { 1,  9,  2,
                                      1,  4,  9,
                                      2,  9,  3,
                                      3,  9,  4,
                                      5,  6, 10,
                                      5, 10,  8,
                                      6,  7, 10,
                                      7,  8, 10,
                                      1,  2, 11,
                                      1, 11,  5,
                                      2,  6, 11,
                                      5, 11,  6,
                                      2,  3, 12,
                                      2, 12,  6,
                                      3,  7, 12,
                                      6, 12,  7,
                                      3,  4, 13,
                                      3, 13,  7,
                                      4,  8, 13,
                                      7, 13,  8,
                                      1, 14,  4,
                                      1,  5, 14,
                                      4, 14,  8,
                                      5,  8, 14 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Generate elements surrounding points
  auto psup = tk::genPsup( inpoel, 3, tk::genEsup(inpoel,3) );

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
    { { 8 }, { 0, 1, 2, 3 } },
    { { 9 }, { 4, 5, 6, 7 } },
    { { 10 }, { 0, 1, 4, 5 } },
    { { 11 }, { 1, 2, 5, 6 } },
    { { 12 }, { 2, 3, 6, 7 } },
    { { 13 }, { 0, 3, 4, 7 } }
  };

  // find out number of points from mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

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

//! \brief Attempt to generate elements surrounding points of elements with
//!   empty connectivity
template<> template<>
void DerivedData_object::test< 13 >() {
  set_test_name( "genEsupel throws with empty inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    std::vector< std::size_t > empty;
    tk::genEsupel( empty, 4, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genEsupel if it throws on non-positive nodes per elements
template<> template<>
void DerivedData_object::test< 14 >() {
  set_test_name( "genEsupel throws on non-positive nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEsupel( inpoel, 0, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genEsupel if it throws with empty element surrounding points . first
template<> template<>
void DerivedData_object::test< 15 >() {
  set_test_name( "genEsupel throws with empty esup.first" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEsupel( inpoel, 4, {} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! Test genEsupel if it throws with empty element surrounding points . second
template<> template<>
void DerivedData_object::test< 16 >() {
  set_test_name( "genEsupel throws with empty esup.second" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEsupel( inpoel, 4, {{1},{}} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! \brief Test genEsupel if it throws on inpoel non-divisible by the number of
//!   nodes per elements
template<> template<>
void DerivedData_object::test< 17 >() {
  set_test_name( "genEsupel throws on inpoel non-div nnpe" );

  try {
    // Partial mesh mesh connectivity
    std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                        14,  4, 13 };
    tk::genEsupel( inpoel, 4, tk::genEsup(inpoel,4) );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, genEsupel still graceful, test ok
  }
}

//! \brief Generate and test elements surrounding points of elements for
//!   tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 18 >() {
  set_test_name( "genEsupel for tetrahedra" );

  // mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12, 14,  9, 11,
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

//! \brief Generate and test elements surrounding points of elements for
//!   triangle-only mesh
template<> template<>
void DerivedData_object::test< 19 >() {
  set_test_name( "genEsupel for triangles" );

  // Mesh connectivity for simple triangle-only mesh
  std::vector< std::size_t > inpoel { 1,  9,  2,
                                      1,  4,  9,
                                      2,  9,  3,
                                      3,  9,  4,
                                      5,  6, 10,
                                      5, 10,  8,
                                      6,  7, 10,
                                      7,  8, 10,
                                      1,  2, 11,
                                      1, 11,  5,
                                      2,  6, 11,
                                      5, 11,  6,
                                      2,  3, 12,
                                      2, 12,  6,
                                      3,  7, 12,
                                      6, 12,  7,
                                      3,  4, 13,
                                      3, 13,  7,
                                      4,  8, 13,
                                      7, 13,  8,
                                      1, 14,  4,
                                      1,  5, 14,
                                      4, 14,  8,
                                      5,  8, 14 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Generate elements surrounding points
  auto esupel = tk::genEsupel( inpoel, 3, tk::genEsup(inpoel,3) );

  // Generate correct solution for elements surrounding points of elements
  std::map< std::size_t, std::vector< std::size_t > > correct_esupel {
    { { 0 }, { 1, 2, 3, 8, 9, 10, 12, 13, 20, 21 } },
    { { 1 }, { 0, 2, 3, 8, 9, 16, 18, 20, 21, 22 } },
    { { 2 }, { 0, 1, 3, 8, 10, 12, 13, 14, 16, 17 } },
    { { 3 }, { 0, 1, 2, 12, 14, 16, 17, 18, 20, 22 } },
    { { 4 }, { 5, 6, 7, 9, 10, 11, 13, 15, 21, 23 } },
    { { 5 }, { 4, 6, 7, 9, 11, 18, 19, 21, 22, 23 } },
    { { 6 }, { 4, 5, 7, 10, 11, 13, 14, 15, 17, 19 } },
    { { 7 }, { 4, 5, 6, 14, 15, 17, 18, 19, 22, 23 } },
    { { 8 }, { 0, 1, 2, 9, 10, 11, 12, 13, 20, 21 } },
    { { 9 }, { 0, 1, 4, 5, 8, 10, 11, 20, 21, 23 } },
    { { 10 }, { 0, 2, 4, 6, 8, 9, 11, 12, 13, 15 } },
    { { 11 }, { 4, 5, 6, 8, 9, 10, 13, 15, 21, 23 } },
    { { 12 }, { 0, 2, 3, 8, 10, 13, 14, 15, 16, 17 } },
    { { 13 }, { 0, 2, 4, 6, 8, 10, 11, 12, 14, 15 } },
    { { 14 }, { 2, 3, 6, 7, 12, 13, 15, 16, 17, 19 } },
    { { 15 }, { 4, 6, 7, 10, 11, 12, 13, 14, 17, 19 } },
    { { 16 }, { 1, 2, 3, 12, 14, 17, 18, 19, 20, 22 } },
    { { 17 }, { 2, 3, 6, 7, 12, 14, 15, 16, 18, 19 } },
    { { 18 }, { 1, 3, 5, 7, 16, 17, 19, 20, 22, 23 } },
    { { 19 }, { 5, 6, 7, 14, 15, 16, 17, 18, 22, 23 } },
    { { 20 }, { 0, 1, 3, 8, 9, 16, 18, 21, 22, 23 } },
    { { 21 }, { 0, 1, 4, 5, 8, 9, 11, 20, 22, 23 } },
    { { 22 }, { 1, 3, 5, 7, 16, 18, 19, 20, 21, 23 } },
    { { 23 }, { 4, 5, 7, 9, 11, 18, 19, 20, 21, 22 } }
  };

  // find out number of elements from mesh connectivity
  auto nelem = inpoel.size()/3;

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

//! Attempt to generate elements surrounding elements with empty connectivity
template<> template<>
void DerivedData_object::test< 20 >() {
  set_test_name( "genEsuel throws with empty inpoel" );

  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    std::vector< std::size_t > empty;
    tk::genEsuel( empty, 4, tk::genEsup(inpoel,4) );
    #ifndef NDEBUG
    fail( "should throw exception in DEBUG mode" );
    #endif
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
    // Assert skipped in RELEASE mode, genEsuel still graceful, test ok
  }
}

//! Test genEsuel if it throws on non-positive nodes per elements
template<> template<>
void DerivedData_object::test< 21 >() {
  set_test_name( "genEsuel throws on non-positive nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEsuel( inpoel, 0, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genEsuel if it throws with empty element surrounding points .first
template<> template<>
void DerivedData_object::test< 22 >() {
  set_test_name( "genEsuel throws with empty esup.first" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEsuel( inpoel, 4, {} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! Test genEsuel if it throws with empty element surrounding points . second
template<> template<>
void DerivedData_object::test< 23 >() {
  set_test_name( "genEsuel throws with empty esup.second" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEsuel( inpoel, 4, {{1},{}} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! \brief Test genEsuel if it throws on inpoel non-divisible by the number of
//!   nodes per elements
template<> template<>
void DerivedData_object::test< 24 >() {
  set_test_name( "genEsuel throws on inpoel non-div nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield invalid read" );
  try {
    // Partial mesh mesh connectivity
    std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                        14,  4, 13 };
    tk::genEsuel( inpoel, 4, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Generate and test elements surrounding elements for tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 25 >() {
  set_test_name( "genEsuel for tetrahedra" );

  // mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12, 14,  9, 11,
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
  auto esuel = tk::genEsuel( inpoel, 4, tk::genEsup(inpoel,4) );

  // Generate correct solution for elements surrounding elements
  std::map< std::size_t, std::vector< std::size_t > > correct_esuel {
    { { 0 }, { 2, 3, 17, 23 } },
    { { 1 }, { 2, 3, 19, 21 } },
    { { 2 }, { 0, 1, 20, 22 } },
    { { 3 }, { 0, 1, 16, 18 } },
    { { 4 }, { 18, 23 } },
    { { 5 }, { 16, 21 } },
    { { 6 }, { 18, 19 } },
    { { 7 }, { 19, 21 } },
    { { 8 }, { 20, 21 } },
    { { 9 }, { 22, 23 } },
    { { 10 }, { 20, 22 } },
    { { 11 }, { 17, 20 } },
    { { 12 }, { 19, 22 } },
    { { 13 }, { 16, 18 } },
    { { 14 }, { 17, 23 } },
    { { 15 }, { 16, 17 } },
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
    std::vector< std::size_t > elements;
    for (auto i=esuel.second[e]+1; i<=esuel.second[e+1]; ++i)
      elements.push_back( esuel.first[i] );
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
                     it->second.size(), elements.size() );
      ensure( "element ids surrounding element '" + std::to_string(e) +
              "' incorrect", elements == it->second );
    }
  }
}

//! Generate and test elements surrounding elements for triangle-only mesh
template<> template<>
void DerivedData_object::test< 26 >() {
  set_test_name( "genEsuel for triangles" );

  // Mesh connectivity for simple triangle-only mesh
  std::vector< std::size_t > inpoel { 1,  9,  2,
                                      1,  4,  9,
                                      2,  9,  3,
                                      3,  9,  4,
                                      5,  6, 10,
                                      5, 10,  8,
                                      6,  7, 10,
                                      7,  8, 10,
                                      1,  2, 11,
                                      1, 11,  5,
                                      2,  6, 11,
                                      5, 11,  6,
                                      2,  3, 12,
                                      2, 12,  6,
                                      3,  7, 12,
                                      6, 12,  7,
                                      3,  4, 13,
                                      3, 13,  7,
                                      4,  8, 13,
                                      7, 13,  8,
                                      1, 14,  4,
                                      1,  5, 14,
                                      4, 14,  8,
                                      5,  8, 14 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Generate elements surrounding points
  auto esuel = tk::genEsuel( inpoel, 3, tk::genEsup(inpoel,3) );

  // Generate correct solution for elements surrounding points of elements
  std::map< std::size_t, std::vector< std::size_t > > correct_esuel {
    { { 0 }, { 1, 2, 8 } },
    { { 1 }, { 0, 3, 20 } },
    { { 2 }, { 0, 3, 12 } },
    { { 3 }, { 1, 2, 16 } },
    { { 4 }, { 5, 6, 11 } },
    { { 5 }, { 4, 7, 23 } },
    { { 6 }, { 4, 7, 15 } },
    { { 7 }, { 5, 6, 19 } },
    { { 8 }, { 0, 9, 10 } },
    { { 9 }, { 8, 11, 21 } },
    { { 10 }, { 8, 11, 13 } },
    { { 11 }, { 4, 9, 10 } },
    { { 12 }, { 2, 13, 14 } },
    { { 13 }, { 10, 12, 15 } },
    { { 14 }, { 12, 15, 17 } },
    { { 15 }, { 6, 13, 14 } },
    { { 16 }, { 3, 17, 18 } },
    { { 17 }, { 14, 16, 19 } },
    { { 18 }, { 16, 19, 22 } },
    { { 19 }, { 7, 17, 18 } },
    { { 20 }, { 1, 21, 22 } },
    { { 21 }, { 9, 20, 23 } },
    { { 22 }, { 18, 20, 23 } },
    { { 23 }, { 5, 21, 22 } }
  };

  // find out number of elements from mesh connectivity
  auto nelem = inpoel.size()/3;

  // this is more of a test on this test
  ensure_equals( "number of elements in esuel incorrect",
                 nelem, correct_esuel.size() );

  // test generated derived data structure, elements surrounding elements
  for (std::size_t e=0; e<nelem; ++e) {
    // extract element ids from generated elements surrounding elements
    std::vector< std::size_t > elements;
    for (auto i=esuel.second[e]+1; i<=esuel.second[e+1]; ++i)
      elements.push_back( esuel.first[i] );
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
                     it->second.size(), elements.size() );
      ensure( "element ids surrounding element '" + std::to_string(e) +
              "' incorrect", elements == it->second );
    }
  }
}

//! Attempt to generate edges surrounding points with empty connectivity
template<> template<>
void DerivedData_object::test< 27 >() {
  set_test_name( "genEdsup throws with empty inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    std::vector< std::size_t > empty;
    tk::genEdsup( empty, 4, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genEdsup if it throws on non-positive number of nodes per elements
template<> template<>
void DerivedData_object::test< 28 >() {
  set_test_name( "genEdsup throws on non-positive nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEdsup( inpoel, 0, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! \brief Test genEdsup if it throws on non-tet or non-tri number of nodes per
//!   elements
template<> template<>
void DerivedData_object::test< 29 >() {
  set_test_name( "genEdsup throws on unsupported nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEdsup( inpoel, 5, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genEdsup if it throws with empty element surrounding points
template<> template<>
void DerivedData_object::test< 30 >() {
  set_test_name( "genEdsup throws with empty esup.first" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEdsup( inpoel, 4, {} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! Test genEdsup if it throws with empty element surrounding points
template<> template<>
void DerivedData_object::test< 31 >() {
  set_test_name( "genEdsup throws with empty esup.second" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEdsup( inpoel, 4, {{1},{}} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! \brief Test genEdsup if it throws on inpoel non-divisible by the number of
//!   nodes per elements
template<> template<>
void DerivedData_object::test< 32 >() {
  set_test_name( "genEdsup throws on inpoel non-div nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield invalid read" );
  #else
  try {
    // Partial mesh mesh connectivity
    std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                        14,  4, 13 };
    tk::genEdsup( inpoel, 4, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Generate and test edges surrounding points for tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 33 >() {
  set_test_name( "genEdsup for tetrahedra" );

  // mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12, 14,  9, 11,
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

  // Generate edges surrounding points
  auto edsup = tk::genEdsup( inpoel, 4, tk::genEsup(inpoel,4) );

  auto& edsup1 = edsup.first;
  auto& edsup2 = edsup.second;

  // Generate correct solution for edges surrounding points
  std::map< std::size_t, std::vector< std::size_t > > correct_edsup {
    { {  0 }, { 1,   3,  4,  8, 10, 13 } },
    { {  1 }, { 2,   5,  8, 10, 11 } },
    { {  2 }, { 3,   6,  8, 11, 12 } },
    { {  3 }, { 7,   8, 12, 13 } },
    { {  4 }, { 5,   7,  9, 10, 13 } },
    { {  5 }, { 6,   9, 10, 11 } },
    { {  6 }, { 7,   9, 11, 12 } },
    { {  7 }, { 9,  12, 13 } },
    { {  8 }, { 10, 11, 12, 13 } },
    { {  9 }, { 10, 11, 12, 13 } },
    { { 10 }, { 11, 13 } },
    { { 11 }, { 12, 13 } },
    { { 12 }, { 13 } },
    // the last star of edsup is always empty, since only edges whose point ids
    // p < q are stored, however, the indices exist in edsup2 to allow client
    // code to be simpler and consistent with using the other derived
    // data structures
    { { 13 }, { } }
  };

  // find out number of points from mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  // this is more of a test on this test
  ensure_equals( "number of points (star centers) in edsup incorrect",
                 npoin, correct_edsup.size() );

  // Test generated derived data structure, edges surrounding points
  for (std::size_t p=0; p<npoin; ++p) {
    // extract edge end-point ids from generated edges surrounding points
    std::vector< std::size_t > edge;
    for (auto i=edsup2[p]+1; i<=edsup2[p+1]; ++i) edge.push_back( edsup1[i] );
    // find correct star-center point id for list of star-end points
    auto it = correct_edsup.find( p );
    // test if star-end point ids exist emanating from star-center p
    ensure( "star-center point id '" + std::to_string(p) + "' generated into "
            "edsup but not in correct edsup",
            it != correct_edsup.end() );
    // test if star-end point ids starting from star-center point p are correct
    if (it != correct_edsup.end()) {
      ensure_equals( "number of star-end points starting from star-center " +
                     std::to_string(p) + " in generated edsup incorrect",
                     it->second.size(), edge.size() );
      ensure( "star-end point ids starting from star-center'" +
              std::to_string(p) + "' incorrect",
              edge == it->second );
    }
  }
}

//! Generate and test edges surrounding points for triangle-only mesh
template<> template<>
void DerivedData_object::test< 34 >() {
  set_test_name( "genEdsup for triangles" );

  // Mesh connectivity for simple triangle-only mesh
  std::vector< std::size_t > inpoel { 1,  9,  2,
                                      1,  4,  9,
                                      2,  9,  3,
                                      3,  9,  4,
                                      5,  6, 10,
                                      5, 10,  8,
                                      6,  7, 10,
                                      7,  8, 10,
                                      1,  2, 11,
                                      1, 11,  5,
                                      2,  6, 11,
                                      5, 11,  6,
                                      2,  3, 12,
                                      2, 12,  6,
                                      3,  7, 12,
                                      6, 12,  7,
                                      3,  4, 13,
                                      3, 13,  7,
                                      4,  8, 13,
                                      7, 13,  8,
                                      1, 14,  4,
                                      1,  5, 14,
                                      4, 14,  8,
                                      5,  8, 14 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Generate edges surrounding points
  auto edsup = tk::genEdsup( inpoel, 3, tk::genEsup(inpoel,3) );

  auto& edsup1 = edsup.first;
  auto& edsup2 = edsup.second;

  // Generate correct solution for edges surrounding points
  std::map< std::size_t, std::vector< std::size_t > > correct_edsup {
    { { 0 }, { 1, 3, 4, 8, 10, 13 } },
    { { 1 }, { 2, 5, 8, 10, 11 } },
    { { 2 }, { 3, 6, 8, 11, 12 } },
    { { 3 }, { 7, 8, 12, 13 } },
    { { 4 }, { 5, 7, 9, 10, 13 } },
    { { 5 }, { 6, 9, 10, 11 } },
    { { 6 }, { 7, 9, 11, 12 } },
    { { 7 }, { 9, 12, 13 } },
    // the last star of edsup is always empty, and, depending on the mesh,
    // sometimes not only the last one but the last few, since only edges whose
    // point ids p < q are stored, however, the indices exist for all points in
    // edsup2 to allow client code to be simpler and consistent with using the
    // other derived data structures
    { { 8 }, { } },
    { { 9 }, { } },
    { { 10 }, { } },
    { { 11 }, { } },
    { { 12 }, { } },
    { { 13 }, { } }
  };

  // find out number of points from mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  // this is more of a test on this test
  ensure_equals( "number of points (star centers) in edsup incorrect",
                 npoin, correct_edsup.size() );

  // Test generated derived data structure, edges surrounding points
  for (std::size_t p=0; p<npoin; ++p) {
    // extract edge end-point ids from generated edges surrounding points
    std::vector< std::size_t > edge;
    for (auto i=edsup2[p]+1; i<=edsup2[p+1]; ++i) edge.push_back( edsup1[i] );
    // find correct star-center point id for list of star-end points
    auto it = correct_edsup.find( p );
    // test if star-end point ids exist emanating from star-center p
    ensure( "star-center point id '" + std::to_string(p) + "' generated into "
            "edsup but not in correct edsup",
            it != correct_edsup.end() );
    // test if star-end point ids starting from star-center point p are correct
    if (it != correct_edsup.end()) {
      ensure_equals( "number of star-end points starting from star-center " +
                     std::to_string(p) + " in generated edsup incorrect",
                     it->second.size(), edge.size() );
      ensure( "star-end point ids starting from star-center'" +
              std::to_string(p) + "' incorrect",
              edge == it->second );
    }
  }
}

//! Attempt to generate edge connectivity with empty connectivity
template<> template<>
void DerivedData_object::test< 35 >() {
  set_test_name( "genInpoed throws with empty inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    std::vector< std::size_t > empty;
    tk::genInpoed( empty, 4, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genInpoed if it throws on non-positive number of nodes per elements
template<> template<>
void DerivedData_object::test< 36 >() {
  set_test_name( "genInpoed throws on non-positive nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genInpoed( inpoel, 0, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! \brief Test genInpoed if it throws on non-tet or non-tri number of nodes per
//!   elements
template<> template<>
void DerivedData_object::test< 37 >() {
  set_test_name( "genInpoed throws on unsupported nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genInpoed( inpoel, 5, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genInpoed if it throws with empty element surrounding points . first
template<> template<>
void DerivedData_object::test< 38 >() {
  set_test_name( "genInpoed throws with empty esup.first" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genInpoed( inpoel, 4, {} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! Test genInpoed if it throws with empty element surrounding points .second
template<> template<>
void DerivedData_object::test< 39 >() {
  set_test_name( "genInpoed throws with empty esup.second" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genInpoed( inpoel, 4, {{1},{}} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! \brief Test genInpoed if it throws on inpoel non-divisible by the number of
//!   nodes per elements
template<> template<>
void DerivedData_object::test< 40 >() {
  set_test_name( "genInpoed throws on inpoel non-div nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield invalid read" );
  #else
  try {
    // Partial mesh mesh connectivity
    std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                        14,  4, 13 };
    tk::genInpoed( inpoel, 4, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Generate and test edge connectivity for tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 41 >() {
  set_test_name( "genInpoed for tetrahedra" );

  // mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12, 14,  9, 11,
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
    12, 13
  };

  // this is more of a test on this test
  ensure_equals(
    "number of edges in correct edge connectivity non-divisble by 2",
    correct_inpoed.size() % 2, 0 );

  // test if edge connectivity is the correct size
  ensure_equals(
    "number of edges in generated edge connectivity non-divisble by 2",
    inpoed.size() % 2, 0 );
  ensure_equals( "number of edges in edge connectivity incorrect",
                 inpoed.size(), correct_inpoed.size() );

  // this if edge connectivity correct
  ensure( "generatd edge connectivity incorrect", inpoed == correct_inpoed );
}

//! Generate and test edge connectivity for triangle-only mesh
template<> template<>
void DerivedData_object::test< 42 >() {
  set_test_name( "genInpoed for triangles" );

  // Mesh connectivity for simple triangle-only mesh
  std::vector< std::size_t > inpoel { 1,  9,  2,
                                      1,  4,  9,
                                      2,  9,  3,
                                      3,  9,  4,
                                      5,  6, 10,
                                      5, 10,  8,
                                      6,  7, 10,
                                      7,  8, 10,
                                      1,  2, 11,
                                      1, 11,  5,
                                      2,  6, 11,
                                      5, 11,  6,
                                      2,  3, 12,
                                      2, 12,  6,
                                      3,  7, 12,
                                      6, 12,  7,
                                      3,  4, 13,
                                      3, 13,  7,
                                      4,  8, 13,
                                      7, 13,  8,
                                      1, 14,  4,
                                      1,  5, 14,
                                      4, 14,  8,
                                      5,  8, 14 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Generate edge connectivity
  auto inpoed = tk::genInpoed( inpoel, 3, tk::genEsup(inpoel,3) );

  // Generate correct solution for edge connectivity
  std::vector< std::size_t > correct_inpoed {
    0, 1, 0, 3, 0, 4, 0, 8, 0, 10, 0, 13,
    1, 2, 1, 5, 1, 8, 1, 10, 1, 11,
    2, 3, 2, 6, 2, 8, 2, 11, 2, 12,
    3, 7, 3, 8, 3, 12, 3, 13,
    4, 5, 4, 7, 4, 9, 4, 10, 4, 13,
    5, 6, 5, 9, 5, 10, 5, 11,
    6, 7, 6, 9, 6, 11, 6, 12,
    7, 9, 7, 12, 7, 13
  };

  // this is more of a test on this test
  ensure_equals(
    "number of edges in correct edge connectivity non-divisble by 2",
    correct_inpoed.size() % 2, 0 );

  // test if edge connectivity is the correct size
  ensure_equals(
    "number of edges in generated edge connectivity non-divisble by 2",
    inpoed.size() % 2, 0 );
  ensure_equals( "number of edges in edge connectivity incorrect",
                 inpoed.size(), correct_inpoed.size() );

  // this if edge connectivity correct
  ensure( "generatd edge connectivity incorrect", inpoed == correct_inpoed );
}

//! Attempt to generate edges of elements with empty connectivity
template<> template<>
void DerivedData_object::test< 43 >() {
  set_test_name( "genInedel throws with empty inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    std::vector< std::size_t > empty;
    tk::genInedel( empty, 4, tk::genInpoed(inpoel,4,tk::genEsup(inpoel,4)) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genInedel if it throws on non-positive number of nodes per elements
template<> template<>
void DerivedData_object::test< 44 >() {
  set_test_name( "genInpoed throws on non-positive nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genInedel( inpoel, 0, tk::genInpoed(inpoel,4,tk::genEsup(inpoel,4)) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! \brief Test genInedel if it throws on non-tet or non-tri number of nodes per
//!   elements
template<> template<>
void DerivedData_object::test< 45 >() {
  set_test_name( "genInpoed throws on unsupported nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genInedel( inpoel, 5, tk::genInpoed(inpoel,4,tk::genEsup(inpoel,4)) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genInedel if it throws with empty edge connectivity
template<> template<>
void DerivedData_object::test< 46 >() {
  set_test_name( "genInpoed throws with empty inpoed" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genInedel( inpoel, 0, {} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! \brief Test genInedel if it throws on inpoel non-divisible by the number of
//!   nodes per elements
template<> template<>
void DerivedData_object::test< 47 >() {
  set_test_name( "genInedel throws on inpoel non-div nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield invalid read" );
  #else
  try {
    // Partial mesh mesh connectivity
    std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                        14,  4, 13 };
    tk::genInedel( inpoel, 4, tk::genInpoed(inpoel,4,tk::genEsup(inpoel,4)) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Generate and test edges of elements for tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 48 >() {
  set_test_name( "genInedel for tetrahedra" );

  // mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12, 14,  9, 11,
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

  // Generate edges of elements
  auto inedel =
    tk::genInedel( inpoel, 4, tk::genInpoed(inpoel,4,tk::genEsup(inpoel,4)) );

  // Generate correct solution for edges of elements
  std::map< std::size_t, std::vector< std::size_t > > correct_inedel {
    { { 0 }, { 47, 36, 37, 39, 44, 45 } },
    { { 1 }, { 41, 42, 43, 48, 46, 47 } },
    { { 2 }, { 48, 46, 47, 37, 38, 39 } },
    { { 3 }, { 40, 41, 43, 47, 44, 45 } },
    { { 4 }, { 2, 4, 5, 23, 24, 45 } },
    { { 5 }, { 30, 31, 25, 26, 28, 41 } },
    { { 6 }, { 33, 35, 21, 22, 24, 43 } },
    { { 7 }, { 33, 34, 29, 30, 32, 42 } },
    { { 8 }, { 31, 32, 12, 14, 15, 46 } },
    { { 9 }, { 1, 3, 5, 17, 19, 39 } },
    { { 10 }, { 17, 18, 11, 13, 15, 38 } },
    { { 11 }, { 13, 14, 6, 8, 10, 37 } },
    { { 12 }, { 16, 18, 19, 34, 35, 48 } },
    { { 13 }, { 26, 27, 20, 22, 23, 40 } },
    { { 14 }, { 0, 3, 4, 8, 9, 36 } },
    { { 15 }, { 7, 9, 10, 27, 28, 44 } },
    { { 16 }, { 26, 27, 28, 40, 41, 44 } },
    { { 17 }, { 8, 9, 10, 36, 37, 44 } },
    { { 18 }, { 22, 23, 24, 40, 43, 45 } },
    { { 19 }, { 33, 34, 35, 42, 43, 48 } },
    { { 20 }, { 13, 14, 15, 46, 37, 38 } },
    { { 21 }, { 30, 31, 32, 41, 42, 46 } },
    { { 22 }, { 17, 18, 19, 48, 38, 39 } },
    { { 23 }, { 3, 4, 5, 36, 39, 45 } }
  };

  // find out number of elements from mesh connectivity
  auto nelem = inpoel.size()/4;

  // this is more of a test on this test
  ensure_equals( "number of elements in inedel incorrect",
                 nelem, correct_inedel.size() );

  // test generated derived data structure, edges of elements
  for (std::size_t e=0; e<nelem; ++e) {
    // extract element ids from generated edges of elements
    std::vector< std::size_t > edges;
    for (std::size_t i=0; i<6; ++i) edges.push_back( inedel[e*6+i] );
    // find correct edge ids of element e
    auto it = correct_inedel.find( e );
    // test if edge ids exist for element e
    ensure( "element id '" + std::to_string(e) + "' referred by inedel but not "
            "in correct inedel",
            it != correct_inedel.end() );
    // test if edge ids of element e are correct
    if (it != correct_inedel.end()) {
      ensure_equals( "number of edges of element " + std::to_string(e) +
                     " in generated inedel incorrect",
                     it->second.size(), 6 );
      ensure( "edge ids surrounding element '" + std::to_string(e) +
              "' incorrect", edges == it->second );
    }
  }
}

//! Generate and test edges of elements for triangle-only mesh
template<> template<>
void DerivedData_object::test< 49 >() {
  set_test_name( "genInedel for triangles" );

  // Mesh connectivity for simple triangle-only mesh
  std::vector< std::size_t > inpoel { 1,  9,  2,
                                      1,  4,  9,
                                      2,  9,  3,
                                      3,  9,  4,
                                      5,  6, 10,
                                      5, 10,  8,
                                      6,  7, 10,
                                      7,  8, 10,
                                      1,  2, 11,
                                      1, 11,  5,
                                      2,  6, 11,
                                      5, 11,  6,
                                      2,  3, 12,
                                      2, 12,  6,
                                      3,  7, 12,
                                      6, 12,  7,
                                      3,  4, 13,
                                      3, 13,  7,
                                      4,  8, 13,
                                      7, 13,  8,
                                      1, 14,  4,
                                      1,  5, 14,
                                      4, 14,  8,
                                      5,  8, 14 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  auto inedel =
    tk::genInedel( inpoel, 3, tk::genInpoed(inpoel,3,tk::genEsup(inpoel,3)) );

  // Generate correct solution for edges of elements
  std::map< std::size_t, std::vector< std::size_t > > correct_inedel {
    { { 0 }, { 0, 3, 8 } },
    { { 1 }, { 1, 3, 17 } },
    { { 2 }, { 6, 8, 13 } },
    { { 3 }, { 11, 13, 17 } },
    { { 4 }, { 20, 22, 26 } },
    { { 5 }, { 21, 22, 33 } },
    { { 6 }, { 25, 26, 30 } },
    { { 7 }, { 29, 30, 33 } },
    { { 8 }, { 0, 4, 9 } },
    { { 9 }, { 2, 4, 23 } },
    { { 10 }, { 7, 9, 27 } },
    { { 11 }, { 20, 23, 27 } },
    { { 12 }, { 6, 10, 14 } },
    { { 13 }, { 7, 10, 28 } },
    { { 14 }, { 12, 14, 31 } },
    { { 15 }, { 25, 28, 31 } },
    { { 16 }, { 11, 15, 18 } },
    { { 17 }, { 12, 15, 32 } },
    { { 18 }, { 16, 18, 34 } },
    { { 19 }, { 29, 32, 34 } },
    { { 20 }, { 1, 5, 19 } },
    { { 21 }, { 2, 5, 24 } },
    { { 22 }, { 16, 19, 35 } },
    { { 23 }, { 21, 24, 35 } }
  };

  // find out number of elements from mesh connectivity
  auto nelem = inpoel.size()/3;

  // this is more of a test on this test
  ensure_equals( "number of elements in inedel incorrect",
                 nelem, correct_inedel.size() );

  // test generated derived data structure, edges of elements
  for (std::size_t e=0; e<nelem; ++e) {
    // extract element ids from generated edges of elements
    std::vector< std::size_t > edges;
    for (std::size_t i=0; i<3; ++i) edges.push_back( inedel[e*3+i] );
    // find correct edge ids of element e
    auto it = correct_inedel.find( e );
    // test if edge ids exist for element e
    ensure( "element id '" + std::to_string(e) + "' referred by inedel but not "
            "in correct inedel",
            it != correct_inedel.end() );
    // test if edge ids of element e are correct
    if (it != correct_inedel.end()) {
      ensure_equals( "number of edges of element " + std::to_string(e) +
                     " in generated inedel incorrect",
                     it->second.size(), 3 );
      ensure( "edge ids surrounding element '" + std::to_string(e) +
              "' incorrect", edges == it->second );
    }
  }
}

//! Attempt to generate elements surrounding edges with empty connectivity
template<> template<>
void DerivedData_object::test< 50 >() {
  set_test_name( "genEsued throws with empty inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    std::vector< std::size_t > empty;
    tk::genEsued( empty, 4, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genEsued if it throws on non-positive number of nodes per elements
template<> template<>
void DerivedData_object::test< 51 >() {
  set_test_name( "genEsued throws on non-positive nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEsued( inpoel, 0, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! \brief Test genEsued if it throws on non-tet or non-tri number of nodes per
//!   elements
template<> template<>
void DerivedData_object::test< 52 >() {
  set_test_name( "genEsued throws on unsupported nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield floating point exception" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEsued( inpoel, 5, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test genEsued if it throws with empty element surrounding points . first
template<> template<>
void DerivedData_object::test< 53 >() {
  set_test_name( "genEsued throws with empty esup.first" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEsued( inpoel, 4, {} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! Test genEsued if it throws with empty element surrounding points . second
template<> template<>
void DerivedData_object::test< 54 >() {
  set_test_name( "genEsued throws with empty esup.second" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > inpoel { 0, 1, 2, 3 };
    tk::genEsued( inpoel, 4, {{1},{}} );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown, test ok
  }
  #endif
}

//! \brief Test genEsued if it throws on inpoel non-divisible by the number of
//!   nodes per elements
template<> template<>
void DerivedData_object::test< 55 >() {
  set_test_name( "genEsued throws on inpoel non-div nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield invalid read" );
  #else
  try {
    // Partial mesh mesh connectivity
    std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                        14,  4, 13 };
    tk::genEsued( inpoel, 4, tk::genEsup(inpoel,4) );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Generate and test elements surrounding edges for tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 56 >() {
  set_test_name( "genEsued for tetrahedra" );

  // mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12, 14,  9, 11,
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

  // Generate edges surrounding points
  auto esup = tk::genEsup( inpoel, 4 );
  auto esued = tk::genEsued( inpoel, 4, esup );

  auto& esued1 = esued.first;
  auto& esued2 = esued.second;

  // Generate correct solution for elements surrounding edges
  std::set< std::vector< std::size_t > > correct_esued {
    { 4, 9, 23 },
    { 4 },
    { 4, 14, 23 },
    { 9 },
    { 9, 14, 23 },
    { 14 },
    { 11 },
    { 11, 15, 17 },
    { 11, 14, 17 },
    { 14, 15, 17 },
    { 15 },
    { 8 },
    { 8, 10, 20 },
    { 8, 11, 20 },
    { 10 },
    { 10, 11, 20 },
    { 9, 12, 22 },
    { 9, 10, 22 },
    { 10, 12, 22 },
    { 12 },
    { 4, 6, 18 },
    { 4, 13, 18 },
    { 6 },
    { 6, 13, 18 },
    { 13 },
    { 5 },
    { 5, 13, 16 },
    { 5, 15, 16 },
    { 13, 15, 16 },
    { 5, 7, 21 },
    { 5, 8, 21 },
    { 7 },
    { 7, 8, 21 },
    { 6, 12, 19 },
    { 6, 7, 19 },
    { 7, 12, 19 },
    { 0, 2, 11, 17, 20 },
    { 0, 2, 9, 22, 23 },
    { 0, 14, 17, 23 },
    { 2, 10, 20, 22 },
    { 1, 3, 6, 18, 19 },
    { 1, 7, 19, 21 },
    { 1, 3, 5, 16, 21 },
    { 3, 13, 16, 18 },
    { 0, 3, 15, 16, 17 },
    { 0, 3, 4, 18, 23 },
    { 0, 1, 2, 3 },
    { 1, 2, 8, 20, 21 },
    { 1, 2, 12, 19, 22 }
  };

  // find out number of edges in mesh
  auto nedge = tk::genInpoed(inpoel,4,esup).size()/2;

  // this is more of a test on this test
  ensure_equals( "number of edges in esued incorrect",
                 nedge, correct_esued.size() );

  // Test generated derived data structure, elements surrounding edges
  for (std::size_t e=0; e<nedge; ++e) {
    // extract element list generated for edge e
    std::vector< std::size_t > elements;
    for (auto i=esued2[e]+1; i<=esued2[e+1]; ++i) elements.push_back(esued1[i]);
    // store element list as string to output it in case test fails
    std::stringstream ss;
    for (auto i : elements) ss << i << " ";
    // attempt to find element list in correct esued
    auto it = correct_esued.find( elements );
    // test if element list can be found among the correct ones
    ensure( "element list { " + ss.str() + "} surrounding edge '" +
            std::to_string(e) + "' generated into esued but not in correct "
            "esued",
            it != correct_esued.end() );
    // remove element list just tested from correct esued, this ensures that the
    // generated esued does not contain edges whose element lists would be
    // exactly the same, as later tests in this loop would fail in that case
    correct_esued.erase( elements );
  }
}

//! Generate and test elements surrounding edges for triangle-only mesh
template<> template<>
void DerivedData_object::test< 57 >() {
  set_test_name( "genEsued for triangles" );

  // Mesh connectivity for simple triangle-only mesh
  std::vector< std::size_t > inpoel { 1,  9,  2,
                                      1,  4,  9,
                                      2,  9,  3,
                                      3,  9,  4,
                                      5,  6, 10,
                                      5, 10,  8,
                                      6,  7, 10,
                                      7,  8, 10,
                                      1,  2, 11,
                                      1, 11,  5,
                                      2,  6, 11,
                                      5, 11,  6,
                                      2,  3, 12,
                                      2, 12,  6,
                                      3,  7, 12,
                                      6, 12,  7,
                                      3,  4, 13,
                                      3, 13,  7,
                                      4,  8, 13,
                                      7, 13,  8,
                                      1, 14,  4,
                                      1,  5, 14,
                                      4, 14,  8,
                                      5,  8, 14 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Generate edges surrounding points
  auto esup = tk::genEsup( inpoel, 3 );
  auto esued = tk::genEsued( inpoel, 3, esup );

  auto& esued1 = esued.first;
  auto& esued2 = esued.second;

  // Generate correct solution for elements surrounding edges
  std::set< std::vector< std::size_t > > correct_esued {
    { 0, 1 },
    { 0, 8 },
    { 1, 20 },
    { 8, 9 },
    { 9, 21 },
    { 20, 21 },
    { 0, 2 },
    { 2, 12 },
    { 8, 10 },
    { 10, 13 },
    { 12, 13 },
    { 2, 3 },
    { 3, 16 },
    { 12, 14 },
    { 14, 17 },
    { 16, 17 },
    { 1, 3 },
    { 16, 18 },
    { 18, 22 },
    { 20, 22 },
    { 4, 11 },
    { 4, 5 },
    { 5, 23 },
    { 9, 11 },
    { 21, 23 },
    { 4, 6 },
    { 6, 15 },
    { 10, 11 },
    { 13, 15 },
    { 6, 7 },
    { 7, 19 },
    { 14, 15 },
    { 17, 19 },
    { 5, 7 },
    { 18, 19 },
    { 22, 23 }
  };

  // find out number of edges in mesh
  auto nedge = tk::genInpoed(inpoel,3,esup).size()/2;

  // this is more of a test on this test
  ensure_equals( "number of edges in esued incorrect",
                 nedge, correct_esued.size() );

  // Test generated derived data structure, elements surrounding edges
  for (std::size_t e=0; e<nedge; ++e) {
    // extract element list generated for edge e
    std::vector< std::size_t > elements;
    for (auto i=esued2[e]+1; i<=esued2[e+1]; ++i) elements.push_back(esued1[i]);
    // store element list as string to output it in case test fails
    std::stringstream ss;
    for (auto i : elements) ss << i << " ";
    // attempt to find element list in correct esued
    auto it = correct_esued.find( elements );
    // test if element list can be found among the correct ones
    ensure( "element list { " + ss.str() + "} surrounding edge '" +
            std::to_string(e) + "' generated into esued but not in correct "
            "esued",
            it != correct_esued.end() );
    // remove element list just tested from correct esued, this ensures that the
    // generated esued does not contain edges whose element lists would be
    // exactly the same, as later tests in this loop would fail in that case
    correct_esued.erase( elements );
  }
}

//! Generate and test number of boundary-faces for tet-only mesh
template<> template<>
void DerivedData_object::test< 65 >() {
  set_test_name( "Boundary faces (genNbfacTet) for tetrahedra" );

  std::map< int, std::vector< std::size_t > > bface, t_bface;

  // Create unstructured-mesh object to read into
  tk::UnsMesh inmesh;
  // Read in mesh from file
  std::string infile( REGRESSION_DIR
                      "/meshconv/gmsh_output/box_24_ss1.exo" );
  tk::ExodusIIMeshReader er( infile );
  auto tnbfac = er.readSidesetFaces( t_bface );

  // Test if the number of boundary faces is correct
  ensure_equals( "total number of boundary faces incorrect",
                 tnbfac, 24 );

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                      10, 14, 13, 12,
                                      14, 13, 12,  9,
                                      10, 14, 12, 11,
                                       1, 14,  5, 11,
                                       7,  6, 10, 12,
                                      14,  8,  5, 10,
                                       8,  7, 10, 13,
                                       7, 13,  3, 12,
                                       1,  4, 14,  9,
                                      13,  4,  3,  9,
                                       3,  2, 12,  9,
                                       4,  8, 14, 13,
                                       6,  5, 10, 11,
                                       1,  2,  9, 11,
                                       2,  6, 12, 11,
                                       6, 10, 12, 11,
                                       2, 12,  9, 11,
                                       5, 14, 10, 11,
                                      14,  8, 10, 13,
                                      13,  3, 12,  9,
                                       7, 10, 13, 12,
                                      14,  4, 13,  9,
                                      14,  1,  9, 11 };

  // Boundary-face node connectivity for the entire mesh
  std::vector< std::size_t > t_triinpoel {  2,  9,  3,
                                            1,  9,  2,
                                            3,  9,  4,
                                            1,  4,  9,
                                            7, 10,  6,
                                            6, 10,  5,
                                            8, 10,  7,
                                           10,  8,  5,
                                            6, 11,  2,
                                            2, 11,  1,
                                           11,  6,  5,
                                           11,  5,  1,
                                            3, 12,  2,
                                           12,  6,  2,
                                            7, 12,  3,
                                           12,  7,  6,
                                           13,  7,  3,
                                            4, 13,  3,
                                           13,  8,  7,
                                            8, 13,  4,
                                            5,  8, 14,
                                            1,  5, 14,
                                            4, 14,  8,
                                            1, 14,  4 };

  // Correct Boundary-face node connectivity
  std::vector< std::size_t > correct_triinpoel {  7, 10,  6,
                                                  6, 10,  5,
                                                  8, 10,  7,
                                                 10,  8,  5,
                                                  5,  8, 14,
                                                  1,  5, 14,
                                                  4, 14,  8,
                                                  1, 14,  4,
                                                  2,  9,  3,
                                                  1,  9,  2,
                                                  3,  9,  4,
                                                  1,  4,  9,
                                                  3, 12,  2,
                                                 12,  6,  2,
                                                  7, 12,  3,
                                                 12,  7,  6,
                                                 13,  7,  3,
                                                  4, 13,  3,
                                                 13,  8,  7,
                                                  8, 13,  4,
                                                  6, 11,  2,
                                                  2, 11,  1,
                                                 11,  6,  5,
                                                 11,  5,  1 };

  std::unordered_map< std::size_t, std::size_t > lid {
          { {0}, {0} },
          { {1}, {1} },
          { {2}, {2} },
          { {3}, {3} },
          { {4}, {4} },
          { {5}, {5} },
          { {6}, {6} },
          { {7}, {7} },
          { {8}, {8} },
          { {9}, {9} },
          { {10}, {10} },
          { {11}, {11} },
          { {12}, {12} },
          { {13}, {13} } };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  tk::shiftToZero( t_triinpoel );
  tk::shiftToZero( correct_triinpoel );

  std::vector< std::size_t > triinpoel;

  auto nbfac = tk::genNbfacTet( tnbfac, inpoel, t_triinpoel, t_bface, lid,
                                triinpoel, bface );

  ensure_equals( "number of boundary-faces is incorrect",
                 nbfac, tnbfac );

  ensure_equals( "total number of entries in triinpoel is incorrect",
                 triinpoel.size(), correct_triinpoel.size() );

  for(std::size_t i=0 ; i<triinpoel.size(); ++i)
  {
    ensure_equals("incorrect entry " + std::to_string(i) + " in triinpoel",
                    triinpoel[i], correct_triinpoel[i]);
  }
}

//! Generate and test face-data structures for tet-only mesh
template<> template<>
void DerivedData_object::test< 66 >() {
  set_test_name( "genEsuf for tetrahedra" );

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12,  9,  1, 24,
                                       1, 12, 24, 20,
                                      12, 24, 20, 23,
                                      24, 20, 23, 13,
                                      20, 23, 13,  5,
                                      23, 13,  5, 16,
                                      23, 24, 13, 27,
                                      12, 24, 23, 28,
                                      12,  9, 24, 28,
                                      23, 13, 16, 27,
                                      24, 23, 28, 27,
                                      23, 12, 28, 19,
                                       9, 24, 28, 29,
                                      24, 13, 27, 17,
                                      16, 23, 27, 30,
                                      12,  9, 28, 21,
                                      13, 16, 27, 22,
                                       9, 24, 29, 17,
                                      16, 23, 30, 19,
                                      23, 28, 27, 30,
                                      28, 24, 27, 29,
                                      28, 23, 19, 30,
                                      27, 24, 17, 29,
                                      12, 28, 19, 21,
                                      13, 27, 17,  6,
                                       9, 28, 21, 29,
                                      27, 16, 30, 22,
                                      13, 27,  6, 14,
                                      19, 12, 21,  4,
                                      21,  9, 29, 10,
                                      29,  9, 17,  2,
                                      29,  9,  2, 10,
                                      30, 16, 19,  8,
                                      30, 16,  8, 15,
                                      27, 13, 22, 14,
                                      16, 30, 22, 15,
                                      28, 27, 30, 31,
                                      27, 28, 29, 31,
                                      28, 19, 21, 25,
                                      17, 27, 29, 26,
                                      28, 21, 29, 31,
                                      19, 28, 30, 25,
                                      27, 17,  6, 14,
                                      21, 19,  4, 25,
                                      19, 30,  8, 15,
                                      17, 29,  2, 10,
                                      30, 27, 22, 31,
                                      21, 28, 25, 31,
                                      28, 30, 25, 31,
                                      27, 17, 14, 26,
                                       4, 21, 25, 11,
                                      29, 27, 31, 26,
                                      30, 19, 25, 15,
                                      29, 21, 10, 31,
                                      17, 29, 10, 26,
                                      22, 27, 14, 31,
                                      30, 22, 15, 31,
                                      27, 14, 31, 26,
                                      21, 25, 11, 31,
                                      25, 30, 15, 31,
                                      21, 10, 31, 11,
                                      10, 29, 31, 26,
                                      22, 15, 31,  7,
                                      14, 22, 31,  7,
                                      15, 25, 31, 18,
                                      25, 11, 31, 18,
                                      10, 31, 11,  3,
                                      31, 10, 26, 18,
                                      14, 31, 26,  7,
                                      15, 31,  7, 18,
                                      11, 31, 18,  3,
                                      31, 10, 18,  3,
                                      31, 26,  7, 18 };

  // boundary elements for this mesh
  std::vector< std::size_t > belem {
    { 1,
      1,
      2,
      2,
      3,
      4,
      5,
      5,
      6,
      6,
     12,
     14,
     16,
     17,
     18,
     19,
     25,
     28,
     29,
     29,
     30,
     31,
     32,
     33,
     34,
     35,
     36,
     43,
     44,
     45,
     46,
     50,
     51,
     51,
     53,
     55,
     61,
     63,
     64,
     65,
     66,
     67,
     68,
     69,
     70,
     71,
     72,
     73 } };

  // Shift elem IDs to start from zero
  tk::shiftToZero( belem );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Generate elements surrounding elements
  auto esuelTet = tk::genEsuelTet( inpoel,
                                   tk::genEsup(inpoel, 4) );

  // Generate number of total and boundary faces
  std::size_t nbfac(48);
  auto ntfac = tk::genNtfac(4, nbfac, esuelTet);

  // Generate esuf
  auto esuf = tk::genEsuf(4, ntfac, nbfac, belem, esuelTet);

  // this is more of a test on this test
  ensure_equals( "total number of elems is incorrect",
                 esuelTet.size()/4, 73 );
  ensure_equals( "total number of faces is incorrect",
                 ntfac, 170 );
  ensure_equals( "number of boundary faces is incorrect",
                 nbfac, 48 );

  // Generate correct solution for elements surrounding faces
  std::vector< int > correct_esuf { 1,  0,
                                    1,  0,
                                    2,  0,
                                    2,  0,
                                    3,  0,
                                    4,  0,
                                    5,  0,
                                    5,  0,
                                    6,  0,
                                    6,  0,
                                   12,  0,
                                   14,  0,
                                   16,  0,
                                   17,  0,
                                   18,  0,
                                   19,  0,
                                   25,  0,
                                   28,  0,
                                   29,  0,
                                   29,  0,
                                   30,  0,
                                   31,  0,
                                   32,  0,
                                   33,  0,
                                   34,  0,
                                   35,  0,
                                   36,  0,
                                   43,  0,
                                   44,  0,
                                   45,  0,
                                   46,  0,
                                   50,  0,
                                   51,  0,
                                   51,  0,
                                   53,  0,
                                   55,  0,
                                   61,  0,
                                   63,  0,
                                   64,  0,
                                   65,  0,
                                   66,  0,
                                   67,  0,
                                   68,  0,
                                   69,  0,
                                   70,  0,
                                   71,  0,
                                   72,  0,
                                   73,  0,
                                    1,  2,
                                    1,  9,
                                    2,  3,
                                    3,  4,
                                    3,  8,
                                    4,  5,
                                    4,  7,
                                    5,  6,
                                    6, 10,
                                    7, 14,
                                    7, 10,
                                    7, 11,
                                    8, 11,
                                    8, 12,
                                    8,  9,
                                    9, 13,
                                    9, 16,
                                   10, 17,
                                   10, 15,
                                   11, 20,
                                   11, 21,
                                   12, 24,
                                   12, 22,
                                   13, 21,
                                   13, 26,
                                   13, 18,
                                   14, 25,
                                   14, 23,
                                   15, 20,
                                   15, 27,
                                   15, 19,
                                   16, 26,
                                   16, 24,
                                   17, 27,
                                   17, 35,
                                   18, 23,
                                   18, 31,
                                   19, 22,
                                   19, 33,
                                   20, 37,
                                   20, 22,
                                   21, 23,
                                   21, 38,
                                   22, 42,
                                   23, 40,
                                   24, 39,
                                   24, 29,
                                   25, 43,
                                   25, 28,
                                   26, 41,
                                   26, 30,
                                   27, 36,
                                   27, 47,
                                   28, 43,
                                   28, 35,
                                   29, 44,
                                   30, 32,
                                   30, 54,
                                   31, 46,
                                   31, 32,
                                   32, 46,
                                   33, 45,
                                   33, 34,
                                   34, 45,
                                   34, 36,
                                   35, 56,
                                   36, 57,
                                   37, 47,
                                   37, 49,
                                   37, 38,
                                   38, 41,
                                   38, 52,
                                   39, 44,
                                   39, 48,
                                   39, 42,
                                   40, 52,
                                   40, 55,
                                   40, 50,
                                   41, 54,
                                   41, 48,
                                   42, 49,
                                   42, 53,
                                   43, 50,
                                   44, 51,
                                   45, 53,
                                   46, 55,
                                   47, 56,
                                   47, 57,
                                   48, 49,
                                   48, 59,
                                   49, 60,
                                   50, 58,
                                   51, 59,
                                   52, 58,
                                   52, 62,
                                   53, 60,
                                   54, 61,
                                   54, 62,
                                   55, 62,
                                   56, 58,
                                   56, 64,
                                   57, 63,
                                   57, 60,
                                   58, 69,
                                   59, 66,
                                   59, 61,
                                   60, 65,
                                   61, 67,
                                   62, 68,
                                   63, 70,
                                   63, 64,
                                   64, 69,
                                   65, 66,
                                   65, 70,
                                   66, 71,
                                   67, 71,
                                   67, 72,
                                   68, 73,
                                   68, 72,
                                   69, 73,
                                   70, 73,
                                   71, 72 };

  ensure_equals( "total number of entries in esuf is incorrect",
                 esuf.size(), correct_esuf.size() );

  for(std::size_t i=0 ; i<esuf.size(); ++i)
  {
          ensure_equals("incorrect entry " + std::to_string(i) + " in esuf",
                          esuf[i], correct_esuf[i]-1);
  }
}

//! Generate and test face-node connectivities for tet-only mesh
template<> template<>
void DerivedData_object::test< 67 >() {
  set_test_name( "genInpofa for tetrahedra" );

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12,  9,  1, 24,
                                       1, 12, 24, 20,
                                      12, 24, 20, 23,
                                      24, 20, 23, 13,
                                      20, 23, 13,  5,
                                      23, 13,  5, 16,
                                      23, 24, 13, 27,
                                      12, 24, 23, 28,
                                      12,  9, 24, 28,
                                      23, 13, 16, 27,
                                      24, 23, 28, 27,
                                      23, 12, 28, 19,
                                       9, 24, 28, 29,
                                      24, 13, 27, 17,
                                      16, 23, 27, 30,
                                      12,  9, 28, 21,
                                      13, 16, 27, 22,
                                       9, 24, 29, 17,
                                      16, 23, 30, 19,
                                      23, 28, 27, 30,
                                      28, 24, 27, 29,
                                      28, 23, 19, 30,
                                      27, 24, 17, 29,
                                      12, 28, 19, 21,
                                      13, 27, 17,  6,
                                       9, 28, 21, 29,
                                      27, 16, 30, 22,
                                      13, 27,  6, 14,
                                      19, 12, 21,  4,
                                      21,  9, 29, 10,
                                      29,  9, 17,  2,
                                      29,  9,  2, 10,
                                      30, 16, 19,  8,
                                      30, 16,  8, 15,
                                      27, 13, 22, 14,
                                      16, 30, 22, 15,
                                      28, 27, 30, 31,
                                      27, 28, 29, 31,
                                      28, 19, 21, 25,
                                      17, 27, 29, 26,
                                      28, 21, 29, 31,
                                      19, 28, 30, 25,
                                      27, 17,  6, 14,
                                      21, 19,  4, 25,
                                      19, 30,  8, 15,
                                      17, 29,  2, 10,
                                      30, 27, 22, 31,
                                      21, 28, 25, 31,
                                      28, 30, 25, 31,
                                      27, 17, 14, 26,
                                       4, 21, 25, 11,
                                      29, 27, 31, 26,
                                      30, 19, 25, 15,
                                      29, 21, 10, 31,
                                      17, 29, 10, 26,
                                      22, 27, 14, 31,
                                      30, 22, 15, 31,
                                      27, 14, 31, 26,
                                      21, 25, 11, 31,
                                      25, 30, 15, 31,
                                      21, 10, 31, 11,
                                      10, 29, 31, 26,
                                      22, 15, 31,  7,
                                      14, 22, 31,  7,
                                      15, 25, 31, 18,
                                      25, 11, 31, 18,
                                      10, 31, 11,  3,
                                      31, 10, 26, 18,
                                      14, 31, 26,  7,
                                      15, 31,  7, 18,
                                      11, 31, 18,  3,
                                      31, 10, 18,  3,
                                      31, 26,  7, 18 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Generate elements surrounding elements
  auto esuelTet = tk::genEsuelTet( inpoel,
                                    tk::genEsup(inpoel, 4) );

  // Generate number of total and boundary faces
  std::size_t nbfac(48);
  auto ntfac = tk::genNtfac(4, nbfac, esuelTet);

  std::map< int, std::vector< std::size_t > > bface {
          { { 0 }, {0,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                    31,
                    32,
                    33,
                    34,
                    35,
                    36,
                    37,
                    38,
                    39,
                    40,
                    41, 
                    42,
                    43,
                    44,
                    45,
                    46,
                    47 } } 
  };

  // Read boundary face-node connectivity
  std::vector< std::size_t > triinpoel { 24,  1,  9,
                                          9,  1, 12,
                                         20,  1, 24,
                                         12,  1, 20,
                                         23, 12, 20,
                                         20, 24, 13,
                                          5, 20, 13,
                                         23, 20,  5,
                                         16,  5, 13,
                                         16, 23,  5,
                                         12, 23, 19,
                                         13, 24, 17,
                                          9, 12, 21,
                                         16, 13, 22,
                                         24,  9, 17,
                                         23, 16, 19,
                                          6, 13, 17,
                                         14, 13,  6,
                                          4, 21, 12,
                                         12, 19,  4,
                                          9, 21, 10,
                                          2, 17,  9,
                                         10,  2,  9,
                                          8, 19, 16,
                                         15,  8, 16,
                                         14, 22, 13,
                                         15, 16, 22,
                                         14,  6, 17,
                                         25,  4, 19,
                                         15, 19,  8,
                                         10, 17,  2,
                                         26, 14, 17,
                                         11,  4, 25,
                                         21,  4, 11,
                                         15, 25, 19,
                                         26, 17, 10,
                                         10, 21, 11,
                                         15, 22,  7,
                                         22, 14,  7,
                                         25, 15, 18,
                                         11, 25, 18,
                                          3, 10, 11,
                                         18, 26, 10,
                                          7, 14, 26,
                                         18, 15,  7,
                                          3, 11, 18,
                                          3, 18, 10,
                                         18,  7, 26 };

  // Shift node IDs to start from zero
  tk::shiftToZero( triinpoel );

  // Generate inpofa
  auto inpofa = tk::genInpofaTet( ntfac, nbfac, inpoel, triinpoel, esuelTet );

  // Generate correct solution for elements surrounding faces
  std::vector< int > correct_inpofa {   9,  1, 24,
                                       12,  1,  9,
                                       24,  1, 20,
                                       20,  1, 12,
                                       20, 12, 23,
                                       13, 24, 20,
                                       13, 20,  5,
                                        5, 20, 23,
                                       13,  5, 16,
                                        5, 23, 16,
                                       19, 23, 12,
                                       17, 24, 13,
                                       21, 12,  9,
                                       22, 13, 16,
                                       17,  9, 24,
                                       19, 16, 23,
                                       17, 13,  6,
                                        6, 13, 14,
                                       12, 21,  4,
                                        4, 19, 12,
                                       10, 21,  9,
                                        9, 17,  2,
                                        9,  2, 10,
                                       16, 19,  8,
                                       16,  8, 15,
                                       13, 22, 14,
                                       22, 16, 15,
                                       17,  6, 14,
                                       19,  4, 25,
                                        8, 19, 15,
                                        2, 17, 10,
                                       17, 14, 26,
                                       25,  4, 11,
                                       11,  4, 21,
                                       19, 25, 15,
                                       10, 17, 26,
                                       11, 21, 10,
                                        7, 22, 15,
                                        7, 14, 22,
                                       18, 15, 25,
                                       18, 25, 11,
                                       11, 10,  3,
                                       10, 26, 18,
                                       26, 14,  7,
                                        7, 15, 18,
                                       18, 11,  3,
                                       10, 18,  3,
                                       26,  7, 18,
                                        1, 12, 24,
                                       24, 12,  9,
                                       12, 24, 20,
                                       24, 20, 23,
                                       23, 12, 24,
                                       20, 23, 13,
                                       23, 24, 13,
                                       23, 13,  5,
                                       16, 23, 13,
                                       24, 13, 27,
                                       13, 23, 27,
                                       27, 23, 24,
                                       24, 23, 28,
                                       23, 12, 28,
                                       28, 12, 24,
                                        9, 24, 28,
                                       28, 12,  9,
                                       13, 16, 27,
                                       16, 23, 27,
                                       23, 28, 27,
                                       28, 24, 27,
                                       12, 28, 19,
                                       28, 23, 19,
                                       24, 28, 29,
                                       28,  9, 29,
                                       29,  9, 24,
                                       13, 27, 17,
                                       27, 24, 17,
                                       23, 27, 30,
                                       27, 16, 30,
                                       30, 16, 23,
                                        9, 28, 21,
                                       28, 12, 21,
                                       16, 27, 22,
                                       27, 13, 22,
                                       24, 29, 17,
                                       29,  9, 17,
                                       23, 30, 19,
                                       30, 16, 19,
                                       28, 27, 30,
                                       30, 23, 28,
                                       24, 27, 29,
                                       27, 28, 29,
                                       19, 28, 30,
                                       17, 27, 29,
                                       28, 19, 21,
                                       19, 12, 21,
                                       27, 17,  6,
                                        6, 13, 27,
                                       28, 21, 29,
                                       21,  9, 29,
                                       16, 30, 22,
                                       30, 27, 22,
                                       27,  6, 14,
                                       14, 13, 27,
                                       21, 19,  4,
                                        9, 29, 10,
                                       29, 21, 10,
                                       17, 29,  2,
                                        2, 29,  9,
                                        2, 29, 10,
                                       19, 30,  8,
                                        8, 30, 16,
                                        8, 30, 15,
                                       15, 30, 16,
                                       22, 27, 14,
                                       30, 22, 15,
                                       27, 30, 31,
                                       30, 28, 31,
                                       31, 28, 27,
                                       28, 29, 31,
                                       29, 27, 31,
                                       19, 21, 25,
                                       21, 28, 25,
                                       25, 28, 19,
                                       27, 29, 26,
                                       29, 17, 26,
                                       26, 17, 27,
                                       21, 29, 31,
                                       31, 28, 21,
                                       28, 30, 25,
                                       30, 19, 25,
                                       14, 27, 17,
                                        4, 21, 25,
                                       15, 19, 30,
                                       10, 17, 29,
                                       27, 22, 31,
                                       22, 30, 31,
                                       28, 25, 31,
                                       25, 21, 31,
                                       30, 25, 31,
                                       14, 27, 26,
                                       21, 25, 11,
                                       27, 31, 26,
                                       31, 29, 26,
                                       25, 30, 15,
                                       21, 10, 31,
                                       10, 29, 31,
                                       29, 10, 26,
                                       27, 14, 31,
                                       14, 22, 31,
                                       22, 15, 31,
                                       15, 30, 31,
                                       14, 31, 26,
                                       25, 11, 31,
                                       11, 21, 31,
                                       15, 25, 31,
                                       10, 31, 11,
                                       31, 10, 26,
                                       15, 31,  7,
                                       31, 22,  7,
                                       31, 14,  7,
                                       25, 31, 18,
                                       31, 15, 18,
                                       11, 31, 18,
                                       31, 11,  3,
                                        3, 10, 31,
                                       26, 31, 18,
                                       18, 31, 10,
                                       31, 26,  7,
                                       31,  7, 18,
                                       31, 18,  3 };

  ensure_equals( "total number of entries in inpofa is incorrect",
                 inpofa.size(), correct_inpofa.size() );

  for(std::size_t i=0 ; i<inpofa.size(); ++i)
  {
          ensure_equals("incorrect entry " + std::to_string(i) + " in inpofa",
                          inpofa[i], correct_inpofa[i]-1);
  }
}

//! Generate and test boundary-element vector for tet-only mesh
template<> template<>
void DerivedData_object::test< 68 >() {
  set_test_name( "genBelemTet for tetrahedra" );

  std::map< int, std::vector< std::size_t > > bface;

  // Create unstructured-mesh object to read into
  tk::UnsMesh inmesh;
  // Read in mesh from file
  std::string infile( REGRESSION_DIR
                      "/meshconv/gmsh_output/box_24_ss1.exo" );
  tk::ExodusIIMeshReader er( infile );
  auto nbfac = er.readSidesetFaces( bface );

  // Test if the number of boundary faces is correct
  ensure_equals( "total number of boundary faces incorrect",
                 nbfac, 24 );

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                      10, 14, 13, 12,
                                      14, 13, 12,  9,
                                      10, 14, 12, 11,
                                       1, 14,  5, 11,
                                       7,  6, 10, 12,
                                      14,  8,  5, 10,
                                       8,  7, 10, 13,
                                       7, 13,  3, 12,
                                       1,  4, 14,  9,
                                      13,  4,  3,  9,
                                       3,  2, 12,  9,
                                       4,  8, 14, 13,
                                       6,  5, 10, 11,
                                       1,  2,  9, 11,
                                       2,  6, 12, 11,
                                       6, 10, 12, 11,
                                       2, 12,  9, 11,
                                       5, 14, 10, 11,
                                      14,  8, 10, 13,
                                      13,  3, 12,  9,
                                       7, 10, 13, 12,
                                      14,  4, 13,  9,
                                      14,  1,  9, 11 };

  // Boundary-face node connectivity
  std::vector< std::size_t > triinpoel {  2,  9,  3,
                                          1,  9,  2,
                                          3,  9,  4,
                                          1,  4,  9,
                                          7, 10,  6,
                                          6, 10,  5,
                                          8, 10,  7,
                                         10,  8,  5,
                                          6, 11,  2,
                                          2, 11,  1,
                                         11,  6,  5,
                                         11,  5,  1,
                                          3, 12,  2,
                                         12,  6,  2,
                                          7, 12,  3,
                                         12,  7,  6,
                                         13,  7,  3,
                                          4, 13,  3,
                                         13,  8,  7,
                                          8, 13,  4,
                                          5,  8, 14,
                                          1,  5, 14,
                                          4, 14,  8,
                                          1, 14,  4 };

  // Correct boundary elements
  std::vector< std::size_t > correct_belem{ 12,
                                            15,
                                            11,
                                            10,
                                             6,
                                            14,
                                             8,
                                             7,
                                            16,
                                            15,
                                            14,
                                             5,
                                            12,
                                            16,
                                             9,
                                             6,
                                             9,
                                            11,
                                             8,
                                            13,
                                             7,
                                             5,
                                            13,
                                            10 };

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  tk::shiftToZero( triinpoel );

  auto esup = tk::genEsup( inpoel, 4 );
  auto esuel = tk::genEsuelTet( inpoel, esup );
  auto ntfac = tk::genNtfac( 4, nbfac, esuel );
  auto inpofa = tk::genInpofaTet( ntfac, nbfac, inpoel, triinpoel, esuel );
  auto belem = tk::genBelemTet( nbfac, inpofa, esup );

  ensure_equals( "total number of entries in belem is incorrect",
                 belem.size(), correct_belem.size() );

  for(std::size_t i=0 ; i<belem.size(); ++i)
  {
    ensure_equals("incorrect entry " + std::to_string(i) + " in belem",
                    belem[i], correct_belem[i]-1);
  }
}

//! Generate and test face-geometry vector for a single tetrahedron
template<> template<>
void DerivedData_object::test< 69 >() {
  set_test_name( "Face-geometry (genGeoFaceTri) for a tetrahedron" );

  // total number of faces
  std::size_t ntfac(4);

  // coordinates of tetrahedron vertices
  tk::UnsMesh::Coords coord {{ {1.0, 0.0, 0.0, 0.0},
                               {0.0, 0.0, 1.0, 0.0},
                               {0.0, 0.0, 0.0, 1.0} }};

  // face-node connectivity
  std::vector< std::size_t > inpofa { 0, 1, 2,
                                      0, 3, 1,
                                      1, 3, 2,
                                      2, 3, 0 };

  // get face-geometries
  auto geoFace = tk::genGeoFaceTri( ntfac, inpofa, coord );

  // correct face-areas
  std::vector< tk::real > correct_farea { 0.5, 0.5, 0.5, 0.8660254037844389 };

  // correct face-normals
  std::array< std::vector< tk::real >, 3 > correct_fnorm {{
                                       { 0.0,  0.0, -1.0, 1.0/sqrt(3.0)},
                                       { 0.0, -1.0,  0.0, 1.0/sqrt(3.0)},
                                       {-1.0,  0.0,  0.0, 1.0/sqrt(3.0)},
                                       }};

  // correct face-centroids
  tk::UnsMesh::Coords correct_fcent {{ {1.0/3.0, 1.0/3.0, 0.0,     1.0/3.0},
                                       {1.0/3.0, 0.0,     1.0/3.0, 1.0/3.0},
                                       {0.0,     1.0/3.0, 1.0/3.0, 1.0/3.0} }};

  tk::real prec = std::numeric_limits< tk::real >::epsilon();

  for(std::size_t f=0 ; f<ntfac; ++f)
  {
    ensure_equals("incorrect entry " + std::to_string(f) + " in geoFace-area",
                    geoFace(f,0,0), correct_farea[f], prec);

    ensure_equals("incorrect entry " + std::to_string(f) + " in geoFace-nx",
                    geoFace(f,1,0), correct_fnorm[0][f], prec);

    ensure_equals("incorrect entry " + std::to_string(f) + " in geoFace-ny",
                    geoFace(f,2,0), correct_fnorm[1][f], prec);

    ensure_equals("incorrect entry " + std::to_string(f) + " in geoFace-nz",
                    geoFace(f,3,0), correct_fnorm[2][f], prec);

    ensure_equals("incorrect entry " + std::to_string(f) + " in geoFace-cx",
                    geoFace(f,4,0), correct_fcent[0][f], prec);

    ensure_equals("incorrect entry " + std::to_string(f) + " in geoFace-cy",
                    geoFace(f,5,0), correct_fcent[1][f], prec);

    ensure_equals("incorrect entry " + std::to_string(f) + " in geoFace-cz",
                    geoFace(f,6,0), correct_fcent[2][f], prec);
  }
}

//! Generate and test element-geometry vector for a single tetrahedron
template<> template<>
void DerivedData_object::test< 70 >() {
  set_test_name( "Element-geometry (genGeoElemTet) for a tetrahedron" );

  // coordinates of tetrahedron vertices
  tk::UnsMesh::Coords coord {{ {1.0, 0.0, 0.0, 0.0},
                               {0.0, 0.0, 1.0, 0.0},
                               {0.0, 0.0, 0.0, 1.0} }};

  // element-node connectivity
  std::vector< std::size_t > inpoel { 0, 3, 2, 1 };

  // get element-geometries
  auto geoElem = tk::genGeoElemTet( inpoel, coord );

  // correct element-volume
  tk::real correct_vole { 1.0/6.0 };

  // correct element-centroid
  tk::UnsMesh::Coords correct_ecent {{ {1.0/4.0},
                                       {1.0/4.0},
                                       {1.0/4.0} }};

  tk::real prec = std::numeric_limits< tk::real >::epsilon();

  ensure_equals("incorrect entry in geoElem-vol",
                  geoElem(0,0,0), correct_vole, prec);

  ensure_equals("incorrect entry in geoElem-cx",
                  geoElem(0,1,0), correct_ecent[0][0], prec);

  ensure_equals("incorrect entry in geoElem-cy",
                  geoElem(0,2,0), correct_ecent[1][0], prec);

  ensure_equals("incorrect entry in geoElem-cz",
                  geoElem(0,3,0), correct_ecent[2][0], prec);
}

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

} // tut::

#endif // test_DerivedData_h
