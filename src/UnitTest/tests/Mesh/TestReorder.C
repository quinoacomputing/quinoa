// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Mesh/TestReorder.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Unit tests for Mesh/Reorder
  \details   Unit tests for Mesh/Reorder. All unit tests start from simple mesh
     connectivities defined in the code. The tetrahedron mesh in Gmsh ASCII
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

#include "NoWarning/tut.h"

#include "TUTConfig.h"
#include "Reorder.h"
#include "DerivedData.h"

namespace tut {

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif

//! All tests in group inherited from this base
struct Reorder_common {

  // Mesh node coordinates
  std::array< std::vector< tk::real >, 3 > tetcoord {{
    {{ 0, 1, 1, 0, 0, 1, 1, 0, 0.5, 0.5, 0.5, 1, 0.5, 0 }},
    {{ 0, 0, 1, 1, 0, 0, 1, 1, 0.5, 0.5, 0, 0.5, 1, 0.5 }},
    {{ 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0.5, 0.5, 0.5, 0.5 }} }};

  // Mesh connectivity for simple tetrahedron-only mesh
  std::vector< std::size_t > tetinpoel { 12, 14,  9, 11,
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
};

// Test group shortcuts
// The 2nd template argument is the max number of tests in this group. If
// omitted, the default is 50, specified in tut/tut.hpp.
using Reorder_group = test_group< Reorder_common, MAX_TESTS_IN_GROUP >;
using Reorder_object = Reorder_group::object;

//! Define test group
static Reorder_group Reorder( "Mesh/Reorder" );

//! Test definitions for group

//! Attempt to shift empty container using shiftToZero
template<> template<>
void Reorder_object::test< 1 >() {
  set_test_name( "shiftToZero graceful with empty inpoel" );

  // Attempt to shift node IDs with empty connectivity. If some error happens or an
  // exception is throw that will go to the screen; no further tests are
  // necessary.
  std::vector< std::size_t > empty;
  tk::shiftToZero( empty );
}

//! Shift node ids to zero in line mesh
template<> template<>
void Reorder_object::test< 2 >() {
  set_test_name( "shiftToZero for lines" );

  // Mesh connectivity for simple line-only mesh
  std::vector< std::size_t > inpoel { 1, 2,
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
  auto min = std::min_element( begin(inpoel), end(inpoel) );
  ensure_equals( "node ids should start from zero", *min, 0 );
}

//! Shift node ids to zero in triangle mesh
template<> template<>
void Reorder_object::test< 3 >() {
  set_test_name( "shiftToZero for triangles" );

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

  // Test new extents of node IDs in element connectivity
  auto min = std::min_element( begin(inpoel), end(inpoel) );
  ensure_equals( "node ids should start from zero", *min, 0 );
}

//! Shift node ids to zero in tetrahedron-only mesh
template<> template<>
void Reorder_object::test< 4 >() {
  set_test_name( "shiftToZero for tetrahedra" );

  // Shift node IDs to start from zero
  auto inpoel = tetinpoel;
  tk::shiftToZero( inpoel );

  // Test new extents of node IDs in element connectivity
  auto min = std::min_element( begin(inpoel), end(inpoel) );
  ensure_equals( "node ids should start from zero", *min, 0 );
}

//! Renumber triangle mesh
template<> template<>
void Reorder_object::test< 5 >() {
  set_test_name( "renumber triangle mesh" );

  // Mesh connectivity for simple triangle mesh
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

  // Renumber triangle mesh
  const auto psup = tk::genPsup( inpoel, 3, tk::genEsup( inpoel, 3 ) );
  auto map = tk::renumber( psup );
  tk::remap( inpoel, map );

  // Test the result of reordering
  std::vector< std::size_t > correct_renumbered_inpoel{ 0, 4, 1,
                                                        0, 2, 4,
                                                        1, 4, 7,
                                                        7, 4, 2,
                                                        3, 8, 12,
                                                        3, 12, 10,
                                                        8, 13, 12,
                                                        13, 10, 12,
                                                        0, 1, 5,
                                                        0, 5, 3,
                                                        1, 8, 5,
                                                        3, 5, 8,
                                                        1, 7, 9,
                                                        1, 9, 8,
                                                        7, 13, 9,
                                                        8, 9, 13,
                                                        7, 2, 11,
                                                        7, 11, 13,
                                                        2, 10, 11,
                                                        13, 11, 10,
                                                        0, 6, 2,
                                                        0, 3, 6,
                                                        2, 6, 10,
                                                        3, 10, 6 };
  ensure( "reordered triangle mesh incorrect",
          inpoel == correct_renumbered_inpoel );
}

//! Renumber tetrahedron mesh
template<> template<>
void Reorder_object::test< 6 >() {
  set_test_name( "renumber tetrahedron mesh" );

  // Shift node IDs to start from zero
  auto inpoel = tetinpoel;
  tk::shiftToZero( inpoel );

  // Renumber triangle mesh
  const auto psup = tk::genPsup( inpoel, 4, tk::genEsup( inpoel, 4 ) );
  auto map = tk::renumber( psup );
  tk::remap( inpoel, map );

  // Test the result of reordering
  std::vector< std::size_t > correct_renumbered_inpoel{ 9, 6, 4, 5,
                                                        12, 6, 11, 9,
                                                        6, 11, 9, 4,
                                                        12, 6, 9, 5,
                                                        0, 6, 3, 5,
                                                        13, 8, 12, 9,
                                                        6, 10, 3, 12,
                                                        10, 13, 12, 11,
                                                        13, 11, 7, 9,
                                                        0, 2, 6, 4,
                                                        11, 2, 7, 4,
                                                        7, 1, 9, 4,
                                                        2, 10, 6, 11,
                                                        8, 3, 12, 5,
                                                        0, 1, 4, 5,
                                                        1, 8, 9, 5,
                                                        8, 12, 9, 5,
                                                        1, 9, 4, 5,
                                                        3, 6, 12, 5,
                                                        6, 10, 12, 11,
                                                        11, 7, 9, 4,
                                                        13, 12, 11, 9,
                                                        6, 2, 11, 4,
                                                        6, 0, 4, 5 };

  ensure( "reordered tetrahedron mesh incorrect",
          inpoel == correct_renumbered_inpoel );
}

//! Test all positive Jacbians in tetrahedron mesh
template<> template<>
void Reorder_object::test< 7 >() {
  set_test_name( "all positive Jacobians" );

  // Shift node IDs to start from zero
  auto inpoel = tetinpoel;
  tk::shiftToZero( inpoel );

  ensure( "Not all Jacobians are positive",
          tk::positiveJacobians( inpoel, tetcoord ) );
}

//! Test not all positive Jacbians in tetrahedron mesh
template<> template<>
void Reorder_object::test< 8 >() {
  set_test_name( "not all positive Jacobians" );

  // Shift node IDs to start from zero
  auto inpoel = tetinpoel;
  tk::shiftToZero( inpoel );

  // Switch two vertices of a tet
  std::swap( inpoel[4*4+0], inpoel[4*4+1] );

  ensure( "All Jacobians are positive",
          !tk::positiveJacobians( inpoel, tetcoord ) );
}

//! \brief Test if positiveJacobians throws on inpoel non-divisible by the
//!   number of nodes per elements
template<> template<>
void Reorder_object::test< 9 >() {
  set_test_name( "positiveJacobians throws on inpoel non-div nnpe" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield invalid read" );
  #else
  try {
    // Partial mesh mesh connectivity
    std::vector< std::size_t > inpoel { 12, 14,  9, 11,
                                        14,  4, 13 };
    tk::positiveJacobians( inpoel, tetcoord );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test if positiveJacobians throws on empty inpoel
template<> template<>
void Reorder_object::test< 10 >() {
  set_test_name( "positivaJacobians throws with empty inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    std::vector< std::size_t > empty;
    tk::positiveJacobians( empty, tetcoord );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test if positiveJacobians throws on empty coord
template<> template<>
void Reorder_object::test< 11 >() {
  set_test_name( "positivaJacobians throws with empty coord" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    decltype(tetcoord) empty;
    tk::positiveJacobians( tetinpoel, empty );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test if positiveJacobians throws on unique(inpoel).size != coord.size
template<> template<>
void Reorder_object::test< 12 >() {
  set_test_name( "positivaJacobians throws w inconsistent inpoel & coord" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    // Shift node IDs to start from zero
    auto inpoel = tetinpoel;
    tk::shiftToZero( inpoel );
    // Remove last coordinate from coord[0]
    auto coord = tetcoord;
    coord[0].pop_back();
    tk::positiveJacobians( tetinpoel, coord );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

//! Test if positiveJacobians throws with non-zero-based inpoel
template<> template<>
void Reorder_object::test< 13 >() {
  set_test_name( "positivaJacobians throws with non-zero-based inpoel" );

  #ifdef NDEBUG        // exception only thrown in DEBUG mode
    skip( "in RELEASE mode, would yield segmentation fault" );
  #else
  try {
    // Do not shift node IDs to start from zero
    auto inpoel = tetinpoel;
    tk::positiveJacobians( inpoel, tetcoord );
    fail( "should throw exception in DEBUG mode" );
  }
  catch ( tk::Exception& ) {
    // exception thrown in DEBUG mode, test ok
  }
  #endif
}

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

} // tut::
