// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Mesh/TestAround.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Unit tests for Mesh/Around
  \details   Unit tests for Mesh/Around.
*/
// *****************************************************************************
#ifndef test_Around_h
#define test_Around_h

#include "NoWarning/tut.h"

#include "Around.h"

namespace tut {

//! All tests in group inherited from this base
struct Around_common {

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

  // Mesh connectivity for simple triangle-only mesh
  std::vector< std::size_t > triinpoel { 1,  9,  2,
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
};

// Test group shortcuts
// The 2nd template argument is the max number of tests in this group. If
// omitted, the default is 50, specified in tut/tut.hpp.
using Around_group = test_group< Around_common, MAX_TESTS_IN_GROUP >;
using Around_object = Around_group::object;

//! Define test group
static Around_group Around( "Mesh/Around" );

//! Test definitions for group

//! Iterate via Around on elements surrounding points for tetrahedron-only mesh
template<> template<>
void Around_object::test< 1 >() {
  set_test_name( "Esup for tetrahedra" );

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
    // Iterate through all elements surrounding point p
    for (auto e : tk::Around(esup,p)) points.push_back(e);
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

//! Iterate via Around on elements surrounding points for triangle-only mesh
template<> template<>
void Around_object::test< 2 >() {
  set_test_name( "Esup for triangles" );

  // Shift node IDs to start from zero
  tk::shiftToZero( triinpoel );

  // Generate elements surrounding points
  auto esup = tk::genEsup( triinpoel, 3 );

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
  auto minmax = std::minmax_element( begin(triinpoel), end(triinpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  // this is more of a test on this test
  ensure_equals( "number of points in 'correct' esup incorrect",
                 npoin, correct_esup.size() );

  // test generated derived data structure, elements surrounding points
  for (std::size_t p=0; p<npoin; ++p) {
    // extract element ids from generated elements surrounding point p
    std::vector< std::size_t > points;
    // Iterate through all elements surrounding point p
    for (auto e : tk::Around(esup,p)) points.push_back(e);
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

//! Iterate via Around on points surrounding points for tetrahedron-only mesh
template<> template<>
void Around_object::test< 3 >() {
  set_test_name( "Psup for tetrahedra" );

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
    // Iterate through all points surrounding point p
    for (auto i : tk::Around(psup,p)) points.push_back(i);
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

//! Iterate via Around on points surrounding points for triangle-only mesh
template<> template<>
void Around_object::test< 4 >() {
  set_test_name( "Psup for triangles" );

  // Shift node IDs to start from zero
  tk::shiftToZero( triinpoel );

  // Generate elements surrounding points
  auto psup = tk::genPsup( triinpoel, 3, tk::genEsup(triinpoel,3) );

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
  auto minmax = std::minmax_element( begin(triinpoel), end(triinpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  // this is more of a test on this test
  ensure_equals( "number of points in psup incorrect",
                 npoin, correct_psup.size() );

  // test generated derived data structure, elements surrounding points
  for (std::size_t p=0; p<npoin; ++p) {
    // extract element ids from generated elements surrounding point p
    std::vector< std::size_t > points;
    // Iterate through all points surrounding point p
    for (auto i : tk::Around(psup,p)) points.push_back(i);
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

//! Iterate via Around on edges surrounding points for tetrahedron-only mesh
template<> template<>
void Around_object::test< 5 >() {
  set_test_name( "Edsup for tetrahedra" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Generate edges surrounding points
  auto edsup = tk::genEdsup( inpoel, 4, tk::genEsup(inpoel,4) );

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
    // Iterate through all edges surrounding point p
    for (auto i : tk::Around(edsup,p)) edge.push_back(i);
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

//! Iterate via Around on edges surrounding points for triangle-only mesh
template<> template<>
void Around_object::test< 6 >() {
  set_test_name( "Edsup for triangles" );

  // Shift node IDs to start from zero
  tk::shiftToZero( triinpoel );

  // Generate edges surrounding points
  auto edsup = tk::genEdsup( triinpoel, 3, tk::genEsup(triinpoel,3) );

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
  auto minmax = std::minmax_element( begin(triinpoel), end(triinpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  // this is more of a test on this test
  ensure_equals( "number of points (star centers) in edsup incorrect",
                 npoin, correct_edsup.size() );

  // Test generated derived data structure, edges surrounding points
  for (std::size_t p=0; p<npoin; ++p) {
    // extract edge end-point ids from generated edges surrounding points
    std::vector< std::size_t > edge;
    // Iterate through all edges surrounding point p
    for (auto i : tk::Around(edsup,p)) edge.push_back(i);
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


//! \brief Iterate via Around on elements surrounding points of elements for
//!    tetrahedron-only mesh
template<> template<>
void Around_object::test< 7 >() {
  set_test_name( "Esupel for tetrahedra" );

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
    // Iterate through all elements surrounding points of elements e
    for (auto i : tk::Around(esupel,e)) elements.push_back(i);
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

//! \brief Iterate via Around on elements surrounding points of elements for
//!    triangle-only mesh
template<> template<>
void Around_object::test< 8 >() {
  set_test_name( "Esupel for triangles" );

  // Shift node IDs to start from zero
  tk::shiftToZero( triinpoel );

  // Generate elements surrounding points
  auto esupel = tk::genEsupel( triinpoel, 3, tk::genEsup(triinpoel,3) );

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
  auto nelem = triinpoel.size()/3;

  // this is more of a test on this test
  ensure_equals( "number of elements in esupel incorrect",
                 nelem, correct_esupel.size() );

  // test generated derived data structure, elements surrounding points of
  // elements
  for (std::size_t e=0; e<nelem; ++e) {
    // extract element ids from generated elements surrounding points of
    // elements
    std::vector< std::size_t > elements;
    // Iterate through all elements surrounding points of elements e
    for (auto i : tk::Around(esupel,e)) elements.push_back(i);
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

//! \brief Iterate via Around on elements surrounding elements for
//!    tetrahedron-only mesh
template<> template<>
void Around_object::test< 9 >() {
  set_test_name( "Esuel for tetrahedra" );

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
    // Iterate through all elements surrounding element e
    for (auto i : tk::Around(esuel,e)) elements.push_back(i);
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

//! \brief Iterate via Around on elements surrounding elements for
//!    triangle-only mesh
template<> template<>
void Around_object::test< 10 >() {
  set_test_name( "Esuel for triangles" );

  // Shift node IDs to start from zero
  tk::shiftToZero( triinpoel );

  // Generate elements surrounding points
  auto esuel = tk::genEsuel( triinpoel, 3, tk::genEsup(triinpoel,3) );

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
  auto nelem = triinpoel.size()/3;

  // this is more of a test on this test
  ensure_equals( "number of elements in esuel incorrect",
                 nelem, correct_esuel.size() );

  // test generated derived data structure, elements surrounding elements
  for (std::size_t e=0; e<nelem; ++e) {
    // extract element ids from generated elements surrounding elements
    std::vector< std::size_t > elements;
    // Iterate through all elements surrounding element e
    for (auto i : tk::Around(esuel,e)) elements.push_back(i);
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

//! Iterate via Around on elements surrounding edges for tetrahedron-only mesh
template<> template<>
void Around_object::test< 11 >() {
  set_test_name( "Esued for tetrahedra" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // Generate edges surrounding points
  auto esup = tk::genEsup( inpoel, 4 );
  auto esued = tk::genEsued( inpoel, 4, esup );

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
    // Iterate through all elements surrounding edges e
    for (auto i : tk::Around(esued,e)) elements.push_back(i);
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

//!Iterate via Around on elements surrounding edges for triangle-only mesh
template<> template<>
void Around_object::test< 12 >() {
  set_test_name( "Esued for triangles" );

  // Shift node IDs to start from zero
  tk::shiftToZero( triinpoel );

  // Generate edges surrounding points
  auto esup = tk::genEsup( triinpoel, 3 );
  auto esued = tk::genEsued( triinpoel, 3, esup );

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
  auto nedge = tk::genInpoed(triinpoel,3,esup).size()/2;

  // this is more of a test on this test
  ensure_equals( "number of edges in esued incorrect",
                 nedge, correct_esued.size() );

  // Test generated derived data structure, elements surrounding edges
  for (std::size_t e=0; e<nedge; ++e) {
    // extract element list generated for edge e
    std::vector< std::size_t > elements;
    // Iterate through all elements surrounding edges e
    for (auto i : tk::Around(esued,e)) elements.push_back(i);
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

} // tut::

#endif // test_Around_h
