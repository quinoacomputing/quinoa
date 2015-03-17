//******************************************************************************
/*!
  \file      src/UnitTest/tests/Mesh/DerivedData.h
  \author    J. Bakosi
  \date      Tue 17 Mar 2015 07:53:56 AM MDT
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

//! Test correction of node ordering for line mesh
template<> template<>
void DerivedData_object::test< 1 >() {
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

//! Test correction of node ordering for triangle mesh
template<> template<>
void DerivedData_object::test< 2 >() {
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

//! Test correction of node ordering for tetrahedron-only mesh
template<> template<>
void DerivedData_object::test< 3 >() {
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

} // tut::

#endif // test_DerivedData_h
