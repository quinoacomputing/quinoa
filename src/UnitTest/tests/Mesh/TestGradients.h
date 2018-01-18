// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Mesh/TestGradients.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Unit tests for Mesh/Gradients
  \details   Unit tests for Mesh/Gradients. All unit tests start from simple
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
*/
// *****************************************************************************
#ifndef test_Gradients_h
#define test_Gradients_h

#include <algorithm>

#include "NoWarning/tut.h"

#include "Gradients.h"
#include "Fields.h"

namespace tut {

//! All tests in group inherited from this base
struct Gradients_common {
  const tk::real pr = 5.0*std::numeric_limits< tk::real >::epsilon();

  // mesh node coordinates
  std::array< std::vector< tk::real >, 3 > coord {{
    {{ 0, 1, 1, 0, 0, 1, 1, 0, 0.5, 0.5, 0.5, 1, 0.5, 0 }},
    {{ 0, 0, 1, 1, 0, 0, 1, 1, 0.5, 0.5, 0, 0.5, 1, 0.5 }},
    {{ 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0.5, 0.5, 0.5, 0.5 }} }};

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
};

// Test group shortcuts
// The 2nd template argument is the max number of tests in this group. If
// omitted, the default is 50, specified in tut/tut.hpp.
using Gradients_group = test_group< Gradients_common, MAX_TESTS_IN_GROUP >;
using Gradients_object = Gradients_group::object;

//! Define test group
static Gradients_group Gradients( "Mesh/Gradients" );

//! Test definitions for group

//! Test nodal gradients for tetrahedron-only mesh
template<> template<>
void Gradients_object::test< 1 >() {
  set_test_name( "node gradients of tetrahedra mesh" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  // Generate elements surrounding points
  auto esup = tk::genEsup( inpoel, 4 );

  // generate a constant scalar field
  tk::Fields uc( npoin, 1 );
  uc.fill( 1.2 );
  // test gradients
  for (std::size_t p=0; p<npoin; ++p) {
    auto g = nodegrad( p, coord, inpoel, esup, uc, 0 );
    ensure_equals( "x-gradient of constant field incorrect", g[0], 0, pr );
    ensure_equals( "y-gradient of constant field incorrect", g[1], 0, pr );
    ensure_equals( "z-gradient of constant field incorrect", g[2], 0, pr );
  }

  // generate a linear scalar field with a slope in only x direction
  tk::Fields ux( npoin, 1 );
  for (std::size_t p=0; p<npoin; ++p) ux(p,0,0) = coord[0][p];
  // test gradients
  for (std::size_t p=0; p<npoin; ++p) {
    auto g = nodegrad( p, coord, inpoel, esup, ux, 0 );
    ensure_equals( "x-gradient of x-linear field incorrect", g[0], 1, pr );
    ensure_equals( "y-gradient of x-linear field incorrect", g[1], 0, pr );
    ensure_equals( "z-gradient of x-linear field incorrect", g[2], 0, pr );
  }

  // generate a linear scalar field with a slope in only y direction
  tk::Fields uy( npoin, 1 );
  for (std::size_t p=0; p<npoin; ++p) uy(p,0,0) = coord[1][p];
  // test gradients
  for (std::size_t p=0; p<npoin; ++p) {
    auto g = nodegrad( p, coord, inpoel, esup, uy, 0 );
    ensure_equals( "x-gradient of y-linear field incorrect", g[0], 0, pr );
    ensure_equals( "y-gradient of y-linear field incorrect", g[1], 1, pr );
    ensure_equals( "z-gradient of y-linear field incorrect", g[2], 0, pr );
  }

  // generate a linear scalar field with a slope in only z direction
  tk::Fields uz( npoin, 1 );
  for (std::size_t p=0; p<npoin; ++p) uz(p,0,0) = coord[2][p];
  // test gradients
  for (std::size_t p=0; p<npoin; ++p) {
    auto g = nodegrad( p, coord, inpoel, esup, uz, 0 );
    ensure_equals( "x-gradient of z-linear field incorrect", g[0], 0, pr );
    ensure_equals( "y-gradient of z-linear field incorrect", g[1], 0, pr );
    ensure_equals( "z-gradient of z-linear field incorrect", g[2], 1, pr );
  }

  // generate linear vector field with different slopes for different components
  tk::Fields u3( npoin, 3 );
  for (std::size_t p=0; p<npoin; ++p) {
     u3(p,0,0) = 2.0*coord[0][p];
     u3(p,1,0) = 1.5*coord[1][p];
     u3(p,2,0) = -0.5*coord[2][p];
  }
  // test gradients
  for (std::size_t p=0; p<npoin; ++p) {
    auto gx = nodegrad( p, coord, inpoel, esup, u3, 0 );
    ensure_equals( "x-gradient of x-linear vector incorrect", gx[0], 2.0, pr );
    ensure_equals( "y-gradient of x-linear vector incorrect", gx[1], 0, pr );
    ensure_equals( "z-gradient of x-linear vector incorrect", gx[2], 0, pr );
    auto gy = nodegrad( p, coord, inpoel, esup, u3, 1 );
    ensure_equals( "x-gradient of y-linear vector incorrect", gy[0], 0, pr );
    ensure_equals( "y-gradient of y-linear vector incorrect", gy[1], 1.5, pr );
    ensure_equals( "z-gradient of y-linear vector incorrect", gy[2], 0, pr );
    auto gz = nodegrad( p, coord, inpoel, esup, u3, 2 );
    ensure_equals( "x-gradient of z-linear vector incorrect", gz[0], 0, pr );
    ensure_equals( "y-gradient of z-linear vector incorrect", gz[1], 0, pr );
    ensure_equals( "z-gradient of z-linear vector incorrect", gz[2], -0.5, pr );
  }
}

//! Test edge gradients for tetrahedron-only mesh
template<> template<>
void Gradients_object::test< 2 >() {
  set_test_name( "edge gradients of tetrahedra mesh" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );

  // find out number of points in mesh connectivity
  auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
  Assert( *minmax.first == 0, "node ids should start from zero" );
  auto npoin = *minmax.second + 1;

  // Generate elements surrounding points
  auto esup = tk::genEsup( inpoel, 4 );
  // Generate elements surrounding edges
  auto esued = tk::genEsued( inpoel, 4, esup );

  // find out number of edges in mesh
  auto nedge = tk::genInpoed(inpoel,4,esup).size()/2;

  // generate a constant scalar field
  tk::Fields uc( npoin, 1 );
  uc.fill( 1.2 );
  // test gradients
  for (std::size_t e=0; e<nedge; ++e) {
    auto g = edgegrad( e, coord, inpoel, esued, uc, 0 );
    ensure_equals( "x-gradient of constant field incorrect", g[0], 0, pr );
    ensure_equals( "y-gradient of constant field incorrect", g[1], 0, pr );
    ensure_equals( "z-gradient of constant field incorrect", g[2], 0, pr );
  }

  // generate a linear scalar field with a slope in only x direction
  tk::Fields ux( npoin, 1 );
  for (std::size_t p=0; p<npoin; ++p) ux(p,0,0) = coord[0][p];
  // test gradients
  for (std::size_t e=0; e<nedge; ++e) {
    auto g = edgegrad( e, coord, inpoel, esued, ux, 0 );
    ensure_equals( "x-gradient of x-linear field incorrect", g[0], 1, pr );
    ensure_equals( "y-gradient of x-linear field incorrect", g[1], 0, pr );
    ensure_equals( "z-gradient of x-linear field incorrect", g[2], 0, pr );
  }

  // generate a linear scalar field with a slope in only y direction
  tk::Fields uy( npoin, 1 );
  for (std::size_t p=0; p<npoin; ++p) uy(p,0,0) = coord[1][p];
  // test gradients
  for (std::size_t e=0; e<nedge; ++e) {
    auto g = edgegrad( e, coord, inpoel, esued, uy, 0 );
    ensure_equals( "x-gradient of y-linear field incorrect", g[0], 0, pr );
    ensure_equals( "y-gradient of y-linear field incorrect", g[1], 1, pr );
    ensure_equals( "z-gradient of y-linear field incorrect", g[2], 0, pr );
  }

  // generate a linear scalar field with a slope in only z direction
  tk::Fields uz( npoin, 1 );
  for (std::size_t p=0; p<npoin; ++p) uz(p,0,0) = coord[2][p];
  // test gradients
  for (std::size_t e=0; e<nedge; ++e) {
    auto g = edgegrad( e, coord, inpoel, esued, uz, 0 );
    ensure_equals( "x-gradient of z-linear field incorrect", g[0], 0, pr );
    ensure_equals( "y-gradient of z-linear field incorrect", g[1], 0, pr );
    ensure_equals( "z-gradient of z-linear field incorrect", g[2], 1, pr );
  }

  // generate linear vector field with different slopes for different components
  tk::Fields u3( npoin, 3 );
  for (std::size_t p=0; p<npoin; ++p) {
     u3(p,0,0) = 2.0*coord[0][p];
     u3(p,1,0) = 1.5*coord[1][p];
     u3(p,2,0) = -0.5*coord[2][p];
  }
  // test gradients
  for (std::size_t e=0; e<nedge; ++e) {
    auto gx = edgegrad( e, coord, inpoel, esued, u3, 0 );
    ensure_equals( "x-gradient of x-linear vector incorrect", gx[0], 2.0, pr );
    ensure_equals( "y-gradient of x-linear vector incorrect", gx[1], 0, pr );
    ensure_equals( "z-gradient of x-linear vector incorrect", gx[2], 0, pr );
    auto gy = edgegrad( e, coord, inpoel, esued, u3, 1 );
    ensure_equals( "x-gradient of y-linear vector incorrect", gy[0], 0, pr );
    ensure_equals( "y-gradient of y-linear vector incorrect", gy[1], 1.5, pr );
    ensure_equals( "z-gradient of y-linear vector incorrect", gy[2], 0, pr );
    auto gz = edgegrad( e, coord, inpoel, esued, u3, 2 );
    ensure_equals( "x-gradient of z-linear vector incorrect", gz[0], 0, pr );
    ensure_equals( "y-gradient of z-linear vector incorrect", gz[1], 0, pr );
    ensure_equals( "z-gradient of z-linear vector incorrect", gz[2], -0.5, pr );
  }
}

} // tut::

#endif // test_DerivedData_h
