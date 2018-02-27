// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Inciter/AMR/TestError.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Unit tests for AMR error indicators in Inciter/AMR/Error.h
  \details   Unit tests for AMR error indicators in Inciter/AMR/Error.h. All
     unit tests start from simple mesh connectivities defined in the code. The
     tetrahedron mesh in Gmsh ASCII format is as follows. Note that ids start
     from zero in the code, but from one in Gmsh.
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
#ifndef test_AMRError_h
#define test_AMRError_h

#include <limits>

#include "NoWarning/tut.h"

#include "Types.h"
#include "Fields.h"
#include "Reorder.h"
#include "DerivedData.h"
#include "AMR/Error.h"

namespace tut {

//! All tests in group inherited from this base
struct AMRError_common {
  // floating point precision tolerance
  const tk::real pr = std::numeric_limits< tk::real >::epsilon();

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

  //! Test jump error indicator (mostly bounds) for tetrahedron mesh
  void TestErrorIndicator( inciter::ctr::AMRErrorType errtype ) {
    // Shift node IDs to start from zero
    tk::shiftToZero( inpoel );
  
    // find out number of points in mesh connectivity
    auto minmax = std::minmax_element( begin(inpoel), end(inpoel) );
    Assert( *minmax.first == 0, "node ids should start from zero" );
    auto npoin = *minmax.second + 1;
  
    // Generate elements surrounding points
    auto esup = tk::genEsup( inpoel, 4 );
  
    // Generate edge connectivity
    auto inpoed = tk::genInpoed( inpoel, 4, esup );
  
    // create error indicator object to get access to error indicators
    AMR::Error err;
  
    // generate a constant scalar field
    tk::Fields uc( npoin, 1 );
    uc.fill( 1.2 );
    // test jump error indicator on all edges
    for (std::size_t e=0; e<inpoed.size()/2; ++e) {
      std::pair< std::size_t, std::size_t > edge{ inpoed[e*2],inpoed[e*2+1] };
      auto r = err.scalar( uc, edge, 0, coord, inpoel, esup, errtype );
      ensure_equals( "edge error of constant field incorrect", r, 0.0, pr );
    }
  
    // generate a linear scalar field with a slope in only x direction
    tk::Fields ux( npoin, 1 );
    for (std::size_t p=0; p<npoin; ++p) ux(p,0,0) = coord[0][p];
    // test jump error indicator on all edges
    for (std::size_t e=0; e<inpoed.size()/2; ++e) {
      std::pair< std::size_t, std::size_t > edge{ inpoed[e*2],inpoed[e*2+1] };
      auto r = err.scalar( ux, edge, 0, coord, inpoel, esup, errtype );
      ensure( "edge error < 0.0", r > -pr );
      ensure( "edge error > 1.0", r < 1.0+pr );
    }
  
    // generate a linear scalar field with a slope in only y direction
    tk::Fields uy( npoin, 1 );
    for (std::size_t p=0; p<npoin; ++p) uy(p,0,0) = coord[1][p];
    // test jump error indicator on all edges
    for (std::size_t e=0; e<inpoed.size()/2; ++e) {
      std::pair< std::size_t, std::size_t > edge{ inpoed[e*2],inpoed[e*2+1] };
      auto r = err.scalar( uy, edge, 0, coord, inpoel, esup, errtype );
      ensure( "edge error < 0.0", r > -pr );
      ensure( "edge error > 1.0", r < 1.0+pr );
    }
  
    // generate a linear scalar field with a slope in only z direction
    tk::Fields uz( npoin, 1 );
    for (std::size_t p=0; p<npoin; ++p) uz(p,0,0) = coord[2][p];
    // test jump error indicator on all edges
    for (std::size_t e=0; e<inpoed.size()/2; ++e) {
      std::pair< std::size_t, std::size_t > edge{ inpoed[e*2],inpoed[e*2+1] };
      auto r = err.scalar( uz, edge, 0, coord, inpoel, esup, errtype );
      ensure( "edge error < 0.0", r > -pr );
      ensure( "edge error > 1.0", r < 1.0+pr );
    }
  
    // generate linear vector field with different slopes for different components
    tk::Fields u3( npoin, 3 );
    for (std::size_t p=0; p<npoin; ++p) {
       u3(p,0,0) = 2.0*coord[0][p];
       u3(p,1,0) = 1.5*coord[1][p];
       u3(p,2,0) = -0.5*coord[2][p];
    }
    // test jump error indicator on all edges
    for (std::size_t e=0; e<inpoed.size()/2; ++e) {
      std::pair< std::size_t, std::size_t > edge{ inpoed[e*2],inpoed[e*2+1] };
      auto rx = err.scalar( u3, edge, 0, coord, inpoel, esup, errtype );
      ensure( "edge error < 0.0", rx > -pr );
      ensure( "edge error > 1.0", rx < 1.0+pr );
      auto ry = err.scalar( u3, edge, 1, coord, inpoel, esup, errtype );
      ensure( "edge error < 0.0", ry > -pr );
      ensure( "edge error > 1.0", ry < 1.0+pr );
      auto rz = err.scalar( u3, edge, 2, coord, inpoel, esup, errtype );
      ensure( "edge error < 0.0", rz > -pr );
      ensure( "edge error > 1.0", rz < 1.0+pr );
    }

    // generate a quadratic scalar field with a slope in only z direction
    tk::Fields u2z( npoin, 1 );
    for (std::size_t p=0; p<npoin; ++p) u2z(p,0,0) = coord[2][p]*coord[2][p];
    // test jump error indicator on all edges
    for (std::size_t e=0; e<inpoed.size()/2; ++e) {
      std::pair< std::size_t, std::size_t > edge{ inpoed[e*2],inpoed[e*2+1] };
      auto r = err.scalar( u2z, edge, 0, coord, inpoel, esup, errtype );
      ensure( "edge error < 0.0", r > -pr );
      ensure( "edge error > 1.0", r < 1.0+pr );
    }
  }
};

//! Test group shortcuts
using AMRError_group = test_group< AMRError_common, MAX_TESTS_IN_GROUP >;
using AMRError_object = AMRError_group::object;

//! Define test group
static AMRError_group AMRError( "Inciter/AMR/Error" );

//! Test definitions for group

//! Test jump error indicator for tetrahedron mesh
template<> template<>
void AMRError_object::test< 1 >() {
  set_test_name( "jump indicator on scalar" );
  TestErrorIndicator( inciter::ctr::AMRErrorType::JUMP );
}

//! Test jump error indicator for tetrahedron mesh
template<> template<>
void AMRError_object::test< 2 >() {
  set_test_name( "Hessian indicator on scalar" );
  TestErrorIndicator( inciter::ctr::AMRErrorType::HESSIAN );
}

} // tut::

#endif // test_AMRError_h
