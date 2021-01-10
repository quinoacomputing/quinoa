// *****************************************************************************
/*!
  \file      tests/unit/LinearSolver/TestCSR.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for LinearSolver/CSR
  \details   Unit tests for LinearSolver/CSR
*/
// *****************************************************************************

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "CSR.hpp"
#include "Reorder.hpp"
#include "DerivedData.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif

//! All tests in group inherited from this base
struct CSR_common {

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
};

//! Test group shortcuts
using CSR_group =
  test_group< CSR_common, MAX_TESTS_IN_GROUP >;
using CSR_object = CSR_group::object;

//! Define test group
static CSR_group CSR( "LinearSolver/CSR" );

//! Test definitions for group

//! Test if constructor does not throw on positive dof and non-empty psup
template<> template<>
void CSR_object::test< 1 >() {
  set_test_name( "ctor doesn't throw on positive DOF and valid psup" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  // Generate points surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );

  tk::CSR c( 1, psup );
}

//! Test if constructor does not throw on positive size
template<> template<>
void CSR_object::test< 2 >() {
  set_test_name( "ctor throws on zero DOF" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  // Generate points surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );

  // exception only thrown in DEBUG mode, but no problem in RELEASE, only
  // exception is not thrown
  try {

    tk::CSR c( 0, psup );

  } catch( tk::Exception& e ) {
    // exception thrown, test ok
    // if any other type of exception is thrown, test fails with except
    // find out if exception was thrown due to the correct reason
    ensure( std::string("wrong exception thrown: ") + e.what(),
            std::string( e.what() ).find( "DOF must be positive" ) !=
              std::string::npos );
  }
}

//! Test if constructor throws on empty psup
template<> template<>
void CSR_object::test< 3 >() {
  set_test_name( "ctor throws on empty psup" );

  // exception thrown in both DEBUG and RELEASE
  try {

    tk::CSR c( 1, {} );

  } catch( tk::Exception& ) {
    // exception thrown, test ok
  }
}

//! Test write_as_stored()
template<> template<>
void CSR_object::test< 4 >() {
  set_test_name("write_as_stored" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  // Generate points surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );

  tk::CSR c( 1, psup );

  std::stringstream ss;
  c.write_as_stored( ss );

  auto correct = R"(size (npoin) = 14
dof = 1
rsize (size*dof) = 14
nnz = 112
rnz[npoin=14] = { 7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 9, 10, 9, 10 }
ia[rsize+1=15] = { 1, 8, 15, 22, 29, 36, 43, 50, 57, 66, 75, 84, 94, 103, 113 }
ja[nnz=112] = { 1, 2, 4, 5, 9, 11, 14, 1, 2, 3, 6, 9, 11, 12, 2, 3, 4, 7, 9, 12, 13, 1, 3, 4, 8, 9, 13, 14, 1, 5, 6, 8, 10, 11, 14, 2, 5, 6, 7, 10, 11, 12, 3, 6, 7, 8, 10, 12, 13, 4, 5, 7, 8, 10, 13, 14, 1, 2, 3, 4, 9, 11, 12, 13, 14, 5, 6, 7, 8, 10, 11, 12, 13, 14, 1, 2, 5, 6, 9, 10, 11, 12, 14, 2, 3, 6, 7, 9, 10, 11, 12, 13, 14, 3, 4, 7, 8, 9, 10, 12, 13, 14, 1, 4, 5, 8, 9, 10, 11, 12, 13, 14 }
a[nnz=112] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
)";

  ensure_equals( "CSR write_as_stored incorrect", ss.str(), correct );
}

//! Test write_as_structure()
template<> template<>
void CSR_object::test< 5 >() {
  set_test_name("write_as_structure" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  // Generate points surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );

  tk::CSR c( 1, psup );

  std::stringstream ss;
  c.write_as_structure( ss );

  auto correct = R"(o o . o o . . . o . o . . o 
o o o . . o . . o . o o . . 
. o o o . . o . o . . o o . 
o . o o . . . o o . . . o o 
o . . . o o . o . o o . . o 
. o . . o o o . . o o o . . 
. . o . . o o o . o . o o . 
. . . o o . o o . o . . o o 
o o o o . . . . o . o o o o 
. . . . o o o o . o o o o o 
o o . . o o . . o o o o . o 
. o o . . o o . o o o o o o 
. . o o . . o o o o . o o o 
o . . o o . . o o o o o o o 
)";

  ensure_equals( "CSR write_as_structure incorrect", ss.str(), correct );
}

//! Test write_as_matrix()
template<> template<>
void CSR_object::test< 6 >() {
  set_test_name("write_as_matrix" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  // Generate points surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );

  tk::CSR c( 1, psup );

  // fill in some matrix entries
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    auto A = inpoel[e*4+0];
    auto B = inpoel[e*4+1];
    auto C = inpoel[e*4+2];
    auto D = inpoel[e*4+3];
    c(A,A) += 1.0;
    c(B,B) += 1.0;
    c(C,C) += 1.0;
    c(D,D) += 1.0;
  }

  std::stringstream ss;
  c.write_as_matrix( ss );
  auto s = ss.str();

  // remove white space for easier comparison
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());

  auto correct = "4000000000000004000000000000004000000000000004000000000000004000000000000004000000000000004000000000000004000000000000001000000000000000100000000000000010000000000000001200000000000000100000000000000012";
  ensure_equals( "CSR write_as_matrix incorrect", s, correct );
}

//! Test write_as_matlab()
template<> template<>
void CSR_object::test< 7 >() {
  set_test_name("write_as_matlab" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  // Generate points surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );

  tk::CSR c( 1, psup );

  // fill in some matrix entries
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    auto A = inpoel[e*4+0];
    auto B = inpoel[e*4+1];
    auto C = inpoel[e*4+2];
    auto D = inpoel[e*4+3];
    c(A,A) += 1.0;
    c(B,B) += 1.0;
    c(C,C) += 1.0;
    c(D,D) += 1.0;
  }

  std::stringstream ss;
  c.write_as_matlab( ss );
  //std::cout << ss.str();

  auto s = ss.str();

  // remove white space for easier comparison
  //s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());

  auto correct = R"(A = [ 4 0 0 0 0 0 0 0 0 0 0 0 0 0 ;
0 4 0 0 0 0 0 0 0 0 0 0 0 0 ;
0 0 4 0 0 0 0 0 0 0 0 0 0 0 ;
0 0 0 4 0 0 0 0 0 0 0 0 0 0 ;
0 0 0 0 4 0 0 0 0 0 0 0 0 0 ;
0 0 0 0 0 4 0 0 0 0 0 0 0 0 ;
0 0 0 0 0 0 4 0 0 0 0 0 0 0 ;
0 0 0 0 0 0 0 4 0 0 0 0 0 0 ;
0 0 0 0 0 0 0 0 10 0 0 0 0 0 ;
0 0 0 0 0 0 0 0 0 10 0 0 0 0 ;
0 0 0 0 0 0 0 0 0 0 10 0 0 0 ;
0 0 0 0 0 0 0 0 0 0 0 12 0 0 ;
0 0 0 0 0 0 0 0 0 0 0 0 10 0 ;
0 0 0 0 0 0 0 0 0 0 0 0 0 12 ;
]
)";

  ensure_equals( "CSR write_as_matlab incorrect", ss.str(), correct );
}

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT