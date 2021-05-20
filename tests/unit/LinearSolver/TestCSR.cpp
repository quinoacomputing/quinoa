// *****************************************************************************
/*!
  \file      tests/unit/LinearSolver/TestCSR.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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
#include "Vector.hpp"

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

  // Mesh node coordinates for simple tet mesh above
  std::array< std::vector< tk::real >, 3 > coord {{
    {{ 0, 1, 1, 0, 0, 1, 1, 0, 0.5, 0.5, 0.5, 1,   0.5, 0 }},
    {{ 0, 0, 1, 1, 0, 0, 1, 1, 0.5, 0.5, 0,   0.5, 1,   0.5 }},
    {{ 0, 0, 0, 0, 1, 1, 1, 1, 0,   1,   0.5, 0.5, 0.5, 0.5 }} }};
};

//! Test group shortcuts
using CSR_group =
  test_group< CSR_common, MAX_TESTS_IN_GROUP >;
using CSR_object = CSR_group::object;

//! Define test group
static CSR_group CSR( "LinearSolver/CSR" );

//! Test definitions for group

//! Test if constructor does not throw on positive ncomp and non-empty psup
template<> template<>
void CSR_object::test< 1 >() {
  set_test_name( "ctor doesn't throw on positive ncomp & valid psup" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  // Generate points surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );

  tk::CSR c( 1, psup );
}

//! Test if constructor does not throw on positive size
template<> template<>
void CSR_object::test< 2 >() {
  set_test_name( "ctor throws on zero ncomp" );

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
            std::string( e.what() ).find( "ncomp must be positive" ) !=
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

//! Test write_stored()
template<> template<>
void CSR_object::test< 4 >() {
  set_test_name("write_stored" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  // Generate points surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );

  tk::CSR c( 1, psup );

  std::stringstream ss;
  c.write_stored( ss );

  auto correct = R"(size (npoin) = 14
ncomp = 1
rsize (size*ncomp) = 14
nnz = 112
rnz[npoin=14] = { 7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 9, 10, 9, 10 }
ia[rsize+1=15] = { 1, 8, 15, 22, 29, 36, 43, 50, 57, 66, 75, 84, 94, 103, 113 }
ja[nnz=112] = { 1, 2, 4, 5, 9, 11, 14, 1, 2, 3, 6, 9, 11, 12, 2, 3, 4, 7, 9, 12, 13, 1, 3, 4, 8, 9, 13, 14, 1, 5, 6, 8, 10, 11, 14, 2, 5, 6, 7, 10, 11, 12, 3, 6, 7, 8, 10, 12, 13, 4, 5, 7, 8, 10, 13, 14, 1, 2, 3, 4, 9, 11, 12, 13, 14, 5, 6, 7, 8, 10, 11, 12, 13, 14, 1, 2, 5, 6, 9, 10, 11, 12, 14, 2, 3, 6, 7, 9, 10, 11, 12, 13, 14, 3, 4, 7, 8, 9, 10, 12, 13, 14, 1, 4, 5, 8, 9, 10, 11, 12, 13, 14 }
a[nnz=112] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
)";

  ensure_equals( "CSR write_stored incorrect", ss.str(), correct );
}

//! Test write_structure()
template<> template<>
void CSR_object::test< 5 >() {
  set_test_name("write_structure" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  // Generate points surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );

  tk::CSR c( 1, psup );

  std::stringstream ss;
  c.write_structure( ss );

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

  ensure_equals( "CSR write_structure incorrect", ss.str(), correct );
}

//! Test write_matrix()
template<> template<>
void CSR_object::test< 6 >() {
  set_test_name("write_matrix" );

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
  c.write_matrix( ss );
  auto s = ss.str();

  // remove white space for easier comparison
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());

  auto correct = "4000000000000004000000000000004000000000000004000000000000004000000000000004000000000000004000000000000004000000000000001000000000000000100000000000000010000000000000001200000000000000100000000000000012";
  ensure_equals( "CSR write_matrix incorrect", s, correct );
}

//! Test write_matlab()
template<> template<>
void CSR_object::test< 7 >() {
  set_test_name("write_matlab" );

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
  c.write_matlab( ss );

  auto s = ss.str();

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

  ensure_equals( "CSR write_matlab incorrect", ss.str(), correct );
}

//! Test rsize
template<> template<>
void CSR_object::test< 8 >() {
  set_test_name( "rsize" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  // Generate points surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );

  tk::CSR c( 1, psup );
  ensure_equals( "CSR::rsize (ncomp=1) incorrect", c.rsize(),
                 psup.second.size()-1 );
  tk::CSR d( 3, psup );
  ensure_equals( "CSR::rsize (ncomp=3) incorrect", d.rsize(),
                 (psup.second.size()-1)*3 );
}

//! Test matrix-vector multiply
template<> template<>
void CSR_object::test< 9 >() {
  set_test_name( "matrix-vector multiply" );

  // Shift node IDs to start from zero
  tk::shiftToZero( inpoel );
  // Generate points surrounding points
  auto psup = tk::genPsup( inpoel, 4, tk::genEsup(inpoel,4) );
  // Query number of nodes in mesh
  auto npoin = psup.second.size()-1;

  tk::CSR A( 1, psup );

  const auto& X = coord[0];
  const auto& Y = coord[1];
  const auto& Z = coord[2];

  // fill matrix with Laplacian
  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    // access node IDs
    const std::array< std::size_t, 4 >
      N{{ inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] }};
    // compute element Jacobi determinant
    const std::array< tk::real, 3 >
      ba{{ X[N[1]]-X[N[0]], Y[N[1]]-Y[N[0]], Z[N[1]]-Z[N[0]] }},
      ca{{ X[N[2]]-X[N[0]], Y[N[2]]-Y[N[0]], Z[N[2]]-Z[N[0]] }},
      da{{ X[N[3]]-X[N[0]], Y[N[3]]-Y[N[0]], Z[N[3]]-Z[N[0]] }};
    const auto J = tk::triple( ba, ca, da );        // J = 6V
    Assert( J > 0, "Element Jacobian non-positive" );

    // shape function derivatives, nnode*ndim [4][3]
    std::array< std::array< tk::real, 3 >, 4 > grad;
    grad[1] = tk::crossdiv( ca, da, J );
    grad[2] = tk::crossdiv( da, ba, J );
    grad[3] = tk::crossdiv( ba, ca, J );
    for (std::size_t i=0; i<3; ++i)
      grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

    for (std::size_t a=0; a<4; ++a)
      for (std::size_t k=0; k<3; ++k)
         for (std::size_t b=0; b<4; ++b)
           A(N[a],N[b]) += J/6 * grad[a][k] * grad[b][k];
  }

  std::vector< tk::real > x( npoin );
  std::iota( begin(x), end(x), 0.0 );
  auto r = x;

  A.mult( x, r );

  std::vector< tk::real > correct{
    -6.124999999999999, -5.25, -5.041666666666666, -5, -4.166666666666666,
    -3.291666666666666, -3.083333333333333, -3.041666666666666,
     2.916666666666672, 1.083333333333338, 5.500000000000007,
    7.833333333333343,6.833333333333341,10.83333333333334 };

  tk::real prec = std::numeric_limits< tk::real >::epsilon()*100;

  for (std::size_t i=0; i<r.size(); ++i)
    ensure_equals( "incorrect matrix-vector product",
                   r[i], correct[i], prec );
}

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
