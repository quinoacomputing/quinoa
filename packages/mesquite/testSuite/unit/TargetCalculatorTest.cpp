/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TargetCalculatorTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "UnitUtil.hpp"
#include "Mesquite_TargetCalculator.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_PatchData.hpp"
#include "Mesquite_IdealElements.hpp"
#include "Mesquite_LinearHexahedron.hpp"
#include "Mesquite_LinearTriangle.hpp"
#include "Mesquite_ReferenceMesh.hpp"
#include <iostream>

using namespace Mesquite;

double EPS = 0.05;
double EPSBIG = 0.5;

class TargetCalculatorTest : public CppUnit::TestFixture
{
private:

  CPPUNIT_TEST_SUITE(TargetCalculatorTest);

  CPPUNIT_TEST (test_factor_2D);
  CPPUNIT_TEST (test_factor_surface);
  CPPUNIT_TEST (test_factor_3D);
  CPPUNIT_TEST (test_factor_2D_zero);
  CPPUNIT_TEST (test_factor_surface_zero);
  CPPUNIT_TEST (test_factor_3D_zero);
  
  CPPUNIT_TEST (test_size_2D);
  CPPUNIT_TEST (test_size_surface);
  CPPUNIT_TEST (test_size_3D);
  
  CPPUNIT_TEST (test_skew_2D);
  CPPUNIT_TEST (test_skew_surface);
  CPPUNIT_TEST (test_skew_3D);
  
  CPPUNIT_TEST (test_aspect_2D);
  CPPUNIT_TEST (test_aspect_surface);
  CPPUNIT_TEST (test_aspect_3D);
  
  CPPUNIT_TEST (test_shape_2D);
  CPPUNIT_TEST (test_shape_surface);
  CPPUNIT_TEST (test_shape_3D);
  
  CPPUNIT_TEST (test_ideal_skew_tri);
  CPPUNIT_TEST (test_ideal_skew_quad);
  CPPUNIT_TEST (test_ideal_skew_tet);
  CPPUNIT_TEST (test_ideal_skew_prism);
  CPPUNIT_TEST (test_ideal_skew_pyramid);
  CPPUNIT_TEST (test_ideal_skew_hex);
  
  CPPUNIT_TEST (test_ideal_shape_tri);
  CPPUNIT_TEST (test_ideal_shape_quad);
  CPPUNIT_TEST (test_ideal_shape_tet);
  CPPUNIT_TEST (test_ideal_shape_prism);
  CPPUNIT_TEST (test_ideal_shape_pyramid);
  CPPUNIT_TEST (test_ideal_shape_hex);
  
  CPPUNIT_TEST (test_new_orientatin_3D);
  CPPUNIT_TEST (test_new_orientatin_2D);
  
  CPPUNIT_TEST (test_new_aspect_3D);
  CPPUNIT_TEST (test_new_aspect_2D);
  
  CPPUNIT_TEST (test_jacobian_3D);
  CPPUNIT_TEST (test_jacobian_2D);
  
  CPPUNIT_TEST (test_get_refmesh_Jacobian_3D);
  CPPUNIT_TEST (test_get_refmesh_Jacobian_2D);

  CPPUNIT_TEST_SUITE_END();
  
    // define some matrices for use in testing.
  MsqMatrix<3,3> V3D_Z45, V3D_X90, Q3D_45, D3D_123;
  MsqMatrix<3,2> V2D_Z45, V2D_X90;
  MsqMatrix<2,2> Q2D_45, D2D_21;
  
public:

  TargetCalculatorTest();
  void setUp();

  void test_factor_2D();
  void test_factor_surface();
  void test_factor_3D();
  void test_factor_2D_zero();
  void test_factor_surface_zero();
  void test_factor_3D_zero();
  
  void test_size_2D();
  void test_size_surface();
  void test_size_3D();
  
  void test_skew_2D();
  void test_skew_surface();
  void test_skew_3D();
  
  void test_aspect_2D();
  void test_aspect_surface();
  void test_aspect_3D();
  
  void test_shape_2D();
  void test_shape_surface();
  void test_shape_3D();
  
  void test_ideal_skew_tri();
  void test_ideal_skew_quad();
  void test_ideal_skew_tet();
  void test_ideal_skew_prism();
  void test_ideal_skew_pyramid();
  void test_ideal_skew_hex();
  
  void test_ideal_shape_tri();
  void test_ideal_shape_quad();
  void test_ideal_shape_tet();
  void test_ideal_shape_prism();
  void test_ideal_shape_pyramid();
  void test_ideal_shape_hex();
  
  void test_new_orientatin_3D();
  void test_new_orientatin_2D();
  
  void test_new_aspect_3D();
  void test_new_aspect_2D();
  
  void test_jacobian_3D();
  void test_jacobian_2D();
  
  void test_get_refmesh_Jacobian_3D();
  void test_get_refmesh_Jacobian_2D();


  template <unsigned R, unsigned C> static inline
  void check_valid_V( MsqMatrix<R,C> V );

  template <unsigned D> static inline
  void check_valid_Q( MsqMatrix<D,D> Q );

  template <unsigned D> static inline
  void check_valid_delta( MsqMatrix<D,D> delta );
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TargetCalculatorTest, "TargetCalculatorTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TargetCalculatorTest, "Unit");

TargetCalculatorTest::TargetCalculatorTest()
{
  const double cos45 = MSQ_SQRT_TWO/2.0;
  const double rotation_3D_Z45[9] = { cos45, -cos45, 0,
                                      cos45,  cos45, 0,
                                          0,      0, 1 };
  const double rotation_2D_Z45[6] = { cos45, -cos45,
                                      cos45,  cos45,
                                          0,      0 };
                                          
  const double rotation_3D_X90[9] = { 1,  0,  0,
                                      0,  0, -1,
                                      0,  1,  0 };
  const double rotation_2D_X90[6] = { 1,  0,
                                      0,  0,
                                      0,  1 };
                                      
  const double rc45 = sqrt(cos45);
  const double skew_2D_45[4] = { 1/rc45, rc45,
                                 0,      rc45 };
  const double skew_3D_45[9] = { 1, cos45, cos45,
                                 0, cos45, 1 - cos45,
                                 0,     0, sqrt(MSQ_SQRT_TWO-1) };
                                 
  const double aspect_2D_2x[4] = { MSQ_SQRT_TWO, 0,
                                   0,            MSQ_SQRT_TWO/2 };
  const double r6 = Mesquite::cbrt(1.0/6.0);
  const double aspect_3D_123[9] = { r6,    0,    0,
                                     0, 2*r6,    0,
                                     0,    0, 3*r6 };
                                     
  V3D_Z45 = MsqMatrix<3,3>(rotation_3D_Z45);
  V3D_X90 = MsqMatrix<3,3>(rotation_3D_X90);
  Q3D_45  = MsqMatrix<3,3>(skew_3D_45);
  Q3D_45  *= 1/Mesquite::cbrt(det(Q3D_45));
  D3D_123 = MsqMatrix<3,3>(aspect_3D_123);
  
  V2D_Z45 = MsqMatrix<3,2>(rotation_2D_Z45);
  V2D_X90 = MsqMatrix<3,2>(rotation_2D_X90);
  Q2D_45  = MsqMatrix<2,2>(skew_2D_45);
  D2D_21  = MsqMatrix<2,2>(aspect_2D_2x);
}

template <unsigned R, unsigned C> inline
void TargetCalculatorTest::check_valid_V( MsqMatrix<R,C> V )
{
  // check that it is a rotation
  MsqMatrix<C,C> I(1.0);
  ASSERT_MATRICES_EQUAL( I, transpose(V) * V, EPS );
}

template <unsigned D> inline
void TargetCalculatorTest::check_valid_Q( MsqMatrix<D,D> Q )
{
  // must have unit determinant
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(Q), EPS );
  
  // must be upper-triangular
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, Q(1,0), EPS );
  if (D == 3) {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, Q(2,0), EPS );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, Q(2,1), EPS );
  }
  
  // columns must be of equal length
  CPPUNIT_ASSERT_DOUBLES_EQUAL( length(Q.column(0)), length(Q.column(1)), EPSBIG );
  if (D == 3) {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( length(Q.column(0)), length(Q.column(2)), EPSBIG );
  }
  
  // diagonal elements must be greater than zero
  CPPUNIT_ASSERT( Q(0,0) - EPS >= 0.0 );
  CPPUNIT_ASSERT( Q(1,1) - EPS >= 0.0 );
  if (D == 3) {
    CPPUNIT_ASSERT( Q(2,2) - EPS >= 0.0 );
  }
}

template <unsigned D> inline
void TargetCalculatorTest::check_valid_delta( MsqMatrix<D,D> delta )
{
  // must have unit determinant
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(delta), EPS );
  
  // must be diagonal
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, delta(1,0), EPS );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, delta(0,1), EPS );
  if (D == 3) {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, delta(2,0), EPS );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, delta(0,2), EPS );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, delta(1,2), EPS );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, delta(2,1), EPS );
  }
}

void TargetCalculatorTest::setUp()
{
    // test the test: make sure we're testing with 
    // valid matrix factors.
    
  check_valid_V( V3D_Z45 );
  check_valid_V( V3D_X90 );
  check_valid_V( V2D_Z45 );
  check_valid_V( V2D_X90 );
  
  check_valid_Q( Q3D_45 );
  check_valid_Q( Q2D_45 );
  
  check_valid_delta( D3D_123 );
  check_valid_delta( D2D_21 );
}

void TargetCalculatorTest::test_factor_2D()
{
  MsqPrintError err(std::cout);
  MsqMatrix<2,2> I(1.0), W, V;
  double lambda;
  MsqMatrix<2,2> Q, delta;
  bool valid;
  
    // first test with I
  valid = TargetCalculator::factor_2D( I, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( I, W, EPS );
  
    // now test with 2*I
  valid = TargetCalculator::factor_2D( 2*I, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( 2*I, W, EPS );
  
    // now rotate the matrix about the X axis by 90 degrees
  MsqMatrix<2,2> I_rot(0.0);
  I_rot(0,1) = 1.0;
  I_rot(1,0) = -1.0;
  valid = TargetCalculator::factor_2D( I_rot, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( I_rot, W, EPS );
  
    // change the aspect ratio
  MsqMatrix<2,2> A(1.0);
  A(1,1) = 1.5;
  valid = TargetCalculator::factor_2D( A, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( pow(det(transpose(A)*A),0.25), lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( A, W, EPS );
  
  // try an arbitrary matrix
  double w[4] = { 3, 1, 
                  4, 2 };
  MsqMatrix<2,2> W2(w);
  valid = TargetCalculator::factor_2D( W2, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( pow(det(transpose(W2)*W2),0.25), lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( W2, W, EPS );
  
  // try a more elaborate test
  const double e = exp((double)1);
  MsqMatrix<2,2> Z45(V2D_Z45.data()); // copy first two rows
  W2 = e * Z45 * Q2D_45 * D2D_21;
  valid = TargetCalculator::factor_2D( W2, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( e, lambda, EPS );
  ASSERT_MATRICES_EQUAL( Z45, V, EPS );
  ASSERT_MATRICES_EQUAL( Q2D_45, Q, EPS );
  ASSERT_MATRICES_EQUAL( D2D_21, delta, EPS );
}

void TargetCalculatorTest::test_factor_surface()
{
  MsqPrintError err(std::cout);
  MsqMatrix<3,2> I(1.0), W, V;
  double lambda;
  MsqMatrix<2,2> Q, delta;
  bool valid;
  
    // first test with I
  valid = TargetCalculator::factor_surface( I, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( I, W, EPS );
  
    // now test with 2*I
  valid = TargetCalculator::factor_surface( 2*I, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( 2*I, W, EPS );
  
    // now rotate the matrix about the X axis by 90 degrees
  MsqMatrix<3,2> I_rot(1.0);
  I_rot(1,1) = 0.0; I_rot(2,1) = 1.0;
  valid = TargetCalculator::factor_surface( I_rot, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( I_rot, W, EPS );
  
    // change the aspect ratio
  MsqMatrix<3,2> A(1.0);
  A(1,1) = 1.5;
  valid = TargetCalculator::factor_surface( A, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( pow(det(transpose(A)*A),0.25), lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( A, W, EPS );
  
  // try an arbitrary matrix
  double w[6] = { 3, 1, 
                  4, 2,
                 -1, 5 };
  MsqMatrix<3,2> W2(w);
  valid = TargetCalculator::factor_surface( W2, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( pow(det(transpose(W2)*W2),0.25), lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( W2, W, EPS );
  
  // try a more elaborate test
  const double e = exp((double)1);
  W2 = e * V2D_Z45 * Q2D_45 * D2D_21;
  valid = TargetCalculator::factor_surface( W2, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( e, lambda, EPS );
  ASSERT_MATRICES_EQUAL( V2D_Z45, V, EPS );
  ASSERT_MATRICES_EQUAL( Q2D_45, Q, EPS );
  ASSERT_MATRICES_EQUAL( D2D_21, delta, EPS );
}

void TargetCalculatorTest::test_factor_3D()
{
  MsqPrintError err(std::cout);
  MsqMatrix<3,3> I(1.0), W, V, Q, delta;
  double lambda;
  bool valid;
  
    // first test with I
  valid = TargetCalculator::factor_3D( I, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( I, W, EPS );
  
    // now test with 2*I
  valid = TargetCalculator::factor_3D( 2*I, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( 2*I, W, EPS );
  
    // now rotate the matrix about the X axis by 90 degrees
  MsqMatrix<3,3> I_rot(1.0);
  I_rot(1,1) = 0.0;  I_rot(2,1) = 1.0;
  I_rot(1,2) = -1.0; I_rot(2,2) = 0.0;
  valid = TargetCalculator::factor_3D( I_rot, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( I_rot, W, EPS );
  
    // change the aspect ratio
  MsqMatrix<3,3> A(1.0);
  A(1,1) = 1.5;
  valid = TargetCalculator::factor_3D( A, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( Mesquite::cbrt(det(A)), lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( A, W, EPS );
  
  // try an arbitrary matrix
  double vals[] = { 9, 8, 7, 
                    1, 5, 4, 
                    3, 2, 6 };
  MsqMatrix<3,3> W2(vals);
  valid = TargetCalculator::factor_3D( W2, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( Mesquite::cbrt(det(W2)), lambda, EPS );
  check_valid_V( V );
  check_valid_Q( Q );
  check_valid_delta( delta );
  W = lambda * V * Q* delta;
  ASSERT_MATRICES_EQUAL( W2, W, EPSBIG );
  
  // try a more elaborate test
  const double e = exp((double)1);
  W2 = e * V3D_Z45 * Q3D_45 * D3D_123;
  valid = TargetCalculator::factor_3D( W2, lambda, V, Q, delta, err );
  CPPUNIT_ASSERT( valid && !err );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( e, lambda, EPS );
  ASSERT_MATRICES_EQUAL( V3D_Z45, V, EPS );
  ASSERT_MATRICES_EQUAL( Q3D_45, Q, EPS );
  ASSERT_MATRICES_EQUAL( D3D_123, delta, EPS );
}

void TargetCalculatorTest::test_factor_2D_zero()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<2,2> Q, delta;
  MsqMatrix<2,2> V, Z(0.0); Z(0,0) = 1.0;
  
  bool valid = TargetCalculator::factor_2D( Z, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!valid);
}

void TargetCalculatorTest::test_factor_surface_zero()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<2,2> Q, delta;
  MsqMatrix<3,2> V, Z(0.0); Z(0,0) = 1.0;
  
  bool valid = TargetCalculator::factor_surface( Z, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!valid);
}

void TargetCalculatorTest::test_factor_3D_zero()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<3,3> V, Q, delta, Z(0.0); Z(0,0) = 1.0;
  
  bool valid = TargetCalculator::factor_3D( Z, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!valid);
}


void TargetCalculatorTest::test_size_2D()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<2,2> V, Q, delta;
  
  MsqMatrix<2,2> I(1.0);
  TargetCalculator::factor_2D( I, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( lambda, TargetCalculator::size(I), EPS );

  double w[] = { 3, 1, 
                 4, 2 };
  MsqMatrix<2,2> W(w);
  TargetCalculator::factor_2D( W, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( lambda, TargetCalculator::size(W), EPS );
}

void TargetCalculatorTest::test_size_surface()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<3,2> V;
  MsqMatrix<2,2> Q, delta;
  
  MsqMatrix<3,2> I(1.0);
  TargetCalculator::factor_surface( I, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( lambda, TargetCalculator::size(I), EPS );

  double w[] = { 3, 1, 
                 4, 2, 
                -1, 5 };
  MsqMatrix<3,2> W(w);
  TargetCalculator::factor_surface( W, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( lambda, TargetCalculator::size(W), EPS );
}

void TargetCalculatorTest::test_size_3D()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<3,3> V, Q, delta;
  
  MsqMatrix<3,3> I(1.0);
  TargetCalculator::factor_3D( I, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( lambda, TargetCalculator::size(I), EPS );

  double w[] = { 9, 8, 7, 
                 1, 5, 4, 
                 3, 2, 6 };
  MsqMatrix<3,3> W(w);
  TargetCalculator::factor_3D( W, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( lambda, TargetCalculator::size(W), EPS );
}
  
void TargetCalculatorTest::test_skew_2D()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<2,2> V, Q, delta;
  
  MsqMatrix<2,2> I(1.0);
  TargetCalculator::factor_2D( I, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q, TargetCalculator::skew(I), EPS );
  
  double r[] = { 0, -1,  
                 2, 0 };
  MsqMatrix<2,2> R(r);
  TargetCalculator::factor_2D( R, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q, TargetCalculator::skew(R), EPS );

  double w[] = { 3, 1, 
                 4, 2 };
  MsqMatrix<2,2> W(w);
  TargetCalculator::factor_2D( W, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q, TargetCalculator::skew(W), EPS );
  
  W = MsqMatrix<2,2>(6) * Q2D_45;
  ASSERT_MATRICES_EQUAL( Q2D_45, TargetCalculator::skew(W), EPS );
}
  
void TargetCalculatorTest::test_skew_surface()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<3,2> V;
  MsqMatrix<2,2> Q, delta;
  
  MsqMatrix<3,2> I(1.0);
  TargetCalculator::factor_surface( I, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q, TargetCalculator::skew(I), EPS );
  
  double r[] = { 1, 0,  
                 0, 0,  
                 0, 2 };
  MsqMatrix<3,2> R(r);
  TargetCalculator::factor_surface( R, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q, TargetCalculator::skew(R), EPS );

  double w[] = { 3, 1, 
                 4, 2, 
                -1, 5 };
  MsqMatrix<3,2> W(w);
  TargetCalculator::factor_surface( W, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q, TargetCalculator::skew(W), EPS );
  
  W = MsqMatrix<3,2>(6) * Q2D_45;
  ASSERT_MATRICES_EQUAL( Q2D_45, TargetCalculator::skew(W), EPS );
}

void TargetCalculatorTest::test_skew_3D()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<3,3> V, Q, delta;
  
  MsqMatrix<3,3> I(1.0);
  TargetCalculator::factor_3D( I, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q, TargetCalculator::skew(I), EPS );
  
  double r[] = { 0.8, 0.0, 0.0,
                 0.0, 0.0, 1.13,
                 0.0,-0.5, 0.0 };
  MsqMatrix<3,3> R(r);
  TargetCalculator::factor_3D( R, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q, TargetCalculator::skew(R), EPS );

  double w[] = { 9, 8, 7, 
                 1, 5, 4, 
                 3, 2, 6 };
  MsqMatrix<3,3> W(w);
  TargetCalculator::factor_3D( W, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q, TargetCalculator::skew(W), 2*EPS );
  
  W = MsqMatrix<3,3>(2.1) * Q3D_45;
  ASSERT_MATRICES_EQUAL( Q3D_45, TargetCalculator::skew(W), EPS );
}
  
void TargetCalculatorTest::test_aspect_2D()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<2,2> V, Q, delta;
  
  MsqMatrix<2,2> I(1.0);
  TargetCalculator::factor_2D( I, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( delta, TargetCalculator::aspect(I), EPS );
  
  double r[] = { 0, -1,  
                 2, 0 };
  MsqMatrix<2,2> R(r);
  TargetCalculator::factor_2D( R, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( delta, TargetCalculator::aspect(R), EPS );

  double w[] = { 3, 1, 
                 4, 2 };
  MsqMatrix<2,2> W(w);
  TargetCalculator::factor_2D( W, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( delta, TargetCalculator::aspect(W), EPS );
  
  W = MsqMatrix<2,2>(6) * Q2D_45;
  ASSERT_MATRICES_EQUAL( I, TargetCalculator::aspect(W), EPS );
}
  
void TargetCalculatorTest::test_aspect_surface()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<3,2> V;
  MsqMatrix<2,2> Q, delta;
  
  MsqMatrix<3,2> I(1.0);
  TargetCalculator::factor_surface( I, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( delta, TargetCalculator::aspect(I), EPS );
  
  double r[] = { 1, 0,  
                 0, 0,  
                 0, 2 };
  MsqMatrix<3,2> R(r);
  TargetCalculator::factor_surface( R, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( delta, TargetCalculator::aspect(R), EPS );

  double w[] = { 3, 1, 
                 4, 2, 
                -1, 5 };
  MsqMatrix<3,2> W(w);
  TargetCalculator::factor_surface( W, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( delta, TargetCalculator::aspect(W), EPS );
  
  W = MsqMatrix<3,2>(6) * Q2D_45;
  ASSERT_MATRICES_EQUAL( (MsqMatrix<2,2>(1.0)), TargetCalculator::aspect(W), EPS );
}

void TargetCalculatorTest::test_aspect_3D()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<3,3> V, Q, delta;
  
  MsqMatrix<3,3> I(1.0);
  TargetCalculator::factor_3D( I, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( delta, TargetCalculator::aspect(I), EPS );
  
  double r[] = { 0.8, 0.0, 0.0,
                 0.0, 0.0, 1.13,
                 0.0,-0.5, 0.0 };
  MsqMatrix<3,3> R(r);
  TargetCalculator::factor_3D( R, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( delta, TargetCalculator::aspect(R), EPS );

  double w[] = { 9, 8, 7, 
                 1, 5, 4, 
                 3, 2, 6 };
  MsqMatrix<3,3> W(w);
  TargetCalculator::factor_3D( W, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( delta, TargetCalculator::aspect(W), 2*EPS );
  
  W = MsqMatrix<3,3>(2.1) * Q3D_45;
  ASSERT_MATRICES_EQUAL( I, TargetCalculator::aspect(W), EPS );
}
  
void TargetCalculatorTest::test_shape_2D()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<2,2> V, Q, delta;
  
  MsqMatrix<2,2> I(1.0);
  TargetCalculator::factor_2D( I, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q*delta, TargetCalculator::shape(I), EPS );
  
  double r[] = { 0, -1,  
                 2, 0 };
  MsqMatrix<2,2> R(r);
  TargetCalculator::factor_2D( R, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q*delta, TargetCalculator::shape(R), EPS );

  double w[] = { 3, 1, 
                 4, 2 };
  MsqMatrix<2,2> W(w);
  TargetCalculator::factor_2D( W, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q*delta, TargetCalculator::shape(W), EPS );
  
  W = MsqMatrix<2,2>(6) * Q2D_45;
  ASSERT_MATRICES_EQUAL( Q2D_45, TargetCalculator::shape(W), EPS );
}
  
void TargetCalculatorTest::test_shape_surface()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<3,2> V;
  MsqMatrix<2,2> Q, delta;
  
  MsqMatrix<3,2> I(1.0);
  TargetCalculator::factor_surface( I, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q*delta, TargetCalculator::shape(I), EPS );
  
  double r[] = { 1, 0,  
                 0, 0,  
                 0, 2 };
  MsqMatrix<3,2> R(r);
  TargetCalculator::factor_surface( R, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q*delta, TargetCalculator::shape(R), EPS );

  double w[] = { 3, 1, 
                 4, 2, 
                -1, 5 };
  MsqMatrix<3,2> W(w);
  TargetCalculator::factor_surface( W, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q*delta, TargetCalculator::shape(W), EPS );
  
  W = MsqMatrix<3,2>(6) * Q2D_45;
  ASSERT_MATRICES_EQUAL( Q2D_45, TargetCalculator::shape(W), EPS );
}

void TargetCalculatorTest::test_shape_3D()
{
  MsqPrintError err(std::cout);
  double lambda;
  MsqMatrix<3,3> V, Q, delta;
  
  MsqMatrix<3,3> I(1.0);
  TargetCalculator::factor_3D( I, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q*delta, TargetCalculator::shape(I), EPS );
  
  double r[] = { 0.8, 0.0, 0.0,
                 0.0, 0.0, 1.13,
                 0.0,-0.5, 0.0 };
  MsqMatrix<3,3> R(r);
  TargetCalculator::factor_3D( R, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q*delta, TargetCalculator::shape(R), EPS );

  double w[] = { 9, 8, 7, 
                 1, 5, 4, 
                 3, 2, 6 };
  MsqMatrix<3,3> W(w);
  TargetCalculator::factor_3D( W, lambda, V, Q, delta, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q*delta, TargetCalculator::shape(W), 2*EPS );
  
  W = MsqMatrix<3,3>(2.1) * Q3D_45;
  ASSERT_MATRICES_EQUAL( Q3D_45, TargetCalculator::shape(W), EPS );
}

void TargetCalculatorTest::test_ideal_skew_tri()
{
  MsqError err;
  PatchData pd;
  MsqMatrix<2,2> Q, Q2;
  
  TargetCalculator::ideal_skew_2D( TRIANGLE, Sample(0,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  check_valid_Q(Q);
  
  TargetCalculator::ideal_skew_2D( TRIANGLE, Sample(2,0), pd, Q2, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q, Q2, EPS );
    
  const Vector3D* coords = unit_edge_element( TRIANGLE );
  MsqMatrix<3,2> W;
  TargetCalculator::jacobian_2D( pd, TRIANGLE, 3, Sample(0,0), coords, W, err );
  ASSERT_NO_ERROR(err);
  Q = TargetCalculator::skew(W);
  ASSERT_MATRICES_EQUAL( Q, Q2, EPS );
}

void TargetCalculatorTest::test_ideal_skew_quad()
{
  MsqError err;
  PatchData pd;
  MsqMatrix<2,2> Q;
  
  TargetCalculator::ideal_skew_2D( QUADRILATERAL, Sample(0,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  ASSERT_IDENTITY_MATRIX( Q );
  
  TargetCalculator::ideal_skew_2D( QUADRILATERAL, Sample(2,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  ASSERT_IDENTITY_MATRIX( Q );
}

void TargetCalculatorTest::test_ideal_skew_tet()
{
  MsqError err;
  PatchData pd;
  MsqMatrix<3,3> Q, Q2;
  
  TargetCalculator::ideal_skew_3D( TETRAHEDRON, Sample(0,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  check_valid_Q(Q);
  
  TargetCalculator::ideal_skew_3D( TETRAHEDRON, Sample(3,0), pd, Q2, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q, Q2, EPS );
    
  const Vector3D* coords = unit_edge_element( TETRAHEDRON );
  MsqMatrix<3,3> W;
  TargetCalculator::jacobian_3D( pd, TETRAHEDRON, 4, Sample(0,0), coords, W, err );
  ASSERT_NO_ERROR(err);
  Q = TargetCalculator::skew(W);
  ASSERT_MATRICES_EQUAL( Q, Q2, EPS );
}

void TargetCalculatorTest::test_ideal_skew_prism()
{
  Sample points[] = { Sample(0,0), Sample(0,1), Sample(0,2), Sample(0,3),
                      Sample(2,0), Sample(3,0) };
  const int num_pts = sizeof(points)/sizeof(points[0]);
  
  MsqError err;
  PatchData pd;
  MsqMatrix<3,3> Q, W;
  const Vector3D* coords = unit_edge_element( PRISM );
  
  for (int i = 0; i < num_pts; ++i) {
  
    TargetCalculator::ideal_skew_3D( PRISM, points[i], pd, Q, err );
    ASSERT_NO_ERROR(err);
    check_valid_Q(Q);
    
    TargetCalculator::jacobian_3D( pd, PRISM, 6, points[i], coords, W, err );
    ASSERT_NO_ERROR(err);
    
    ASSERT_MATRICES_EQUAL( TargetCalculator::skew(W), Q, EPS );
  }
}

void TargetCalculatorTest::test_ideal_skew_pyramid()
{
  Sample points[] = { Sample(0,0), Sample(0,1), Sample(0,2), Sample(0,3),
                      Sample(2,0), Sample(3,0) };
  const int num_pts = sizeof(points)/sizeof(points[0]);
  
  MsqError err;
  PatchData pd;
  MsqMatrix<3,3> Q, W;
  const Vector3D* coords = unit_edge_element( PYRAMID, true );
  
  for (int i = 0; i < num_pts; ++i) {
  
    TargetCalculator::ideal_skew_3D( PYRAMID, points[i], pd, Q, err );
    ASSERT_NO_ERROR(err);
    check_valid_Q(Q);
    
    TargetCalculator::jacobian_3D( pd, PYRAMID, 5, points[i], coords, W, err );
    ASSERT_NO_ERROR(err);
    
    ASSERT_MATRICES_EQUAL( TargetCalculator::skew(W), Q, EPS );
  }
}

void TargetCalculatorTest::test_ideal_skew_hex()
{
  MsqError err;
  PatchData pd;
  MsqMatrix<3,3> Q;
  
  TargetCalculator::ideal_skew_3D( HEXAHEDRON, Sample(0,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  ASSERT_IDENTITY_MATRIX( Q );
  
  TargetCalculator::ideal_skew_3D( HEXAHEDRON, Sample(3,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  ASSERT_IDENTITY_MATRIX( Q );
}


void TargetCalculatorTest::test_ideal_shape_tri()
{
    // Ideal triangle should have Aspect (i.e. delta) == identity,
    // so shape should be equal to skew.

  MsqError err;
  PatchData pd;
  MsqMatrix<2,2> Q, Q2;
  
  TargetCalculator::ideal_shape_2D( TRIANGLE, Sample(0,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  TargetCalculator::ideal_skew_2D( TRIANGLE, Sample(0,0), pd, Q2, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q2, Q, EPS );
  
  TargetCalculator::ideal_shape_2D( TRIANGLE, Sample(2,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  TargetCalculator::ideal_skew_2D( TRIANGLE, Sample(2,0), pd, Q2, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q2, Q, EPS );
}

void TargetCalculatorTest::test_ideal_shape_quad()
{
  MsqError err;
  PatchData pd;
  MsqMatrix<2,2> Q;
  
  TargetCalculator::ideal_shape_2D( QUADRILATERAL, Sample(0,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  ASSERT_IDENTITY_MATRIX( Q );
  
  TargetCalculator::ideal_shape_2D( QUADRILATERAL, Sample(2,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  ASSERT_IDENTITY_MATRIX( Q );
}

void TargetCalculatorTest::test_ideal_shape_tet()
{
    // Ideal tetrahedron should have Aspect (i.e. delta) == identity,
    // so shape should be equal to skew.

  MsqError err;
  PatchData pd;
  MsqMatrix<3,3> Q, Q2;
  
  TargetCalculator::ideal_shape_3D( TETRAHEDRON, Sample(0,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  TargetCalculator::ideal_skew_3D( TETRAHEDRON, Sample(0,0), pd, Q2, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q2, Q, EPS );
  
  TargetCalculator::ideal_shape_3D( TETRAHEDRON, Sample(3,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  TargetCalculator::ideal_skew_3D( TETRAHEDRON, Sample(3,0), pd, Q2, err );
  ASSERT_NO_ERROR(err);
  ASSERT_MATRICES_EQUAL( Q2, Q, EPS );
}

void TargetCalculatorTest::test_ideal_shape_prism()
{
    // Ideal wedge should have Aspect (i.e. delta) == identity,
    // so shape should be equal to skew.

  Sample points[] = { Sample(0,0), Sample(0,1), Sample(0,2), Sample(0,3),
                      Sample(2,0), Sample(3,0) };
  const int num_pts = sizeof(points)/sizeof(points[0]);
  
  MsqError err;
  PatchData pd;
  MsqMatrix<3,3> Q, W;
 
  for (int i = 0; i < num_pts; ++i) {
    TargetCalculator::ideal_shape_3D( PRISM, points[i], pd, Q, err );
    ASSERT_NO_ERROR(err);
    TargetCalculator::ideal_skew_3D( PRISM, points[i], pd, W, err );
    ASSERT_NO_ERROR(err);
    ASSERT_MATRICES_EQUAL( W, Q, EPS );
  }
}

void TargetCalculatorTest::test_ideal_shape_pyramid()
{
  Sample points[] = { Sample(0,0), Sample(0,1), Sample(0,2), Sample(0,3),
                      Sample(2,0), Sample(3,0) };
  const int num_pts = sizeof(points)/sizeof(points[0]);
  
  MsqError err;
  PatchData pd;
  MsqMatrix<3,3> Q, W;
  const Vector3D* coords = unit_edge_element( PYRAMID, true );
  
  for (int i = 0; i < num_pts; ++i) {
  
    TargetCalculator::ideal_shape_3D( PYRAMID, points[i], pd, Q, err );
    ASSERT_NO_ERROR(err);
    
    TargetCalculator::jacobian_3D( pd, PYRAMID, 5, points[i], coords, W, err );
    ASSERT_NO_ERROR(err);
    
    ASSERT_MATRICES_EQUAL( TargetCalculator::shape(W), Q, EPS );
  }
}

void TargetCalculatorTest::test_ideal_shape_hex()
{
  MsqError err;
  PatchData pd;
  MsqMatrix<3,3> Q;
  
  TargetCalculator::ideal_shape_3D( HEXAHEDRON, Sample(0,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  ASSERT_IDENTITY_MATRIX( Q );
  
  TargetCalculator::ideal_shape_3D( HEXAHEDRON, Sample(3,0), pd, Q, err );
  ASSERT_NO_ERROR(err);
  ASSERT_IDENTITY_MATRIX( Q );
}
  
void TargetCalculatorTest::test_new_orientatin_3D()
{
  MsqVector<3> x(0.0), y(0.0), z(0.0);
  x[0] = y[1] = z[2] = 1;
  
  MsqMatrix<3,3> V = TargetCalculator::new_orientation_3D( x, y );
  check_valid_V(V);
  ASSERT_MATRICES_EQUAL( z, V.column(0) * V.column(1), EPS );
  ASSERT_MATRICES_EQUAL( z, V.column(2), EPS );
}

void TargetCalculatorTest::test_new_orientatin_2D()
{
  MsqVector<3> x(0.0), y(0.0), z(0.0);
  x[0] = y[1] = z[2] = 1;
  
  MsqMatrix<3,2> V = TargetCalculator::new_orientation_2D( x, y );
  check_valid_V(V);
  ASSERT_MATRICES_EQUAL( z, V.column(0) * V.column(1), EPS );
}
  
void TargetCalculatorTest::test_new_aspect_3D()
{
  MsqVector<3> r;
  r[0] = 2; r[1] = 4; r[2] = 6;
  MsqMatrix<3,3> delta = TargetCalculator::new_aspect_3D(r);
  check_valid_delta(delta);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( r[0]/r[1], delta(0,0)/delta(1,1), EPS );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( r[0]/r[2], delta(0,0)/delta(2,2), EPS );
}

void TargetCalculatorTest::test_new_aspect_2D()
{
  MsqVector<2> r;
  r[0] = 7; r[1] = 3;
  MsqMatrix<2,2> delta = TargetCalculator::new_aspect_2D(r);
  check_valid_delta(delta);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( r[0]/r[1], delta(0,0)/delta(1,1), EPS );
}
    
void TargetCalculatorTest::test_jacobian_3D()
{
  MsqError err;
  PatchData pd;
  LinearHexahedron map;
  const Vector3D hex_coords[] = { Vector3D(0,0,0),
                                  Vector3D(2,0,0),
                                  Vector3D(3,2,0),
                                  Vector3D(1,2,0),
                                  Vector3D(1,0,1),
                                  Vector3D(3,0,2),
                                  Vector3D(3,2,2),
                                  Vector3D(2,2,2) };
  const MsqVector<3>* coords = reinterpret_cast<const MsqVector<3>*>(hex_coords);
  
  Sample pts[] = { Sample(0,0), Sample(1,6), Sample(3,0) };
  const int num_pts = sizeof(pts)/sizeof(pts[0]);
  for (int i = 0; i < num_pts; ++i) {
    MsqMatrix<3,3> J(0.0), W;
    size_t indices[8], n;
    MsqVector<3> derivs[8];
    map.derivatives( pts[i], NodeSet(), indices, derivs, n, err );
    ASSERT_NO_ERROR(err);
    for (size_t j = 0; j < n; ++j)
      J += outer( coords[indices[j]], derivs[j] );
    
    TargetCalculator::jacobian_3D( pd, HEXAHEDRON, 8, pts[i], hex_coords, W, err );
    ASSERT_NO_ERROR(err);
    
    ASSERT_MATRICES_EQUAL( J, W, EPS );
  }
}

void TargetCalculatorTest::test_jacobian_2D()
{
  MsqError err;
  PatchData pd;
  LinearTriangle map;
  const Vector3D tri_coords[] = { Vector3D(0,0,0),
                                  Vector3D(2,0,0),
                                  Vector3D(0,0,1) };
  const MsqVector<3>* coords = reinterpret_cast<const MsqVector<3>*>(tri_coords);
  
  Sample pts[] = { Sample(0,0), Sample(1,1), Sample(2,0) };
  const int num_pts = sizeof(pts)/sizeof(pts[0]);
  for (int i = 0; i < num_pts; ++i) {
    MsqMatrix<3,2> J(0.0), W;
    size_t indices[3], n;
    MsqVector<2> derivs[3];
    map.derivatives( pts[i], NodeSet(), indices, derivs, n, err );
    ASSERT_NO_ERROR(err);
    for (size_t j = 0; j < n; ++j)
      J += outer( coords[indices[j]], derivs[j] );
    
    TargetCalculator::jacobian_2D( pd, TRIANGLE, 3, pts[i], tri_coords, W, err );
    ASSERT_NO_ERROR(err);
    
    ASSERT_MATRICES_EQUAL( J, W, EPS );
  }
}


class DummyRefMesh : public ReferenceMeshInterface
{
  std::vector<Vector3D> mCoords;
public:
  DummyRefMesh( const Vector3D* coords, int num_verts )
    : mCoords(coords, coords+num_verts) {}
    
  void get_reference_vertex_coordinates( const Mesh::VertexHandle* handles,
                                         const size_t num_vertices,
                                         Vector3D* coords, 
                                         MsqError& err )
  {
    const size_t* indices = reinterpret_cast<const size_t*>(handles);
    for (size_t i = 0; i < num_vertices; ++i) {
      if (i >= mCoords.size()) {
        MSQ_SETERR(err)(MsqError::INVALID_ARG);
        return;
      }
      coords[i] = mCoords[indices[i]];
    }
  }
};
    

void TargetCalculatorTest::test_get_refmesh_Jacobian_3D()
{
  MsqError err;

  const Vector3D hex_coords[] = { Vector3D(0,0,0),
                                  Vector3D(2,0,0),
                                  Vector3D(3,2,0),
                                  Vector3D(1,2,0),
                                  Vector3D(1,0,1),
                                  Vector3D(3,0,2),
                                  Vector3D(3,2,2),
                                  Vector3D(2,2,2) };
  DummyRefMesh ref_mesh( hex_coords, 8 );                                
                                  
  const Vector3D rect_coords[] = { Vector3D(0,0,0),
                                   Vector3D(1,0,0),
                                   Vector3D(1,1,0),
                                   Vector3D(0,1,0),
                                   Vector3D(0,0,5),
                                   Vector3D(1,0,5),
                                   Vector3D(1,1,5),
                                   Vector3D(0,1,5) };
  size_t conn[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  bool fixed[8];
  std::fill( fixed, fixed + sizeof(fixed)/sizeof(fixed[0]), false );
  PatchData pd;
  pd.fill( 8, rect_coords[0].to_array(), 1, HEXAHEDRON, conn, fixed, err );
  ASSERT_NO_ERROR(err);
  
  Sample pts[] = { Sample(0,0), Sample(1,6), Sample(3,0) };
  const int num_pts = sizeof(pts)/sizeof(pts[0]);
  for (int i = 0; i < num_pts; ++i) {
    MsqMatrix<3,3> J, W;
    TargetCalculator::jacobian_3D( pd, HEXAHEDRON, 8, pts[i], hex_coords, J, err );
    ASSERT_NO_ERROR(err);
    
    TargetCalculator::get_refmesh_Jacobian_3D( &ref_mesh, pd, 0, pts[i], W, err );
    ASSERT_NO_ERROR(err);
    
    ASSERT_MATRICES_EQUAL( J, W, EPS );
  }
}

void TargetCalculatorTest::test_get_refmesh_Jacobian_2D()
{
  MsqError err;

  const Vector3D tri_coords[] = { Vector3D(0,0,0),
                                  Vector3D(3,0,1),
                                  Vector3D(1,0,3) };
  DummyRefMesh ref_mesh( tri_coords, 3 );                                
                                  
  const Vector3D right_coords[] = { Vector3D(0,0,0),
                                    Vector3D(5,0,0),
                                    Vector3D(5,1,0) };
  size_t conn[] = { 0, 1, 2 };
  bool fixed[3];
  std::fill( fixed, fixed + sizeof(fixed)/sizeof(fixed[0]), false );
  PatchData pd;
  pd.fill( 3, right_coords[0].to_array(), 1, TRIANGLE, conn, fixed, err );
  ASSERT_NO_ERROR(err);
  
  Sample pts[] = { Sample(0,0), Sample(1,1), Sample(2,0) };
  const int num_pts = sizeof(pts)/sizeof(pts[0]);
  for (int i = 0; i < num_pts; ++i) {
    MsqMatrix<3,2> J, W;
    TargetCalculator::jacobian_2D( pd, TRIANGLE, 3, pts[i], tri_coords, J, err );
    ASSERT_NO_ERROR(err);
    
    TargetCalculator::get_refmesh_Jacobian_2D( &ref_mesh, pd, 0, pts[i], W, err );
    ASSERT_NO_ERROR(err);
    
    ASSERT_MATRICES_EQUAL( J, W, EPS );
  }
}

