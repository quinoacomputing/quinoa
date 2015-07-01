/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retian certain rights to this software.

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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file MsqMatrixTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "MsqMatrix.hpp"
#include <sstream>
#include "cppunit/extensions/HelperMacros.h"
#include "UnitUtil.hpp"

using namespace Mesquite;
using namespace std;

class MsqMatrixTest : public CppUnit::TestFixture
{
  private:
    CPPUNIT_TEST_SUITE( MsqMatrixTest );
    
    CPPUNIT_TEST(test_initialize);
    CPPUNIT_TEST(test_assign);
    CPPUNIT_TEST(test_matrix_multiply);
    CPPUNIT_TEST(test_scalar_multiply);
    CPPUNIT_TEST(test_matrix_add);
    CPPUNIT_TEST(test_scalar_add);
    CPPUNIT_TEST(test_matrix_subtract);
    CPPUNIT_TEST(test_scalar_subtract);
    CPPUNIT_TEST(test_get_set_row);
    CPPUNIT_TEST(test_get_set_column);
    CPPUNIT_TEST(test_inverse);
    CPPUNIT_TEST(test_qr);
    CPPUNIT_TEST(test_frobenius);
    CPPUNIT_TEST(test_vec_length);
    CPPUNIT_TEST(test_determinant);
    CPPUNIT_TEST(test_vec_outer_product);
    
    CPPUNIT_TEST_SUITE_END();
    
  public:

    void test_initialize();
    void test_assign();
    void test_matrix_multiply();
    void test_scalar_multiply();
    void test_matrix_add();
    void test_scalar_add();
    void test_matrix_subtract();
    void test_scalar_subtract();
    void test_get_set_row();
    void test_get_set_column();
    void test_inverse();
    void test_qr();
    void test_frobenius();
    void test_vec_length();
    void test_determinant();
    void test_vec_outer_product();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqMatrixTest, "MsqMatrixTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqMatrixTest, "Unit");

void MsqMatrixTest::test_initialize()
{
  MsqMatrix<4,4> I44(1.0);
  ASSERT_IDENTITY_MATRIX( I44 );
  
  double values[] = { 1.0, 0.0, 0.0, 
                      0.0, 1.0, 0.0,
                      0.0, 0.0, 1.0 };
  MsqMatrix<3,3> I33(values);
  ASSERT_IDENTITY_MATRIX( I33 );
  
  double c0v[] = { 1.0, 0.0 }, c1v[] = { 0.0, 1.0 };
  MsqMatrix<2,1> cols[] = { c0v, c1v };
  MsqMatrix<2,2> I22( cols );
  ASSERT_IDENTITY_MATRIX( I22 );
  
  double r0v[] = { 1.0, 2.0 }, r1v[] = { 3.0, 4.0 };
  MsqMatrix<1,2> rows[] = { r0v, r1v };
  MsqMatrix<2,2> I( rows );
  CPPUNIT_ASSERT_DOUBLES_EQUAL(r0v[0], I(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL(r0v[1], I(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL(r1v[0], I(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL(r1v[1], I(1,1), 1e-6 );
  
  std::string str = "-1 -2 1 2 3 4";
  MsqMatrix<2,3> S1( str );
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1, S1(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-2, S1(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1, S1(0,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2, S1(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 3, S1(1,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 4, S1(1,2), 1e-6 );

  MsqMatrix<2,3> S2( str.c_str() );
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1, S2(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-2, S2(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1, S2(0,2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2, S2(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 3, S2(1,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 4, S2(1,2), 1e-6 );
  
  MsqMatrix<2,2> II( I33, 1, 1 );
  ASSERT_IDENTITY_MATRIX( II );
}
  
  
void MsqMatrixTest::test_assign()
{
  MsqMatrix<3,2> m;
  for (int i = 0; i < 6; ++i) {
    m(i/2,i%2) = (double)i;
    CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)i, m(i/2,i%2), DBL_EPSILON );
  }
}

void MsqMatrixTest::test_matrix_multiply()
{
  double v1[] = { 1,  2,  3, 
                 -3, -2, -1 };
  double v2[] = { 4,   3,  
                  5,  10, 
                  2, 0.5 };

  MsqMatrix<2,3> A(v1);
  MsqMatrix<3,2> B(v2);
  MsqMatrix<2,2> R = A * B;
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 20.0, R(0,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 24.5, R(0,1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-24.0, R(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-29.5, R(1,1), 1e-6 );
}

void MsqMatrixTest::test_scalar_multiply()
{
  MsqMatrix<4,4> A(33.0);
  MsqMatrix<4,4> B(1.0/33 * A);
  A *= 1.0/33.0;
  ASSERT_IDENTITY_MATRIX( A );
  ASSERT_IDENTITY_MATRIX( B );
}

void MsqMatrixTest::test_matrix_add()
{
  double v1[] = { 3, 4, 5,
                  2, 1, 0,
                  7, 8, 9 };
  double v2[] = {-2,-4,-5,
                 -2, 0, 0,
                 -7,-8,-8 };
  MsqMatrix<3,3> A(v1), B(v2);
  MsqMatrix<3,3> R = A + B;
  ASSERT_IDENTITY_MATRIX( R );
  A += B;
  ASSERT_IDENTITY_MATRIX( A );
}

void MsqMatrixTest::test_scalar_add()
{
  double v[] = { -6, -7, -7, -6 };
  MsqMatrix<2,2> A(v), B(A+7);
  A += 7;
  ASSERT_IDENTITY_MATRIX( A );
  ASSERT_IDENTITY_MATRIX( B );
}

void MsqMatrixTest::test_matrix_subtract()
{
  double v1[] = { 3, 4, 5,
                  2, 1, 0,
                  7, 8, 9 };
  double v2[] = { 2, 4, 5,
                  2, 0, 0,
                  7, 8, 8 };
  MsqMatrix<3,3> A(v1), B(v2);
  MsqMatrix<3,3> R = A - B;
  ASSERT_IDENTITY_MATRIX( R );
  A -= B;
  ASSERT_IDENTITY_MATRIX( A );
}

void MsqMatrixTest::test_scalar_subtract()
{
  double v[] = { -6, -7, -7, -6 };
  MsqMatrix<2,2> A(v), B(A - -7.0);
  A -= -7.0;
  ASSERT_IDENTITY_MATRIX( A );
  ASSERT_IDENTITY_MATRIX( B );
}

void MsqMatrixTest::test_get_set_row()
{
  double values[] = { -1, 1, -2, 
                       2, -3, 3, 
                      -4, 4, -5 };
  MsqMatrix<3,3> A(values);
  MsqMatrix<1,3> r0 = A.row(0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[0], r0(0,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[1], r0(0,1), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[2], r0(0,2), DBL_EPSILON );
  MsqMatrix<1,3> r1 = A.row(1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[3], r1(0,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[4], r1(0,1), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[5], r1(0,2), DBL_EPSILON );
  MsqMatrix<1,3> r2 = A.row(2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[6], r2(0,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[7], r2(0,1), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[8], r2(0,2), DBL_EPSILON );
  
  r0(0,0) = 20;
  r0(0,1) = 30;
  r0(0,2) = 10;
  A.set_row( 0, r0 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A(0,0), 20, DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A(0,1), 30, DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A(0,2), 10, DBL_EPSILON );
}

void MsqMatrixTest::test_get_set_column()
{
  double values[] = { -1, 1, -2, 
                       2, -3, 3, 
                      -4, 4, -5 };
  MsqMatrix<3,3> A(values);
  MsqMatrix<3,1> c0 = A.column(0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[0], c0(0,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[3], c0(1,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[6], c0(2,0), DBL_EPSILON );
  MsqMatrix<3,1> c1 = A.column(1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[1], c1(0,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[4], c1(1,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[7], c1(2,0), DBL_EPSILON );
  MsqMatrix<3,1> c2 = A.column(2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[2], c2(0,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[5], c2(1,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( values[8], c2(2,0), DBL_EPSILON );
  
  c0(0,0) = 20;
  c0(1,0) = 30;
  c0(2,0) = 10;
  A.set_column( 0, c0 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A(0,0), 20, DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A(1,0), 30, DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A(2,0), 10, DBL_EPSILON );
}

void MsqMatrixTest::test_inverse()
{
  double v2[] = { 4, 5, 6, 7 };
  MsqMatrix<2,2> M2(v2);
  MsqMatrix<2,2> M2I = inverse(M2);
  MsqMatrix<2,2> M2P = M2 * M2I;
  ASSERT_IDENTITY_MATRIX( M2P );
  
  double v3[] = { 0.1, 0.5, 0.3, 1, 2, 3, 0.5, 0.1, 10 };
  MsqMatrix<3,3> M3(v3);
  MsqMatrix<3,3> M3I = inverse(M3);
  MsqMatrix<3,3> M3P = M3 * M3I;
  ASSERT_IDENTITY_MATRIX( M3P );
}

void MsqMatrixTest::test_qr()
{
  double v2[] = { 4, 5, 6, 7 };
  MsqMatrix<2,2> M2(v2), Q2, R2;
  QR( M2, Q2, R2 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, R2(1,0), 1e-6 );
  MsqMatrix<2,2> P2 = Q2 * R2;
  ASSERT_MATRICES_EQUAL( M2, P2, 1e-6 );
  MsqMatrix<2,2> I2 = Q2 * transpose(Q2);
  ASSERT_IDENTITY_MATRIX( I2 );

  double v3[] = { 0.1, 0.5, 0.3, 1, 2, 3, 0.5, 0.1, 10 };
  MsqMatrix<3,3> M3(v3), Q3, R3;
  QR( M3, Q3, R3 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, R3(1,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, R3(2,0), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, R3(2,1), 1e-6 );
  MsqMatrix<3,3> P3 = Q3 * R3;
  ASSERT_MATRICES_EQUAL( M3, P3, 1e-6 );
  MsqMatrix<3,3> I3 = Q3 * transpose(Q3);
  ASSERT_IDENTITY_MATRIX( I3 );
}

void MsqMatrixTest::test_frobenius()
{
  double v3[] = { 0.1, 0.5, 0.3, 1, 2, 3, 0.5, 0.1, 10 };
  MsqMatrix<3,3> M3(v3);
  double exp = 0.0;
  for (unsigned i = 0; i < 9; ++i)
    exp += v3[i]*v3[i];
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, sqr_Frobenius( M3 ), 1e-10 );
}

void MsqMatrixTest::test_vec_length()
{
  double values[] = { 2, 4, 4 };
  MsqMatrix<3,1> c( values );
  MsqMatrix<1,3> r( values );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 36.0, sqr_length(c), 1e-12 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 36.0, sqr_length(r), 1e-12 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 6.0, length(c), 1e-12 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 6.0, length(r), 1e-12 );
}

void MsqMatrixTest::test_determinant()
{
  MsqMatrix<5,5> m5(1.0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(m5), 1e-8 );
  
  double m4vals[] = { 1, 2, 3, 4,
                      5, 6, 7, 8, 
                      9,10,11,12,
                     13,14,15,16 };
  MsqMatrix<4,4> m4(m4vals);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, det(m4), 1e-8 );
  
  double m4vals2[] = { 2, 5, 0, 0,
                       1, 1, 0, 0,
                       0, 0, 1, 0,
                       0, 0, 0, 1 };
  MsqMatrix<4,4> m4b( m4vals2 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( -3.0, det(m4b), 1e-8 );
  
  double m3vals[] = { 1, 2, -1,
                     -1, 0,  1,
                      1, 1,  3 };
  MsqMatrix<3,3> m3(m3vals);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.0, det(m3), 1e-8 );
  
  double m2vals[] = { sqrt(2.0), 1, 1, sqrt(2.0) };
  MsqMatrix<2,2> m2(m2vals);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(m2), 1e-8 );
}

void MsqMatrixTest::test_vec_outer_product()
{
  double v1[] = { 1, 3, -1, 12 };
  double v2[] = { -10, 4, 0.5, -9, 2 };
  MsqMatrix<4,1> V1(v1);
  MsqMatrix<5,1> V2(v2);
  
  MsqMatrix<4,5> M1 = V1 * transpose(V2);
  MsqMatrix<4,5> M2 = outer( V1, V2 );
  ASSERT_MATRICES_EQUAL( M1, M2, 1e-8 );
}
  
