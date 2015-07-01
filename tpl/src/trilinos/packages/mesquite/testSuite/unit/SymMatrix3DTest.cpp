/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


#include "SymMatrix3D.hpp"

#include <math.h>

#include "UnitUtil.hpp"

using namespace Mesquite;

class SymMatrix3DTest : public CppUnit::TestFixture
{

private:
  CPPUNIT_TEST_SUITE(SymMatrix3DTest);
  CPPUNIT_TEST (test_init_diag);
  CPPUNIT_TEST (test_indices);
  CPPUNIT_TEST (test_plus_eq);
  CPPUNIT_TEST (test_minus_eq);
  CPPUNIT_TEST (test_times_eq);
  CPPUNIT_TEST (test_divide_eq);
  CPPUNIT_TEST (test_plus);
  CPPUNIT_TEST (test_minus);
  CPPUNIT_TEST (test_plus_nonsym);
  CPPUNIT_TEST (test_minus_nonsym);
  CPPUNIT_TEST (test_multiply);
  CPPUNIT_TEST (test_scalar_multiply);
  CPPUNIT_TEST (test_divide);
  CPPUNIT_TEST (test_vector_multiply);
  CPPUNIT_TEST (test_nonsym_multiply);
  CPPUNIT_TEST (test_determinant);
  CPPUNIT_TEST (test_inverse);
  CPPUNIT_TEST_SUITE_END();

public:

  void test_init_diag();
  void test_indices();
  void test_plus_eq();
  void test_minus_eq();
  void test_times_eq();
  void test_divide_eq();
  void test_plus();
  void test_minus();
  void test_plus_nonsym();
  void test_minus_nonsym();
  void test_multiply();
  void test_scalar_multiply();
  void test_divide();
  void test_vector_multiply();
  void test_nonsym_multiply();
  void test_determinant();
  void test_inverse();

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(SymMatrix3DTest, "SymMatrix3DTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(SymMatrix3DTest, "Unit");

void SymMatrix3DTest::test_init_diag()
{
  const SymMatrix3D A( 2.0 );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, A[SymMatrix3D::T00], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A[SymMatrix3D::T01], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A[SymMatrix3D::T02], DBL_EPSILON );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A[SymMatrix3D::T10], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, A[SymMatrix3D::T11], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A[SymMatrix3D::T12], DBL_EPSILON );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A[SymMatrix3D::T20], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A[SymMatrix3D::T21], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, A[SymMatrix3D::T22], DBL_EPSILON );
}

void SymMatrix3DTest::test_indices()
{
  const SymMatrix3D A( 1, 2, 3, 4, 5, 6 );
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, A(0,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, A(0,1), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.0, A(0,2), DBL_EPSILON );
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, A(1,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 4.0, A(1,1), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, A(1,2), DBL_EPSILON );
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.0, A(2,0), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 5.0, A(2,1), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 6.0, A(2,2), DBL_EPSILON );
}

void SymMatrix3DTest::test_plus_eq()
{
  const SymMatrix3D A( 1, 2, 3, 4, 5, 6 );
  const SymMatrix3D B( 101, 102, 103, 104, 105, 106 );
  SymMatrix3D C(A);
  C += B;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[0]+B[0], C[0], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[1]+B[1], C[1], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[2]+B[2], C[2], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[3]+B[3], C[3], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[4]+B[4], C[4], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[5]+B[5], C[5], DBL_EPSILON );
}

void SymMatrix3DTest::test_minus_eq()
{
  const SymMatrix3D A( 1, 2, 3, 4, 5, 6 );
  const SymMatrix3D B( 101, 102, 103, 104, 105, 106 );
  SymMatrix3D C(A);
  C -= B;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[0]-B[0], C[0], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[1]-B[1], C[1], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[2]-B[2], C[2], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[3]-B[3], C[3], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[4]-B[4], C[4], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[5]-B[5], C[5], DBL_EPSILON );
}

void SymMatrix3DTest::test_times_eq()
{
  const SymMatrix3D A( 0.5, 2, 3, 4, 5, 6 );
  const double s = 0.25;
  SymMatrix3D C(A);
  C *= s;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( s*A[0], C[0], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( s*A[1], C[1], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( s*A[2], C[2], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( s*A[3], C[3], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( s*A[4], C[4], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( s*A[5], C[5], DBL_EPSILON );
}

void SymMatrix3DTest::test_divide_eq()
{
  const SymMatrix3D A( 0.5, 2, 3, 4, 5, 6 );
  const double s = 0.25;
  SymMatrix3D C(A);
  C /= s;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[0]/s, C[0], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[1]/s, C[1], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[2]/s, C[2], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[3]/s, C[3], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[4]/s, C[4], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[5]/s, C[5], DBL_EPSILON );
}

void SymMatrix3DTest::test_plus()
{
  const SymMatrix3D A( 1, 2, 3, 4, 5, 6 );
  const SymMatrix3D B( 101, 102, 103, 104, 105, 106 );
  const SymMatrix3D C(A+B);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[0]+B[0], C[0], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[1]+B[1], C[1], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[2]+B[2], C[2], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[3]+B[3], C[3], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[4]+B[4], C[4], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[5]+B[5], C[5], DBL_EPSILON );
}

void SymMatrix3DTest::test_minus()
{
  const SymMatrix3D A( 1, 2, 3, 4, 5, 6 );
  const SymMatrix3D B( 101, 102, 103, 104, 105, 106 );
  const SymMatrix3D C(A-B);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[0]-B[0], C[0], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[1]-B[1], C[1], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[2]-B[2], C[2], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[3]-B[3], C[3], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[4]-B[4], C[4], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[5]-B[5], C[5], DBL_EPSILON );
}

void SymMatrix3DTest::test_plus_nonsym()
{
  const SymMatrix3D A( 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625 );
  const Matrix3D B( 1, 2, 3, 4, 5, 6, 7, 8, 9 );
  const Matrix3D C( A + B );
  const Matrix3D D( B + A );
  const Matrix3D E( A(0,0)+B(0,0), A(0,1)+B(0,1), A(0,2)+B(0,2),
                    A(1,0)+B(1,0), A(1,1)+B(1,1), A(1,2)+B(1,2),
                    A(2,0)+B(2,0), A(2,1)+B(2,1), A(2,2)+B(2,2) );
  CPPUNIT_ASSERT_MATRICES_EQUAL( E, C, DBL_EPSILON );
  CPPUNIT_ASSERT_MATRICES_EQUAL( E, D, DBL_EPSILON );
}

void SymMatrix3DTest::test_minus_nonsym()
{
  const SymMatrix3D A( 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625 );
  const Matrix3D B( 1, 2, 3, 4, 5, 6, 7, 8, 9 );

  const Matrix3D C( A - B );
  const Matrix3D E( A(0,0)-B(0,0), A(0,1)-B(0,1), A(0,2)-B(0,2),
                    A(1,0)-B(1,0), A(1,1)-B(1,1), A(1,2)-B(1,2),
                    A(2,0)-B(2,0), A(2,1)-B(2,1), A(2,2)-B(2,2) );
  CPPUNIT_ASSERT_MATRICES_EQUAL( E, C, DBL_EPSILON );
  
  const Matrix3D D( B - A );
  const Matrix3D F( B(0,0)-A(0,0), B(0,1)-A(0,1), B(0,2)-A(0,2),
                    B(1,0)-A(1,0), B(1,1)-A(1,1), B(1,2)-A(1,2),
                    B(2,0)-A(2,0), B(2,1)-A(2,1), B(2,2)-A(2,2) );
  CPPUNIT_ASSERT_MATRICES_EQUAL( F, D, DBL_EPSILON );
}

void SymMatrix3DTest::test_multiply()
{
  const SymMatrix3D A( 1, 2, 3, 4, 5, 6 );
  const SymMatrix3D B( 101, 102, 103, 104, 105, 106 );
  const Matrix3D C(A), D(B);
  CPPUNIT_ASSERT_MATRICES_EQUAL( C*D, A*B, DBL_EPSILON );
  CPPUNIT_ASSERT_MATRICES_EQUAL( D*C, B*A, DBL_EPSILON );
}

void SymMatrix3DTest::test_scalar_multiply()
{
  const SymMatrix3D A( 1, 2, 3, 4, 5, 6 );
  const double s = 0.25;
  const SymMatrix3D B(A*s), C(s*A);

  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[0]*s, B[0], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[1]*s, B[1], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[2]*s, B[2], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[3]*s, B[3], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[4]*s, B[4], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[5]*s, B[5], DBL_EPSILON );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[0]*s, C[0], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[1]*s, C[1], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[2]*s, C[2], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[3]*s, C[3], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[4]*s, C[4], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[5]*s, C[5], DBL_EPSILON );
}

void SymMatrix3DTest::test_divide()
{
  const SymMatrix3D A( 1, 2, 3, 4, 5, 6 );
  const double s = 0.25;
  const SymMatrix3D B(A/s);

  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[0]/s, B[0], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[1]/s, B[1], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[2]/s, B[2], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[3]/s, B[3], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[4]/s, B[4], DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( A[5]/s, B[5], DBL_EPSILON );
}

void SymMatrix3DTest::test_vector_multiply()
{
  const SymMatrix3D A( 1, 2, 3, 4, 5, 6 );
  const Vector3D v(101, 102, 103);
  const Matrix3D B(A);
  CPPUNIT_ASSERT_VECTORS_EQUAL( B*v, A*v, DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( v*B, v*A, DBL_EPSILON );
}

void SymMatrix3DTest::test_nonsym_multiply()
{
  const SymMatrix3D A( 1, 2, 3, 4, 5, 6 );
  const Matrix3D B( 101, 102, 103, 104, 105, 106, 107, 108, 109 );
  const Matrix3D C(A);
  CPPUNIT_ASSERT_MATRICES_EQUAL( C*B, A*B, DBL_EPSILON );
  CPPUNIT_ASSERT_MATRICES_EQUAL( B*C, B*A, DBL_EPSILON );
}

void SymMatrix3DTest::test_determinant()
{
  const SymMatrix3D A( 1, 2, 3, 4, 5, 6 );
  const Matrix3D B(A);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( det(B), det(A), DBL_EPSILON );
}

void SymMatrix3DTest::test_inverse()
{
  const SymMatrix3D A( 1, 2, 3, 4, 5, 6 );
  const Matrix3D B( A * inverse(A) );
  const Matrix3D I( 1, 0, 0, 0, 1, 0, 0, 0, 1 );
  CPPUNIT_ASSERT_MATRICES_EQUAL( I, B, 1e-8 );
}

