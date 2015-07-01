/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

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


/** \file UnitUtil.hpp
 *  \brief Utility functions for use in unit tests
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_UNIT_UTIL_HPP
#define MSQ_UNIT_UTIL_HPP

#include "Mesquite.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"

#include <string>
#include <stdio.h>

#include "cppunit/extensions/HelperMacros.h"

#define ASSERT_MESSAGE( MSG, COND ) \
  CPPUNIT_NS::Asserter::failIf( !(COND), (MSG), CPPUNIT_SOURCELINE() )

/** Assert that Mesquite API has not flagged an error */
#define ASSERT_NO_ERROR( A ) \
  ASSERT_MESSAGE( (A).error_message(), ! MSQ_CHKERR( (A) ) )

/**\brief compare two vectors (Vector3D)
 *
 * Ensure that the test result Vector3D \a v2 is within \a eps of the
 * expected vector \a v1 .  
 */
#define CPPUNIT_ASSERT_VECTORS_EQUAL( v1, v2, eps ) \
  ASSERT_MESSAGE( utest_vect_message((v1),(v2)), \
                          utest_vect_equal((v1),(v2),(eps)) )

/**\brief compare two matrices (Matrix3D)
 *
 * Ensure that the test result Matrix3D \a m2 is within \a eps of the
 * expected matrix \a m1
 */
#define CPPUNIT_ASSERT_MATRICES_EQUAL( m1, m2, eps ) \
  ASSERT_MESSAGE( utest_mat_message((m1),(m2)), \
                          utest_mat_equal((m1),(m2),(eps)) )


/** compare matrix (MsqMatrix) with the identity matrix  */
#define ASSERT_IDENTITY_MATRIX( M ) \
  ASSERT_MESSAGE( ident_check_msg(M), ident_check(M) )

/** compare two matrices (MsqMatrix) */
#define ASSERT_MATRICES_EQUAL( A, B, E ) \
  ASSERT_MESSAGE( mat_equal_check_msg( A, B), mat_equal_check(A,B,E) )

/** compare two matrices (MsqMatrix) */
#define ASSERT_MATRICES_DIFFERENT( A, B, E ) \
  ASSERT_MESSAGE( mat_not_equal_check_msg(A, B), !mat_equal_check(A,B,E) )

/** compare two arrays of values */
#define ASSERT_ARRAYS_EQUAL( A, B, LEN ) \
  CPPUNIT_NS::Asserter::failIf( !(arrays_equal((A),(B),(LEN))), arrays_not_equal_msg((A),(LEN),(B),(LEN)), CPPUNIT_SOURCELINE() )

#define ASSERT_STD_VECTORS_EQUAL( A, B ) \
  CPPUNIT_NS::Asserter::failIf( ((A) != (B)), arrays_not_equal_msg(&(A)[0],(A).size(),&(B)[0],(B).size()), CPPUNIT_SOURCELINE() )


/** make string representation of cartesian vector */
inline std::string utest_vect_str( const Mesquite::Vector3D& v )
{
  char buffer[128];
  sprintf(buffer, "[%f, %f, %f]", v[0], v[1], v[2]);
  return buffer;
}

/** make string representation of 3x3 matrix */
inline std::string utest_mat_str( const Mesquite::Matrix3D& m )
{
  char buffer[256];
  sprintf(buffer, "[%f, %f, %f] [%f, %f, %f] [%f, %f, %f]", 
          m[0][0], m[0][1], m[0][2],
          m[1][0], m[1][1], m[1][2],
          m[2][0], m[2][1], m[2][2] );
  return buffer;
}

/** make string representation of 3x3 symetric matrix */
inline std::string utest_mat_str( const Mesquite::SymMatrix3D& m )
{
  char buffer[256];
  sprintf(buffer, "[%f, %f, %f] [%f, %f, %f] [%f, %f, %f]", 
          m(0,0), m(0,1), m(0,2),
          m(1,0), m(1,1), m(1,2),
          m(2,0), m(2,1), m(2,2) );
  return buffer;
}

/** make error message for failed vector copmarison */
inline CppUnit::Message utest_vect_message( const Mesquite::Vector3D& v1,
                                      const Mesquite::Vector3D& v2 )
{
  CppUnit::Message m( "equality assertion failed" );
  m.addDetail( std::string("Expected: ") + utest_vect_str(v1) );
  m.addDetail( std::string("Actual  : ") + utest_vect_str(v2) );
  return m;
}

/** make error message for failed matrix copmarison */
inline CppUnit::Message utest_mat_message( const Mesquite::Matrix3D& m1,
                                     const Mesquite::Matrix3D& m2 )
{
  CppUnit::Message m( "equality assertion failed" );
  m.addDetail( std::string("Expected: ") + utest_mat_str(m1) );
  m.addDetail( std::string("Actual  : ") + utest_mat_str(m2) );
  return m;
}

/** make error message for failed symmetric matrix copmarison */
inline CppUnit::Message utest_mat_message( const Mesquite::SymMatrix3D& m1,
                                     const Mesquite::SymMatrix3D& m2 )
{
  CppUnit::Message m( "equality assertion failed" );
  m.addDetail( std::string("Expected: ") + utest_mat_str(m1) );
  m.addDetail( std::string("Actual  : ") + utest_mat_str(m2) );
  return m;
}

/** compare vectors */
inline bool utest_vect_equal( const Mesquite::Vector3D& v1, const Mesquite::Vector3D& v2, double eps )
{
  return (fabs(v1[0] - v2[0]) < eps) &&
         (fabs(v1[1] - v2[1]) < eps) &&
         (fabs(v1[2] - v2[2]) < eps);
}

/** compare matrices */
inline bool utest_mat_equal( const Mesquite::Matrix3D& m1, const Mesquite::Matrix3D& m2, double eps )
{
  return (fabs(m1[0][0] - m2[0][0]) < eps) &&
         (fabs(m1[0][1] - m2[0][1]) < eps) &&
         (fabs(m1[0][2] - m2[0][2]) < eps) &&
         (fabs(m1[1][0] - m2[1][0]) < eps) &&
         (fabs(m1[1][1] - m2[1][1]) < eps) &&
         (fabs(m1[1][2] - m2[1][2]) < eps) &&
         (fabs(m1[2][0] - m2[2][0]) < eps) &&
         (fabs(m1[2][1] - m2[2][1]) < eps) &&
         (fabs(m1[2][2] - m2[2][2]) < eps);
}

/** compare matrices */
inline bool utest_mat_equal( const Mesquite::SymMatrix3D& m1, const Mesquite::SymMatrix3D& m2, double eps )
{
  return (fabs(m1(0,0) - m2(0,0)) < eps) &&
         (fabs(m1(0,1) - m2(0,1)) < eps) &&
         (fabs(m1(0,2) - m2(0,2)) < eps) &&
         (fabs(m1(1,1) - m2(1,1)) < eps) &&
         (fabs(m1(1,2) - m2(1,2)) < eps) &&
         (fabs(m1(2,2) - m2(2,2)) < eps);
}

template <unsigned R, unsigned C>
inline std::string msq_mat_str( const Mesquite::MsqMatrix<R,C>& m )
{
  std::ostringstream os;
  for (unsigned i = 0; i < R; ++i) {
    os << "[" << m(i,0);
    for (unsigned j = 1; j < C; ++j)
      os << ", " << m(i,j);
    os << "]";
  }
  return os.str();
}

template <unsigned R, unsigned C>
inline CppUnit::Message ident_check_msg( const Mesquite::MsqMatrix<R,C>& m )
{
  CppUnit::Message mes( "Identity Assertion Failed" );
  mes.addDetail( std::string("Actual: ") + msq_mat_str(m) );
  return mes;
}

template <unsigned R, unsigned C>
inline bool ident_check( const Mesquite::MsqMatrix<R,C>& m )
{
  for (unsigned i = 0; i < R; ++i)
    for (unsigned j = 0; j < C; ++j)
      if (i == j && fabs(m(i,j) - 1.0) > 1e-6)
        return false;
      else if (i != j && fabs(m(i,j)) > 1e-6)
        return false;
  return true;
}

template <unsigned R, unsigned C>
inline CppUnit::Message mat_equal_check_msg( const Mesquite::MsqMatrix<R,C>& A,
                                             const Mesquite::MsqMatrix<R,C>& B )
{
  CppUnit::Message mes( "Matrix Equality Assertion Failed" );
  mes.addDetail( std::string( "Expected: ") + msq_mat_str(A) );
  mes.addDetail( std::string( "Actual:   ") + msq_mat_str(B) );
  return mes;
}

template <unsigned R, unsigned C>
inline CppUnit::Message mat_not_equal_check_msg( const Mesquite::MsqMatrix<R,C>& A,
                                             const Mesquite::MsqMatrix<R,C>& B )
{
  CppUnit::Message mes( "Matrix Inequality Assertion Failed" );
  mes.addDetail( std::string( "Expected: ") + msq_mat_str(A) );
  mes.addDetail( std::string( "Actual:   ") + msq_mat_str(B) );
  return mes;
}

template <unsigned R, unsigned C>
inline bool mat_equal_check( const Mesquite::MsqMatrix<R,C>& A, 
                             const Mesquite::MsqMatrix<R,C>& B, 
                             double eps  )
{
  for (unsigned i = 0; i < R; ++i)
    for (unsigned j = 0; j < C; ++j)
      if (fabs(A(i,j)-B(i,j)) > eps)
        return false;
  return true;
}

template <typename T1, typename T2>
inline bool arrays_equal( const T1* A, const T2* B, size_t len )
{
  for (size_t i = 0; i < len; ++i)
    if (A[i] != B[i])
      return false;
  return true;
}

template <typename T1, typename T2>
inline CppUnit::Message arrays_not_equal_msg( const T1* A, size_t A_len,
                                              const T2* B, size_t B_len )
{
  CppUnit::Message mes( "Equality Assertion Failed for Arrays" );
  
  std::ostringstream strA;
  if (A_len == 0)
    strA << "(empty)";
  else {
    strA << '[' << A[0];
    for (size_t i = 1; i < A_len; ++i)
      strA << ',' << A[i];
    strA << ']';
  }
  mes.addDetail( std::string( "Expected: ") + strA.str() );
  
  std::ostringstream strB;
  if (B_len == 0)
    strB << "(empty)";
  else {
    strB << '[' << B[0];
    for (size_t i = 1; i < B_len; ++i)
      strB << ',' << B[i];
    strB << ']';
  }
  mes.addDetail( std::string( "Actual: ") + strB.str() );
  
  return mes;
}

#endif
