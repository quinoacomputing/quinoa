/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file LineDomainTest.cpp
 *  \brief UnitTests for LineDomain class
 *  \author Jason Kraftcheck 
 */

#include "UnitUtil.hpp"
#include "Mesquite_MeshDomain1D.hpp"

namespace MESQUITE_NS {

class LineDomainTest: public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(LineDomainTest);
  CPPUNIT_TEST (test_snap_to);
  CPPUNIT_TEST (test_arc_length);
  CPPUNIT_TEST (test_position_from_length);
  CPPUNIT_TEST_SUITE_END();
 
public:

  void test_snap_to();
  void test_arc_length();
  void test_position_from_length();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(LineDomainTest, "LineDomainTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(LineDomainTest, "Unit");

void LineDomainTest::test_snap_to()
{
  const Vector3D base1(1,1,1);
  const Vector3D dir1(-1,-2,-3);
  LineDomain dom1( base1, dir1 );

  Vector3D point(0,0,0), result(point);
  dom1.snap_to( 0, result );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, ((result - base1) * dir1).length(), 1e-6 ); // on the line
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, ((point - result) % dir1), 1e-6 ); // moved perp to line
  
  point = Vector3D(10,11,12);
  result = point;
  dom1.snap_to( 0, result );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, ((result - base1) * dir1).length(), 1e-6 ); // on the line
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, ((point - result) % dir1), 1e-6 ); // moved perp to line
}

void LineDomainTest::test_arc_length()
{
  MsqPrintError err(std::cerr);
  
  const Vector3D base1(1,1,1);
  const Vector3D dir1(-1,-2,-3);
  LineDomain dom1( base1, dir1 );

  double l1 = 1.0, l2 = M_PI;
  Vector3D p1 = base1 + l1/dir1.length() * dir1;
  Vector3D p2 = base1 + l2/dir1.length() * dir1;
  double len = dom1.arc_length( p1.to_array(), p2.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( l2 - l1, len, 1e-6 );
  len = dom1.arc_length( p2.to_array(), p1.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( l1 - l2, len, 1e-6 );
}

void LineDomainTest::test_position_from_length()
{
  MsqPrintError err(std::cerr);
  
  const Vector3D base1(1,1,1);
  const Vector3D dir1(-1,-2,-3);
  LineDomain dom1( base1, dir1 );
  const Vector3D unit = dir1 / dir1.length();
  double l1 = 1.5;
  Vector3D result;
  dom1.position_from_length( base1.to_array(), l1, result.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( base1 + l1 * unit, result, 1e-6 );
  
  double l2 = 3.333;
  Vector3D result2;
  dom1.position_from_length( result.to_array(), l2, result2.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( result + l2 * unit, result2, 1e-6 );
}

} // namespace MESQUITE_NS
