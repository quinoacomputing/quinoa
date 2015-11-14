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


/** \file CircleDomainTest.cpp
 *  \brief UnitTests for CircleDomain class
 *  \author Jason Kraftcheck 
 */

#include "UnitUtil.hpp"
#include "Mesquite_MeshDomain1D.hpp"

namespace MESQUITE_NS {

class CircleDomainTest: public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CircleDomainTest);
  CPPUNIT_TEST (test_snap_to);
  CPPUNIT_TEST (test_arc_length);
  CPPUNIT_TEST (test_position_from_length);
  CPPUNIT_TEST_SUITE_END();
 
public:

  void test_snap_to();
  void test_arc_length();
  void test_position_from_length();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(CircleDomainTest, "CircleDomainTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(CircleDomainTest, "Unit");

void CircleDomainTest::test_snap_to()
{
  Vector3D origin(0,0,0);
  Vector3D z(0,0,1);
  double rad1 = 4.0/3.0;
  CircleDomain dom1( origin, z, rad1 );
  Vector3D pt( 1, 0, 0 );
  dom1.snap_to( 0, pt );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(rad1,0,0), pt, 1e-6 );
  Vector3D a = Vector3D( 1, 2, 3 );
  pt = a;
  dom1.snap_to( 0, pt );
  a = Vector3D( 1, 2, 0 );
  a *= rad1 / a.length();
  CPPUNIT_ASSERT_VECTORS_EQUAL( a, pt, 1e-6 );
  
  Vector3D some_pt(5,-1,6);
  Vector3D some_dir(-5,-4,1);
  double rad2 = 1.0;
  CircleDomain dom2( some_pt, some_dir, rad2 );
  
  a = Vector3D( 0, 0, 0);
  pt = a;
  dom2.snap_to( 0, pt );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( rad2, (pt - some_pt).length(), 1e-6 ); // rad from center
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, (pt - some_pt) % some_dir, 1e-6 );// in plane
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, ((pt - some_pt) * (a - some_pt)) % some_dir, 1e-6 ); // correct direction from center

  a = Vector3D( 0, -1, -2 );
  pt = a;
  dom2.snap_to( 0, pt );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( rad2, (pt - some_pt).length(), 1e-6 ); // rad from center
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, (pt - some_pt) % some_dir, 1e-6 );// in plane
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, ((pt - some_pt) * (a - some_pt)) % some_dir, 1e-6 ); // correct direction from center
}

void CircleDomainTest::test_arc_length()
{
  MsqPrintError err(std::cerr);
  
  Vector3D origin(0,0,0);
  Vector3D z(0,0,1);
  double rad1 = 4.0/3.0;
  CircleDomain dom1( origin, z, rad1 );
  
  Vector3D xp = Vector3D(1,0,0);
  Vector3D xn = Vector3D(-1,0,0);
  Vector3D yp = Vector3D(0,1,0);
  
  double len = dom1.arc_length( xp.to_array(), yp.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( rad1 * M_PI * 0.5, len, 1e-6 );
  
  len = dom1.arc_length( xp.to_array(), xn.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( rad1 * M_PI, len, 1e-6 );
  
  len = dom1.arc_length( yp.to_array(), xp.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( -rad1 * M_PI * 0.5, len, 1e-6 );
  
  len = dom1.arc_length( (xp+z).to_array(), yp.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( rad1 * M_PI * 0.5, len, 1e-6 );
  
  len = dom1.arc_length( xp.to_array(), (yp+z).to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( rad1 * M_PI * 0.5, len, 1e-6 );
  
  Vector3D center(-1,-2,-1);
  Vector3D normal(-2,-1,-2);
  double rad = 1.5;
  CircleDomain dom2( center, normal, rad );
  
  Vector3D v(1,0,0);
  v = v * normal;
  v *= rad / v.length();
  v += center;
  len = dom2.arc_length( v.to_array(), v.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, len, 1e-6 );
}

void CircleDomainTest::test_position_from_length()
{
  MsqPrintError err(std::cerr);
  
  Vector3D origin(0,0,0);
  Vector3D z(0,0,1);
  double rad1 = 4.0/3.0;
  CircleDomain dom1( origin, z, rad1 );
  
  Vector3D xp = Vector3D(1,0,0);
  Vector3D xn = Vector3D(-1,0,0);
  Vector3D yp = Vector3D(0,1,0);
  
  const double qc = 0.5 * rad1 * M_PI;
  
  Vector3D result;
  dom1.position_from_length( xp.to_array(), qc, result.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( yp*rad1, result, 1e-6 );
  
  dom1.position_from_length( xp.to_array(), 2*qc, result.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( xn*rad1, result, 1e-6 );
  
  dom1.position_from_length( yp.to_array(), -qc, result.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( xp*rad1, result, 1e-6 );
  
  dom1.position_from_length( (xp+z).to_array(), qc, result.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( yp*rad1, result, 1e-6 );
  
  Vector3D center(-1,-2,-1);
  Vector3D normal(-2,-1,-2);
  double rad2 = 1.5;
  CircleDomain dom2( center, normal, rad2 );
  
  Vector3D v(1,0,0);
  v = v * normal;
  v *= rad2 / v.length();
  v += center;
  dom2.position_from_length( v.to_array(), 0.0, result.to_array(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( v, result, 1e-6 );
}

} // namespace MESQUITE_NS
