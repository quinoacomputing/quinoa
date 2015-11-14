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


/** \file ConicDomainTest.cpp
 *  \brief Unit tests for ConicDomain class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "UnitUtil.hpp"
#include "Mesquite_ConicDomain.hpp"

using namespace Mesquite;

class ConicDomainTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(ConicDomainTest);
  CPPUNIT_TEST (test_construct);
  CPPUNIT_TEST (test_snap_to);
  CPPUNIT_TEST (test_normal_at);
  CPPUNIT_TEST (test_closest_point);
  CPPUNIT_TEST (test_domain_DoF);
  CPPUNIT_TEST_SUITE_END();

public:

  void test_construct();
  void test_snap_to();
  void test_normal_at();
  void test_closest_point();
  void test_domain_DoF();
   
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ConicDomainTest, "ConicDomainTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ConicDomainTest, "Unit");

void ConicDomainTest::test_construct()
{
  double rad = 1.5;
  double hei = 20.1;
  Vector3D axis( 1, 2, 3 );
  Vector3D point( -1, -1, 1 );
  ConicDomain dom( rad, hei, axis, point );
  
  axis /= axis.length();
  CPPUNIT_ASSERT_VECTORS_EQUAL( axis, dom.axis(), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( point, dom.point(), 1e-18 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( rad, dom.point_radius(), 1e-18 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( hei, dom.height_from_point(), 1e-18 );
}

void ConicDomainTest::test_snap_to()
{
  const double a = 3.0, b = 4.0;
  const double f = a*b/(a*a + b*b);
  ConicDomain cone( a, b );
  
    // test some points
  Vector3D pt( b, 0, a );
  Vector3D close( pt );
  cone.snap_to( 0, close );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt * f, close, 1e-6 );
  
  pt = Vector3D( 0, b, a );
  close = pt;
  cone.snap_to( 0, close );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt * f, close, 1e-6 );
  
  pt = Vector3D( 0, -b, a );
  close = pt;
  cone.snap_to( 0, close );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt * f, close, 1e-6 );
  
  pt = Vector3D( -b, 0, a );
  close = pt;
  cone.snap_to( 0, close );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt * f, close, 1e-6 );
  
  // test point at apex
  pt = cone.point() + cone.height_from_point()*cone.axis();
  close = pt;
  cone.snap_to( 0, close );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt, close, 1e-6 );  
}

void ConicDomainTest::test_normal_at()
{
  const double a = 3.0, b = 4.0;
  ConicDomain cone( a, b );
  
  Vector3D pt( b, 0, a );
  Vector3D norm( pt );
  cone.vertex_normal_at( 0, norm );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt/pt.length(), norm, 1e-6 );
  norm = pt;
  cone.element_normal_at( 0, norm );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt/pt.length(), norm, 1e-6 );
  
  pt = Vector3D( 0, b, a );
  norm = pt;
  cone.vertex_normal_at( 0, norm );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt/pt.length(), norm, 1e-6 );
  norm = pt;
  cone.element_normal_at( 0, norm );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt/pt.length(), norm, 1e-6 );
  
  pt = Vector3D( 0, -b, a );
  norm = pt;
  cone.vertex_normal_at( 0, norm );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt/pt.length(), norm, 1e-6 );
  norm = pt;
  cone.element_normal_at( 0, norm );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt/pt.length(), norm, 1e-6 );
  
  pt = Vector3D( -b, 0, a );
  norm = pt;
  cone.vertex_normal_at( 0, norm );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt/pt.length(), norm, 1e-6 );
  norm = pt;
  cone.element_normal_at( 0, norm );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt/pt.length(), norm, 1e-6 );
}

void ConicDomainTest::test_closest_point()
{
  const double a = 3.0, b = 4.0;
  const double f = a*b/(a*a + b*b);
  ConicDomain cone( a, b );
  MsqError err;
  
  Vector3D pt( b, 0, a ), close, norm;
  cone.closest_point( 0, pt, close, norm, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt * f, close, 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt/pt.length(), norm, 1e-6 );
  
  pt = Vector3D( 0, b, a );
  cone.closest_point( 0, pt, close, norm, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt * f, close, 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt/pt.length(), norm, 1e-6 );
  
  pt = Vector3D( 0, -b, a );
  cone.closest_point( 0, pt, close, norm, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt * f, close, 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt/pt.length(), norm, 1e-6 );
  
  pt = Vector3D( -b, 0, a );
  cone.closest_point( 0, pt, close, norm, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt * f, close, 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( pt/pt.length(), norm, 1e-6 );
}

void ConicDomainTest::test_domain_DoF()
{
  std::vector<Mesh::VertexHandle> junk(10);
  std::vector<unsigned short> dof(junk.size());
  std::vector<unsigned short> expected(dof.size(), 2);
  ConicDomain dom;
  MsqPrintError err(std::cout);
  dom.domain_DoF( arrptr(junk), arrptr(dof), junk.size(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT( expected == dof );  
}
