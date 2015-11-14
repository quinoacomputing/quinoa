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


/** \file SphericalDomainTest.cpp
 *  \brief Unit tests for SphericalDomain class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "UnitUtil.hpp"
#include "Mesquite_SphericalDomain.hpp"
#include "Mesquite_ArrayMesh.hpp"

using namespace Mesquite;

class SphericalDomainTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(SphericalDomainTest);
  CPPUNIT_TEST (test_construct);
  CPPUNIT_TEST (test_fit_vertices);
  CPPUNIT_TEST (test_snap_to);
  CPPUNIT_TEST (test_normal_at);
  CPPUNIT_TEST (test_closest_point);
  CPPUNIT_TEST (test_domain_DoF);
  CPPUNIT_TEST_SUITE_END();

public:

  void test_construct();
  void test_fit_vertices();
  void test_snap_to();
  void test_normal_at();
  void test_closest_point();
  void test_domain_DoF();

private:

  void check_closest_pt( const SphericalDomain& dom, 
                         const Vector3D& input_pt,
                         const Vector3D& output_pt );

  void check_normal( const SphericalDomain& dom, 
                     const Vector3D& point,
                     const Vector3D& normal );

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(SphericalDomainTest, "SphericalDomainTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(SphericalDomainTest, "Unit");

void SphericalDomainTest::check_closest_pt( const SphericalDomain& dom, 
                                            const Vector3D& input_pt,
                                            const Vector3D& output_pt )
{
  Vector3D vo = output_pt - dom.center();
  Vector3D vi = input_pt  - dom.center();
  CPPUNIT_ASSERT_DOUBLES_EQUAL( dom.radius(), vo.length(), 1e-6 );
  vi *= dom.radius() / vi.length();
  CPPUNIT_ASSERT_VECTORS_EQUAL( vo, vi, 1e-6 );
}
void SphericalDomainTest::check_normal( const SphericalDomain& dom, 
                                        const Vector3D& point,
                                        const Vector3D& normal )
{
  Vector3D vi = point - dom.center();
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, normal.length(), 1e-6 );
  vi /= vi.length();
  CPPUNIT_ASSERT_VECTORS_EQUAL( vi, normal, 1e-6 );
}

void SphericalDomainTest::test_construct()
{
  Vector3D cen( 1, 2, -1 );
  double rad = 3.14159;
  SphericalDomain sph( cen, rad );
  CPPUNIT_ASSERT_VECTORS_EQUAL( cen, sph.center(), 1e-18 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( rad, sph.radius(), 1e-18 );
  
  cen = Vector3D( 5, 6, 1.14 );
  rad = 1/rad;
  sph.set_sphere( cen, rad );
  CPPUNIT_ASSERT_VECTORS_EQUAL( cen, sph.center(), 1e-18 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( rad, sph.radius(), 1e-18 );
}

static void sphere_point( Vector3D cen, double rad, double theta, double phi, double pt[3] )
{
  Vector3D pdir( cos(phi)*cos(theta), cos(phi)*sin(theta), sin(phi) );
  Vector3D rslt = cen + rad*pdir;
  pt[0] = rslt[0];
  pt[1] = rslt[1];
  pt[2] = rslt[2];
}
  

void SphericalDomainTest::test_fit_vertices()
{
  Vector3D cen( 1, 2, 3 );
  double rad = 2.5;
  
    // point locations on a expected sphere
  const int num_pt = 6;
  double coords[3*num_pt];
  sphere_point( cen, rad, 0.1, 0.2, coords     );
  sphere_point( cen, rad, 2.0, 2.0, coords + 3 );
  sphere_point( cen, rad,-1.0, 0.0, coords + 6 );
  sphere_point( cen, rad, 3.1,-0.5, coords + 9 );
  sphere_point( cen, rad,-1.5,-1.0, coords +12 );
  sphere_point( cen, rad, 0.2, 0.1, coords +15 );
  
    // make sure our setup is valid
  for (int i = 0; i < num_pt; ++i) {
    Vector3D pt(coords + 3*i);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( rad, (pt-cen).length(), 1e-6 );
  }
  
  std::vector<int> fixed(num_pt, 0);
  ArrayMesh mesh( 3, num_pt, coords, arrptr(fixed), 0, TRIANGLE, 0 );
  SphericalDomain sph;
  MsqError err;
  sph.fit_vertices( &mesh, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( cen, sph.center(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( rad, sph.radius(), 1e-6 );
}

void SphericalDomainTest::test_snap_to()
{
  Vector3D cen( 3.14, 15, 0.91 );
  double rad = 6.02;
  SphericalDomain dom(cen, rad);
  
  const int num_pts = 5;
  double points[num_pts][3] = { { 0, 0, 0 },
                                { 10, 11, 8 },
                                { 1, 2, 3 },
                                { -5, 1, 1 },
                                { -1, 0, -2 } };
  for (int i = 0; i < num_pts; ++i) {
    Vector3D v(points[i]);
    dom.snap_to( 0, v );
    check_closest_pt( dom, Vector3D(points[i]), v );
  }
}

void SphericalDomainTest::test_normal_at()
{
  Vector3D cen( -3.14, 0, 0.91 );
  double rad = 2;
  SphericalDomain dom(cen, rad);
  
  const int num_pts = 5;
  double points[num_pts][3] = { { 0, 0, 0 },
                                { 10, 11, 8 },
                                { 1, 2, 3 },
                                { -5, 1, 1 },
                                { -1, 0, -2 } };
  for (int i = 0; i < num_pts; ++i) {
    Vector3D v(points[i]);
    dom.vertex_normal_at( 0, v );
    check_normal( dom, Vector3D(points[i]), v );
    
    v = Vector3D(points[i]);
    dom.element_normal_at( 0, v );
    check_normal( dom, Vector3D(points[i]), v );
  }
}

void SphericalDomainTest::test_closest_point()
{
  Vector3D cen( -1, -1, -2 );
  double rad = 1.4;
  SphericalDomain dom(cen, rad);
  MsqPrintError err(std::cout);
  
  const int num_pts = 5;
  double points[num_pts][3] = { { 0, 0, 0 },
                                { 10, 11, 8 },
                                { 1, 2, 3 },
                                { -5, 1, 1 },
                                { -1, 0, -2 } };
  for (int i = 0; i < num_pts; ++i) {
    Vector3D p(points[i]), c, n;
    dom.closest_point( 0, p, c, n, err );
    ASSERT_NO_ERROR(err);
    check_closest_pt( dom, p, c );
    check_normal( dom, p, n );
  }
}

void SphericalDomainTest::test_domain_DoF()
{
  std::vector<Mesh::VertexHandle> junk(10);
  std::vector<unsigned short> dof(junk.size());
  std::vector<unsigned short> expected(dof.size(), 2);
  SphericalDomain dom;
  MsqPrintError err(std::cout);
  dom.domain_DoF( arrptr(junk), arrptr(dof), junk.size(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT( expected == dof );  
}
