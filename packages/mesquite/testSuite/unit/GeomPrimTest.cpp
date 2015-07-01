/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2008 Sandia National Laboratories.  Developed at the
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


/** \file GeomPrimTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "MsqGeomPrim.hpp"
#include "UnitUtil.hpp"

using namespace Mesquite;
using namespace std;


class GeomPrimTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(GeomPrimTest);
  CPPUNIT_TEST (test_line_basic);
  CPPUNIT_TEST (test_line_closest_to_point);
  CPPUNIT_TEST (test_line_intersect);
  CPPUNIT_TEST (test_line_closest_to_line);
  CPPUNIT_TEST (test_circle_basic);
  CPPUNIT_TEST (test_circle_from_three_points);
  CPPUNIT_TEST (test_circle_from_two_points);
  CPPUNIT_TEST (test_circle_closest_to_point);
  CPPUNIT_TEST (test_circle_closest_with_tangent);
  CPPUNIT_TEST (test_plane_basic);
  CPPUNIT_TEST (test_plane_distance);
  CPPUNIT_TEST (test_plane_closest_to_point);
  CPPUNIT_TEST (test_plane_intersect_plane);
  CPPUNIT_TEST (test_plane_intersect_line);
  CPPUNIT_TEST (test_sphere_basic);
  CPPUNIT_TEST (test_sphere_closest_to_point);
  CPPUNIT_TEST (test_sphere_intersect_plane);
  CPPUNIT_TEST (test_sphere_intersect_sphere);
  CPPUNIT_TEST_SUITE_END();
  
public:

  void test_line_basic();
  void test_line_closest_to_point();
  void test_line_intersect();
  void test_line_closest_to_line();
  void test_circle_basic();
  void test_circle_from_three_points();
  void test_circle_from_two_points();
  void test_circle_closest_to_point();
  void test_circle_closest_with_tangent();
  void test_plane_basic();
  void test_plane_distance();
  void test_plane_closest_to_point();
  void test_plane_intersect_plane();
  void test_plane_intersect_line();
  void test_sphere_basic();
  void test_sphere_closest_to_point();
  void test_sphere_intersect_plane();
  void test_sphere_intersect_sphere();

};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(GeomPrimTest, "Unit");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(GeomPrimTest, "GeomPrimTest");

void GeomPrimTest::test_line_basic()
{
  const Vector3D point(1,2,3), direction(-1,-2,-3);
  const MsqLine line( point, direction );
  CPPUNIT_ASSERT_VECTORS_EQUAL( point, line.point(), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, line.direction().length(), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( direction/direction.length(), line.direction(), 1e-6 );
  double param = direction.length();
  CPPUNIT_ASSERT_VECTORS_EQUAL( point+direction, line.point(param), 1e-6 );

  const MsqLine line2( MsqLine::two_point( point, point+direction ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( line.point(), line2.point(), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( line.direction(), line2.direction(), 1e-6 );
}

void GeomPrimTest::test_line_closest_to_point()
{
  const MsqLine line( Vector3D(1,2,3), Vector3D(3,2,1) );

  const Vector3D p1( 0, 0, 0), p2( 1, 1, 1 ), p3( -5, 0, 0 );
  
  double p = line.closest( p1 );
  Vector3D diff = line.point(p) - p1;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, diff % line.direction(), 1e-6 );
  
  p = line.closest( p2 );
  diff = line.point(p) - p2;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, diff % line.direction(), 1e-6 );
  
  p = line.closest( p3 );
  diff = line.point(p) - p3;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, diff % line.direction(), 1e-6 );
}

void GeomPrimTest::test_line_intersect()
{
  const Vector3D xsect( 2, 0, 0 );
  const MsqLine line1( Vector3D(0,0,0), Vector3D(1,0,0) );
  const MsqLine line2( Vector3D(2,-1,0), Vector3D(0,1,0) );
  double param;
  
  CPPUNIT_ASSERT( line1.intersect(line2, param, 1e-6) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( xsect, line1.point(param), 1e-6 );
  CPPUNIT_ASSERT( line2.intersect(line1, param, 1e-6) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( xsect, line2.point(param), 1e-6 );
}
  
void GeomPrimTest::test_line_closest_to_line()
{
  const MsqLine line1( Vector3D( 3, 2, 1), Vector3D(5, 6, 5) );
  const MsqLine line2( Vector3D( -1, -2, -3), Vector3D(1, 0, 0) );
  double param1, param2;
  Vector3D diff;
  
  CPPUNIT_ASSERT( line1.closest( line2, param1 ) );
  CPPUNIT_ASSERT( line2.closest( line1, param2 ) );
  diff = line1.point(param1) - line2.point(param2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, diff % line1.direction(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, diff % line2.direction(), 1e-6 );
}


void GeomPrimTest::test_circle_basic()
{
  const Vector3D center(0,0,2), normal(1,0,0);
  const double radius = 1.5;
  const MsqCircle circle( center, normal, radius );
  CPPUNIT_ASSERT_VECTORS_EQUAL( center, circle.center(), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( normal, circle.normal(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( radius, circle.radius(), 1e-6 );
}


void GeomPrimTest::test_circle_from_three_points()
{
  const Vector3D p1(1,1,1), p2(2,2,2), p3(3,0,0);
  MsqCircle circ;
  CPPUNIT_ASSERT( MsqCircle::three_point(p1,p2,p3,circ) );
  const Vector3D r1(p1 - circ.center());
  const Vector3D r2(p2 - circ.center());
  const Vector3D r3(p3 - circ.center());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( circ.radius(), r1.length(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( circ.radius(), r2.length(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( circ.radius(), r3.length(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, r1 % circ.normal(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, r2 % circ.normal(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, r3 % circ.normal(), 1e-6 );
}

void GeomPrimTest::test_circle_from_two_points() 
{
  const Vector3D center(1,2,3);
  const Vector3D normal(1,1,1);
  const Vector3D somedir(0,0,1);
  const double radius = 0.33;
  Vector3D p1dir = normal * somedir;
  Vector3D p2dir = p1dir + normal * p1dir;
  p1dir *= radius / p1dir.length();
  p2dir *= radius / p2dir.length();
  
  MsqCircle circ;
  const Vector3D p1(center + p1dir);
  const Vector3D p2(center + p2dir);
  CPPUNIT_ASSERT( MsqCircle::two_point( center, p1, p2, circ ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( center, circ.center(), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( normal / normal.length(), circ.normal(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( radius, circ.radius(), 1e-6 );
}

void GeomPrimTest::test_circle_closest_to_point()
{
  const MsqCircle circ( Vector3D( 1, 1, 2 ), Vector3D( 0, 1, 0 ), 0.5 );
  Vector3D v1 = circ.radial_vector();
  Vector3D v2 = circ.normal() * v1;
  v1 /= v1.length();
  v2 /= v2.length();
  const Vector3D q1 = circ.center() + 0.2 * v1;
  const Vector3D r1 = circ.center() + circ.radius() * v1;
  const Vector3D q2 = circ.center() + 100 * v2;
  const Vector3D r2 = circ.center() + circ.radius() * v2;
  
  CPPUNIT_ASSERT_VECTORS_EQUAL( r1, circ.closest(q1), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( r2, circ.closest(q2), 1e-6 );
}

void GeomPrimTest::test_circle_closest_with_tangent()
{
  const Vector3D center( 1, 3, 2 );
  const Vector3D x(1,0,0), y(0,1,0), z(0, 0, 1);
  const double radius = 2.81;
  const MsqCircle circ( center, z, radius);

  Vector3D in, out, tan;

    // at x minimum
  in = center; in[0] -= radius;
  CPPUNIT_ASSERT( circ.closest( in, out, tan ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( in, out, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, tan % z, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, tan % x, 1e-6 );

    // at x maximum
  in = center; in[0] += radius;
  CPPUNIT_ASSERT( circ.closest( in, out, tan ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( in, out, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, tan % z, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, tan % x, 1e-6 );

    // at y minimum
  in = center; in[1] -= radius;
  CPPUNIT_ASSERT( circ.closest( in, out, tan ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( in, out, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, tan % z, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, tan % y, 1e-6 );

    // at y maximum
  in = center; in[1] += radius;
  CPPUNIT_ASSERT( circ.closest( in, out, tan ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( in, out, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, tan % z, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, tan % y, 1e-6 );
}

void GeomPrimTest::test_plane_basic()
{
  const Vector3D norm(1,1,2);
  const double coeff = -1;
  const MsqPlane plane( norm, coeff );
  CPPUNIT_ASSERT_VECTORS_EQUAL( norm/norm.length(),  plane.normal(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( coeff/norm.length(), plane.coefficient(), 1e-6 );
  
  const MsqPlane plane2( norm, plane.point() );
  CPPUNIT_ASSERT_VECTORS_EQUAL( plane.normal(), plane2.normal(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( plane.coefficient(), plane2.coefficient(), 1e-6 );
}

void GeomPrimTest::test_plane_distance()
{
  double offset = 2.3;
  const Vector3D x(1,0,0), y(0,1,0), z(0,0,1);
  const MsqPlane xyplane( z, z*offset );
  const MsqPlane yzplane( x, x*offset );
  const MsqPlane xzplane( y, y*offset );
  
  const Vector3D p1(1,2,3), p2(-1,-5,0), p3(0.5,0.5,100);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( fabs(p1[2] - offset), xyplane.distance(p1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( fabs(p2[2] - offset), xyplane.distance(p2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( fabs(p3[2] - offset), xyplane.distance(p3), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( fabs(p1[0] - offset), yzplane.distance(p1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( fabs(p2[0] - offset), yzplane.distance(p2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( fabs(p3[0] - offset), yzplane.distance(p3), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( fabs(p1[1] - offset), xzplane.distance(p1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( fabs(p2[1] - offset), xzplane.distance(p2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( fabs(p3[1] - offset), xzplane.distance(p3), 1e-6 );
}

void GeomPrimTest::test_plane_closest_to_point()
{
  double offset = 2.3;
  const Vector3D x(1,0,0), y(0,1,0), z(0,0,1);
  const MsqPlane xyplane( z, z*offset );
  const MsqPlane yzplane( x, x*offset );
  const MsqPlane xzplane( y, y*offset );
  
  const Vector3D p1(1,2,3), p2(-1,-5,0), p3(0.5,0.5,100);
  const Vector3D p1xy(p1[0],p1[1],offset);
  
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(p1[0],p1[1],offset), xyplane.closest(p1), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(p2[0],p2[1],offset), xyplane.closest(p2), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(p3[0],p3[1],offset), xyplane.closest(p3), 1e-6 );
  
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(p1[0],offset,p1[2]), xzplane.closest(p1), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(p2[0],offset,p2[2]), xzplane.closest(p2), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(p3[0],offset,p3[2]), xzplane.closest(p3), 1e-6 );
  
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(offset,p1[1],p1[2]), yzplane.closest(p1), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(offset,p2[1],p2[2]), yzplane.closest(p2), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(offset,p3[1],p3[2]), yzplane.closest(p3), 1e-6 );
}

void GeomPrimTest::test_plane_intersect_plane()
{
  const double zoffset = 2.0;
  const double yoffset = -3.0;
  const MsqPlane xy( Vector3D(0,0,1), Vector3D(0,0,zoffset) );
  const MsqPlane xz( Vector3D(0,1,0), Vector3D(0,yoffset,0) );
  MsqLine line;
  
  CPPUNIT_ASSERT( xy.intersect(xz, line) );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, xy.distance(line.point()), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, xz.distance(line.point()), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, line.direction() % xy.normal(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, line.direction() % xz.normal(), 1e-6 );
}

void GeomPrimTest::test_plane_intersect_line()
{
  const Vector3D normal( 1, 1, 1 );
  const double coeff = -3;
  const Vector3D direction( -1, -1, -1 );
  const Vector3D point( 0, 0, 0);
  const MsqPlane plane( normal, coeff );
  const MsqLine line( point, direction );
  
  double param;
  CPPUNIT_ASSERT( plane.intersect( line, param ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(1,1,1), line.point(param), 1e-6 );
}

void GeomPrimTest::test_sphere_basic()
{
  const Vector3D center( -1, 2, 0.5 );
  const double radius = 3.5;
  const MsqSphere sphere( center, radius );
  CPPUNIT_ASSERT_VECTORS_EQUAL( center, sphere.center(), 1e-10 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( radius, sphere.radius(), 1e-10 );
}

void GeomPrimTest::test_sphere_closest_to_point()
{
  const Vector3D center( -1, 2, 0.5 );
  const double radius = 3.5;
  MsqSphere sphere( center, radius );
  
  const Vector3D xpt = center + Vector3D(1,0,0);
  const Vector3D ypt = center + Vector3D(0,-100,0);
  const Vector3D zpt = center + Vector3D(0,0,0.5);
  
  Vector3D closest, normal;
  closest = sphere.closest( xpt );
  CPPUNIT_ASSERT_VECTORS_EQUAL( center + Vector3D(radius,0,0), closest, 1e-6 );
  closest = Vector3D(0,0,0);
  CPPUNIT_ASSERT( sphere.closest( xpt, closest, normal ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( center + Vector3D(radius,0,0), closest, 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(1,0,0), normal, 1e-6 );  
  
  closest = sphere.closest( ypt );
  CPPUNIT_ASSERT_VECTORS_EQUAL( center + Vector3D(0,-radius,0), closest, 1e-6 );
  closest = Vector3D(0,0,0);
  CPPUNIT_ASSERT( sphere.closest( ypt, closest, normal ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( center + Vector3D(0,-radius,0), closest, 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,-1,0), normal, 1e-6 );  
  
  closest = sphere.closest( zpt );
  CPPUNIT_ASSERT_VECTORS_EQUAL( center + Vector3D(0,0,radius), closest, 1e-6 );
  closest = Vector3D(0,0,0);
  CPPUNIT_ASSERT( sphere.closest( zpt, closest, normal ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( center + Vector3D(0,0,radius), closest, 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,0,1), normal, 1e-6 );  
}

void GeomPrimTest::test_sphere_intersect_plane()
{
  const MsqSphere sphere( Vector3D(1,1,1), 2.0 );
  MsqCircle result;
  
  const MsqPlane noxsect( Vector3D(1,1,1), 5 );
  CPPUNIT_ASSERT(!sphere.intersect( noxsect, result ));
  
  const MsqPlane bisect( Vector3D(0,0,1), -1.0);
  CPPUNIT_ASSERT(sphere.intersect( bisect, result ) );
  CPPUNIT_ASSERT_VECTORS_EQUAL( sphere.center(), result.center(), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( bisect.normal(), result.normal(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( sphere.radius(), result.radius(), 1e-6 );
  
  const MsqPlane tangent( Vector3D(0,1,0), 1.0 );
  CPPUNIT_ASSERT(sphere.intersect( tangent, result ) );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, result.radius(), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( sphere.center() - Vector3D( 0, sphere.radius(), 0), 
                               result.center(), 1e-6 );
                               
  const MsqPlane xyz0( Vector3D(1,1,1), 0.0 );
  CPPUNIT_ASSERT(sphere.intersect( xyz0, result ) );
  if (xyz0.normal() % result.normal() > 0.0)
    CPPUNIT_ASSERT_VECTORS_EQUAL( xyz0.normal(), result.normal(), 1e-6 );
  else
    CPPUNIT_ASSERT_VECTORS_EQUAL( -xyz0.normal(), result.normal(), 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,0,0), result.center(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( sqrt(sphere.radius()*sphere.radius()-3), result.radius(), 1e-6 );
}

void GeomPrimTest::test_sphere_intersect_sphere()
{
  const MsqSphere orig( Vector3D(0,0,0), 5.0 );
  const MsqSphere zsph( Vector3D(0,0,3), 3.0 );
  const MsqSphere ysph( Vector3D(0,6,0), 2.0 );
  const MsqSphere noxsct( Vector3D(7,0,0), 1.0 );
  Vector3D p1, p2, v1, v2;
  MsqCircle result;

  CPPUNIT_ASSERT( !orig.intersect( noxsct, result ) );
  
  CPPUNIT_ASSERT( orig.intersect( zsph, result ) );
    // choose two points on resulting circle
  v1 = result.radial_vector();
  v2 = result.normal() * result.radial_vector();
  p1 = result.center() + v1;
  p2 = result.center() + v2;
    // verify points are on both spheres
  CPPUNIT_ASSERT_DOUBLES_EQUAL( (p1 - orig.center()).length(), orig.radius(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( (p1 - zsph.center()).length(), zsph.radius(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( (p2 - orig.center()).length(), orig.radius(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( (p2 - zsph.center()).length(), zsph.radius(), 1e-6 );
  
  CPPUNIT_ASSERT( orig.intersect( ysph, result ) );
    // choose two points on resulting circle
  v1 = result.radial_vector();
  v2 = result.normal() * result.radial_vector();
  p1 = result.center() + v1;
  p2 = result.center() + v2;
    // verify points are on both spheres
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, orig.distance(p1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, ysph.distance(p1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, orig.distance(p2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, ysph.distance(p2), 1e-6 );
}

