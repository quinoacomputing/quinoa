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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file XYRectangleTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Mesquite_XYRectangle.hpp"
#include "Mesquite_ArrayMesh.hpp"
#include "Mesquite_MsqVertex.hpp"
#include "UnitUtil.hpp"
#include "Mesquite_MsqError.hpp"


using namespace Mesquite;
using namespace std;


//const double WIDTH = 2, HEIGHT = 3, XMIN = -2, YMIN = 1;
const double WIDTH = 4, HEIGHT = 2, XMIN = 0, YMIN = 0;


class XYRectangleTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(XYRectangleTest);
  CPPUNIT_TEST (test_snap_to);
  CPPUNIT_TEST (test_normal_at);
  CPPUNIT_TEST (test_closest_point);
  CPPUNIT_TEST (test_domain_DoF);
  CPPUNIT_TEST_SUITE_END();
  
  vector<double> vertCoords,invertCoords;
  vector<int> fixedFlags;
  vector<unsigned long> triConn, invertConn;
  
  ArrayMesh myMesh;
  XYRectangle myDomain;
  std::vector<double> mCoords;
  std::vector<int> mFlags;
  
public:
  void setUp();
  void tearDown();
  
public:
  
  XYRectangleTest() : myDomain( WIDTH, HEIGHT, XMIN, YMIN ) {}

  void test_snap_to();
  void test_normal_at();
  void test_closest_point();
  void test_domain_DoF();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(XYRectangleTest, "Unit");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(XYRectangleTest, "XYRectangleTest");

/*   (7)-------(8)-------(9)
 *    |\       / \       /|
 *    | \  8  /   \  9  / |
 *    |  \   /  6  \   /  |
 *    | 5 \ /       \ / 7 |
 *   (3)--(4)-------(5)--(6)
 *    | 2 / \       / \ 4 |
 *    |  /   \  3  /   \  |
 *    | /  0  \   /  1  \ |
 *    |/       \ /       \|
 *   (0)-------(1)-------(2)
 */
double TEST_MESH_COORDS[] = { 
  XMIN            , YMIN             , 0,
  XMIN + 0.5*WIDTH, YMIN             , 0,
  XMIN +     WIDTH, YMIN             , 0,
  XMIN            , YMIN + 0.5*HEIGHT, 0,
  XMIN + .25*WIDTH, YMIN + 0.5*HEIGHT, 0,
  XMIN + .75*WIDTH, YMIN + 0.5*HEIGHT, 0,
  XMIN +     WIDTH, YMIN + 0.5*HEIGHT, 0,
  XMIN            , YMIN +     HEIGHT, 0,
  XMIN + 0.5*WIDTH, YMIN +     HEIGHT, 0,
  XMIN +     WIDTH, YMIN +     HEIGHT, 0 };
const int NUM_TEST_MESH_VERTS = sizeof(TEST_MESH_COORDS)/sizeof(double)/3;
const unsigned long TEST_MESH_CONN[] = {
  0, 1, 4,
  1, 2, 5,
  0, 4, 3,
  1, 5, 4,
  2, 6, 5, 
  3, 4, 7,
  4, 5, 8,
  5, 6, 9,
  4, 8, 7,
  5, 9, 8 };
const int NUM_TEST_MESH_TRIS = sizeof(TEST_MESH_CONN)/sizeof(TEST_MESH_CONN[0])/3;

void XYRectangleTest::setUp()
{
  MsqError err;
  mCoords.resize( 3*NUM_TEST_MESH_VERTS );
  std::copy( TEST_MESH_COORDS, TEST_MESH_COORDS + 3*NUM_TEST_MESH_VERTS, mCoords.begin() );
  mFlags.clear();
  mFlags.resize( NUM_TEST_MESH_VERTS, 0 );
  myMesh.set_mesh( 3, NUM_TEST_MESH_VERTS, arrptr(mCoords), arrptr(mFlags),
                   NUM_TEST_MESH_TRIS, TRIANGLE, TEST_MESH_CONN );
  myDomain.setup( &myMesh, err );
  ASSERT_NO_ERROR(err);
}

void XYRectangleTest::tearDown()
{
}

Vector3D snap_to( const Vector3D& vertex, const Vector3D& point )
{
  Vector3D result;
  result[2] = 0.0;
  
  if (fabs(vertex[0] - XMIN) < 1e-6)
    result[0] = XMIN;
  else if (fabs(vertex[0] - XMIN - WIDTH) < 1e-6)
    result[0] = XMIN + WIDTH;
  else
    result[0] = point[0];
  
  if (fabs(vertex[1] - YMIN) < 1e-6)
    result[1] = YMIN;
  else if (fabs(vertex[1] - YMIN - HEIGHT) < 1e-6)
    result[1] = YMIN + HEIGHT;
  else
    result[1] = point[1];
  
  return result;
}

void XYRectangleTest::test_snap_to()
{
  MsqError err;
  Vector3D off, exp, act;

  const Vector3D d1( 0.1, 0.1, 0.0 );
  const Vector3D d2( 0.0, 0.2, 0.0 );
  const Vector3D d3( 0.3, 0.0, 0.0 );
  const Vector3D d4( 0.0, 0.0, 5.0 );
  const Vector3D d5( 0.5, 0.5, 0.5 );

  std::vector<Mesh::VertexHandle> verts;
  myMesh.get_all_vertices( verts, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!verts.empty());
  std::vector<MsqVertex> coords( verts.size() );
  myMesh.vertices_get_coordinates( arrptr(verts), arrptr(coords), verts.size(), err );
  ASSERT_NO_ERROR(err);
  
  for (size_t i = 0; i < coords.size(); ++i) {

    off = coords[i] + d1;
    exp = snap_to( coords[i], off );
    act = off;
    myDomain.snap_to( verts[i], act );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp, act, 1e-6 );

    off = coords[i] + d2;
    exp = snap_to( coords[i], off );
    act = off;
    myDomain.snap_to( verts[i], act );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp, act, 1e-6 );

    off = coords[i] + d3;
    exp = snap_to( coords[i], off );
    act = off;
    myDomain.snap_to( verts[i], act );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp, act, 1e-6 );

    off = coords[i] + d4;
    exp = snap_to( coords[i], off );
    act = off;
    myDomain.snap_to( verts[i], act );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp, act, 1e-6 );

    off = coords[i] + d5;
    exp = snap_to( coords[i], off );
    act = off;
    myDomain.snap_to( verts[i], act );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp, act, 1e-6 );
  }
}

void XYRectangleTest::test_normal_at()
{
  MsqError err;
  std::vector<Mesh::VertexHandle> vertices;
  myMesh.get_all_vertices( vertices, err );
  ASSERT_NO_ERROR(err);
  
  std::vector<MsqVertex> coords(vertices.size());
  myMesh.vertices_get_coordinates( arrptr(vertices), arrptr(coords), vertices.size(), err );
  ASSERT_NO_ERROR(err);
  
  std::vector<Vector3D> normals(vertices.size());
  std::copy( coords.begin(), coords.end(), normals.begin() );
  myDomain.vertex_normal_at( arrptr(vertices), arrptr(normals), vertices.size(), err );
  ASSERT_NO_ERROR(err);
  
  for (size_t i = 0; i < normals.size(); ++i) {
    CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,0,1), normals[i], 1e-6 );
  }    
}

void XYRectangleTest::test_closest_point()
{
  MsqError err;
  std::vector<Mesh::ElementHandle> elems;
  myMesh.get_all_elements( elems, err );
  ASSERT_NO_ERROR(err);
  
  for (size_t i = 0; i < elems.size(); ++i) {
    std::vector<Mesh::VertexHandle> verts;
    std::vector<size_t> junk;
    MsqVertex coords;
    myMesh.elements_get_attached_vertices( &elems[i], 1, verts, junk, err );
    ASSERT_NO_ERROR(err);
    myMesh.vertices_get_coordinates( arrptr(verts), &coords, 1, err );
    ASSERT_NO_ERROR(err);
    
    Vector3D offset(coords + Vector3D(0,0,3)), closest, norm;
    myDomain.closest_point( elems[i], offset, closest, norm, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_VECTORS_EQUAL( coords, closest, 1e-6 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,0,1), norm, 1e-6 );
  }    
}

unsigned short dof( const Vector3D& point )
{
  unsigned short result = 2;
  if ((fabs(point[0] - XMIN) < 1e-6) || (fabs(point[0] - XMIN - WIDTH) < 1e-6))
    --result;
  if ((fabs(point[1] - YMIN) < 1e-6) || (fabs(point[1] - YMIN - HEIGHT) < 1e-6))
    --result;
  return result;
}

void XYRectangleTest::test_domain_DoF()
{
  MsqError err;

  std::vector<Mesh::VertexHandle> verts;
  myMesh.get_all_vertices( verts, err );
  ASSERT_NO_ERROR(err);
  std::vector<MsqVertex> coords( verts.size() );
  myMesh.vertices_get_coordinates( arrptr(verts), arrptr(coords), verts.size(), err );
  ASSERT_NO_ERROR(err);
  
  for (size_t i = 0; i < coords.size(); ++i) {
    unsigned short exp = dof( coords[i] );
    unsigned short act = 100;
    myDomain.domain_DoF( &verts[i], &act, 1, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_EQUAL( exp, act );
  }
}
