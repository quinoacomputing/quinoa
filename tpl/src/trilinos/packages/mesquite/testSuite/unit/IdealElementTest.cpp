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


/** \file IdealElementTest.cpp
 *  \brief Test the ExtraData functionality of the PatchData class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "IdealElements.hpp"
#include "TopologyInfo.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "cppunit/extensions/HelperMacros.h"

using namespace Mesquite;

class IdealElementTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(IdealElementTest);
  CPPUNIT_TEST (test_unit_tri);
  CPPUNIT_TEST (test_unit_quad);
  CPPUNIT_TEST (test_unit_tet);
  CPPUNIT_TEST (test_unit_pyr);
  CPPUNIT_TEST (test_unit_wdg);
  CPPUNIT_TEST (test_unit_hex);
  CPPUNIT_TEST (test_side_height_pyr);
  CPPUNIT_TEST (test_unit_edge_tri);
  CPPUNIT_TEST (test_unit_edge_quad);
  CPPUNIT_TEST (test_unit_edge_tet);
  CPPUNIT_TEST (test_unit_edge_pyr);
  CPPUNIT_TEST (test_unit_edge_wdg);
  CPPUNIT_TEST (test_unit_edge_hex);
  CPPUNIT_TEST (test_unit_height_pyr);
  CPPUNIT_TEST_SUITE_END();

public:
  void test_unit_tri();
  void test_unit_quad();
  void test_unit_tet();
  void test_unit_pyr();
  void test_unit_wdg();
  void test_unit_hex();
  void test_unit_edge_tri();
  void test_unit_edge_quad();
  void test_unit_edge_tet();
  void test_unit_edge_pyr();
  void test_unit_edge_wdg();
  void test_unit_edge_hex();
  void test_side_height_pyr();
  void test_unit_height_pyr();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(IdealElementTest, "IdealElementTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(IdealElementTest, "Unit");

static double distance_from_origin( unsigned n, const Vector3D* vtx );
static bool unit_edge_lengths( EntityTopology type, const Vector3D* coords );
static bool equal_edge_lengths( EntityTopology type, const Vector3D* coords );

void IdealElementTest::test_unit_tri()
{
  const Vector3D* coords = unit_element( TRIANGLE );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(3, coords), 1e-6);
  CPPUNIT_ASSERT( equal_edge_lengths( TRIANGLE, coords ) );
  
  Vector3D v1 = coords[1] - coords[0];
  Vector3D v2 = coords[2] - coords[0];
  double area = 0.5 * (v1 * v2).length();
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, area, 1e-6 );
}

void IdealElementTest::test_unit_quad()
{
  const Vector3D* coords = unit_element( QUADRILATERAL );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(4, coords), 1e-6);
  CPPUNIT_ASSERT( equal_edge_lengths( QUADRILATERAL, coords ) );
  
  Vector3D p1 = 0.5 * (coords[0] + coords[1]);
  Vector3D p2 = 0.5 * (coords[1] + coords[2]);
  Vector3D p3 = 0.5 * (coords[2] + coords[3]);
  Vector3D p4 = 0.5 * (coords[3] + coords[0]);
  Vector3D v1 = p1 - p3;
  Vector3D v2 = p2 - p4;
  double area = (v1 * v2).length();
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, area, 1e-6 );
}

void IdealElementTest::test_unit_tet()
{
  const Vector3D* coords = unit_element( TETRAHEDRON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(4, coords), 1e-6);
  CPPUNIT_ASSERT( equal_edge_lengths( TETRAHEDRON, coords ) );
  
  Vector3D v1 = coords[1] - coords[0];
  Vector3D v2 = coords[2] - coords[0];
  Vector3D v3 = coords[3] - coords[0];
  double vol = (v3 % (v1 * v2 ))/6.0;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, vol, 1e-6 );
}

void IdealElementTest::test_unit_pyr()
{
  const Vector3D* coords = unit_element( PYRAMID );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(5, coords), 1e-6);
  CPPUNIT_ASSERT( equal_edge_lengths( PYRAMID, coords ) );
  
  Vector3D p1 = 0.5 * (coords[0] + coords[1]);
  Vector3D p2 = 0.5 * (coords[1] + coords[2]);
  Vector3D p3 = 0.5 * (coords[2] + coords[3]);
  Vector3D p4 = 0.5 * (coords[3] + coords[0]);
  Vector3D v1 = p1 - p3;
  Vector3D v2 = p2 - p4;
  double area = (v1 * v2).length();
  
  Vector3D n = (v1 * v2);
  n.normalize();
  Vector3D c = 0.25*(coords[0]+coords[1]+coords[2]+coords[3]);
  double height = n % (coords[4] - c);
  double vol = area * height / 3.0;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, vol, 1e-6 );
}

void IdealElementTest::test_unit_wdg()
{
  const Vector3D* coords = unit_element( PRISM );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(6, coords), 1e-6);
  CPPUNIT_ASSERT( equal_edge_lengths( PRISM, coords ) );
  
  Vector3D v1 = coords[1] - coords[0];
  Vector3D v2 = coords[2] - coords[0];
  Vector3D v3 = coords[3] - coords[0];
  double vol = 0.5 * v3 % (v1 * v2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, vol, 1e-6 );
}

void IdealElementTest::test_unit_hex()
{
  const Vector3D* coords = unit_element( HEXAHEDRON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(8, coords), 1e-6);
  CPPUNIT_ASSERT( equal_edge_lengths( HEXAHEDRON, coords ) );
  
  Vector3D p1 = 0.25*(coords[0]+coords[1]+coords[2]+coords[3]);
  Vector3D p2 = 0.25*(coords[4]+coords[5]+coords[6]+coords[7]);
  Vector3D p3 = 0.25*(coords[0]+coords[1]+coords[5]+coords[4]);
  Vector3D p4 = 0.25*(coords[7]+coords[6]+coords[6]+coords[7]);
  Vector3D p5 = 0.25*(coords[0]+coords[3]+coords[4]+coords[7]);
  Vector3D p6 = 0.25*(coords[1]+coords[2]+coords[6]+coords[5]);
  Vector3D v1 = p1 - p2;
  Vector3D v2 = p3 - p4;
  Vector3D v3 = p5 - p6;
  double vol = v3 % (v1 * v2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, vol, 1e-6 );
}

void IdealElementTest::test_side_height_pyr()
{
  const Vector3D* coords = unit_element( PYRAMID, true );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(5, coords), 1e-6);
  CPPUNIT_ASSERT( equal_edge_lengths( QUADRILATERAL, coords ) );
  
  Vector3D p1 = 0.5 * (coords[0] + coords[1]);
  Vector3D p2 = 0.5 * (coords[1] + coords[2]);
  Vector3D p3 = 0.5 * (coords[2] + coords[3]);
  Vector3D p4 = 0.5 * (coords[3] + coords[0]);
  Vector3D v1 = p1 - p3;
  Vector3D v2 = p2 - p4;
  double area = (v1 * v2).length();
  
  Vector3D n = (v1 * v2);
  n.normalize();
  Vector3D c = 0.25*(coords[0]+coords[1]+coords[2]+coords[3]);
  double height = n % (coords[4] - c);
  double vol = area * height / 3.0;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, vol, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( height, (coords[0] - coords[1]).length(), 1e-6 );
}

void IdealElementTest::test_unit_edge_tri()
{
  const Vector3D* coords = unit_edge_element( TRIANGLE );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(3, coords), 1e-6);
  CPPUNIT_ASSERT( unit_edge_lengths( TRIANGLE, coords ) );
}

void IdealElementTest::test_unit_edge_quad()
{
  const Vector3D* coords = unit_edge_element( QUADRILATERAL );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(4, coords), 1e-6);
  CPPUNIT_ASSERT( unit_edge_lengths( QUADRILATERAL, coords ) );
}

void IdealElementTest::test_unit_edge_tet()
{
  const Vector3D* coords = unit_edge_element( TETRAHEDRON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(4, coords), 1e-6);
  CPPUNIT_ASSERT( unit_edge_lengths( TETRAHEDRON, coords ) );
}

void IdealElementTest::test_unit_edge_pyr()
{
  const Vector3D* coords = unit_edge_element( PYRAMID );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(5, coords), 1e-6);
  CPPUNIT_ASSERT( unit_edge_lengths( PYRAMID, coords ) );
}

void IdealElementTest::test_unit_edge_wdg()
{
  const Vector3D* coords = unit_edge_element( PRISM );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(6, coords), 1e-6);
  CPPUNIT_ASSERT( unit_edge_lengths( PRISM, coords ) );
}

void IdealElementTest::test_unit_edge_hex()
{
  const Vector3D* coords = unit_edge_element( HEXAHEDRON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(8, coords), 1e-6);
  CPPUNIT_ASSERT( unit_edge_lengths( HEXAHEDRON, coords ) );
}

void IdealElementTest::test_unit_height_pyr()
{
  const Vector3D* coords = unit_edge_element( PYRAMID, true );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, distance_from_origin(5, coords), 1e-6);
  CPPUNIT_ASSERT( unit_edge_lengths( QUADRILATERAL, coords ) );
  
  Vector3D p1 = 0.5 * (coords[0] + coords[1]);
  Vector3D p2 = 0.5 * (coords[1] + coords[2]);
  Vector3D p3 = 0.5 * (coords[2] + coords[3]);
  Vector3D p4 = 0.5 * (coords[3] + coords[0]);
  Vector3D v1 = p1 - p3;
  Vector3D v2 = p2 - p4;
 
  Vector3D n = (v1 * v2);
  n.normalize();
  Vector3D c = 0.25*(coords[0]+coords[1]+coords[2]+coords[3]);
  double height = n % (coords[4] - c);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( height, (coords[0] - coords[1]).length(), 1e-6 );
}

static void get_edge_lengths( EntityTopology type, 
                              const Vector3D* coords,
                              double& min, double& max )
{
  MsqError err;
  unsigned n = TopologyInfo::edges( type );
  min = HUGE_VAL;
  max = -HUGE_VAL;
  for (unsigned j = 0; j < n; ++j) {
    const unsigned* edge = TopologyInfo::edge_vertices( type, j, err );
    CPPUNIT_ASSERT(!err);
    double len = (coords[edge[0]] - coords[edge[1]]).length();
    if (min > len) min = len;
    if (max < len) max = len;
  }
}

static bool unit_edge_lengths( EntityTopology type, const Vector3D* coords )
{
  double min, max;
  get_edge_lengths( type, coords, min, max );
  return fabs(1.0 - min) < 1e-6 && fabs(1.0 - max) < 1e-6;
}

static bool equal_edge_lengths( EntityTopology type, const Vector3D* coords )
{
  double min, max;
  get_edge_lengths( type, coords, min, max );
  return (max - min) < 1e-6;
}

static double distance_from_origin( unsigned n, const Vector3D* vtx )
{
  Vector3D sum(0,0,0);
  for (unsigned i = 0; i < n; ++i)
    sum += vtx[i];
  if (n == 5) // pyrmid centroid is not mean of corners
    sum += (vtx[4]/3.0);

  return sum.length();
}
