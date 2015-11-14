/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 12-Nov-02 at 18:05:56
//  LAST-MOD:  7-May-03 at 13:28:47 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqMeshEntityTest.cpp

Unit testing of various functions in the MsqMeshEntity class.

\author Michael Brewer
\author Thomas Leurent

 */
// DESCRIP-END.
//


#include "Mesquite_MsqMeshEntity.hpp"
#include "Mesquite_Vector3D.hpp"
#include "Mesquite_PatchData.hpp"
#include "PatchDataInstances.hpp"
#include <math.h>
#include <iostream>
#include <sstream>
#include "UnitUtil.hpp"

using namespace Mesquite;
using std::cout;
using std::endl;

class MsqMeshEntityTest : public CppUnit::TestFixture
{

private:
  CPPUNIT_TEST_SUITE(MsqMeshEntityTest);
  CPPUNIT_TEST (test_hex_vertices);
  CPPUNIT_TEST (test_centroid_tri);
  CPPUNIT_TEST (test_centroid_quad);
  CPPUNIT_TEST (test_centroid_hex);
  CPPUNIT_TEST (test_unsigned_area);
  CPPUNIT_TEST (test_unsigned_area_poly);
  CPPUNIT_TEST (test_unsigned_area_tet);
  CPPUNIT_TEST (test_unsigned_area_pyr);
  CPPUNIT_TEST (test_unsigned_area_pri);
  CPPUNIT_TEST (test_unsigned_area_hex);
  CPPUNIT_TEST (test_all_nodes);
  CPPUNIT_TEST (test_check_element_orientation_linear);
  CPPUNIT_TEST (test_check_element_orientation_quadratic);
  CPPUNIT_TEST_SUITE_END();

  void test_all_nodes( EntityTopology type, unsigned num_nodes );

private:
  PatchData oneHexPatch;
  PatchData oneTetPatch;
  PatchData oneQuadPatch;
  PatchData oneTriPatch;
  Vector3D e1, e2, e3;
  double tolEps;
public:
  void setUp()
  {
    tolEps=1.e-12;
    
    // sets up the unit vectors
    e1.set(1,0,0);
    e2.set(0,1,0);
    e3.set(0,0,1);

    MsqPrintError err(cout);

    // creates empty Patch
    create_one_hex_patch(oneHexPatch, err); CPPUNIT_ASSERT(!err);
    create_one_tet_patch(oneTetPatch, err); CPPUNIT_ASSERT(!err);
    create_one_tri_patch(oneTriPatch, err); CPPUNIT_ASSERT(!err);
    create_one_quad_patch(oneQuadPatch, err); CPPUNIT_ASSERT(!err);
  }

  void tearDown()
  {
    destroy_patch_with_domain( oneTriPatch );
    destroy_patch_with_domain( oneQuadPatch );
  }
  
public:
  MsqMeshEntityTest()
  {  }
  
  void test_hex_vertices()
  {
    MsqPrintError err(cout);
    // prints out the vertices.
    const MsqVertex* ideal_vertices = oneHexPatch.get_vertex_array(err); CPPUNIT_ASSERT(!err);
    size_t num_vtx = oneHexPatch.num_nodes();
    CPPUNIT_ASSERT_EQUAL(size_t(8), num_vtx);
    
    MsqVertex vtx;

    vtx.set(1,1,1);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[0]);
    
    vtx.set(2,2,2);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[6]);
    
    vtx.set(1,2,2);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[7]);
  }

  //! test the centroid of the first element in the Patch
  void test_centroid(PatchData &pd, Vector3D &correct)
  {
    MsqPrintError err(cout);
    double eps = 1e-6;
    Vector3D centroid;

    MsqMeshEntity* elem = pd.get_element_array(err); CPPUNIT_ASSERT(!err);
    elem->get_centroid(centroid, pd, err); CPPUNIT_ASSERT(!err);

//     cout << "centroid: "<< centroid <<endl; 
//     cout << "correct: "<< correct <<endl; 

    for (int i=0; i<3; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL( centroid[i], correct[i], eps);
  }

  void test_centroid_tri()
  {
    Vector3D correct(1.5, 1+1/(2.0*sqrt(3.0)), 1.0);    
    test_centroid(oneTriPatch, correct);
  }

  void test_centroid_quad()
  {
    Vector3D correct(1.5, 1.5, 1.0);    
    test_centroid(oneQuadPatch, correct);
  }

  void test_centroid_hex()
  {
    Vector3D correct(1.5, 1.5, 1.5);    
    test_centroid(oneHexPatch, correct);
  }


  void test_unsigned_area()
     {
       MsqPrintError err(cout);
       MsqMeshEntity* tri = oneTriPatch.get_element_array(err);
       CPPUNIT_ASSERT(!err);
       CPPUNIT_ASSERT(fabs(tri->compute_unsigned_area(oneTriPatch,err)
                           -(sqrt(3.0)/4.0)) < tolEps);
       MsqMeshEntity* quad = oneQuadPatch.get_element_array(err);
       CPPUNIT_ASSERT(!err);
       CPPUNIT_ASSERT(fabs(quad->compute_unsigned_area(oneQuadPatch,err)
                           -1.0) < tolEps);
     }

  void test_unsigned_area_poly();
  void test_unsigned_area_tet();
  void test_unsigned_area_pyr();
  void test_unsigned_area_pri();
  void test_unsigned_area_hex();
  void test_all_nodes();
  
  void test_unsigned_area_common( EntityTopology type,
                                  const double* coords,
                                  double expected );
                                  
  void test_check_element_orientation_linear();
  void test_check_element_orientation_quadratic();
  void test_check_element_orientation( EntityTopology type, int nodes );
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqMeshEntityTest, "MsqMeshEntityTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqMeshEntityTest, "Unit");
 
const size_t conn[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
const bool fixed[] = { false, false, false, false, false, false, false, false };

void MsqMeshEntityTest::test_unsigned_area_poly()
{
  const double coords[] = { 0, 0, 0,
                            1, 0, 0,
                            1, 1, 0,
                          0.5, 1.5, 0,
                            0, 1, 0 };
  size_t n_vtx = 5;
  EntityTopology type = POLYGON;
  MsqError err;
  
  PatchData pd;
  pd.fill( n_vtx, coords, 1, &type, &n_vtx, conn, fixed, err );
  ASSERT_NO_ERROR(err);
  
  double a = pd.element_by_index(0).compute_unsigned_area( pd, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.25, a, 1e-8 );
}

void MsqMeshEntityTest::test_unsigned_area_common( EntityTopology type,
                                                   const double* coords,
                                                   double expected )
{
  MsqError err;
  PatchData pd;

  pd.fill( TopologyInfo::corners( type ), coords, 1, type, conn, fixed, err );
  ASSERT_NO_ERROR(err);
  
  double a = pd.element_by_index(0).compute_unsigned_area( pd, err );
  ASSERT_NO_ERROR(err);
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, a, 1e-8 );
}


void MsqMeshEntityTest::test_unsigned_area_tet()
{
  const double coords[] = { 0, 0, 0,
                            1, 0, 0,
                            0, 1, 0,
                            0, 0, 1 };
  test_unsigned_area_common( TETRAHEDRON, coords, 1.0/6.0 );
}

void MsqMeshEntityTest::test_unsigned_area_pyr()
{
  const double coords[] = { 0, 0, 0,
                            1, 0, 0,
                            1, 1, 0,
                            0, 1, 0,
                            0, 0, 1 };
  test_unsigned_area_common( PYRAMID, coords, 1.0/3.0 );

  const double pyr_coords[] = {-1, -1, -1,
                                1, -1, -1,
                                1,  1, -1,
                               -1,  1, -1,
                                0,  0,  0 };
  test_unsigned_area_common( PYRAMID, pyr_coords, 4.0/3.0 );
}

void MsqMeshEntityTest::test_unsigned_area_pri()
{
  const double coords[] = { 0, 0, 0,
                            1, 0, 0,
                            0, 1, 0,
                            0, 0, 1,
                            1, 0, 1,
                            0, 1, 1 };
  test_unsigned_area_common( PRISM, coords, 0.5 );
  
  const double tet_coords[] = { 0, 0, 0,
                                1, 0, 0,
                                0, 1, 0,
                                0, 0, 1,
                                0, 0, 1,
                                0, 0, 1 };
  test_unsigned_area_common( PRISM, tet_coords, 1.0/6.0 );
}

void MsqMeshEntityTest::test_unsigned_area_hex()
{
  const double coords[] = { 0, 0, 0,
                            1, 0, 0,
                            1, 1, 0,
                            0, 1, 0,
                            0, 0, 1,
                            1, 0, 1,
                            1, 1, 1,
                            0, 1, 1 };
  test_unsigned_area_common( HEXAHEDRON, coords, 1.0 );

  const double coords2[] = { 0, 0, 0,
                             2, 0, 0,
                             2, 2, 0,
                             0, 2, 0,
                             0, 0, 2,
                             2, 0, 2,
                             2, 2, 2,
                             0, 2, 2 };
  test_unsigned_area_common( HEXAHEDRON, coords2, 8.0 );
  
  const double pyr_coords[] = {-1,-1, 0,
                                1,-1, 0,
                                1, 1, 0,
                               -1, 1, 0,
                                0, 0, 1,
                                0, 0, 1,
                                0, 0, 1,
                                0, 0, 1 };
  test_unsigned_area_common( HEXAHEDRON, pyr_coords, 4.0/3.0 );
}

void MsqMeshEntityTest::test_all_nodes( EntityTopology type, unsigned num_nodes )
{
  const unsigned num_vtx = 27;
  double coords[3*num_vtx] = {0.0};
  size_t conn[num_vtx];
  for (size_t i = 0; i < num_vtx; ++i)
    conn[i] = i;
  bool fixed[num_vtx] = {false};
  CPPUNIT_ASSERT(num_nodes <= num_vtx);
  
  MsqError err;
  PatchData pd;
  size_t n = num_nodes;
  pd.fill( num_nodes, coords, 1, &type, &n, conn, fixed, err );
  ASSERT_NO_ERROR(err);
  
  MsqMeshEntity& elem = pd.element_by_index(0);
  NodeSet all = elem.all_nodes(err);
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( num_nodes, all.num_nodes() );
  CPPUNIT_ASSERT( all.have_any_corner_node() );
  bool mid_edge, mid_face, mid_reg;
  TopologyInfo::higher_order( type, num_nodes, mid_edge, mid_face, mid_reg, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( mid_edge, !!all.have_any_mid_edge_node() );
  CPPUNIT_ASSERT_EQUAL( mid_face, !!all.have_any_mid_face_node() );
  CPPUNIT_ASSERT_EQUAL( mid_reg,  !!all.have_any_mid_region_node() );
  
}

void MsqMeshEntityTest::test_all_nodes() 
{
  test_all_nodes( TRIANGLE, 3 );
  test_all_nodes( TRIANGLE, 4 );
  test_all_nodes( TRIANGLE, 6 );
  test_all_nodes( TRIANGLE, 7 );

  test_all_nodes( QUADRILATERAL, 4 );
  test_all_nodes( QUADRILATERAL, 5 );
  test_all_nodes( QUADRILATERAL, 8 );
  test_all_nodes( QUADRILATERAL, 9 );
  
  test_all_nodes( TETRAHEDRON, 4 );
  test_all_nodes( TETRAHEDRON, 5 );
  test_all_nodes( TETRAHEDRON, 10 );
  test_all_nodes( TETRAHEDRON, 11 );
  test_all_nodes( TETRAHEDRON, 8 );
  test_all_nodes( TETRAHEDRON, 9 );
  test_all_nodes( TETRAHEDRON, 14 );
  test_all_nodes( TETRAHEDRON, 15 );

  test_all_nodes( HEXAHEDRON, 8 );
  test_all_nodes( HEXAHEDRON, 9 );
  test_all_nodes( HEXAHEDRON, 20 );
  test_all_nodes( HEXAHEDRON, 21 );
  test_all_nodes( HEXAHEDRON, 14 );
  test_all_nodes( HEXAHEDRON, 15 );
  test_all_nodes( HEXAHEDRON, 26 );
  test_all_nodes( HEXAHEDRON, 27 );
}

void MsqMeshEntityTest::test_check_element_orientation_linear()
{
  const EntityTopology types[] = { TRIANGLE,
                                   QUADRILATERAL,
                                   TETRAHEDRON,
                                   PYRAMID,
                                   PRISM,
                                   HEXAHEDRON };
  const int num_types = sizeof(types)/sizeof(types[0]);

  for (int i = 0; i < num_types; ++i) {
    test_check_element_orientation( types[i], TopologyInfo::corners(types[i]) );
  }
}

void MsqMeshEntityTest::test_check_element_orientation_quadratic()
{
  struct ElemType {
    EntityTopology topo;
    unsigned nodes;
  };

  const ElemType types[] = {
    { TRIANGLE, 6 },
    { QUADRILATERAL, 8 },
    { QUADRILATERAL, 9 },
    { TETRAHEDRON, 10 } };
  const int num_types = sizeof(types)/sizeof(types[0]);
  
  for (int i = 0; i < num_types; ++i) {
    test_check_element_orientation( types[i].topo, types[i].nodes );
  }
}

void MsqMeshEntityTest::test_check_element_orientation( EntityTopology type, 
                                                        int nodes )
{
    // get an ideal element
  MsqError err;
  PatchData pd;
  create_ideal_element_patch( pd, type, nodes, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( (size_t)1, pd.num_elements() );
  CPPUNIT_ASSERT_EQUAL( (size_t)nodes, pd.num_nodes() );
  MsqMeshEntity& elem = pd.element_by_index(0);
  CPPUNIT_ASSERT_EQUAL( (size_t)nodes, elem.node_count() );
  CPPUNIT_ASSERT_EQUAL( type, elem.get_element_type() );
  const size_t* conn = elem.get_vertex_index_array();
  
    // test that ideal element is not reported as inverted
  int inverted, tested;
  elem.check_element_orientation( pd, inverted, tested, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( 0, inverted );
  CPPUNIT_ASSERT( tested > 0 );

  bool mids[4] = {false};
  TopologyInfo::higher_order( type, nodes, mids[1], mids[2], mids[3], err );
  MSQ_ERRRTN(err);
  
    // invert element at each vertex and test
  Vector3D centroid;
  elem.get_centroid( centroid, pd, err );
  ASSERT_NO_ERROR(err);
  for (int i = 0; i < nodes; ++i) {
    unsigned dim, num;
    TopologyInfo::side_from_higher_order( type, nodes, i, dim, num, err );
    ASSERT_NO_ERROR(err);
    const Vector3D old_pos = pd.vertex_by_index( conn[i] );
    Vector3D new_pos = old_pos;
    if (dim == TopologyInfo::dimension(type)) { 
      // move mid-element node 3/4 of the way to corner 0
      new_pos += 3*pd.vertex_by_index( conn[0] );
      new_pos *= 0.25;
    }
    else if (dim == 0) { // if a corner vertex
      if (type == TRIANGLE || type == TETRAHEDRON) {
        // move tri/tet vertex past opposite side of element
        new_pos += 2*(centroid - old_pos);
      }
      else if (mids[1]) {
        // if have mid-edge nodes move 3/4 of the way to center vertex
        new_pos += 3*centroid;
        new_pos *= 0.25;
      }
      else {
        // move vertex past centroid
        new_pos += 1.5*(centroid - old_pos);
      }
    }
    else {
      // otherwise move vertex past centroid
      new_pos += 2.5*(centroid - old_pos);
    }
    
    pd.set_vertex_coordinates( new_pos, conn[i], err );
    ASSERT_NO_ERROR(err);
    
      // test that element is inverted
    inverted = tested = 0;
    elem.check_element_orientation( pd, inverted, tested, err );
    ASSERT_NO_ERROR(err);
    std::ostringstream str;
    str << TopologyInfo::short_name(type) << nodes
        << " Vertex " << i 
        << " (Dimension " << dim
        << " Index " << num << ")";
    CppUnit::Message m( "MsqMeshEntity failed to detect inverted element" );
    m.addDetail( str.str() );
    ASSERT_MESSAGE( m, inverted > 0 );
    
      // move vertex back to ideal position
    pd.set_vertex_coordinates( old_pos, conn[i], err );
    ASSERT_NO_ERROR(err);
  }
}
