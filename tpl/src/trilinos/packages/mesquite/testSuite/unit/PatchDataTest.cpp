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
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 13-Nov-02 at 18:05:56
//  LAST-MOD:  9-Jun-04 at 14:50:51 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file PatchDataTest.cpp

Unit testing of various functions in the PatchData class. 

 */
// DESCRIP-END.
//



#include "Mesquite.hpp"
#include "Mesquite_PatchData.hpp"
#include "PatchDataInstances.hpp"
#include "UnitUtil.hpp"
#include "Mesquite_Settings.hpp"
#include "Mesquite_Instruction.hpp"

#include "Mesquite_ArrayMesh.hpp"
#include "Mesquite_DomainClassifier.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_MeshDomain1D.hpp"
#include "Mesquite_MeshDecorator.hpp"

//#include "Mesquite_TriLagrangeShape.hpp"

#include "Mesquite_cppunit/extensions/HelperMacros.h"

#include <algorithm>
#include <set>
#include <map>

using namespace Mesquite;

using std::cout;
using std::cerr;
using std::endl;

class PatchDataTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(PatchDataTest);
  CPPUNIT_TEST (test_num_corners);
  CPPUNIT_TEST (test_get_element_vertex_indices);
  CPPUNIT_TEST (test_get_vertex_element_indices);
  CPPUNIT_TEST (test_get_element_vertex_coordinates);
  CPPUNIT_TEST (test_move_free_vertices_constrained);
  CPPUNIT_TEST (test_movement_function);
  CPPUNIT_TEST (test_get_adj_elems_2d);
  CPPUNIT_TEST (test_get_minmax_element_area);
  CPPUNIT_TEST (test_sub_patch);
  CPPUNIT_TEST (test_fill);
  CPPUNIT_TEST (test_reorder);
  CPPUNIT_TEST (test_update_slave_node_coords);
  CPPUNIT_TEST (test_patch_data_fill_slaved_ho_nodes);
  CPPUNIT_TEST (test_patch_reorder_ho_nodes);
  CPPUNIT_TEST (test_patch_data_fill_free_ho_nodes);
  CPPUNIT_TEST (test_patch_data_mesh_slaved_ho_nodes);
  CPPUNIT_TEST (test_patch_data_mesh_free_ho_nodes);
  CPPUNIT_TEST (test_patch_data_mesh_calcualted_ho_nodes);
  CPPUNIT_TEST (test_patch_data_mesh_flagged_ho_nodes);
  CPPUNIT_TEST (test_vertex_verts_fixed);
  CPPUNIT_TEST (test_curve_verts_fixed);
  CPPUNIT_TEST (test_surf_verts_fixed);
  CPPUNIT_TEST_SUITE_END();
   
private:

   MsqVertex vtx_0_0;
   MsqVertex vtx_0_1;
   MsqVertex vtx_1_0;
   MsqVertex vtx_1_1;
   MsqVertex vtx_2_0;
   MsqVertex vtx_2_1;

   MsqMeshEntity tri1;
   MsqMeshEntity tri2;
   MsqMeshEntity quad1;   
   
   PatchData mPatch2D;
   
   void test_quad8_patch( bool reorder, bool ho_nodes_slaved );
   void get_quad8_mesh( Mesh*& mesh_out );
   void get_quad8_mesh_and_domain( Mesh*& mesh_out, MeshDomain*& domain_out );
   void get_higher_order_vertices( Mesh* mesh, 
                                   std::map<Mesh::VertexHandle,bool>& ho_verts,
                                   bool initial_value = false,
                                   bool non_fixed_only = true );
   void check_higher_order_vertices_slaved( Mesh* mesh,
                   Settings::HigherOrderSlaveMode mode,
                   const std::map<Mesh::VertexHandle,bool>& expected );
   

public:
  void setUp()
  {
     MsqPrintError err(cout);

     /* our 2D set up: 2 triangles and one quad are available
       1___3___5
        |\1|   |
        |0\| 2 |
       0---2---4
     */
     vtx_0_0.set(0,0,0);
     vtx_0_1.set(0,1,0);
     vtx_1_0.set(1,0,0);
     vtx_1_1.set(1,1,0);
     vtx_2_0.set(2,0,0);
     vtx_2_1.set(2,1,0);
     
     double coords[] = { 0,0,0,
                         0,1,0,
                         1,0,0,
                         1,1,0,
                         2,0,0,
                         2,1,0 };
     
     EntityTopology types[] = { TRIANGLE, TRIANGLE, QUADRILATERAL };
     
     size_t connectivity[] = { 0, 2, 1,
                               1, 2, 3,
                               3, 2, 4, 5 };
     
     size_t counts[] = { 3, 3, 4 };
     
     mPatch2D.fill( 6, coords,
                    3, types,
                    counts, connectivity,
                    0, err );
  }
  
  void tearDown()
  {
  }
  
public:
  PatchDataTest()
    {}
  
   void test_num_corners()
   {
     MsqPrintError err(cout);
     size_t n = mPatch2D.num_corners();
     CPPUNIT_ASSERT(n==10);
   }

   void test_get_element_vertex_indices()
   {

      MsqPrintError err(cout);
      
      std::vector<size_t> vtx_ind;
      std::vector<size_t> res;

      // test we get the right vertices for element 1 (tri)
      mPatch2D.get_element_vertex_indices(1, vtx_ind, err); CPPUNIT_ASSERT(!err);
      res.push_back(1); res.push_back(2); res.push_back(3);
      CPPUNIT_ASSERT( vtx_ind==res );

      // test we get the right vertices for element 2 (quad)
      vtx_ind.clear(); res.clear();
      mPatch2D.get_element_vertex_indices(2, vtx_ind, err); CPPUNIT_ASSERT(!err);
      res.push_back(3); res.push_back(2); res.push_back(4); res.push_back(5);
      CPPUNIT_ASSERT( vtx_ind==res );
   }

   void test_get_vertex_element_indices()
   {
     /*  1___3___5
         |\1|   |
         |0\| 2 |
         0---2---4   */
     MsqPrintError err(cout);
     
     std::vector<size_t> elem_ind;
     std::vector<size_t> res;
     
     mPatch2D.generate_vertex_to_element_data();
     
     // test we get the elements contiguous to vertex 3
     mPatch2D.get_vertex_element_indices(3, elem_ind,err); CPPUNIT_ASSERT(!err);
     res.push_back(2); res.push_back(1);
     CPPUNIT_ASSERT(res==elem_ind);
     
     // test we get the elements contiguous to vertex 2
     elem_ind.clear(); res.clear();
     mPatch2D.get_vertex_element_indices(2, elem_ind,err); CPPUNIT_ASSERT(!err);
     res.push_back(2); res.push_back(1); res.push_back(0);
     CPPUNIT_ASSERT(res==elem_ind);
   }

   void test_get_element_vertex_coordinates()
   {
      MsqPrintError err(cout);

      std::vector< Vector3D > coords;
      mPatch2D.get_element_vertex_coordinates(1, coords,err); CPPUNIT_ASSERT(!err);
      
      CPPUNIT_ASSERT( coords[0]==vtx_0_1 );
      CPPUNIT_ASSERT( coords[1]==vtx_1_0 );
      CPPUNIT_ASSERT( coords[2]==vtx_1_1 );
   }

   /* This tests the move_vertices() function as well as the
      PatchDataCoordsMemento functionality
      */
   void test_move_free_vertices_constrained()
   {
      MsqPrintError err(cout);

      // gets a memento of the patch coordinates.
      PatchDataVerticesMemento* coords_mem = mPatch2D.create_vertices_memento(err);
      CPPUNIT_ASSERT(!err);
      
      // Move the two first vertices in direction dk by step size s;
      Vector3D dk[6];
      dk[0].set(-1,-2,0);
      dk[1].set(-1, 2,0);
      double s = 0.3;
      mPatch2D.move_free_vertices_constrained(dk, 6, s, err); CPPUNIT_ASSERT(!err);

      // gets the new coordinates and  checks the vertices were displaced as expected.
      std::vector< Vector3D > coords;
      mPatch2D.get_element_vertex_coordinates(0, coords,err);
      Vector3D new_vtx_0_0 = vtx_0_0 + s*dk[0];
      Vector3D new_vtx_0_1 = vtx_0_1 + s*dk[1];
      CPPUNIT_ASSERT(coords[0] == new_vtx_0_0);
      CPPUNIT_ASSERT(coords[2] == new_vtx_0_1);

      // restore the PatchData to previous coords.
      mPatch2D.set_to_vertices_memento(coords_mem, err); CPPUNIT_ASSERT(!err);

      // gets the new coordinates and  checks the vertices are back to original.
      coords.clear();
      mPatch2D.get_element_vertex_coordinates(0, coords,err);
      CPPUNIT_ASSERT(coords[0] == vtx_0_0);
      CPPUNIT_ASSERT(coords[2] == vtx_0_1);

      delete coords_mem;
   }

   void test_movement_function()
   {
      MsqPrintError err(cout);
      // gets a memento of the patch coordinates.
      PatchDataVerticesMemento* coords_mem = mPatch2D.create_vertices_memento(err);
      CPPUNIT_ASSERT(!err);
      
      // Move the two first vertices in direction dk by step size s;
      Vector3D dk[6];
      dk[0].set(0,-2,0);
      dk[1].set(-1,0,0);
      double s = 1;
      mPatch2D.move_free_vertices_constrained(dk, 6, 1, err); CPPUNIT_ASSERT(!err);
      // gets the new coordinates and  checks the vertices were displaced as expected.
      std::vector< Vector3D > coords;
      mPatch2D.get_element_vertex_coordinates(0, coords,err);
      Vector3D new_vtx_0_0 = vtx_0_0 + s*dk[0];
      Vector3D new_vtx_0_1 = vtx_0_1 + s*dk[1];
      CPPUNIT_ASSERT(coords[0] == new_vtx_0_0);
      CPPUNIT_ASSERT(coords[2] == new_vtx_0_1);
      double m_dist=mPatch2D.get_max_vertex_movement_squared(coords_mem,err);
      CPPUNIT_ASSERT(m_dist==4.0);
      // restore the PatchData to previous coords.
      mPatch2D.set_to_vertices_memento(coords_mem, err); CPPUNIT_ASSERT(!err);
      // gets the new coordinates and  checks the vertices are back to original.
      coords.clear();
      mPatch2D.get_element_vertex_coordinates(0, coords,err);
      CPPUNIT_ASSERT(coords[0] == vtx_0_0);
      CPPUNIT_ASSERT(coords[2] == vtx_0_1);

      delete coords_mem;
   }
  
/*Tests the function PatchData::get_adjacent_entities_via_n_dim()
  which finds the elements adjacent to a given element.  If 'n'
  equals 0 the elements must share a vertex; if 'n' equals 1 the
  elements must share an edge; and if 'n' equals 2 the elements
  must share a face.*/
   void test_get_adj_elems_2d()
   {
     MsqPrintError err(cout);
     std::vector<size_t> elems_0;
       //find elements sharing an edge with oth elem (should be 1)
     mPatch2D.get_adjacent_entities_via_n_dim(1, 0, elems_0, err);
     CPPUNIT_ASSERT(!err);
     CPPUNIT_ASSERT(elems_0.back() == 1);
     std::vector<size_t> elems_1;
       //find elements sharing an edge with 1st elem (should be 0 and 2)
     mPatch2D.get_adjacent_entities_via_n_dim(1, 1, elems_1, err);
     CPPUNIT_ASSERT(!err);
     CPPUNIT_ASSERT(elems_1.size() == 2);
     std::vector<size_t> elems_2;
       //find elements sharing an vert with 0th elem (should be 1 and 2).
     mPatch2D.get_adjacent_entities_via_n_dim(0, 0, elems_2, err);
     CPPUNIT_ASSERT(!err);
     CPPUNIT_ASSERT(elems_2.size() == 2);
     std::vector<size_t> elems_3;
     //find elements sharing an face with 0th elem (should be empty).
     mPatch2D.get_adjacent_entities_via_n_dim(2, 0, elems_3, err);
     CPPUNIT_ASSERT(!err);
     CPPUNIT_ASSERT(elems_3.size() == 0);
   }

  
   void test_get_minmax_element_area()
   {
     MsqPrintError err(cout);
     double min, max;
     mPatch2D.get_minmax_element_unsigned_area(min, max, err); CPPUNIT_ASSERT(!err);

     CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, min, 0.0001 );
     CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, max, 0.0001 );
   }        
     
  void check_sub_patch( unsigned vtx, unsigned layers, PatchData& pd, PatchData& sub );
  
  void test_sub_patch( );
  
  void test_patch_contents( bool reorder );
  
  void test_fill() { test_patch_contents(false); }
  void test_reorder() { test_patch_contents(true); }
  
  void test_update_slave_node_coords();
  void test_patch_data_fill_slaved_ho_nodes() { test_quad8_patch(false, true ); }
  void test_patch_reorder_ho_nodes()          { test_quad8_patch(true,  true ); }
  void test_patch_data_fill_free_ho_nodes()   { test_quad8_patch(false, false); }

  void test_patch_data_mesh_slaved_ho_nodes();
  void test_patch_data_mesh_free_ho_nodes();
  void test_patch_data_mesh_calcualted_ho_nodes();
  void test_patch_data_mesh_flagged_ho_nodes();

  void test_fixed_by_geom_dim( unsigned dim );
  void test_vertex_verts_fixed() { test_fixed_by_geom_dim(0); }
  void test_curve_verts_fixed()  { test_fixed_by_geom_dim(1); }
  void test_surf_verts_fixed()   { test_fixed_by_geom_dim(2); }
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PatchDataTest, "PatchDataTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PatchDataTest, "Unit");

void PatchDataTest::check_sub_patch( unsigned vtx, unsigned layers, PatchData& pd, PatchData& sub )
{
  unsigned i, j;
  std::set<size_t> seen;
  std::vector<size_t> vtx_map, elem_map;

    // test vertex list consistency 
  vtx_map.resize( sub.num_nodes() );
  for (i = 0; i < sub.num_nodes(); ++i) {
      // get index in old patch for this vertex
    Mesh::VertexHandle h = sub.get_vertex_handles_array()[i];
    Mesh::VertexHandle* end = pd.get_vertex_handles_array() + pd.num_nodes();
    Mesh::VertexHandle* ptr = std::find( pd.get_vertex_handles_array(), end, h );
    CPPUNIT_ASSERT(ptr != end);
    size_t idx = ptr - pd.get_vertex_handles_array();
    CPPUNIT_ASSERT( idx < pd.num_nodes() );
// put handle in map
    vtx_map[i] = idx;
      // make sure we don't have duplicates of vertices
    CPPUNIT_ASSERT( seen.insert(idx).second );
      // make sure vertices have same coords
    CPPUNIT_ASSERT_VECTORS_EQUAL( pd.vertex_by_index(idx), sub.vertex_by_index(i), 1e-12 );
  }

    // test element list consistency 
  seen.clear();
  elem_map.resize( sub.num_elements() );
  for (i = 0; i < sub.num_elements(); ++i) {
      // get index in old patch for element
    Mesh::ElementHandle h = sub.get_element_handles_array()[i];
    Mesh::ElementHandle* end = pd.get_element_handles_array() + pd.num_nodes();
    Mesh::ElementHandle* ptr = std::find( pd.get_element_handles_array(), end, h );
    CPPUNIT_ASSERT(ptr != end);
    size_t idx = ptr - pd.get_element_handles_array();
    CPPUNIT_ASSERT( idx < pd.num_elements() );
// put handle in map
    elem_map[i] = idx;
      // make sure we don't have duplicate elements
    CPPUNIT_ASSERT( seen.insert(idx).second );
      // get elements
    MsqMeshEntity& elem1 = pd.element_by_index(idx);
    MsqMeshEntity& elem2 = sub.element_by_index(i);
      // compare element data
    CPPUNIT_ASSERT_EQUAL( elem1.get_element_type(), elem2.get_element_type() );
    CPPUNIT_ASSERT_EQUAL( elem1.node_count(), elem2.node_count() );
      // get connectivity for elements
    std::vector<size_t> vtx1, vtx2;
    elem1.get_node_indices( vtx1 );
    elem2.get_node_indices( vtx2 );
    CPPUNIT_ASSERT_EQUAL( vtx1.size(), vtx2.size() );
      // compare connectivity
    for (j = 0; j < vtx1.size(); ++j) {
      CPPUNIT_ASSERT( vtx1[j] < pd.num_nodes() );
      CPPUNIT_ASSERT( vtx2[j] < sub.num_nodes() );
      CPPUNIT_ASSERT_EQUAL( vtx1[j], vtx_map[vtx2[j]] );
    }
  }

    // test that the subpatch has the elements adjacent to the specified 
    // vertex.

    // first get list of adjacent elements in original patch
  seen.clear();
  for (i = 0; i < pd.num_elements(); ++i) {
    std::vector<size_t> vtx_list;
    pd.element_by_index(i).get_node_indices( vtx_list );
    if (std::find( vtx_list.begin(), vtx_list.end(), vtx ) != vtx_list.end())
      seen.insert(i);
  }

    // if 1 layer, then should match element count
  if (layers == 1) {
    CPPUNIT_ASSERT_EQUAL( seen.size(), sub.num_elements() );
  }

    // remove from the set each element in the subpatch
  for (i = 0; i < sub.num_elements(); ++i) {
    size_t idx = elem_map[i];
    std::set<size_t>::iterator it = seen.find( idx );
    if (it != seen.end()) {
      seen.erase(it);
    }
    else {
      CPPUNIT_ASSERT( layers > 1 );
    }
  }
  CPPUNIT_ASSERT( seen.empty() );
}

void PatchDataTest::test_sub_patch( )
{
  MsqPrintError err( std::cout );
  PatchData pd, sub;
  create_twelve_hex_patch( pd, err );
  CPPUNIT_ASSERT(!err);

  for (unsigned i = 0; i < pd.num_free_vertices(); ++i) {
    unsigned layers = i % 2 ? 2 : 1;
    pd.get_subpatch( i, layers, sub, err );
    CPPUNIT_ASSERT(!err);

    check_sub_patch( i, layers, pd, sub );
  }
}



void PatchDataTest::test_patch_contents( bool reorder )
{
  const unsigned NUM_VERTEX = 15;
  const unsigned NUM_ELEMENT = 9;

    // Mesh data used to populate patch
    // Use a relatively randomized order for verices
    // so patch reordering will result in a changed
    // vertex ordering.
  double coords[3*NUM_VERTEX] = { 6, 6, 3, // 0
                                  0, 0, 0,
                                  0, 6, 3,
                                  4, 2, 2, 
                                  2, 4, 2,
                                  4, 4, 2, // 5
                                  0, 6, 3,
                                  2, 2, 1, 
                                  2, 6, 3, 
                                  4, 0, 2, 
                                  6, 3, 3, // 10
                                  0, 4, 2,
                                  2, 0, 1,
                                  6, 2, 2,
                                  0, 2, 1 }; // 14
  size_t conn[] = { 3, 5, 4, 7,
                    7, 4, 11, 14,
                    5, 0, 8,
                    11, 5, 8,
                    13, 3, 9, 6, 
                    10, 5, 3, 13,
                    12, 7, 14, 1,
                    10, 0, 5,
                    7, 12, 9, 3 };
  size_t conn_len[NUM_ELEMENT] = { 4, 4, 3, 3, 4, 4, 4, 3, 4 };
  EntityTopology types[NUM_ELEMENT] = { QUADRILATERAL,
                                        QUADRILATERAL,
                                        TRIANGLE,
                                        TRIANGLE,
                                        QUADRILATERAL,
                                        QUADRILATERAL,
                                        QUADRILATERAL,
                                        TRIANGLE,
                                        QUADRILATERAL };
    // mark vertices along X and Y axis as fixed
  bool fixed[NUM_VERTEX] = { false, true, true, false, false, false, 
                             true, false, false, true, false, true,
                             true, false, false };

    // populate patch data
  PatchData pd;
  MsqPrintError err(std::cout);
  pd.fill( NUM_VERTEX, coords, NUM_ELEMENT, types, conn_len, conn, fixed, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  
  if (reorder)
    pd.reorder();
    
    // count free vertices
  unsigned i, j;
  size_t num_free = 0;
  for (i = 0; i < NUM_VERTEX; ++i)
    if (!fixed[i])
      ++num_free;
  CPPUNIT_ASSERT_EQUAL( num_free, pd.num_free_vertices() );
  
    // NOTE: PatchData will reorder contents either because reorder()
    //       was called or to group vertices by fixed/free status.
    //       We will assume that the handles arrays for vertices and
    //       elements contain the initial positions in the input
    //       arrays used to populate the patch data.
  
    // Test vertex handles
  std::vector<bool> seen( NUM_VERTEX, false );
  for (i = 0; i < pd.num_nodes(); ++i) {
    size_t h = (size_t)(pd.get_vertex_handles_array()[i]);
    CPPUNIT_ASSERT( h < NUM_VERTEX );
    CPPUNIT_ASSERT( !seen[h] );
    seen[h] = true;
  }
  
    // Test vertex coordinates
  for (i = 0; i < pd.num_nodes(); ++i) {
    size_t h = (size_t)(pd.get_vertex_handles_array()[i]);
    MsqVertex vtx = pd.vertex_by_index(i);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( vtx[0], coords[3*h  ], DBL_EPSILON );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( vtx[1], coords[3*h+1], DBL_EPSILON );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( vtx[2], coords[3*h+2], DBL_EPSILON );
  }
  
    // Test vertex fixed flags
  for (i = 0; i < pd.num_nodes(); ++i) {
    size_t h = (size_t)(pd.get_vertex_handles_array()[i]);
    if (fixed[h]) {
      CPPUNIT_ASSERT( i >= pd.num_free_vertices() );
      CPPUNIT_ASSERT( !pd.vertex_by_index(i).is_free_vertex() );
    }
    else {
      CPPUNIT_ASSERT( i < pd.num_free_vertices() );
      CPPUNIT_ASSERT( pd.vertex_by_index(i).is_free_vertex() );
    }
  }
  
    // Test element handles
  seen.clear();
  seen.resize( NUM_ELEMENT, false );
  for (i = 0; i < pd.num_elements(); ++i) {
    size_t h = (size_t)(pd.get_element_handles_array()[i]);
    CPPUNIT_ASSERT( h < NUM_ELEMENT );
    CPPUNIT_ASSERT( !seen[h] );
    seen[h] = true;
  }
  
    // Test element types 
  for (i = 0; i < pd.num_elements(); ++i) {
    size_t h = (size_t)(pd.get_element_handles_array()[i]);
    CPPUNIT_ASSERT_EQUAL( types[h], pd.element_by_index(i).get_element_type() );
  }
  
    // Test element connectivity 
  for (i = 0; i < pd.num_elements(); ++i) {
    size_t h = (size_t)(pd.get_element_handles_array()[i]);
    MsqMeshEntity& elem = pd.element_by_index(i);
    CPPUNIT_ASSERT_EQUAL( conn_len[h], elem.vertex_count() );
    CPPUNIT_ASSERT_EQUAL( conn_len[h], elem.node_count() );
    
      // calculate offset in input list for element connectivity
    unsigned conn_pos = 0;
    for (j = 0; j < h; ++j)
      conn_pos += conn_len[j];

    const size_t* elem_conn = elem.get_vertex_index_array();
    for (unsigned j = 0; j < elem.vertex_count(); ++j) {
      size_t vh = (size_t)(pd.get_vertex_handles_array()[elem_conn[j]]);
      CPPUNIT_ASSERT_EQUAL( vh, conn[conn_pos] );
      ++conn_pos;
    }
  }
}

void PatchDataTest::test_update_slave_node_coords()
{
  MsqPrintError err(std::cerr);

    // create a patch containing a single 6-node triangle
    // with a) two mid-edge nodes marked as slave vertices and
    // b) with all mid-edge nodes moved away from the center
    // of their corresponding edge.
  const double coords[] = { 0, 0, 0, 
                            1, 0, 0,
                            0, 1, 0,
                            0.5, -0.1, -0.1,
                            0.6,  0.6,  0.1,
                           -0.1,  0.5,  0.0 };
  const size_t init_conn[] = { 0, 1, 2, 3, 4, 5 };
  const bool fixed[] = { false, false, false, false, true, false };
  PatchData pd;
  EntityTopology type = TRIANGLE;
  size_t node_per_tri = 6;
  pd.fill( 6, coords, 1, &type, &node_per_tri, init_conn, fixed, err );
  ASSERT_NO_ERROR(err);
  
    // update_slave_node_coords requires a mapping function
  Settings settings;
//  TriLagrangeShape tri_func;
//  settings.set_mapping_function( &tri_func );
  pd.attach_settings( &settings );
  
    // call the function we're trying to test.
  pd.update_slave_node_coordinates( err );
  ASSERT_NO_ERROR(err);
  
    // get vertex coordinates in the same order that we passed them in
  MsqMeshEntity elem = pd.element_by_index(0);
  const size_t* conn = elem.get_vertex_index_array();
  Vector3D vtx_coords[6];
  for (size_t i = 0; i < 6; ++i)
    vtx_coords[i] = pd.vertex_by_index(conn[i]);
  
    // check that corner vertex coordinates are unchanged
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D( coords   ), vtx_coords[0], 1e-12 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D( coords+3 ), vtx_coords[1], 1e-12 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D( coords+6 ), vtx_coords[2], 1e-12 );
    // check that fixed HO node is unchanged
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D( coords+12 ), vtx_coords[4], 1e-12 );
    // check that slave HO nodes were updated
  const Vector3D mid1 = 0.5 * (Vector3D(coords) + Vector3D(coords+3));
  const Vector3D mid2 = 0.5 * (Vector3D(coords) + Vector3D(coords+6));
  CPPUNIT_ASSERT_VECTORS_EQUAL( mid1, vtx_coords[3], 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( mid2, vtx_coords[5], 1e-6 );
}

/* This is the input mesh topology
     (0)------(16)-----(1)------(17)-----(2)------(18)-----(3)
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
     (19)      0       (20)      1       (21)      2       (22)
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
     (4)------(23)-----(5)------(24)-----(6)------(25)-----(7)
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
     (26)      3       (27)      4       (28)      5       (29)
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
     (8)------(30)-----(9)------(31)-----(10)-----(32)-----(11)
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
     (33)      6       (34)      7       (35)      8       (36)
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
     (12)-----(37)-----(13)-----(38)-----(14)-----(39)-----(15)
*/
    // input mesh definition
  const int NUM_VTX = 40, NUM_ELEM = 9, NUM_CORNER = 16;
  const double input_coords[3*NUM_VTX] = {
        -3.0,  3.0, 0,
        -1.0,  3.0, 0,
         1.0,  3.0, 0,
         3.0,  3.0, 0,
        -3.0,  1.0, 0,
        -1.0,  1.0, 0,
         1.0,  1.0, 0,
         3.0,  1.0, 0,
        -3.0, -1.0, 0,
        -1.0, -1.0, 0,
         1.0, -1.0, 0,
         3.0, -1.0, 0,
        -3.0, -3.0, 0,
        -1.0, -3.0, 0,
         1.0, -3.0, 0,
         3.0, -3.0, 0,
        -2.0,  3.0, 0,
         0.0,  3.0, 0,
         2.0,  3.0, 0,
        -3.0,  2.0, 0,
        -1.0,  2.0, 0,
         1.0,  2.0, 0,
         3.0,  2.0, 0,
        -2.0,  1.0, 0,
         0.0,  1.0, 0,
         2.0,  1.0, 0,
        -3.0,  0.0, 0,
        -1.0,  0.0, 0,
         1.0,  0.0, 0,
         3.0,  0.0, 0,
        -2.0, -1.0, 0,
         0.0, -1.0, 0,
         2.0, -1.0, 0,
        -3.0, -2.0, 0,
        -1.0, -2.0, 0,
         1.0, -2.0, 0,
         3.0, -2.0, 0,
        -2.0, -3.0, 0,
         0.0, -3.0, 0,
         2.0, -3.0, 0 };
  const bool fixed[NUM_VTX] = {
         true,  true,  true, true, 
         true, false, false, true,
         true, false, false, true,
         true,  true,  true, true, 
         true,  true,  true,
         true, false, false, true,
        false, false, false,
         true, false, false, true,
        false, false, false,
         true, false, false, true };
  const size_t input_conn[8*NUM_ELEM] = {
        1,  0,  4,  5, 16, 19, 23, 20,
        2,  1,  5,  6, 17, 20, 24, 21,
        3,  2,  6,  7, 18, 21, 25, 22,
        5,  4,  8,  9, 23, 26, 30, 27,
        6,  5,  9, 10, 24, 27, 31, 28,
        7,  6, 10, 11, 25, 28, 32, 29,
        9,  8, 12, 13, 30, 33, 37, 34,
       10,  9, 13, 14, 31, 34, 38, 35,
       11, 10, 14, 15, 32, 35, 39, 36
  };

  const size_t node_per_elem[NUM_ELEM] = { 8, 8, 8, 8, 8, 8, 8, 8, 8 };
  const EntityTopology elem_types[NUM_ELEM] = {  QUADRILATERAL,
    QUADRILATERAL, QUADRILATERAL, QUADRILATERAL, QUADRILATERAL, 
    QUADRILATERAL, QUADRILATERAL, QUADRILATERAL, QUADRILATERAL };
 
void PatchDataTest::test_quad8_patch( bool reorder, bool slaved )
{
    // create PatchData
  MsqPrintError err( std::cerr );
  PatchData pd;
  Settings settings;
  settings.set_slaved_ho_node_mode(Settings::SLAVE_NONE);
  if (!slaved) // default is slaved, so only change settings if no slaved
    pd.attach_settings(&settings);
  pd.fill( NUM_VTX, input_coords, NUM_ELEM, elem_types, node_per_elem, input_conn, fixed, err );
  ASSERT_NO_ERROR(err);
  
  
    // reorder if testing that
  if (reorder)
    pd.reorder();

    // Check sizes.  Assume that all non-fixed HO nodes are slave vertices
    // unless 'slaved' is false.
    
  CPPUNIT_ASSERT_EQUAL( NUM_VTX, (int)pd.num_nodes() );
  CPPUNIT_ASSERT_EQUAL( NUM_ELEM, (int)pd.num_elements() );
  int num_free = 0, num_slave = 0, num_fixed = 0;
  for (int i = 0; i < NUM_VTX; ++i) {
    if (fixed[i])
      ++num_fixed;
    else if (i < NUM_CORNER)
      ++num_free;
    else if (slaved)
      ++num_slave;
    else
      ++num_free;
  }
  CPPUNIT_ASSERT_EQUAL( num_free , (int)pd.num_free_vertices() );
  CPPUNIT_ASSERT_EQUAL( num_slave, (int)pd.num_slave_vertices() );
  CPPUNIT_ASSERT_EQUAL( num_fixed, (int)pd.num_fixed_vertices() );
  
    // Check that vertex handles and vertex coords are correct.
    // Assume that handles array contains input vertex indices.
    
  for (int i = 0; i < NUM_VTX; ++i) {
    const MsqVertex& vtx = pd.vertex_by_index(i);
    Mesh::VertexHandle hdl = pd.get_vertex_handles_array()[i];
    size_t idx = (size_t)hdl;
    Vector3D exp_coords( input_coords + 3*idx );
    if ((exp_coords - vtx).length_squared() > 1e-16) {
      std::cerr << "Input Index: " << idx << std::endl;
      std::cerr << "Patch Index: " << i << std::endl;
    }
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp_coords, vtx, 1e-16 );
  }
  
    // Check that vertex flags are correct.  
    // Assume all non-fixed HO noes are slave vertices unless 'slaved' is false.
    // Assume that handles array contains input vertex indices.
    
  for (int i = 0; i < NUM_VTX; ++i) {
    const MsqVertex& vtx = pd.vertex_by_index(i);
    Mesh::VertexHandle hdl = pd.get_vertex_handles_array()[i];
    size_t idx = (size_t)hdl;
    if (fixed[idx]) {
      CPPUNIT_ASSERT( vtx.is_flag_set( MsqVertex::MSQ_HARD_FIXED ) );
    }
    else if (slaved && idx >= (size_t)NUM_CORNER) {
      CPPUNIT_ASSERT( vtx.is_flag_set( MsqVertex::MSQ_DEPENDENT ) );
    }
    else {
      CPPUNIT_ASSERT( !vtx.is_flag_set( MsqVertex::MSQ_HARD_FIXED ) );
      CPPUNIT_ASSERT( !vtx.is_flag_set( MsqVertex::MSQ_DEPENDENT  ) );
    }
  }
  
    // Check that element connectivity is correct.
    // Assume that handles array contains input vertex and element indices.
  
  for (int i = 0; i < NUM_ELEM; ++i) {
    MsqMeshEntity& elem = pd.element_by_index(i);
    CPPUNIT_ASSERT_EQUAL( QUADRILATERAL, elem.get_element_type() );
    CPPUNIT_ASSERT_EQUAL( (size_t)4, elem.vertex_count() );
    CPPUNIT_ASSERT_EQUAL( (size_t)8,elem.node_count() );
    CPPUNIT_ASSERT_EQUAL( (size_t)4, elem.corner_count() );
    std::vector<std::size_t> conn;
    elem.get_node_indices( conn );
    CPPUNIT_ASSERT_EQUAL( elem.node_count(), conn.size() );
    for (int j = 0; j < 8; ++j)
      conn[j] = (size_t)pd.get_vertex_handles_array()[conn[j]];
    size_t idx = (size_t)pd.get_element_handles_array()[i];
    ASSERT_ARRAYS_EQUAL( input_conn + 8*idx, arrptr(conn), 8 );
  }
}

void PatchDataTest::get_quad8_mesh( Mesh*& mesh_out )
{
  static std::vector<int> fixed_flags(fixed, fixed+NUM_VTX);
  static std::vector<double> coords(input_coords, input_coords+3*NUM_VTX);
  static std::vector<unsigned long> conn( input_conn, input_conn+8*NUM_ELEM );
  mesh_out = new ArrayMesh( 3, NUM_VTX, arrptr(coords), arrptr(fixed_flags), 
                  NUM_ELEM, QUADRILATERAL, arrptr(conn), false,
                  8 );
}

void PatchDataTest::get_quad8_mesh_and_domain( Mesh*& mesh_out,
                                               MeshDomain*& domain_out )
{
  MsqPrintError err(std::cerr);

  get_quad8_mesh( mesh_out );
  DomainClassifier geom;
  Vector3D corners[4] = {  Vector3D(-3, 3, 0),
                           Vector3D( 3, 3, 0),
                           Vector3D(-3,-3, 0),
                           Vector3D( 3,-3, 0) };
  MeshDomain* geomarr[] = { new PointDomain( corners[0] ),
                            new PointDomain( corners[1] ),
                            new PointDomain( corners[2] ),
                            new PointDomain( corners[3] ),
                            new LineDomain( corners[0], corners[1] - corners[0] ),
                            new LineDomain( corners[1], corners[2] - corners[1] ),
                            new LineDomain( corners[2], corners[3] - corners[2] ),
                            new LineDomain( corners[3], corners[0] - corners[3] ),
                            new PlanarDomain( PlanarDomain::XY ) };
  int dimarr[] = { 0, 0, 0, 0, 1, 1, 1, 1, 2 };
  
  DomainClassifier* domain;
  domain_out = domain = new DomainClassifier;
  DomainClassifier::classify_geometrically( *domain, mesh_out, 1e-6, geomarr, dimarr, 9, err );
  domain->delete_sub_domains(true);
  ASSERT_NO_ERROR(err);
}

void PatchDataTest::test_fixed_by_geom_dim( unsigned dim )
{
  MsqPrintError err(std::cerr);

  Settings settings;
  switch (dim) {
    case 0: settings.set_fixed_vertex_mode( Settings::FIXED_VERTEX  ); break;
    case 1: settings.set_fixed_vertex_mode( Settings::FIXED_CURVE   ); break;
    case 2: settings.set_fixed_vertex_mode( Settings::FIXED_SURFACE ); break;
    default: CPPUNIT_ASSERT(false);
  }
  
  Mesh* mesh = 0;
  MeshDomain* domain = 0;
  get_quad8_mesh_and_domain( mesh, domain );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(mesh, domain);
  Instruction::initialize_vertex_byte( &mesh_and_domain, &settings, err ); 
  ASSERT_NO_ERROR(err);
  
  PatchData pd;
  pd.attach_settings( &settings );
  pd.set_mesh( mesh );
  pd.set_domain( domain );
  
  std::vector<Mesh::ElementHandle> elems;
  std::vector<Mesh::VertexHandle> verts;
  mesh->get_all_elements( elems, err ); ASSERT_NO_ERROR(err);
  mesh->get_all_vertices( verts, err ); ASSERT_NO_ERROR(err);
  pd.set_mesh_entities( elems, verts, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!elems.empty());
  
  std::vector<unsigned short> dims(verts.size());
  domain->domain_DoF( arrptr(verts), arrptr(dims), verts.size(), err );
  ASSERT_NO_ERROR(err);
  
  for (size_t i = 0; i < pd.num_free_vertices(); ++i) {
    Mesh::VertexHandle handle = pd.get_vertex_handles_array()[i];
    unsigned short d;
    domain->domain_DoF( &handle, &d, 1, err ); ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT( d > dim );
  }
  for (size_t i = 0; i < pd.num_fixed_vertices(); ++i) {
    size_t j = i + pd.num_free_vertices() + pd.num_slave_vertices();
    Mesh::VertexHandle handle = pd.get_vertex_handles_array()[j];
    unsigned short d;
    domain->domain_DoF( &handle, &d, 1, err ); ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT( d <= dim );
  }
  
  delete mesh;
  delete domain;
}

void PatchDataTest::get_higher_order_vertices( Mesh* mesh, 
                                           std::map<Mesh::VertexHandle,bool>& ho_verts,
                                           bool initial_value,
                                           bool non_fixed_only )
{
    // get mesh data
  MsqPrintError err(std::cerr);
  std::vector<Mesh::ElementHandle> elems;
  std::vector<Mesh::VertexHandle> verts;
  std::vector<size_t> offsets;
  mesh->get_all_elements( elems, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!elems.empty());
  mesh->elements_get_attached_vertices( arrptr(elems), elems.size(),
                                        verts, offsets, err );
  CPPUNIT_ASSERT_EQUAL( elems.size()+1, offsets.size() );
  ASSERT_NO_ERROR(err);
  std::vector<EntityTopology> types(elems.size());
  mesh->elements_get_topologies( arrptr(elems), arrptr(types), elems.size(), err );
  ASSERT_NO_ERROR(err);
  
    // clear initial state
  ho_verts.clear();
  
    // for each element, add ho nodes
  for (size_t i = 0; i < elems.size(); ++i) 
    for (size_t j = offsets[i] + TopologyInfo::corners(types[i]); 
         j < offsets[i+1]; ++j) 
      ho_verts[verts[j]] = initial_value;
      
  if (non_fixed_only) {
    std::map<Mesh::VertexHandle,bool>::iterator p;
    std::sort( verts.begin(), verts.end() );
    verts.erase( std::unique( verts.begin(), verts.end() ), verts.end() );
    std::vector<bool> fixed;
    mesh->vertices_get_fixed_flag( arrptr(verts), fixed, verts.size(), err );
    ASSERT_NO_ERROR(err);
    for (size_t i = 0; i < verts.size(); ++i) {
      if (fixed[i]) {
        p = ho_verts.find(verts[i]);
        if (p != ho_verts.end())
          ho_verts.erase(p);
      }
    }
  }
}

void PatchDataTest::check_higher_order_vertices_slaved( 
                   Mesh* mesh,
                   Settings::HigherOrderSlaveMode mode,
                   const std::map<Mesh::VertexHandle,bool>& expected )
{
  MsqPrintError err(std::cerr);
  
  Settings settings;
  settings.set_slaved_ho_node_mode( mode );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(mesh, 0);
  Instruction::initialize_vertex_byte( &mesh_and_domain, &settings, err ); 
  ASSERT_NO_ERROR(err);
  
  PatchData pd;
  pd.attach_settings( &settings );
  pd.set_mesh( mesh );
  
  std::vector<Mesh::ElementHandle> elements;
  std::vector<Mesh::VertexHandle> vertices;
  mesh->get_all_elements( elements, err ); 
  ASSERT_NO_ERROR(err);
  pd.set_mesh_entities( elements, vertices, err ); 
  ASSERT_NO_ERROR(err);
  
  std::map<Mesh::VertexHandle,bool>::const_iterator p;
  for (size_t i = 0; i < pd.num_nodes(); ++i) {
    p = expected.find( pd.get_vertex_handles_array()[i] );
    bool found = (p != expected.end());
    bool exp = found && p->second;
    bool act = pd.vertex_by_index(i).is_flag_set(MsqVertex::MSQ_DEPENDENT);
    CPPUNIT_ASSERT_EQUAL( exp, act );
  }
}                

void PatchDataTest::test_patch_data_mesh_slaved_ho_nodes()
{
  Mesh* mesh = 0;
  get_quad8_mesh( mesh );
  
  std::map<Mesh::VertexHandle,bool> ho_verts;
  get_higher_order_vertices( mesh, ho_verts, true, false );
  
  check_higher_order_vertices_slaved( mesh, Settings::SLAVE_ALL, ho_verts );
  delete mesh;
}

void PatchDataTest::test_patch_data_mesh_free_ho_nodes()
{
  Mesh* mesh = 0;
  get_quad8_mesh( mesh );
  
  std::map<Mesh::VertexHandle,bool> ho_verts;
  get_higher_order_vertices( mesh, ho_verts, false );
  
  check_higher_order_vertices_slaved( mesh, Settings::SLAVE_NONE, ho_verts );
  delete mesh;
}

void PatchDataTest::test_patch_data_mesh_calcualted_ho_nodes()
{
  Mesh* mesh = 0;
  get_quad8_mesh( mesh );
  
  std::map<Mesh::VertexHandle,bool> ho_verts;
  get_higher_order_vertices( mesh, ho_verts, false );
  
    // set bit on every other higher-order vertex
  std::map<Mesh::VertexHandle,bool>::iterator i = ho_verts.end();
  std::vector<Mesh::VertexHandle> slaved;
  std::vector<unsigned char> bytes;
  while (i != ho_verts.end()) {
    slaved.push_back(i->first);
    bytes.push_back( MsqVertex::MSQ_DEPENDENT );
    i->second = true;
    if (++i == ho_verts.end()); 
      break;
    ++i;
  }
  
  if (!slaved.empty()) {
    MsqPrintError err(std::cerr);
    mesh->vertices_set_byte( arrptr(slaved), arrptr(bytes), slaved.size(), err );
    ASSERT_NO_ERROR(err);
  }
  
  check_higher_order_vertices_slaved( mesh, Settings::SLAVE_CALCULATED, ho_verts );
  delete mesh;
}

//! create a wrapper around the real mesh that returns what we
//! want from vertices_get_slaved_flag.
class HoSlavedMesh : public MeshDecorator
{
  public:
    typedef std::map<Mesh::VertexHandle,bool> SMap;
    HoSlavedMesh( Mesh* real_mesh, SMap& slaved ) 
      : MeshDecorator( real_mesh ),
        slavedVerts(slaved)
      {}
      
    virtual void vertices_get_slaved_flag( const VertexHandle vert_array[], 
                                           std::vector<bool>& slaved_flag_array,
                                           size_t num_vtx, 
                                           MsqError &err );
  private:
    SMap slavedVerts;
};
void HoSlavedMesh::vertices_get_slaved_flag( const VertexHandle vert_array[], 
                                             std::vector<bool>& slaved_flag_array,
                                             size_t num_vtx, 
                                             MsqError &err )
{
  slaved_flag_array.resize( num_vtx );
  for (size_t i = 0; i < num_vtx; ++i) {
    SMap::iterator j = slavedVerts.find( vert_array[i] );
    slaved_flag_array[i] = (j != slavedVerts.end()) && j->second;
  }
}

void PatchDataTest::test_patch_data_mesh_flagged_ho_nodes()
{
  Mesh* mesh = 0;
  get_quad8_mesh( mesh );
  
  std::map<Mesh::VertexHandle,bool> ho_verts;
  get_higher_order_vertices( mesh, ho_verts, false );
  
    // set every other higher-order vertex as slaved
  std::map<Mesh::VertexHandle,bool>::iterator i = ho_verts.end();
  while (i != ho_verts.end()) {
    if (++i == ho_verts.end()); 
      break;
    i->second = true;
    ++i;
  }
  
    // create a wrapper mesh that returns what we want 
    // from vertices_get_slaved_flag
  HoSlavedMesh wrapper( mesh, ho_verts );
  
  check_higher_order_vertices_slaved( &wrapper, Settings::SLAVE_FLAG, ho_verts );
  delete mesh;
}
