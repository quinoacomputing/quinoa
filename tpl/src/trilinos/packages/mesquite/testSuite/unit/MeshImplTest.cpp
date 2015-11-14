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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov,
    kraftche@cae.wisc.edu
   
  ***************************************************************** */
/**
 *\file   MeshImplTest.hpp
 *\brief  Test MeshImpl class
 *
 * Test MeshImpl functionality.  
 *
 * Note: there are additional tests for MeshImpl in the mis-named 
 *       MeshInterfaceTest.cpp, VtkTest.cpp, and presumably ExodusTest.cpp.
 *
 *\author Jason Kraftcheck
 *\date   2007-10-18
*/

#include "Mesquite.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_MsqVertex.hpp"

#include "UnitUtil.hpp"

#include <iostream>
#include <stdio.h>

using namespace Mesquite;

const bool DUMP_MESH = false;

class MeshImplTest : public CppUnit::TestFixture
{
public:

  CPPUNIT_TEST_SUITE( MeshImplTest );
  CPPUNIT_TEST( test_zero_length_data );
  CPPUNIT_TEST( skin_mesh_2D );
  CPPUNIT_TEST( skin_mesh_3D );
  CPPUNIT_TEST( skin_mesh_mixed );
  CPPUNIT_TEST( skin_mesh_higher_order );
  CPPUNIT_TEST_SUITE_END();

  void test_zero_length_data();
  void skin_mesh_2D();
  void skin_mesh_3D();
  void skin_mesh_mixed();
  void skin_mesh_higher_order();
  
  static void load_vtk( const char* file_data, MeshImpl& mesh, MsqError& err );
  static void dump_mesh( const char* filename, MeshImpl& data, MsqError& err );
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MeshImplTest, "MeshImplTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MeshImplTest, "Unit");


void MeshImplTest::load_vtk( const char* file_data, MeshImpl& mesh, MsqError& err )
{
  const char* fname = "MeshImplTest.vtk";
  FILE* f = fopen( fname, "w" );
  if (!f) {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS, "Cannot create temp file: %s\n", fname);
    return;
  }
  
  if (!fwrite(file_data, strlen(file_data), 1, f)) {
    MSQ_SETERR(err)("I/O error while writing temporary file\n", MsqError::IO_ERROR);
    fclose( f );
    remove( fname );
    return;
  }
  
  fclose( f );
  mesh.read_vtk( fname, err ); MSQ_CHKERR(err);
  remove( fname );  
  
  mesh.mark_skin_fixed( err );
  MSQ_CHKERR(err);
}

void MeshImplTest::dump_mesh( const char* filename, MeshImpl& mesh, MsqError& err )
{
  if (DUMP_MESH) {
    mesh.write_vtk( filename, err );
    if (MSQ_CHKERR(err))
      std::cerr << err << std::endl;
  }
}


void MeshImplTest::test_zero_length_data()
{
  const size_t num_vtx = 2, zero = 0;
  const double coords[3*num_vtx] = { 0, 0, 0, 1, 1, 1 };
  const bool fixed[num_vtx] = { false, false };
  std::vector<Mesh::ElementHandle> elems;
  std::vector<Mesh::VertexHandle> verts;
  MsqError err;
  const int conn[6] = { 0, 1, 0, 1, 0, 1 };
  EntityTopology type = TRIANGLE;

  MeshImpl no_elem1( 2, 0, TRIANGLE, fixed, coords, 0 );
  verts.clear();
  no_elem1.get_all_vertices( verts, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( num_vtx, verts.size() );
  elems.clear();
  no_elem1.get_all_elements( elems, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( zero, elems.size() );

  MeshImpl no_elem2( 2, 0, 0, fixed, coords, 0 );
  verts.clear();
  no_elem2.get_all_vertices( verts, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( num_vtx, verts.size() );
  elems.clear();
  no_elem2.get_all_elements( elems, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( zero, elems.size() );

  MeshImpl no_elem3( 2, 0, TRIANGLE, fixed, coords, conn );
  verts.clear();
  no_elem3.get_all_vertices( verts, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( num_vtx, verts.size() );
  elems.clear();
  no_elem3.get_all_elements( elems, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( zero, elems.size() );

  MeshImpl no_elem4( 2, 0, &type, fixed, coords, conn );
  verts.clear();
  no_elem4.get_all_vertices( verts, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( num_vtx, verts.size() );
  elems.clear();
  no_elem4.get_all_elements( elems, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( zero, elems.size() );
  
  MeshImpl no_vert1( 0, 0, TRIANGLE, 0, 0, 0 );
  verts.clear();
  no_vert1.get_all_vertices( verts, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( zero, verts.size() );
  elems.clear();
  no_vert1.get_all_elements( elems, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( zero, elems.size() );

  MeshImpl no_vert2( 0, 0, 0, 0, 0, 0 );
  verts.clear();
  no_vert2.get_all_vertices( verts, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( zero, verts.size() );
  elems.clear();
  no_vert2.get_all_elements( elems, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( zero, elems.size() );
  
  MeshImpl no_vert3( 0, 0, TRIANGLE, fixed, coords, 0 );
  verts.clear();
  no_vert3.get_all_vertices( verts, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( zero, verts.size() );
  elems.clear();
  no_vert3.get_all_elements( elems, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( zero, elems.size() );

  MeshImpl no_vert4( 0, 0, 0, fixed, coords, 0 );
  verts.clear();
  no_vert4.get_all_vertices( verts, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( zero, verts.size() );
  elems.clear();
  no_vert4.get_all_elements( elems, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( zero, elems.size() );
}

void MeshImplTest::skin_mesh_2D()
{
  MsqPrintError err(std::cerr);
  const char vtk_file[] = 
    "#vtk DataFile Version 2.0\n"
    "test data for MeshImplTest::skin_mesh_2D\n"
    "ASCII\n"
    "DATASET STRUCTURED_POINTS\n"
    "DIMENSIONS 10 10 1\n"
    "ORIGIN 0 0 0\n"
    "SPACING 1 1 1\n";
  
  MeshImpl mesh;
  load_vtk( vtk_file, mesh, err );
  CPPUNIT_ASSERT(!err);
  dump_mesh( "MeshSkin2D.vtk", mesh, err );
  
  std::vector<Mesh::VertexHandle> verts;
  mesh.get_all_vertices( verts, err );
  CPPUNIT_ASSERT(!err);
  
  std::vector<Mesh::VertexHandle> elems;
  std::vector<size_t> offsets;
  std::vector<bool> fixed;
  
  
  for (unsigned i = 0; i < verts.size(); ++i) {
    elems.clear();
    offsets.clear();
    mesh.vertices_get_attached_elements( &verts[i], 1, elems, offsets, err );
    CPPUNIT_ASSERT(!err);
    
    mesh.vertices_get_fixed_flag( &verts[i], fixed, 1, err );
    CPPUNIT_ASSERT(!err);
    
    if (elems.size() == 4)
      CPPUNIT_ASSERT(!fixed[0]);
    else 
      CPPUNIT_ASSERT(fixed[0]);
  }
}

void MeshImplTest::skin_mesh_3D()
{
  MsqPrintError err(std::cerr);
  const char vtk_file[] = 
    "#vtk DataFile Version 2.0\n"
    "test data for MeshImplTest::skin_mesh_3D\n"
    "ASCII\n"
    "DATASET STRUCTURED_POINTS\n"
    "DIMENSIONS 10 10 10\n"
    "ORIGIN 0 0 0\n"
    "SPACING 1 1 1\n";
  
  MeshImpl mesh;
  load_vtk( vtk_file, mesh, err );
  CPPUNIT_ASSERT(!err);
  dump_mesh( "MeshSkin3D.vtk", mesh, err );
  
  std::vector<Mesh::VertexHandle> verts;
  mesh.get_all_vertices( verts, err );
  CPPUNIT_ASSERT(!err);
  
  std::vector<Mesh::VertexHandle> elems;
  std::vector<size_t> offsets;
  std::vector<bool> fixed;
  
  for (unsigned i = 0; i < verts.size(); ++i) {
    elems.clear();
    offsets.clear();
    mesh.vertices_get_attached_elements( &verts[i], 1, elems, offsets, err );
    CPPUNIT_ASSERT(!err);
    
    mesh.vertices_get_fixed_flag( &verts[i], fixed, 1, err );
    CPPUNIT_ASSERT(!err);
    
    if (elems.size() == 8)
      CPPUNIT_ASSERT(!fixed[0]);
    else 
      CPPUNIT_ASSERT(fixed[0]);
  }
}

// Make sure that if we have a volume mesh with one interior
// face that the vertices on that face are not marked as fixed
// due soley to the fact that there are no adjacent faces.
// That is, the iterior face doesn't 'count' because it is the
// side of a volume element.
void MeshImplTest::skin_mesh_mixed()
{
  MsqPrintError err(std::cerr);
    // define the mesh of a tetrahedron centered
    // at the origin as four tetrahedral elements
    // sharing a vertex at the origin.
  const char vtk_file[] = 
    "#vtk DataFile Version 2.0\n"
    "test data for MeshImplTest::skin_mesh_mixed\n"
    "ASCII\n"
    "DATASET UNSTRUCTURED_GRID\n"
    "POINTS 5 float\n"
    " 0  0  0\n"  // center vertex
    " 2 -1 -1\n"  
    " 0  2 -1\n"
    "-2 -1 -1\n"
    " 0  0  2\n"
    "CELLS 5 24\n"
    "4 1 2 3 0\n"
    "4 1 3 4 0\n"
    "4 4 2 1 0\n"
    "4 4 3 2 0\n"
    "3 0 1 2\n"   // iterior triangle
    "CELL_TYPES 5\n"
    "10 10 10 10 5\n";
  
  MeshImpl mesh;
  load_vtk( vtk_file, mesh, err );
  CPPUNIT_ASSERT(!err);
  dump_mesh( "MeshSkinMixed.vtk", mesh, err );
  
  std::vector<Mesh::VertexHandle> verts;
  mesh.get_all_vertices( verts, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL(5u, (unsigned)verts.size());
  MsqVertex coords[5];
  mesh.vertices_get_coordinates( arrptr(verts), coords, 5, err );
  CPPUNIT_ASSERT(!err);
  std::vector<bool> fixed;
  mesh.vertices_get_fixed_flag( arrptr(verts), fixed, 5, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( (size_t)5, fixed.size() );
  
  int free_idx = -1;
  for (int i = 0; i < 5; ++i) {
    if (!fixed[i]) {
      CPPUNIT_ASSERT_EQUAL( -1, free_idx ); // at most 1 free vertex
      free_idx = i;
    }
  }
  CPPUNIT_ASSERT( -1 != free_idx ); // at least one free vertex
  
  // free vertex must be at origin
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,0,0), coords[free_idx], DBL_EPSILON );
}


void MeshImplTest::skin_mesh_higher_order()
{
  MsqPrintError err(std::cerr);
    // define a simple 2-quad mesh:
    //
    //  0    1    2    3    4
    //  o----o----o----o----o
    //  |         |         |
    //  |         |         | 
    //  o5        o6        o7
    //  |         |         |
    //  |         |         | 
    //  o----o----o----o----o
    //  8    9    10   11   12
    //
    // Expect all vertices but 6 to be fixed.
    // Position mesh such that vertex 6 is at the origin.
    
  const char vtk_file[] = 
    "#vtk DataFile Version 2.0\n"
    "test data for MeshImplTest::skin_mesh_mixed\n"
    "ASCII\n"
    "DATASET UNSTRUCTURED_GRID\n"
    "POINTS 13 float\n"
    "-2  1  0\n"
    "-1  1  0\n"
    " 0  1  0\n"
    " 1  1  0\n"
    " 2  1  0\n"
    "-2  0  0\n"
    " 0  0  0\n"
    " 2  0  0\n"
    "-2 -1  0\n"
    "-1 -1  0\n"
    " 0 -1  0\n"
    " 1 -1  0\n"
    " 2 -1  0\n"
    "CELLS 2 18\n"
    "8 0 8 10 2 5 9 6 1\n"
    "8 2 10 12 4 6 11 7 3\n"
    "CELL_TYPES 2\n"
    "23 23\n";
  
  MeshImpl mesh;
  load_vtk( vtk_file, mesh, err );
  CPPUNIT_ASSERT(!err);
  dump_mesh( "MeshSkinHO.vtk", mesh, err );
  
  std::vector<Mesh::VertexHandle> verts;
  mesh.get_all_vertices( verts, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL(13u, (unsigned)verts.size());
  MsqVertex coords[13];
  mesh.vertices_get_coordinates( arrptr(verts), coords, 13, err );
  CPPUNIT_ASSERT(!err);
  std::vector<bool> fixed;
  mesh.vertices_get_fixed_flag( arrptr(verts), fixed, 13, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( (size_t)13, fixed.size() );
  
  int free_idx = -1;
  for (int i = 0; i < 13; ++i) {
    if (!fixed[i]) {
      CPPUNIT_ASSERT_EQUAL( -1, free_idx ); // at most 1 free vertex
      free_idx = i;
    }
  }
  CPPUNIT_ASSERT( -1 != free_idx ); // at least one free vertex
  
  // free vertex must be at origin
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,0,0), coords[free_idx], DBL_EPSILON );
}


