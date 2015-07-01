/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TagVertexMeshTest.cpp
 *  \brief unit tests for TagVertexMesh class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TagVertexMesh.hpp"
#include "MsqError.hpp"
#include "UnitUtil.hpp"
#include "MeshImpl.hpp"
#include "MsqVertex.hpp"
#include "InstructionQueue.hpp"
#include <cppunit/extensions/HelperMacros.h>
#include <stdio.h>

using namespace Mesquite;

const char TEMP_FILE_NAME[] = "TagVertexMesh.vtk";

class TagVertexMeshTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(TagVertexMeshTest);
  CPPUNIT_TEST (test_vertex_coordinates);
  CPPUNIT_TEST (test_save_coordinates);
  CPPUNIT_TEST (test_cleanup);
  CPPUNIT_TEST (test_alternate_name);
  CPPUNIT_TEST (test_reference_mesh);
  CPPUNIT_TEST_SUITE_END();

  MeshImpl* realMesh;

public:

  TagVertexMeshTest() : realMesh(0) {}

  void setUp();
  void tearDown();
  
  void test_vertex_coordinates();
  void test_save_coordinates();
  void test_cleanup();
  void test_alternate_name();
  void test_reference_mesh();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TagVertexMeshTest, "TagVertexMeshTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TagVertexMeshTest, "Unit");

void TagVertexMeshTest::setUp()
{
  const char vtk_data[] = 
    "# vtk DataFile Version 2.0\n"
    "test mesh\n"
    "ASCII\n"
    "DATASET UNSTRUCTURED_GRID\n"
    "POINTS 3 float\n"
    "0 0 0\n"
    "1 0 0\n"
    "0 1 0\n"
    "CELLS 1 4\n"
    "3 0 1 2\n"
    "CELL_TYPES 1\n"
    "5\n";
  
  FILE* file = fopen( TEMP_FILE_NAME, "w" );
  CPPUNIT_ASSERT( !!file );
  size_t r = fwrite( vtk_data, sizeof(vtk_data)-1, 1, file );
  fclose( file );
  CPPUNIT_ASSERT( r == 1 );
  
  MsqPrintError err( std::cerr );
  realMesh = new MeshImpl;
  realMesh->read_vtk( TEMP_FILE_NAME, err );
  remove( TEMP_FILE_NAME );
  ASSERT_NO_ERROR(err);
}

void TagVertexMeshTest::tearDown()
{
  delete realMesh;
  realMesh = 0;
}

void TagVertexMeshTest::test_vertex_coordinates()
{
  MsqPrintError err( std::cerr );
  TagVertexMesh tag_mesh( err, realMesh, true );
  ASSERT_NO_ERROR(err);
  
  std::vector<Mesh::VertexHandle> vertices;
  realMesh->get_all_vertices( vertices, err );
  ASSERT_NO_ERROR(err);
  
    // Check that initial position for vertex matches that of real mesh
  Mesh::VertexHandle vertex = vertices[0];
  MsqVertex get_coords;
  Vector3D orig_coords, real_coords, tag_coords;
  realMesh->vertices_get_coordinates( &vertex, &get_coords, 1, err );
  ASSERT_NO_ERROR(err);
  orig_coords = get_coords;
  tag_mesh.vertices_get_coordinates( &vertex, &get_coords, 1, err );
  ASSERT_NO_ERROR(err);
  tag_coords = get_coords;
  CPPUNIT_ASSERT_VECTORS_EQUAL( orig_coords, tag_coords, DBL_EPSILON );
  
    // Check that modified vertex coords show up in tag mesh but not
    // real mesh.
  Vector3D new_coords(5,5,5);
  tag_mesh.vertex_set_coordinates( vertex, new_coords, err );
  ASSERT_NO_ERROR(err);
  tag_mesh.vertices_get_coordinates( &vertex, &get_coords, 1, err );
  ASSERT_NO_ERROR(err);
  tag_coords = get_coords;
  CPPUNIT_ASSERT_VECTORS_EQUAL( new_coords, tag_coords, DBL_EPSILON );
  realMesh->vertices_get_coordinates( &vertex, &get_coords, 1, err );
  ASSERT_NO_ERROR(err);
  real_coords = get_coords;
  CPPUNIT_ASSERT_VECTORS_EQUAL( orig_coords, real_coords, DBL_EPSILON );
}


void TagVertexMeshTest::test_save_coordinates()
{
  MsqPrintError err( std::cerr );
  Vector3D new_coords(5, 5, 5);
  MsqVertex get_coords;
  
  std::vector<Mesh::VertexHandle> vertices;
  realMesh->get_all_vertices( vertices, err );
  ASSERT_NO_ERROR(err);
  Mesh::VertexHandle vertex = vertices[0];

  realMesh->vertices_get_coordinates( &vertex, &get_coords, 1, err );
  ASSERT_NO_ERROR(err);
  Vector3D orig_coords = get_coords;

    // modify a vertex in the tag interface
  {
    TagVertexMesh tag_mesh( err, realMesh, false );
    ASSERT_NO_ERROR(err);
    tag_mesh.vertex_set_coordinates( vertex, new_coords, err );
    ASSERT_NO_ERROR(err);
  }

    // check that it exists in a new TagVertexMesh
  {
    TagVertexMesh tag_mesh( err, realMesh, false );
    ASSERT_NO_ERROR(err);
    tag_mesh.vertices_get_coordinates( &vertex, &get_coords, 1, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_VECTORS_EQUAL( new_coords, get_coords, DBL_EPSILON );
  }
}

void TagVertexMeshTest::test_cleanup()
{
  MsqPrintError err( std::cerr );
  Vector3D new_coords(5, 5, 5);
  MsqVertex get_coords;
  
  std::vector<Mesh::VertexHandle> vertices;
  realMesh->get_all_vertices( vertices, err );
  ASSERT_NO_ERROR(err);
  Mesh::VertexHandle vertex = vertices[0];

  realMesh->vertices_get_coordinates( &vertex, &get_coords, 1, err );
  ASSERT_NO_ERROR(err);
  Vector3D orig_coords = get_coords;

    // modify a vertex in the tag interface
  {
    TagVertexMesh tag_mesh( err, realMesh, true );
    ASSERT_NO_ERROR(err);
    tag_mesh.vertex_set_coordinates( vertex, new_coords, err );
    ASSERT_NO_ERROR(err);
  }

    // check that values were cleaned up when previous instance was destroyed
  {
    TagVertexMesh tag_mesh( err, realMesh, false );
    ASSERT_NO_ERROR(err);
    tag_mesh.vertices_get_coordinates( &vertex, &get_coords, 1, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_VECTORS_EQUAL( orig_coords, get_coords, DBL_EPSILON );
  }
}

void TagVertexMeshTest::test_alternate_name()
{
  MsqPrintError err( std::cerr );
  Vector3D new_coords(5, 5, 5);
  MsqVertex get_coords;
  
  std::vector<Mesh::VertexHandle> vertices;
  realMesh->get_all_vertices( vertices, err );
  ASSERT_NO_ERROR(err);
  Mesh::VertexHandle vertex = vertices[0];

  realMesh->vertices_get_coordinates( &vertex, &get_coords, 1, err );
  ASSERT_NO_ERROR(err);
  Vector3D orig_coords = get_coords;

    // modify a vertex in the tag interface and save it
  {
    TagVertexMesh tag_mesh( err, realMesh, false, "foobar" );
    ASSERT_NO_ERROR(err);
    tag_mesh.vertex_set_coordinates( vertex, new_coords, err );
    ASSERT_NO_ERROR(err);
  }
  
    // verify that it is modified in the new interface
  {
    TagVertexMesh tag_mesh( err, realMesh, false, "foobar" );
    ASSERT_NO_ERROR(err);
    tag_mesh.vertices_get_coordinates( &vertex, &get_coords, 1, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_VECTORS_EQUAL( new_coords, get_coords, DBL_EPSILON );
  }
}


void TagVertexMeshTest::test_reference_mesh()
{
  MsqPrintError err( std::cerr );
  TagVertexMesh tag_mesh( err, realMesh, true );
  ASSERT_NO_ERROR(err);
  
  std::vector<Mesh::VertexHandle> vertices;
  realMesh->get_all_vertices( vertices, err );
  ASSERT_NO_ERROR(err);
  
    // copy real mesh coordinates into tag data in TagVertexMesh
  InstructionQueue q;
  q.add_tag_vertex_mesh( &tag_mesh, err );
  ASSERT_NO_ERROR(err);
  q.run_instructions( realMesh, err );
  ASSERT_NO_ERROR(err);
  
    // Check that initial position for vertex matches that of real mesh
  Mesh::VertexHandle vertex = vertices[0];
  MsqVertex get_coords;
  Vector3D orig_coords, real_coords, tag_coords;
  realMesh->vertices_get_coordinates( &vertex, &get_coords, 1, err );
  ASSERT_NO_ERROR(err);
  orig_coords = get_coords;
  tag_mesh.vertices_get_coordinates( &vertex, &get_coords, 1, err );
  ASSERT_NO_ERROR(err);
  tag_coords = get_coords;
  CPPUNIT_ASSERT_VECTORS_EQUAL( orig_coords, tag_coords, DBL_EPSILON );
  
    // Check that modified vertex coords show up in real mesh but not
    // tag mesh.
  realMesh->vertices_get_coordinates( &vertex, &get_coords, 1, err );
  ASSERT_NO_ERROR(err);
  orig_coords = get_coords;
  Vector3D new_coords(5,5,5);
  realMesh->vertex_set_coordinates( vertex, new_coords, err );
  ASSERT_NO_ERROR(err);
  tag_mesh.vertices_get_coordinates( &vertex, &get_coords, 1, err );
  ASSERT_NO_ERROR(err);
  tag_coords = get_coords;
  CPPUNIT_ASSERT_VECTORS_EQUAL( orig_coords, tag_coords, DBL_EPSILON );
    // restore realMesh to initial state
  realMesh->vertex_set_coordinates( vertex, orig_coords, err );
  ASSERT_NO_ERROR(err);
}
