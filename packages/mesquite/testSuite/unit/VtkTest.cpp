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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//

#include "meshfiles.h"

#include <iostream>
using std::cout;

#ifndef DEBUG
#  include "Mesquite.hpp"
#  include "PatchData.hpp"
#  include "MeshImpl.hpp"
#  include "MeshInterface.hpp"
#  include "VertexPatches.hpp"
#  include "PatchIterator.hpp"
#  include "UnitUtil.hpp"
#  include "Instruction.hpp"
   using Mesquite::Mesh;
   using Mesquite::MeshImpl;
   using Mesquite::Vector3D;
   using Mesquite::MsqVertex;
   using Mesquite::MsqPrintError;
   using Mesquite::EntityTopology;
   using Mesquite::VertexPatches;
   using Mesquite::PatchIterator;
   using Mesquite::arrptr;
#else
#  include <stdio.h>
#endif

#include <algorithm>


extern const char temp_file_name[] = "VtkTest.vtk";

  // VTK file for 2x2x2 block of hexes as structured-point 
extern const char structured_3d_points_data[] = 
"# vtk DataFile Version 2.0\n"
"Mesquite Mesh\n"
"ASCII\n"
"DATASET STRUCTURED_POINTS\n"
"DIMENSIONS 3 3 3\n"
"ORIGIN 0 0 0\n"
"SPACING 1.5 1.5 1.5\n";

  // VTK file for 2x2 block of quads as structured-point 
extern const char structured_2d_points_data[] = 
"# vtk DataFile Version 2.0\n"
"Mesquite Mesh\n"
"ASCII\n"
"DATASET STRUCTURED_POINTS\n"
"DIMENSIONS 3 3 1\n"
"ORIGIN 0 0 0\n"
"SPACING 1.5 1.5 0.0\n";

  // VTK file for 2x2x2 block of hexes as structured-grid
extern const char structured_grid_data[] = 
"# vtk DataFile Version 2.0\n"
"Mesquite Mesh\n"
"ASCII\n"
"DATASET STRUCTURED_GRID\n"
"DIMENSIONS 3 3 3\n"
"POINTS 27 float\n"
"0.0 0.0 0.0\n1.5 0.0 0.0\n3.0 0.0 0.0\n"
"0.0 1.5 0.0\n1.5 1.5 0.0\n3.0 1.5 0.0\n"
"0.0 3.0 0.0\n1.5 3.0 0.0\n3.0 3.0 0.0\n"
"0.0 0.0 1.5\n1.5 0.0 1.5\n3.0 0.0 1.5\n"
"0.0 1.5 1.5\n1.5 1.5 1.5\n3.0 1.5 1.5\n"
"0.0 3.0 1.5\n1.5 3.0 1.5\n3.0 3.0 1.5\n"
"0.0 0.0 3.0\n1.5 0.0 3.0\n3.0 0.0 3.0\n"
"0.0 1.5 3.0\n1.5 1.5 3.0\n3.0 1.5 3.0\n"
"0.0 3.0 3.0\n1.5 3.0 3.0\n3.0 3.0 3.0\n";

  // VTK file for 2x2x2 block of hexes as rectilinear-grid
extern const char rectilinear_grid_data[] =
"# vtk DataFile Version 2.0\n"
"Mesquite Mesh\n"
"ASCII\n"
"DATASET RECTILINEAR_GRID\n"
"DIMENSIONS 3 3 3\n"
"X_COORDINATES 3 float\n"
"0.0 1.5 3.0\n"
"Y_COORDINATES 3 float\n"
"0.0 1.5 3.0\n"
"Z_COORDINATES 3 float\n"
"0.0 1.5 3.0\n";

  // VTK file containing mixed-element unstructured mesh
  // First 8 elems result in same 2x2x2 block of hexes as
  // structured cases above.  The next 4 elems are the same
  // as the quads from the 2D structured-point file above.
extern const char mixed_unstructured_data[] = 
"# vtk DataFile Version 2.0\n"
"Mesquite Mesh\n"
"ASCII\n"
"DATASET UNSTRUCTURED_GRID\n"
"POINTS 35 float\n"
"\n" // points for an 2x2x2 brick of hexes (same geom/topo as above structured meshes)
"0.0 0.0 0.0\n1.5 0.0 0.0\n3.0 0.0 0.0\n"
"0.0 1.5 0.0\n1.5 1.5 0.0\n3.0 1.5 0.0\n"
"0.0 3.0 0.0\n1.5 3.0 0.0\n3.0 3.0 0.0\n"
"0.0 0.0 1.5\n1.5 0.0 1.5\n3.0 0.0 1.5\n"
"0.0 1.5 1.5\n1.5 1.5 1.5\n3.0 1.5 1.5\n"
"0.0 3.0 1.5\n1.5 3.0 1.5\n3.0 3.0 1.5\n"
"0.0 0.0 3.0\n1.5 0.0 3.0\n3.0 0.0 3.0\n"
"0.0 1.5 3.0\n1.5 1.5 3.0\n3.0 1.5 3.0\n"
"0.0 3.0 3.0\n1.5 3.0 3.0\n3.0 3.0 3.0\n"
"\n" // more points on +x side of brick for pyramids and tets
"4 0.75 0.75\n4 2.25 0.75\n4 0.75 2.25\n4 2.25 2.25\n"
"\n" // more points for two prisms/wedges
"6 0.75 0.75\n6 2.25 0.75\n6 0.75 2.25\n6 2.25 2.25\n"
"\n"
"CELLS 38 216\n"
"\n" // 8 hexes in 2x2x2 block
"8  0  1  4  3  9 10 13 12\n"
"8  1  2  5  4 10 11 14 13\n"
"8  3  4  7  6 12 13 16 15\n"
"8  4  5  8  7 13 14 17 16\n"
"8  9 10 13 12 18 19 22 21\n"
"8 10 11 14 13 19 20 23 22\n"
"8 12 13 16 15 21 22 25 24\n"
"8 13 14 17 16 22 23 26 25\n"
"\n"// Quads on -z face of hex (inverted to match structured data)
"4 0 1 4 3\n" 
"4 1 2 5 4\n"
"4 3 4 7 6\n"
"4 4 5 8 7\n"
"\n" // some pyramids on the +x side of the block
"5  2  5 14 11 27\n" 
"5  5  8 17 14 28\n"
"5 11 14 23 20 29\n"
"5 14 17 26 23 30\n"
"\n" // Some tetrahedrons around the pyramids
"4  5 14 27 28\n"
"4 14 28 17 30\n"
"4 29 14 23 30\n"
"4 11 27 14 29\n"
"4 27 29 30 14\n"
"4 28 27 30 14\n"
"\n" // Triangles bounding the pyramid/tet region
"3  2  5 27\n"
"3  5 28 27\n"
"3  5  8 28\n"
"3  8 17 28\n"
"3 17 30 28\n"
"3 17 26 30\n"
"3 26 23 30\n"
"3 23 29 30\n"
"3 23 20 29\n"
"3 20 11 29\n"
"3 11 27 29\n"
"3  2 27 11\n"
"3 27 28 30\n"
"3 27 30 29\n"
"\n" // A couple wedges/prisms
"6 27 30 28 31 34 32\n"
"6 30 27 29 34 31 33\n"
"\n" 
"CELL_TYPES 38\n"
"12 12 12 12 12 12 12 12\n"  //  8 hexes
" 9  9  9  9\n"              //  4 quads
"14 14 14 14\n"              //  4 pyramids
"10 10 10 10 10 10\n"        //  6 tets
" 5  5  5  5  5  5  5 "
" 5  5  5  5  5  5  5\n"     // 14 tri
"13 13\n";                   // 2 wedges

extern const char quadratic_unstructured_data[] = 
"# vtk DataFile Version 2.0\n"
"Mesquite Mesh\n"
"ASCII\n"
"DATASET UNSTRUCTURED_GRID\n"
"POINTS 41 float\n"
"\n" // points for a single quadtratic hex, in ExodusII order
" 1.0 -1.0 -1.0 \n" //  0  bottom corners
" 1.0  1.0 -1.0 \n" //  1
"-1.0  1.0 -1.0 \n" //  2
"-1.0 -1.0 -1.0 \n" //  3
" 1.0 -1.0  1.0 \n" //  4  top corners
" 1.0  1.0  1.0 \n" //  5
"-1.0  1.0  1.0 \n" //  6
"-1.0 -1.0  1.0 \n" //  7
" 1.0  0.0 -1.0 \n" //  8  bottom mid-nodes
" 0.0  1.0 -1.0 \n" //  9
"-1.0  0.0 -1.0 \n" // 10
" 0.0 -1.0 -1.0 \n" // 11
" 1.0 -1.0  0.0 \n" // 12  side mid-nodes 
" 1.0  1.0  0.0 \n" // 13
"-1.0  1.0  0.0 \n" // 14
"-1.0 -1.0  0.0 \n" // 15
" 1.0  0.0  1.0 \n" // 16  top mid-nodes
" 0.0  1.0  1.0 \n" // 17
"-1.0  0.0  1.0 \n" // 18
" 0.0 -1.0  1.0 \n" // 19
" 1.0  0.0  0.0 \n" // 20  mid-nodes for side faces
" 0.0  1.0  0.0 \n" // 21
"-1.0  0.0  0.0 \n" // 22
" 0.0 -1.0  0.0 \n" // 23
" 0.0  0.0 -1.0 \n" // 24  mid-nodes for top/bottom faces
" 0.0  0.0  1.0 \n" // 25
" 0.0  0.0  0.0 \n" // 26  mid-region
"\n" // points for a single quadtatic tet, in ExodusII order
" 1.0 -1.0 -1.0 \n" // 27  base triangle
" 1.0  1.0 -1.0 \n" // 28
"-1.0  0.0 -1.0 \n" // 29
" 0.0  0.0  1.0 \n" // 30  apex                            
" 1.0  0.0 -1.0 \n" // 31  base mid-nodes
" 0.0  0.5 -1.0 \n" // 32
" 0.0 -0.5 -1.0 \n" // 33
" 0.5 -0.5  0.0 \n" // 34  side mid-nodes
" 0.5  0.5  0.0 \n" // 35
"-0.5  0.0  0.0 \n" // 36
"\n" // mid-edge nodes for apex edges of quadratic pyramid
" 0.5 -0.5 0 \n" // 37
" 0.5  0.5 0 \n" // 38
"-0.5  0.5 0 \n" // 39
"-0.5 -0.5 0 \n" // 40
"\n"
"CELLS 8 116\n"
"20 0 1 2 3 4 5 6 7 8 9 10 11 16 17 18 19 12 13 14 15\n" // Hex20
"27 0 1 2 3 4 5 6 7 8 9 10 11 16 17 18 19 12 13 14 15 23 21 20 22 24 25 26\n" // Hex27
"10 27 28 29 30 31 32 33 34 35 36\n" // Tet10
"8 0 1 2 3 8 9 10 11\n" // Quad8 (base of hex)
"9 0 1 2 3 8 9 10 11 24\n" // Quad9 (base of hex)
"6 27 28 29 31 32 33\n" // Tri6 (base of tet)
"15 3 2 1 7 6 5 10 9 24 18 17 25 15 14 13\n" // quadratic prism as half of hex
"13 0 1 2 3 25 8 9 10 11 37 38 39 40\n" // quadratic pyramid with same base as hex
"\n"
"CELL_TYPES 8\n"
"25 29 24 23 28 22 26 27\n";


  // A simple scalar attribute specifying point and hex
  // IDs.  May be appended to any of the above 3D structured
  // mesh files.
extern const char simple_scalar_attrib[] = 
"CELL_DATA 8\n"
"SCALARS global_id int 1\n"
"LOOKUP_TABLE default\n"
"1 2 3 4 5 6 7 8\n"
"\n"
"POINT_DATA 27\n"
"SCALARS global_id int\n"
"LOOKUP_TABLE default\n"
" 1  2  3  4  5  6  7  8  9\n"
"10 11 12 13 14 15 16 17 18\n"
"19 20 21 22 23 24 25 26 27\n";

  // A VTK vector attribute.  May be appended to any of the 
  // above 3D structured mesh files.
extern const char simple_vector_attrib[] =
"CELL_DATA 8\n"
"VECTORS hexvect float\n"
"1 1 1\n"
"2 2 2\n"
"3 3 3\n"
"4 4 4\n"
"5 5 5\n"
"6 6 6\n"
"7 7 7\n"
"8 8 8\n";

extern const char simple_field_attrib[] = 
"CELL_DATA 8\n"
"FIELD test_field 2\n"
"elem_ids 1 8 int\n"
"1 2 3 4 5 6 7 8\n"
"elem_vects 3 8 float\n"
"1 1 1\n"
"2 2 2\n"
"3 3 3\n"
"4 4 4\n"
"5 5 5\n"
"6 6 6\n"
"7 7 7\n"
"8 8 8\n"
"FIELD field1 1\n"
"values 1 8 int\n"
"8 7 6 5 4 3 2 1\n";


  // A scalar VTK attribute with the name and datatype
  // expected by MeshImpl for specifying boundary vertices.
  // May be appended to any of the above 3D structured
  // mesh files
extern const char fixed_vertex_attrib[] = 
"POINT_DATA 27\n"
"SCALARS fixed float\n"
"LOOKUP_TABLE default\n"
"1 1 1 1 1 1 1 1 1\n"
"1 1 1 1 0 1 1 1 1\n"
"1 1 1 1 1 1 1 1 1\n";

#ifndef DEBUG

class VtkTest : public CppUnit::TestFixture
{
private:
   CPPUNIT_TEST_SUITE(VtkTest);

    // Original test for old Vtk parser
   CPPUNIT_TEST (test_elements);

    // Additional tests for new Vtk parser - J.Kraftcheck, 2004-10-12
   CPPUNIT_TEST (test_read_unstructured);
   CPPUNIT_TEST (test_read_structured_2d_points);
   CPPUNIT_TEST (test_read_structured_3d_points);
   CPPUNIT_TEST (test_read_structured_grid);
   CPPUNIT_TEST (test_read_rectilinear_grid);
   CPPUNIT_TEST (test_read_simple_scalar_attrib);
   CPPUNIT_TEST (test_read_vector_attrib);
   CPPUNIT_TEST (test_read_field_attrib);
   CPPUNIT_TEST (test_write_field_attrib);
   CPPUNIT_TEST (test_read_fixed_attrib);
   CPPUNIT_TEST (test_read_quadratic);
   CPPUNIT_TEST (test_write_quadratic);
   
    // Test writer - J.Kraftcheck, 2004-10-12
   CPPUNIT_TEST (test_write);

   CPPUNIT_TEST_SUITE_END();
   
public:
   /* Automatically called by CppUnit before each test function. */
  void setUp()
  {
  }
  
    // Automatically called by CppUnit after each test function.
  void tearDown()
  {
  }
  
public:
  VtkTest()
    {}
  
    // Check if the 2x2x2 brick of structured mesh
    // read from file is as expected.
  void check_8hex_structured( Mesh& mesh );
  
    // Check if the 2x2x2 brick of hexes
    // read from file is as expected.
  void check_8hex_block( Mesh& mesh, std::vector<Mesh::VertexHandle>::iterator connectivity );
  
    // Check if the 2x2 brick of structured mesh
    // read from file is as expected.
  void check_4quad_structured( Mesh& mesh );
  
    // Check if the 2x2 brick of quads
    // read from file is as expected.
  void check_4quad_block( Mesh& mesh, std::vector<Mesh::VertexHandle>::iterator connectivity );
        
    
    // Test reading VTK unstructured mesh
  void test_read_unstructured();
  
  void test_read_unstructured( const char* filename );
  
  
    // Test reading 2D Vtk structured-points mesh
  void test_read_structured_2d_points();
  
  
    // Test reading 3D Vtk structured-points mesh
  void test_read_structured_3d_points();
  
  
    // Test reading 3D Vtk structured-grid mesh
  void test_read_structured_grid();
  
  
    // Test reading 3D Vtk rectilinear-grid mesh
  void test_read_rectilinear_grid();
  
  
    // Test reading Vtk simple (one-component) scalar attribute
  void test_read_simple_scalar_attrib();
   
  
    // Test reading Vtk vector attribute
  void test_read_vector_attrib();
  

    // Test reading VTK FIELD attribute data
  void test_read_field_attrib();
  void test_write_field_attrib();
  void check_field_attrib( const char* file );

    // Test reading MeshImpl boundary-vertex bit
    // from Vtk scalar attribute.
  void test_read_fixed_attrib();
 
    
    // Test writing VTK unstructured mesh
  void test_write();
  
  void test_read_quadratic();
  void test_read_quadratic(const char* filename);
  void test_write_quadratic();
  
  void test_elements();
  
  int tri_check_validity(const Mesquite::MsqMeshEntity* element_array,
                         size_t num_elements,
                         const Mesquite::MsqVertex* vtx_array,
                         size_t num_vertices);
  
  int tet_validity_check(const Mesquite::MsqMeshEntity* element_array,
                         size_t num_elements,
                         const Mesquite::MsqVertex *vtx_array);
};
    // Check if the 2x2x2 brick of structured mesh
    // read from file is as expected.
  void VtkTest::check_8hex_structured( Mesh& mesh )
  {
    MsqPrintError err(cout);
    
    std::vector<Mesh::ElementHandle> elems(8);
    std::vector<Mesh::VertexHandle> verts(64);
    std::vector<size_t> offsets(9);
    
    mesh.get_all_elements( elems, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(elems.size(), (size_t)8);
    
    mesh.elements_get_attached_vertices( arrptr(elems), elems.size(),
                                         verts, offsets, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(verts.size() == 64);
    CPPUNIT_ASSERT(offsets.size() == 9);
    
    check_8hex_block( mesh, verts.begin() );
  }
  
    // Check if the 2x2x2 brick of hexes
    // read from file is as expected.
  void VtkTest::check_8hex_block( Mesh& mesh, std::vector<Mesh::VertexHandle>::iterator connectivity )
  {
    MsqPrintError err(cout);
    const int base_corners[] =   { 0, 1, 3, 4, 9, 10, 12, 13 };
    const int corner_offsets[] = { 0, 1, 4, 3, 9, 10, 13, 12 };

    for (int hex = 0; hex < 8; ++hex)
    {
      for (int node = 0; node < 8; ++node)
      {
        const int index = base_corners[hex] + corner_offsets[node];
        const int x = index % 3;
        const int y = (index / 3) % 3;
        const int z = index / 9;
        const Vector3D expected_coords( 1.5*x, 1.5*y, 1.5*z );
        MsqVertex actual_coords;
        Mesh::VertexHandle* conn_ptr = &*connectivity;
        ++connectivity;
        mesh.vertices_get_coordinates( conn_ptr, &actual_coords, 1, err );
        CPPUNIT_ASSERT( !err );
        CPPUNIT_ASSERT( expected_coords.within_tolerance_box( actual_coords, DBL_EPSILON ) );
      }
    }
  }
  
    // Check if the 2x2 brick of structured mesh
    // read from file is as expected.
  void VtkTest::check_4quad_structured( Mesh& mesh )
  {
    MsqPrintError err(cout);
    
    std::vector<Mesh::ElementHandle> elems(4);
    std::vector<Mesh::VertexHandle> verts(9);
    std::vector<size_t> offsets(5);
    
    mesh.get_all_elements( elems, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(elems.size(), (size_t)4);
    
    mesh.elements_get_attached_vertices( arrptr(elems), elems.size(),
                                         verts, offsets, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(verts.size() == 16);
    CPPUNIT_ASSERT(offsets.size() == 5);
    
    check_4quad_block( mesh, verts.begin() );
  }
  
    // Check if the 2x2 brick of quads
    // read from file is as expected.
  void VtkTest::check_4quad_block( Mesh& mesh, std::vector<Mesh::VertexHandle>::iterator connectivity )
  {
    MsqPrintError err(cout);
    const int base_corners[]   = { 0, 1, 3, 4 };
    const int corner_offsets[] = { 0, 1, 4, 3 };
    for (int quad = 0; quad < 4; ++quad)
    {
      for (int node = 0; node < 4; ++node)
      {
        const int index = base_corners[quad] + corner_offsets[node];
        const int x = index % 3;
        const int y = index / 3;
        const Vector3D expected_coords( 1.5*x, 1.5*y, 0.0 );
        MsqVertex actual_coords;
        Mesh::VertexHandle* conn_ptr = &*connectivity;
        ++connectivity;
        mesh.vertices_get_coordinates( conn_ptr, &actual_coords, 1, err );
        CPPUNIT_ASSERT( !err );
        CPPUNIT_ASSERT( expected_coords.within_tolerance_box( actual_coords, DBL_EPSILON ) );
      }
    }
  }
        
    
    // Test reading VTK unstructured mesh
  void VtkTest::test_read_unstructured()
  {
    FILE* file = fopen( temp_file_name, "w+" );
    CPPUNIT_ASSERT( file );
    int rval = fputs( mixed_unstructured_data, file );
    fclose( file );
    if (rval == EOF) remove( temp_file_name );
    CPPUNIT_ASSERT( rval != EOF );

    test_read_unstructured( temp_file_name );
  }
  
  void VtkTest::test_read_unstructured( const char* filename )
  {
    MeshImpl mesh;
    MsqPrintError err(cout);
    
    mesh.read_vtk( filename, err );
    ASSERT_NO_ERROR(err);
    
      // Get mesh data
    std::vector<Mesh::VertexHandle> conn;
    std::vector<Mesh::ElementHandle> elems(38);
    std::vector<size_t> offsets(39);
    mesh.get_all_elements( elems, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_EQUAL( elems.size(), (size_t)38 );
    mesh.elements_get_attached_vertices( arrptr(elems), elems.size(),
                                         conn, offsets, err );
    ASSERT_NO_ERROR(err);

    unsigned i;
    struct meshdata { EntityTopology type; size_t nodes; size_t count; };
    meshdata list[] = {
      { Mesquite::HEXAHEDRON,    8,  8 },
      { Mesquite::QUADRILATERAL, 4,  4 },
      { Mesquite::PYRAMID,       5,  4 },
      { Mesquite::TETRAHEDRON,   4,  6 },
      { Mesquite::TRIANGLE,      3, 14 },
      { Mesquite::PRISM,         6,  2 }, 
      { Mesquite::MIXED,         0,  0 } };
      
      // Count expected lenght of connectivity list
    size_t conn_len = 0;
    for (i = 0; list[i].nodes; ++i)
      conn_len += list[i].nodes * list[i].count;
    CPPUNIT_ASSERT_EQUAL( conn_len, conn.size() );
    
    check_8hex_block( mesh, conn.begin() );
    check_4quad_block( mesh, conn.begin() + 64 );
  }
  
  
    // Test reading 2D Vtk structured-points mesh
  void VtkTest::test_read_structured_2d_points()
  {
    MeshImpl mesh;
    MsqPrintError err(cout);
   
    FILE* file = fopen( temp_file_name, "w+" );
    fputs( structured_2d_points_data, file );
    fclose( file );
    
    mesh.read_vtk( temp_file_name, err );
    remove( temp_file_name );
    ASSERT_NO_ERROR(err);
    
    check_4quad_structured( mesh );
  }
  
  
    // Test reading 3D Vtk structured-points mesh
  void VtkTest::test_read_structured_3d_points()
  {
    MeshImpl mesh;
    MsqPrintError err(cout);
    
    FILE* file = fopen( temp_file_name, "w+" );
    fputs( structured_3d_points_data, file );
    fclose( file );
    
    mesh.read_vtk( temp_file_name, err );
    remove( temp_file_name );
    ASSERT_NO_ERROR(err);
    
    check_8hex_structured( mesh );
  }
  
  
    // Test reading 3D Vtk structured-grid mesh
  void VtkTest::test_read_structured_grid()
  {
    MeshImpl mesh;
    MsqPrintError err(cout);
    
    FILE* file = fopen( temp_file_name, "w+" );
    fputs( structured_grid_data, file );
    fclose( file );
    
    mesh.read_vtk( temp_file_name, err );
    remove( temp_file_name );
    ASSERT_NO_ERROR(err);
    
    check_8hex_structured( mesh );
  }
  
  
    // Test reading 3D Vtk rectilinear-grid mesh
  void VtkTest::test_read_rectilinear_grid()
  {
    MeshImpl mesh;
    MsqPrintError err(cout);
    
    FILE* file = fopen( temp_file_name, "w+" );
    fputs( rectilinear_grid_data, file );
    fclose( file );
    
    mesh.read_vtk( temp_file_name, err );
    remove( temp_file_name );
    ASSERT_NO_ERROR(err);
    
    check_8hex_structured( mesh );
  }
  
  
    // Test reading Vtk simple (one-component) scalar attribute
  void VtkTest::test_read_simple_scalar_attrib()
  {
    MeshImpl mesh;
    MsqPrintError err(cout);
   
    FILE* file = fopen( temp_file_name, "w+" );
    fputs( structured_3d_points_data, file );
    fputs( simple_scalar_attrib, file );
    fclose( file );
    
    mesh.read_vtk( temp_file_name, err );
    remove( temp_file_name );
    ASSERT_NO_ERROR(err);
    
    std::vector<Mesh::ElementHandle> elems;
    mesh.get_all_elements( elems, err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( elems.size(), (size_t)8 );
   
    void* th = mesh.tag_get( "global_id", err );
    CPPUNIT_ASSERT( !err );

    std::string name;
    Mesh::TagType type;
    unsigned tagsize;
    mesh.tag_properties( th, name, type, tagsize, err );
    CPPUNIT_ASSERT( !err && type == Mesh::INT && tagsize == 1 );
    
    int elem_data[8];
    mesh.tag_get_element_data( th, 8, arrptr(elems), elem_data, err );
    CPPUNIT_ASSERT( !err );
    
    for (int i = 0; i < 8; ++i)
      CPPUNIT_ASSERT( elem_data[i] == (1+i) );
  }
   
  
    // Test reading Vtk vector attribute
  void VtkTest::test_read_vector_attrib()
  {
    MeshImpl mesh;
    MsqPrintError err(cout);
    
    FILE* file = fopen( temp_file_name, "w+" );
    fputs( structured_3d_points_data, file );
    fputs( simple_vector_attrib, file );
    fclose( file );
    
    mesh.read_vtk( temp_file_name, err );
    remove( temp_file_name );
    ASSERT_NO_ERROR(err);
    
    std::vector<Mesh::ElementHandle> elems;
    mesh.get_all_elements( elems, err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( elems.size(), (size_t)8 );
   
    void* th = mesh.tag_get( "hexvect", err );
    CPPUNIT_ASSERT( !err );

    std::string name;
    Mesh::TagType type;
    unsigned tagsize;
    mesh.tag_properties( th, name, type, tagsize, err );
    CPPUNIT_ASSERT( !err && type == Mesh::DOUBLE && tagsize == 3 );
    
    double elem_data[24];
    mesh.tag_get_element_data( th, 8, arrptr(elems), elem_data, err );
    CPPUNIT_ASSERT( !err );
    
    for (int i = 0; i < 8; ++i)
      CPPUNIT_ASSERT( Vector3D( elem_data+3*i ) == Vector3D( i+1, i+1, i+1 ) );
  }
  
  
    // Test reading Vtk field attribute
  void VtkTest::test_read_field_attrib()
  {
    FILE* file = fopen( temp_file_name, "w+" );
    fputs( structured_3d_points_data, file );
    fputs( simple_field_attrib, file );
    fclose( file );
    check_field_attrib( temp_file_name );
  }
  
  void VtkTest::check_field_attrib( const char* temp_file_name )
  {
    MeshImpl mesh;
    MsqPrintError err(cout);
    
    mesh.read_vtk( temp_file_name, err );
    remove( temp_file_name );
    ASSERT_NO_ERROR(err);
    
    std::vector<Mesh::ElementHandle> elems;
    mesh.get_all_elements( elems, err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( elems.size(), (size_t)8 );
   
    std::string name;
    Mesh::TagType type;
    unsigned tagsize;


    void* th = mesh.tag_get( "test_field elem_vects", err );
    CPPUNIT_ASSERT( !err );

    mesh.tag_properties( th, name, type, tagsize, err );
    CPPUNIT_ASSERT( !err && type == Mesh::DOUBLE && tagsize == 3 );
    
    double elem_data[24];
    mesh.tag_get_element_data( th, 8, arrptr(elems), elem_data, err );
    CPPUNIT_ASSERT( !err );
    
    for (int i = 0; i < 8; ++i)
      CPPUNIT_ASSERT( Vector3D( elem_data+3*i ) == Vector3D( i+1, i+1, i+1 ) );
   
    
    th = mesh.tag_get( "test_field elem_ids", err );
    CPPUNIT_ASSERT( !err );

    mesh.tag_properties( th, name, type, tagsize, err );
    CPPUNIT_ASSERT( !err && type == Mesh::INT && tagsize == 1 );
    
    int elem_ids[8];
    mesh.tag_get_element_data( th, 8, arrptr(elems), elem_ids, err );
    CPPUNIT_ASSERT( !err );
    
    for (int i = 0; i < 8; ++i)
      CPPUNIT_ASSERT( elem_ids[i] == i+1 );
   
    
    th = mesh.tag_get( "field1", err );
    CPPUNIT_ASSERT( !err );

    mesh.tag_properties( th, name, type, tagsize, err );
    CPPUNIT_ASSERT( !err && type == Mesh::INT && tagsize == 1 );
    
    int values[8];
    mesh.tag_get_element_data( th, 8, arrptr(elems), values, err );
    CPPUNIT_ASSERT( !err );
    
    for (int i = 0; i < 8; ++i)
      CPPUNIT_ASSERT( values[i] == 8-i );
  }

    
    // Test writing quadtratic elements
  void VtkTest::test_write_field_attrib()
  {
    MeshImpl mesh;
    MsqPrintError err(cout);
    
      // Create file containing unstructured mesh test case
    FILE* file = fopen( temp_file_name, "w+" );
    fputs( structured_3d_points_data, file );
    fputs( simple_field_attrib, file );
    fclose( file );
    
      // Read unstructured mesh file
    mesh.read_vtk( temp_file_name, err );
    remove(temp_file_name);
    ASSERT_NO_ERROR(err);
    
      // Write unstructured mesh file back out
    mesh.write_vtk( temp_file_name, err );
    if (err)  remove( temp_file_name );
    ASSERT_NO_ERROR(err);
    
      // Check if file contained expected mesh
    check_field_attrib( temp_file_name );
  }

    // Test reading MeshImpl boundary-vertex bit
    // from Vtk scalar attribute.
  void VtkTest::test_read_fixed_attrib()
  {
    MeshImpl mesh;
    MsqPrintError err(cout);
    
    FILE* file = fopen( temp_file_name, "w+" );
    fputs( structured_3d_points_data, file );
    fputs( fixed_vertex_attrib, file );
    fclose( file );
    
    mesh.read_vtk( temp_file_name, err );
    remove( temp_file_name );
    ASSERT_NO_ERROR(err);

    std::vector<Mesh::ElementHandle> elems;
    mesh.get_all_elements( elems, err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( elems.size(), (size_t)8 );
    
    std::vector<Mesh::VertexHandle> verts;
    std::vector<size_t> offsets;
    mesh.elements_get_attached_vertices( arrptr(elems), elems.size(), verts, offsets, err );
    ASSERT_NO_ERROR(err);
    
    // get unique list of vertices
    std::vector<Mesh::VertexHandle>::iterator new_end;
    std::sort( verts.begin(), verts.end() );
    new_end = std::unique( verts.begin(), verts.end() );
    verts.resize( new_end - verts.begin() );
    CPPUNIT_ASSERT_EQUAL( verts.size(), (size_t)27 );

    // get fixed flag
    std::vector<bool> fixed;
    mesh.vertices_get_fixed_flag( arrptr(verts), fixed, verts.size(), err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_EQUAL( verts.size(), fixed.size() );
    
    for (int i = 0; i < 27; ++i)
    {
      bool should_be_fixed = (i != 13);
      CPPUNIT_ASSERT_EQUAL( should_be_fixed, (bool)fixed[i] );
    }
  }
 
    
    // Test writing VTK unstructured mesh
  void VtkTest::test_write()
  {
    MeshImpl mesh;
    MsqPrintError err(cout);
    
      // Create file containing unstructured mesh test case
    FILE* file = fopen( temp_file_name, "w+" );
    CPPUNIT_ASSERT(file);
    int rval = fputs( mixed_unstructured_data, file );
    fclose( file );
    if (rval == EOF) remove(temp_file_name);
    CPPUNIT_ASSERT(rval != EOF);
    
      // Read unstructured mesh file
    mesh.read_vtk( temp_file_name, err );
    remove(temp_file_name);
    ASSERT_NO_ERROR(err);
    
      // Write unstructured mesh file back out
    mesh.write_vtk( temp_file_name, err );
    if (err)  remove( temp_file_name );
    ASSERT_NO_ERROR(err);
    
      // Check if file contained expected mesh
    test_read_unstructured( temp_file_name );
    remove( temp_file_name );
  }
  
    
    // Test reading VTK unstructured mesh quadratic elements
  void VtkTest::test_read_quadratic()
  {
    FILE* file = fopen( temp_file_name, "w+" );
    CPPUNIT_ASSERT( file );
    int rval = fputs( quadratic_unstructured_data, file );
    fclose( file );
    if (rval == EOF) remove( temp_file_name );
    CPPUNIT_ASSERT( rval != EOF );

    test_read_quadratic( temp_file_name );
    remove( temp_file_name );
  }
  
  void VtkTest::test_read_quadratic( const char* filename )
  {
    const size_t NUM_ELEM = 8;
  
    MeshImpl mesh;
    MsqPrintError err(cout);
    
    mesh.read_vtk( filename, err );
    ASSERT_NO_ERROR(err);
    
    std::vector<Mesh::ElementHandle> elems(NUM_ELEM);
    mesh.get_all_elements( elems, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_EQUAL(elems.size(), NUM_ELEM );
    
    std::vector<Mesh::VertexHandle> conn;
    std::vector<size_t> offsets;
    mesh.elements_get_attached_vertices( arrptr(elems), elems.size(), conn, offsets, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT_EQUAL( conn.size(), (size_t)108 );
    
    EntityTopology types[NUM_ELEM];
    mesh.elements_get_topologies( arrptr(elems), types, NUM_ELEM, err );
    ASSERT_NO_ERROR(err);

    static const double hex_corners[] = 
     {  1.0, -1.0, -1.0, 
        1.0,  1.0, -1.0, 
       -1.0,  1.0, -1.0, 
       -1.0, -1.0, -1.0,
        1.0, -1.0,  1.0, 
        1.0,  1.0,  1.0, 
       -1.0,  1.0,  1.0, 
       -1.0, -1.0,  1.0 };
    static const double tet_corners[] = 
     {  1.0, -1.0, -1.0,
        1.0,  1.0, -1.0,
       -1.0,  0.0, -1.0,
        0.0,  0.0,  1.0 };
    static const double pyr_corners[] = 
     {  1.0, -1.0, -1.0, 
        1.0,  1.0, -1.0, 
       -1.0,  1.0, -1.0, 
       -1.0, -1.0, -1.0,
        0.0,  0.0,  1.0 };
    static const double pri_corners[] = 
      { -1.0, -1.0, -1.0,
         1.0,  1.0, -1.0,
        -1.0,  1.0, -1.0,
        -1.0, -1.0,  1.0,
         1.0,  1.0,  1.0,
        -1.0,  1.0,  1.0 };
    static const unsigned hex_edges[] =
     { 0, 1, 
       1, 2, 
       2, 3,
       3, 0,
       0, 4,
       1, 5,
       2, 6,
       3, 7,
       4, 5,
       5, 6,
       6, 7,
       7, 4 };
    static const unsigned tet_edges[] = 
     { 0, 1,
       1, 2, 
       2, 0,
       0, 3,
       1, 3, 
       2, 3 };
    static const unsigned pri_edges[] =
     { 0, 1, 
       1, 2, 
       2, 0,
       0, 3,
       1, 4,
       2, 5,
       3, 4,
       4, 5,
       5, 3 };
    static const unsigned pyr_edges[] =
     { 0, 1, 
       1, 2, 
       2, 3,
       3, 0,
       0, 4,
       1, 4,
       2, 4,
       3, 4 };
    static const unsigned hex_faces[] = 
    { 4, 0, 1, 5, 4,
      4, 1, 2, 6, 5,
      4, 2, 3, 7, 6,
      4, 3, 0, 4, 7,
      4, 3, 2, 1, 0,
      4, 4, 5, 6, 7
    };
    static const struct {
      EntityTopology topology;
      unsigned num_corners;
      unsigned num_edges;
      unsigned num_faces; // if non-zero expect mid-face nodes
      unsigned num_region; // if non-zero expect mid-region node
      const double* corners;
      const unsigned* edges;
      const unsigned* faces;
    } expected_elems[NUM_ELEM] = {
      { Mesquite::HEXAHEDRON,    8, 12, 0, 0, hex_corners, hex_edges, hex_faces },
      { Mesquite::HEXAHEDRON,    8, 12, 6, 1, hex_corners, hex_edges, hex_faces },
      { Mesquite::TETRAHEDRON,   4,  6, 0, 0, tet_corners, tet_edges, 0 },
      { Mesquite::QUADRILATERAL, 4,  4, 0, 0, hex_corners, hex_edges, 0 },
      { Mesquite::QUADRILATERAL, 4,  4, 0, 1, hex_corners, hex_edges, 0 },
      { Mesquite::TRIANGLE,      3,  3, 0, 0, tet_corners, tet_edges, 0 },
      { Mesquite::PRISM,         6,  9, 0, 0, pri_corners, pri_edges, 0 },
      { Mesquite::PYRAMID,       5,  8, 0, 0, pyr_corners, pyr_edges, 0 } };
    
    MsqVertex have;
    std::vector<Mesh::VertexHandle>::iterator v_it = conn.begin();
    for (unsigned i = 0; i < NUM_ELEM; ++i)
    {
      CPPUNIT_ASSERT_EQUAL( expected_elems[i].topology, types[i] );

      size_t vtx_start = offsets[i];
      size_t vtx_end = offsets[i+1];
      size_t conn_len = expected_elems[i].num_corners 
                      + expected_elems[i].num_edges
                      + expected_elems[i].num_faces
                      + expected_elems[i].num_region;
      CPPUNIT_ASSERT_EQUAL( conn_len, vtx_end - vtx_start );
      
      for (unsigned c = 0; c < expected_elems[i].num_corners; ++c, ++v_it)
      {
        Vector3D expected(expected_elems[i].corners + 3*c);
        mesh.vertices_get_coordinates( &*v_it, &have, 1, err );
        ASSERT_NO_ERROR(err);
        expected -= have;
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, expected.length(), DBL_EPSILON );
      }
      
      for (unsigned m = 0; m < expected_elems[i].num_edges; ++m, ++v_it)
      {
        unsigned start_idx = expected_elems[i].edges[2*m];
        unsigned end_idx = expected_elems[i].edges[2*m+1];
        Vector3D start( expected_elems[i].corners + 3*start_idx );
        Vector3D end( expected_elems[i].corners + 3*end_idx );
        Vector3D expected = 0.5 * (start + end);
        
        mesh.vertices_get_coordinates( &*v_it, &have, 1, err );
        ASSERT_NO_ERROR(err);
        
        expected -= have;
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, expected.length(), DBL_EPSILON );
      }

      const unsigned* f_it = expected_elems[i].faces;
      for (unsigned m = 0; m < expected_elems[i].num_faces; ++m, ++v_it)
      {
        Vector3D expected(0,0,0);
        const unsigned face_size = *f_it; ++f_it;
        CPPUNIT_ASSERT( face_size == 3u || face_size == 4u );
        for (unsigned f = 0; f < face_size; ++f, ++f_it) 
          expected += Vector3D( expected_elems[i].corners + 3 * *f_it );
        expected /= face_size;
        
        mesh.vertices_get_coordinates( &*v_it, &have, 1, err );
        ASSERT_NO_ERROR(err);
        
        expected -= have;
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, expected.length(), DBL_EPSILON );
      }
      
      if (expected_elems[i].num_region) {
        CPPUNIT_ASSERT_EQUAL( 1u, expected_elems[i].num_region );

        Vector3D expected(0,0,0);
        for (unsigned m = 0; m < expected_elems[i].num_corners; ++m)
          expected += Vector3D( expected_elems[i].corners + 3*m );
        expected /= expected_elems[i].num_corners;

        mesh.vertices_get_coordinates( &*v_it, &have, 1, err );
        ASSERT_NO_ERROR(err);

        expected -= have;
        CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, expected.length(), DBL_EPSILON );

        ++v_it;
      }
    }
  }
    
    // Test writing quadtratic elements
  void VtkTest::test_write_quadratic()
  {
    MeshImpl mesh;
    MsqPrintError err(cout);
    
      // Create file containing unstructured mesh test case
    FILE* file = fopen( temp_file_name, "w+" );
    CPPUNIT_ASSERT(file);
    int rval = fputs( quadratic_unstructured_data, file );
    fclose( file );
    if (rval == EOF) remove(temp_file_name);
    CPPUNIT_ASSERT(rval != EOF);
    
      // Read unstructured mesh file
    mesh.read_vtk( temp_file_name, err );
    remove(temp_file_name);
    ASSERT_NO_ERROR(err);
    
      // Write unstructured mesh file back out
    mesh.write_vtk( temp_file_name, err );
    if (err)  remove( temp_file_name );
    ASSERT_NO_ERROR(err);
    
      // Check if file contained expected mesh
    test_read_quadratic( temp_file_name );
    remove( temp_file_name );
  }

  void VtkTest::test_elements()
  {
    Mesquite::MsqPrintError err(cout);
    MeshImpl mMesh;
    mMesh.read_vtk(MESH_FILES_DIR "2D/vtk/tris/untangled/equil_tri2.vtk", err);
    ASSERT_NO_ERROR(err);
    Mesquite::MeshDomainAssoc mesh_and_domain = Mesquite::MeshDomainAssoc(&mMesh, 0);
    Mesquite::Instruction::initialize_vertex_byte( &mesh_and_domain, 0, err );
    ASSERT_NO_ERROR(err);
    
      // Retrieve a patch
    Mesquite::PatchData pd;
    pd.set_mesh( &mMesh );
    VertexPatches patch_set;
    patch_set.set_mesh( &mMesh );
    PatchIterator patches( &patch_set );
    patches.get_next_patch( pd, err );
    ASSERT_NO_ERROR(err);
    
    int free_vtx = pd.num_free_vertices(); 
//    std::cout << "nb of free vertices: " << free_vtx << std::endl;
    CPPUNIT_ASSERT( free_vtx == 1 );
    
    Mesquite::MsqMeshEntity* element_array =  pd.get_element_array(err); ASSERT_NO_ERROR(err);
    size_t num_elements = pd.num_elements();
    CPPUNIT_ASSERT( num_elements == 6 );
    
    const Mesquite::MsqVertex* vtx_array = pd.get_vertex_array(err); ASSERT_NO_ERROR(err);
    size_t num_vertices = pd.num_nodes();
    CPPUNIT_ASSERT( num_vertices == 7 );
    
    CPPUNIT_ASSERT( tri_check_validity(element_array, num_elements, vtx_array, num_vertices) == 1 );
    
    patches.get_next_patch( pd, err ); ASSERT_NO_ERROR(err);
    
    element_array =  pd.get_element_array(err); ASSERT_NO_ERROR(err);
    num_elements = pd.num_elements();
    CPPUNIT_ASSERT( num_elements == 6 );
    
    vtx_array = pd.get_vertex_array(err); ASSERT_NO_ERROR(err);
    num_vertices = pd.num_nodes();
    CPPUNIT_ASSERT( num_vertices == 7 );
    
    CPPUNIT_ASSERT( tri_check_validity(element_array, num_elements, vtx_array, num_vertices) == 1 );
  }
  
  int VtkTest::tri_check_validity(const Mesquite::MsqMeshEntity* element_array,
                         size_t num_elements,
                         const Mesquite::MsqVertex* vtx_array,
                         size_t num_vertices)
  {
       
      /* check that the simplicial mesh is still valid, 
         based on right handedness. Returns a 1 or a 0 */
    int valid = 1;
    double dEps = 1.e-13;
    
    double x1, x2, x3, y1, y2, y3;// z1, z2, z3;
    std::vector<size_t> vertex_indices;
    
    for (size_t i=0;i<num_elements;i++)
    {
      element_array[i].get_vertex_indices(vertex_indices);
      
      x1 = vtx_array[vertex_indices[0]][0];
      y1 = vtx_array[vertex_indices[0]][1];
      x2 = vtx_array[vertex_indices[1]][0];
      y2 = vtx_array[vertex_indices[1]][1];
      x3 = vtx_array[vertex_indices[2]][0];
      y3 = vtx_array[vertex_indices[2]][1];
      
      double a = x2*y3 - x3*y2;
      double b = y2 - y3;
      double c = x3 - x2;
      
      double area = .5*(a+b*x1+c*y1);
      if (area < dEps) {
          //          printf("x1 y1 = %f %f\n",x1,y1);
          //          printf("x2 y3 = %f %f\n",x2,y2);
          //          printf("x3 y3 = %f %f\n",x3,y3);
          //          printf("area = %f\n",area);
        valid=0;
      }
    }
    
    return(valid);
  }
  
  int VtkTest::tet_validity_check(const Mesquite::MsqMeshEntity* element_array,
                         size_t num_elements,
                         const Mesquite::MsqVertex *vtx_array)
  {
    int valid = 1;
    double dEps = 1.e-13;
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
    std::vector<size_t> vertex_indices;
    
    for (size_t i=0;i<num_elements;i++)
    {
      element_array[i].get_vertex_indices(vertex_indices);
      
      x1=vtx_array[vertex_indices[0]][0];
      y1=vtx_array[vertex_indices[0]][1];
      z1=vtx_array[vertex_indices[0]][2];
      
      x2=vtx_array[vertex_indices[1]][0];
      y2=vtx_array[vertex_indices[1]][1];
      z2=vtx_array[vertex_indices[1]][2];
      
      x3=vtx_array[vertex_indices[2]][0];
      y3=vtx_array[vertex_indices[2]][1];
      z3=vtx_array[vertex_indices[2]][2];
      
      x4=vtx_array[vertex_indices[3]][0];
      y4=vtx_array[vertex_indices[3]][1];
      z4=vtx_array[vertex_indices[3]][2];
      
      double dDX2 = x2 - x1;
      double dDX3 = x3 - x1;
      double dDX4 = x4 - x1;
      
      double dDY2 = y2 - y1;
      double dDY3 = y3 - y1;
      double dDY4 = y4 - y1;
      
      double dDZ2 = z2 - z1;
      double dDZ3 = z3 - z1;
      double dDZ4 = z4 - z1;
      
        /* dDet is proportional to the cell volume */
      double dDet = dDX2*dDY3*dDZ4 + dDX3*dDY4*dDZ2 + dDX4*dDY2*dDZ3
        - dDZ2*dDY3*dDX4 - dDZ3*dDY4*dDX2 - dDZ4*dDY2*dDX3 ;
      
        /* Compute a length scale based on edge lengths. */
      double dScale = ( sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) +
                             (z1-z2)*(z1-z2)) +
                        sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3) +
                             (z1-z3)*(z1-z3)) +
                        sqrt((x1-x4)*(x1-x4) + (y1-y4)*(y1-y4) +
                             (z1-z4)*(z1-z4)) +
                        sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) +
                             (z2-z3)*(z2-z3)) +
                        sqrt((x2-x4)*(x2-x4) + (y2-y4)*(y2-y4) +
                             (z2-z4)*(z2-z4)) +
                        sqrt((x3-x4)*(x3-x4) + (y3-y4)*(y3-y4) +
                             (z3-z4)*(z3-z4)) ) / 6.;
      
        /* Use the length scale to get a better idea if the tet is flat or
           just really small. */
      if (fabs(dScale) < dEps)
      {
        return(valid = 0);
      }
      else
      {
        dDet /= (dScale*dScale*dScale);
      }
      
      if (dDet > dEps)
      {
        valid = 1;
      }
      else if (dDet < -dEps)
      {
        valid = -1;
      }
      else
      {
        valid = 0;
      }
    }  // end for i=1,numElements
    
    return(valid);
  }


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(VtkTest, "VtkTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(VtkTest, "Unit");

#else /* ifndef DEBUG */

int  main()
{
  FILE* file;
  
  file = fopen( "2d_structured_points.vtk", "w" );
  fputs( structured_2d_points_data, file );
  fclose( file );
  
  file = fopen( "3d_structured_points.vtk", "w" );
  fputs( structured_3d_points_data, file );
  fputs( fixed_vertex_attrib, file );
  fclose( file );
  
  file = fopen( "structured_grid.vtk", "w" );
  fputs( structured_grid_data, file );
  fputs( fixed_vertex_attrib, file );
  fclose( file );
  
  file = fopen( "rectilinear_grid.vtk", "w" );
  fputs( rectilinear_grid_data, file );
  fputs( fixed_vertex_attrib, file );
  fclose( file );
  
  file = fopen( "mixed_unstructured.vtk", "w" );
  fputs( mixed_unstructured_data, file );
  fclose( file );
  
  file = fopen( "quadratic.vtk", "w" );
  fputs( quadratic_unstructured_data, file );
  fclose( file );
  
  return 0;
}
#endif
