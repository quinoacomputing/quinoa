/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
 
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */

#include <iostream>
using std::cout;
using std::endl;
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
using std::set;
using std::map;
using std::vector;
using std::string;

#include <stdlib.h>


#include "Mesquite.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_Vector3D.hpp"
#include "Mesquite_MsqVertex.hpp"
#include "Mesquite_MsqIMesh.hpp"
#include "UnitUtil.hpp"

using namespace Mesquite;

/* This is the mesh defined in the two arrays below.  This test
 * should work with any mesh composed entirely of triangle
 * elements.  The mesh used for the test can be changed by modifying
 * the two arrays below.
 *                                             
 *            3-------------------2            
 *           / \                 / \           
 *          /   \               /   \          
 *         /     \             /     \         
 *        /       \    (1)    /       \        
 *       /         \         /         \       
 *      /    (2)    \       /    (0)    \      
 *     /             \     /             \     
 *    /               \   /               \    
 *   /                 \ /                 \   
 *  4-------------------0-------------------1  
 *   \                 / \                 /   
 *    \               /   \               /    
 *     \             /     \             /     
 *      \    (3)    /       \    (5)    /      
 *       \         /         \         /       
 *        \       /    (4)    \       /        
 *         \     /             \     /         
 *          \   /               \   /          
 *           \ /                 \ /           
 *            5-------------------6            
 *                                             
 */

extern const double vertexCoords[] = {
   0.00,  0.00,  0.00,  // 0
   2.00,  0.00,  0.00,  // 1
   1.00,  1.73,  0.00,  // 2
  -1.00,  1.73,  0.00,  // 3
  -2.00,  0.00,  0.00,  // 4
  -1.00, -1.73,  0.00,  // 5
   1.00, -1.73,  0.00   // 6
};

extern const int triangleConnectivity[] = {
 0, 1, 2, // 0
 0, 2, 3, // 1
 0, 3, 4, // 2
 0, 4, 5, // 3
 0, 5, 6, // 4
 0, 6, 1  // 5
};



class iMeshTest : public CppUnit::TestFixture
{
  private:
    
    CPPUNIT_TEST_SUITE( iMeshTest );
//    CPPUNIT_TEST( testVertexIterator );
    CPPUNIT_TEST( testVertexByte );
    CPPUNIT_TEST( testFixedFlag );
    CPPUNIT_TEST( testSlavedFlag );
    CPPUNIT_TEST( testVertexAdjacency );
    CPPUNIT_TEST( testElementConnectivity );
    CPPUNIT_TEST( testElementTopology );
    CPPUNIT_TEST( testIntTag );
    CPPUNIT_TEST( testDoubleTag );
    CPPUNIT_TEST_SUITE_END();

    MsqIMesh* myMesh;
    iMesh_Instance myIMesh;
    
    Mesh::VertexHandle vtxIndexToHandle[7];
    Mesh::ElementHandle triIndexToHandle[7];
    map<Mesh::VertexHandle,int> vtxHandleToIndex;
    map<Mesh::ElementHandle,int> triHandleToIndex;

   bool match_triangles( const int* tri1, const Mesh::VertexHandle* tri2_handles );
   bool writeVtkFile( const char* name );

  public:
  
    iMeshTest()
     : myMesh(0)
      {}
  
    void setUp();
    void tearDown();
    
    void testVertexFlag(bool fixed, iBase_TagValueType type);
    void testVertexFlagNone( bool fixed );
    
    void matchVertexCoordinates();
    void matchElementConnectivity();
    void testVertexIterator();
    void testVertexByte();
    void testFixedFlag() { testVertexFlag(true,iBase_INTEGER); }
    void testSlavedFlag() { testVertexFlag(false,iBase_INTEGER); }
    void testFixedFlagByte() { testVertexFlag(true,iBase_BYTES); }
    void testSlavedFlagByte() { testVertexFlag(false,iBase_BYTES); }
    void testFixedFlagNone() { testVertexFlagNone(true); }
    void testSlavedFlagNone() { testVertexFlagNone(false); }
    void testVertexAdjacency();
    void testElementConnectivity();
    void testElementTopology();
    void testIntTag();
    void testDoubleTag();
    
};

bool iMeshTest::writeVtkFile( const char* filename )
{
  int i;
  
  FILE* file = fopen( filename, "w" );
  if (!file) {
    perror( filename );
    return false;
  }
  
  fputs( "# vtk DataFile Version 2.0\n"
         "Mesquite Mesh\n"
         "ASCII\n"
         "DATASET UNSTRUCTURED_GRID\n",
         file );
  
  const int num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  fprintf( file, "POINTS %d float\n", num_pts );
  for (i = 0; i < num_pts; ++i)
    fprintf( file, "%f %f %f\n",
                   vertexCoords[3*i    ],
                   vertexCoords[3*i + 1],
                   vertexCoords[3*i + 2] );
  
  int num_tris = sizeof(triangleConnectivity) / (3*sizeof(int));
  fprintf( file, "CELLS %d %d\n", num_tris, 4*num_tris );
  for (i = 0; i < num_tris; ++i)
    fprintf( file, "3 %d %d %d\n",
                   triangleConnectivity[3*i    ],
                   triangleConnectivity[3*i + 1],
                   triangleConnectivity[3*i + 2] );
  
  fprintf( file, "CELL_TYPES %d\n", num_tris );
  for (i = 0; i < num_tris; ++i)
    fprintf( file, "5\n" );
    
  fclose( file );
  return true;
}


void iMeshTest::setUp()
{
  const char* TMP_FILE_NAME = "hexagon.vtk";
  MsqPrintError err(cout);
  int ierr;
  
  CPPUNIT_ASSERT(writeVtkFile(TMP_FILE_NAME));
  
  iMesh_newMesh( NULL, &myIMesh, &ierr, 0 ); 
  CPPUNIT_ASSERT_EQUAL( (int)iBase_SUCCESS, ierr );
  
  iMesh_load( myIMesh, 0, TMP_FILE_NAME, NULL, &ierr, strlen(TMP_FILE_NAME), 0 );
  CPPUNIT_ASSERT_EQUAL( (int)iBase_SUCCESS, ierr );
  
  iBase_EntitySetHandle root_set;
  iMesh_getRootSet( myIMesh, &root_set, &ierr );
  CPPUNIT_ASSERT_EQUAL( (int)iBase_SUCCESS, ierr );
          
  myMesh = new MsqIMesh( myIMesh, root_set, iBase_ALL_TYPES, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT( myMesh != NULL );  
  
  matchVertexCoordinates();
  matchElementConnectivity();
}

void iMeshTest::tearDown()
{
  int err;
  delete myMesh;
  myMesh = 0;
  if (myIMesh) {
    iMesh_dtor( myIMesh, &err );
    CPPUNIT_ASSERT_EQUAL( (int)iBase_SUCCESS, err );
  }
}

void iMeshTest::matchVertexCoordinates()
{
  MsqPrintError err(cout);
  const size_t num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );
    // make sure this hasn't been called yet
  CPPUNIT_ASSERT( 0 == vtxHandleToIndex.size() );
    // initialize data
  memset( vtxIndexToHandle, 0, sizeof(vtxIndexToHandle) );
  
    // Get vertex handles
  vector<Mesh::VertexHandle> vertices;
  myMesh->get_all_vertices( vertices, err ); 
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT_EQUAL( num_pts, vertices.size() );
  
    // get vertex coordinates
  vector<MsqVertex> coordinates( num_pts );
  myMesh->vertices_get_coordinates( arrptr(vertices),
                                    arrptr(coordinates),
                                    num_pts,
                                    err );
  CPPUNIT_ASSERT( !err );
  
    // match vertex coordiantes
  for (size_t i = 0; i < num_pts; ++i)
  {
    Mesquite::Vector3D coord( vertexCoords[3*i], vertexCoords[3*i+1], vertexCoords[3*i+2] );
    size_t j;
    for (j = 0; j < vertices.size(); ++j)
    {
      if (((coordinates[j]) - coord).length() < DBL_EPSILON)
      {
        vtxIndexToHandle[i] = vertices[j];
        CPPUNIT_ASSERT(vtxHandleToIndex.find(vertices[j]) == vtxHandleToIndex.end());
        vtxHandleToIndex[vertices[j]] = i;
        break;
      }
    }
    
    CPPUNIT_ASSERT(j < vertices.size()); // found a match
  }
}

bool iMeshTest::match_triangles( const int* tri1, 
                                 const Mesh::VertexHandle* tri2_handles )
{
  int tri2[3] = { vtxHandleToIndex[tri2_handles[0]],
                  vtxHandleToIndex[tri2_handles[1]],
                  vtxHandleToIndex[tri2_handles[2]] };
  
  return (tri1[0] == tri2[0] && tri1[1] == tri2[1] && tri1[2] == tri2[2]) ||
         (tri1[0] == tri2[1] && tri1[1] == tri2[2] && tri1[2] == tri2[0]) ||
         (tri1[0] == tri2[2] && tri1[1] == tri2[0] && tri1[2] == tri2[1]);
}

void iMeshTest::matchElementConnectivity()
{
  MsqPrintError err(cout);
  const size_t num_tri = sizeof(triangleConnectivity) / (3*sizeof(int));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );
    // make sure this hasn't been called yet
  CPPUNIT_ASSERT( 0 == triHandleToIndex.size() );
    // initialize data
  memset( triIndexToHandle, 0, sizeof(triIndexToHandle) );
  
  vector<Mesh::VertexHandle> vertices;
  vector<Mesh::ElementHandle> elements;
  vector<size_t> offsets;
  myMesh->get_all_elements( elements, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( num_tri, elements.size() );
  myMesh->elements_get_attached_vertices( arrptr(elements),
                                          elements.size(),
                                          vertices,
                                          offsets,
                                          err );
  CPPUNIT_ASSERT(!err);
                                          
    // Make sure all are triangles
  size_t i;
  for (i = 0; i < elements.size(); ++i)
    CPPUNIT_ASSERT_EQUAL( offsets[i] + 3, offsets[i+1] );
  
    // Match triangles
  for (size_t i = 0; i < num_tri; ++i)
  {
    size_t j;
    for (j = 0; j < elements.size(); ++j)
    {
      Mesh::VertexHandle verts[3] = {
        vertices[offsets[j]  ],
        vertices[offsets[j]+1],
        vertices[offsets[j]+2] };
      
      if (match_triangles( triangleConnectivity + 3*i, verts ))
      {
        triIndexToHandle[i] = elements[j];
        CPPUNIT_ASSERT(triHandleToIndex.find(elements[j]) == triHandleToIndex.end());
        triHandleToIndex[elements[j]] = i;
        break;
      }
    }
    
    CPPUNIT_ASSERT(j < elements.size()); // found a match
  }
}

/*
void iMeshTest::testVertexIterator() 
{
  MsqPrintError err(cout);
  const size_t num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

    // mark each vertex as it is encountered
  std::vector<int> marks(num_pts);
  memset( arrptr(marks), 0, num_pts * sizeof(int) );
  
    // iterate over vertices
  size_t count = 0;
  VertexIterator* iter = myMesh->vertex_iterator(err);
  CPPUNIT_ASSERT(!err && iter);
  while (!iter->is_at_end())
  {
    Mesh::VertexHandle handle = iter->operator*();
    iter->operator++();
    
    map<Mesh::VertexHandle,int>::iterator f = vtxHandleToIndex.find(handle);
    CPPUNIT_ASSERT( f != vtxHandleToIndex.end() );
    
    unsigned index = f->second;
    CPPUNIT_ASSERT(index < num_pts);
    CPPUNIT_ASSERT(marks[index] == 0);
    marks[index] = 1;
    ++count;
  }
  CPPUNIT_ASSERT_EQUAL(count , num_pts);
}  
*/
void iMeshTest::testVertexByte()
{
  MsqPrintError err(cout);
  const size_t num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

    // set each one individually
  unsigned char bytes[num_pts];
  size_t i;
  for (i = 0; i < num_pts; ++i)
  {
    bytes[i] = (unsigned char)(rand() % 128);
    myMesh->vertex_set_byte( vtxIndexToHandle[i], bytes[i], err );
    CPPUNIT_ASSERT( !err );
  }
  for (i = 0; i < num_pts; ++i)
  {
    unsigned char byte = 0;
    myMesh->vertex_get_byte( vtxIndexToHandle[i], &byte, err );
    CPPUNIT_ASSERT( !err );
    CPPUNIT_ASSERT_EQUAL( byte , bytes[i] );
  }
  
    // set all at once
  for (i = 0; i < num_pts; ++i)
    bytes[i] = (unsigned char)(rand() % 128);
  myMesh->vertices_set_byte( vtxIndexToHandle, bytes, num_pts, err );
  CPPUNIT_ASSERT( !err );
  unsigned char bytes2[num_pts];
  myMesh->vertices_get_byte( vtxIndexToHandle, bytes2, num_pts, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT( !memcmp( bytes, bytes2, num_pts ) );
}  

void iMeshTest::testVertexFlagNone( bool fixed )
{
  if (fixed) {
    myMesh->clear_fixed_tag();
    CPPUNIT_ASSERT( NULL == myMesh->get_fixed_tag() );
  }
  else {
    myMesh->clear_slaved_tag();
    CPPUNIT_ASSERT( NULL == myMesh->get_slaved_tag() );
  }
  
    // get all vertices
  MsqPrintError err(cout);
  std::vector<Mesh::VertexHandle> handles;
  myMesh->get_all_vertices( handles, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!handles.empty());

  std::vector<bool> b(handles.size(), true );
  if (fixed)
    myMesh->vertices_get_fixed_flag( arrptr(handles), b, handles.size(), err );
  else
    myMesh->vertices_get_slaved_flag( arrptr(handles), b, handles.size(), err );
  CPPUNIT_ASSERT_EQUAL( handles.size(), b.size() );
  size_t first_true = std::find( b.begin(), b.end(), true ) - b.begin();
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(  handles.size(), first_true );
}


void iMeshTest::testVertexFlag( bool fixed, iBase_TagValueType type )
{
  int ierr;
  MsqPrintError err(cout);
  iBase_TagHandle tag;

  const char* name = fixed ? "fixed" : "slaved";
  iMesh_createTag( myIMesh, name, 1, type, &tag, &ierr, strlen(name) );
  CPPUNIT_ASSERT_EQUAL( (int)iBase_SUCCESS, ierr );

  if (fixed) {
    myMesh->set_fixed_tag( tag, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT( NULL != myMesh->get_fixed_tag() );
    CPPUNIT_ASSERT_EQUAL( tag, *(myMesh->get_fixed_tag()) );
  }
  else {
    myMesh->set_slaved_tag( tag, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT( NULL != myMesh->get_slaved_tag() );
    CPPUNIT_ASSERT_EQUAL( tag, *(myMesh->get_slaved_tag()) );
  }
  
    // get all vertices
  std::vector<Mesh::VertexHandle> handles;
  myMesh->get_all_vertices( handles, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!handles.empty());
  iBase_EntityHandle* ihandles = reinterpret_cast<iBase_EntityHandle*>(arrptr(handles));
  
  if (type == iBase_INTEGER) {
      // define alternating values for flag
    std::vector<int> values(handles.size(),!fixed);
    for (size_t i = 0; i < handles.size(); i+=2)
      values[i] = fixed;

      // set flag on vertices
    iMesh_setIntArrData( myIMesh, ihandles, handles.size(), tag, arrptr(values), values.size(), &ierr );
    CPPUNIT_ASSERT_EQUAL( (int)iBase_SUCCESS, ierr );
  }
  else if (type == iBase_BYTES ) {
       // define alternating values for flag
    std::vector<char> values(handles.size(),!fixed);
    for (size_t i = 0; i < handles.size(); i+=2)
      values[i] = fixed;

      // set flag on vertices
    iMesh_setArrData( myIMesh, ihandles, handles.size(), tag, arrptr(values), values.size(), &ierr );
    CPPUNIT_ASSERT_EQUAL( (int)iBase_SUCCESS, ierr );
  }
  else {
    CPPUNIT_ASSERT(!"Unexpected tag data type in test code"); 
  }
  
    // get flag through MsqIMesh
  std::vector<bool> flags;
  if (fixed)
    myMesh->vertices_get_fixed_flag( arrptr(handles), flags, handles.size(), err );
  else
    myMesh->vertices_get_slaved_flag( arrptr(handles), flags, handles.size(), err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( handles.size(), flags.size() );
  
    // check flag values
  for (size_t i = 0; i < handles.size(); ++i)
    CPPUNIT_ASSERT_EQUAL( ((!(i%2)) == fixed), (bool)flags[i] );
}

void iMeshTest::testVertexAdjacency()
{
  MsqPrintError err(cout);
  const size_t num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  const size_t num_tri = sizeof(triangleConnectivity) / (3*sizeof(int));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

    // construct adjacency list from input data to compare against
  std::set<int> adjset[num_pts];
  size_t i;
  for (i = 0; i < num_tri; ++i)
    for (size_t j = 3*i; j < 3*i+3; ++j)
      adjset[triangleConnectivity[j]].insert(i);
  
  std::vector<Mesh::ElementHandle> elements;
  std::vector<size_t> offsets;
  myMesh->vertices_get_attached_elements( vtxIndexToHandle, num_pts,
                                          elements, offsets, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL(offsets.size() , num_pts + 1);
  
    // compare connectivity for each vertex
  for (i = 0; i < num_pts; ++i)
  {
    size_t count = offsets[i+1] - offsets[i];
    CPPUNIT_ASSERT_EQUAL(adjset[i].size() , count);
    Mesh::ElementHandle* elems = &elements[offsets[i]];
    
    for (size_t j = 0; j < count; ++j)
    {
        // Get element index from handle
      Mesh::ElementHandle handle = elems[j];
      std::map<Mesh::ElementHandle,int>::iterator idx = triHandleToIndex.find(handle);
      CPPUNIT_ASSERT(idx != triHandleToIndex.end());

        // Find element index in vertex adjacency set
      std::set<int>::iterator iter = adjset[i].find(idx->second);
      CPPUNIT_ASSERT(iter != adjset[i].end());
      adjset[i].erase(iter);
    }
      // Make sure we got all the adjacent elements
    CPPUNIT_ASSERT_EQUAL( (size_t)0, adjset[i].size() );
  }
}


void iMeshTest::testElementConnectivity()
{
  MsqPrintError err(cout);
  //const size_t num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  const size_t num_tri = sizeof(triangleConnectivity) / (3*sizeof(int));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

    // get element connectivity list
  std::vector<Mesh::VertexHandle> vertices;
  std::vector<size_t> offsets;
  myMesh->elements_get_attached_vertices( triIndexToHandle, num_tri,
                                          vertices, offsets, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL(offsets.size() , num_tri + 1);
  CPPUNIT_ASSERT_EQUAL(vertices.size() , 3*num_tri);
  
    // check each element's connectivity
  Mesh::VertexHandle elem_vertices[3];
  for (size_t i = 0; i < num_tri; ++i)
  {
      // check that connectivity list contains three vertices
    CPPUNIT_ASSERT_EQUAL(offsets[i] + 3 , offsets[i+1]);
    
      // get list of vertex indices from connectivity data
    for (size_t j = 0; j < 3; j++)
    {
      size_t offset = offsets[i] + j;
      CPPUNIT_ASSERT( offset < vertices.size() );
      elem_vertices[j] = vertices[offset];
    }
    
      // compare connectivity
    CPPUNIT_ASSERT( match_triangles( triangleConnectivity + 3*i, elem_vertices ) );
  }
}
    
    
void iMeshTest::testElementTopology()
{
  MsqPrintError err(cout);
  const size_t num_tri = sizeof(triangleConnectivity) / (3*sizeof(int));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

  EntityTopology topo[num_tri];
  myMesh->elements_get_topologies( triIndexToHandle, topo, num_tri, err );
  CPPUNIT_ASSERT(!err);
  for (size_t i = 0; i < num_tri; ++i)
    CPPUNIT_ASSERT_EQUAL( topo[i] , Mesquite::TRIANGLE );
}


void iMeshTest::testIntTag()
{
  const char* tagname = "TEST_TEST_INT_TAG";
  MsqPrintError err(cout);
  const size_t num_tri = sizeof(triangleConnectivity) / (3*sizeof(int));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

    // create a tag
  TagHandle tag = myMesh->tag_create( tagname, Mesh::INT, 2, NULL, err );
  CPPUNIT_ASSERT( !err );
  
    // get the tag
  TagHandle tag2 = myMesh->tag_get( tagname, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT_EQUAL( tag , tag2 );
  
    // check tag metadata
  string name;
  Mesh::TagType type;
  unsigned length;
  myMesh->tag_properties( tag, name, type, length, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT_EQUAL( name , std::string(tagname) );
  CPPUNIT_ASSERT_EQUAL( type , Mesh::INT );
  CPPUNIT_ASSERT_EQUAL( length , 2u );
  
    // set the tag on all triangles
  std::vector<int> data1(2*num_tri), data2(2*num_tri);
  std::vector<int>::iterator iter1, iter2;
  for (iter1 = data1.begin(); iter1 != data1.end(); ++iter1)
    *iter1 = rand();
  myMesh->tag_set_element_data( tag, num_tri, triIndexToHandle, arrptr(data1), err );
  CPPUNIT_ASSERT( !err );
  
    // get tag data from all triangles and compare
  myMesh->tag_get_element_data( tag, num_tri, triIndexToHandle, arrptr(data2), err );
  CPPUNIT_ASSERT( !err );
  for (iter1 = data1.begin(), iter2 = data2.begin(); 
       iter1 != data1.end(); ++iter1, ++iter2)
    CPPUNIT_ASSERT_EQUAL( *iter1 , *iter2 );
    
    // destroy the tag
  myMesh->tag_delete( tag, err );
  CPPUNIT_ASSERT(!err);
  tag2 = myMesh->tag_get( tagname, err );
  CPPUNIT_ASSERT_EQUAL(err.error_code() , MsqError::TAG_NOT_FOUND);
  err.clear();
}

  
void iMeshTest::testDoubleTag()
{
  const char* tagname = "TEST_TEST_DOUBLE_TAG";
  MsqPrintError err(cout);
  const size_t num_pts = sizeof(vertexCoords) / (3*sizeof(double));
  
    // check if initialized properly
  CPPUNIT_ASSERT( myMesh );

    // create a tag
  TagHandle tag = myMesh->tag_create( tagname, Mesh::DOUBLE, 1, NULL, err );
  CPPUNIT_ASSERT( !err );
  
    // get the tag
  TagHandle tag2 = myMesh->tag_get( tagname, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT_EQUAL( tag , tag2 );
  
    // check tag metadata
  string name;
  Mesh::TagType type;
  unsigned length;
  myMesh->tag_properties( tag, name, type, length, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT_EQUAL( name , std::string(tagname) );
  CPPUNIT_ASSERT_EQUAL( type , Mesh::DOUBLE );
  CPPUNIT_ASSERT_EQUAL( length , 1u );
  
    // set the tag on all vertices
  std::vector<double> data1(num_pts), data2(num_pts);
  std::vector<double>::iterator iter1, iter2;
  for (iter1 = data1.begin(); iter1 != data1.end(); ++iter1)
    *iter1 = sqrt(abs(rand()));
  myMesh->tag_set_vertex_data( tag, num_pts, vtxIndexToHandle, arrptr(data1), err );
  CPPUNIT_ASSERT( !err );
  
    // get tag data from all vertices and compare
  myMesh->tag_get_vertex_data( tag, num_pts, vtxIndexToHandle, arrptr(data2), err );
  CPPUNIT_ASSERT( !err );
  for (iter1 = data1.begin(), iter2 = data2.begin(); 
       iter1 != data1.end(); ++iter1, ++iter2)
    CPPUNIT_ASSERT_EQUAL( *iter1 , *iter2 );
    
    // destroy the tag
  myMesh->tag_delete( tag, err );
  CPPUNIT_ASSERT(!err);
  tag2 = myMesh->tag_get( tagname, err );
  CPPUNIT_ASSERT_EQUAL(err.error_code() , MsqError::TAG_NOT_FOUND);
  err.clear();
}

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(iMeshTest, "iMeshTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(iMeshTest, "Unit");
