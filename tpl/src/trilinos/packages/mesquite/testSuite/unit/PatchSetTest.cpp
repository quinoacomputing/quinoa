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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */

#include "Mesquite.hpp"
#include "Mesquite_VertexPatches.hpp"
#include "Mesquite_GlobalPatch.hpp"
#include "Mesquite_MeshInterface.hpp"
#include "Mesquite_Instruction.hpp"
#include "Mesquite_MsqError.hpp"
#include "UnitUtil.hpp"

#include <assert.h>
#include <vector>
#include <algorithm>
#include <set>
#include <numeric>
#include <iostream>
using namespace std;


// Define a fake implementation of Mesquite::Mesh to pass
// to the PatchSet implementations.  This implementation 
// provides only the four methods required by the PatchSet
// implementations. 

#define NI(A) MSQ_SETERR(A)(MsqError::NOT_IMPLEMENTED)

using namespace Mesquite;

class FakeMesh : public Mesquite::Mesh
{
public:
  FakeMesh( size_t num_verts );
  ~FakeMesh();
  
  // Specify whether or not calls should fail, to test error behavior
  void should_fail( bool yesno ) { doError = yesno; }
  
  
  void get_all_elements( vector<ElementHandle>& elems, MsqError& err );
  
  void get_all_vertices( vector<VertexHandle>& verts, MsqError& err );
  
  void vertices_get_fixed_flag( const VertexHandle* verts, std::vector<bool>& fixed, size_t n, MsqError& err );
  
  void vertices_get_slaved_flag( const VertexHandle* verts, std::vector<bool>& fixed, size_t n, MsqError& err );
  
  void vertices_get_attached_elements( const VertexHandle* array, size_t len,
                                       vector<ElementHandle>& elems,
                                       vector<size_t>& offsets,
                                       MsqError& err );
  
  int get_geometric_dimension(MsqError& err) { NI(err); return 3; }
  VertexIterator* vertex_iterator(MsqError& err) { NI(err); return 0; }
  ElementIterator* element_iterator(MsqError& err) { NI(err); return 0; }
  void vertices_get_coordinates( const VertexHandle*, MsqVertex*, size_t, MsqError& err ) { NI(err); }
  void vertex_set_coordinates( VertexHandle, const Vector3D&, MsqError& err ) { NI(err); }
  void vertex_set_byte( VertexHandle h, unsigned char b, MsqError& err ) { vertices_get_byte( &h, &b, 1, err ); }
  void vertices_set_byte( const VertexHandle*, const unsigned char*, size_t, MsqError& err );
  void vertex_get_byte( const VertexHandle h, unsigned char* b, MsqError& err) { vertices_get_byte( &h, b, 1, err );}
  void vertices_get_byte( const VertexHandle*, unsigned char*, size_t, MsqError& err );
  
  void elements_get_attached_vertices(const ElementHandle*, size_t,
                                      vector<VertexHandle>&,
                                      vector<size_t>&,
                                      MsqError& err ) { NI(err);}
  void elements_get_topologies(const ElementHandle*, EntityTopology*, size_t, MsqError& err ) { NI(err); }
  TagHandle tag_create( const string&, TagType, unsigned, const void*, MsqError& err ) { NI(err); return 0; }
  void tag_delete( TagHandle, MsqError& err ) { NI(err); }
  TagHandle tag_get( const string&, MsqError& err ) { NI(err); return 0; }
  void tag_properties( TagHandle, string&, TagType&, unsigned&, MsqError& err ) { NI(err); }
  void tag_set_element_data( TagHandle, size_t, const ElementHandle*, const void*, MsqError& err ) { NI(err); }
  void tag_set_vertex_data( TagHandle, size_t, const VertexHandle*, const void*, MsqError& err ) { NI(err); }
  void tag_get_element_data( TagHandle, size_t, const ElementHandle*, void*, MsqError& err ) { NI(err); }
  void tag_get_vertex_data( TagHandle, size_t, const VertexHandle*, void*, MsqError& err ) { NI(err); }
  void release_entity_handles( const EntityHandle*, size_t, MsqError& err ) { NI(err); }
  void release() {}
private:
  std::vector<VertexHandle> vertHandles;
  std::vector<ElementHandle> elemHandles;
  std::vector<size_t> vertOffsets;
  std::vector<bool> fixedFlags;
  std::vector<unsigned char> vertexBytes;
  bool doError;
};

FakeMesh::FakeMesh( size_t num_vtx )
  : doError(false)
{
  vertHandles.resize(num_vtx);
  vertOffsets.resize(num_vtx+1);
  fixedFlags.resize(num_vtx);
  elemHandles.clear();
  for (size_t i = 0; i < num_vtx; ++i)
  {
    vertHandles[i] = (Mesh::VertexHandle)i;
    vertOffsets[i] = elemHandles.size();
    for (size_t j = 0; j < (num_vtx%5); ++j)
      elemHandles.push_back((Mesh::ElementHandle)(i*num_vtx+j));
    fixedFlags[i] = !(i%2);
  }
  vertOffsets[vertOffsets.size()-1] = elemHandles.size();
  vertexBytes.resize( num_vtx, 0 );
}

FakeMesh::~FakeMesh() {}

void FakeMesh::get_all_elements( vector<ElementHandle>& elems, MsqError& err )
{
  if (doError) {
    MSQ_SETERR(err)(MsqError::UNKNOWN_ERROR, "Expected error");
    return;
  }
  
  elems = elemHandles;
}
  
void FakeMesh::get_all_vertices( vector<Mesh::VertexHandle>& verts, MsqError& err )
{
  if (doError) {
    MSQ_SETERR(err)(MsqError::UNKNOWN_ERROR, "Expected error");
    return;
  }
  
  verts = vertHandles;
}
  
void FakeMesh::vertices_get_fixed_flag( const VertexHandle* verts, 
                              std::vector<bool>& fixed, size_t n, MsqError& err )
{
  if (doError) {
    MSQ_SETERR(err)(MsqError::UNKNOWN_ERROR, "Expected error");
    return;
  }
  
  if (!verts) {
    MSQ_SETERR(err)(MsqError::INVALID_STATE, "NULL array pointer");
    return;
  }
  
  fixed.resize(n);
  for (size_t i = 0; i < n; ++i)
  {
    size_t vert = (size_t)verts[i];
    if (vert >= vertHandles.size()) {
      MSQ_SETERR(err)(MsqError::INVALID_STATE, "Vertex handle out of range");
      return;
    }
    fixed[i] = fixedFlags[vert];
  }
}
   
void FakeMesh::vertices_get_slaved_flag( const VertexHandle* , 
                              std::vector<bool>& , size_t , MsqError&  )
{
  CPPUNIT_ASSERT(false);
}
 
void FakeMesh::vertices_get_attached_elements( const VertexHandle* verts, 
                                     size_t n,
                                     vector<ElementHandle>& elems,
                                     vector<size_t>& offsets,
                                     MsqError& err )
{
  if (doError) {
    MSQ_SETERR(err)(MsqError::UNKNOWN_ERROR, "Expected error");
    return;
  }
  
  if (!verts) {
    MSQ_SETERR(err)(MsqError::INVALID_STATE, "NULL array pointer");
    return;
  }
  
  elems.clear();
  offsets.clear();
  for (size_t i = 0; i < n; ++i)
  {
    size_t vert = (size_t)verts[i];
    if (vert >= vertHandles.size()) {
      MSQ_SETERR(err)(MsqError::INVALID_STATE, "Vertex handle out of range");
      return;
    }
    offsets.push_back(elems.size());
    size_t s = vertOffsets[vert];
    size_t e = vertOffsets[vert+1];
    for (size_t j = s; j < e; ++j)
      elems.push_back( elemHandles[j] );
  }
  offsets.push_back( elems.size() );
}
  
void FakeMesh::vertices_get_byte( const VertexHandle* handles, 
                                  unsigned char* bytes, 
                                  size_t count, 
                                  MsqError& err )
{
  for (size_t i = 0; i < count; ++i) {
    size_t vert = (size_t)handles[i];
    if (vert >= vertHandles.size()) {
      MSQ_SETERR(err)(MsqError::INVALID_STATE, "Vertex handle out of range");
      return;
    }
    bytes[i] = vertexBytes[vert];
  }
}
  
void FakeMesh::vertices_set_byte( const VertexHandle* handles, 
                                  const unsigned char* bytes, 
                                  size_t count, 
                                  MsqError& err )
{
  for (size_t i = 0; i < count; ++i) {
    size_t vert = (size_t)handles[i];
    if (vert >= vertHandles.size()) {
      MSQ_SETERR(err)(MsqError::INVALID_STATE, "Vertex handle out of range");
      return;
    }
    vertexBytes[vert] = bytes[i];
  }
}
    


class PatchSetTest : public CppUnit::TestFixture
{
  private:
    CPPUNIT_TEST_SUITE( PatchSetTest );
    
    CPPUNIT_TEST(test_vertex_patches);
    CPPUNIT_TEST(test_vertex_patches_bad_handle);
    CPPUNIT_TEST(test_vertex_patches_fail_handles);
    CPPUNIT_TEST(test_vertex_patches_fail_entities);
    
    CPPUNIT_TEST(test_global_patch);
    CPPUNIT_TEST(test_global_patch_fail_entities);
    
    CPPUNIT_TEST_SUITE_END();
    
    FakeMesh myMesh;
    
    // make sure methods fail if passed a bad handle
    void test_bad_handle( PatchSet& ps );
    
    // test that get_patch_handles propogates failure from Mesh
    void test_fail_handles( PatchSet& ps );
    
    // test that get_patch propogates failures from Mesh
    void test_fail_entities( PatchSet& ps );
  
  public:
  
    PatchSetTest() : myMesh(10) {}
  
    void setUp() { myMesh.should_fail(false); }
    void tearDown() {}
    
    // test that VertexPatches returns expected data
    void test_vertex_patches();
    
    // call failure tests for VertexPatches
    void test_vertex_patches_bad_handle();
    void test_vertex_patches_fail_handles();
    void test_vertex_patches_fail_entities();
    
    // test that GlobalPatch returns expected data
    void test_global_patch();
    
    // call failure tests for GlobalPatch
    void test_global_patch_fail_entities();
};

void PatchSetTest::test_vertex_patches()
{
  size_t i;
  VertexPatches vp;
  MsqPrintError err(std::cout);
  vp.set_mesh( &myMesh );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, 0);
  Instruction::initialize_vertex_byte( &mesh_and_domain, 0, err );
  ASSERT_NO_ERROR(err);
  
    // Get data from myMesh to compare to
  
  vector<Mesh::VertexHandle> vertex_handles, patch_verts;
  vector<Mesh::ElementHandle> element_handles, patch_elems;
  myMesh.get_all_vertices( vertex_handles, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!vertex_handles.empty());
  
  std::vector<bool> fixed;
  myMesh.vertices_get_fixed_flag(arrptr(vertex_handles), fixed, vertex_handles.size(), err );
  ASSERT_NO_ERROR(err);
  
  set<Mesh::VertexHandle> free_verts;
  for (i = 0; i < vertex_handles.size(); ++i)
    if (!fixed[i])
      free_verts.insert( vertex_handles[i] );
  
    // Get list of patch handles
  
  vector<PatchSet::PatchHandle> patch_handles;
  vp.get_patch_handles( patch_handles, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(free_verts.size(), patch_handles.size());
  
  
    // Check each patch handle
  vector<size_t> offsets;
  for (i = 0; i < patch_handles.size(); ++i)
  {
    vp.get_patch(patch_handles[i], patch_elems, patch_verts, err );
    ASSERT_NO_ERROR(err);

      // Check that each patch contains exactly 1 free vertex
      // and that it is always a different free vertex.
    CPPUNIT_ASSERT(patch_verts.size() == 1);
    set<Mesh::VertexHandle>::iterator i = free_verts.find( patch_verts[0] );
    CPPUNIT_ASSERT(i != free_verts.end());
    free_verts.erase(i);
    
      // Get adjacent elements from myMesh to compare with
    element_handles.clear();
    myMesh.vertices_get_attached_elements( arrptr(patch_verts), 1, 
                                           element_handles, offsets, err );
    ASSERT_NO_ERROR(err);
    
      // Compare element handle lists
    sort( element_handles.begin(), element_handles.end() );
    sort( patch_elems.begin(), patch_elems.end() );
    CPPUNIT_ASSERT( element_handles == patch_elems );
  }
}


void PatchSetTest::test_global_patch()
{
  GlobalPatch gp;
  MsqPrintError err(std::cout);
  gp.set_mesh( &myMesh );
  
    // Get data from myMesh to compare to
  vector<Mesh::VertexHandle> vertex_handles, patch_verts;
  vector<Mesh::ElementHandle> element_handles, patch_elems;
  myMesh.get_all_vertices( vertex_handles, err );
  CPPUNIT_ASSERT(!err);
  myMesh.get_all_elements( element_handles, err );
  CPPUNIT_ASSERT(!err);

    // Get list of patch handles
  vector<PatchSet::PatchHandle> patch_handles;
  gp.get_patch_handles( patch_handles, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(1 == patch_handles.size());
  
    // Get mesh data from GlobalPatch
  gp.get_patch( patch_handles[0], patch_elems, patch_verts, err );
  CPPUNIT_ASSERT(!err);
  
    // compare element list
  sort( element_handles.begin(), element_handles.end() );
  sort( patch_elems.begin(), patch_elems.end() );
  CPPUNIT_ASSERT( patch_elems == element_handles );
  
    // compare vertex list
  sort( vertex_handles.begin(), vertex_handles.end() );
  sort( patch_verts.begin(), patch_verts.end() );
  CPPUNIT_ASSERT( patch_verts.empty() || patch_verts == vertex_handles );
}


void PatchSetTest::test_bad_handle( PatchSet& ps )
{
  MsqPrintError err(std::cout);
  ps.set_mesh( &myMesh );
  
    // Get list of patch handles
  vector<PatchSet::PatchHandle> patch_handles;
  ps.get_patch_handles( patch_handles, err );
  CPPUNIT_ASSERT(!err);
  
    // create an invalid handle
  size_t max_handle = (size_t)*max_element( patch_handles.begin(), patch_handles.end() );
  size_t bad_handle = max_handle + 1;
  
    // try to get patch for invalid handle
  vector<Mesh::VertexHandle> patch_verts;
  vector<Mesh::ElementHandle> patch_elems;
  ps.get_patch( (PatchSet::PatchHandle)bad_handle, patch_elems, patch_verts, err );
  CPPUNIT_ASSERT(err);
  err.clear();
}
  
void PatchSetTest::test_fail_handles( PatchSet& ps )
{
  MsqPrintError err(std::cout);
  ps.set_mesh( &myMesh );
  
  myMesh.should_fail(true);
  vector<PatchSet::PatchHandle> patch_handles;
  ps.get_patch_handles( patch_handles, err );
  myMesh.should_fail(false);
  CPPUNIT_ASSERT(err);
  err.clear();
}

void PatchSetTest::test_fail_entities( PatchSet& ps )
{
  MsqPrintError err(std::cout);
  ps.set_mesh( &myMesh );
  
    // Get list of patch handles
  vector<PatchSet::PatchHandle> patch_handles;
  ps.get_patch_handles( patch_handles, err );
  CPPUNIT_ASSERT(!err);
  
    // try to get patch for invalid handle
  vector<Mesh::VertexHandle> patch_verts;
  vector<Mesh::ElementHandle> patch_elems;
  myMesh.should_fail(true);
  ps.get_patch( patch_handles[0], patch_elems, patch_verts, err );
  myMesh.should_fail(false);
  CPPUNIT_ASSERT(err);
  err.clear();
}


void PatchSetTest::test_vertex_patches_bad_handle()
{
  VertexPatches ps;
  test_bad_handle(ps);
}

void PatchSetTest::test_vertex_patches_fail_handles()
{
  VertexPatches ps;
  test_fail_handles(ps);
}

void PatchSetTest::test_vertex_patches_fail_entities()
{
  VertexPatches ps;
  test_fail_entities(ps);
}

void PatchSetTest::test_global_patch_fail_entities()
{
  GlobalPatch ps;
  test_fail_entities(ps);
}

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PatchSetTest, "PatchSetTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PatchSetTest, "Unit");

