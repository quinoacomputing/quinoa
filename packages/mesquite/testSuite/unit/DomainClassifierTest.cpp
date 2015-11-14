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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file DomainClassifierTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "UnitUtil.hpp"
#include "Mesquite_DomainClassifier.hpp"
#include "Mesquite_MeshDomain1D.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_MsqVertex.hpp"
#include "Mesquite_PlanarDomain.hpp"

#include <stdio.h>
#include <iostream>
#include <algorithm>

using namespace Mesquite;
using namespace std;


class DomainClassifierTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(DomainClassifierTest);
  CPPUNIT_TEST (test_classify_by_handle);
  CPPUNIT_TEST (test_valid_classification);
  CPPUNIT_TEST (test_classify_by_tag);
  CPPUNIT_TEST (test_classify_skin);
  CPPUNIT_TEST (test_classify_by_geometry);
  CPPUNIT_TEST_SUITE_END();
  
public:

  typedef DomainClassifier::DomainSet DomSet;
  typedef std::vector<DomSet> DomSetList;
  MeshImpl myMesh;
  DomSetList myDomains;
  std::vector<int> domainDims;
  
  void setUp();
  void tearDown();

  void test_classify_by_handle();
  void test_valid_classification();
  void test_classify_by_tag();
  void test_classify_skin();
  void test_classify_by_geometry();

  void check_domain( DomainClassifier& dom );
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(DomainClassifierTest, "Unit");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(DomainClassifierTest, "DomainClassifierTest");

static void print_domain( int i, DomainClassifier::DomainSet& set )
{
  if (PointDomain* pd = dynamic_cast<PointDomain*>(set.domain)) 
    printf("%d: PointDomain( %f %f %f ) @ %p\n", i, pd->geom()[0], pd->geom()[1], pd->geom()[2], pd );
  else if (LineDomain* ld = dynamic_cast<LineDomain*>(set.domain)) {
    if (fabs(ld->geom().direction()[0]) < 1e-6 && fabs(ld->geom().direction()[1]) < 1e-6)
      printf("%d: LineDomain( x = %f, y = %f ) @ %p\n", i, ld->geom().point()[0], ld->geom().point()[1], ld );
    else if (fabs(ld->geom().direction()[1]) < 1e-6 && fabs(ld->geom().direction()[2]) < 1e-6)
      printf("%d: LineDomain( y = %f, z = %f ) @ %p\n", i, ld->geom().point()[1], ld->geom().point()[2], ld );
    else if (fabs(ld->geom().direction()[0]) < 1e-6 && fabs(ld->geom().direction()[2]) < 1e-6)
      printf("%d: LineDomain( x = %f, z = %f ) @ %p\n", i, ld->geom().point()[0], ld->geom().point()[2], ld );
    else
      printf("%d: LineDomain( ? ) @ %p\n", i, ld );
  }
  else if (PlanarDomain* pd = dynamic_cast<PlanarDomain*>(set.domain)) {
    if (fabs(pd->get_normal()[0]) < 1e-6 && fabs(pd->get_normal()[1]) < 1e-6)
      printf("%d: PlanarDomain( z = %f ) @ %p\n", i, pd->get_origin()[2], pd );
    else if (fabs(pd->get_normal()[1]) < 1e-6 && fabs(pd->get_normal()[2]) < 1e-6)
      printf("%d: PlanarDomain( x = %f ) @ %p\n", i, pd->get_origin()[1], pd );
    else if (fabs(pd->get_normal()[0]) < 1e-6 && fabs(pd->get_normal()[2]) < 1e-6)
      printf("%d: PlanarDomain( y = %f ) @ %p\n", i, pd->get_origin()[0], pd );
    else
      printf("%d: PlanarDomain( ? ) @ %p\n", i, pd );
  }
  else {
    printf("%d: unknown domain type @ %p\n", i, set.domain );
  }
  
  if (!set.vertices.empty()) {
    printf("  vertices: ");
    std::vector<Mesh::VertexHandle>::iterator vi = set.vertices.begin();
    for (; vi != set.vertices.end(); ++vi)
      printf( "%lu, ", (unsigned long)*vi);
    printf("\n");
  }
  if (!set.elements.empty()) {
    printf("  elements: ");
    std::vector<Mesh::ElementHandle>::iterator ei = set.elements.begin();
    for (; ei != set.elements.end(); ++ei)
      printf( "%lu, ", (unsigned long)*ei);
    printf("\n");
  }
}
    


void DomainClassifierTest::setUp()
{
  myMesh.clear();
  myDomains.clear();
  domainDims.clear();
    // vertex coodinates
  const char vertex_data[] =
  "POINTS 64 float\n"
  "0 0 0  1 0 0  2 0 0  3 0 0\n"
  "0 1 0  1 1 0  2 1 0  3 1 0\n"
  "0 2 0  1 2 0  2 2 0  3 2 0\n"
  "0 3 0  1 3 0  2 3 0  3 3 0\n"
  "\n"
  "0 0 1  1 0 1  2 0 1  3 0 1\n"
  "0 1 1  1 1 1  2 1 1  3 1 1\n"
  "0 2 1  1 2 1  2 2 1  3 2 1\n"
  "0 3 1  1 3 1  2 3 1  3 3 1\n"
  "\n"
  "0 0 2  1 0 2  2 0 2  3 0 2\n"
  "0 1 2  1 1 2  2 1 2  3 1 2\n"
  "0 2 2  1 2 2  2 2 2  3 2 2\n"
  "0 3 2  1 3 2  2 3 2  3 3 2\n"
  "\n"
  "0 0 3  1 0 3  2 0 3  3 0 3\n"
  "0 1 3  1 1 3  2 1 3  3 1 3\n"
  "0 2 3  1 2 3  2 2 3  3 2 3\n"
  "0 3 3  1 3 3  2 3 3  3 3 3\n"
  "\n";
    // quad connectivity for quads on mesh skin
  const int num_quads = 9*6; // nine per side
  const char quad_data[] = 
  "4  1  0  4  5\n" // -z face (z == 0)
  "4  2  1  5  6\n"
  "4  3  2  6  7\n"
  "4  5  4  8  9\n"
  "4  6  5  9 10\n"
  "4  7  6 10 11\n"
  "4  9  8 12 13\n"
  "4 10  9 13 14\n"
  "4 11 10 14 15\n" 
  "\n"
  "4 48 49 53 52\n" // +z face (z == 3)
  "4 49 50 54 53\n"
  "4 50 51 55 54\n"
  "4 52 53 57 56\n"
  "4 53 54 58 57\n"
  "4 54 55 59 58\n"
  "4 56 57 61 60\n"
  "4 57 58 62 61\n"
  "4 58 59 63 62\n" 
  "\n"
  "4  0  1 17 16\n" // -y face (y == 0)
  "4  1  2 18 17\n"
  "4  2  3 19 18\n"
  "4 16 17 33 32\n"
  "4 17 18 34 33\n"
  "4 18 19 35 34\n"
  "4 32 33 49 48\n"
  "4 33 34 50 49\n"
  "4 34 35 51 50\n"
  "\n"
  "4 13 12 28 29\n" // +y face (y == 3)
  "4 14 13 29 30\n"
  "4 15 14 30 31\n"
  "4 29 28 44 45\n"
  "4 30 29 45 46\n"
  "4 31 30 46 47\n"
  "4 45 44 60 61\n"
  "4 46 45 61 62\n"
  "4 47 46 62 63\n"
  "\n"
  "4  4  0 16 20\n" // -x face (x == 0)
  "4  8  4 20 24\n"
  "4 12  8 24 28\n"
  "4 20 16 32 36\n"
  "4 24 20 36 40\n"
  "4 28 24 40 44\n"
  "4 36 32 48 52\n"
  "4 40 36 52 56\n"
  "4 44 40 56 60\n"
  "\n"
  "4  3  7 23 19\n" // +x face (x == 3)
  "4  7 11 27 23\n"
  "4 11 15 31 27\n"
  "4 19 23 39 35\n"
  "4 23 27 43 39\n"
  "4 27 31 47 43\n"
  "4 35 39 55 51\n"
  "4 39 43 59 55\n"
  "4 43 47 63 59\n"
  "\n";
    // hexahedron connectivity
  const int num_hexes = 3*3*3;
  const char hex_data[] =
  "8  0  1  5  4 16 17 21 20\n"
  "8  1  2  6  5 17 18 22 21\n"
  "8  2  3  7  6 18 19 23 22\n"
  "8  4  5  9  8 20 21 25 24\n"
  "8  5  6 10  9 21 22 26 25\n"
  "8  6  7 11 10 22 23 27 26\n"
  "8  8  9 13 12 24 25 29 28\n"
  "8  9 10 14 13 25 26 30 29\n"
  "8 10 11 15 14 26 27 31 30\n"
  "\n"
  "8 16 17 21 20 32 33 37 36\n"
  "8 17 18 22 21 33 34 38 37\n"
  "8 18 19 23 22 34 35 39 38\n"
  "8 20 21 25 24 36 37 41 40\n"
  "8 21 22 26 25 37 38 42 41\n"
  "8 22 23 27 26 38 39 43 42\n"
  "8 24 25 29 28 40 41 45 44\n"
  "8 25 26 30 29 41 42 46 45\n"
  "8 26 27 31 30 42 43 47 46\n"
  "\n"
  "8 32 33 37 36 48 49 53 52\n"
  "8 33 34 38 37 49 50 54 53\n"
  "8 34 35 39 38 50 51 55 54\n"
  "8 36 37 41 40 52 53 57 56\n"
  "8 37 38 42 41 53 54 58 57\n"
  "8 38 39 43 42 54 55 59 58\n"
  "8 40 41 45 44 56 57 61 60\n"
  "8 41 42 46 45 57 58 62 61\n"
  "8 42 43 47 46 58 59 63 62\n"
  "\n";
    // a few interior quads
  const int num_interior_quads = 3;
  const char interior_quad_data[] = 
  "4  1  5 25 17\n"
  "4  4  5 25 24\n"
  "4 16 17 25 24\n"
  "\n";
  
  
  const char filename[] = "dctest.vtk";
  FILE* file = fopen( filename, "w" );
  fputs( "# vtk DataFile Version 2.0\n", file );
  fputs( "Mesquite Mesh\n", file );
  fputs( "ASCII\n", file );
  fputs( "DATASET UNSTRUCTURED_GRID\n", file );
  fputs( vertex_data, file );
  
  int num_elem = num_quads + num_hexes + num_interior_quads;
  int num_elem_data = 5*num_quads + 9*num_hexes * 5*num_interior_quads;
  fprintf( file, "CELLS %d %d\n", num_elem, num_elem_data );
  fputs( quad_data, file );
  fputs( hex_data, file );
  fputs( interior_quad_data, file );
  fprintf( file, "CELL_TYPES %d\n", num_elem );
  for (int i = 0; i < num_quads; ++i)
    fputs( "9\n", file );
  for (int i = 0; i < num_hexes; ++i)
    fputs( "12\n", file );
  for (int i = 0; i < num_interior_quads; ++i)
    fputs( "9\n", file );
  
  fclose( file );
  MsqPrintError err(std::cerr);
  myMesh.read_vtk( filename, err );
  remove( filename );
  CPPUNIT_ASSERT(!err);

  std::vector<Mesh::VertexHandle> verts;
  std::vector<Mesh::ElementHandle> elems;
  myMesh.get_all_vertices(verts, err);
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( (size_t)64, verts.size() );
  myMesh.get_all_elements(elems, err);
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( (size_t)num_elem, elems.size() );
  
    // define point domains
  PointDomain* pdom[8];
  pdom[0] = new PointDomain( Vector3D(0,0,0) );
  pdom[1] = new PointDomain( Vector3D(3,0,0) );
  pdom[2] = new PointDomain( Vector3D(0,3,0) );
  pdom[3] = new PointDomain( Vector3D(3,3,0) );
  pdom[4] = new PointDomain( Vector3D(0,0,3) );
  pdom[5] = new PointDomain( Vector3D(3,0,3) );
  pdom[6] = new PointDomain( Vector3D(0,3,3) );
  pdom[7] = new PointDomain( Vector3D(3,3,3) );
  size_t pdidx[8] = { 0, 3, 12, 15, 48, 51, 60, 63 };
  for (unsigned i = 0; i < 8; ++i) {
    MsqVertex coords;
    Mesh::VertexHandle h = verts[pdidx[i]];
    myMesh.vertices_get_coordinates( &h, &coords, 1, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_VECTORS_EQUAL( pdom[i]->geom(), coords, 1e-6 );
    DomSet set;
    set.domain = pdom[i];
    set.vertices.push_back( h );
    myDomains.push_back( set );
    domainDims.push_back( 0 );
  }
  
    // define line domains
  LineDomain* ldom[12];
  ldom[0] = new LineDomain( Vector3D(0,0,0), Vector3D(1,0,0) ); // y=0,z=0
  ldom[1] = new LineDomain( Vector3D(0,3,0), Vector3D(1,0,0) ); // y=3,z=0
  ldom[2] = new LineDomain( Vector3D(0,0,3), Vector3D(1,0,0) ); // y=0,z=3
  ldom[3] = new LineDomain( Vector3D(0,3,3), Vector3D(1,0,0) ); // y=3,z=3
  ldom[4] = new LineDomain( Vector3D(0,0,0), Vector3D(0,1,0) ); // x=0,z=0
  ldom[5] = new LineDomain( Vector3D(3,0,0), Vector3D(0,1,0) ); // x=3,z=0
  ldom[6] = new LineDomain( Vector3D(0,0,3), Vector3D(0,1,0) ); // x=0,z=3
  ldom[7] = new LineDomain( Vector3D(3,0,3), Vector3D(0,1,0) ); // x=3,z=3
  ldom[8] = new LineDomain( Vector3D(0,0,0), Vector3D(0,0,1) ); // x=0,y=0
  ldom[9] = new LineDomain( Vector3D(3,0,0), Vector3D(0,0,1) ); // x=3,y=0
  ldom[10]= new LineDomain( Vector3D(0,3,0), Vector3D(0,0,1) ); // x=0,y=3
  ldom[11]= new LineDomain( Vector3D(3,3,0), Vector3D(0,0,1) ); // x=3,y=3
  size_t ldidx[12][2] = { {  1,  2 }, { 13, 14 }, { 49, 50 }, { 61, 62 },
                          {  4,  8 }, {  7, 11 }, { 52, 56 }, { 55, 59 },
                          { 16, 32 }, { 19, 35 }, { 28, 44 }, { 31, 47 } };
  for (unsigned i = 0; i < 12; ++i) {
    Mesh::VertexHandle v[2];
    v[0] = verts[ldidx[i][0]];
    v[1] = verts[ldidx[i][1]];
    MsqVertex coords[2];
    myMesh.vertices_get_coordinates( v, coords, 2, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, ldom[i]->geom().distance( coords[0] ), 1e-6 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, ldom[i]->geom().distance( coords[1] ), 1e-6 );
    DomSet set;
    set.domain = ldom[i];
    set.vertices.push_back(v[0]);
    set.vertices.push_back(v[1]);
    myDomains.push_back( set );
    domainDims.push_back( 1 );
  }
  
    // define planar domains
  PlanarDomain* sdom[6];
  sdom[0] = new PlanarDomain( Vector3D( 0, 0,-1), Vector3D(0,0,0) );
  sdom[1] = new PlanarDomain( Vector3D( 0, 0, 1), Vector3D(0,0,3) );
  sdom[2] = new PlanarDomain( Vector3D( 0,-1, 0), Vector3D(0,0,0) );
  sdom[3] = new PlanarDomain( Vector3D( 0, 1, 0), Vector3D(0,3,0) );
  sdom[4] = new PlanarDomain( Vector3D(-1, 0, 0), Vector3D(0,0,0) );
  sdom[5] = new PlanarDomain( Vector3D( 1, 0, 0), Vector3D(3,0,0) );
  size_t sdidx[6][4] = { {  5,  6,  9, 10 }, { 53, 54, 57, 58 },
                         { 17, 18, 33, 34 }, { 29, 30, 45, 46 },
                         { 20, 24, 36, 40 }, { 23, 27, 39, 43 } };
  for (unsigned i = 0; i < 6; ++i) {
    DomSet set;
    set.domain = sdom[i];
    for (unsigned j = 0; j < 4; ++j)
      set.vertices.push_back( verts[sdidx[i][j]] );
    for (unsigned j = 0; j < 9; ++j)
      set.elements.push_back( elems[9*i+j] );
    myDomains.push_back( set );
    domainDims.push_back( 2 );
  }
  
  
//  for (unsigned i = 0; i < myDomains.size(); ++i) 
//    print_domain( i, myDomains[i] );
}

void DomainClassifierTest::tearDown()
{
  myMesh.clear();
  while (!myDomains.empty()) {
    delete myDomains.back().domain;
    myDomains.pop_back();
  }
  domainDims.clear();
}

void DomainClassifierTest::check_domain( DomainClassifier& domain )
{
  std::vector<Mesh::VertexHandle> vertices, cverts;
  std::vector<Mesh::ElementHandle> elements, celems;
  
    // Check that, for each entity with a domain, the 
    // DomainClassifier instance returns that domain.
    // Also, put all entities with domains into cverts and
    // celems for later.
  for (unsigned i = 0; i < myDomains.size(); ++i) {
    for (unsigned j = 0; j < myDomains[i].vertices.size(); ++j) {
      Mesh::VertexHandle v = myDomains[i].vertices[j];
      const MeshDomain* ptr = domain.find_vertex_domain( v );
      CPPUNIT_ASSERT( myDomains[i].domain == ptr );
      cverts.push_back( v );
    }
    for (unsigned k = 0; k < myDomains[i].elements.size(); ++k) {
      Mesh::ElementHandle e = myDomains[i].elements[k];
      const MeshDomain* ptr = domain.find_element_domain( e );
      CPPUNIT_ASSERT( myDomains[i].domain == ptr );
      celems.push_back( e );
    }
  }
  
    // sort cverts and celems so we can do binary_search later
  std::sort( cverts.begin(), cverts.end() );
  std::sort( celems.begin(), celems.end() );
    // get all vertices and elements in mesh
  MsqPrintError err(std::cerr);
  myMesh.get_all_vertices( vertices, err );
  CPPUNIT_ASSERT(!err);
  myMesh.get_all_elements( elements, err );
  CPPUNIT_ASSERT(!err);
  
    // For each vertex not in a domain (not in cverts), make sure
    // that the domain is NULL.
  for (size_t i = 0; i < vertices.size(); ++i) {
    if (std::binary_search( cverts.begin(), cverts.end(), vertices[i] ))
      continue;
    
    const MeshDomain* ptr = domain.find_vertex_domain( vertices[i] );
    CPPUNIT_ASSERT( NULL == ptr );
  }
    // For each element not in a domain (not in celems), make sure
    // that the domain is NULL.
  for (size_t i = 0; i < elements.size(); ++i) {
    if (std::binary_search( celems.begin(), celems.end(), elements[i] ))
      continue;
    
    const MeshDomain* ptr = domain.find_element_domain( elements[i] );
    CPPUNIT_ASSERT( NULL == ptr );
  }
}

void DomainClassifierTest::test_classify_by_handle()
{
  MsqPrintError err(std::cerr);
  DomainClassifier domain;
  DomainClassifier::classify_by_handle( domain, &myMesh,
                                        arrptr(myDomains), myDomains.size(),
                                        err );
  CPPUNIT_ASSERT(!err);
  
  check_domain( domain );
}

void DomainClassifierTest::test_valid_classification()
{
  MsqPrintError err(std::cerr);
  DomainClassifier domain;
  DomainClassifier::classify_by_handle( domain, &myMesh,
                                        arrptr(myDomains), myDomains.size(),
                                        err );
  CPPUNIT_ASSERT(!err);

  domain.test_valid_classification( &myMesh, err );
  CPPUNIT_ASSERT(!err);
}

void DomainClassifierTest::test_classify_by_tag()
{
  CPPUNIT_ASSERT( !myDomains.empty() );
  
  MsqPrintError err(std::cerr);
  int def = myDomains.size();
  TagHandle tag = myMesh.tag_create( "domain", Mesh::INT, 1, &def, err );
  CPPUNIT_ASSERT(!err);
  
  std::vector<MeshDomain*> dom_list;
  std::vector<int> id_list;
  for (unsigned i = 0; i < myDomains.size(); ++i) {
    std::vector<int> vtx_data( myDomains[i].vertices.size(), i );
    std::vector<int> elm_data( myDomains[i].elements.size(), i );
    if (!vtx_data.empty()) {
      myMesh.tag_set_vertex_data( tag, vtx_data.size(), 
                                  &(myDomains[i].vertices[0]),
                                  arrptr(vtx_data), err );
      CPPUNIT_ASSERT(!err);
    }
    if (!elm_data.empty()) {
      myMesh.tag_set_element_data( tag, elm_data.size(), 
                                   &(myDomains[i].elements[0]),
                                   arrptr(elm_data), err );
      CPPUNIT_ASSERT(!err);
    }
    
    dom_list.push_back( myDomains[i].domain );
    id_list.push_back( i );
  }

  DomainClassifier domain;
  DomainClassifier::classify_by_tag( domain, &myMesh,
                                     "domain", 
                                     arrptr(dom_list),
                                     arrptr(id_list),
                                     myDomains.size(),
                                     err );
  CPPUNIT_ASSERT(!err);
  
  check_domain( domain );
}

void DomainClassifierTest::test_classify_skin()
{
  CPPUNIT_ASSERT( !myDomains.empty() );

  MsqPrintError err(std::cerr);
  std::vector<MeshDomain*> arr( myDomains.size() );
  for (size_t i = 0; i < myDomains.size(); ++i)
    arr[i] = myDomains[i].domain;
    
  DomainClassifier domain;
  DomainClassifier::classify_skin_geometrically( domain, &myMesh,
                                        1e-6,
                                        arrptr(arr), arrptr(domainDims),
                                        arr.size(), err );
  CPPUNIT_ASSERT(!err);

  check_domain( domain );
}

void DomainClassifierTest::test_classify_by_geometry()
{
  CPPUNIT_ASSERT( !myDomains.empty() );

  MsqPrintError err(std::cerr);
  std::vector<MeshDomain*> arr( myDomains.size() );
  for (size_t i = 0; i < myDomains.size(); ++i)
    arr[i] = myDomains[i].domain;
    
  DomainClassifier domain;
  DomainClassifier::classify_geometrically( domain, &myMesh,
                                        1e-6,
                                        arrptr(arr), arrptr(domainDims),
                                        arr.size(), err );
  CPPUNIT_ASSERT(!err);

  check_domain( domain );
}

