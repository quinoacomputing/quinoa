#include "Mesquite_BoundedCylinderDomain.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_MsqVertex.hpp"
#include <cppunit/extensions/HelperMacros.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

const double EPSILON = 1e-6;
#define ASSERT_VECTORS_EQUAL( A, B ) \
  CPPUNIT_ASSERT( ((A)-(B)).length() < EPSILON )

using namespace Mesquite;

class BoundedCylinderDomainTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE( BoundedCylinderDomainTest );
  CPPUNIT_TEST( test_domain_DoF );
  CPPUNIT_TEST( test_z_snap_to );
  CPPUNIT_TEST( test_x_snap_to );
  CPPUNIT_TEST( test_create_curve_from_mesh );
  CPPUNIT_TEST_SUITE_END();

public:

  BoundedCylinderDomainTest() {}
  
  void setUp() {}
  void tearDown() {}
  
  void test_domain_DoF();
  
  void test_z_snap_to();
  void test_x_snap_to();
  
  void test_create_curve_from_mesh();
};
 
void BoundedCylinderDomainTest::test_domain_DoF()
{
  MsqError err;
  const Mesh::EntityHandle list[] = { (void*)5, (void*)9, (void*)1, (void*)2, (void*)3 };
  const size_t len = sizeof(list)/sizeof(list[0]); 
  
  BoundedCylinderDomain d( 1000 );
  std::vector<Mesh::EntityHandle> handles;
  std::copy( list, list+len, std::back_inserter( handles ) );
  d.create_curve( 10, handles );
  
  const Mesh::EntityHandle list2[] = 
    { (void*)1, (void*)2, (void*)3, (void*)4, (void*)5, (void*)6, (void*)7, (void*)8,(void*) 9, (void*)10 };
  const size_t len2 = sizeof(list2)/sizeof(list2[0]);
  unsigned short dof_list[len2];
  d.domain_DoF( list2, dof_list, len2, err );
  CPPUNIT_ASSERT(!err);
  
  for (size_t i = 0; i < len2; ++i)
  {
    Mesh::EntityHandle handle = list2[i];
    const Mesh::EntityHandle* f = std::find( list, list + len, handle );
    unsigned short expected_dim = (f < list+len) ? 1 : 2;
    CPPUNIT_ASSERT_EQUAL( expected_dim, dof_list[i] );
  }
}

void BoundedCylinderDomainTest::test_z_snap_to()
{
  const Mesh::EntityHandle HANDLE = (void*)20;
  const double ZPLANE = 5;
  
    // Define a circle in the Z=5 plane, centered at (0,0,5)
    // with a radius of 10.
  std::vector<Mesh::EntityHandle> handles(1);
  handles[0] = HANDLE;
  BoundedCylinderDomain d( 10 );
  d.create_curve( ZPLANE, handles );
  
  Vector3D point;
  
    // try a point inside and below the circle
  point.set( 5, 5, 2 );
  d.snap_to( HANDLE, point );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( ZPLANE, point.z(), EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( point.x(), point.y(), EPSILON );
  point[2] = 0.0;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( d.radius(), point.length(), EPSILON );
  
    // try a point outside and above the circle
  point.set( 50, 50, 50 );
  d.snap_to( HANDLE, point );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( ZPLANE, point.z(), EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( point.x(), point.y(), EPSILON );
  point[2] = 0.0;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( d.radius(), point.length(), EPSILON );

    // try a point on the circle
  Vector3D on_circ(d.radius(), 0, ZPLANE);
  point = on_circ;
  d.snap_to( HANDLE, point );
  ASSERT_VECTORS_EQUAL( on_circ, point );
  
    // try a point at the center of the circle
  point.set( d.radius(), 0, ZPLANE );
  d.snap_to( HANDLE, point );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( ZPLANE, point.z(), EPSILON );
  point[2] = 0.0;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( d.radius(), point.length(), EPSILON);
  
    // try a different handle to make sure a handle not
    // on the the circle results in a point not on the circle
  double pplane = ZPLANE - 5;
  point.set( 0, 0, pplane );
  d.snap_to( (void*)((size_t)HANDLE+1), point );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( pplane, point.z(), EPSILON );
} 
     
  

void BoundedCylinderDomainTest::test_x_snap_to()
{
  const Mesh::EntityHandle HANDLE = 0;
  const double XPLANE = -3;
  
    // Define a circle in the X=-3 plane, centered at (-3,1,1)
    // with a radius of 1
  std::vector<Mesh::EntityHandle> handles(1);
  handles[0] = HANDLE;
  BoundedCylinderDomain d( 1, Vector3D(1,0,0), Vector3D(0,1,1) );
  d.create_curve( XPLANE, handles );
  
  Vector3D point, center(-3,1,1);
  
    // try a point inside and below the circle
  point.set( XPLANE-2, 1.5, 1.5 );
  d.snap_to( HANDLE, point );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( XPLANE, point.x(), EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( point.y(), point.z(), EPSILON );
  point -= center;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( d.radius(), point.length(), EPSILON );
  
    // try a point outside and above the circle
  point.set( XPLANE+20, 50, 50 );
  d.snap_to( HANDLE, point );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( XPLANE, point.x(), EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( point.y(), point.z(), EPSILON );
  point -= center;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( d.radius(), point.length(), EPSILON );

    // try a point on the circle
  Vector3D on_circ(XPLANE, d.radius(), 0);
  point = on_circ;
  d.snap_to( HANDLE, point );
  ASSERT_VECTORS_EQUAL( on_circ, point );
  
    // try a point at the center of the circle
  point = center;
  d.snap_to( HANDLE, point );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( XPLANE, point.x(), EPSILON );
  point -= center;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( d.radius(), point.length(), EPSILON);
  
    // try a different handle to make sure a handle not
    // on the the circle results in a point not on the circle
  double pplane = XPLANE + 5;
  point.set( pplane, 0, 0 );
  d.snap_to( (void*)((size_t)HANDLE+1), point );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( pplane, point.x(), EPSILON );
} 


// define a VDK file containing four quads quartering a
// square with corners at ( +/-1, +/-1, 0 )
const char vtk_file[] = 
"# vtk DataFile Version 2.0\n"
"Mesquite Mesh\n"
"ASCII\n"
"DATASET UNSTRUCTURED_GRID\n"
"POINTS 9 float\n"
" 0  0 0\n"
" 1  0 0\n"
" 1  1 0\n"
" 0  1 0\n"
"-1  1 0\n"
"-1  0 0\n"
"-1 -1 0\n"
" 0 -1 0\n"
" 1 -1 0\n"
"CELLS 4 20\n"
"4 0 1 2 3\n"
"4 0 3 4 5\n"
"4 0 5 6 7\n"
"4 0 7 8 1\n"
"CELL_TYPES 4\n"
"9 9 9 9\n"
;

void BoundedCylinderDomainTest::test_create_curve_from_mesh()
{
    // create the input file to read
  char filename[] = "BCDTestTMP";
  FILE* file = fopen( filename, "w" );
  size_t s = fwrite( vtk_file, sizeof(vtk_file), 1, file );
  if (!s)
  {
    remove( filename );
    CPPUNIT_ASSERT(s);
  }
  fclose( file );
  
    // read the file into a Mesh object
  MsqPrintError err(std::cout);
  MeshImpl mesh;
  mesh.read_vtk( filename, err );
  remove( filename );
  CPPUNIT_ASSERT( !err );
  
    // create a domain such that the four corner vertices
    // of the 4-quad square mesh lie on the bounding curves
    // of the cylinder.
  BoundedCylinderDomain d( 1, Vector3D( 0, 1, 0 ) );
  d.create_curve( 1, &mesh, EPSILON, err );
  CPPUNIT_ASSERT(!err);
  d.create_curve( -1, &mesh, EPSILON, err );
  CPPUNIT_ASSERT(!err);
  
    // Get list of vertex coordinates from mesh
  std::vector<Mesh::VertexHandle> handles;
  mesh.get_all_vertices( handles, err ); 
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( handles.size(), (size_t)9 );
  std::vector<MsqVertex> vertices( handles.size() );
  mesh.vertices_get_coordinates( arrptr(handles), arrptr(vertices), handles.size(), err );
  CPPUNIT_ASSERT(!err);
  
    // Get the dimension of the geometry each vertex is bound to
  std::vector<unsigned short> dims( handles.size() );
  d.domain_DoF( arrptr(handles), arrptr(dims), handles.size(), err );
  CPPUNIT_ASSERT(!err);
  
    // Check that the four corner vertices are on the curves
    // and all others are not
  for (size_t i = 0; i < handles.size(); ++i)
  {
    if ( fabs(vertices[i].x()) < EPSILON ||
         fabs(vertices[i].y()) < EPSILON )
      CPPUNIT_ASSERT_EQUAL( (unsigned short)2, dims[i] );
    else 
      CPPUNIT_ASSERT_EQUAL( (unsigned short)1, dims[i] );
  }
}  

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( BoundedCylinderDomainTest, "BoundedCylinderDomainTest" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( BoundedCylinderDomainTest, "Unit" );
  
