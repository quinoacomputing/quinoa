

#include "Mesquite.hpp"
#include "PatchData.hpp"
#include "SphericalDomain.hpp"
#include "BoundedCylinderDomain.hpp"
#include <cppunit/extensions/HelperMacros.h>
#include <iostream>

const double EPSILON = 1e-6;
#define ASSERT_VECTORS_EQUAL( A, B ) \
  CPPUNIT_ASSERT( ((A)-(B)).length() < EPSILON )

using namespace Mesquite;

using std::cout;
using std::cerr;
using std::endl;

class PatchDataTestNormals : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(PatchDataTestNormals);
  //CPPUNIT_TEST (test_get_vertex_normals_infinite_domain);
  //CPPUNIT_TEST (test_get_vertex_normals_bounded_domain);
  CPPUNIT_TEST (test_get_corner_normals_infinite_domain);
  CPPUNIT_TEST (test_get_corner_normals_bounded_domain);
  CPPUNIT_TEST (test_get_element_normals_infinite_domain);
  CPPUNIT_TEST (test_get_element_normals_bounded_domain);
  CPPUNIT_TEST_SUITE_END();
   
  PatchData boundedMesh, unboundedMesh;
  BoundedCylinderDomain boundedDomain;
  SphericalDomain unboundedDomain;
  

public:

  PatchDataTestNormals()
    : boundedDomain( 1 ), 
      unboundedDomain( Vector3D(0,0,0), 1 )
    {}

  void setUp();
  
  void tearDown() {}

  //void test_get_vertex_normals_infinite_domain();
  //void test_get_vertex_normals_bounded_domain();
  void test_get_corner_normals_infinite_domain();
  void test_get_corner_normals_bounded_domain();
  void test_get_element_normals_infinite_domain();
  void test_get_element_normals_bounded_domain();
  
};

void PatchDataTestNormals::setUp()
{
  MsqPrintError err(cout);
  
    // Define a mesh on the unit sphere
    // Make six quads corresponding to the six faces
    // of a cube inscribed in the sphere
  const double T = 1.0 / sqrt(3.0);
  double ucoords[] = {  T, -T, -T,
                        T,  T, -T,
                       -T,  T, -T,
                       -T, -T, -T,
                        T, -T,  T,
                        T,  T,  T,
                       -T,  T,  T,
                       -T, -T,  T };
  size_t uconn[] = { 3, 2, 1, 0, // -Z face
                     4, 5, 6, 7, // +Z face
                     0, 1, 5, 4, // +X face
                     1, 2, 6, 5, // +Y face
                     2, 3, 7, 6, // -X face
                     3, 0, 4, 7  // -Y face
                   };
  unboundedMesh.fill( 8, ucoords, 6, QUADRILATERAL, uconn, 0, err );
  CPPUNIT_ASSERT( !err );
  unboundedMesh.set_domain( &unboundedDomain );
  
    // Define a mesh on a cylinder with a radius of
    // one that is capped at z = +/- 2.  Define the
    // mesh as the 8 quads defining the sides of a pair of cubes
    // stacked axially in the cylinder
  const double V = 1.0 / sqrt(2.0);
  double bcoords[] = {  V, -V, -2,
                        V,  V, -2,
                       -V,  V, -2,
                       -V, -V, -2,
                       
                        V, -V,  0,
                        V,  V,  0,
                       -V,  V,  0,
                       -V, -V,  0,
                       
                        V, -V,  2,
                        V,  V,  2,
                       -V,  V,  2,
                       -V, -V,  2};
  size_t bconn[] = { // lower cube side faces
                      0,  1,  5,  4, // +X face
                      1,  2,  6,  5, // +Y face
                      2,  3,  7,  6, // -X face
                      3,  0,  4,  7, // -Y face
                     // upper cube side faces
                      4,  5,  9,  8, // +X face
                      5,  6, 10,  9, // +Y face
                      6,  7, 11, 10, // -X face
                      7,  4,  8, 11, // -Y face
                    };
  boundedMesh.fill( 12, bcoords, 8, QUADRILATERAL, bconn, 0, err );
  CPPUNIT_ASSERT( !err );
  boundedMesh.set_domain( &boundedDomain );
  
    // set element and vertex handles arrays
  size_t i = 0;
  for (i = 0; i < 12; ++i)
    boundedMesh.get_vertex_handles_array()[i] = (Mesh::VertexHandle)i;
  for (i = 0; i < 8; ++i)
    boundedMesh.get_element_handles_array()[i] = (Mesh::ElementHandle)i;
  
    // Bound the unit cylinder at +/- 1 on the z axis  
  std::vector<Mesh::VertexHandle> upper_curve(4), lower_curve(4);
  for (i = 0; i < 4; ++i)
  {
    lower_curve[i] = (Mesh::VertexHandle)i;
    upper_curve[i] = (Mesh::VertexHandle)(i+8);
  }
  boundedDomain.create_curve( -2, lower_curve );
  boundedDomain.create_curve(  2, upper_curve );
}
/*
void PatchDataTestNormals::test_get_vertex_normals_infinite_domain()
{
  MsqPrintError err(cout);
  
  for (size_t i = 0; i < unboundedMesh.num_nodes(); ++i)
  {
    Vector3D pos = unboundedMesh.vertex_by_index( i );
    Vector3D norm;
    unboundedMesh.get_domain_normal_at_vertex( i, false, norm, err );
    CPPUNIT_ASSERT(!err);
      // all points are on unit sphere centered at origin, so
      // the normal should be the same as the point.
    ASSERT_VECTORS_EQUAL( pos, norm );
  }
}

void PatchDataTestNormals::test_get_vertex_normals_bounded_domain()
{
  MsqPrintError err(cout);
  Vector3D norm;
  
    // Vertices 0 to 3 and 8 to 11 should lie on the end
    // curves of the cylinder.  There is no single valid normal
    // for a vertex on a curve, so it should fail for each
    // of these.  
  size_t indices[] = { 0, 1, 2, 3, 8, 9, 10, 11 };
  for (size_t i = 0; i < (sizeof(indices)/sizeof(indices[0])); ++i)
  {
    boundedMesh.get_domain_normal_at_vertex( indices[i], false, norm, err );
    CPPUNIT_ASSERT( err );
    err.clear();
  }
    
    // The remaining for vertices lie in the Z plane and the
    // cylinder's axis is the Z axis, so the normal should
    // be the same as the point coordinates.
  for (size_t j = 4; j < 8; ++j)
  {
    Vector3D pos = unboundedMesh.vertex_by_index( j );
    unboundedMesh.get_domain_normal_at_vertex( j, false, norm, err );
    CPPUNIT_ASSERT(!err);
      // all points are on unit sphere centered at origin, so
      // the normal should be the same as the point.
    ASSERT_VECTORS_EQUAL( pos, norm );
  }
}
*/ 
    
void PatchDataTestNormals::test_get_corner_normals_infinite_domain()
{
  MsqPrintError err(cout);
  
    // Element 0 is a quad parallel to and below the Z plane.
    // All corners of the element lie on the unit sphere and
    // thus the normal should be the same as the location.
  const size_t elem_index = 0;
  std::vector<Vector3D> coords;
  Vector3D normals[4];
  unboundedMesh.get_element_vertex_coordinates( elem_index, coords, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT( coords.size() == 4 /*quad*/ );
  unboundedMesh.get_domain_normals_at_corners( elem_index, normals, err );
  CPPUNIT_ASSERT( !err );
  for (size_t i = 0; i < 4; ++i)
  {
    ASSERT_VECTORS_EQUAL( coords[i], normals[i] );
  }
}


void PatchDataTestNormals::test_get_corner_normals_bounded_domain()
{
  MsqPrintError err(cout);
  std::vector<Vector3D> coords;
  Vector3D normals[4];
  
    // Element 0 is a quad in the plane X=1/sqrt(2).  Two of
    // the vertices of this element lie on the lower bounding
    // curve of the cylinder, and the other two lie in the
    // Z=0 plane
  const size_t elem_index = 0;
  boundedMesh.get_element_vertex_coordinates( elem_index, coords, err );
  CPPUNIT_ASSERT( !err );
  CPPUNIT_ASSERT( coords.size() == 4 /*quad*/ );
  boundedMesh.get_domain_normals_at_corners( elem_index, normals, err );
  CPPUNIT_ASSERT( !err );
  for (size_t i = 0; i < 4; ++i)
  {
    coords[i][2] = 0; // project into Z plane
    ASSERT_VECTORS_EQUAL( coords[i], normals[i] );
  }
} 


void PatchDataTestNormals::test_get_element_normals_infinite_domain()
{
  MsqPrintError err(cout);
  Vector3D expected_normals[] = { Vector3D(  0,  0, -1),
                                  Vector3D(  0,  0,  1),
                                  Vector3D(  1,  0,  0),
                                  Vector3D(  0,  1,  0),
                                  Vector3D( -1,  0,  0),
                                  Vector3D(  0, -1,  0) };
                                   
  CPPUNIT_ASSERT( unboundedMesh.num_elements() == 6u );
  for (size_t i = 0; i < 6u; ++i)
  {
    Vector3D norm;
    unboundedMesh.get_domain_normal_at_element( i, norm, err );
    CPPUNIT_ASSERT(!err);
    ASSERT_VECTORS_EQUAL( expected_normals[i], norm );
  }    
}

void PatchDataTestNormals::test_get_element_normals_bounded_domain()
{
  MsqPrintError err(cout);
  Vector3D expected_normals[] = { Vector3D(  1,  0,  0),
                                  Vector3D(  0,  1,  0),
                                  Vector3D( -1,  0,  0),
                                  Vector3D(  0, -1,  0),
                                  Vector3D(  1,  0,  0),
                                  Vector3D(  0,  1,  0),
                                  Vector3D( -1,  0,  0),
                                  Vector3D(  0, -1,  0)
                                };
                                   
  CPPUNIT_ASSERT( boundedMesh.num_elements() == 8u );
  for (size_t i = 0; i < 8u; ++i)
  {
    Vector3D norm;
    boundedMesh.get_domain_normal_at_element( i, norm, err );
    CPPUNIT_ASSERT(!err);
    ASSERT_VECTORS_EQUAL( expected_normals[i], norm );
  }    
}
  


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PatchDataTestNormals, "PatchDataTestNormals");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PatchDataTestNormals, "Unit");
