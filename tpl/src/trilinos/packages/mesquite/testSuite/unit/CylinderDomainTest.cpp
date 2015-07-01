#include "CylinderDomain.hpp"
#include "MsqError.hpp"
#include <cppunit/extensions/HelperMacros.h>

const double EPSILON = 1e-6;
#define ASSERT_VECTORS_EQUAL( A, B ) \
  CPPUNIT_ASSERT( ((A)-(B)).length() < EPSILON )

using namespace Mesquite;

class CylinderDomainTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE( CylinderDomainTest );
  CPPUNIT_TEST( test_z_basic );
  CPPUNIT_TEST( test_z_snap_to );
  CPPUNIT_TEST( test_z_normal_at );
  CPPUNIT_TEST( test_z_closest_point );
  CPPUNIT_TEST( test_x_basic );
  CPPUNIT_TEST( test_x_snap_to );
  CPPUNIT_TEST( test_x_normal_at );
  CPPUNIT_TEST( test_x_closest_point );
  CPPUNIT_TEST( test_domain_DoF );
  CPPUNIT_TEST_SUITE_END();

  CylinderDomain z, x;

public:

  CylinderDomainTest() 
    : z( 1 ), x( 2, Vector3D(1,0,0) ) {}

  void setUp() {}
  void tearDown() {}
  
  void test_z_basic();
  void test_z_snap_to();
  void test_z_normal_at();
  void test_z_closest_point();
  
  void test_x_basic();
  void test_x_snap_to();
  void test_x_normal_at();
  void test_x_closest_point();
  
  void test_domain_DoF();
};

void CylinderDomainTest::test_z_basic()
{
  CPPUNIT_ASSERT_EQUAL( Vector3D(0,0,1), z.axis()   );
  CPPUNIT_ASSERT_EQUAL(            1.0 , z.radius() );
  CPPUNIT_ASSERT_EQUAL( Vector3D(0,0,0), z.center() );
}

void CylinderDomainTest::test_x_basic()
{
  CPPUNIT_ASSERT_EQUAL( Vector3D(1,0,0), x.axis()   );
  CPPUNIT_ASSERT_EQUAL(            2.0 , x.radius() );
  CPPUNIT_ASSERT_EQUAL( Vector3D(0,0,0), x.center() );
}

void CylinderDomainTest::test_z_snap_to()
{
  Vector3D vect;
  
  vect.set( 0.5, 0, 1 );
  z.snap_to( 0, vect );
  ASSERT_VECTORS_EQUAL( Vector3D( z.radius(), 0, 1 ), vect );
  
  vect.set( 0, 100, -5 );
  z.snap_to( 0, vect );
  ASSERT_VECTORS_EQUAL( Vector3D( 0, z.radius(), -5 ), vect );
  
  vect = z.center();
  z.snap_to( 0, vect );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, vect.z(), EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( z.radius(), (vect - z.center()).length(), EPSILON );
}

void CylinderDomainTest::test_x_snap_to()
{
  Vector3D vect;
  
  vect.set( 1, 0, 0.5 );
  x.snap_to( 0, vect );
  ASSERT_VECTORS_EQUAL( Vector3D( 1, 0, x.radius() ), vect );
  
  vect.set( -5, 100, 0 );
  x.snap_to( 0, vect );
  ASSERT_VECTORS_EQUAL( Vector3D( -5, x.radius(), 0 ), vect );
  
  vect = x.center();
  x.snap_to( 0, vect );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, vect.x(), EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( x.radius(), (vect - x.center()).length(), EPSILON );
}

void CylinderDomainTest::test_z_normal_at()
{
  Vector3D vect;
  
  vect.set( 0.5, 0, 1 );
  z.vertex_normal_at( 0, vect );
  ASSERT_VECTORS_EQUAL( Vector3D( 1, 0, 0 ), vect );
  
  vect.set( 0, 100, -5 );
  z.vertex_normal_at( 0, vect );
  ASSERT_VECTORS_EQUAL( Vector3D( 0, 1, 0 ), vect );
  
  vect = z.center();
  z.vertex_normal_at( 0, vect );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, vect.z(), EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, vect.length(), EPSILON );
  
  vect.set( 0.5, 0, 1 );
  z.element_normal_at( 0, vect );
  ASSERT_VECTORS_EQUAL( Vector3D( 1, 0, 0 ), vect );
  
  vect.set( 0, 100, -5 );
  z.element_normal_at( 0, vect );
  ASSERT_VECTORS_EQUAL( Vector3D( 0, 1, 0 ), vect );
  
  vect = z.center();
  z.element_normal_at( 0, vect );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, vect.z(), EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, vect.length(), EPSILON );
}

void CylinderDomainTest::test_x_normal_at()
{
  Vector3D vect;
  
  vect.set( 1, 0, 0.5 );
  x.vertex_normal_at( 0, vect );
  ASSERT_VECTORS_EQUAL( Vector3D( 0, 0, 1 ), vect );
  
  vect.set( -5, 100, 0 );
  x.vertex_normal_at( 0, vect );
  ASSERT_VECTORS_EQUAL( Vector3D( 0, 1, 0 ), vect );
  
  vect = x.center();
  x.vertex_normal_at( 0, vect );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, vect.x(), EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, vect.length(), EPSILON );
  
  vect.set( 1, 0, 0.5 );
  x.element_normal_at( 0, vect );
  ASSERT_VECTORS_EQUAL( Vector3D( 0, 0, 1 ), vect );
  
  vect.set( -5, 100, 0 );
  x.element_normal_at( 0, vect );
  ASSERT_VECTORS_EQUAL( Vector3D( 0, 1, 0 ), vect );
  
  vect = x.center();
  x.element_normal_at( 0, vect );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, vect.x(), EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, vect.length(), EPSILON );
}

void CylinderDomainTest::test_z_closest_point()
{
  MsqError err;
  Vector3D point( z.radius(), 0, 20 );
  Vector3D close, normal;
  z.closest_point( 0, point, close, normal, err );
  CPPUNIT_ASSERT(!err);
  ASSERT_VECTORS_EQUAL( point, close );
  ASSERT_VECTORS_EQUAL( Vector3D(1, 0, 0), normal );
}

void CylinderDomainTest::test_x_closest_point()
{
  MsqError err;
  Vector3D point( 20, 0, x.radius() );
  Vector3D close, normal;
  x.closest_point( 0, point, close, normal, err );
  CPPUNIT_ASSERT(!err);
  ASSERT_VECTORS_EQUAL( point, close );
  ASSERT_VECTORS_EQUAL( Vector3D(0, 0, 1), normal );
}

void CylinderDomainTest::test_domain_DoF()
{
  const size_t count = 3;
  unsigned short dof_vals[count];
  MsqError err;
  z.domain_DoF( 0, dof_vals, count, err );
  CPPUNIT_ASSERT(!err);
  for (size_t i = 0; i < count; ++i)
    CPPUNIT_ASSERT_EQUAL( (unsigned short)2, dof_vals[i] );
}

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( CylinderDomainTest, "CylinderDomainTest" );
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( CylinderDomainTest, "Unit" );
  
