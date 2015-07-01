/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retian certain rights to this software.

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


/** \file TetLagrangeShapeTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "TetLagrangeShape.hpp"
#include "TopologyInfo.hpp"
#include "MsqError.hpp"
#include "IdealElements.hpp"
#include "JacobianCalculator.hpp"
#include <math.h>

#include "UnitUtil.hpp"

#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>

using namespace Mesquite;
using namespace std;

const double epsilon = 1e-6;
#define ASSERT_VALUES_EQUAL( v1, v2, location, bits ) \
  ASSERT_MESSAGE( value_message( (location), (bits), (v1), (v2) ), \
                          (fabs((v1) - (v2)) < epsilon) )

static inline CppUnit::Message value_message( unsigned location, NodeSet bits, double v1, double v2 )
{
  CppUnit::Message m( "equality assertion failed" );

  std::ostringstream buffer1;
  buffer1 << "Expected : " << v1;
   m.addDetail( buffer1.str() );

  std::ostringstream buffer2;
  buffer2 << "Actual   : " << v2;
  m.addDetail( buffer2.str() );

  std::ostringstream buffer3;
  buffer3 << "Location : ";
  if (location < 4) 
    buffer3 << "Corner " << location;
  else if (location < 10)
    buffer3 << "Edge " << location-4;
  else if (location < 14)
    buffer3 << "Face " << location-10;
  else if (location == 14)
    buffer3 << "Mid-element";
  else
    buffer3 << "INVALID!!";
  m.addDetail( buffer3.str() );

  std::ostringstream buffer4;
  buffer4 << "Node Bits: " << bits;
  m.addDetail( buffer4.str() );
  return m;
}

class TetLagrangeShapeTest : public CppUnit::TestFixture
{
  private:
    CPPUNIT_TEST_SUITE(TetLagrangeShapeTest);

    CPPUNIT_TEST(test_coeff_corners);
    CPPUNIT_TEST(test_coeff_edges);
    CPPUNIT_TEST(test_coeff_faces);
    CPPUNIT_TEST(test_coeff_center);

    CPPUNIT_TEST(test_deriv_corners);
    CPPUNIT_TEST(test_deriv_edges);
    CPPUNIT_TEST(test_deriv_faces);
    CPPUNIT_TEST(test_deriv_center);

    CPPUNIT_TEST(test_mid_elem_node_coeff);
    CPPUNIT_TEST(test_mid_elem_node_deriv);
    CPPUNIT_TEST(test_mid_face_node_coeff);
    CPPUNIT_TEST(test_mid_face_node_deriv);
    
    CPPUNIT_TEST(test_ideal_jacobian);
    
    CPPUNIT_TEST_SUITE_END();
  
    TetLagrangeShape sf;
    
    void test_corner_coeff( int corner, NodeSet nodeset );
    void test_edge_coeff( int edge, NodeSet nodeset );
    void test_face_coeff( int face, NodeSet nodeset );
    void test_mid_coeff( NodeSet nodebits );
    
    void test_corner_derivs( int corner, NodeSet nodeset );
    void test_edge_derivs( int edge, NodeSet nodeset );
    void test_face_derivs( int face, NodeSet nodeset );
    void test_mid_derivs( NodeSet nodeset );
    
    void test_invalid_nodebits_coeff( NodeSet nodeset );
    void test_invalid_nodebits_deriv( NodeSet nodeset );
    
  public:

    void test_coeff_corners();
    void test_coeff_edges();
    void test_coeff_faces();
    void test_coeff_center();

    void test_deriv_corners();
    void test_deriv_edges();
    void test_deriv_faces();
    void test_deriv_center();
    
    void test_mid_elem_node_coeff();
    void test_mid_elem_node_deriv();
    void test_mid_face_node_coeff();
    void test_mid_face_node_deriv();
    
    void test_ideal_jacobian();
};



CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TetLagrangeShapeTest, "TetLagrangeShapeTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TetLagrangeShapeTest, "Unit");

static double N0( double r, double s, double t ) { double u = 1 - r - s - t; return u*(2*u - 1); }
static double N1( double r, double  , double   ) { return r*(2*r - 1); }
static double N2( double  , double s, double   ) { return s*(2*s - 1); }
static double N3( double  , double  , double t ) { return t*(2*t - 1);  }
static double N4( double r, double s, double t ) { double u = 1 - r - s - t; return 4*r*u; }
static double N5( double r, double s, double   ) { return 4*r*s; }
static double N6( double r, double s, double t ) { double u = 1 - r - s - t; return 4*s*u; }
static double N7( double r, double s, double t ) { double u = 1 - r - s - t; return 4*u*t; }
static double N8( double r, double  , double t ) { return 4*r*t; }
static double N9( double  , double s, double t ) { return 4*s*t; }

static double dN0dr( double r, double s, double t ) { double u = 1 - r - s - t; return 1 - 4*u; }
static double dN0ds( double r, double s, double t ) { double u = 1 - r - s - t; return 1 - 4*u; }
static double dN0dt( double r, double s, double t ) { double u = 1 - r - s - t; return 1 - 4*u; }

static double dN1dr( double r, double  , double   ) { return 4*r - 1; }
static double dN1ds( double  , double  , double   ) { return 0; }
static double dN1dt( double  , double  , double   ) { return 0; }

static double dN2dr( double  , double  , double   ) { return 0; }
static double dN2ds( double  , double s, double   ) { return 4*s - 1; }
static double dN2dt( double  , double  , double   ) { return 0; }

static double dN3dr( double  , double  , double   ) { return 0; }
static double dN3ds( double  , double  , double   ) { return 0; }
static double dN3dt( double  , double  , double t ) { return 4*t - 1; }

static double dN4dr( double r, double s, double t ) { double u = 1 - r - s - t; return 4*(u - r); }
static double dN4ds( double r, double  , double   ) { return -4*r; }
static double dN4dt( double r, double  , double   ) { return -4*r; }

static double dN5dr( double  , double s, double   ) { return 4*s; }
static double dN5ds( double r, double  , double   ) { return 4*r; }
static double dN5dt( double  , double  , double   ) { return 0; }

static double dN6dr( double  , double s, double   ) { return -4*s; }
static double dN6ds( double r, double s, double t ) { double u = 1 - r - s - t; return 4*(u - s); }
static double dN6dt( double  , double s, double   ) { return -4*s; }

static double dN7dr( double  , double  , double t ) { return -4*t; }
static double dN7ds( double  , double  , double t ) { return -4*t; }
static double dN7dt( double r, double s, double t ) { double u = 1 - r - s - t; return 4*(u - t); }

static double dN8dr( double  , double  , double t ) { return 4*t; }
static double dN8ds( double  , double  , double   ) { return 0; }
static double dN8dt( double r, double  , double   ) { return 4*r; }

static double dN9dr( double  , double  , double   ) { return 0; }
static double dN9ds( double  , double  , double t ) { return 4*t; }
static double dN9dt( double  , double s, double   ) { return 4*s; }

typedef double (*N_t)(double,double,double);
static const N_t N[] = { &N0, &N1, &N2, &N3, &N4, &N5, &N6, &N7, &N8, &N9 };
static const N_t dNdr[] = { &dN0dr, &dN1dr, &dN2dr, &dN3dr, &dN4dr, 
                     &dN5dr, &dN6dr, &dN7dr, &dN8dr, &dN9dr };
static const N_t dNds[] = { &dN0ds, &dN1ds, &dN2ds, &dN3ds, &dN4ds, 
                     &dN5ds, &dN6ds, &dN7ds, &dN8ds, &dN9ds };
static const N_t dNdt[] = { &dN0dt, &dN1dt, &dN2dt, &dN3dt, &dN4dt, 
                     &dN5dt, &dN6dt, &dN7dt, &dN8dt, &dN9dt };

static const double rst_corner[][3] = { {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
static const double rst_edge[][3] = { {0.5, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.0, 0.5, 0.0},
                                      {0.0, 0.0, 0.5}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5} };
static const double rst_face[][3] = { {1./3, 0.00, 1./3}, {1./3, 1./3, 1./3}, 
                                      {0.00, 1./3, 1./3}, {1./3, 1./3, 0.00} };
static const double rst_mid[3] = { 0.25, 0.25, 0.25 };

static unsigned edges[][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 },
                               { 0, 3 }, { 1, 3 }, { 2, 3 } };

static void get_coeff( NodeSet nodeset, const double* rst, double* coeffs )
{
  for (int i = 0; i < 10; ++i) 
    coeffs[i] = (*N[i])(rst[0], rst[1], rst[2]);
  for (int i = 0; i < 6; ++i) 
    if (!nodeset.mid_edge_node(i)) {
      coeffs[edges[i][0]] += 0.5 * coeffs[i+4];
      coeffs[edges[i][1]] += 0.5 * coeffs[i+4];
      coeffs[i+4] = 0;
    }
}

static void get_derivs( NodeSet nodeset, const double* rst, MsqVector<3>* derivs )
{
  for (int i = 0; i < 10; ++i) {
    derivs[i][0] = (*dNdr[i])(rst[0], rst[1], rst[2]);
    derivs[i][1] = (*dNds[i])(rst[0], rst[1], rst[2]);
    derivs[i][2] = (*dNdt[i])(rst[0], rst[1], rst[2]);
  }
  for (int i = 0; i < 6; ++i) 
    if (!nodeset.mid_edge_node(i)) {
      int j = edges[i][0];
      derivs[j][0] += 0.5 * derivs[i+4][0];
      derivs[j][1] += 0.5 * derivs[i+4][1];
      derivs[j][2] += 0.5 * derivs[i+4][2];
      j = edges[i][1];
      derivs[j][0] += 0.5 * derivs[i+4][0];
      derivs[j][1] += 0.5 * derivs[i+4][1];
      derivs[j][2] += 0.5 * derivs[i+4][2];
      derivs[i+4][0] = 0.0;
      derivs[i+4][1] = 0.0;
      derivs[i+4][2] = 0.0;
    }
}

static void check_valid_indices( const size_t* vertices, size_t num_vtx )
{
  CPPUNIT_ASSERT( num_vtx < 11 );
  CPPUNIT_ASSERT( num_vtx > 3 );
  size_t vtxcopy[10];
  std::copy( vertices, vertices+num_vtx, vtxcopy );
  std::sort( vtxcopy, vtxcopy+num_vtx );
  for (unsigned i = 1; i < num_vtx; ++i) {
    CPPUNIT_ASSERT( vtxcopy[i] != vtxcopy[i-1] );
    CPPUNIT_ASSERT( vtxcopy[i] < 10 );
  }
}

static void check_no_zeros( const MsqVector<3>* derivs, size_t num_vtx )
{
  for (unsigned i = 0; i < num_vtx; ++i) {
    double dr = derivs[i][0];
    double ds = derivs[i][1]; 
    double dt = derivs[i][2];
    CPPUNIT_ASSERT( (fabs(dr) > 1e-6) || (fabs(ds) > 1e-6) || (fabs(dt) > 1e-6) );
  }
}

static void compare_coefficients( const double* coeffs,
                                  const size_t* indices,
                                  const double* expected_coeffs,
                                  size_t num_coeff,
                                  unsigned loc, NodeSet bits )
{
    // find the location in the returned list for each node
  size_t revidx[10];
  double test_vals[10];
  for (size_t i = 0; i < 10; ++i) {
    revidx[i] = std::find( indices, indices+num_coeff, i ) - indices;
    test_vals[i] = (revidx[i] == num_coeff) ? 0.0 : coeffs[revidx[i]];
  }

    // Check that index list doesn't contain any nodes not actually
    // present in the element.
  CPPUNIT_ASSERT( bits.mid_edge_node(0) || (revidx[4] == num_coeff) );
  CPPUNIT_ASSERT( bits.mid_edge_node(1) || (revidx[5] == num_coeff) );
  CPPUNIT_ASSERT( bits.mid_edge_node(2) || (revidx[6] == num_coeff) );
  CPPUNIT_ASSERT( bits.mid_edge_node(3) || (revidx[7] == num_coeff) );
  CPPUNIT_ASSERT( bits.mid_edge_node(4) || (revidx[8] == num_coeff) );
  CPPUNIT_ASSERT( bits.mid_edge_node(5) || (revidx[9] == num_coeff) );
    
    // compare expected and actual coefficient values
  ASSERT_VALUES_EQUAL( expected_coeffs[0], test_vals[0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[1], test_vals[1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[2], test_vals[2], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[3], test_vals[3], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[4], test_vals[4], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[5], test_vals[5], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[6], test_vals[6], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[7], test_vals[7], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[8], test_vals[8], loc, bits );
  ASSERT_VALUES_EQUAL( expected_coeffs[9], test_vals[9], loc, bits );
}

static void compare_derivatives( const size_t* vertices,
                                 size_t num_vtx,
                                 const MsqVector<3>* derivs,
                                 const MsqVector<3>* expected_derivs,
                                 unsigned loc, NodeSet bits )
{
  check_valid_indices( vertices, num_vtx );
  check_no_zeros( derivs, num_vtx );
  MsqVector<3> expanded_derivs[30];
  memset( expanded_derivs, 0, sizeof(expanded_derivs) );
  for (unsigned i = 0; i < num_vtx; ++i) 
    expanded_derivs[vertices[i]] = derivs[i];
  
  ASSERT_VALUES_EQUAL( expected_derivs[0][0], expanded_derivs[0][0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[0][1], expanded_derivs[0][1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[0][2], expanded_derivs[0][2], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[1][0], expanded_derivs[1][0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[1][1], expanded_derivs[1][1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[1][2], expanded_derivs[1][2], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[2][0], expanded_derivs[2][0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[2][1], expanded_derivs[2][1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[2][2], expanded_derivs[2][2], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[3][0], expanded_derivs[3][0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[3][1], expanded_derivs[3][1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[3][2], expanded_derivs[3][2], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[4][0], expanded_derivs[4][0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[4][1], expanded_derivs[4][1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[4][2], expanded_derivs[4][2], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[5][0], expanded_derivs[5][0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[5][1], expanded_derivs[5][1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[5][2], expanded_derivs[5][2], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[6][0], expanded_derivs[6][0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[6][1], expanded_derivs[6][1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[6][2], expanded_derivs[6][2], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[7][0], expanded_derivs[7][0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[7][1], expanded_derivs[7][1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[7][2], expanded_derivs[7][2], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[8][0], expanded_derivs[8][0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[8][1], expanded_derivs[8][1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[8][2], expanded_derivs[8][2], loc, bits );
  
  ASSERT_VALUES_EQUAL( expected_derivs[9][0], expanded_derivs[9][0], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[9][1], expanded_derivs[9][1], loc, bits );
  ASSERT_VALUES_EQUAL( expected_derivs[9][2], expanded_derivs[9][2], loc, bits );
}

void TetLagrangeShapeTest::test_corner_coeff( int corner, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[10];
  get_coeff( nodebits, rst_corner[corner], expected );
  
  double coeff[100];
  size_t n = 29, indices[100];
  sf.coefficients( Sample(0, corner), nodebits, coeff, indices, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, expected, n, corner, nodebits );
}

void TetLagrangeShapeTest::test_edge_coeff( int edge, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[10];
  get_coeff( nodebits, rst_edge[edge], expected );
  
  double coeff[100];
  size_t n = 29, indices[100];
  sf.coefficients( Sample(1, edge), nodebits, coeff, indices, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, expected, n, edge+4, nodebits );
}

void TetLagrangeShapeTest::test_face_coeff( int face, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[10];
  get_coeff( nodebits, rst_face[face], expected );
  
  double coeff[100];
  size_t n = 29, indices[100];
  sf.coefficients( Sample(2, face), nodebits, coeff, indices, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, expected, n, face+10, nodebits );
}

void TetLagrangeShapeTest::test_mid_coeff( NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  double expected[10];
  get_coeff( nodebits, rst_mid, expected );
  
  double coeff[100];
  size_t n = 29, indices[100];
  sf.coefficients( Sample(3, 0), nodebits, coeff, indices, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, expected, n, 14, nodebits );
}

void TetLagrangeShapeTest::test_corner_derivs( int corner, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  MsqVector<3> expected[10];
  get_derivs( nodebits, rst_corner[corner], expected );
  
  MsqVector<3> derivs[100];
  size_t vertices[100], n = 29;
  sf.derivatives( Sample(0, corner), nodebits, vertices, derivs, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, n, derivs, expected, corner, nodebits );
}

void TetLagrangeShapeTest::test_edge_derivs( int edge, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  MsqVector<3> expected[10];
  get_derivs( nodebits, rst_edge[edge], expected );
  
  MsqVector<3> derivs[100];
  size_t vertices[100], n = 29;
  sf.derivatives( Sample(1, edge), nodebits, vertices, derivs, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, n, derivs, expected, edge+4, nodebits );
}

void TetLagrangeShapeTest::test_face_derivs( int face, NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  MsqVector<3> expected[10];
  get_derivs( nodebits, rst_face[face], expected );
  
  MsqVector<3> derivs[100];
  size_t vertices[100], n = 29;
  sf.derivatives( Sample(2, face), nodebits, vertices, derivs, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, n, derivs, expected, face+10, nodebits );
}

void TetLagrangeShapeTest::test_mid_derivs( NodeSet nodebits )
{
  MsqPrintError err(std::cout);
  
  MsqVector<3> expected[10];
  get_derivs( nodebits, rst_mid, expected );
  
  MsqVector<3> derivs[100];
  size_t vertices[100], n = 29;
  sf.derivatives( Sample(3, 0), nodebits, vertices, derivs, n, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, n, derivs, expected, 14, nodebits );
}


static NodeSet nodeset_from_bits( unsigned bits )
{
  NodeSet result;
  for (unsigned i = 0; i < 6; ++i)
    if (bits & (1<<i))
      result.set_mid_edge_node(i);
  for (unsigned i = 6; i < 10; ++i)
    if (bits & (1<<i))
      result.set_mid_face_node(i-6);
  if (bits & (1<<10))
    result.set_mid_region_node();
  return result;
}

void TetLagrangeShapeTest::test_coeff_corners()
{
  for (unsigned i = 0; i < 0x40; ++i) {
    test_corner_coeff( 0, nodeset_from_bits(i) );
    test_corner_coeff( 1, nodeset_from_bits(i) );
    test_corner_coeff( 2, nodeset_from_bits(i) );
    test_corner_coeff( 3, nodeset_from_bits(i) );
  }
}

void TetLagrangeShapeTest::test_coeff_edges()
{
  for (unsigned i = 0; i < 0x40; ++i) {
    test_edge_coeff( 0, nodeset_from_bits(i) );
    test_edge_coeff( 1, nodeset_from_bits(i) );
    test_edge_coeff( 2, nodeset_from_bits(i) );
    test_edge_coeff( 3, nodeset_from_bits(i) );
    test_edge_coeff( 4, nodeset_from_bits(i) );
    test_edge_coeff( 5, nodeset_from_bits(i) );
  }
}

void TetLagrangeShapeTest::test_coeff_faces()
{
  for (unsigned i = 0; i < 0x40; ++i) {
    test_face_coeff( 0, nodeset_from_bits(i) );
    test_face_coeff( 1, nodeset_from_bits(i) );
    test_face_coeff( 2, nodeset_from_bits(i) );
    test_face_coeff( 3, nodeset_from_bits(i) );
  }
}

void TetLagrangeShapeTest::test_coeff_center()
{
  for (unsigned i = 0; i < 0x40; ++i) {
    test_mid_coeff( nodeset_from_bits(i) );
  }
}

void TetLagrangeShapeTest::test_deriv_corners()
{
  for (unsigned i = 0; i < 0x40; ++i) {
    test_corner_derivs( 0, nodeset_from_bits(i) );
    test_corner_derivs( 1, nodeset_from_bits(i) );
    test_corner_derivs( 2, nodeset_from_bits(i) );
    test_corner_derivs( 3, nodeset_from_bits(i) );
  }
}

void TetLagrangeShapeTest::test_deriv_edges()
{
  for (unsigned i = 0; i < 0x40; ++i) {
    test_edge_derivs( 0, nodeset_from_bits(i) );
    test_edge_derivs( 1, nodeset_from_bits(i) );
    test_edge_derivs( 2, nodeset_from_bits(i) );
    test_edge_derivs( 3, nodeset_from_bits(i) );
    test_edge_derivs( 4, nodeset_from_bits(i) );
    test_edge_derivs( 5, nodeset_from_bits(i) );
  }
}

void TetLagrangeShapeTest::test_deriv_faces()
{
  for (unsigned i = 0; i < 0x40; ++i) {
    test_face_derivs( 0, nodeset_from_bits(i) );
    test_face_derivs( 1, nodeset_from_bits(i) );
    test_face_derivs( 2, nodeset_from_bits(i) );
    test_face_derivs( 3, nodeset_from_bits(i) );
  }
}

void TetLagrangeShapeTest::test_deriv_center()
{
  for (unsigned i = 0; i < 0x40; ++i) {
    test_mid_derivs( nodeset_from_bits(i) );
  }
}

void TetLagrangeShapeTest::test_invalid_nodebits_coeff( NodeSet bits )
{
  MsqError err;
  double coeff[100];
  size_t n, indices[100];
  
  sf.coefficients( Sample( 0, 0 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients( Sample( 0, 1 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients( Sample( 0, 2 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients( Sample( 0, 3 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  
  sf.coefficients( Sample( 1, 0 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients( Sample( 1, 1 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients( Sample( 1, 2 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients( Sample( 1, 3 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients( Sample( 1, 4 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients( Sample( 1, 5 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  
  sf.coefficients( Sample( 2, 0 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients( Sample( 2, 1 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients( Sample( 2, 2 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  sf.coefficients( Sample( 2, 3 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
  
  sf.coefficients( Sample( 3, 0 ), bits, coeff, indices, n, err );
  CPPUNIT_ASSERT( err );
}

void TetLagrangeShapeTest::test_invalid_nodebits_deriv( NodeSet bits )
{
  MsqError err;
  size_t verts[100], n;
  MsqVector<3> derivs[100];
  
  sf.derivatives( Sample( 0, 0 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives( Sample( 0, 1 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives( Sample( 0, 2 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives( Sample( 0, 3 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  
  sf.derivatives( Sample( 1, 0 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives( Sample( 1, 1 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives( Sample( 1, 2 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives( Sample( 1, 3 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives( Sample( 1, 4 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives( Sample( 1, 5 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  
  sf.derivatives( Sample( 2, 0 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives( Sample( 2, 1 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives( Sample( 2, 2 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  sf.derivatives( Sample( 2, 3 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
  
  sf.derivatives( Sample( 3, 0 ), bits, verts, derivs, n, err );
  CPPUNIT_ASSERT( err );
}


void TetLagrangeShapeTest::test_mid_elem_node_coeff()
{
  NodeSet nodeset;
  nodeset.set_mid_region_node();
  test_invalid_nodebits_coeff( nodeset );
}

void TetLagrangeShapeTest::test_mid_elem_node_deriv()
{
  NodeSet nodeset;
  nodeset.set_mid_region_node();
  test_invalid_nodebits_deriv( nodeset );
}

void TetLagrangeShapeTest::test_mid_face_node_coeff()
{
  NodeSet nodeset1, nodeset2;
  nodeset1.set_mid_face_node(0);
  test_invalid_nodebits_coeff( nodeset1 );
  nodeset2.set_mid_face_node(2);
  test_invalid_nodebits_coeff( nodeset2 );
}

void TetLagrangeShapeTest::test_mid_face_node_deriv()
{
  NodeSet nodeset1, nodeset2;
  nodeset1.set_mid_face_node(0);
  test_invalid_nodebits_deriv( nodeset1 );
  nodeset2.set_mid_face_node(2);
  test_invalid_nodebits_deriv( nodeset2 );
}


void TetLagrangeShapeTest::test_ideal_jacobian()
{
  MsqError err;
  MsqMatrix<3,3> J_act, J_exp;
  sf.ideal( Sample(3,0), J_act, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(J_act), 1e-6 );

  const Vector3D* verts = unit_edge_element( TETRAHEDRON );
  CPPUNIT_ASSERT(verts);
  
  JacobianCalculator jc;
  jc.get_Jacobian_3D( &sf, NodeSet(), Sample(2,0), verts, 4, J_exp, err );
  ASSERT_NO_ERROR(err);
  J_exp /= Mesquite::cbrt(det(J_exp));
  
    // Matrices should be a rotation of each other.
    // First, calculate tentative rotation matrix
  MsqMatrix<3,3> R = inverse(J_exp) * J_act;
    // next check that it is a rotation
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(R), 1e-6 ); // no scaling
  ASSERT_MATRICES_EQUAL( transpose(R), inverse(R), 1e-6 ); // orthogonal
}
