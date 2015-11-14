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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
 
    (2010) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file HexLagrangeShapeTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Mesquite_HexLagrangeShape.hpp"
#include "Mesquite_TopologyInfo.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_IdealElements.hpp"
#include "Mesquite_JacobianCalculator.hpp"

#include "UnitUtil.hpp"

#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>

using namespace Mesquite;
using namespace std;
const double epsilon = 1e-6;
#define ASSERT_VALUES_EQUAL( v1, v2, location ) \
  ASSERT_MESSAGE( value_message( (location), (v1), (v2) ), test_value((v1), (v2)) )

static bool test_value( double v1, double v2 )
  { return fabs(v1 - v2) < epsilon; }

static inline CppUnit::Message value_message( unsigned location, double v1, double v2 )
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
  if (location < 9) 
    buffer3 << "Corner " << location;
  else if (location < 20)
    buffer3 << "Edge " << location-8;
  else if (location < 26)
    buffer3 << "Face " << location-20;
  else if (location == 26)
    buffer3 << "Mid-element";
  else
    buffer3 << "INVALID!!";
  m.addDetail( buffer3.str() );
  return m;
}

static bool test_value( MsqVector<3> v1, MsqVector<3> v2 )
  { return length(v1 - v2) < epsilon; }

static inline CppUnit::Message value_message( unsigned location, 
                                              MsqVector<3> v1, 
                                              MsqVector<3> v2 )
{
  CppUnit::Message m( "equality assertion failed" );

  std::ostringstream buffer1;
  buffer1 << "Expected : [" << v1[0] << ", " << v1[1] << ", " << v1[2] << "]";
   m.addDetail( buffer1.str() );

  std::ostringstream buffer2;
  buffer2 << "Actual   : [" << v2[0] << ", " << v2[1] << ", " << v2[2] << "]";
  m.addDetail( buffer2.str() );

  std::ostringstream buffer3;
  buffer3 << "Location : ";
  if (location < 9) 
    buffer3 << "Corner " << location;
  else if (location < 20)
    buffer3 << "Edge " << location-8;
  else if (location < 26)
    buffer3 << "Face " << location-20;
  else if (location == 26)
    buffer3 << "Mid-element";
  else
    buffer3 << "INVALID!!";
  m.addDetail( buffer3.str() );
  return m;
}

class HexLagrangeShapeTest : public CppUnit::TestFixture
{
  private:
    CPPUNIT_TEST_SUITE(HexLagrangeShapeTest);

    CPPUNIT_TEST(test_corners_coeff);
    CPPUNIT_TEST(test_edges_coeff);
    CPPUNIT_TEST(test_faces_coeff);
    CPPUNIT_TEST(test_mid_coeff);

    CPPUNIT_TEST(test_corners_derivs);
    CPPUNIT_TEST(test_edges_derivs);
    CPPUNIT_TEST(test_faces_derivs);
    CPPUNIT_TEST(test_mid_derivs);
    
    CPPUNIT_TEST(test_ideal_jacobian);
    
    CPPUNIT_TEST_SUITE_END();
  
    HexLagrangeShape sf;
    
  public:
    
    void test_corner_coeff( int corner );
    void test_edge_coeff( int edge );
    void test_face_coeff( int face );
    void test_mid_coeff();
    
    void test_corner_derivs( int corner );
    void test_edge_derivs( int edge );
    void test_face_derivs( int face );
    void test_mid_derivs();

    NodeSet allNodes;
    void setUp() { allNodes.set_all_nodes( HEXAHEDRON ); }

    void test_corners_coeff();
    void test_edges_coeff();
    void test_faces_coeff();

    void test_corners_derivs();
    void test_edges_derivs();
    void test_faces_derivs();
    
    void test_ideal_jacobian();
};



CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(HexLagrangeShapeTest, "HexLagrangeShapeTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(HexLagrangeShapeTest, "Unit");

// Indices
enum { XI = 0, ETA = 1, ZETA = 2 };
// Parameter range
const double min_xi =  0;
const double max_xi =  1;
const double mid_xi = 0.5*(min_xi + max_xi);

// Some lagrange polynomial terms
static double l21( double xi ) { return (2 * xi - 1) * (xi - 1); }
static double l22( double xi ) { return 4 * xi * (1 - xi); }
static double l23( double xi ) { return xi * (2 * xi - 1); }
// First derivatives of lagrange polynomial terms
static double dl21( double xi ) { return 4 * xi - 3; }
static double dl22( double xi ) { return  4 - 8 * xi; }
static double dl23( double xi ) { return  4 * xi - 1; }
// Indexible polynomail terms
typedef double (*f_t)(double);
f_t l[4] = { 0, &l21, &l22, &l23 };
f_t dl[4] = { 0, &dl21, &dl22, &dl23 };
const int C[][3] = {
   { 1, 1, 1 }, //   0 
   { 3, 1, 1 }, //   1 
   { 3, 3, 1 }, //   2 
   { 1, 3, 1 }, //   3 
   
   { 1, 1, 3 }, //   4 
   { 3, 1, 3 }, //   5 
   { 3, 3, 3 }, //   6 
   { 1, 3, 3 }, //   7 
   
   { 2, 1, 1 }, //   8 
   { 3, 2, 1 }, //   9 
   { 2, 3, 1 }, //  10
   { 1, 2, 1 }, //  11
   
   { 1, 1, 2 }, //  12
   { 3, 1, 2 }, //  13
   { 3, 3, 2 }, //  14
   { 1, 3, 2 }, //  15
   
   { 2, 1, 3 }, //  16
   { 3, 2, 3 }, //  17
   { 2, 3, 3 }, //  18
   { 1, 2, 3 }, //  19
   
   { 2, 1, 2 }, //  20
   { 3, 2, 2 }, //  21
   { 2, 3, 2 }, //  22
   { 1, 2, 2 }, //  23
   
   { 2, 2, 1 }, //  24
   { 2, 2, 3 }, //  25
   
   { 2, 2, 2 }  //  26
};

const double corners[8][3] = {
  { min_xi, min_xi, min_xi },
  { max_xi, min_xi, min_xi },
  { max_xi, max_xi, min_xi },
  { min_xi, max_xi, min_xi },
  { min_xi, min_xi, max_xi },
  { max_xi, min_xi, max_xi },
  { max_xi, max_xi, max_xi },
  { min_xi, max_xi, max_xi } };
  
static MsqVector<3> XI_corner( int i )
  { return MsqVector<3>(corners[i]); }

static MsqVector<3> XI_edge( int i )
  { 
    const unsigned* idx = TopologyInfo::edge_vertices( HEXAHEDRON, i );
    return 0.5 * (XI_corner(idx[0]) + XI_corner(idx[1]));
  }

static MsqVector<3> XI_face( int i )
  { 
    unsigned n;
    const unsigned* idx = TopologyInfo::face_vertices( HEXAHEDRON, i, n );
    MsqVector<3> result = XI_corner(idx[0]);
    for (unsigned j = 1; j < n; ++j)
      result += XI_corner(idx[j]);
    result *= 1.0/n;
    return result;
  }

static MsqVector<3> XI_elem()
  { 
    const double vals[] = { mid_xi, mid_xi, mid_xi };
    return MsqVector<3>(vals);
  }

static double N( int i, MsqVector<3> xi ) 
{
  const int* c = C[i];
  return l[c[XI]](xi[XI]) * l[c[ETA]](xi[ETA]) * l[c[ZETA]](xi[ZETA]);
}

static MsqVector<3> dN( int i, MsqVector<3> xi ) 
{
  MsqVector<3> result;
  const int* c = C[i];
  result[XI  ] = dl[c[XI]](xi[XI]) *  l[c[ETA]](xi[ETA]) *  l[c[ZETA]](xi[ZETA]);
  result[ETA ] =  l[c[XI]](xi[XI]) * dl[c[ETA]](xi[ETA]) *  l[c[ZETA]](xi[ZETA]);
  result[ZETA] =  l[c[XI]](xi[XI]) *  l[c[ETA]](xi[ETA]) * dl[c[ZETA]](xi[ZETA]);
  return result;
}

static void get_coeffs( MsqVector<3> xi, double coeffs_out[27] )
{
  for (int i = 0; i < 27; ++i)
    coeffs_out[i] = N(i, xi);
}

static void get_derivs( MsqVector<3> xi, MsqVector<3> derivs_out[27] )
{
  for (int i = 0; i < 27; ++i)
    derivs_out[i] = dN(i, xi);
}



static void check_valid_indices( const size_t* vertices, size_t num_vtx )
{
    // check valid size of list (at least fout, at most all nodes)
  CPPUNIT_ASSERT( num_vtx <= 27 );
  CPPUNIT_ASSERT( num_vtx >= 4 );
    // make sure vertex indices are valid (in [0,26])
  size_t vertcopy[27];
  std::copy( vertices, vertices+num_vtx, vertcopy );
  std::sort( vertcopy, vertcopy+num_vtx );
  CPPUNIT_ASSERT( vertcopy[num_vtx-1] <= 26 ); // max value less than 27
    // make sure there are no duplicates in the list
  const size_t* iter = std::unique( vertcopy, vertcopy+num_vtx );
  CPPUNIT_ASSERT( iter == vertcopy+num_vtx );
}

static void check_no_zeros( const MsqVector<3>* derivs, size_t num_vtx )
{
  for (unsigned i = 0; i < num_vtx; ++i) {
    CPPUNIT_ASSERT( length(derivs[i]) > 1e-6 );
  }
}

static void compare_coefficients( const double* coeffs,
                                  const size_t* indices,
                                  size_t num_coeff,
                                  const double* expected_coeffs,
                                  unsigned loc )
{
    // find the location in the returned list for each node
  size_t revidx[27];
  double test_vals[27];
  for (size_t i = 0; i < 27; ++i) {
    revidx[i] = std::find( indices, indices+num_coeff, i ) - indices;
    test_vals[i] = (revidx[i] == num_coeff) ? 0.0 : coeffs[revidx[i]];
  }
    
    // compare expected and actual coefficient values
  ASSERT_VALUES_EQUAL( expected_coeffs[0], test_vals[0], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[1], test_vals[1], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[2], test_vals[2], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[3], test_vals[3], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[4], test_vals[4], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[5], test_vals[5], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[6], test_vals[6], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[7], test_vals[7], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[8], test_vals[8], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[9], test_vals[9], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[10], test_vals[10], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[11], test_vals[11], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[12], test_vals[12], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[13], test_vals[13], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[14], test_vals[14], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[15], test_vals[15], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[16], test_vals[16], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[17], test_vals[17], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[18], test_vals[18], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[19], test_vals[19], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[20], test_vals[20], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[21], test_vals[21], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[22], test_vals[22], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[23], test_vals[23], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[24], test_vals[24], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[25], test_vals[25], loc );
  ASSERT_VALUES_EQUAL( expected_coeffs[26], test_vals[26], loc );
}

static void compare_derivatives( const size_t* vertices,
                                 size_t num_vtx,
                                 const MsqVector<3>* actual,
                                 const MsqVector<3>* expected,
                                 unsigned loc )
{
  check_valid_indices( vertices, num_vtx );
  check_no_zeros( actual, num_vtx );
  
    // Input has values in dxi & deta & dzeta only for nodes in 'vertices'
    // Convert to values for every possible node, with zero's for
    // nodes that are not present.
  CPPUNIT_ASSERT(num_vtx <= 9);
  MsqVector<3> expanded[27];
  std::fill( expanded, expanded+27,MsqVector<3>(0.0) );
  for (unsigned i = 0; i < num_vtx; ++i) {
    CPPUNIT_ASSERT(vertices[i] <= 26);
    expanded[vertices[i]] = actual[i];
  }
  
  ASSERT_VALUES_EQUAL( expected[0], expanded[0], loc );
  ASSERT_VALUES_EQUAL( expected[1], expanded[1], loc );
  ASSERT_VALUES_EQUAL( expected[2], expanded[2], loc );
  ASSERT_VALUES_EQUAL( expected[3], expanded[3], loc );
  ASSERT_VALUES_EQUAL( expected[4], expanded[4], loc );
  ASSERT_VALUES_EQUAL( expected[5], expanded[5], loc );
  ASSERT_VALUES_EQUAL( expected[6], expanded[6], loc );
  ASSERT_VALUES_EQUAL( expected[7], expanded[7], loc );
  ASSERT_VALUES_EQUAL( expected[8], expanded[8], loc );
  ASSERT_VALUES_EQUAL( expected[9], expanded[9], loc );
  ASSERT_VALUES_EQUAL( expected[10], expanded[10], loc );
  ASSERT_VALUES_EQUAL( expected[11], expanded[11], loc );
  ASSERT_VALUES_EQUAL( expected[12], expanded[12], loc );
  ASSERT_VALUES_EQUAL( expected[13], expanded[13], loc );
  ASSERT_VALUES_EQUAL( expected[14], expanded[14], loc );
  ASSERT_VALUES_EQUAL( expected[15], expanded[15], loc );
  ASSERT_VALUES_EQUAL( expected[16], expanded[16], loc );
  ASSERT_VALUES_EQUAL( expected[17], expanded[17], loc );
  ASSERT_VALUES_EQUAL( expected[18], expanded[18], loc );
  ASSERT_VALUES_EQUAL( expected[19], expanded[19], loc );
  ASSERT_VALUES_EQUAL( expected[20], expanded[20], loc );
  ASSERT_VALUES_EQUAL( expected[21], expanded[21], loc );
  ASSERT_VALUES_EQUAL( expected[22], expanded[22], loc );
  ASSERT_VALUES_EQUAL( expected[23], expanded[23], loc );
  ASSERT_VALUES_EQUAL( expected[24], expanded[24], loc );
  ASSERT_VALUES_EQUAL( expected[25], expanded[25], loc );
  ASSERT_VALUES_EQUAL( expected[26], expanded[26], loc );
}

void HexLagrangeShapeTest::test_corner_coeff( int corner )
{
  MsqPrintError err(std::cout);
  
  double expected[27];
  get_coeffs( XI_corner(corner), expected );
  
  double coeff[100];
  size_t num_coeff = 11, indices[100];
  sf.coefficients( Sample(0, corner), allNodes, coeff, indices, num_coeff, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, num_coeff, expected, corner );
}

void HexLagrangeShapeTest::test_edge_coeff( int edge )
{
  MsqPrintError err(std::cout);
  
  double expected[27];
  get_coeffs( XI_edge(edge), expected );
  
  double coeff[100];
  size_t num_coeff = 11, indices[100];
  sf.coefficients( Sample(1, edge), allNodes, coeff, indices, num_coeff, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, num_coeff, expected, edge+8 );
}

void HexLagrangeShapeTest::test_face_coeff( int face )
{
  MsqPrintError err(std::cout);
  
  double expected[27];
  get_coeffs( XI_face(face), expected );
  
  double coeff[100];
  size_t num_coeff = 11, indices[100];
  sf.coefficients( Sample(2, face), allNodes, coeff, indices, num_coeff, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, num_coeff, expected, face+20 );
}

void HexLagrangeShapeTest::test_mid_coeff()
{
  MsqPrintError err(std::cout);
  
  double expected[27];
  get_coeffs( XI_elem(), expected );
  
  double coeff[100];
  size_t num_coeff = 11, indices[100];
  sf.coefficients( Sample(3, 0), allNodes, coeff, indices, num_coeff, err );
  CPPUNIT_ASSERT( !err );
  
  compare_coefficients( coeff, indices, num_coeff, expected, 26 );
}

void HexLagrangeShapeTest::test_corner_derivs( int corner )
{
  MsqPrintError err(std::cout);
  
  MsqVector<3> expected[27];
  get_derivs( XI_corner(corner), expected );
  
  size_t vertices[100], num_vtx = 23;
  MsqVector<3> derivs[100];
  sf.derivatives( Sample(0, corner), allNodes, vertices, derivs, num_vtx, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, num_vtx, derivs, expected, corner );
}

void HexLagrangeShapeTest::test_edge_derivs( int edge )
{
  MsqPrintError err(std::cout);
  
  MsqVector<3> expected[27];
  get_derivs( XI_edge(edge), expected );
  
  size_t vertices[100], num_vtx = 23;
  MsqVector<3> derivs[100];
  sf.derivatives( Sample(1, edge), allNodes, vertices, derivs, num_vtx, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, num_vtx, derivs, expected, edge+8 );
}

void HexLagrangeShapeTest::test_face_derivs( int face )
{
  MsqPrintError err(std::cout);
  
  MsqVector<3> expected[27];
  get_derivs( XI_face(face), expected );
  
  size_t vertices[100], num_vtx = 23;
  MsqVector<3> derivs[100];
  sf.derivatives( Sample(2, face), allNodes, vertices, derivs, num_vtx, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, num_vtx, derivs, expected, face+20 );
}

void HexLagrangeShapeTest::test_mid_derivs( )
{
  MsqPrintError err(std::cout);
  
  MsqVector<3> expected[27];
  get_derivs( XI_elem(), expected );
  
  size_t vertices[100], num_vtx = 23;
  MsqVector<3> derivs[100];
  sf.derivatives( Sample(3, 0), allNodes, vertices, derivs, num_vtx, err );
  CPPUNIT_ASSERT( !err );
  
  compare_derivatives( vertices, num_vtx, derivs, expected, 26 );
}

void HexLagrangeShapeTest::test_corners_coeff()
{
  for (unsigned i = 0; i < 8; ++i) 
    test_corner_coeff( i );
}

void HexLagrangeShapeTest::test_edges_coeff()
{
  for (unsigned i = 0; i < 12; ++i) 
    test_edge_coeff( i );
}

void HexLagrangeShapeTest::test_faces_coeff()
{
  for (unsigned i = 0; i < 6; ++i) 
    test_face_coeff( i );
}

void HexLagrangeShapeTest::test_corners_derivs()
{
  for (unsigned i = 0; i < 8; ++i) 
    test_corner_derivs( i );
}

void HexLagrangeShapeTest::test_edges_derivs()
{
  for (unsigned i = 0; i < 12; ++i) 
    test_edge_derivs( i );
}

void HexLagrangeShapeTest::test_faces_derivs()
{
  for (unsigned i = 0; i < 6; ++i) 
    test_face_derivs( i );
}

void HexLagrangeShapeTest::test_ideal_jacobian()
{
  MsqError err;
  MsqMatrix<3,3> J_act, J_exp(1.0);
  sf.ideal( Sample(3,0), J_act, err );
  ASSERT_NO_ERROR(err);

    // Matrices should be a rotation of each other.
    // First, calculate tentative rotation matrix
  MsqMatrix<3,3> R = inverse(J_exp) * J_act;
    // next check that it is a rotation
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(R), 1e-6 ); // no scaling
  ASSERT_MATRICES_EQUAL( transpose(R), inverse(R), 1e-6 ); // orthogonal
}
