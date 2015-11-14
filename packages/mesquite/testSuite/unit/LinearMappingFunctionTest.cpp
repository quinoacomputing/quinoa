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


/** \file LinearMappingFunctionTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Mesquite_MappingFunction.hpp"
#include "Mesquite_LinearHexahedron.hpp"
#include "Mesquite_LinearQuadrilateral.hpp"
#include "Mesquite_LinearTetrahedron.hpp"
#include "Mesquite_LinearTriangle.hpp"
#include "Mesquite_LinearPrism.hpp"
#include "Mesquite_LinearPyramid.hpp"
#include "Mesquite_TopologyInfo.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_JacobianCalculator.hpp"
#include "Mesquite_IdealElements.hpp"

#include "UnitUtil.hpp"

#include <vector>
#include <algorithm>

using namespace Mesquite;
using namespace std;

class LinearMappingFunctionTest : public CppUnit::TestFixture
{
  private:
    CPPUNIT_TEST_SUITE(LinearMappingFunctionTest);

    CPPUNIT_TEST(test_linear_hex_coeff_corners);
    CPPUNIT_TEST(test_linear_hex_coeff_edges);
    CPPUNIT_TEST(test_linear_hex_coeff_faces);
    CPPUNIT_TEST(test_linear_hex_coeff_center);

    CPPUNIT_TEST(test_linear_quad_coeff_corners);
    CPPUNIT_TEST(test_linear_quad_coeff_edges);
//    CPPUNIT_TEST(test_linear_quad_coeff_faces);
    CPPUNIT_TEST(test_linear_quad_coeff_center);

    CPPUNIT_TEST(test_linear_tet_coeff_corners);
    CPPUNIT_TEST(test_linear_tet_coeff_edges);
    CPPUNIT_TEST(test_linear_tet_coeff_faces);
    CPPUNIT_TEST(test_linear_tet_coeff_center);

    CPPUNIT_TEST(test_linear_tri_coeff_corners);
    CPPUNIT_TEST(test_linear_tri_coeff_edges);
    CPPUNIT_TEST(test_linear_tri_coeff_faces);
    CPPUNIT_TEST(test_linear_tri_coeff_center);

    CPPUNIT_TEST(test_linear_prism_coeff_corners);
    CPPUNIT_TEST(test_linear_prism_coeff_edges);
    CPPUNIT_TEST(test_linear_prism_coeff_faces);
    CPPUNIT_TEST(test_linear_prism_coeff_center);

    CPPUNIT_TEST(test_linear_pyr_coeff_corners);
    CPPUNIT_TEST(test_linear_pyr_coeff_edges);
    CPPUNIT_TEST(test_linear_pyr_coeff_faces);
    CPPUNIT_TEST(test_linear_pyr_coeff_center);

    CPPUNIT_TEST(test_linear_hex_deriv_corners);
    CPPUNIT_TEST(test_linear_hex_deriv_edges);
    CPPUNIT_TEST(test_linear_hex_deriv_faces);
    CPPUNIT_TEST(test_linear_hex_deriv_center);

    CPPUNIT_TEST(test_linear_quad_deriv_corners);
    CPPUNIT_TEST(test_linear_quad_deriv_edges);
//    CPPUNIT_TEST(test_linear_quad_deriv_faces);
    CPPUNIT_TEST(test_linear_quad_deriv_center);

    CPPUNIT_TEST(test_linear_tet_deriv_corners);
    CPPUNIT_TEST(test_linear_tet_deriv_edges);
    CPPUNIT_TEST(test_linear_tet_deriv_faces);
    CPPUNIT_TEST(test_linear_tet_deriv_center);

    CPPUNIT_TEST(test_linear_tri_deriv_corners);
    CPPUNIT_TEST(test_linear_tri_deriv_edges);
//    CPPUNIT_TEST(test_linear_tri_deriv_faces);
    CPPUNIT_TEST(test_linear_tri_deriv_center);

    CPPUNIT_TEST(test_linear_prism_deriv_corners);
    CPPUNIT_TEST(test_linear_prism_deriv_edges);
    CPPUNIT_TEST(test_linear_prism_deriv_faces);
    CPPUNIT_TEST(test_linear_prism_deriv_center);

    CPPUNIT_TEST(test_linear_pyr_deriv_corners);
    CPPUNIT_TEST(test_linear_pyr_deriv_edges);
    CPPUNIT_TEST(test_linear_pyr_deriv_faces);
    CPPUNIT_TEST(test_linear_pyr_deriv_center);

    CPPUNIT_TEST(test_linear_hex_ideal);
    CPPUNIT_TEST(test_linear_quad_ideal);
    CPPUNIT_TEST(test_linear_tet_ideal);
    CPPUNIT_TEST(test_linear_tri_ideal);
    CPPUNIT_TEST(test_linear_prism_ideal);
   
    CPPUNIT_TEST_SUITE_END();
  
    LinearHexahedron hex;
    LinearQuadrilateral quad;
    LinearTetrahedron tet;
    LinearTriangle tri;
    LinearPrism prism;
    LinearPyramid pyr;
    
    static void hex_coeff( double xi[3], double coeff[8] );
    static void tet_coeff( double xi[3], double coeff[4] );
    static void quad_coeff( double xi[2], double coeff[4] );
    static void tri_coeff( double xi[2], double coeff[3] );
    static void prism_coeff( double xi[3], double coeff[6] );
    static void pyr_coeff( double xi[3], double coeff[5] );
    
    static void hex_deriv( double xi[3], double coeff_deriv[24] );
    static void tet_deriv( double xi[3], double coeff_deriv[12] );
    static void quad_deriv( double xi[2], double coeff_deriv[8] );
    static void tri_deriv( double xi[2], double coeff_deriv[6] );
    static void prism_deriv( double xi[3], double coeff_deriv[18] );
    static void pyr_deriv( double xi[3], double coeff_deriv[15] );
    
    typedef void (*map_func)(double*, double* );
    
    void do_coeff_test( MappingFunction& mf, 
                        unsigned subdim,
                        map_func mf2,
                        unsigned count,
                        double* xi );
    void do_deriv_test( MappingFunction2D& mf, 
                        unsigned subdim,
                        map_func mf2,
                        unsigned count,
                        double* xi );
    void do_deriv_test( MappingFunction3D& mf, 
                        unsigned subdim,
                        map_func mf2,
                        unsigned count,
                        double* xi );
    void do_ideal_test( MappingFunction2D& mf );
    void do_ideal_test( MappingFunction3D& mf );
                 
                 
    void test_coeff_fail( MappingFunction& mf, unsigned subdim );
    void test_deriv_fail( MappingFunction2D& mf, unsigned subdim );
    void test_deriv_fail( MappingFunction3D& mf, unsigned subdim );
    
    
    void xi_at_corners( EntityTopology type, double* xi, const int* corners );
    void xi_at_edges( EntityTopology type, double* xi, const int* corners );
    void xi_at_faces( EntityTopology type, double* xi, const int* corners );
    
  public:
  
    void setUp();
    void tearDown();
    

    void test_linear_hex_coeff_corners();
    void test_linear_hex_coeff_edges();
    void test_linear_hex_coeff_faces();
    void test_linear_hex_coeff_center();

    void test_linear_quad_coeff_corners();
    void test_linear_quad_coeff_edges();
//    void test_linear_quad_coeff_faces();
    void test_linear_quad_coeff_center();

    void test_linear_tet_coeff_corners();
    void test_linear_tet_coeff_edges();
    void test_linear_tet_coeff_faces();
    void test_linear_tet_coeff_center();

    void test_linear_tri_coeff_corners();
    void test_linear_tri_coeff_edges();
    void test_linear_tri_coeff_faces();
    void test_linear_tri_coeff_center();

    void test_linear_prism_coeff_corners();
    void test_linear_prism_coeff_edges();
    void test_linear_prism_coeff_faces();
    void test_linear_prism_coeff_center();

    void test_linear_pyr_coeff_corners();
    void test_linear_pyr_coeff_edges();
    void test_linear_pyr_coeff_faces();
    void test_linear_pyr_coeff_center();

    void test_linear_hex_deriv_corners();
    void test_linear_hex_deriv_edges();
    void test_linear_hex_deriv_faces();
    void test_linear_hex_deriv_center();

    void test_linear_quad_deriv_corners();
    void test_linear_quad_deriv_edges();
//    void test_linear_quad_deriv_faces();
    void test_linear_quad_deriv_center();

    void test_linear_tet_deriv_corners();
    void test_linear_tet_deriv_edges();
    void test_linear_tet_deriv_faces();
    void test_linear_tet_deriv_center();

    void test_linear_tri_deriv_corners();
    void test_linear_tri_deriv_edges();
//    void test_linear_tri_deriv_faces();
    void test_linear_tri_deriv_center();

    void test_linear_prism_deriv_corners();
    void test_linear_prism_deriv_edges();
    void test_linear_prism_deriv_faces();
    void test_linear_prism_deriv_center();

    void test_linear_pyr_deriv_corners();
    void test_linear_pyr_deriv_edges();
    void test_linear_pyr_deriv_faces();
    void test_linear_pyr_deriv_center();

    void test_linear_hex_ideal()   { do_ideal_test( hex   ); }
    void test_linear_quad_ideal()  { do_ideal_test( quad  ); }
    void test_linear_tet_ideal()   { do_ideal_test( tet   ); }
    void test_linear_tri_ideal()   { do_ideal_test( tri   ); }
    void test_linear_prism_ideal() { do_ideal_test( prism ); }
};



CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(LinearMappingFunctionTest, "LinearMappingFunctionTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(LinearMappingFunctionTest, "Unit");

void LinearMappingFunctionTest::setUp() {}
void LinearMappingFunctionTest::tearDown() {}



/*******************************************************************************
 *                             Xi values at element corners
 *******************************************************************************/

const int HexCorners[8][3] = { {  0,  0,  0 },
                               {  1,  0,  0 },
                               {  1,  1,  0 },
                               {  0,  1,  0 },
                               {  0,  0,  1 },
                               {  1,  0,  1 },
                               {  1,  1,  1 },
                               {  0,  1,  1 } };

const int QuadCorners[4][2] = { {  0,  0 },
                                {  1,  0 },
                                {  1,  1 },
                                {  0,  1 } };


const int TetCorners[12] = { 0, 0, 0, 
                             1, 0, 0, 
                             0, 1, 0,
                             0, 0, 1 };
  
const int TriCorners[6] = { 0, 0,
                            1, 0, 
                            0, 1 };

const int PrismCorners[18] = {  0, 0, 0,
                                0, 1, 0,
                                0, 0, 1,
                                1, 0, 0,
                                1, 1, 0,
                                1, 0, 1 };
                                
/*******************************************************************************
 * Functions to calculate element xi at different locations, given the xi
 * at the corners.
 *******************************************************************************/

void LinearMappingFunctionTest::xi_at_corners( EntityTopology type, double* xi, const int* corners )
{
  unsigned d = TopologyInfo::dimension( type );
  unsigned c = TopologyInfo::corners( type );
  for (unsigned i = 0; i < c*d; ++i)
    xi[i] = corners[i];
}

void LinearMappingFunctionTest::xi_at_edges( EntityTopology type, double* xi, const int* corners )
{
  MsqError err;
  unsigned d = TopologyInfo::dimension( type );
  unsigned e = TopologyInfo::edges( type );
  for (unsigned i = 0; i < e; ++i)
  {
    const unsigned* vtx = TopologyInfo::edge_vertices( type, i, err );
    CPPUNIT_ASSERT(!err);
    for (unsigned j = 0; j < d; ++j)
      xi[d*i+j] = (corners[d*vtx[0]+j] + corners[d*vtx[1]+j])/2.0;
  }
}


void LinearMappingFunctionTest::xi_at_faces( EntityTopology type, double* xi, const int* corners )
{
  MsqError err;
  unsigned d = TopologyInfo::dimension( type );
  unsigned f = TopologyInfo::faces( type );
  for (unsigned i = 0; i < f; ++i)
  {
    unsigned c;
    const unsigned* vtx = TopologyInfo::face_vertices( type, i, c, err );
    CPPUNIT_ASSERT(!err);
    for (unsigned j = 0; j < d; ++j)
    {
      int sum = 0;
      for (unsigned k = 0; k < c; ++k)
        sum += corners[d*vtx[k]+j];
      xi[d*i+j] = (double)sum / c;
    }
  }
}




/*******************************************************************************
 *                 Test functions for mapping function coefficients
 *******************************************************************************/

void LinearMappingFunctionTest::test_linear_hex_coeff_corners()
{ 
  double xi[24];
  xi_at_corners( HEXAHEDRON, xi, &HexCorners[0][0] );
  do_coeff_test( hex, 0, hex_coeff, 8, xi );
}
  
void LinearMappingFunctionTest::test_linear_hex_coeff_edges()
{
  double xi[36];
  xi_at_edges( HEXAHEDRON, xi, &HexCorners[0][0] );
  do_coeff_test( hex, 1, hex_coeff, 12, xi );
}

void LinearMappingFunctionTest::test_linear_hex_coeff_faces()
{
  double xi[18];
  xi_at_faces( HEXAHEDRON, xi, &HexCorners[0][0] );
  do_coeff_test( hex, 2, hex_coeff, 6, xi );
}

void LinearMappingFunctionTest::test_linear_hex_coeff_center()
{
  double xi[3] = { 0.5, 0.5, 0.5 };
  do_coeff_test( hex, 3, hex_coeff, 1, xi );
}

void LinearMappingFunctionTest::test_linear_quad_coeff_corners()
{ 
  double xi[8];
  xi_at_corners( QUADRILATERAL, xi, &QuadCorners[0][0] );
  do_coeff_test( quad, 0, quad_coeff, 4, xi );
}

void LinearMappingFunctionTest::test_linear_quad_coeff_edges()
{
  double xi[8];
  xi_at_edges( QUADRILATERAL, xi, &QuadCorners[0][0] );
  do_coeff_test( quad, 1, quad_coeff, 4, xi );
}

//void LinearMappingFunctionTest::test_linear_quad_coeff_faces()
//{
//  test_coeff_fail( quad, 3 );
//}

void LinearMappingFunctionTest::test_linear_quad_coeff_center()
{
  double xi[2] = { 0.5, 0.5 };
  do_coeff_test( quad, 2, quad_coeff, 1, xi );
}

void LinearMappingFunctionTest::test_linear_tet_coeff_corners()
{
  double xi[12];
  xi_at_corners( TETRAHEDRON, xi, TetCorners );
  do_coeff_test( tet, 0, tet_coeff, 4, xi );
}

void LinearMappingFunctionTest::test_linear_tet_coeff_edges()
{
  double xi[18];
  xi_at_edges( TETRAHEDRON, xi, TetCorners );
  do_coeff_test( tet, 1, tet_coeff, 6, xi );
}

void LinearMappingFunctionTest::test_linear_tet_coeff_faces()
{
  double xi[12];
  xi_at_faces( TETRAHEDRON, xi, TetCorners );
  do_coeff_test( tet, 2, tet_coeff, 4, xi );
}

void LinearMappingFunctionTest::test_linear_tet_coeff_center()
{
  double xi[3] = { 0.25, 0.25, 0.25 };
  do_coeff_test( tet, 3, tet_coeff, 1, xi );
}

void LinearMappingFunctionTest::test_linear_tri_coeff_corners()
{
  double xi[12];
  xi_at_corners( TRIANGLE, xi, TriCorners );
  do_coeff_test( tri, 0, tri_coeff, 3, xi );
}

void LinearMappingFunctionTest::test_linear_tri_coeff_edges()
{
  double xi[18];
  xi_at_edges( TRIANGLE, xi, TriCorners );
  do_coeff_test( tri, 1, tri_coeff, 3, xi );
}

void LinearMappingFunctionTest::test_linear_tri_coeff_faces()
{
  test_coeff_fail( tri, 3 );
}

void LinearMappingFunctionTest::test_linear_tri_coeff_center()
{
  double xi[2] = { 1./3, 1./3 };
  do_coeff_test( tri, 2, tri_coeff, 1, xi );
}

void LinearMappingFunctionTest::test_linear_prism_coeff_corners()
{
  double xi[18];
  xi_at_corners( PRISM, xi, PrismCorners );
  do_coeff_test( prism, 0, prism_coeff, 6, xi );
}

void LinearMappingFunctionTest::test_linear_prism_coeff_edges()
{
  double xi[27];
  xi_at_edges( PRISM, xi, PrismCorners );
  do_coeff_test( prism, 1, prism_coeff, 9, xi );
}

void LinearMappingFunctionTest::test_linear_prism_coeff_faces()
{
  double xi[15];
  xi_at_faces( PRISM, xi, PrismCorners );
  do_coeff_test( prism, 2, prism_coeff, 5, xi );
}

void LinearMappingFunctionTest::test_linear_prism_coeff_center()
{
  double xi[3] = { 0.5, 1./3, 1./3 };
  do_coeff_test( prism, 3, prism_coeff, 1, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_coeff_corners()
{
  double xi[15] = {  0,  0,  0,
                     1,  0,  0,
                     1,  1,  0,
                     0,  1,  0,
                     0,  0,  1 };
  do_coeff_test( pyr, 0, pyr_coeff, 5, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_coeff_edges()
{
  double xi[24] = { 0.5, 0.0, 0.0,
                    1.0, 0.5, 0.0,
                    0.5, 1.0, 0.0,
                    0.0, 0.5, 0.0,
                    0.0, 0.0, 0.5,
                    1.0, 0.0, 0.5,
                    1.0, 1.0, 0.5,
                    0.0, 1.0, 0.5 };
  do_coeff_test( pyr, 1, pyr_coeff, 8, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_coeff_faces()
{
  double xi[15] = { 0.5, 0.0, 0.5,
                    1.0, 0.5, 0.5,
                    0.5, 1.0, 0.5,
                    0.0, 0.5, 0.5,
                    0.5, 0.5, 0.0 };
  do_coeff_test( pyr, 2, pyr_coeff, 5, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_coeff_center()
{
  double xi[3] = { 0.5, 0.5, 0.5 };
  do_coeff_test( pyr, 3, pyr_coeff, 1, xi );
}



/*******************************************************************************
 *                 Test functions for mapping function derivatives
 *******************************************************************************/

void LinearMappingFunctionTest::test_linear_hex_deriv_corners()
{ 
  double xi[24];
  xi_at_corners( HEXAHEDRON, xi, &HexCorners[0][0] );
  do_deriv_test( hex, 0, hex_deriv, 8, xi );
}
  
void LinearMappingFunctionTest::test_linear_hex_deriv_edges()
{
  double xi[36];
  xi_at_edges( HEXAHEDRON, xi, &HexCorners[0][0] );
  do_deriv_test( hex, 1, hex_deriv, 12, xi );
}

void LinearMappingFunctionTest::test_linear_hex_deriv_faces()
{
  double xi[18];
  xi_at_faces( HEXAHEDRON, xi, &HexCorners[0][0] );
  do_deriv_test( hex, 2, hex_deriv, 6, xi );
}

void LinearMappingFunctionTest::test_linear_hex_deriv_center()
{
  double xi[3] = { 0.5, 0.5, 0.5 };
  do_deriv_test( hex, 3, hex_deriv, 1, xi );
}

void LinearMappingFunctionTest::test_linear_quad_deriv_corners()
{ 
  double xi[8];
  xi_at_corners( QUADRILATERAL, xi, &QuadCorners[0][0] );
  do_deriv_test( quad, 0, quad_deriv, 4, xi );
}

void LinearMappingFunctionTest::test_linear_quad_deriv_edges()
{
  double xi[8];
  xi_at_edges( QUADRILATERAL, xi, &QuadCorners[0][0] );
  do_deriv_test( quad, 1, quad_deriv, 4, xi );
}

//void LinearMappingFunctionTest::test_linear_quad_deriv_faces()
//{
//  test_deriv_fail( quad, 3 );
//}

void LinearMappingFunctionTest::test_linear_quad_deriv_center()
{
  double xi[2] = { 0.5, 0.5 };
  do_deriv_test( quad, 2, quad_deriv, 1, xi );
}

void LinearMappingFunctionTest::test_linear_tet_deriv_corners()
{
  double xi[12];
  xi_at_corners( TETRAHEDRON, xi, TetCorners );
  do_deriv_test( tet, 0, tet_deriv, 4, xi );
}

void LinearMappingFunctionTest::test_linear_tet_deriv_edges()
{
  double xi[18];
  xi_at_edges( TETRAHEDRON, xi, TetCorners );
  do_deriv_test( tet, 1, tet_deriv, 6, xi );
}

void LinearMappingFunctionTest::test_linear_tet_deriv_faces()
{
  double xi[12];
  xi_at_faces( TETRAHEDRON, xi, TetCorners );
  do_deriv_test( tet, 2, tet_deriv, 4, xi );
}

void LinearMappingFunctionTest::test_linear_tet_deriv_center()
{
  double xi[3] = { 0.25, 0.25, 0.25 };
  do_deriv_test( tet, 3, tet_deriv, 1, xi );
}
  
void LinearMappingFunctionTest::test_linear_tri_deriv_corners()
{
  double xi[12];
  xi_at_corners( TRIANGLE, xi, TriCorners );
  do_deriv_test( tri, 0, tri_deriv, 3, xi );
}

void LinearMappingFunctionTest::test_linear_tri_deriv_edges()
{
  double xi[18];
  xi_at_edges( TRIANGLE, xi, TriCorners );
  do_deriv_test( tri, 1, tri_deriv, 3, xi );
}

//void LinearMappingFunctionTest::test_linear_tri_deriv_faces()
//{
//  test_deriv_fail( tri, 3 );
//}

void LinearMappingFunctionTest::test_linear_tri_deriv_center()
{
  double xi[2] = { 1./3, 1./3 };
  do_deriv_test( tri, 2, tri_deriv, 1, xi );
}
   
void LinearMappingFunctionTest::test_linear_prism_deriv_corners()
{
  double xi[18];
  xi_at_corners( PRISM, xi, PrismCorners );
  do_deriv_test( prism, 0, prism_deriv, 6, xi );
}

void LinearMappingFunctionTest::test_linear_prism_deriv_edges()
{
  double xi[27];
  xi_at_edges( PRISM, xi, PrismCorners );
  do_deriv_test( prism, 1, prism_deriv, 9, xi );
}

void LinearMappingFunctionTest::test_linear_prism_deriv_faces()
{
  double xi[15];
  xi_at_faces( PRISM, xi, PrismCorners );
  do_deriv_test( prism, 2, prism_deriv, 5, xi );
}

void LinearMappingFunctionTest::test_linear_prism_deriv_center()
{
  double xi[3] = { 0.5, 1./3, 1./3 };
  do_deriv_test( prism, 3, prism_deriv, 1, xi );
}
  


void LinearMappingFunctionTest::test_linear_pyr_deriv_corners()
{
  double xi[15] = {  0,  0,  0,
                     1,  0,  0,
                     1,  1,  0,
                     0,  1,  0,
                     0.5,  0.5,  1 };
  do_deriv_test( pyr, 0, pyr_deriv, 5, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_deriv_edges()
{
  double xi[24] = { 0.5, 0.0, 0.0,
                    1.0, 0.5, 0.0,
                    0.5, 1.0, 0.0,
                    0.0, 0.5, 0.0,
                    0.0, 0.0, 0.5,
                    1.0, 0.0, 0.5,
                    1.0, 1.0, 0.5,
                    0.0, 1.0, 0.5 };
  do_deriv_test( pyr, 1, pyr_deriv, 8, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_deriv_faces()
{
  double xi[15] = { 0.5, 0.0, 0.5,
                    1.0, 0.5, 0.5,
                    0.5, 1.0, 0.5,
                    0.0, 0.5, 0.5,
                    0.5, 0.5, 0.0 };
  do_deriv_test( pyr, 2, pyr_deriv, 5, xi );
}

void LinearMappingFunctionTest::test_linear_pyr_deriv_center()
{
  double xi[3] = { 0.5, 0.5, 0.5 };
  do_deriv_test( pyr, 3, pyr_deriv, 1, xi );
}


/*******************************************************************************
 *        Mapping function implementations to compare with
 *******************************************************************************/
    
void LinearMappingFunctionTest::hex_coeff( double xi[3], double coeff[8] )
{
  coeff[0] = (1-xi[0]) * (1-xi[1]) * (1-xi[2]);
  coeff[1] =    xi[0]  * (1-xi[1]) * (1-xi[2]);
  coeff[2] =    xi[0]  *    xi[1]  * (1-xi[2]);
  coeff[3] = (1-xi[0]) *    xi[1]  * (1-xi[2]);
  coeff[4] = (1-xi[0]) * (1-xi[1]) *    xi[2] ;
  coeff[5] =    xi[0]  * (1-xi[1]) *    xi[2] ;
  coeff[6] =    xi[0]  *    xi[1]  *    xi[2] ;
  coeff[7] = (1-xi[0]) *    xi[1]  *    xi[2] ;
}

void LinearMappingFunctionTest::tet_coeff( double xi[3], double coeff[4] )
{
  coeff[0] = 1 - xi[0] - xi[1] - xi[2];
  coeff[1] = xi[0];
  coeff[2] = xi[1];
  coeff[3] = xi[2];
}

void LinearMappingFunctionTest::quad_coeff( double xi[2], double coeff[4] )
{
  coeff[0] = (1-xi[0]) * (1-xi[1]);
  coeff[1] =    xi[0]  * (1-xi[1]);
  coeff[2] =    xi[0]  *    xi[1] ;
  coeff[3] = (1-xi[0]) *    xi[1] ;
}

void LinearMappingFunctionTest::tri_coeff( double xi[2], double coeff[3] )
{
  coeff[0] = 1 - xi[0] - xi[1];
  coeff[1] = xi[0];
  coeff[2] = xi[1];
}

void LinearMappingFunctionTest::prism_coeff( double xi[3], double coeff[6] )
{
  coeff[0] = (1 - xi[0]) * (1 - xi[1] - xi[2]);
  coeff[1] = (1 - xi[0]) * xi[1];
  coeff[2] = (1 - xi[0]) * xi[2];
  coeff[3] =      xi[0]  * (1 - xi[1] - xi[2]);
  coeff[4] =      xi[0]  * xi[1];
  coeff[5] =      xi[0]  * xi[2];
}

void LinearMappingFunctionTest::pyr_coeff( double xi[3], double coeff[5] )
{
  coeff[0] = (1 - xi[0]) * (1 - xi[1]) * (1 - xi[2]);
  coeff[1] =      xi[0]  * (1 - xi[1]) * (1 - xi[2]);
  coeff[2] =      xi[0]  *      xi[1]  * (1 - xi[2]);
  coeff[3] = (1 - xi[0]) *      xi[1]  * (1 - xi[2]);
  coeff[4] =                                  xi[2] ;
}

/*******************************************************************************
 *        Mapping function derivatives to compare with
 *******************************************************************************/

void LinearMappingFunctionTest::hex_deriv( double xi[3], double coeff[24] )
{
  coeff[3*0+0] = - (1-xi[1]) * (1-xi[2]);
  coeff[3*0+1] = - (1-xi[0]) * (1-xi[2]);
  coeff[3*0+2] = - (1-xi[0]) * (1-xi[1]);

  coeff[3*1+0] =   (1-xi[1]) * (1-xi[2]);
  coeff[3*1+1] = -    xi[0]  * (1-xi[2]);
  coeff[3*1+2] = -    xi[0]  * (1-xi[1]);

  coeff[3*2+0] =      xi[1]  * (1-xi[2]);
  coeff[3*2+1] =      xi[0]  * (1-xi[2]);
  coeff[3*2+2] = -    xi[0]  *    xi[1] ;

  coeff[3*3+0] = -    xi[1]  * (1-xi[2]);
  coeff[3*3+1] =   (1-xi[0]) * (1-xi[2]);
  coeff[3*3+2] = - (1-xi[0]) *    xi[1] ;

  coeff[3*4+0] = - (1-xi[1]) *    xi[2] ;
  coeff[3*4+1] = - (1-xi[0]) *    xi[2] ;
  coeff[3*4+2] =   (1-xi[0]) * (1-xi[1]);

  coeff[3*5+0] =   (1-xi[1]) *    xi[2] ;
  coeff[3*5+1] = -    xi[0]  *    xi[2] ;
  coeff[3*5+2] =      xi[0]  * (1-xi[1]);

  coeff[3*6+0] =      xi[1]  *    xi[2] ;
  coeff[3*6+1] =      xi[0]  *    xi[2] ;
  coeff[3*6+2] =      xi[0]  *    xi[1] ;

  coeff[3*7+0] = -    xi[1]  *    xi[2] ;
  coeff[3*7+1] =   (1-xi[0]) *    xi[2] ;
  coeff[3*7+2] =   (1-xi[0]) *    xi[1] ;
}

void LinearMappingFunctionTest::tet_deriv( double*, double coeff[12] )
{
  static const double derivs[] = {-1,-1,-1,
                                   1, 0, 0, 
                                   0, 1, 0,
                                   0, 0, 1 };
  memcpy( coeff, derivs, sizeof(derivs) );
}

void LinearMappingFunctionTest::quad_deriv( double xi[2], double coeff[8] )
{
  coeff[2*0+0] = xi[1] - 1;
  coeff[2*0+1] = xi[0] - 1;
  
  coeff[2*1+0] = 1 - xi[1];
  coeff[2*1+1] = -xi[0];
  
  coeff[2*2+0] = xi[1];
  coeff[2*2+1] = xi[0];
  
  coeff[2*3+0] = -xi[1];
  coeff[2*3+1] = 1 - xi[0];
}

void LinearMappingFunctionTest::tri_deriv( double*, double coeff[6] )
{
  static const double derivs[] = {-1,-1,
                                   1, 0, 
                                   0, 1 };
  memcpy( coeff, derivs, sizeof(derivs) );
}

void LinearMappingFunctionTest::prism_deriv( double xi[3], double coeff[18] )
{
  coeff[ 0] = xi[1] + xi[2] - 1.0;;
  coeff[ 1] = xi[0] - 1.0;
  coeff[ 2] = xi[0] - 1.0;
  
  coeff[ 3] = -xi[1];
  coeff[ 4] = 1.0 - xi[0];
  coeff[ 5] = 0.0;
  
  coeff[ 6] = -xi[2];
  coeff[ 7] = 0.0;
  coeff[ 8] = 1.0 - xi[0];
  
  coeff[ 9] =  1.0 - xi[1] - xi[2];
  coeff[10] = -xi[0];
  coeff[11] = -xi[0];
  
  coeff[12] =  xi[1];
  coeff[13] =  xi[0];
  coeff[14] =  0.0;
  
  coeff[15] =  xi[2];
  coeff[16] =  0.0;
  coeff[17] =  xi[0];
}

void LinearMappingFunctionTest::pyr_deriv( double xi[3], double coeff[15] )
{
  coeff[3*0+0] = -(1-xi[1])*(1-xi[2]);
  coeff[3*0+1] = -(1-xi[0])*(1-xi[2]);
  coeff[3*0+2] = -(1-xi[0])*(1-xi[1]);
  
  coeff[3*1+0] =  (1-xi[1])*(1-xi[2]);
  coeff[3*1+1] =    -xi[0] *(1-xi[2]);
  coeff[3*1+2] =    -xi[0] *(1-xi[1]);
  
  coeff[3*2+0] =     xi[1] *(1-xi[2]);
  coeff[3*2+1] =     xi[0] *(1-xi[2]);
  coeff[3*2+2] =    -xi[0] *   xi[1] ;
  
  coeff[3*3+0] =    -xi[1] *(1-xi[2]);
  coeff[3*3+1] =  (1-xi[0])*(1-xi[2]);
  coeff[3*3+2] =    -xi[1] *(1-xi[0]);
  
  coeff[3*4+0] = 0.0;
  coeff[3*4+1] = 0.0;
  coeff[3*4+2] = 1.0;
}

/*******************************************************************************
 *         Some utlity functions for creating CppUnit error messages
 *******************************************************************************/

static string itostr( int i )
{
  char buffer[32];
  sprintf(buffer, "%d", i );
  return buffer;
}

static string dtostr( double i )
{
  char buffer[32];
  sprintf(buffer, "%g", i );
  return buffer;
}



/*******************************************************************************
 *         Actual test imlplementation (common code for many tests)
 *******************************************************************************/
    
void LinearMappingFunctionTest::do_coeff_test( MappingFunction& mf, 
                                               unsigned subdim,
                                               map_func mf2,
                                               unsigned count,
                                               double* xi )
{
    // make sure it fails if passed a nonlinear element
  MsqError err;
  double coeff[100];
  size_t indices[100];
  size_t num_coeff = 100;
  NodeSet tmp_set;
  tmp_set.set_mid_edge_node(1);
  mf.coefficients( Sample(0, 1), tmp_set, coeff, indices, num_coeff, err );
  CPPUNIT_ASSERT(err);
  err.clear();
  
    // get number of vertices in element
  const unsigned n = TopologyInfo::corners( mf.element_topology() );
  const unsigned d = TopologyInfo::dimension( mf.element_topology() );
  
    // compare coefficients at each location
  vector<double> comp(n);
  for (unsigned i = 0; i < count; ++i)
  {
    num_coeff = 101;
    mf.coefficients( Sample(subdim, i), NodeSet(), coeff, indices, num_coeff, err );
    CPPUNIT_ASSERT(!err);
    
    mf2( xi, &comp[0] );
    string xi_str;
    for (unsigned j = 0; j < d; ++j) {
      xi_str += !j ? "(" : ", ";
      xi_str += dtostr(xi[j]);
    }
    xi_str += ")";
    xi += d;
    
    for (unsigned j = 0; j < n; ++j)
    {
      double mf_val = 0.0;
      size_t idx = std::find( indices, indices+num_coeff, j ) - indices;
      if (idx < num_coeff)
	mf_val = coeff[idx];

      CppUnit::Message message( "Coefficients do not match." );
      message.addDetail( string("Entity:             ") + itostr( i ) );
      message.addDetail( string("Coefficient number: ") + itostr( j ) );
      message.addDetail( string("Xi:             ") + xi_str );
      message.addDetail( string("Expected value: ") + dtostr( comp[j] ) );
      message.addDetail( string("Actual value:   ") + dtostr( mf_val ) );
      ASSERT_MESSAGE( message, fabs(comp[j]-mf_val) < DBL_EPSILON );
    }
  }
}


void LinearMappingFunctionTest::do_deriv_test( MappingFunction2D& mf, 
                                               unsigned subdim,
                                               map_func mf2,
                                               unsigned count,
                                               double* xi )
{
    // make sure it fails if passed a nonlinear element
  MsqError err;
  MsqVector<2> derivs[100];
  size_t verts[100], num_vtx = 37;
  NodeSet tmp_set;
  tmp_set.set_mid_edge_node(1);
  mf.derivatives( Sample(subdim, 0), tmp_set, verts, derivs, num_vtx, err );
  CPPUNIT_ASSERT(err);
  err.clear();
  
    // get number of vertices in element
  const unsigned n = TopologyInfo::corners( mf.element_topology() );
  
    // compare coefficients at each location
  vector<double> comp(2*n);
  for (unsigned i = 0; i < count; ++i)
  {
    num_vtx = 33;
    mf.derivatives( Sample(subdim, i), NodeSet(), verts, derivs, num_vtx, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( num_vtx > 0 );
    CPPUNIT_ASSERT( num_vtx <= n );
    
    mf2( xi, &comp[0] );
    string xi_str;
    for (unsigned j = 0; j < 2; ++j) {
      xi_str += !j ? "(" : ", ";
      xi_str += dtostr(xi[j]);
    }
    xi_str += ")";
    xi += 2;
    
    for (unsigned j = 0; j < num_vtx; ++j)
    {
      bool all_zero = true;
      for (unsigned k = 0; k < 2; ++k)
      {
        CppUnit::Message message( "Coefficient derivatives do not match." );
        message.addDetail( string("Entity:             ") + itostr( i ) );
        message.addDetail( string("Coefficient number: ") + itostr( j ) );
        message.addDetail( string("Xi:             ") + xi_str );
        message.addDetail( string("Axis:           ") + itostr( k ) );
        message.addDetail( string("Expected value: ") + dtostr( comp[2*verts[j]+k] ) );
        message.addDetail( string("Actual value:   ") + dtostr( derivs[j][k] ) );
        ASSERT_MESSAGE( message, fabs(comp[2*verts[j]+k]-derivs[j][k]) < DBL_EPSILON );
        if (fabs(derivs[j][k]) > DBL_EPSILON)
          all_zero = false;
      }

        // if vertex has all zero values, it shouldn't have been in the
        // vertex list at all, as the Jacobian will not depend on that vertex.
      CPPUNIT_ASSERT( !all_zero );
    }
    
      // If any vertex is not in the list, then its values must be zero.
    sort( verts, verts + num_vtx );
    for (unsigned j = 0; j < num_vtx; ++j) {
      if (!binary_search( verts, verts+num_vtx, j )) {
        for (unsigned k = 0; k < 2; ++k)
        {
          CppUnit::Message message( "Missing coefficient derivatives." );
          message.addDetail( string("Entity:              ") + itostr( i ) );
          message.addDetail( string("Coefficient number:  ") + itostr( j ) );
          message.addDetail( string("Axis:                ") + itostr( k ) );
          message.addDetail( string("Expected derivative: ") + dtostr( comp[2*j+k] ) );
          ASSERT_MESSAGE( message, fabs(comp[2*j+k]) < DBL_EPSILON );
        }
      }
    }
  }
}

void LinearMappingFunctionTest::do_deriv_test( MappingFunction3D& mf, 
                                               unsigned subdim,
                                               map_func mf2,
                                               unsigned count,
                                               double* xi )
{
    // make sure it fails if passed a nonlinear element
  MsqError err;
  MsqVector<3> derivs[100];
  size_t verts[100], num_vtx = 37;
  NodeSet tmp_set;
  tmp_set.set_mid_edge_node(1);
  mf.derivatives( Sample(subdim, 0), tmp_set, verts, derivs, num_vtx, err );
  CPPUNIT_ASSERT(err);
  err.clear();
  
    // get number of vertices in element
  const unsigned n = TopologyInfo::corners( mf.element_topology() );
  
    // compare coefficients at each location
  vector<double> comp(3*n);
  for (unsigned i = 0; i < count; ++i)
  {
    num_vtx = 33;
    mf.derivatives( Sample(subdim, i), NodeSet(), verts, derivs, num_vtx, err );
    CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT( num_vtx > 0 );
    CPPUNIT_ASSERT( num_vtx <= n );
    
    mf2( xi, &comp[0] );
    string xi_str;
    for (unsigned j = 0; j < 3; ++j) {
      xi_str += !j ? "(" : ", ";
      xi_str += dtostr(xi[j]);
    }
    xi_str += ")";
    xi += 3;
    
    for (unsigned j = 0; j < num_vtx; ++j)
    {
      bool all_zero = true;
      for (unsigned k = 0; k < 3; ++k)
      {
        CppUnit::Message message( "Coefficient derivatives do not match." );
        message.addDetail( string("Entity:             ") + itostr( i ) );
        message.addDetail( string("Coefficient number: ") + itostr( j ) );
        message.addDetail( string("Xi:             ") + xi_str );
        message.addDetail( string("Axis:           ") + itostr( k ) );
        message.addDetail( string("Expected value: ") + dtostr( comp[3*verts[j]+k] ) );
        message.addDetail( string("Actual value:   ") + dtostr( derivs[j][k] ) );
        ASSERT_MESSAGE( message, fabs(comp[3*verts[j]+k]-derivs[j][k]) < DBL_EPSILON );
        if (fabs(derivs[j][k]) > DBL_EPSILON)
          all_zero = false;
      }

        // if vertex has all zero values, it shouldn't have been in the
        // vertex list at all, as the Jacobian will not depend on that vertex.
      CPPUNIT_ASSERT( !all_zero );
    }
    
      // If any vertex is not in the list, then its values must be zero.
    sort( verts, verts + num_vtx );
    for (unsigned j = 0; j < num_vtx; ++j) {
      if (!binary_search( verts, verts+num_vtx, j )) {
        for (unsigned k = 0; k < 3; ++k)
        {
          CppUnit::Message message( "Missing coefficient derivatives." );
          message.addDetail( string("Entity:              ") + itostr( i ) );
          message.addDetail( string("Coefficient number:  ") + itostr( j ) );
          message.addDetail( string("Axis:                ") + itostr( k ) );
          message.addDetail( string("Expected derivative: ") + dtostr( comp[3*j+k] ) );
          ASSERT_MESSAGE( message, fabs(comp[3*j+k]) < DBL_EPSILON );
        }
      }
    }
  }
}

 
void LinearMappingFunctionTest::do_ideal_test( MappingFunction2D& mf )
{
  MsqError err;
  MsqMatrix<3,2> W_prime;
  mf.ideal( Sample(2,0), W_prime, err );
  ASSERT_NO_ERROR(err);
  
    // for this test that everything is in the xy-plane
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, W_prime(2,0), 1e-12 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, W_prime(2,1), 1e-12 );
  MsqMatrix<2,2> W( W_prime.data() );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(W), 1e-6 );

  const Vector3D* verts = unit_edge_element( mf.element_topology() );
  CPPUNIT_ASSERT(verts);
  
  JacobianCalculator jc;
  jc.get_Jacobian_2D( &mf, NodeSet(), Sample(2,0), verts, TopologyInfo::corners(mf.element_topology()), W_prime, err );
  ASSERT_NO_ERROR(err);
  
    // for this test that everything is in the xy-plane
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, W_prime(2,0), 1e-12 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, W_prime(2,1), 1e-12 );
  MsqMatrix<2,2> W_exp( W_prime.data() );
  W_exp /= sqrt(det(W_exp));
  
    // Matrices should be a rotation of each other.
    // First, calculate tentative rotation matrix
  MsqMatrix<2,2> R = inverse(W_exp) * W;
    // next check that it is a rotation
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(R), 1e-6 ); // no scaling
  ASSERT_MATRICES_EQUAL( transpose(R), inverse(R), 1e-6 ); // orthogonal
}

void LinearMappingFunctionTest::do_ideal_test( MappingFunction3D& mf )
{
  MsqError err;
  MsqMatrix<3,3> W, I(1.0);
  mf.ideal( Sample(3,0), W, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(W), 1e-6 );

  const Vector3D* verts = unit_edge_element( mf.element_topology() );
  CPPUNIT_ASSERT(verts);
  
  JacobianCalculator jc;
  MsqMatrix<3,3> W_exp;
  jc.get_Jacobian_3D( &mf, NodeSet(), Sample(3,0), verts, TopologyInfo::corners(mf.element_topology()), W_exp, err );
  ASSERT_NO_ERROR(err);
  W_exp /= Mesquite::cbrt(det(W_exp));
  
    // Matrices should be a rotation of each other.
    // First, calculate tentative rotation matrix
  MsqMatrix<3,3> R = W * inverse(W_exp);
    // next check that it is a rotation
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(R), 1e-6 ); // no scaling
  ASSERT_MATRICES_EQUAL( I, transpose(R) * R, 1e-6 ); // orthogonal
}

                 
void LinearMappingFunctionTest::test_coeff_fail( MappingFunction& mf, unsigned subdim )
{
    // make sure it fails if called
  MsqError err;
  double coeff[100];
  size_t num_coeff, indices[100];
  mf.coefficients( Sample(subdim, 0), NodeSet(), coeff, indices, num_coeff, err );
  CPPUNIT_ASSERT(err);
  err.clear();
}  

void LinearMappingFunctionTest::test_deriv_fail( MappingFunction2D& mf, unsigned subdim )
{
    // make sure it fails if called
  MsqError err;
  MsqVector<2> coeff[100];
  size_t verts[100], num_coeff;
  mf.derivatives( Sample(subdim, 0), NodeSet(), verts, coeff, num_coeff, err );
  CPPUNIT_ASSERT(err);
  err.clear();
}

void LinearMappingFunctionTest::test_deriv_fail( MappingFunction3D& mf, unsigned subdim )
{
    // make sure it fails if called
  MsqError err;
  MsqVector<3> coeff[100];
  size_t verts[100], num_coeff;
  mf.derivatives( Sample(subdim, 0), NodeSet(), verts, coeff, num_coeff, err );
  CPPUNIT_ASSERT(err);
  err.clear();
}
