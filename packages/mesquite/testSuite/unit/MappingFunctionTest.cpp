/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file MappingFunctionTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "UnitUtil.hpp"
#include "LinearTetrahedron.hpp"
#include "LinearTriangle.hpp"
#include "PatchData.hpp"
#include "IdealElements.hpp"
#include <algorithm>

using namespace Mesquite;

class MappingFunctionTest : public CppUnit::TestFixture
{
  private:
    CPPUNIT_TEST_SUITE(MappingFunctionTest);
    CPPUNIT_TEST(test_jacobian_2d);
    CPPUNIT_TEST(test_jacobian_3d);
    CPPUNIT_TEST(test_ideal_2d);
    CPPUNIT_TEST(test_ideal_3d);
    CPPUNIT_TEST_SUITE_END();
  
    LinearTriangle tri;
    LinearTetrahedron tet;
  public:
  
    void test_jacobian_2d();
    void test_jacobian_3d();
    void test_ideal_2d();
    void test_ideal_3d();
    
    MsqMatrix<3,2> calculate_jacobian_2d( Sample s, const Vector3D* coords );
    MsqMatrix<3,3> calculate_jacobian_3d( Sample s, const Vector3D* coords );
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MappingFunctionTest, "MappingFunctionTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MappingFunctionTest, "Unit");
   
void MappingFunctionTest::test_jacobian_2d()
{
  MsqError err;
  PatchData pd;
  const size_t conn[3] = { 0, 1, 2 };
  const bool fixed[3] = { false, false, false };
  const double coords[3][3] = { {0.5, 0.5, 0.0},
                                {2.0, 0.2, 0.0},
                                {0.7, 0.8, 0.0} };
                                
  size_t indices_exp[3], indices_act[3], num_exp, num_act;
  MsqVector<2> coeff_exp[3], coeff_act[3];
  MsqMatrix<3,2> J_exp, J_act;
  for (int r = 0; r < 3; ++r ) {
    const double verts[] = { coords[ r     ][0], coords[ r     ][1], coords[ r     ][2],
                             coords[(r+1)%3][0], coords[(r+1)%3][1], coords[(r+1)%3][2],
                             coords[(r+2)%3][0], coords[(r+2)%3][1], coords[(r+2)%3][2] };
    pd.fill( 3, verts, 1, TRIANGLE, conn, fixed, err ); 
    ASSERT_NO_ERROR(err);
    
    for (int d = 0; d < 3; ++d) {
      int n = (d == 2) ? 1 : 3;
      for (int s = 0; s < n; ++s) {
        tri.MappingFunction2D::jacobian( pd, 0, NodeSet(), Sample(d,s), indices_act, coeff_act, num_act, J_act, err );
        ASSERT_NO_ERROR(err);
        CPPUNIT_ASSERT(num_act <= 3);
        CPPUNIT_ASSERT(num_act > 0);
        
        tri.derivatives( Sample(d,s), NodeSet(), indices_exp, coeff_exp, num_exp, err );
        ASSERT_NO_ERROR(err);
        CPPUNIT_ASSERT_EQUAL(num_exp, num_act);
        
        J_exp = MsqMatrix<3,2>(0.0); // zero
        for (size_t j = 0; j < num_exp; ++j) {
          size_t idx = pd.element_by_index(0).get_vertex_index_array()[j];
          size_t k = std::find( indices_act, indices_act+num_act, idx ) - indices_act;
          CPPUNIT_ASSERT( k < num_act ); // found
          ASSERT_MATRICES_EQUAL( coeff_exp[j], coeff_act[k], 1e-10 );
          
          J_exp += MsqVector<3>( pd.vertex_by_index( idx ).to_array() ) * transpose( coeff_exp[j] );
        }
        
        ASSERT_MATRICES_EQUAL( J_exp, J_act, 1e-6 );
      }
    }
  }
}
   
void MappingFunctionTest::test_jacobian_3d()
{
  MsqError err;
  PatchData pd;
  const size_t conn[4] = { 0, 1, 2, 3 };
  const bool fixed[4] = { false, false, false, false };
  const double coords[4][3] = { {0.5, 0.5, 0.0},
                                {2.0, 0.2, 0.0},
                                {0.7, 0.8, 0.0},
                                {0.6, 1.1, 13.} };
                                
  size_t indices_exp[4], indices_act[4], num_exp, num_act;
  MsqVector<3> coeff_exp[4], coeff_act[4];
  MsqMatrix<3,3> J_exp, J_act;
  for (int r = 0; r < 3; ++r ) {
    const double verts[] = { coords[ r     ][0], coords[ r     ][1], coords[ r     ][2],
                             coords[(r+1)%3][0], coords[(r+1)%3][1], coords[(r+1)%3][2],
                             coords[(r+2)%3][0], coords[(r+2)%3][1], coords[(r+2)%3][2],
                             coords[ 3     ][0], coords[ 3     ][1], coords[ 3     ][2]};
    pd.fill( 4, verts, 1, TETRAHEDRON, conn, fixed, err ); 
    ASSERT_NO_ERROR(err);
    
    for (int d = 0; d < 3; ++d) {
      int n = (d == 3) ? 1 : (d ==1 ) ? 6 : 4;
      for (int s = 0; s < n; ++s) {
        tet.MappingFunction3D::jacobian( pd, 0, NodeSet(), Sample(d,s), indices_act, coeff_act, num_act, J_act, err );
        ASSERT_NO_ERROR(err);
        CPPUNIT_ASSERT(num_act <= 4);
        CPPUNIT_ASSERT(num_act > 0);
        
        tet.derivatives( Sample(d,s), NodeSet(), indices_exp, coeff_exp, num_exp, err );
        ASSERT_NO_ERROR(err);
        CPPUNIT_ASSERT_EQUAL(num_exp, num_act);
        
        J_exp = MsqMatrix<3,3>(0.0); // zero
        for (size_t j = 0; j < num_exp; ++j) {
          size_t idx = pd.element_by_index(0).get_vertex_index_array()[j];
          size_t k = std::find( indices_act, indices_act+num_act, idx ) - indices_act;
          CPPUNIT_ASSERT( k < num_act ); // found
          ASSERT_MATRICES_EQUAL( coeff_exp[j], coeff_act[k], 1e-10 );
          
          J_exp += MsqVector<3>( pd.vertex_by_index( idx ).to_array() ) * transpose( coeff_exp[j] );
        }
        
        ASSERT_MATRICES_EQUAL( J_exp, J_act, 1e-6 );
      }
    }
  }
}
   
void MappingFunctionTest::test_ideal_2d()
{
  MsqError err;
  MsqMatrix<3,2> W_prime;
  tri.MappingFunction2D::ideal( Sample(2,0), W_prime, err );
  ASSERT_NO_ERROR(err);
  
    // for this test that everything is in the xy-plane
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, W_prime(2,0), 1e-12 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, W_prime(2,1), 1e-12 );
  MsqMatrix<2,2> W( W_prime.data() );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(W), 1e-6 );

  const Vector3D* verts = unit_edge_element( TRIANGLE );
  CPPUNIT_ASSERT(verts);
  
  size_t indices[3], num;
  MsqVector<2> coeff[3];
  tri.derivatives( Sample(2,0), NodeSet(), indices, coeff, num, err );
  ASSERT_NO_ERROR(err);

  MsqMatrix<3,2> J_exp(0.0); // zero
  for (size_t j = 0; j < num; ++j) 
    J_exp += MsqVector<3>( verts[j].to_array() ) * transpose( coeff[j] );
  
    // for this test that everything is in the xy-plane
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, W_prime(2,0), 1e-12 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, W_prime(2,1), 1e-12 );
  MsqMatrix<2,2> W_exp( J_exp.data() );
  W_exp /= sqrt(det(W_exp));
  
    // Matrices should be a rotation of each other.
    // First, calculate tentative rotation matrix
  MsqMatrix<2,2> R = inverse(W_exp) * W;
    // next check that it is a rotation
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(R), 1e-6 ); // no scaling
  ASSERT_MATRICES_EQUAL( transpose(R), inverse(R), 1e-6 ); // orthogonal
}

void MappingFunctionTest::test_ideal_3d()
{
  MsqError err;
  MsqMatrix<3,3> J;
  tet.MappingFunction3D::ideal( Sample(3,0), J, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(J), 1e-6 );

  const Vector3D* verts = unit_edge_element( TETRAHEDRON );
  CPPUNIT_ASSERT(verts);
  
  size_t indices[4], num;
  MsqVector<3> coeff[4];
  tet.derivatives( Sample(2,0), NodeSet(), indices, coeff, num, err );
  ASSERT_NO_ERROR(err);

  MsqMatrix<3,3> J_exp(0.0); // zero
  for (size_t j = 0; j < num; ++j) 
    J_exp += MsqVector<3>( verts[j].to_array() ) * transpose( coeff[j] );
  J_exp /= Mesquite::cbrt(det(J_exp));
  
    // Matrices should be a rotation of each other.
    // First, calculate tentative rotation matrix
  MsqMatrix<3,3> R = inverse(J_exp) * J;
    // next check that it is a rotation
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, det(R), 1e-6 ); // no scaling
  ASSERT_MATRICES_EQUAL( transpose(R), inverse(R), 1e-6 ); // orthogonal
}

