/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */

/*! \file QualityMetricTest.cpp

Unit testing for the QualityMetric class
\author Jason Kraftcheck
*/
#include "Mesquite.hpp"
#include "Mesquite_VertexQM.hpp"
#include "Mesquite_ElementQM.hpp"
#include "Mesquite_IdealElements.hpp"
#include "UnitUtil.hpp"
#include "Mesquite_TopologyInfo.hpp"
#include "Mesquite_PatchData.hpp"
#include "Mesquite_cppunit/extensions/HelperMacros.h"
#include "Mesquite_TMPDerivs.hpp"

#include <string>

using namespace Mesquite;

class QualityMetricTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QualityMetricTest);
  CPPUNIT_TEST (test_fixed_vertex_list);
  CPPUNIT_TEST (test_remove_fixed_gradients);
  CPPUNIT_TEST (test_remove_fixed_hessians);
  CPPUNIT_TEST( test_gradient_constant );
  CPPUNIT_TEST( test_gradient_linear );
  CPPUNIT_TEST( test_gradient_parabolic );
  CPPUNIT_TEST( test_gradient_tau );
  CPPUNIT_TEST( test_Hessian_constant );
  CPPUNIT_TEST( test_Hessian_linear );
  CPPUNIT_TEST( test_Hessian_parabolic );
  CPPUNIT_TEST( test_Hessian_tau );
  CPPUNIT_TEST( test_diagonal_constant );
  CPPUNIT_TEST( test_diagonal_linear );
  CPPUNIT_TEST( test_diagonal_parabolic );
  CPPUNIT_TEST( test_diagonal_tau );
  CPPUNIT_TEST_SUITE_END();
  PatchData tri_pd;

public:
  void setUp();
  
  void test_fixed_vertex_list();
  void test_remove_fixed_gradients();
  void test_remove_fixed_hessians();
  void test_gradient_constant();
  void test_gradient_linear();
  void test_gradient_parabolic();
  void test_gradient_tau();
  void test_Hessian_constant();
  void test_Hessian_linear();
  void test_Hessian_parabolic();
  void test_Hessian_tau();
  void test_diagonal_constant();
  void test_diagonal_linear();
  void test_diagonal_parabolic();
  void test_diagonal_tau();
  
  void compare_indices( QualityMetric& qm, PatchData& pd, size_t sample,
                        double value,
                        const std::vector<size_t>& indices );
  void compare_gradient( QualityMetric& qm, PatchData& pd, size_t sample,
                         double value,
                         const std::vector<size_t>& indices,
                         const std::vector<Vector3D>& grad );
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QualityMetricTest, "QualityMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QualityMetricTest, "Unit");


// Create a single-triangle patch with one fixed vertex
void QualityMetricTest::setUp()
{
  MsqError err;
  const double vtx_coords[] = { 1, 0, 0,
                                0, 1, 0,
                                3, 0, 0,
                                0, 0, 0,
                                2, 1, 0 };
  const size_t connectivity[] = { 0, 1, 3, 1, 2, 4 };
  const bool fixed[] = { false, false, false, true, true };
  tri_pd.fill( 5, vtx_coords, 2, TRIANGLE, connectivity, fixed, err );
  CPPUNIT_ASSERT(!err);
}

// tolerance on numerical gradient/hessian values
const double EPSILON = 5e-3;

/**\brief Fake quality metric for testing numerial gradient
 *
 * Implement a vertex metric for which the "quality" at a given
 * vertex is the value of the x-coordinate of that vertex.  Thus
 * the gradient of the "quality" at every vertex should be {1,0,0}.
 */
class LinearVertexMetric : public VertexQM
{
public:
  std::string get_name() const { return "Fake metric for testing numerical gradient"; }
  int get_negate_flag() const { return 1; }
  bool evaluate( PatchData& pd, size_t vtx_idx, double& value, MsqError& )
  {
    value = pd.vertex_by_index( vtx_idx )[0];
    return true;
  }
  bool evaluate_with_indices( PatchData& pd, size_t vtx_idx, double& value, std::vector<size_t>& indices, MsqError& )
  {
    value = pd.vertex_by_index( vtx_idx )[0];
    indices.resize(1);
    indices[0] = vtx_idx;
    return true;
  }
};

/**\brief Fake quality metric for testing numerial gradient
 *
 * Implement an element metric for which the "quality" is always 1.0.
 * The resulting gradient and Hessian should always be zero.
 */
class ConstantElementMetric : public ElementQM
{
public:
  std::string get_name() const { return "Fake metric for testing numerical gradient"; }
  int get_negate_flag() const { return 1; }
  bool evaluate( PatchData& pd, size_t , double& value, MsqError& )
    { value = 1.0; return true; }
  bool evaluate_with_indices( PatchData& pd, size_t elem_idx, double& value, std::vector<size_t>& indices, MsqError& )
  { 
    MsqMeshEntity& elem = pd.element_by_index( elem_idx );
    unsigned nv = elem.node_count();
    const size_t* conn = elem.get_vertex_index_array();
    indices.clear();
    for (unsigned i = 0; i < nv; ++i)
      if (conn[i] < pd.num_free_vertices())
        indices.push_back( conn[i] );

    value = 1.0; 
    return true; 
  }
};

/**\brief Fake quality metric for testing numerial gradient
 *
 * Implement a vertex metric for which the "quality" at a given
 * vertex is the square of the value of the y-coordinate of that vertex.  Thus
 * the Hessian of the "quality" at every vertex should be {2,0,0}.
 */
class ParabolicVertexMetric : public VertexQM
{
public:
  std::string get_name() const { return "Fake metric for testing numerical gradient"; }
  int get_negate_flag() const { return 1; }
  bool evaluate( PatchData& pd, size_t vtx_idx, double& value, MsqError& )
  {
    value = pd.vertex_by_index( vtx_idx )[1];
    value *= value;
    return true;
  }
  bool evaluate_with_indices( PatchData& pd, size_t vtx_idx, double& value, std::vector<size_t>& indices, MsqError& )
  {
    value = pd.vertex_by_index( vtx_idx )[1];
    value *= value;
    indices.resize(1);
    indices[0] = vtx_idx;
    return true;
  }
};

/**\brief More complex fake metric for testing finite difference 
 *
 * Calculate /f$ (det(M)-1)^2 /f$ where the rows of M are the
 * coordinates of the three vertices in the triangle.
 */
class TriTauMetric : public ElementQM
{
public:
  std::string get_name() const { return "Triangle Tau Metric"; }
  int get_negate_flag() const { return 1; }
  MsqMatrix<3,3> matrix( PatchData& pd, size_t elem_idx )
  {
    MsqMeshEntity& elem = pd.element_by_index( elem_idx );
    unsigned nv = elem.node_count();
    const size_t* conn = elem.get_vertex_index_array();
    CPPUNIT_ASSERT(nv == 3);
    
    MsqMatrix<3,3> M;
    M.set_row( 0, pd.vertex_by_index( conn[0] ).to_array() );
    M.set_row( 1, pd.vertex_by_index( conn[1] ).to_array() );
    M.set_row( 2, pd.vertex_by_index( conn[2] ).to_array() );
    
    return M;
  }
  
  bool evaluate( PatchData& pd, size_t element_idx, double& value, MsqError& )
  {
    MsqMatrix<3,3> M = matrix( pd, element_idx );
    value = det(M) - 1;
    value *= value;
    return true;
  }
  
  bool evaluate_with_indices( PatchData& pd, size_t elem_idx, double& value, std::vector<size_t>& indices, MsqError& )
  {
    MsqMatrix<3,3> M = matrix( pd, elem_idx );
    value = det(M) - 1;
    value *= value;
   
    indices.clear();
    unsigned nv = pd.element_by_index(elem_idx).node_count();
    const size_t* conn = pd.element_by_index(elem_idx).get_vertex_index_array();
    for (unsigned i = 0; i < nv; ++i) 
      if (pd.is_vertex_free( conn[i] ))
        indices.push_back( conn[i] );
    
    return true;
  }
};

void QualityMetricTest::test_fixed_vertex_list()
{
  // define a patch of four quads such that
  // the number of fixed vertices in each quad
  // varies from 0 to 3
  /*   2------1-----8
   *   |      |     |\
   *   | (0)  | (3) | \
   *   |      |     |  \
   *   3------0-----7   \
   *   |      |     |\   \
   *   | (1)  | (2) | \   \
   *   |      |     |  \   \
   *   4------5-----6   \   \
   *           \____\____\___\__ fixed
   */
  const double coords[] = { 0, 0, 0,
                            0, 1, 0,
                           -1, 1, 0,
                           -1, 0, 0,
                           -1,-1, 0,
                            0,-1, 0,
                            1,-1, 0,
                            1, 0, 0,
                            1, 1, 0 };
  const size_t conn[] = { 0, 1, 2, 3,
                          3, 4, 5, 0, 
                          6, 7, 0, 5,
                          0, 7, 8, 1 };
  const bool fixed[] = { false, false, false, 
                         false, false, true,
                         true,  true,  true };
  
  MsqPrintError err(std::cout);
  PatchData pd;
  pd.fill( 9, coords, 4, QUADRILATERAL, conn, fixed, err );
  ASSERT_NO_ERROR(err);
  
  uint32_t bits;
  std::vector<size_t> indices;
  std::vector<size_t>::iterator it;
  const size_t* verts;
  unsigned i;
  
  indices.clear();
  bits = QualityMetric::fixed_vertex_bitmap( pd, &pd.element_by_index(0), indices );
  CPPUNIT_ASSERT_EQUAL( (size_t)4, indices.size() );
  CPPUNIT_ASSERT_EQUAL( (uint32_t)0, bits&0xF );
  CPPUNIT_ASSERT( pd.num_free_vertices() > *std::max_element(indices.begin(), indices.end()) );
  verts = pd.element_by_index(0).get_vertex_index_array();
  for (i = 0; i < 4; ++i) 
    CPPUNIT_ASSERT( std::find( verts, verts+4, indices[i] ) != verts+4 );
  
  indices.clear();
  bits = QualityMetric::fixed_vertex_bitmap( pd, &pd.element_by_index(1), indices );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, indices.size() );
  verts = pd.element_by_index(1).get_vertex_index_array();
  for (i = 0; i < 4; ++i) {
    it = std::find( indices.begin(), indices.end(), verts[i] );
    if (verts[i] < pd.num_free_vertices()) {
      CPPUNIT_ASSERT( it != indices.end() );
      CPPUNIT_ASSERT_EQUAL( 0u, bits & (1<<i) );
    }
    else {
      CPPUNIT_ASSERT( it == indices.end() );
      CPPUNIT_ASSERT_EQUAL( 1u, (bits>>i) & 1 );
    }
  }
  
  indices.clear();
  bits = QualityMetric::fixed_vertex_bitmap( pd, &pd.element_by_index(2), indices );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, indices.size() );
  verts = pd.element_by_index(2).get_vertex_index_array();
  for (i = 0; i < 4; ++i) {
    it = std::find( indices.begin(), indices.end(), verts[i] );
    if (verts[i] < pd.num_free_vertices()) {
      CPPUNIT_ASSERT( it != indices.end() );
      CPPUNIT_ASSERT_EQUAL( 0u, bits & (1<<i) );
    }
    else {
      CPPUNIT_ASSERT( it == indices.end() );
      CPPUNIT_ASSERT_EQUAL( 1u, (bits>>i) & 1 );
    }
  }
  
  indices.clear();
  bits = QualityMetric::fixed_vertex_bitmap( pd, &pd.element_by_index(3), indices );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, indices.size() );
  verts = pd.element_by_index(3).get_vertex_index_array();
  for (i = 0; i < 4; ++i) {
    it = std::find( indices.begin(), indices.end(), verts[i] );
    if (verts[i] < pd.num_free_vertices()) {
      CPPUNIT_ASSERT( it != indices.end() );
      CPPUNIT_ASSERT_EQUAL( 0u, bits & (1<<i) );
    }
    else {
      CPPUNIT_ASSERT( it == indices.end() );
      CPPUNIT_ASSERT_EQUAL( 1u, (bits>>i) & 1 );
    }
  }
}  

void QualityMetricTest::test_remove_fixed_gradients()
{
    // define a list of vectors
  std::vector<Vector3D> grads(6);
  grads[0] = Vector3D(0,0,0);
  grads[1] = Vector3D(1,1,1);
  grads[2] = Vector3D(2,2,2);
  grads[3] = Vector3D(3,3,3);
  grads[4] = Vector3D(4,4,4);
  grads[5] = Vector3D(5,5,5);
    // remove the first, third, and last
  uint32_t flags = 1u | 4u | 32u;
    // call function, choose element type w/ correct number of corners
  QualityMetric::remove_fixed_gradients( PRISM, flags, grads );
    // check results
  CPPUNIT_ASSERT_EQUAL( (size_t)3, grads.size() );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(1,1,1), grads[0], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(3,3,3), grads[1], DBL_EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(4,4,4), grads[2], DBL_EPSILON );
}

void QualityMetricTest::test_remove_fixed_hessians()
{
    // define Hessian matrix for a quadrilateral
  Matrix3D m[10] = {
    Matrix3D(0.0), Matrix3D(1.0), Matrix3D(2.0), Matrix3D(3.0),
                   Matrix3D(4.0), Matrix3D(5.0), Matrix3D(6.0),
                                  Matrix3D(7.0), Matrix3D(8.0),
                                                 Matrix3D(9.0)
  };
    // convert to std::vector
  std::vector<Matrix3D> hess(10);
  std::copy( m, m+10, hess.begin() );
    // mark fist and third vertices as fixed
  uint32_t flags = 1u | 4u;
    // call function to remove grads for fixed vertices
  QualityMetric::remove_fixed_hessians( QUADRILATERAL, flags, hess );
    // the submatrix with the first and third rows/columns removed should be
    // { 4, 6,
    //      9 }
  CPPUNIT_ASSERT_EQUAL( (size_t)3, hess.size() );
  CPPUNIT_ASSERT_MATRICES_EQUAL( Matrix3D(4.0), hess[0], DBL_EPSILON );
  CPPUNIT_ASSERT_MATRICES_EQUAL( Matrix3D(6.0), hess[1], DBL_EPSILON );
  CPPUNIT_ASSERT_MATRICES_EQUAL( Matrix3D(9.0), hess[2], DBL_EPSILON );
}


void QualityMetricTest::test_gradient_constant()
{
  MsqError err;
  std::vector<size_t> indices;
  std::vector<Vector3D> gradient;
  double value;
  size_t ELEMENT = 0;

  ConstantElementMetric constant;
  QualityMetric& qm = constant;

  bool valid = qm.evaluate_with_gradient( tri_pd, ELEMENT, value, indices, gradient, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(valid);
  compare_indices( qm, tri_pd, ELEMENT, value, indices );
  CPPUNIT_ASSERT_EQUAL(indices.size(), gradient.size());
  CPPUNIT_ASSERT_EQUAL((size_t)2, indices.size());  // two free vertices in triangle
  
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,0,0), gradient[0], EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,0,0), gradient[1], EPSILON );
}


void QualityMetricTest::test_gradient_linear()
{
  MsqError err;
  std::vector<size_t> indices;
  std::vector<Vector3D> gradient;
  double value;
  size_t VERTEX = 0;
  
  LinearVertexMetric linear;
  QualityMetric& qm = linear;

  bool valid = qm.evaluate_with_gradient( tri_pd, VERTEX, value, indices, gradient, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(valid);
  compare_indices( qm, tri_pd, VERTEX, value, indices );
  CPPUNIT_ASSERT_EQUAL(indices.size(), gradient.size());
  CPPUNIT_ASSERT_EQUAL((size_t)1, indices.size());
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(1,0,0), gradient[0], EPSILON );
}


void QualityMetricTest::test_gradient_parabolic()
{
  const size_t VERTEX = 1;  // pick vertex with non-zero Y-coordinate
  MsqError err;
  std::vector<size_t> indices;
  std::vector<Vector3D> gradient;
  double value;

  ParabolicVertexMetric parab;
  QualityMetric& qm = parab;

  bool valid = qm.evaluate_with_gradient( tri_pd, VERTEX, value, indices, gradient, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(valid);
  compare_indices( qm, tri_pd, VERTEX, value, indices );
  CPPUNIT_ASSERT_EQUAL(indices.size(), gradient.size());
  CPPUNIT_ASSERT_EQUAL((size_t)1, indices.size());  // two free vertices in triangle
  
  const double expected_y = 2*tri_pd.vertex_by_index(VERTEX)[1];
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,expected_y,0), gradient[0], EPSILON );
}


void QualityMetricTest::test_gradient_tau()
{
  MsqError err;
  std::vector<size_t> indices;
  std::vector<Vector3D> gradient;
  double value;
  size_t ELEMENT = 1;

  TriTauMetric qm;

  bool valid = qm.evaluate_with_gradient( tri_pd, ELEMENT, value, indices, gradient, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(valid);
  compare_indices( qm, tri_pd, ELEMENT, value, indices );
  CPPUNIT_ASSERT_EQUAL(indices.size(), gradient.size());
  CPPUNIT_ASSERT_EQUAL((size_t)2, indices.size());  // two free vertices in triangle
  
  MsqMatrix<3,3> M = qm.matrix( tri_pd, ELEMENT );
  MsqMatrix<3,3> dM = transpose(adj(M));
  dM *= 2 * (det(M) - 1);
  
  if (indices[0] == 2) {
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices[1] );
    std::swap( indices[0], indices[1] );
    std::swap( gradient[0], gradient[1] );
  }
  else {
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices[0] );
    CPPUNIT_ASSERT_EQUAL( (size_t)2, indices[1] );
  }  
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(dM.row(0).data()), gradient[0], EPSILON );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(dM.row(1).data()), gradient[1], EPSILON );
}

void QualityMetricTest::test_Hessian_constant()
{
  MsqError err;
  std::vector<size_t> indices;
  std::vector<Vector3D> gradient;
  std::vector<Matrix3D> Hessian;
  double value;
  size_t ELEMENT = 0;

  ConstantElementMetric constant;
  QualityMetric& qm = constant;
  
  bool valid = qm.evaluate_with_Hessian( tri_pd, ELEMENT, value, indices, gradient, Hessian, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(valid);
  compare_gradient( qm, tri_pd, ELEMENT, value, indices, gradient );
  CPPUNIT_ASSERT_EQUAL(indices.size()*(indices.size()+1)/2, Hessian.size());
  CPPUNIT_ASSERT_EQUAL((size_t)2, indices.size());
  
  CPPUNIT_ASSERT_MATRICES_EQUAL( Matrix3D(0.0), Hessian[0], EPSILON );
  CPPUNIT_ASSERT_MATRICES_EQUAL( Matrix3D(0.0), Hessian[1], EPSILON );
}


void QualityMetricTest::test_Hessian_linear()
{
  MsqError err;
  std::vector<size_t> indices;
  std::vector<Vector3D> gradient;
  std::vector<Matrix3D> Hessian;
  double value;
  size_t VERTEX = 0;
  
  LinearVertexMetric linear;
  QualityMetric& qm = linear;

  bool valid = qm.evaluate_with_Hessian( tri_pd, VERTEX, value, indices, gradient, Hessian, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(valid);
  compare_gradient( qm, tri_pd, VERTEX, value, indices, gradient );
  CPPUNIT_ASSERT_EQUAL(indices.size()*(indices.size()+1)/2, Hessian.size());
  CPPUNIT_ASSERT_EQUAL((size_t)1, indices.size());
  CPPUNIT_ASSERT_MATRICES_EQUAL( Matrix3D(0.0), Hessian[0], EPSILON );
}


void QualityMetricTest::test_Hessian_parabolic()
{
  MsqError err;
  std::vector<size_t> indices;
  std::vector<Vector3D> gradient;
  std::vector<Matrix3D> Hessian;
  double value;
  size_t VERTEX = 0;

  ParabolicVertexMetric parab;
  QualityMetric& qm = parab;
  
  bool valid = qm.evaluate_with_Hessian( tri_pd, VERTEX, value, indices, gradient, Hessian, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(valid);
  compare_gradient( qm, tri_pd, VERTEX, value, indices, gradient );
  CPPUNIT_ASSERT_EQUAL(indices.size()*(indices.size()+1)/2, Hessian.size());
  CPPUNIT_ASSERT_EQUAL((size_t)1, indices.size());
  
  Matrix3D expected( 0.0, 0.0, 0.0,
                     0.0, 2.0, 0.0,
                     0.0, 0.0, 0.0 );
  CPPUNIT_ASSERT_MATRICES_EQUAL( expected, Hessian[0], EPSILON );
}


void QualityMetricTest::test_Hessian_tau()
{
  MsqError err;
  std::vector<size_t> indices;
  std::vector<Vector3D> gradient;
  std::vector<Matrix3D> Hessian;
  double value;
  size_t ELEMENT = 1;

  TriTauMetric qm;

  bool valid = qm.evaluate_with_Hessian( tri_pd, ELEMENT, value, indices, gradient, Hessian, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(valid);
  compare_gradient( qm, tri_pd, ELEMENT, value, indices, gradient );
  CPPUNIT_ASSERT_EQUAL(indices.size()*(indices.size()+1)/2, Hessian.size());
  CPPUNIT_ASSERT_EQUAL((size_t)2, indices.size());  // two free vertices in triangle
  
  MsqMatrix<3,3> M = qm.matrix( tri_pd, ELEMENT );
  MsqMatrix<3,3> d2M[6];
  set_scaled_2nd_deriv_of_det( d2M, 2 * (det(M) - 1), M );
  pluseq_scaled_outer_product( d2M, 2, transpose(adj(M)) );
  
  if (indices[0] == 2) {
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices[1] );
    std::swap( indices[0], indices[1] );
    std::swap( gradient[0], gradient[1] );
    std::swap( Hessian[0], Hessian[2] );
    Hessian[1] = transpose(Hessian[1]);
  }
  else {
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices[0] );
    CPPUNIT_ASSERT_EQUAL( (size_t)2, indices[1] );
  }  
  
  CPPUNIT_ASSERT_MATRICES_EQUAL( Matrix3D(d2M[0].data()), Hessian[0], EPSILON );
  CPPUNIT_ASSERT_MATRICES_EQUAL( Matrix3D(d2M[1].data()), Hessian[1], EPSILON );
  CPPUNIT_ASSERT_MATRICES_EQUAL( Matrix3D(d2M[3].data()), Hessian[2], EPSILON );
}

  
void QualityMetricTest::test_diagonal_constant()
{
  MsqError err;
  std::vector<size_t> indices;
  std::vector<Vector3D> gradient;
  std::vector<SymMatrix3D> Hessian;
  double value;
  size_t ELEMENT = 0;

  ConstantElementMetric constant;
  QualityMetric& qm = constant;
  
  bool valid = qm.evaluate_with_Hessian_diagonal( tri_pd, ELEMENT, value, indices, gradient, Hessian, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(valid);
  compare_gradient( qm, tri_pd, ELEMENT, value, indices, gradient );
  CPPUNIT_ASSERT_EQUAL(indices.size(), Hessian.size());
  CPPUNIT_ASSERT_EQUAL((size_t)2, indices.size());
  
  CPPUNIT_ASSERT_MATRICES_EQUAL( SymMatrix3D(0.0), Hessian[0], EPSILON );
  CPPUNIT_ASSERT_MATRICES_EQUAL( SymMatrix3D(0.0), Hessian[1], EPSILON );
}


void QualityMetricTest::test_diagonal_linear()
{
  MsqError err;
  std::vector<size_t> indices;
  std::vector<Vector3D> gradient;
  std::vector<SymMatrix3D> Hessian;
  double value;
  size_t VERTEX = 0;
  
  LinearVertexMetric linear;
  QualityMetric& qm = linear;

  bool valid = qm.evaluate_with_Hessian_diagonal( tri_pd, VERTEX, value, indices, gradient, Hessian, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(valid);
  compare_gradient( qm, tri_pd, VERTEX, value, indices, gradient );
  CPPUNIT_ASSERT_EQUAL(indices.size(), Hessian.size());
  CPPUNIT_ASSERT_EQUAL((size_t)1, indices.size());
  CPPUNIT_ASSERT_MATRICES_EQUAL( SymMatrix3D(0.0), Hessian[0], EPSILON );
}


void QualityMetricTest::test_diagonal_parabolic()
{
  MsqError err;
  std::vector<size_t> indices;
  std::vector<Vector3D> gradient;
  std::vector<SymMatrix3D> Hessian;
  double value;
  size_t VERTEX = 0;

  ParabolicVertexMetric parab;
  QualityMetric& qm = parab;
  
  bool valid = qm.evaluate_with_Hessian_diagonal( tri_pd, VERTEX, value, indices, gradient, Hessian, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(valid);
  compare_gradient( qm, tri_pd, VERTEX, value, indices, gradient );
  CPPUNIT_ASSERT_EQUAL(indices.size(), Hessian.size());
  CPPUNIT_ASSERT_EQUAL((size_t)1, indices.size());
  
  SymMatrix3D expected( 0.0, 0.0, 0.0,
                             2.0, 0.0,
                                  0.0 );
  CPPUNIT_ASSERT_MATRICES_EQUAL( expected, Hessian[0], EPSILON );
}


void QualityMetricTest::test_diagonal_tau()
{
  MsqError err;
  std::vector<size_t> indices;
  std::vector<Vector3D> gradient;
  std::vector<SymMatrix3D> Hessian;
  double value;
  size_t ELEMENT = 1;

  TriTauMetric qm;

  bool valid = qm.evaluate_with_Hessian_diagonal( tri_pd, ELEMENT, value, indices, gradient, Hessian, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(valid);
  compare_gradient( qm, tri_pd, ELEMENT, value, indices, gradient );
  CPPUNIT_ASSERT_EQUAL(indices.size(), Hessian.size());
  CPPUNIT_ASSERT_EQUAL((size_t)2, indices.size());  // two free vertices in triangle
  
  MsqMatrix<3,3> M = qm.matrix( tri_pd, ELEMENT );
  MsqMatrix<3,3> d2M[6];
  set_scaled_2nd_deriv_of_det( d2M, 2 * (det(M) - 1), M );
  pluseq_scaled_outer_product( d2M, 2, transpose(adj(M)) );
  
  if (indices[0] == 2) {
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices[1] );
    std::swap( indices[0], indices[1] );
    std::swap( gradient[0], gradient[1] );
    std::swap( Hessian[0], Hessian[1] );
  }
  else {
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices[0] );
    CPPUNIT_ASSERT_EQUAL( (size_t)2, indices[1] );
  }  
  
  SymMatrix3D exp0( d2M[0](0,0), d2M[0](0,1), d2M[0](0,2),
                                 d2M[0](1,1), d2M[0](1,2),
                                              d2M[0](2,2) );
  SymMatrix3D exp1( d2M[3](0,0), d2M[3](0,1), d2M[3](0,2),
                                 d2M[3](1,1), d2M[3](1,2),
                                              d2M[3](2,2) );
  
  CPPUNIT_ASSERT_MATRICES_EQUAL( exp0, Hessian[0], EPSILON );
  CPPUNIT_ASSERT_MATRICES_EQUAL( exp1, Hessian[1], EPSILON );
}

void QualityMetricTest::compare_indices( QualityMetric& qm, PatchData& pd, 
                                         size_t sample, double value,
                                         const std::vector<size_t>& indices )
{
  double value2;
  std::vector<size_t> indices2;
  MsqError err;
  bool valid = qm.evaluate_with_indices( pd, sample, value2, indices2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( value2, value, EPSILON );
  std::vector<size_t> indices1(indices);
  std::sort( indices1.begin(), indices1.end() );
  std::sort( indices2.begin(), indices2.end() );
  ASSERT_STD_VECTORS_EQUAL( indices2, indices1 );
}

void QualityMetricTest::compare_gradient( QualityMetric& qm, PatchData& pd, 
                                          size_t sample, double value,
                                          const std::vector<size_t>& indices,
                                          const std::vector<Vector3D>& grad )
{
  double value2;
  std::vector<size_t> indices2;
  std::vector<Vector3D> grad2;
  MsqError err;
  bool valid = qm.evaluate_with_gradient( pd, sample, value2, indices2, grad2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( value2, value, EPSILON );
  CPPUNIT_ASSERT_EQUAL( indices2.size(), indices.size() );
  CPPUNIT_ASSERT_EQUAL( indices.size(), grad.size() );
  
  for (size_t i = 0; i < indices.size(); ++i) {
    size_t j = std::find( indices2.begin(), indices2.end(), indices[i] ) - indices2.begin();
    CPPUNIT_ASSERT( j < indices2.size() ); // not found->indices don't match
    CPPUNIT_ASSERT_VECTORS_EQUAL( grad2[j], grad[i], EPSILON );
  }
}
