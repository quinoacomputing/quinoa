/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TQualityMetricTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Mesquite_TQualityMetric.hpp"
#include "Mesquite_TInverseMeanRatio.hpp"
#include "Mesquite_IdealWeightInverseMeanRatio.hpp"
#include "TMPQualityMetricTest.hpp"
#include "Mesquite_TShapeNB1.hpp"

template <> class TMPTypes<TQualityMetric> {
public:
  typedef TMetric MetricType;
  typedef TShapeNB1 TestType;
};

class TQualityMetricTest : public TMPQualityMetricTest<TQualityMetric>
{
  CPPUNIT_TEST_SUITE(TQualityMetricTest);
  
  REGISTER_TMP_TESTS
  
  CPPUNIT_TEST (test_inverse_mean_ratio_grad);
  CPPUNIT_TEST (test_inverse_mean_ratio_hess);
  CPPUNIT_TEST (test_inverse_mean_ratio_hess_diag);
  CPPUNIT_TEST (regression_inverse_mean_ratio_grad);
  CPPUNIT_TEST (regression_inverse_mean_ratio_hess);
  
  CPPUNIT_TEST_SUITE_END();
  
public:
  
  void test_inverse_mean_ratio_grad();
  void test_inverse_mean_ratio_hess();
  void test_inverse_mean_ratio_hess_diag();
  void regression_inverse_mean_ratio_grad();
  void regression_inverse_mean_ratio_hess();

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TQualityMetricTest, "TQualityMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TQualityMetricTest, "Unit");

  
void TQualityMetricTest::test_inverse_mean_ratio_grad()
{
  TInverseMeanRatio tm;
  IdealShapeTarget target;
  TQualityMetric metric( &target, &tm );
  ElementPMeanP avg( 1.0, &metric );
  
  tester.test_gradient_reflects_quality( &metric );
  compare_analytical_and_numerical_gradients( &metric );
  tester.test_gradient_with_fixed_vertex( &avg );
}


void TQualityMetricTest::test_inverse_mean_ratio_hess()
{
  TInverseMeanRatio tm;
  IdealShapeTarget target;
  TQualityMetric metric( &target, &tm );
  ElementPMeanP avg( 1.0, &metric );
 
  compare_analytical_and_numerical_hessians( &metric );
  tester.test_symmetric_Hessian_diagonal_blocks( &metric );
  tester.test_hessian_with_fixed_vertex( &avg );
}

void TQualityMetricTest::test_inverse_mean_ratio_hess_diag()
{
  TInverseMeanRatio tm;
  IdealShapeTarget target;
  TQualityMetric metric( &target, &tm );
  
  compare_analytical_and_numerical_diagonals( &metric );
  tester.compare_eval_with_diag_and_eval_with_hessian( &metric );
}
  
void TQualityMetricTest::regression_inverse_mean_ratio_grad()
{
  MsqError err;
  TInverseMeanRatio tm;
  IdealShapeTarget target;
  TQualityMetric metric( &target, &tm );
  const double coords[] = { -0.80000000000000004, -0.80000000000000004, 0,
                             0.00000000000000000,  2.00000000000000000, 0,
                            -1.73205079999999990,  1.00000000000000000, 0 };
  const size_t indices[] = { 0, 1, 2 };
  PatchData pd;
  pd.fill( 3, coords, 1, TRIANGLE, indices, 0, err );
  pd.attach_settings( &settings );
  PlanarDomain dom( PlanarDomain::XY, coords[0] );
  pd.set_domain( &dom );
  
  IdealWeightInverseMeanRatio ref_metric;
  
  double exp_val, act_val;
  std::vector<size_t> exp_idx, act_idx, handles;
  std::vector<Vector3D> exp_grad, act_grad;
  
  handles.clear();
  ref_metric.get_evaluations( pd, handles, false, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( (size_t)1, handles.size() );
  const size_t hand1 = handles.front();
  handles.clear();
  metric.get_evaluations( pd, handles, false, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( (size_t)1, handles.size() );
  const size_t hand2 = handles.front();
  
  bool exp_rval, act_rval;
  exp_rval = ref_metric.evaluate_with_gradient( pd, hand1, exp_val, exp_idx, exp_grad, err );
  ASSERT_NO_ERROR(err);
  act_rval = metric.evaluate_with_gradient( pd, hand2, act_val, act_idx, act_grad, err );
  ASSERT_NO_ERROR(err);
  
  CPPUNIT_ASSERT( exp_rval );
  CPPUNIT_ASSERT( act_rval );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_val - 1.0, act_val, 1e-5 );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, exp_idx.size() );
  CPPUNIT_ASSERT_EQUAL( (size_t)3, act_idx.size() );
  
  std::vector<size_t> sorted(exp_idx);
  std::sort( sorted.begin(), sorted.end() );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, sorted[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, sorted[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, sorted[2] );
  
  sorted = act_idx;
  std::sort( sorted.begin(), sorted.end() );
  CPPUNIT_ASSERT_EQUAL( (size_t)0, sorted[0] );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, sorted[1] );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, sorted[2] );
  
  const size_t idx_map[] = { 
    std::find(act_idx.begin(),act_idx.end(),exp_idx[0]) - act_idx.begin(),
    std::find(act_idx.begin(),act_idx.end(),exp_idx[1]) - act_idx.begin(),
    std::find(act_idx.begin(),act_idx.end(),exp_idx[2]) - act_idx.begin() };
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[0], act_grad[idx_map[0]], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[1], act_grad[idx_map[1]], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[2], act_grad[idx_map[2]], 1e-5 );
}

  
void TQualityMetricTest::regression_inverse_mean_ratio_hess()
{
  MsqError err;
  TInverseMeanRatio tm;
  IdealShapeTarget target;
  TQualityMetric metric( &target, &tm );
  const double coords[] = { 4.158984727, 4.6570859130000004, 5,
                            4.51742825, 4.51742825, 5,
                            4.3103448279999999, 5, 5 };
  const bool fixed[] = { false, false, true };
  const size_t indices[] = { 0, 1, 2 };
  PatchData pd;
  pd.fill( 3, coords, 1, TRIANGLE, indices, fixed, err );
  pd.attach_settings( &settings );
  PlanarDomain dom( PlanarDomain::XY, coords[2] );
  pd.set_domain( &dom );
  
  IdealWeightInverseMeanRatio ref_metric;
  
  double exp_val, act_val;
  std::vector<size_t> exp_idx, act_idx, handles;
  std::vector<Vector3D> exp_grad, act_grad;
  std::vector<Matrix3D> exp_hess, act_hess;
  
  handles.clear();
  ref_metric.get_evaluations( pd, handles, false, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( (size_t)1, handles.size() );
  const size_t hand1 = handles.front();
  handles.clear();
  metric.get_evaluations( pd, handles, false, err ); ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( (size_t)1, handles.size() );
  const size_t hand2 = handles.front();
  
    // first make sure that non-TMP IMR metric works correctly
  bool exp_rval, act_rval;
  exp_rval = ref_metric.evaluate_with_Hessian( pd, hand1, exp_val, exp_idx, exp_grad, exp_hess, err );
  ASSERT_NO_ERROR(err);
  act_rval = ref_metric.QualityMetric::evaluate_with_Hessian( pd, hand1, act_val, act_idx, act_grad, act_hess, err );
  CPPUNIT_ASSERT( exp_rval );
  CPPUNIT_ASSERT( act_rval );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_val, act_val, 1e-5 );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, exp_idx.size() );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, act_idx.size() );

  if (act_idx[0] == exp_idx[0]) {
    CPPUNIT_ASSERT_EQUAL( exp_idx[1], act_idx[1] );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[0], act_grad[0], 1e-5 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[1], act_grad[1], 1e-5 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[0], act_hess[0], 5e-3 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[1], act_hess[1], 5e-3 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[2], act_hess[2], 5e-3 );
  }
  else {
    CPPUNIT_ASSERT_EQUAL( exp_idx[0], act_idx[1] );
    CPPUNIT_ASSERT_EQUAL( exp_idx[1], act_idx[0] );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[0], act_grad[1], 1e-5 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[1], act_grad[0], 1e-5 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[0], act_hess[2], 5e-3 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[1], transpose(act_hess[1]), 5e-3 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[2], act_hess[0], 5e-3 );
  }
  
    // now compare TMP metric with non-TMP metric
  act_rval = metric.evaluate_with_Hessian( pd, hand2, act_val, act_idx, act_grad, act_hess, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT( act_rval );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_val - 1.0, act_val, 1e-5 );
  CPPUNIT_ASSERT_EQUAL( (size_t)2, act_idx.size() );

#ifdef PLANAR_HESSIAN
  // zero derivatives with respect to Z
  for (int i = 0; i < 2; ++i) 
    exp_grad[i][2] = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) 
      exp_hess[i][j][2] = exp_hess[i][2][j] = 0.0;
  }
#else
  // don't compare double out-of-plane term because metrics
  // make varying assumptions about behavior of Hessian
  for (int i = 0; i < 3; ++i)
    exp_hess[i][2][2] = act_hess[i][2][2] = 0.0;
#endif
  
  if (act_idx[0] == exp_idx[0]) {
    CPPUNIT_ASSERT_EQUAL( exp_idx[1], act_idx[1] );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[0], act_grad[0], 1e-5 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[1], act_grad[1], 1e-5 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[0], act_hess[0], 5e-3 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[1], act_hess[1], 5e-3 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[2], act_hess[2], 5e-3 );
  }
  else {
    CPPUNIT_ASSERT_EQUAL( exp_idx[0], act_idx[1] );
    CPPUNIT_ASSERT_EQUAL( exp_idx[1], act_idx[0] );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[0], act_grad[1], 1e-5 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[1], act_grad[0], 1e-5 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[0], act_hess[2], 5e-3 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[1], transpose(act_hess[1]), 5e-3 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[2], act_hess[0], 5e-3 );
  }
}

