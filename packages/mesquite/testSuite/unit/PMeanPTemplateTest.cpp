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


/** \file PMeanPTemplateTest.cpp
 *  \brief previous name: PowerMeanPTest.cpp
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Mesquite_PMeanPTemplate.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_PatchData.hpp"
#include "UnitUtil.hpp"
#include "Mesquite_VertexQM.hpp"
#include "Mesquite_MsqHessian.hpp"

#include "ObjectiveFunctionTests.hpp"


using namespace Mesquite;
using namespace std;

const double EPSILON = 1e-4;

class PMeanPTemplateTest : public CppUnit::TestFixture, public ObjectiveFunctionTests
{
private:
  void test_eval_type( OFTestMode eval_func, ObjectiveFunction::EvalType type );
  void test_evaluate( double power );
  void test_gradient( double power );
  void test_diagonal( double power );
  void test_Hessian( double power );
  
  void check_result( PatchData& pd, double power, double value, 
                     Vector3D* gradient = 0, Matrix3D* Hessian = 0 );

  CPPUNIT_TEST_SUITE( PMeanPTemplateTest );

  CPPUNIT_TEST( test_eval_calc );
  CPPUNIT_TEST( test_eval_accum );
  CPPUNIT_TEST( test_eval_save );
  CPPUNIT_TEST( test_eval_update );
  CPPUNIT_TEST( test_eval_temp );

  CPPUNIT_TEST( test_grad_calc );
  CPPUNIT_TEST( test_grad_save );
  CPPUNIT_TEST( test_grad_update );
  CPPUNIT_TEST( test_grad_temp );

  CPPUNIT_TEST( test_diag_calc );
  CPPUNIT_TEST( test_diag_save );
  CPPUNIT_TEST( test_diag_update );
  CPPUNIT_TEST( test_diag_temp );

  CPPUNIT_TEST( test_Hess_calc );
  CPPUNIT_TEST( test_Hess_save );
  CPPUNIT_TEST( test_Hess_update );
  CPPUNIT_TEST( test_Hess_temp );
  
  CPPUNIT_TEST( test_clone );
  
  CPPUNIT_TEST( test_failed_metric_in_eval );
  CPPUNIT_TEST( test_failed_metric_in_grad );
  CPPUNIT_TEST( test_failed_metric_in_diag );
  CPPUNIT_TEST( test_failed_metric_in_Hess );
  
  CPPUNIT_TEST( test_false_metric_in_eval );
  CPPUNIT_TEST( test_false_metric_in_grad );
  CPPUNIT_TEST( test_false_metric_in_diag );
  CPPUNIT_TEST( test_false_metric_in_Hess );
  
  CPPUNIT_TEST( test_evaluate_arithmatic );
  CPPUNIT_TEST( test_evaluate_rms );
  
  CPPUNIT_TEST( test_gradient_arithmatic );
  CPPUNIT_TEST( test_gradient_rms );

  CPPUNIT_TEST( test_Hessian_arithmatic );
  CPPUNIT_TEST( test_Hessian_rms );
  
  CPPUNIT_TEST( compare_gradient_arithmatic );
  CPPUNIT_TEST( compare_gradient_rms );
 
  CPPUNIT_TEST( compare_diagonal_gradient_arithmatic );
  CPPUNIT_TEST( compare_diagonal_gradient_rms );
 
  CPPUNIT_TEST( compare_hessian_gradient_arithmatic );
  CPPUNIT_TEST( compare_hessian_gradient_rms );
 
  CPPUNIT_TEST( compare_hessian_diagonal_arithmatic );
  CPPUNIT_TEST( compare_hessian_diagonal_rms );
  
  CPPUNIT_TEST( compare_hessian_arithmatic );
  CPPUNIT_TEST( compare_hessian_rms );
  
  CPPUNIT_TEST( compare_hessian_diag_arithmatic );
  CPPUNIT_TEST( compare_hessian_diag_rms );
  
  CPPUNIT_TEST( test_negate_eval );
  CPPUNIT_TEST( test_negate_grad );
  CPPUNIT_TEST( test_negate_diag );
  CPPUNIT_TEST( test_negate_hess );

  CPPUNIT_TEST_SUITE_END();
  PatchData mPatch;

public:

  void setUp();

  void test_eval_calc()   { test_eval_type( EVAL, ObjectiveFunction::CALCULATE ); }
  void test_eval_accum()  { test_eval_type( EVAL, ObjectiveFunction::ACCUMULATE ); }
  void test_eval_save()   { test_eval_type( EVAL, ObjectiveFunction::SAVE ); }
  void test_eval_update() { test_eval_type( EVAL, ObjectiveFunction::UPDATE ); }
  void test_eval_temp()   { test_eval_type( EVAL, ObjectiveFunction::TEMPORARY ); }

  void test_grad_calc()   { test_eval_type( GRAD, ObjectiveFunction::CALCULATE ); }
  void test_grad_save()   { test_eval_type( GRAD, ObjectiveFunction::SAVE ); }
  void test_grad_update() { test_eval_type( GRAD, ObjectiveFunction::UPDATE ); }
  void test_grad_temp()   { test_eval_type( GRAD, ObjectiveFunction::TEMPORARY ); }

  void test_diag_calc()   { test_eval_type( DIAG, ObjectiveFunction::CALCULATE ); }
  void test_diag_save()   { test_eval_type( DIAG, ObjectiveFunction::SAVE ); }
  void test_diag_update() { test_eval_type( DIAG, ObjectiveFunction::UPDATE ); }
  void test_diag_temp()   { test_eval_type( DIAG, ObjectiveFunction::TEMPORARY ); }

  void test_Hess_calc()   { test_eval_type( HESS, ObjectiveFunction::CALCULATE ); }
  void test_Hess_save()   { test_eval_type( HESS, ObjectiveFunction::SAVE ); }
  void test_Hess_update() { test_eval_type( HESS, ObjectiveFunction::UPDATE ); }
  void test_Hess_temp()   { test_eval_type( HESS, ObjectiveFunction::TEMPORARY ); }

  void test_clone() { PMeanPTemplate of( 1, NULL ); 
                      ObjectiveFunctionTests::test_clone(&of); }
  
  void test_failed_metric_in_eval() 
    { PMeanPTemplate of( 1, NULL ); test_handles_qm_error( EVAL, &of); }
  void test_failed_metric_in_grad() 
    { PMeanPTemplate of( 1, NULL ); test_handles_qm_error( GRAD, &of); }
  void test_failed_metric_in_diag() 
    { PMeanPTemplate of( 1, NULL ); test_handles_qm_error( DIAG, &of); }
  void test_failed_metric_in_Hess() 
    { PMeanPTemplate of( 1, NULL ); test_handles_qm_error( HESS, &of); }
  
  void test_false_metric_in_eval() 
    { PMeanPTemplate of( 1, NULL ); test_handles_invalid_qm( EVAL, &of); }
  void test_false_metric_in_grad() 
    { PMeanPTemplate of( 1, NULL ); test_handles_invalid_qm( GRAD, &of); }
  void test_false_metric_in_diag() 
    { PMeanPTemplate of( 1, NULL ); test_handles_invalid_qm( DIAG, &of); }
  void test_false_metric_in_Hess() 
    { PMeanPTemplate of( 1, NULL ); test_handles_invalid_qm( HESS, &of); }
  
  void test_evaluate_arithmatic() { test_evaluate( 1 ); }
  void test_evaluate_rms()        { test_evaluate( 2 ); }
  
  void test_gradient_arithmatic() { test_gradient( 1 ); }
  void test_gradient_rms()        { test_gradient( 2 ); }
  
  void test_diagonal_arithmatic() { test_diagonal( 1 ); }
  void test_diagonal_rms()        { test_diagonal( 2 ); }

  void test_Hessian_arithmatic()  { test_Hessian( 1 ); }
  void test_Hessian_rms()         { test_Hessian( 2 ); }
  
  void compare_gradient_arithmatic() 
    { PMeanPTemplate of( 1, NULL ); compare_numerical_gradient( &of ); }
  void compare_gradient_rms()
    { PMeanPTemplate of( 2, NULL ); compare_numerical_gradient( &of ); }
  
  void compare_diagonal_gradient_arithmatic() 
    { PMeanPTemplate of( 1, NULL ); compare_diagonal_gradient( &of ); }
  void compare_diagonal_gradient_rms()
    { PMeanPTemplate of( 2, NULL ); compare_diagonal_gradient( &of ); }
  
  void compare_hessian_gradient_arithmatic() 
    { PMeanPTemplate of( 1, NULL ); compare_hessian_gradient( &of ); }
  void compare_hessian_gradient_rms()
    { PMeanPTemplate of( 2, NULL ); compare_hessian_gradient( &of ); }
  
  void compare_hessian_diagonal_arithmatic() 
    { PMeanPTemplate of( 1, NULL ); compare_hessian_diagonal( &of ); }
  void compare_hessian_diagonal_rms()
    { PMeanPTemplate of( 2, NULL ); compare_hessian_diagonal( &of ); }
  
  void compare_hessian_arithmatic() 
    { PMeanPTemplate of( 1, NULL ); compare_numerical_hessian( &of ); }
  void compare_hessian_rms()
    { PMeanPTemplate of( 2, NULL ); compare_numerical_hessian( &of ); }
  
  void compare_hessian_diag_arithmatic() 
    { PMeanPTemplate of( 1, NULL ); compare_numerical_hessian_diagonal( &of ); }
  void compare_hessian_diag_rms()
    { PMeanPTemplate of( 2, NULL ); compare_numerical_hessian_diagonal( &of ); }
    
  void test_negate_eval()
    { PMeanPTemplate of( 2, NULL ); test_negate_flag( EVAL, &of ); }
  void test_negate_grad()
    { PMeanPTemplate of( 2, NULL ); test_negate_flag( GRAD, &of ); }
  void test_negate_diag()
    { PMeanPTemplate of( 2, NULL ); test_negate_flag( DIAG, &of ); }
  void test_negate_hess()
    { PMeanPTemplate of( 2, NULL ); test_negate_flag( HESS, &of ); }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PMeanPTemplateTest, "PMeanPTemplateTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PMeanPTemplateTest, "Unit");

void PMeanPTemplateTest::setUp()
{
  MsqPrintError err(std::cout);
  
  // Create a triangle mesh with three free vertices
  const double coords[] = { 1, 1, 0,
                            2, 1, 0,
                            3, 1, 0,
                            2, 2, 0,
                            1, 3, 0,
                            1, 2, 0 };
  const bool fixed_vtx[] = { true,
                             false,
                             true,
                             false,
                             true,
                             false };
  const size_t tri_conn[] = { 0, 1, 5,
                              1, 2, 3,
                              3, 4, 5,
                              1, 3, 5 };
  mPatch.fill( 6, coords, 
               4, TRIANGLE, tri_conn,
               fixed_vtx, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));

}

/** Define a fake quality metric for testing the objective function
 *
 * Returns, for each vertex in the patch, the distance of that
 * vertex from the origin as the quality.  Each evaluation depends
 * only on a single vertex.
 */
class DistTestMetric : public VertexQM
{
public:
  DistTestMetric() : falseEval(false), failEval(false) {}
  string get_name() const { return "Fake metric for testing objective function"; }
  int get_negate_flag() const { return 1; }
  bool evaluate( PatchData& pd, size_t vtx_idx, double &value, MsqError& err );
  bool evaluate_with_indices( PatchData& pd, size_t vtx_idx, double &value, vector<size_t>& indices, MsqError& err );
  //bool evaluate_with_gradient( PatchData& pd, size_t vtx_idx, double &value, vector<size_t>& indices, vector<Vector3D>& grad, MsqError& err );
  bool falseEval;
  bool failEval;
};

bool DistTestMetric::evaluate( PatchData& pd, size_t vtx_idx, 
                               double &value, MsqError& err )
{
  if (failEval) {
    MSQ_SETERR(err)(MsqError::INVALID_STATE);
    return true;
  }

  const MsqVertex& vtx = pd.vertex_by_index( vtx_idx );
  value = vtx.length_squared();
  return !falseEval;
}

bool DistTestMetric::evaluate_with_indices( PatchData& pd, size_t vtx_idx, 
                               double &value, vector<size_t>& indices, 
                               MsqError& err )
{
  if (failEval) {
    MSQ_SETERR(err)(MsqError::INVALID_STATE);
    return true;
  }

  indices.clear();
  if (vtx_idx < pd.num_free_vertices())
    indices.push_back( vtx_idx );
  
  const MsqVertex& vtx = pd.vertex_by_index( vtx_idx );
  value = vtx.length_squared();
  return !falseEval;
}


void PMeanPTemplateTest::test_eval_type( OFTestMode eval_func, ObjectiveFunction::EvalType type )
{
  PMeanPTemplate func( 1, NULL );
  ObjectiveFunctionTests::test_eval_type( type, eval_func, &func );
}


void PMeanPTemplateTest::test_evaluate( double power )
{
  MsqPrintError err(cout);
  double value;
  bool rval;
  
  DistTestMetric metric;
  PMeanPTemplate func( power, &metric );
  rval = func.evaluate( ObjectiveFunction::CALCULATE, mPatch, value, OF_FREE_EVALS_ONLY, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  
  check_result( mPatch, power, value );
}

void PMeanPTemplateTest::test_gradient( double power )
{
  MsqPrintError err(cout);
  double value;
  bool rval;
  vector<Vector3D> grad;
  
  DistTestMetric metric;
  PMeanPTemplate func( power, &metric );
  rval = func.evaluate_with_gradient( ObjectiveFunction::CALCULATE, mPatch, value, grad, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_EQUAL(mPatch.num_free_vertices(), grad.size());
  
  if (!grad.empty())
    check_result( mPatch, power, value, arrptr(grad) );
}

void PMeanPTemplateTest::test_diagonal( double power )
{
  MsqPrintError err(cout);
  double value;
  bool rval;
  vector<Vector3D> grad;
  vector<SymMatrix3D> Hess;
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  
  DistTestMetric metric;
  PMeanPTemplate func( power, &metric );
  rval = func.evaluate_with_Hessian_diagonal( ObjectiveFunction::CALCULATE, mPatch, value, grad, Hess, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  size_t n = mPatch.num_free_vertices();
  CPPUNIT_ASSERT_EQUAL( n, grad.size() );
  CPPUNIT_ASSERT_EQUAL( n, Hess.size() );
  
  vector<Matrix3D> Hessians(n);
  for (size_t r = 0; r < n; ++r) 
    Hessians[r] = Hess[r];
 
  if (!grad.empty())
    check_result( mPatch, power, value, arrptr(grad), arrptr(Hessians) );
}


void PMeanPTemplateTest::test_Hessian( double power )
{
  MsqPrintError err(cout);
  double value;
  bool rval;
  vector<Vector3D> grad;
  MsqHessian Hess;
  Hess.initialize( mPatch, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  
  DistTestMetric metric;
  PMeanPTemplate func( power, &metric );
  rval = func.evaluate_with_Hessian( ObjectiveFunction::CALCULATE, mPatch, value, grad, Hess, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  size_t n = mPatch.num_free_vertices();
  CPPUNIT_ASSERT_EQUAL( n, grad.size() );
  CPPUNIT_ASSERT_EQUAL( n, Hess.size() );
  
  Matrix3D zero( 0, 0, 0, 0, 0, 0, 0, 0, 0 );
  vector<Matrix3D> Hessians(n);
  for (size_t r = 0; r < n; ++r) {
    Matrix3D* mat = Hess.get_block( r, r );
    CPPUNIT_ASSERT( mat != 0 );
    Hessians[r] = *mat;
    
    for (size_t c = r+1; c < n; ++c) {
      mat = Hess.get_block( r, c );
      if (mat)
        CPPUNIT_ASSERT_MATRICES_EQUAL( zero, *mat, EPSILON );
    }
  }
 
  if (!grad.empty())
    check_result( mPatch, power, value, arrptr(grad), arrptr(Hessians) );
}

void PMeanPTemplateTest::check_result( PatchData& pd, double power, double value, 
                                   Vector3D* gradient, Matrix3D* Hessian )
{
  MsqPrintError err(cout);
  double mvalue, sum = 0;
  bool rval;
  vector<Vector3D> grads;
  vector<Matrix3D> Hess;
  vector<size_t> indices;
  
  
  DistTestMetric metric;
  vector<size_t> handles;
  metric.get_evaluations( pd, handles, OF_FREE_EVALS_ONLY, err );
  ASSERT_NO_ERROR(err);
  
  for (size_t i = 0; i < handles.size(); ++i)
  {
    rval = metric.evaluate_with_Hessian( pd, handles[i], mvalue, indices, grads, Hess, err );
    CPPUNIT_ASSERT( !MSQ_CHKERR(err) && rval );
    sum += pow( mvalue, power );
    
    //if (!OF_FREE_EVALS_ONLY && indices.empty())
    //  continue;
      
    CPPUNIT_ASSERT_EQUAL( (size_t)1, indices.size() );
    CPPUNIT_ASSERT_EQUAL( handles[i], indices[0] );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, grads.size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, Hess.size() );
    
    if (gradient)
    {
      double f = power * pow( mvalue, power - 1 ) / handles.size();
      CPPUNIT_ASSERT_VECTORS_EQUAL( f * grads[0], gradient[i], EPSILON );
    }
    if (Hessian)
    {
      double f = power / handles.size();
      double p2 = (power - 1) * pow( mvalue, power - 2 );
      double p1 = pow( mvalue, power - 1 );
      Matrix3D m;
      m.outer_product( grads[0], grads[0] );
      m *= p2;
      m += p1 * Hess[0];
      m *= f;
      CPPUNIT_ASSERT_MATRICES_EQUAL( m, Hessian[i], EPSILON );
    }  
  }
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( sum / handles.size(), value, EPSILON );
}

