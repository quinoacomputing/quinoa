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


/** \file StdDevTemplateTest.cpp
 *  \brief Unit tests for StdDevTemplate and VarianceTemplate
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Mesquite_StdDevTemplate.hpp"
#include "Mesquite_VarianceTemplate.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_PatchData.hpp"
#include "ObjectiveFunctionTests.hpp"
#include "Mesquite_MsqHessian.hpp"

using namespace Mesquite;
using namespace std;

const double EPSILON = 1e-4;

class StdDevTemplateTest : public CppUnit::TestFixture, public ObjectiveFunctionTests
{
private:
  CPPUNIT_TEST_SUITE( StdDevTemplateTest );

  CPPUNIT_TEST( test_eval_calc );
  CPPUNIT_TEST( test_eval_accum );
  CPPUNIT_TEST( test_eval_save );
  CPPUNIT_TEST( test_eval_update );
  CPPUNIT_TEST( test_eval_temp );

  CPPUNIT_TEST( test_eval_calc_sqr );
  CPPUNIT_TEST( test_eval_accum_sqr );
  CPPUNIT_TEST( test_eval_save_sqr );
  CPPUNIT_TEST( test_eval_update_sqr );
  CPPUNIT_TEST( test_eval_temp_sqr );

  CPPUNIT_TEST( test_grad_calc );
  CPPUNIT_TEST( test_grad_save );
  CPPUNIT_TEST( test_grad_update );
  CPPUNIT_TEST( test_grad_temp );

  CPPUNIT_TEST( test_grad_calc_sqr );
  CPPUNIT_TEST( test_grad_save_sqr );
  CPPUNIT_TEST( test_grad_update_sqr );
  CPPUNIT_TEST( test_grad_temp_sqr );

  CPPUNIT_TEST( test_diag_calc );
  CPPUNIT_TEST( test_diag_save );
  CPPUNIT_TEST( test_diag_update );
  CPPUNIT_TEST( test_diag_temp );

  CPPUNIT_TEST( test_diag_calc_sqr );
  CPPUNIT_TEST( test_diag_save_sqr );
  CPPUNIT_TEST( test_diag_update_sqr );
  CPPUNIT_TEST( test_diag_temp_sqr );

  CPPUNIT_TEST( test_hessian_fails );
  CPPUNIT_TEST( test_hessian_fails_sqr );
  
  CPPUNIT_TEST( test_failed_metric_in_eval );
  CPPUNIT_TEST( test_failed_metric_in_grad );
  CPPUNIT_TEST( test_failed_metric_in_eval_sqr  );
  CPPUNIT_TEST( test_failed_metric_in_grad_sqr  );
  
  CPPUNIT_TEST( test_false_metric_in_eval );
  CPPUNIT_TEST( test_false_metric_in_grad );
  CPPUNIT_TEST( test_false_metric_in_eval_sqr  );
  CPPUNIT_TEST( test_false_metric_in_grad_sqr  );
  
  CPPUNIT_TEST( test_evaluate );
  CPPUNIT_TEST( test_evaluate_sqr  );
  
  CPPUNIT_TEST( test_numerical_gradient );
  CPPUNIT_TEST( test_numerical_gradient_sqr  );
  
  CPPUNIT_TEST( test_diagonal_gradient );
  CPPUNIT_TEST( test_diagonal_gradient_sqr );
  
  CPPUNIT_TEST( test_hessian_diag );
  CPPUNIT_TEST( test_hessian_diag_sqr );
 
  CPPUNIT_TEST( test_clone );
  CPPUNIT_TEST( test_clone_sqr );
  
  CPPUNIT_TEST( test_eval_negate );
  CPPUNIT_TEST( test_eval_negate_sqr );
  CPPUNIT_TEST( test_grad_negate );
  CPPUNIT_TEST( test_grad_negate_sqr );
  CPPUNIT_TEST( test_diag_negate );
  CPPUNIT_TEST( test_diag_negate_sqr );

  CPPUNIT_TEST_SUITE_END();

public:

  void test_eval_calc()   
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::CALCULATE,  EVAL, &of ); }
  void test_eval_accum() 
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::ACCUMULATE, EVAL, &of ); }
  void test_eval_save()  
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::SAVE,       EVAL, &of ); }
  void test_eval_update()
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::UPDATE,     EVAL, &of ); }
  void test_eval_temp()  
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::TEMPORARY,  EVAL, &of ); }

  void test_eval_calc_sqr()  
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::CALCULATE,  EVAL, &of ); }
  void test_eval_accum_sqr() 
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::ACCUMULATE, EVAL, &of ); }
  void test_eval_save_sqr()  
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::SAVE,       EVAL, &of ); }
  void test_eval_update_sqr()
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::UPDATE,     EVAL, &of ); }
  void test_eval_temp_sqr()  
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::TEMPORARY,  EVAL, &of ); }

  void test_grad_calc()   
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::CALCULATE,  GRAD, &of ); }
  void test_grad_accum() 
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::ACCUMULATE, GRAD, &of ); }
  void test_grad_save()  
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::SAVE,       GRAD, &of ); }
  void test_grad_update()
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::UPDATE,     GRAD, &of ); }
  void test_grad_temp()  
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::TEMPORARY,  GRAD, &of ); }

  void test_grad_calc_sqr()  
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::CALCULATE,  GRAD, &of ); }
  void test_grad_accum_sqr() 
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::ACCUMULATE, GRAD, &of ); }
  void test_grad_save_sqr()  
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::SAVE,       GRAD, &of ); }
  void test_grad_update_sqr()
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::UPDATE,     GRAD, &of ); }
  void test_grad_temp_sqr()  
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::TEMPORARY,  GRAD, &of ); }

  void test_diag_calc()   
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::CALCULATE,  DIAG, &of ); }
  void test_diag_accum() 
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::ACCUMULATE, DIAG, &of ); }
  void test_diag_save()  
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::SAVE,       DIAG, &of ); }
  void test_diag_update()
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::UPDATE,     DIAG, &of ); }
  void test_diag_temp()  
    { StdDevTemplate of(NULL); test_eval_type( ObjectiveFunction::TEMPORARY,  DIAG, &of ); }

  void test_diag_calc_sqr()  
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::CALCULATE,  DIAG, &of ); }
  void test_diag_accum_sqr() 
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::ACCUMULATE, DIAG, &of ); }
  void test_diag_save_sqr()  
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::SAVE,       DIAG, &of ); }
  void test_diag_update_sqr()
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::UPDATE,     DIAG, &of ); }
  void test_diag_temp_sqr()  
    { VarianceTemplate of(NULL); test_eval_type( ObjectiveFunction::TEMPORARY,  DIAG, &of ); }

  void test_evaluate( ) ;
  void test_evaluate_sqr( );
  
  void test_numerical_gradient( )      
    { StdDevTemplate of(NULL); compare_numerical_gradient(&of); }
  void test_numerical_gradient_sqr( )  
    { VarianceTemplate of(NULL); compare_numerical_gradient(&of); }
  
  void test_diagonal_gradient( )      
    { StdDevTemplate of(NULL); compare_diagonal_gradient(&of); }
  void test_diagonal_gradient_sqr( )  
    { VarianceTemplate of(NULL); compare_diagonal_gradient(&of); }

  void test_hessian_diag()
    { StdDevTemplate of(NULL); compare_numerical_hessian_diagonal(&of); }
  void test_hessian_diag_sqr()
    { VarianceTemplate of(NULL); compare_numerical_hessian_diagonal(&of); }

  void test_hessian_fails();
  void test_hessian_fails_sqr();
  
  void test_failed_metric_in_eval()   
    { StdDevTemplate of(NULL); test_handles_qm_error(EVAL, &of); }
  void test_failed_metric_in_grad()     
    { StdDevTemplate of(NULL); test_handles_qm_error(GRAD, &of); }
  void test_failed_metric_in_eval_sqr() 
    { VarianceTemplate of(NULL); test_handles_qm_error(EVAL, &of); }
  void test_failed_metric_in_grad_sqr() 
    { VarianceTemplate of(NULL); test_handles_qm_error(GRAD, &of); }
  
  void test_false_metric_in_eval()     
    { StdDevTemplate of(NULL); test_handles_invalid_qm(EVAL, &of); }
  void test_false_metric_in_grad()     
    { StdDevTemplate of(NULL); test_handles_invalid_qm(GRAD, &of); }
  void test_false_metric_in_eval_sqr() 
    { VarianceTemplate of(NULL); test_handles_invalid_qm(EVAL, &of); }
  void test_false_metric_in_grad_sqr() 
    { VarianceTemplate of(NULL); test_handles_invalid_qm(GRAD, &of); }
    
  void test_clone()
    { StdDevTemplate of(NULL); ObjectiveFunctionTests::test_clone(&of); }
  void test_clone_sqr()
    { VarianceTemplate of(NULL); ObjectiveFunctionTests::test_clone(&of); }
  
  void test_eval_negate()
    { StdDevTemplate of(NULL); test_negate_flag( EVAL, &of ); }
  void test_eval_negate_sqr()
    { VarianceTemplate of(NULL); test_negate_flag( EVAL, &of ); }
  void test_grad_negate()
    { StdDevTemplate of(NULL); test_negate_flag( GRAD, &of ); }
  void test_grad_negate_sqr()
    { VarianceTemplate of(NULL); test_negate_flag( GRAD, &of ); }
  void test_diag_negate()
    { StdDevTemplate of(NULL); test_negate_flag( DIAG, &of ); }
  void test_diag_negate_sqr()
    { VarianceTemplate of(NULL); test_negate_flag( DIAG, &of ); }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(StdDevTemplateTest, "StdDevTemplateTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(StdDevTemplateTest, "Unit");

static double std_dev_sqr( const double* array, unsigned len )
{
  double sum = 0, sqr_sum = 0;
  for (size_t i = 0; i < len; ++i)
  {
    sum += array[i];
    sqr_sum += array[i]*array[i];
  }
  
  return sqr_sum/len - (sum/len)*(sum/len);
}


void StdDevTemplateTest::test_evaluate()
{
  StdDevTemplate OF(NULL);
  
  const double list1[] = { 1.0, 0.0, -5.0, 0.2, 6.0, 3 };
  const unsigned len1 = sizeof(list1)/sizeof(list1[0]);
  test_value( list1, len1, sqrt(std_dev_sqr(list1,len1)), EVAL, &OF );
  
  const double list2[] = { 20, -30, 40, -50, 60, -70, 80, -90 };
  const unsigned len2 = sizeof(list2)/sizeof(list2[0]);
  test_value( list2, len2, sqrt(std_dev_sqr(list2,len2)), EVAL, &OF );
}

void StdDevTemplateTest::test_evaluate_sqr()
{
  VarianceTemplate OF(NULL);
  
  const double list1[] = { 1.0, 0.0, -5.0, 0.2, 6.0, 3 };
  const unsigned len1 = sizeof(list1)/sizeof(list1[0]);
  test_value( list1, len1, std_dev_sqr(list1,len1), EVAL, &OF );
  
  const double list2[] = { 20, -30, 40, -50, 60, -70, 80, -90 };
  const unsigned len2 = sizeof(list2)/sizeof(list2[0]);
  test_value( list2, len2, std_dev_sqr(list2,len2), EVAL, &OF );
}


void StdDevTemplateTest::test_hessian_fails()
{
  MsqError err;
  double value;
  bool rval;
  vector<Vector3D> grad;
  MsqHessian Hess;
  Hess.initialize( patch(), err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  
  OFTestQM metric( &value, 1 );
  StdDevTemplate func( &metric );
  rval = func.evaluate_with_Hessian( ObjectiveFunction::CALCULATE, patch(), value, grad, Hess, err );
  CPPUNIT_ASSERT(err);
}

void StdDevTemplateTest::test_hessian_fails_sqr()
{
  MsqError err;
  double value;
  bool rval;
  vector<Vector3D> grad;
  MsqHessian Hess;
  Hess.initialize( patch(), err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  
  OFTestQM metric( &value, 1 );
  VarianceTemplate func( &metric );
  rval = func.evaluate_with_Hessian( ObjectiveFunction::CALCULATE, patch(), value, grad, Hess, err );
  CPPUNIT_ASSERT(err);
}
