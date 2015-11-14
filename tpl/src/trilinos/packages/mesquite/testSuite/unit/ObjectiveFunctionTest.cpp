/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 13-Nov-02 at 18:05:56
//  LAST-MOD: 23-Jul-03 at 17:40:28 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file ObjectiveFunctionTest.cpp

Unit testing of various functions in the ObjectiveFunction class. 
*/
// DESCRIP-END.
//



#include "Mesquite.hpp"
#include "Mesquite_LPtoPTemplate.hpp"
#include "Mesquite_MaxTemplate.hpp"
#include "Mesquite_LInfTemplate.hpp"
#include "Mesquite_IdealWeightInverseMeanRatio.hpp"
#include "Mesquite_MsqHessian.hpp"

#include "PatchDataInstances.hpp"
#include "ObjectiveFunctionTests.hpp"
#include "Mesquite_cppunit/extensions/HelperMacros.h"
#include "UnitUtil.hpp"
#include <iterator>
#include <algorithm>

using namespace Mesquite;
using std::cout;
using std::endl;
using std::cerr;

class ObjectiveFunctionTest : public CppUnit::TestFixture, public ObjectiveFunctionTests
{
private:
  CPPUNIT_TEST_SUITE(ObjectiveFunctionTest);
  
  CPPUNIT_TEST (test_eval_OF_value_LPtoP_L1);
  CPPUNIT_TEST (test_grad_OF_value_LPtoP_L1);
  CPPUNIT_TEST (test_diag_OF_value_LPtoP_L1);
  CPPUNIT_TEST (test_hess_OF_value_LPtoP_L1);
  CPPUNIT_TEST (test_eval_OF_value_LPtoP_L2);
  CPPUNIT_TEST (test_grad_OF_value_LPtoP_L2);
  CPPUNIT_TEST (test_diag_OF_value_LPtoP_L2);
  CPPUNIT_TEST (test_hess_OF_value_LPtoP_L2);
  CPPUNIT_TEST (test_eval_OF_value_LPtoP_scaled);
  CPPUNIT_TEST (test_grad_OF_value_LPtoP_scaled);
  CPPUNIT_TEST (test_diag_OF_value_LPtoP_scaled);
  CPPUNIT_TEST (test_hess_OF_value_LPtoP_scaled);
  CPPUNIT_TEST (test_eval_OF_value_LInf);
  CPPUNIT_TEST (test_eval_OF_value_max);
  
  CPPUNIT_TEST (test_compute_gradient_LPtoPTemplate_L1);
  CPPUNIT_TEST (test_compute_gradient_LPtoPTemplate_L2);
  CPPUNIT_TEST (test_compute_gradient_LPtoPTemplate_L2_scaled);
  CPPUNIT_TEST (test_hessian_gradient_LPtoPTemplate_L1);
  CPPUNIT_TEST (test_hessian_gradient_LPtoPTemplate_L2);
  CPPUNIT_TEST (test_hessian_gradient_LPtoPTemplate_L2_scaled);
  CPPUNIT_TEST (test_diagonal_gradient_LPtoPTemplate_L1);
  CPPUNIT_TEST (test_diagonal_gradient_LPtoPTemplate_L2);
  CPPUNIT_TEST (test_diagonal_gradient_LPtoPTemplate_L2_scaled);
  CPPUNIT_TEST (test_hessian_diagonal_LPtoPTemplate_L1);
  CPPUNIT_TEST (test_hessian_diagonal_LPtoPTemplate_L2);
  CPPUNIT_TEST (test_hessian_diagonal_LPtoPTemplate_L2_scaled);
  CPPUNIT_TEST (test_compute_ana_hessian_tet);
  CPPUNIT_TEST (test_compute_ana_hessian_tet_scaled);
  
  CPPUNIT_TEST (test_LPtoP_invalid_qm_eval);
  CPPUNIT_TEST (test_LPtoP_invalid_qm_grad);
  CPPUNIT_TEST (test_LPtoP_invalid_qm_diag);
  CPPUNIT_TEST (test_LPtoP_invalid_qm_hess);
  CPPUNIT_TEST (test_LInf_invalid_qm_eval);
  CPPUNIT_TEST (test_max_invalid_qm_eval);
  
  CPPUNIT_TEST (test_LPtoP_qm_error_eval);
  CPPUNIT_TEST (test_LPtoP_qm_error_grad);
  CPPUNIT_TEST (test_LPtoP_qm_error_diag);
  CPPUNIT_TEST (test_LPtoP_qm_error_hess);
  CPPUNIT_TEST (test_LInf_qm_error_eval);
  CPPUNIT_TEST (test_max_qm_error_eval);
  
  CPPUNIT_TEST (test_LPtoP_eval_calc);
  CPPUNIT_TEST (test_LPtoP_eval_accum);
  CPPUNIT_TEST (test_LPtoP_eval_save);
  CPPUNIT_TEST (test_LPtoP_eval_update);
  CPPUNIT_TEST (test_LPtoP_eval_temp);
  
  CPPUNIT_TEST (test_LPtoP_grad_calc);
  CPPUNIT_TEST (test_LPtoP_grad_save);
  CPPUNIT_TEST (test_LPtoP_grad_update);
  CPPUNIT_TEST (test_LPtoP_grad_temp);
  
  CPPUNIT_TEST (test_LPtoP_diag_calc);
  CPPUNIT_TEST (test_LPtoP_diag_save);
  CPPUNIT_TEST (test_LPtoP_diag_update);
  CPPUNIT_TEST (test_LPtoP_diag_temp);
  
  CPPUNIT_TEST (test_LPtoP_hess_calc);
  CPPUNIT_TEST (test_LPtoP_hess_save);
  CPPUNIT_TEST (test_LPtoP_hess_update);
  CPPUNIT_TEST (test_LPtoP_hess_temp);
  
  CPPUNIT_TEST (test_LPtoP_clone_L1);
  CPPUNIT_TEST (test_LPtoP_clone_L2);
  CPPUNIT_TEST (test_LInf_clone);
  CPPUNIT_TEST (test_max_clone);

  CPPUNIT_TEST (test_LPtoP_negate_flag_eval);
  CPPUNIT_TEST (test_LPtoP_negate_flag_grad);
  CPPUNIT_TEST (test_LPtoP_negate_flag_diag);
  CPPUNIT_TEST (test_LPtoP_negate_flag_hess);
  CPPUNIT_TEST (test_LInf_negate_flag_eval);
  CPPUNIT_TEST (test_max_negate_flag_eval);

  CPPUNIT_TEST_SUITE_END();
  
  void test_LPtoP_value( short P, bool scale,
                         const std::vector<double>& values, OFTestMode mode );
  
  void test_LInf_value( const std::vector<double>& values );
  
  void test_max_value( const std::vector<double>& values );
  
  void test_max_negate_flag( ObjectiveFunctionTemplate& OF );
  
  std::vector<double> values1, values2;
   
public:
  ObjectiveFunctionTest()
  {}
  
  void setUp() 
  {
    const double v1[] = { 0, 1, 2, 3, 4, 5 };
    const double v2[] = { -2.5, 0.5, 3.0, M_PI, -25.0, 1.0 };
    std::copy( v1, v1+sizeof(v1)/sizeof(v1[0]), std::back_inserter(values1) );
    std::copy( v2, v2+sizeof(v2)/sizeof(v2[0]), std::back_inserter(values2) );
  }
  
  void test_eval_OF_value_LPtoP_L1() 
  {
    test_LPtoP_value( 1, false, values1, EVAL );
    test_LPtoP_value( 1, false, values2, EVAL );
  }
  
  void test_grad_OF_value_LPtoP_L1()
  {
    test_LPtoP_value( 1, false, values1, GRAD );
    test_LPtoP_value( 1, false, values2, GRAD );
  }
  
  void test_diag_OF_value_LPtoP_L1()
  {
    test_LPtoP_value( 1, false, values1, DIAG );
    test_LPtoP_value( 1, false, values2, DIAG );
  }
  
  void test_hess_OF_value_LPtoP_L1()
  {
    test_LPtoP_value( 1, false, values1, HESS );
    test_LPtoP_value( 1, false, values2, HESS );
  }
  
  void test_eval_OF_value_LPtoP_L2()
  {
    test_LPtoP_value( 2, false, values1, EVAL );
    test_LPtoP_value( 2, false, values2, EVAL );
  }
  
  void test_grad_OF_value_LPtoP_L2()
  {
    test_LPtoP_value( 2, false, values1, GRAD );
    test_LPtoP_value( 2, false, values2, GRAD );
  }
  
  void test_diag_OF_value_LPtoP_L2()
  {
    test_LPtoP_value( 2, false, values1, DIAG );
    test_LPtoP_value( 2, false, values2, DIAG );
  }
  
  void test_hess_OF_value_LPtoP_L2()
  {
    test_LPtoP_value( 2, false, values1, HESS );
    test_LPtoP_value( 2, false, values2, HESS );
  }
  
  void test_eval_OF_value_LPtoP_scaled()
  {
    test_LPtoP_value( 2, true, values1, EVAL );
    test_LPtoP_value( 2, true, values2, EVAL );
  }
  
  void test_grad_OF_value_LPtoP_scaled()
  {
    test_LPtoP_value( 2, true, values1, GRAD );
    test_LPtoP_value( 2, true, values2, GRAD );
  }
  
  void test_diag_OF_value_LPtoP_scaled()
  {
    test_LPtoP_value( 2, true, values1, DIAG );
    test_LPtoP_value( 2, true, values2, DIAG );
  }
  
  void test_hess_OF_value_LPtoP_scaled()
  {
    test_LPtoP_value( 2, true, values1, EVAL );
    test_LPtoP_value( 2, true, values2, EVAL );
  }
  
  void test_eval_OF_value_LInf()
  {
    test_LInf_value( values1 );
    test_LInf_value( values2 );
  }
  
  void test_eval_OF_value_max()
  {
    test_max_value( values1 );
    test_max_value( values2 );
  }

  void test_compute_gradient_LPtoPTemplate_L1()
    { LPtoPTemplate LP1( 1, NULL ); compare_numerical_gradient( &LP1 ); }

  void test_compute_gradient_LPtoPTemplate_L2()
    { LPtoPTemplate LP1( 2, NULL ); compare_numerical_gradient( &LP1 ); }
  
  void test_compute_gradient_LPtoPTemplate_L2_scaled()
    { 
      LPtoPTemplate LP1( 2, NULL ); 
      LP1.set_dividing_by_n(true);
      compare_numerical_gradient( &LP1 ); 
    }

  void test_hessian_gradient_LPtoPTemplate_L1()
    { LPtoPTemplate LP1( 1, NULL ); compare_hessian_gradient( &LP1 ); }

  void test_hessian_gradient_LPtoPTemplate_L2()
    { LPtoPTemplate LP1( 2, NULL ); compare_hessian_gradient( &LP1 ); }
  
  void test_hessian_gradient_LPtoPTemplate_L2_scaled()
    { 
      LPtoPTemplate LP1( 2, NULL ); 
      LP1.set_dividing_by_n(true);
      compare_hessian_gradient( &LP1 ); 
    }

  void test_diagonal_gradient_LPtoPTemplate_L1()
    { LPtoPTemplate LP1( 1, NULL ); compare_diagonal_gradient( &LP1 ); }

  void test_diagonal_gradient_LPtoPTemplate_L2()
    { LPtoPTemplate LP1( 2, NULL ); compare_diagonal_gradient( &LP1 ); }
  
  void test_diagonal_gradient_LPtoPTemplate_L2_scaled()
    { 
      LPtoPTemplate LP1( 2, NULL ); 
      LP1.set_dividing_by_n(true);
      compare_diagonal_gradient( &LP1 ); 
    }

  void test_hessian_diagonal_LPtoPTemplate_L1()
    { LPtoPTemplate LP1( 1, NULL ); compare_hessian_diagonal( &LP1 ); }

  void test_hessian_diagonal_LPtoPTemplate_L2()
    { LPtoPTemplate LP1( 2, NULL ); compare_hessian_diagonal( &LP1 ); }
  
  void test_hessian_diagonal_LPtoPTemplate_L2_scaled()
    { 
      LPtoPTemplate LP1( 2, NULL ); 
      LP1.set_dividing_by_n(true);
      compare_hessian_diagonal( &LP1 ); 
    }

  void test_compute_ana_hessian_tet();
  
  void test_compute_ana_hessian_tet_scaled();

  void test_LPtoP_invalid_qm_eval()
    { LPtoPTemplate LP1( 1, NULL ); test_handles_invalid_qm( EVAL, &LP1 ); }
  void test_LPtoP_invalid_qm_grad()
    { LPtoPTemplate LP1( 1, NULL ); test_handles_invalid_qm( GRAD, &LP1 ); }
  void test_LPtoP_invalid_qm_diag()
    { LPtoPTemplate LP1( 1, NULL ); test_handles_invalid_qm( DIAG, &LP1 ); }
  void test_LPtoP_invalid_qm_hess()
    { LPtoPTemplate LP1( 1, NULL ); test_handles_invalid_qm( HESS, &LP1 ); }
  void test_LInf_invalid_qm_eval()
    { LInfTemplate OF( NULL ); test_handles_invalid_qm( EVAL, &OF ); }
  void test_max_invalid_qm_eval()
    { MaxTemplate OF( NULL ); test_handles_invalid_qm( EVAL, &OF ); }
  
  void test_LPtoP_qm_error_eval()
    { LPtoPTemplate LP1( 1, NULL ); test_handles_qm_error( EVAL, &LP1 ); }
  void test_LPtoP_qm_error_grad()
    { LPtoPTemplate LP1( 1, NULL ); test_handles_qm_error( GRAD, &LP1 ); }
  void test_LPtoP_qm_error_diag()
    { LPtoPTemplate LP1( 1, NULL ); test_handles_qm_error( DIAG, &LP1 ); }
  void test_LPtoP_qm_error_hess()
    { LPtoPTemplate LP1( 1, NULL ); test_handles_qm_error( HESS, &LP1 ); }
  void test_LInf_qm_error_eval()
    { LInfTemplate OF( NULL ); test_handles_qm_error( EVAL, &OF ); }
  void test_max_qm_error_eval()
    { MaxTemplate OF( NULL ); test_handles_qm_error( EVAL, &OF ); }
  
  void test_LPtoP_eval_calc()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::CALCULATE, EVAL, &OF );
  }
  void test_LPtoP_eval_accum()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::ACCUMULATE, EVAL, &OF );
  }
  void test_LPtoP_eval_save()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::SAVE, EVAL, &OF );
  }
  void test_LPtoP_eval_update()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::UPDATE, EVAL, &OF );
  }
  void test_LPtoP_eval_temp()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::TEMPORARY, EVAL, &OF );
  }
  
  void test_LPtoP_grad_calc()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::CALCULATE, GRAD, &OF );
  }
  void test_LPtoP_grad_save()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::SAVE, GRAD, &OF );
  }
  void test_LPtoP_grad_update()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::UPDATE, GRAD, &OF );
  }
  void test_LPtoP_grad_temp()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::TEMPORARY, GRAD, &OF );
  }
  
  void test_LPtoP_diag_calc()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::CALCULATE, DIAG, &OF );
  }
  void test_LPtoP_diag_save()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::SAVE, DIAG, &OF );
  }
  void test_LPtoP_diag_update()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::UPDATE, DIAG, &OF );
  }
  void test_LPtoP_diag_temp()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::TEMPORARY, DIAG, &OF );
  }
  
  void test_LPtoP_hess_calc()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::CALCULATE, HESS, &OF );
  }
  void test_LPtoP_hess_save()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::SAVE, HESS, &OF );
  }
  void test_LPtoP_hess_update()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::UPDATE, HESS, &OF );
  }
  void test_LPtoP_hess_temp()
  { 
    LPtoPTemplate OF( 1, NULL ); 
    test_eval_type( ObjectiveFunction::TEMPORARY, HESS, &OF );
  }
  
  void test_LPtoP_clone_L1()
    { LPtoPTemplate OF( 1, NULL ); test_clone(&OF); }
  void test_LPtoP_clone_L2()
    { LPtoPTemplate OF( 2, NULL ); test_clone(&OF); }
  void test_LInf_clone()
    { LInfTemplate OF( NULL ); test_clone(&OF); }
  void test_max_clone()
    { MaxTemplate OF( NULL ); test_clone(&OF); }

    //This function tests to make sure that LPtoP handles the negate
    // flag correctly for the analytical Hessian, gradient, and evaluate.
    // It does this by creating two objective functions that are
    // identical except that the metric in one has negate flag of -1
    // and the metric in the other has a negate flag of 1.  Thus,
    // the Hessians, gradients, and function values should be the
    // the same except for a negative sign.
  void test_LPtoP_negate_flag_eval()
   { LPtoPTemplate LP( 2, NULL ); test_negate_flag( EVAL, &LP ); }
  void test_LPtoP_negate_flag_grad()
   { LPtoPTemplate LP( 2, NULL ); test_negate_flag( GRAD, &LP ); }
  void test_LPtoP_negate_flag_diag()
   { LPtoPTemplate LP( 2, NULL ); test_negate_flag( DIAG, &LP ); }
  void test_LPtoP_negate_flag_hess()
   { LPtoPTemplate LP( 2, NULL ); test_negate_flag( HESS, &LP ); }
  void test_LInf_negate_flag_eval()
   { LInfTemplate OF( NULL ); test_max_negate_flag( OF ); }
  void test_max_negate_flag_eval()
   { MaxTemplate OF( NULL ); test_max_negate_flag( OF ); }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ObjectiveFunctionTest, "ObjectiveFunctionTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ObjectiveFunctionTest, "Unit");

void ObjectiveFunctionTest::test_compute_ana_hessian_tet()
  {
    MsqPrintError err(cout);
    PatchData tetPatch;
    create_qm_two_tet_patch(tetPatch,err);
    ASSERT_NO_ERROR(err);
    
    
    // creates a mean ratio quality metric ...
    IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio(err);CPPUNIT_ASSERT(!err);
    mean_ratio->set_averaging_method(QualityMetric::SUM);
    
    // ... and builds an objective function with it
    LPtoPTemplate LP2(mean_ratio, 2, err);
    
    MsqHessian H;
    std::vector<Vector3D> g;
    double dummy;
    H.initialize(tetPatch, err); CPPUNIT_ASSERT(!err);
    LP2.evaluate_with_Hessian(ObjectiveFunction::CALCULATE, tetPatch, dummy, g, H, err); CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL( tetPatch.num_free_vertices(), g.size() );

    Matrix3D mat00(" 2.44444  0.2566   0.181444 "
		   " 0.2566   2.14815  0.104757 "
		   " 0.181444 0.104757 2.07407 ");


    Matrix3D mat13(" 5.47514 3.16659    9.83479 "
		   " -1.11704 -5.29718 -3.67406 "
		   " 10.3635 -13.5358  -15.5638 ");

    CPPUNIT_ASSERT_MATRICES_EQUAL( mat00, *H.get_block(0,0), 1e-4 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( mat13, *H.get_block(1,3), 1e-4 );

    delete mean_ratio;
  }
  
void ObjectiveFunctionTest::test_compute_ana_hessian_tet_scaled()
  {
    MsqPrintError err(cout);
    PatchData tetPatch;
    create_qm_two_tet_patch(tetPatch,err);
    ASSERT_NO_ERROR(err);
    
    // creates a mean ratio quality metric ...
    IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio(err);CPPUNIT_ASSERT(!err);
    mean_ratio->set_averaging_method(QualityMetric::SUM);
    
    // ... and builds an objective function with it
    LPtoPTemplate LP2(mean_ratio, 2, err);
    LP2.set_dividing_by_n(true);
    
    MsqHessian H;
    std::vector<Vector3D> g;
    double dummy;
    H.initialize(tetPatch, err); CPPUNIT_ASSERT(!err);
    LP2.evaluate_with_Hessian(ObjectiveFunction::CALCULATE, tetPatch, dummy, g, H, err); CPPUNIT_ASSERT(!err);
    CPPUNIT_ASSERT_EQUAL( tetPatch.num_free_vertices(), g.size() );

    Matrix3D mat00(" 2.44444  0.2566   0.181444 "
		   " 0.2566   2.14815  0.104757 "
		   " 0.181444 0.104757 2.07407 ");

    mat00*=.5;
    
    Matrix3D mat13(" 5.47514 3.16659    9.83479 "
		   " -1.11704 -5.29718 -3.67406 "
		   " 10.3635 -13.5358  -15.5638 ");

    mat13*=.5;
    
    CPPUNIT_ASSERT_MATRICES_EQUAL( mat00, *H.get_block(0,0), 1e-4 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( mat13, *H.get_block(1,3), 1e-4 );
    
//    cout << H <<endl;
    delete mean_ratio;
  }
  
void ObjectiveFunctionTest::test_LPtoP_value( short P, bool scale,
                         const std::vector<double>& values, OFTestMode mode )
{
  CPPUNIT_ASSERT(!values.empty());

  LPtoPTemplate OF( P, NULL );
  OF.set_dividing_by_n(scale);
  
  double expected = 0.0;
  for (std::vector<double>::const_iterator i = values.begin(); i != values.end(); ++i)
    expected += std::pow( fabs(*i), P );
  if (scale) 
    expected /= values.size();
    
  test_value( arrptr(values), values.size(), expected, mode, &OF );
}
  
void ObjectiveFunctionTest::test_LInf_value( const std::vector<double>& values )
{
  CPPUNIT_ASSERT(!values.empty());

  LInfTemplate OF( NULL );
  
  double expected = fabs(values[0]);
  for (unsigned i = 1; i < values.size(); ++i)
    if (fabs(values[i]) > expected)
      expected = fabs(values[i]);
    
  test_value( arrptr(values), values.size(), expected, EVAL, &OF );
}
  
void ObjectiveFunctionTest::test_max_value( const std::vector<double>& values )
{
  CPPUNIT_ASSERT(!values.empty());

  MaxTemplate OF( NULL );
  
  double expected = fabs(values[0]);
  for (unsigned i = 1; i < values.size(); ++i)
    if (values[i] > expected)
      expected = values[i];
    
  test_value( arrptr(values), values.size(), expected, EVAL, &OF );
}

void ObjectiveFunctionTest::test_max_negate_flag( ObjectiveFunctionTemplate& of )
{
    // use only positive values so this test works for both
    // MaxTemplate and LInfTemplate
  const double some_vals[] = { 1, 2, 3, 4, 5, 6, 0.5 };
  const unsigned num_vals = sizeof(some_vals)/sizeof(some_vals[0]);
  OFTestQM metric( some_vals, num_vals );
  metric.set_negate_flag(-1);
  of.set_quality_metric(&metric);

    // get OF value
  MsqPrintError err(cout);
  bool rval;
  double value;
  rval = of.evaluate( ObjectiveFunction::CALCULATE, patch(), value, false, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT(rval);
  
    // find min value
  double expected = *std::min_element( some_vals, some_vals+num_vals );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( -expected, value, 1e-6 );
}
