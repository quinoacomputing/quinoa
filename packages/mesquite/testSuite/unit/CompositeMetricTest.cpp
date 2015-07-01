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


/** \file CompositeMetricTest.cpp
 *  \brief unit tests for IdealWeightMeanRatio quality metric
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "IdealWeightInverseMeanRatio.hpp"
#include "LocalSizeQualityMetric.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include "QualityMetricTester.hpp"
#include "AddQualityMetric.hpp"
#include "MultiplyQualityMetric.hpp"
#include "PowerQualityMetric.hpp"
#include "ScalarAddQualityMetric.hpp"
#include "ScalarMultiplyQualityMetric.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include "UnitUtil.hpp"

using namespace Mesquite;

class CompositeMetricTestBase : public CppUnit::TestFixture
{
private:  

  void test_evaluate( bool ideal, EntityTopology type );

protected:
  QualityMetricTester tester;
  QualityMetric* mMetric;
  virtual bool evaluate( PatchData&, size_t, double&, MsqError& ) = 0;
  
public:

  CompositeMetricTestBase( ) 
    : tester(QualityMetricTester::ALL_FE_EXCEPT_SEPTAHEDRON),
      mMetric(0)
    { tester.ideal_pyramid_base_equals_height( true ); }

  void test_supported_types()
    { tester.test_supported_element_types( mMetric ); }
    
  void test_ideal_element_eval()
  {
    test_evaluate( true, TRIANGLE );
    test_evaluate( true, QUADRILATERAL );
    test_evaluate( true, TETRAHEDRON );
    test_evaluate( true, HEXAHEDRON );
    test_evaluate( true, PRISM );
    test_evaluate( true, PYRAMID );
  }
  
  void test_non_ideal_eval()
  {
    test_evaluate( false, TRIANGLE );
    test_evaluate( false, QUADRILATERAL );
    test_evaluate( false, TETRAHEDRON );
    test_evaluate( false, HEXAHEDRON );
    test_evaluate( false, PRISM );
    test_evaluate( false, PYRAMID );
  }
  
  void test_ideal_element_grad()
    { tester.test_ideal_element_zero_gradient( mMetric, false ); }
  
  void test_ideal_element_hess()
    { tester.test_ideal_element_positive_definite_Hessian( mMetric, false ); }
    
  void test_valid_hessian()
    { tester.test_symmetric_Hessian_diagonal_blocks( mMetric ); }
  
  void test_measures_quality()
    { tester.test_measures_quality( mMetric ); }
  
  void test_gradient_reflects_quality()
    { tester.test_gradient_reflects_quality( mMetric ); }
  
  void test_domain_deviation()
  {
    tester.test_domain_deviation_quality( mMetric );
    tester.test_domain_deviation_gradient( mMetric );
  }
  
  void test_inverted_elements()
    { tester.test_evaluate_inverted_element( mMetric, false ); }
    
  void test_degenerate_elements()
    { tester.test_evaluate_degenerate_element( mMetric, false ); }
    
  void test_get_evaluations()
    { tester.test_get_element_evaluations( mMetric ); }
    
  void test_get_element_indices()
    { tester.test_get_element_indices( mMetric ); }
  
  void test_get_fixed_indices()
    { tester.test_get_indices_fixed( mMetric ); }
  
  void test_eval_with_indices()
    { tester.compare_eval_and_eval_with_indices( mMetric ); }
    
  void test_eval_with_gradient()
  {
    tester.compare_eval_with_indices_and_eval_with_gradient( mMetric );
    tester.compare_analytical_and_numerical_gradients( mMetric );
  }
  
  void test_eval_with_hessian()
  {
    tester.compare_eval_with_indices_and_eval_with_hessian( mMetric );
    tester.compare_eval_with_grad_and_eval_with_hessian( mMetric );
    tester.compare_analytical_and_numerical_hessians( mMetric );
  }
  
  void test_location_invariant()
  {
    tester.test_location_invariant( mMetric );
    tester.test_grad_location_invariant( mMetric );
    tester.test_hessian_location_invariant( mMetric );
  }
  
  void test_scale_invariant()
  {
    tester.test_scale_invariant( mMetric );
  }
  
  void test_orient_invariant()
  {
    tester.test_orient_invariant( mMetric );
    tester.test_grad_orient_invariant( mMetric );
  }
};

void CompositeMetricTestBase::test_evaluate( bool ideal, EntityTopology type )
{
  MsqPrintError err(std::cout);
  PatchData pd;
  double act, ex;
  bool rval;
  
  if (ideal)
    tester.get_ideal_element( type, false, pd, false );
  else
    tester.get_nonideal_element( type, pd );
  
  rval = mMetric->evaluate( pd, 0, act, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
  
  rval = this->evaluate( pd, 0, ex, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( ex, act, 1e-6 );
}

class AddQualityMetricTest : public CompositeMetricTestBase
{

  CPPUNIT_TEST_SUITE(AddQualityMetricTest);
  CPPUNIT_TEST (test_supported_types);
  CPPUNIT_TEST (test_ideal_element_eval);
  CPPUNIT_TEST (test_ideal_element_grad);
  CPPUNIT_TEST (test_ideal_element_hess);
  CPPUNIT_TEST (test_non_ideal_eval);
  CPPUNIT_TEST (test_valid_hessian);
  CPPUNIT_TEST (test_inverted_elements);
  CPPUNIT_TEST (test_degenerate_elements);
  CPPUNIT_TEST (test_get_evaluations);
  CPPUNIT_TEST (test_get_element_indices);
  CPPUNIT_TEST (test_get_fixed_indices);
  CPPUNIT_TEST (test_eval_with_indices);
  CPPUNIT_TEST (test_eval_with_gradient);
  CPPUNIT_TEST (test_eval_with_hessian);
  CPPUNIT_TEST_SUITE_END();

private:
  MsqError mErr;
  IdealWeightInverseMeanRatio m1;
  IdealWeightInverseMeanRatio m2;
  AddQualityMetric m;
protected:
  virtual bool evaluate( PatchData&, size_t, double&, MsqError& );
public:
  AddQualityMetricTest() 
    : m2(mErr, 2.0), m( &m1, &m2, mErr )
    { mMetric = &m; }
  void setUp() { CPPUNIT_ASSERT(!mErr); }
};

bool AddQualityMetricTest::evaluate( PatchData& pd, size_t h, double& val, MsqError& err )
{
  double v1, v2;
  bool rval = true, rval1;
  
  rval1 = m1.evaluate( pd, h, v1, err );
  MSQ_ERRFALSE(err);
  if (!rval1) rval = false;
  
  rval1 = m2.evaluate( pd, h, v2, err );
  MSQ_ERRFALSE(err);
  if (!rval1) rval = false;
  
  val = v1+v2;
  return rval;
}
  
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(AddQualityMetricTest, "CompositeMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(AddQualityMetricTest, "AddQualityMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(AddQualityMetricTest, "Unit");

class MultiplyQualityMetricTest : public CompositeMetricTestBase
{
private:

  CPPUNIT_TEST_SUITE(MultiplyQualityMetricTest);
  CPPUNIT_TEST (test_supported_types);
  CPPUNIT_TEST (test_ideal_element_eval);
  CPPUNIT_TEST (test_ideal_element_grad);
  CPPUNIT_TEST (test_ideal_element_hess);
  CPPUNIT_TEST (test_non_ideal_eval);
  CPPUNIT_TEST (test_valid_hessian);
  CPPUNIT_TEST (test_inverted_elements);
  CPPUNIT_TEST (test_degenerate_elements);
  CPPUNIT_TEST (test_get_evaluations);
  CPPUNIT_TEST (test_get_element_indices);
  CPPUNIT_TEST (test_get_fixed_indices);
  CPPUNIT_TEST (test_eval_with_indices);
  CPPUNIT_TEST (test_eval_with_gradient);
  CPPUNIT_TEST (test_eval_with_hessian);
  CPPUNIT_TEST_SUITE_END();

  MsqError mErr;
  IdealWeightInverseMeanRatio m1;
  IdealWeightInverseMeanRatio m2;
  MultiplyQualityMetric m;
protected:
  virtual bool evaluate( PatchData&, size_t, double&, MsqError& );
public:
  MultiplyQualityMetricTest() 
    : m2(mErr, 2.0), m( &m1, &m2, mErr )
    { mMetric = &m; }
  void setUp() { CPPUNIT_ASSERT(!mErr); }
};

bool MultiplyQualityMetricTest::evaluate( PatchData& pd, size_t h, double& val, MsqError& err )
{
  double v1, v2;
  bool rval = true, rval1;
  
  rval1 = m1.evaluate( pd, h, v1, err );
  MSQ_ERRFALSE(err);
  if (!rval1) rval = false;
  
  rval1 = m2.evaluate( pd, h, v2, err );
  MSQ_ERRFALSE(err);
  if (!rval1) rval = false;
  
  val = v1*v2;
  return rval;
}

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MultiplyQualityMetricTest, "CompositeMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MultiplyQualityMetricTest, "MultiplyQualityMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MultiplyQualityMetricTest, "Unit");

template <int POWER>
class PowerQualityMetricTest : public CompositeMetricTestBase
{
private:

  CPPUNIT_TEST_SUITE(PowerQualityMetricTest<POWER>);
  CPPUNIT_TEST (test_supported_types);
  CPPUNIT_TEST (test_ideal_element_eval);
  CPPUNIT_TEST (test_ideal_element_grad);
  //CPPUNIT_TEST (test_ideal_element_hess);
  CPPUNIT_TEST (test_non_ideal_eval);
  CPPUNIT_TEST (test_valid_hessian);
  CPPUNIT_TEST (test_inverted_elements);
  CPPUNIT_TEST (test_degenerate_elements);
  CPPUNIT_TEST (test_get_evaluations);
  CPPUNIT_TEST (test_get_element_indices);
  CPPUNIT_TEST (test_get_fixed_indices);
  CPPUNIT_TEST (test_eval_with_indices);
  CPPUNIT_TEST (test_eval_with_gradient);
  CPPUNIT_TEST (test_eval_with_hessian);
  CPPUNIT_TEST_SUITE_END();

  IdealWeightInverseMeanRatio m1;
  PowerQualityMetric m;
protected:
  virtual bool evaluate( PatchData&, size_t, double&, MsqError& );
public:
  PowerQualityMetricTest() 
    : m( &m1, POWER )
    { mMetric = &m; }
};

template <int POWER>
bool PowerQualityMetricTest<POWER>::evaluate( PatchData& pd, size_t h, double& val, MsqError& err )
{
  bool rval = m1.evaluate( pd, h, val, err );
  MSQ_ERRFALSE(err);
  val = std::pow( val, POWER );
  return rval;
}

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PowerQualityMetricTest<1>, "CompositeMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PowerQualityMetricTest<2>, "CompositeMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PowerQualityMetricTest<1>, "PowerQualityMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PowerQualityMetricTest<2>, "PowerQualityMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PowerQualityMetricTest<1>, "Unit");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PowerQualityMetricTest<2>, "Unit");

template <int OFFSET>
class ScalarAddMetricTest : public CompositeMetricTestBase
{
private:

  CPPUNIT_TEST_SUITE(ScalarAddMetricTest<OFFSET>);
  CPPUNIT_TEST (test_supported_types);
  CPPUNIT_TEST (test_ideal_element_eval);
  CPPUNIT_TEST (test_ideal_element_grad);
  CPPUNIT_TEST (test_ideal_element_hess);
  CPPUNIT_TEST (test_non_ideal_eval);
  CPPUNIT_TEST (test_valid_hessian);
  CPPUNIT_TEST (test_inverted_elements);
  CPPUNIT_TEST (test_degenerate_elements);
  CPPUNIT_TEST (test_get_evaluations);
  CPPUNIT_TEST (test_get_element_indices);
  CPPUNIT_TEST (test_get_fixed_indices);
  CPPUNIT_TEST (test_eval_with_indices);
  CPPUNIT_TEST (test_eval_with_gradient);
  CPPUNIT_TEST (test_eval_with_hessian);
  CPPUNIT_TEST_SUITE_END();

  IdealWeightInverseMeanRatio m1;
  ScalarAddQualityMetric m;
protected:
  virtual bool evaluate( PatchData&, size_t, double&, MsqError& );
public:
  ScalarAddMetricTest() 
    : m( &m1, OFFSET )
    { mMetric = &m; }
};

template <int OFFSET>
bool ScalarAddMetricTest<OFFSET>::evaluate( PatchData& pd, size_t h, double& val, MsqError& err )
{
  bool rval = m1.evaluate( pd, h, val, err );
  MSQ_ERRFALSE(err);
  val += OFFSET;
  return rval;
}
  
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ScalarAddMetricTest<0>, "CompositeMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ScalarAddMetricTest<2>, "CompositeMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ScalarAddMetricTest<0>, "ScalarAddMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ScalarAddMetricTest<2>, "ScalarAddMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ScalarAddMetricTest<0>, "Unit");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ScalarAddMetricTest<2>, "Unit");

template <int SCALE>
class ScalarMultiplyMetricTest : public CompositeMetricTestBase
{
private:

  CPPUNIT_TEST_SUITE(ScalarMultiplyMetricTest<SCALE>);
  CPPUNIT_TEST (test_supported_types);
  CPPUNIT_TEST (test_ideal_element_eval);
  CPPUNIT_TEST (test_ideal_element_grad);
  CPPUNIT_TEST (test_ideal_element_hess);
  CPPUNIT_TEST (test_non_ideal_eval);
  CPPUNIT_TEST (test_valid_hessian);
  CPPUNIT_TEST (test_inverted_elements);
  CPPUNIT_TEST (test_degenerate_elements);
  CPPUNIT_TEST (test_get_evaluations);
  CPPUNIT_TEST (test_get_element_indices);
  CPPUNIT_TEST (test_get_fixed_indices);
  CPPUNIT_TEST (test_eval_with_indices);
  CPPUNIT_TEST (test_eval_with_gradient);
  CPPUNIT_TEST (test_eval_with_hessian);
  CPPUNIT_TEST_SUITE_END();

  IdealWeightInverseMeanRatio m1;
  ScalarMultiplyQualityMetric m;
protected:
  virtual bool evaluate( PatchData&, size_t, double&, MsqError& );
public:
  ScalarMultiplyMetricTest() 
    : m( &m1, SCALE )
    { mMetric = &m; }
};

template <int SCALE>
bool ScalarMultiplyMetricTest<SCALE>::evaluate( PatchData& pd, size_t h, double& val, MsqError& err )
{
  bool rval = m1.evaluate( pd, h, val, err );
  MSQ_ERRFALSE(err);
  val *= SCALE;
  return rval;
}
  
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ScalarMultiplyMetricTest<1>, "CompositeMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ScalarMultiplyMetricTest<3>, "CompositeMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ScalarMultiplyMetricTest<1>, "ScalarMultiplyMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ScalarMultiplyMetricTest<3>, "ScalarMultiplyMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ScalarMultiplyMetricTest<1>, "Unit");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ScalarMultiplyMetricTest<3>, "Unit");

