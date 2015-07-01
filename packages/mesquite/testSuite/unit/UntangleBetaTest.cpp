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


/** \file UntangleBetaTest.cpp
 *  \brief unit tests for UntangleBetaQualityMetric
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include "QualityMetricTester.hpp"
#include "PatchData.hpp"
#include "UnitUtil.hpp"

using namespace Mesquite;

class UntangleBetaTest : public CppUnit::TestFixture
{
private:  

  CPPUNIT_TEST_SUITE(UntangleBetaTest);

  CPPUNIT_TEST (test_supported_types);
  CPPUNIT_TEST (test_ideal_element_eval);
  CPPUNIT_TEST (test_inverted_elements);
  CPPUNIT_TEST (test_degenerate_elements);
  CPPUNIT_TEST (test_get_evaluations);
  CPPUNIT_TEST (test_get_element_indices);
  CPPUNIT_TEST (test_get_fixed_indices);
  CPPUNIT_TEST (test_eval_with_indices);
  CPPUNIT_TEST (test_location_invariant);
  CPPUNIT_TEST (test_orient_invariant);
  
  CPPUNIT_TEST_SUITE_END();
  
  UntangleBetaQualityMetric mMetric;
  QualityMetricTester tester;
  
public:

  UntangleBetaTest() 
    : tester(QualityMetricTester::ALL_FE_EXCEPT_SEPTAHEDRON)
    { tester.ideal_pyramid_base_equals_height( true ); }

  void test_supported_types()
    { tester.test_supported_element_types( &mMetric ); }
    
  void test_ideal_element_eval()
  {
    tester.test_evaluate_unit_edge_element( &mMetric, TRIANGLE, 0.0 );
    tester.test_evaluate_unit_edge_element( &mMetric, QUADRILATERAL, 0.0 );
    tester.test_evaluate_unit_edge_element( &mMetric, TETRAHEDRON, 0.0 );
    tester.test_evaluate_unit_edge_element( &mMetric, HEXAHEDRON, 0.0 );
    tester.test_evaluate_unit_edge_element( &mMetric, PRISM, 0.0 );
    tester.test_evaluate_unit_edge_element( &mMetric, PYRAMID, 0.0 );
  }
  
  void test_inverted_elements()
  { 
    MsqPrintError err(std::cout);
    double value;
    PatchData pd;
    char val_str[128];
    
    tester.get_inverted_element( TRIANGLE, pd );
    mMetric.evaluate( pd, 0, value, err );
    ASSERT_NO_ERROR(err);
    sprintf(val_str, "value: %f", value );
    CPPUNIT_ASSERT_MESSAGE(val_str, value > 1e-6);
    
    tester.get_inverted_element( QUADRILATERAL, pd );
    mMetric.evaluate( pd, 0, value, err );
    ASSERT_NO_ERROR(err);
    sprintf(val_str, "value: %f", value );
    CPPUNIT_ASSERT_MESSAGE(val_str, value > 1e-6);
    
    tester.get_inverted_element( TETRAHEDRON, pd );
    mMetric.evaluate( pd, 0, value, err );
    ASSERT_NO_ERROR(err);
    sprintf(val_str, "value: %f", value );
    CPPUNIT_ASSERT_MESSAGE(val_str, value > 1e-6);
    
    tester.get_inverted_element( HEXAHEDRON, pd );
    mMetric.evaluate( pd, 0, value, err );
    ASSERT_NO_ERROR(err);
    sprintf(val_str, "value: %f", value );
    CPPUNIT_ASSERT_MESSAGE(val_str, value > 1e-6);
    
    tester.get_inverted_element( PRISM, pd );
    mMetric.evaluate( pd, 0, value, err );
    ASSERT_NO_ERROR(err);
    sprintf(val_str, "value: %f", value );
    CPPUNIT_ASSERT_MESSAGE(val_str, value > 1e-6);
    
    tester.get_inverted_element( PYRAMID, pd );
    mMetric.evaluate( pd, 0, value, err );
    ASSERT_NO_ERROR(err);
    sprintf(val_str, "value: %f", value );
    CPPUNIT_ASSERT_MESSAGE(val_str, value > 1e-6);
  }
    
  void test_degenerate_elements()
  { 
    MsqPrintError err(std::cout);
    double value;
    PatchData pd;
    char val_str[128];
    
    tester.get_zero_element( TRIANGLE, pd );
    mMetric.evaluate( pd, 0, value, err );
    ASSERT_NO_ERROR(err);
    sprintf(val_str, "value: %f", value );
    CPPUNIT_ASSERT_MESSAGE(val_str, value > 1e-6);
    
    tester.get_zero_element( QUADRILATERAL, pd );
    mMetric.evaluate( pd, 0, value, err );
    ASSERT_NO_ERROR(err);
    sprintf(val_str, "value: %f", value );
    CPPUNIT_ASSERT_MESSAGE(val_str, value > 1e-6);
    
    tester.get_zero_element( TETRAHEDRON, pd );
    mMetric.evaluate( pd, 0, value, err );
    ASSERT_NO_ERROR(err);
    sprintf(val_str, "value: %f", value );
    CPPUNIT_ASSERT_MESSAGE(val_str, value > 1e-6);
    
    tester.get_zero_element( HEXAHEDRON, pd );
    mMetric.evaluate( pd, 0, value, err );
    ASSERT_NO_ERROR(err);
    sprintf(val_str, "value: %f", value );
    CPPUNIT_ASSERT_MESSAGE(val_str, value > 1e-6);
    
    tester.get_zero_element( PRISM, pd );
    mMetric.evaluate( pd, 0, value, err );
    ASSERT_NO_ERROR(err);
    sprintf(val_str, "value: %f", value );
    CPPUNIT_ASSERT_MESSAGE(val_str, value > 1e-6);
    
    tester.get_zero_element( PYRAMID, pd );
    mMetric.evaluate( pd, 0, value, err );
    ASSERT_NO_ERROR(err);
    sprintf(val_str, "value: %f", value );
    CPPUNIT_ASSERT_MESSAGE(val_str, value > 1e-6);
  }
    
  void test_get_evaluations()
    { tester.test_get_element_evaluations( &mMetric ); }
    
  void test_get_element_indices()
    { tester.test_get_element_indices( &mMetric ); }
  
  void test_get_fixed_indices()
    { tester.test_get_indices_fixed( &mMetric ); }
  
  void test_eval_with_indices()
    { tester.compare_eval_and_eval_with_indices( &mMetric ); }
  
  void test_location_invariant()
  {
    tester.test_location_invariant( &mMetric, true );
    tester.test_grad_location_invariant( &mMetric, true );
    tester.test_hessian_location_invariant( &mMetric, true );
  }
  
  void test_orient_invariant()
  {
    tester.test_orient_invariant( &mMetric, true );
    tester.test_grad_orient_invariant( &mMetric, true );
  }
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(UntangleBetaTest, "UntangleBetaTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(UntangleBetaTest, "Unit");
