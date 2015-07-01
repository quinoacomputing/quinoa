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


/** \file EdgeLengthMetricTest.cpp
 *  \brief unit tests for EdgeLengthQualityMetric
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "EdgeLengthMetric.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include "QualityMetricTester.hpp"

using namespace Mesquite;

class EdgeLengthMetricTest : public CppUnit::TestFixture
{
private:  

  CPPUNIT_TEST_SUITE(EdgeLengthMetricTest);

  CPPUNIT_TEST (test_supported_types);
  CPPUNIT_TEST (test_ideal_element_eval);
  CPPUNIT_TEST (test_valid_hessian);
  CPPUNIT_TEST (test_get_evaluations);
  CPPUNIT_TEST (test_get_edge_indices);
  CPPUNIT_TEST (test_eval_with_indices);
  CPPUNIT_TEST (test_eval_with_gradient);
  CPPUNIT_TEST (test_location_invariant);
  CPPUNIT_TEST (test_orient_invariant);
  
  CPPUNIT_TEST_SUITE_END();
  
  EdgeLengthMetric mMetric;
  QualityMetricTester tester;
  
public:

  EdgeLengthMetricTest() 
    : mMetric(0.0), tester(QualityMetricTester::ALL_FE_EXCEPT_SEPTAHEDRON)
    { }

  void test_supported_types()
    { tester.test_supported_element_types( &mMetric ); }
    
  void test_ideal_element_eval()
  {
    tester.test_evaluate_unit_edge_element( &mMetric, TRIANGLE,      1.0 );
    tester.test_evaluate_unit_edge_element( &mMetric, QUADRILATERAL, 1.0 );
    tester.test_evaluate_unit_edge_element( &mMetric, TETRAHEDRON,   1.0 );
    tester.test_evaluate_unit_edge_element( &mMetric, PYRAMID,       1.0 );
    tester.test_evaluate_unit_edge_element( &mMetric, PRISM,         1.0 );
    tester.test_evaluate_unit_edge_element( &mMetric, HEXAHEDRON,    1.0 );
  }
    
  void test_valid_hessian()
    { tester.test_symmetric_Hessian_diagonal_blocks( &mMetric ); }
    
  void test_get_evaluations()
    { tester.test_get_edge_evaluations( &mMetric ); }
    
  void test_get_edge_indices()
    { tester.test_get_edge_indices( &mMetric ); }
  
  void test_eval_with_indices()
    { tester.compare_eval_and_eval_with_indices( &mMetric ); }
    
  void test_eval_with_gradient()
  {
    tester.compare_eval_with_indices_and_eval_with_gradient( &mMetric );
    tester.compare_analytical_and_numerical_gradients( &mMetric );
  }
  
  void test_location_invariant()
  {
    tester.test_location_invariant( &mMetric );
    tester.test_grad_location_invariant( &mMetric );
  }
  
  void test_orient_invariant()
  {
    tester.test_orient_invariant( &mMetric );
    tester.test_grad_orient_invariant( &mMetric );
  }
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(EdgeLengthMetricTest, "EdgeLengthMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(EdgeLengthMetricTest, "Unit");
