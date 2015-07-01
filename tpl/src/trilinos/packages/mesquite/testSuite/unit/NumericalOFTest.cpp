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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file NumericalQMTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ObjectiveFunctionTemplate.hpp"
#include "UnitUtil.hpp"
#include "PatchData.hpp"
#include "MsqHessian.hpp"

#include "cppunit/extensions/HelperMacros.h"

using namespace Mesquite;
using namespace std;

const double EPSILON = 1e-6;

class NumericalOFTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE( NumericalOFTest );
  CPPUNIT_TEST( test_gradient_constant );
  CPPUNIT_TEST( test_gradient_linear );
  CPPUNIT_TEST( test_handles_eval_failure );
  CPPUNIT_TEST( test_handles_eval_false );
  CPPUNIT_TEST( test_changed );
  CPPUNIT_TEST( test_unchanged );
  CPPUNIT_TEST( test_Hessian_fails );
  CPPUNIT_TEST_SUITE_END();
  
  PatchData pd;

public:
  
  void setUp();
  void test_gradient_values( bool constant );

  void test_gradient_constant( ) { test_gradient_values( true ); }
  void test_gradient_linear( )   { test_gradient_values( false ); }
  
  void test_handles_eval_failure( );
  void test_handles_eval_false( );
  void test_changed( );
  void test_unchanged( );
  void test_Hessian_fails( );
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(NumericalOFTest, "NumericalOFTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(NumericalOFTest, "Unit");

/** Define a dummy ObjectiveFunction on which to do numerical gradient calculation
 *
 * If 'constant' flag is true, returns a constant value.
 * If 'constant' flag is false, returns a linear function of vertex coordinates.
 */
class NumericalTestOF : public ObjectiveFunctionTemplate 
{
public:
  NumericalTestOF( bool should_fail, bool should_return_false, bool constant_func )
    : changed(false), 
      fail(should_fail), 
      return_false(should_return_false), 
      constant(constant_func),
      linearGrad(1,2,3)
    {}
  
  bool evaluate( EvalType type, PatchData& pd, double& val, bool free, MsqError& );
  
  ObjectiveFunction* clone() const { return new NumericalTestOF( *this ); }
  
  void clear() { changed = true; }
  
  bool changed;      // "accumulated value" has changed (tester should initialize)
  bool fail;         // if true, calls to accumulate will unconditionally fail
  bool return_false; // if true, evaluate will return false.
  bool constant;     // if true, OF value is constant, otherwise linear.
  
  const Vector3D linearGrad;
};


bool NumericalTestOF::evaluate( EvalType type, PatchData& pd, double& val, bool free, MsqError& err )
{
  if (fail) {
    MSQ_SETERR(err)( "Expected failure of OF::evaluate", MsqError::INVALID_ARG );
    return false;
  }
  
  if (pd.num_free_vertices() < 1) {
    MSQ_SETERR(err)("PatchData without free vertices.", MsqError::INVALID_ARG);
    return false;
  }
  
  if (type != ObjectiveFunction::CALCULATE &&
      type != ObjectiveFunction::TEMPORARY)
    changed = true;
  
  if (constant) {
    val = 1.0;
    return !return_false;
  }
  
  val = 0.0;
  for (size_t i = 0; i < pd.num_nodes(); ++i)
  {
    const MsqVertex& v = pd.vertex_by_index(i);
    val += linearGrad[0]*v[0] + linearGrad[1]*v[1] + linearGrad[2]*v[2];
  }
  
  return !return_false;
}


void NumericalOFTest::setUp()
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
  pd.fill( 6, coords, 
           4, TRIANGLE, tri_conn,
           fixed_vtx, err );
  ASSERT_NO_ERROR(err);

}


void NumericalOFTest::test_gradient_values( bool constant )
{
  NumericalTestOF func( false, false, constant );
  MsqPrintError err(std::cout);
  double value;
  vector<Vector3D> gradient;
  
  bool rval = func.evaluate_with_gradient( ObjectiveFunction::CALCULATE,
                                           pd, value, gradient, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
  
  Vector3D expected = constant ? Vector3D( 0, 0, 0 ) : func.linearGrad;
  CPPUNIT_ASSERT( gradient.size() == pd.num_free_vertices() );
  for (vector<Vector3D>::iterator i = gradient.begin(); i != gradient.end(); ++i)
    CPPUNIT_ASSERT_VECTORS_EQUAL( expected, *i, EPSILON );
}


void NumericalOFTest::test_handles_eval_failure( )
{
  MsqError err;
  NumericalTestOF func( true, false, true );
  double value;
  vector<Vector3D> gradient;
  
  func.evaluate_with_gradient( ObjectiveFunction::CALCULATE,
                               pd, value, gradient, err );
  CPPUNIT_ASSERT( err );  
}


void NumericalOFTest::test_handles_eval_false( )
{
  MsqError err;
  NumericalTestOF func( false, true, true );
  double value;
  vector<Vector3D> gradient;
  
  bool rval = func.evaluate_with_gradient( ObjectiveFunction::CALCULATE,
                                           pd, value, gradient, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT( !rval );  
}


void NumericalOFTest::test_changed( )
{
  MsqPrintError err(cout);
  NumericalTestOF func( false, false, true );
  double value;
  vector<Vector3D> gradient;
  bool rval;
  
  func.changed = false;
  rval = func.evaluate_with_gradient( ObjectiveFunction::SAVE,
                                      pd, value, gradient, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT( rval );  
  CPPUNIT_ASSERT( func.changed );
  
  func.changed = false;
  rval = func.evaluate_with_gradient( ObjectiveFunction::UPDATE,
                                      pd, value, gradient, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT( rval );  
  CPPUNIT_ASSERT( func.changed );
}



void NumericalOFTest::test_unchanged( )
{
  MsqPrintError err(cout);
  NumericalTestOF func( false, false, true );
  double value;
  vector<Vector3D> gradient;
  bool rval;
  
  func.changed = false;
  rval = func.evaluate_with_gradient( ObjectiveFunction::TEMPORARY,
                                      pd, value, gradient, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT( rval );  
  CPPUNIT_ASSERT( !func.changed );
}


void NumericalOFTest::test_Hessian_fails( )
{
  MsqError err;
  NumericalTestOF func( false, false, true );
  double value;
  vector<Vector3D> gradient;
  MsqHessian Hessian;
  
  func.evaluate_with_Hessian( ObjectiveFunction::CALCULATE,
                               pd, value, gradient, Hessian, err );
  CPPUNIT_ASSERT( err );
}

  
  


