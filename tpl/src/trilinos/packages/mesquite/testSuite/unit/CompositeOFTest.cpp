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


/** \file CompositeOFTest.cpp
 *  \brief Unit tests for Composite Objective Functions
 *  \author Jason Kraftcheck 
 */



#include "Mesquite.hpp"
#include "Mesquite_CompositeOFAdd.hpp"
#include "Mesquite_CompositeOFMultiply.hpp"
#include "Mesquite_CompositeOFScalarAdd.hpp"
#include "Mesquite_CompositeOFScalarMultiply.hpp"
#include "Mesquite_MsqHessian.hpp"
#include "Mesquite_IdealWeightInverseMeanRatio.hpp"
#include "Mesquite_LPtoPTemplate.hpp"

#include "ObjectiveFunctionTests.hpp"
#include "PatchDataInstances.hpp"
#include "Mesquite_cppunit/extensions/HelperMacros.h"
#include "UnitUtil.hpp"
#include "Mesquite_Matrix3D.hpp"
#include <iterator>

using namespace Mesquite;
using std::cout;
using std::endl;
using std::cerr;

/** Fake ObjectiveFunction to pass to Composite OFs */
class FauxObjectiveFunction : public ObjectiveFunction
{
  public:
    FauxObjectiveFunction( double value, bool invalid = false, bool error = false) 
      : mValue(value), mInvalid(invalid), mError(error)
      { ++instanceCount; }
    ~FauxObjectiveFunction() { --instanceCount; }
    bool initialize_block_coordinate_descent( MeshDomainAssoc*, const Settings*, PatchSet*, MsqError& )
      { CPPUNIT_ASSERT_MESSAGE("This shouldn't ever get called", false ); return false; }
    bool evaluate( EvalType, PatchData&, double& value_out, bool, MsqError& err )
      { 
        if (mError)
          MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
        value_out = mValue; 
        return !mInvalid;
      }
    ObjectiveFunction* clone() const 
      { 
        ++instanceCount;
        return new FauxObjectiveFunction(*this); 
      }
    void clear() {}
    int min_patch_layers() const { return 0; }
    
    void initialize_queue( MeshDomainAssoc* , 
                           const Settings* ,
                           MsqError&  ) {}
  
    double get_value() const { return mValue; }
    static int get_instance_count() { return instanceCount; }
  private:
    double mValue;
    bool mInvalid, mError;
    static int instanceCount;
};
int FauxObjectiveFunction::instanceCount = 0;

class CompositeOFTest : public CppUnit::TestFixture, public ObjectiveFunctionTests
{
private:
  CPPUNIT_TEST_SUITE(CompositeOFTest);
  
  CPPUNIT_TEST (test_add_value);
  CPPUNIT_TEST (test_multiply_value);
  CPPUNIT_TEST (test_scalar_add_value);
  CPPUNIT_TEST (test_scalar_multiply_value);

  CPPUNIT_TEST (test_add_gradient);
  CPPUNIT_TEST (test_multiply_gradient);
  CPPUNIT_TEST (test_scalar_add_gradient);
  CPPUNIT_TEST (test_scalar_multiply_gradient);

  CPPUNIT_TEST (test_add_hess_diagonal);
  CPPUNIT_TEST (test_multiply_hess_diagonal);
  CPPUNIT_TEST (test_scalar_add_hess_diagonal);
  CPPUNIT_TEST (test_scalar_multiply_hess_diagonal);

  CPPUNIT_TEST (test_add_hessian);
  CPPUNIT_TEST (test_multiply_hessian);
  CPPUNIT_TEST (test_scalar_add_hessian);
  CPPUNIT_TEST (test_scalar_multiply_hessian);

  CPPUNIT_TEST (test_clone_add);
  CPPUNIT_TEST (test_clone_multiply);
  CPPUNIT_TEST (test_clone_scalar_add);
  CPPUNIT_TEST (test_clone_scalar_multiply);
  
  CPPUNIT_TEST (test_add_invalid);
  CPPUNIT_TEST (test_multiply_invalid);
  CPPUNIT_TEST (test_scalar_add_invalid);
  CPPUNIT_TEST (test_scalar_multiply_invalid);
  
  CPPUNIT_TEST (test_add_error);
  CPPUNIT_TEST (test_multiply_error);
  CPPUNIT_TEST (test_scalar_add_error);
  CPPUNIT_TEST (test_scalar_multiply_error);

  CPPUNIT_TEST_SUITE_END();
  
  FauxObjectiveFunction OF1, OF2, OF3, OF4, invalidOF, errorOF;
  IdealWeightInverseMeanRatio metric;
  LPtoPTemplate LP1, LP2;
  
public:
  
  CompositeOFTest() 
    : OF1(1.0), OF2(3.0), OF3(-7.0), OF4(M_PI),
      invalidOF(1.0,true,false), errorOF(1.0,false,true),
      LP1( 1, &metric ), LP2( 2, &metric )
   {}
  
  void test_add_value();
  void test_multiply_value();
  void test_scalar_add_value();
  void test_scalar_multiply_value();

  void test_add_gradient();
  void test_multiply_gradient();
  void test_scalar_add_gradient();
  void test_scalar_multiply_gradient();

  void test_add_hess_diagonal();
  void test_multiply_hess_diagonal();
  void test_scalar_add_hess_diagonal();
  void test_scalar_multiply_hess_diagonal();

  void test_add_hessian();
  void test_multiply_hessian();
  void test_scalar_add_hessian();
  void test_scalar_multiply_hessian();
  
  void test_clone_add();
  void test_clone_multiply();
  void test_clone_scalar_add();
  void test_clone_scalar_multiply();
  
  void test_add_invalid();
  void test_multiply_invalid();
  void test_scalar_add_invalid();
  void test_scalar_multiply_invalid();
  
  void test_add_error();
  void test_multiply_error();
  void test_scalar_add_error();
  void test_scalar_multiply_error();

  void test_evaluate( double expected_value, ObjectiveFunction& of );
  void get_hessians( MsqHessian& LP1_hess, MsqHessian& LP2_hess,
                     ObjectiveFunction& OF, MsqHessian& OF_hess );
  void test_composite_clone( ObjectiveFunction& of );
  void test_invalid_eval( ObjectiveFunction& of );
  void test_eval_fails( ObjectiveFunction& of );
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(CompositeOFTest, "CompositeOFTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(CompositeOFTest, "Unit");

void CompositeOFTest::test_evaluate( double expected, ObjectiveFunction& OF )
{
  MsqPrintError err(cout);
  double value;
  bool rval = OF.evaluate( ObjectiveFunction::CALCULATE, patch(), value, false, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, value, 1e-6 );
}

void CompositeOFTest::test_add_value()
{
  CompositeOFAdd add1( &OF1, &OF2 );
  test_evaluate( OF1.get_value() + OF2.get_value(), add1 );

  CompositeOFAdd add2( &OF3, &OF4 );
  test_evaluate( OF3.get_value() + OF4.get_value(), add2 );
}

void CompositeOFTest::test_multiply_value()
{
  CompositeOFMultiply mult1( &OF1, &OF2 );
  test_evaluate( OF1.get_value() * OF2.get_value(), mult1 );

  CompositeOFMultiply mult2( &OF3, &OF4 );
  test_evaluate( OF3.get_value() * OF4.get_value(), mult2 );
}

void CompositeOFTest::test_scalar_add_value()
{
  CompositeOFScalarAdd add1( sqrt(2.0), &OF1 );
  test_evaluate( OF1.get_value() + sqrt(2.0), add1 );

  CompositeOFScalarAdd add2( -1.0, &OF4 );
  test_evaluate( OF4.get_value() - 1, add2 );
}

void CompositeOFTest::test_scalar_multiply_value()
{
  CompositeOFScalarMultiply mult1( sqrt(2.0), &OF1 );
  test_evaluate( OF1.get_value() * sqrt(2.0), mult1 );

  CompositeOFScalarMultiply mult2( -1.0, &OF4 );
  test_evaluate( -OF4.get_value(), mult2 );
}

void CompositeOFTest::test_add_gradient()
  { CompositeOFAdd OF( &LP1, &LP2 ); compare_numerical_gradient( &OF ); }
void CompositeOFTest::test_multiply_gradient()
  { CompositeOFMultiply OF( &LP1, &LP2 ); compare_numerical_gradient( &OF ); }
void CompositeOFTest::test_scalar_add_gradient()
  { CompositeOFScalarAdd OF( M_PI, &LP1 ); compare_numerical_gradient( &OF ); }
void CompositeOFTest::test_scalar_multiply_gradient()
  { CompositeOFScalarMultiply OF( M_PI, &LP1 ); compare_numerical_gradient( &OF ); }

void CompositeOFTest::get_hessians( MsqHessian& LP1_hess, 
                                    MsqHessian& LP2_hess,
                                    ObjectiveFunction& OF, 
                                    MsqHessian& OF_hess )
{
  MsqPrintError err(cout);
  PatchData pd;
  create_twelve_hex_patch( pd, err ); 
  ASSERT_NO_ERROR( err );

  LP1_hess.initialize( pd, err ); ASSERT_NO_ERROR(err);
  LP2_hess.initialize( pd, err ); ASSERT_NO_ERROR(err);
  OF_hess .initialize( pd, err ); ASSERT_NO_ERROR(err);
  
  std::vector<Vector3D> grad;
  bool rval;
  double value;
  rval = LP1.evaluate_with_Hessian( ObjectiveFunction::CALCULATE, pd, value, grad, LP1_hess, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
  rval = LP2.evaluate_with_Hessian( ObjectiveFunction::CALCULATE, pd, value, grad, LP2_hess, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
  rval = OF .evaluate_with_Hessian( ObjectiveFunction::CALCULATE, pd, value, grad, OF_hess , err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
}

void CompositeOFTest::test_add_hess_diagonal()
{
  CompositeOFAdd OF( &LP1, &LP2 );
  compare_hessian_diagonal( &OF );
}

void CompositeOFTest::test_multiply_hess_diagonal()
{
  CompositeOFMultiply OF( &LP1, &LP2 );
  std::vector<SymMatrix3D> hess1, hess2, hess;
  
  MsqPrintError err(cout);
  PatchData pd;
  create_twelve_hex_patch( pd, err ); 
  ASSERT_NO_ERROR( err );

  std::vector<Vector3D> grad1, grad2, grad;
  bool rval;
  double value1, value2, value;
  rval = LP1.evaluate_with_Hessian_diagonal( ObjectiveFunction::CALCULATE, pd, value1, grad1, hess1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
  rval = LP2.evaluate_with_Hessian_diagonal( ObjectiveFunction::CALCULATE, pd, value2, grad2, hess2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
  rval = OF .evaluate_with_Hessian_diagonal( ObjectiveFunction::CALCULATE, pd, value, grad, hess , err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);

  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), grad1.size() );
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), grad2.size() );
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), grad .size() );

  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), hess1.size() );
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), hess2.size() );
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), hess .size() );
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( value1 * value2, value, 1e-6 );
  
  for (size_t i = 0; i < pd.num_free_vertices(); ++i) {
    const Vector3D expected_grad = value2 * grad1[i] + value1 * grad2[i];
    CPPUNIT_ASSERT_VECTORS_EQUAL( expected_grad, grad[i], 1e-6 );
    
    Matrix3D o;
    o.outer_product( grad1[i], grad2[i] );
    Matrix3D expect = o + transpose(o);
    expect += value2 * hess1[i];
    expect += value1 * hess2[i];
    CPPUNIT_ASSERT_MATRICES_EQUAL( expect, Matrix3D(hess[i]), 1e-6 );
  }
}

void CompositeOFTest::test_scalar_add_hess_diagonal()
{
  CompositeOFScalarAdd OF( 1111.1, &LP1 );
  compare_hessian_diagonal( &OF );
}

void CompositeOFTest::test_scalar_multiply_hess_diagonal()
{
  const double scale = 2.5;
  CompositeOFScalarMultiply OF( scale, &LP1 );
  compare_hessian_diagonal( &OF );
}

void CompositeOFTest::test_add_hessian()
{
    // test value and gradient 
  CompositeOFAdd OF( &LP1, &LP2 );
  compare_hessian_gradient( &OF );
  
    // test actual hessian values
  MsqHessian hess1, hess2, hess;
  get_hessians( hess1, hess2, OF, hess );
  Matrix3D *b1, *b2, *b;
  for (unsigned r = 0; r < hess.size(); ++r) {
    for (unsigned c = r; c < hess.size(); ++c) {
      b1 = hess1.get_block(r,c);
      b2 = hess2.get_block(r,c);
      b  = hess .get_block(r,c);
      if (b) {
        CPPUNIT_ASSERT_MATRICES_EQUAL( *b1 + *b2, *b, 1e-6 );
      }
    }
  }
}

void CompositeOFTest::test_multiply_hessian()
{
  MsqError err;
  PatchData pd;
  create_twelve_hex_patch( pd, err ); 
  ASSERT_NO_ERROR( err );
  
  // this should always fail because the Hessian is not sparse
  CompositeOFMultiply OF( &LP1, &LP2 );
  double value;
  MsqHessian hess;
  hess.initialize( pd, err );
  ASSERT_NO_ERROR(err);
  std::vector<Vector3D> grad;
  OF.evaluate_with_Hessian( ObjectiveFunction::CALCULATE, pd, value, grad, hess, err );
  CPPUNIT_ASSERT(err);
}

void CompositeOFTest::test_scalar_add_hessian()
{
    // test value and gradient 
  CompositeOFScalarAdd OF( 1111.1, &LP1 );
  compare_hessian_gradient( &OF );
  
    // test actual hessian values
  MsqHessian hess1, hess2, hess;
  get_hessians( hess1, hess2, OF, hess );
  Matrix3D *b1, *b;
  for (unsigned r = 0; r < hess.size(); ++r) {
    for (unsigned c = r; c < hess.size(); ++c) {
      b1 = hess1.get_block(r,c);
      b  = hess .get_block(r,c);
      if (b) {
        CPPUNIT_ASSERT_MATRICES_EQUAL( *b1, *b, 1e-6 );
      }
    }
  }
}

void CompositeOFTest::test_scalar_multiply_hessian()
{
    // test value and gradient 
  const double scale = 2.5;
  CompositeOFScalarMultiply OF( scale, &LP1 );
  compare_hessian_gradient( &OF );
  
    // test actual hessian values
  MsqHessian hess1, hess2, hess;
  get_hessians( hess1, hess2, OF, hess );
  Matrix3D *b1, *b;
  for (unsigned r = 0; r < hess.size(); ++r) {
    for (unsigned c = r; c < hess.size(); ++c) {
      b1 = hess1.get_block(r,c);
      b  = hess .get_block(r,c);
      if (b) {
        CPPUNIT_ASSERT_MATRICES_EQUAL( scale * *b1, *b, 1e-6 );
      }
    }
  }
}


void CompositeOFTest::test_composite_clone( ObjectiveFunction& OF )
{
  // save current count of instances of underlying OFs for later
  const int init_count = FauxObjectiveFunction::get_instance_count();

  // clone the objective function
  ObjectiveFunction* clone = OF.clone();
  
  // check that the underlying OFs were also cloned
  CPPUNIT_ASSERT( init_count < FauxObjectiveFunction::get_instance_count() );
  
  // check that the value is the same
  MsqPrintError err(cout);
  double orig_val, clone_val;
  bool rval;
  rval = OF.evaluate( ObjectiveFunction::CALCULATE, patch(), orig_val, false, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
  rval = clone->evaluate( ObjectiveFunction::CALCULATE, patch(), clone_val, false, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( orig_val, clone_val, 1e-6 );
  
  // check that cloned instances of underlying OFs are deleted
  delete clone;
  CPPUNIT_ASSERT_EQUAL( init_count, FauxObjectiveFunction::get_instance_count() );
}
  

void CompositeOFTest::test_clone_add()
  { CompositeOFAdd OF( &OF1, &OF2 ); test_composite_clone( OF ); }
void CompositeOFTest::test_clone_multiply()
  { CompositeOFMultiply OF( &OF1, &OF2 ); test_composite_clone( OF ); }
void CompositeOFTest::test_clone_scalar_add()
  { CompositeOFScalarAdd OF( 2.1, &OF2 ); test_composite_clone( OF ); }
void CompositeOFTest::test_clone_scalar_multiply()
  { CompositeOFScalarMultiply OF( 0.333, &OF2 ); test_composite_clone( OF ); }

void CompositeOFTest::test_invalid_eval( ObjectiveFunction& OF )
{
  MsqPrintError err(cout);
  bool rval;
  double value;
  rval = OF.evaluate( ObjectiveFunction::CALCULATE, patch(), value, false, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(rval == false);
}

void CompositeOFTest::test_eval_fails( ObjectiveFunction& OF )
{
  MsqError err;
  double value;
  OF.evaluate( ObjectiveFunction::CALCULATE, patch(), value, false, err );
  CPPUNIT_ASSERT_EQUAL(MsqError::INTERNAL_ERROR, err.error_code());
}

void CompositeOFTest::test_add_invalid()
{
  CompositeOFAdd add1( &OF1, &invalidOF );
  test_invalid_eval( add1 );
  
  CompositeOFAdd add2( &invalidOF, &OF3 );
  test_invalid_eval( add2 );
}

void CompositeOFTest::test_multiply_invalid()
{
  CompositeOFMultiply mult1( &OF1, &invalidOF );
  test_invalid_eval( mult1 );
  
  CompositeOFMultiply mult2( &invalidOF, &OF3 );
  test_invalid_eval( mult2 );
}

void CompositeOFTest::test_scalar_add_invalid()
{
  CompositeOFScalarAdd OF( 2.0, &invalidOF );
  test_invalid_eval( OF );
}

void CompositeOFTest::test_scalar_multiply_invalid()
{
  CompositeOFScalarMultiply OF( 2.0, &invalidOF );
  test_invalid_eval( OF );
}
  
void CompositeOFTest::test_add_error()
{
  CompositeOFAdd add1( &OF1, &errorOF );
  test_eval_fails( add1 );
  
  CompositeOFAdd add2( &errorOF, &OF3 );
  test_eval_fails( add2 );
}

void CompositeOFTest::test_multiply_error()
{
  CompositeOFMultiply mult1( &OF1, &errorOF );
  test_eval_fails( mult1 );
  
  CompositeOFMultiply mult2( &errorOF, &OF3 );
  test_eval_fails( mult2 );
}

void CompositeOFTest::test_scalar_add_error()
{
  CompositeOFScalarAdd OF( 2.0, &errorOF );
  test_eval_fails( OF );
}

void CompositeOFTest::test_scalar_multiply_error()
{
  CompositeOFScalarMultiply OF( 2.0, &errorOF );
  test_eval_fails( OF );
}
  
