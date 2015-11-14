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


/** \file ObjectiveFunctionTests.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "ObjectiveFunctionTests.hpp"


/** get a patchdata to use for testing */
static bool init_pd( PatchData& pd );
PatchData& ObjectiveFunctionTests::patch() {
  static PatchData the_pd;
  static bool did_init = init_pd( the_pd );
  CPPUNIT_ASSERT(did_init);
  return the_pd;
}
static bool init_pd( PatchData& pd ) {
  MsqError err;
  const double coords[] = { 0,0,0, 1,1,1, 2,2,2 };
  const bool fixed[] = { false, true, true };
  const size_t indices[] = { 0, 1, 2 };
  pd.fill( 3, coords, 1, TRIANGLE, indices, fixed, err );
  return !err;
}

/** Internal helper function for test_eval_type */
double ObjectiveFunctionTests::evaluate_internal( 
                                 ObjectiveFunction::EvalType type, 
                                 OFTestMode test_mode,
                                 ObjectiveFunction* of )
{
  MsqPrintError err(cout);
  vector<Vector3D> grad;
  vector<SymMatrix3D> diag;
  MsqHessian hess;
  bool valid = false;
  double result;
  
  switch (test_mode) {
    case EVAL:
      valid = of->evaluate( type, patch(), result, OF_FREE_EVALS_ONLY, err );
      break;
    case GRAD:
      valid = of->evaluate_with_gradient( type, patch(), result, grad, err );
      break;
    case DIAG:
      valid = of->evaluate_with_Hessian_diagonal( type, patch(), result, grad, diag, err );
      break;
    case HESS:
      hess.initialize( patch(), err );
      ASSERT_NO_ERROR( err );
      valid = of->evaluate_with_Hessian( type, patch(), result, grad, hess, err );
      break;
    default:
      CPPUNIT_ASSERT(false);
  }
  
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT(valid);
  return result;
}

void ObjectiveFunctionTests::test_eval_type( 
                            ObjectiveFunction::EvalType type,
                            OFTestMode test_mode,
                            ObjectiveFunctionTemplate* of )
{
    // define two sets of quality metric values
  const double vals1[] = { 1.0, 3.0, 4.0, 8.0 };
  const size_t vals1_len = sizeof(vals1)/sizeof(vals1[0]);
  const double vals2[] = { 21.5, 11.1, 30.0, 0.5 };
  const size_t vals2_len = sizeof(vals2)/sizeof(vals2[0]);
  
    // create a quality metric to use
  OFTestQM metric;
  of->set_quality_metric( &metric );
  
    // get some initial values to compare to
  of->clear();
  metric.set_values( vals1, vals1_len );
  const double init1 = evaluate_internal( ObjectiveFunction::CALCULATE, test_mode, of );
  of->clear();
  metric.set_values( vals2, vals2_len );
  const double init2 = evaluate_internal( ObjectiveFunction::CALCULATE, test_mode, of );
  of->clear();
  metric.append_values( vals1, vals1_len );
  const double inits = evaluate_internal( ObjectiveFunction::CALCULATE, test_mode, of );
  of->clear();

  double val, expected;  
  switch (type) {
    case ObjectiveFunction::CALCULATE:
    
        // first make sure we get back the same values
      metric.set_values( vals1, vals1_len );
      val = evaluate_internal( ObjectiveFunction::CALCULATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( init1, val, 1e-6 );
      metric.set_values( vals2, vals2_len );
      val = evaluate_internal( ObjectiveFunction::CALCULATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( init2, val, 1e-6 );
      
        // now do something that should modify the accumulated value of the OF
      evaluate_internal( ObjectiveFunction::ACCUMULATE, test_mode, of );
      
        // check that the values are unchanged
      metric.set_values( vals1, vals1_len );
      val = evaluate_internal( ObjectiveFunction::CALCULATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( init1, val, 1e-6 );
      metric.set_values( vals2, vals2_len );
      val = evaluate_internal( ObjectiveFunction::CALCULATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( init2, val, 1e-6 );
      
      break;
      
    case ObjectiveFunction::ACCUMULATE:
      
        // begin with first set
      metric.set_values( vals1, vals1_len );
      val = evaluate_internal( ObjectiveFunction::ACCUMULATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( init1, val, 1e-6 );
        
        // add in second set
      metric.set_values( vals2, vals2_len );
      val = evaluate_internal( ObjectiveFunction::ACCUMULATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( inits, val, 1e-6 );
      
        // clear 
      of->clear();
      
        // begin with second set
      metric.set_values( vals2, vals2_len );
      val = evaluate_internal( ObjectiveFunction::ACCUMULATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( init2, val, 1e-6 );
      
        // add in first set
      metric.set_values( vals1, vals1_len );
      val = evaluate_internal( ObjectiveFunction::ACCUMULATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( inits, val, 1e-6 );
      
      break;
      
    case ObjectiveFunction::SAVE:
      
        // calculate value for first set twice
      metric.set_values( vals1, vals1_len );
      metric.append_values( vals1, vals1_len );
      expected = evaluate_internal( ObjectiveFunction::CALCULATE, test_mode, of );
      
        // begin with the first set
      metric.set_values( vals1, vals1_len );
      val = evaluate_internal( ObjectiveFunction::ACCUMULATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( init1, val, 1e-6 );
      
        // add the second set
      metric.set_values( vals2, vals2_len );
      val = evaluate_internal( ObjectiveFunction::ACCUMULATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( inits, val, 1e-6 );
      
        // now save the second set of values - OF value should not change
      metric.set_values( vals2, vals2_len );
      val = evaluate_internal( ObjectiveFunction::SAVE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( inits, val, 1e-6 );
      
        // now replace the second set with the first
      metric.set_values( vals1, vals1_len );
      val = evaluate_internal( ObjectiveFunction::UPDATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, val, 1e-6 );
      
        // check that saved values are cleared
      of->clear();
      metric.set_values( vals2, vals2_len );
      val = evaluate_internal( ObjectiveFunction::UPDATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( init2, val, 1e-6 );
      
      break;
 
    case ObjectiveFunction::UPDATE:
      
        // With no saved data, an update should produce the same
        // result as CALCULATE
      metric.set_values( vals1, vals1_len );
      val = evaluate_internal( ObjectiveFunction::UPDATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( init1, val, 1e-6 );

        // Doing an update with a second patch should change
        // the global value to that of the second set
      metric.set_values( vals2, vals2_len );
      val = evaluate_internal( ObjectiveFunction::UPDATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( init2, val, 1e-6 );
      
        // Now add in the second set again
      metric.set_values( vals2, vals2_len );
      evaluate_internal( ObjectiveFunction::ACCUMULATE, test_mode, of );
      
        // Now replace one accumulation of the second set with the first
        // Should result in the first set + the second set
      metric.set_values( vals1, vals1_len );
      val = evaluate_internal( ObjectiveFunction::UPDATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( inits, val, 1e-6 );
     
      break;

    case ObjectiveFunction::TEMPORARY:
      
        // With no saved data, an TEMPORARY should produce the same
        // result as CALCULATE
      metric.set_values( vals1, vals1_len );
      val = evaluate_internal( ObjectiveFunction::TEMPORARY, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( init1, val, 1e-6 );
      
        // Begin with two instances of the second set in the 
        // accumulated value, with one instance saved so it 
        // can be removed later
      metric.set_values( vals2, vals2_len );
      evaluate_internal( ObjectiveFunction::ACCUMULATE, test_mode, of );
      evaluate_internal( ObjectiveFunction::UPDATE, test_mode, of );
      
        // Now do a temporary eval, replacing one instance of the
        // second set wiht the first set
      metric.set_values( vals1, vals1_len );
      val = evaluate_internal( ObjectiveFunction::TEMPORARY, test_mode, of );
      
        // TEMPORARY should produce the same value as UPDATE, but without
        // modifying any internal state.  
      expected = evaluate_internal( ObjectiveFunction::UPDATE, test_mode, of );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, val, 1e-6 );
      
      break;
      
    default:
      CPPUNIT_ASSERT_MESSAGE("No test for specified evaluation type", false);
      break;
  }
}

void ObjectiveFunctionTests::test_value( 
                        const double* input_values,
                        unsigned num_input_values,
                        double expected_value,
                        OFTestMode test_mode,
                        ObjectiveFunctionTemplate* of )
{
  OFTestQM metric( input_values, num_input_values );
  of->set_quality_metric( &metric );
  
  double val = evaluate_internal( ObjectiveFunction::CALCULATE, test_mode, of );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( expected_value, val, 1e-6 );
}

void ObjectiveFunctionTests::test_clone( ObjectiveFunctionTemplate* of )
{
  const double some_vals[] = { 1, 2, 3, 4, 5, 6, 0.5 };
  const unsigned num_vals = sizeof(some_vals)/sizeof(some_vals[0]);
  OFTestQM metric( some_vals, num_vals );
  of->set_quality_metric(&metric);
  
    // test that we get the same value from both
  auto_ptr<ObjectiveFunction> of2( of->clone() );
  double exp_val, val;
  exp_val = evaluate_internal( ObjectiveFunction::CALCULATE, EVAL, of );
  val = evaluate_internal( ObjectiveFunction::CALCULATE, EVAL, of2.get() );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_val, val, 1e-12 );
  
    // check if OF supports BCD -- if not then done
  MsqError err;
  of->evaluate( ObjectiveFunction::UPDATE, patch(), val, false, err );
  if (err) {
    err.clear();
    return;
  }
  
    // build up some saved state in the objective function
  of->clear();
  const double vals1[] = { 1.0, 3.0, 4.0, 8.0 };
  const size_t vals1_len = sizeof(vals1)/sizeof(vals1[0]);
  const double vals2[] = { 21.5, 11.1, 30.0, 0.5 };
  const size_t vals2_len = sizeof(vals2)/sizeof(vals2[0]);
  metric.set_values( vals1, vals1_len );
  evaluate_internal( ObjectiveFunction::SAVE, EVAL, of );
  metric.append_values( vals2, vals2_len );
  evaluate_internal( ObjectiveFunction::ACCUMULATE, EVAL, of );
  
    // check that clone has same accumulated data
  of2 = auto_ptr<ObjectiveFunction>(of->clone());
  metric.set_values( some_vals, num_vals );
  exp_val = evaluate_internal( ObjectiveFunction::UPDATE, EVAL, of );
  val = evaluate_internal( ObjectiveFunction::UPDATE, EVAL, of2.get() );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_val, val, 1e-12 );
}

/** Test correct handling of QM negate flag */
void ObjectiveFunctionTests::test_negate_flag( 
                                       OFTestMode test_mode, 
                                       ObjectiveFunctionTemplate* of )
{
  const double some_vals[] = { 1, 2, 3, 4, 5, 6, 0.5 };
  const unsigned num_vals = sizeof(some_vals)/sizeof(some_vals[0]);
  OFTestQM metric( some_vals, num_vals );
  of->set_quality_metric(&metric);

  MsqPrintError err(cout);
  bool rval = false;
  double value[2];
  vector<Vector3D> grad[2];
  vector<SymMatrix3D> diag[2];
  MsqHessian hess[2];
  
    // Do twice, once w/out negate flag set and then once
    // with negate flag == -1.
  ObjectiveFunction::EvalType type = ObjectiveFunction::CALCULATE;
  for (unsigned i = 0; i < 2; ++i ) {
    switch (test_mode) {
      case EVAL:
        rval = of->evaluate( type, patch(), value[i], false, err );
        break;
      case GRAD:
        rval = of->evaluate_with_gradient( type, patch(), value[i], grad[i], err );
        break;
      case DIAG:
        rval = of->evaluate_with_Hessian_diagonal( type, patch(), value[i], grad[i], diag[i], err );
        break;
      case HESS:
        hess[i].initialize(patch(),err);
        ASSERT_NO_ERROR( err );
        rval = of->evaluate_with_Hessian( type, patch(), value[i], grad[i], hess[i], err );
        break;
      default:
        CPPUNIT_ASSERT_MESSAGE("Invalid enum value in test code",false);
        break;
    }
    ASSERT_NO_ERROR( err );
    CPPUNIT_ASSERT(rval);
    metric.set_negate_flag(-1);
  }
  
  switch (test_mode) {
    case HESS:
      CPPUNIT_ASSERT_EQUAL( hess[0].size(), hess[1].size() );
      for (size_t r = 0; r < hess[0].size(); ++r)
        for (size_t c = r; c < hess[0].size(); ++c)
          if (hess[0].get_block(r,c))
            CPPUNIT_ASSERT_MATRICES_EQUAL( -*hess[0].get_block(r,c),
                                           *hess[1].get_block(r,c),
                                           1e-6 );
    case DIAG:
      // NOTE: When case HESS: falls through to here, diag[0] and diag[1]
      // will be empty, making this a no-op.
      CPPUNIT_ASSERT_EQUAL( diag[0].size(), diag[1].size() );
      for (size_t j = 0; j < diag[0].size(); ++j) 
        CPPUNIT_ASSERT_MATRICES_EQUAL( -diag[0][j], diag[1][j], 1e-6 );
    case GRAD:
      CPPUNIT_ASSERT_EQUAL( grad[0].size(), grad[1].size() );
      for (size_t j = 0; j < grad[0].size(); ++j)
        CPPUNIT_ASSERT_VECTORS_EQUAL( -grad[0][j], grad[1][j], 1e-6 );
    default:
      CPPUNIT_ASSERT_DOUBLES_EQUAL( -value[0], value[1], 1e-6 );
  }
}

void ObjectiveFunctionTests::compare_numerical_gradient( ObjectiveFunction* of )
{
  MsqPrintError err(std::cout);
  PatchData pd;
  create_twelve_hex_patch( pd, err ); 
  ASSERT_NO_ERROR( err );
  
  std::vector<Vector3D> num_grad, ana_grad;
  double num_val, ana_val;
  bool valid;
  
  valid = of->evaluate_with_gradient( ObjectiveFunction::CALCULATE, pd, ana_val, ana_grad, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), ana_grad.size() );
  
  valid = of->ObjectiveFunction::evaluate_with_gradient( ObjectiveFunction::CALCULATE, pd, num_val, num_grad, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), num_grad.size() );
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( ana_val, num_val, 1e-6 );
  for (size_t i = 0; i < pd.num_free_vertices(); ++i) {
    CPPUNIT_ASSERT_VECTORS_EQUAL( num_grad[i], ana_grad[i], 1e-3 );
  }
}

void ObjectiveFunctionTests::compare_hessian_gradient( ObjectiveFunction* of )
{
  MsqPrintError err(std::cout);
  PatchData pd;
  create_twelve_hex_patch( pd, err ); 
  ASSERT_NO_ERROR( err );
  
  std::vector<Vector3D> grad, hess_grad;
  MsqHessian hess;
  double grad_val, hess_val;
  bool valid;
  
  valid = of->evaluate_with_gradient( ObjectiveFunction::CALCULATE, pd, grad_val, grad, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), grad.size() );
  
  hess.initialize( pd, err );
  ASSERT_NO_ERROR( err );
  valid = of->evaluate_with_Hessian( ObjectiveFunction::CALCULATE, pd, hess_val, hess_grad, hess, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), hess_grad.size() );
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( grad_val, hess_val, 1e-6 );
  for (size_t i = 0; i < pd.num_free_vertices(); ++i) {
    CPPUNIT_ASSERT_VECTORS_EQUAL( grad[i], hess_grad[i], 1e-6 );
  }
}

void ObjectiveFunctionTests::compare_diagonal_gradient( ObjectiveFunction* of )
{
  MsqPrintError err(std::cout);
  PatchData pd;
  create_twelve_hex_patch( pd, err ); 
  ASSERT_NO_ERROR( err );
  
  std::vector<Vector3D> grad, hess_grad;
  std::vector<SymMatrix3D> hess;
  double grad_val, hess_val;
  bool valid;
  
  valid = of->evaluate_with_gradient( ObjectiveFunction::CALCULATE, pd, grad_val, grad, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), grad.size() );
  
  valid = of->evaluate_with_Hessian_diagonal( ObjectiveFunction::CALCULATE, pd, hess_val, hess_grad, hess, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), hess_grad.size() );
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( grad_val, hess_val, 1e-6 );
  for (size_t i = 0; i < pd.num_free_vertices(); ++i) {
    CPPUNIT_ASSERT_VECTORS_EQUAL( grad[i], hess_grad[i], 1e-6 );
  }
}

void ObjectiveFunctionTests::compare_hessian_diagonal( ObjectiveFunction* of )
{
  MsqPrintError err(std::cout);
  PatchData pd;
  create_twelve_hex_patch( pd, err ); 
  ASSERT_NO_ERROR( err );
  
  std::vector<Vector3D> diag_grad, hess_grad;
  std::vector<SymMatrix3D> diag;
  MsqHessian hess;
  double diag_val, hess_val;
  bool valid;
  
  valid = of->evaluate_with_Hessian_diagonal( ObjectiveFunction::CALCULATE, pd, diag_val, diag_grad, diag, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), diag_grad.size() );
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), diag.size() );
  
  hess.initialize( pd, err );
  ASSERT_NO_ERROR( err );
  valid = of->evaluate_with_Hessian( ObjectiveFunction::CALCULATE, pd, hess_val, hess_grad, hess, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), hess_grad.size() );
  CPPUNIT_ASSERT_EQUAL( pd.num_free_vertices(), hess.size() );
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( hess_val, diag_val, 1e-6 );
  for (size_t i = 0; i < pd.num_free_vertices(); ++i) {
    CPPUNIT_ASSERT_VECTORS_EQUAL( hess_grad[i], diag_grad[i], 1e-6 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( *hess.get_block(i,i), diag[i], 1e-6 );
  }
}

double MAT_EPS(const Matrix3D& A) 
{ return std::max( 0.001 * sqrt(Frobenius_2(A)), 0.001 ); }

#define CHECK_EQUAL_MATRICES(A,B) \
  CPPUNIT_ASSERT_MATRICES_EQUAL( (A), (B), MAT_EPS(A) )

void ObjectiveFunctionTests::compare_numerical_hessian( ObjectiveFunction* of,
                                                        bool diagonal_only )
{
  const double delta = 0.0001;

  MsqPrintError err(std::cout);
  PatchData pd;
  create_qm_two_tet_patch( pd, err ); 
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT( pd.num_free_vertices() != 0 );
  
    // get analytical Hessian from objective function
  std::vector<Vector3D> grad;
  std::vector<SymMatrix3D> diag;
  MsqHessian hess;
  hess.initialize( pd, err );
  ASSERT_NO_ERROR( err );
  double value;
  bool valid;
  if (diagonal_only)
    valid = of->evaluate_with_Hessian_diagonal( ObjectiveFunction::CALCULATE, pd, value, grad, diag, err );
  else
    valid = of->evaluate_with_Hessian( ObjectiveFunction::CALCULATE, pd, value, grad, hess, err );
  ASSERT_NO_ERROR(err); CPPUNIT_ASSERT(valid);

  
    // do numerical approximation of each block and compare to analytical value
  for (size_t i = 0; i < pd.num_free_vertices(); ++i) {
    const size_t j_end = diagonal_only ? i+1 : pd.num_free_vertices();
    for (size_t j = i; j < j_end; ++j) {
        // do numerical approximation for block corresponding to
        // coorindates for ith and jth vertices.
      Matrix3D block;    
      for (int k = 0; k < 3; ++k) {
        for (int m = 0; m < 3; ++m) {
          double dk, dm, dkm;
          Vector3D ik = pd.vertex_by_index(i);
          Vector3D im = pd.vertex_by_index(j);
          
          Vector3D delta_k(0.0); delta_k[k] = delta;
          pd.move_vertex( delta_k, i, err ); ASSERT_NO_ERROR(err);
          valid = of->evaluate( ObjectiveFunction::CALCULATE, pd, dk, true, err );
          ASSERT_NO_ERROR(err); CPPUNIT_ASSERT(valid);
          
          Vector3D delta_m(0.0); delta_m[m] = delta;
          pd.move_vertex( delta_m, j, err ); ASSERT_NO_ERROR(err);
          valid = of->evaluate( ObjectiveFunction::CALCULATE, pd, dkm, true, err );
          ASSERT_NO_ERROR(err); CPPUNIT_ASSERT(valid);
          
            // be careful here that we do the right thing if i==j
          pd.set_vertex_coordinates( ik, i, err ); ASSERT_NO_ERROR(err);
          pd.set_vertex_coordinates( im, j, err ); ASSERT_NO_ERROR(err);
          pd.move_vertex( delta_m, j, err ); ASSERT_NO_ERROR(err);
          valid = of->evaluate( ObjectiveFunction::CALCULATE, pd, dm, true, err );
          ASSERT_NO_ERROR(err); CPPUNIT_ASSERT(valid);
          
          pd.set_vertex_coordinates( ik, i, err ); ASSERT_NO_ERROR(err);
          pd.set_vertex_coordinates( im, j, err ); ASSERT_NO_ERROR(err);
          
          block[k][m] = (dkm - dk - dm + value)/(delta*delta);
        }
      }
        // compare to analytical value
      if (diagonal_only) {
        CPPUNIT_ASSERT(i == j); // see j_end above
        CPPUNIT_ASSERT(i < diag.size());
        CHECK_EQUAL_MATRICES( block, Matrix3D(diag[i]) );
      }
      else {
        Matrix3D* m = hess.get_block( i, j );
        Matrix3D* mt = hess.get_block( j, i );
        if (NULL != m) {
          CHECK_EQUAL_MATRICES( block, *m );
        }
        if (NULL != mt) {
          CHECK_EQUAL_MATRICES( transpose(block), *m );
        }
        if (NULL == mt && NULL == m) {
          CHECK_EQUAL_MATRICES( Matrix3D(0.0), block );
        }
      }
    }
  }
}
  


const size_t HANDLE_VAL = 0xDEADBEEF;
class OFTestBadQM : public QualityMetric
{
  public:
    
    OFTestBadQM( bool return_error )  : mError(return_error) {}
     
    virtual MetricType get_metric_type() const 
      { return ELEMENT_BASED; }
    
    virtual string get_name() const
      { return "ObjectiveFunctionTests"; }
    
    virtual int get_negate_flag() const
      { return 1; }
    
    virtual void get_evaluations( PatchData&, vector<size_t>& h, bool, MsqError& )
      {
        h.clear();
        h.push_back(HANDLE_VAL);
      }
    
    virtual bool evaluate( PatchData&, size_t h, double&, MsqError& err )
      {
        CPPUNIT_ASSERT_EQUAL( HANDLE_VAL, h );
        if (mError) 
          MSQ_SETERR(err)(MsqError::INVALID_STATE);
        return false;
      }
    
    virtual bool evaluate_with_indices( PatchData&, size_t h, double&, vector<size_t>&, MsqError& err )
      {
        CPPUNIT_ASSERT_EQUAL( HANDLE_VAL, h );
        if (mError) 
          MSQ_SETERR(err)(MsqError::INVALID_STATE);
        return false;
      }
    
    virtual bool evaluate_with_gradient( PatchData&, size_t h, double&, vector<size_t>&, vector<Vector3D>&, MsqError& err )
      {
        CPPUNIT_ASSERT_EQUAL( HANDLE_VAL, h );
        if (mError) 
          MSQ_SETERR(err)(MsqError::INVALID_STATE);
        return false;
      }

    virtual bool evaluate_with_Hessian( PatchData&, size_t h, double&, vector<size_t>&, vector<Vector3D>&, vector<Matrix3D>&, MsqError& err )
      {
        CPPUNIT_ASSERT_EQUAL( HANDLE_VAL, h );
        if (mError) 
          MSQ_SETERR(err)(MsqError::INVALID_STATE);
        return false;
      }
  
  private:
    
    bool mError;
};

void ObjectiveFunctionTests::test_handles_invalid_qm( 
                                     OFTestMode test_mode,
                                     ObjectiveFunctionTemplate* of )
{
  OFTestBadQM metric(false);
  of->set_quality_metric( &metric );
  
  MsqPrintError err(cout);
  vector<Vector3D> grad;
  vector<SymMatrix3D> diag;
  MsqHessian hess;
  double result;
  bool valid = false;
  
  switch (test_mode) {
    case EVAL:
      valid = of->evaluate( ObjectiveFunction::CALCULATE, patch(), result, OF_FREE_EVALS_ONLY, err );
      break;
    case GRAD:
      valid = of->evaluate_with_gradient( ObjectiveFunction::CALCULATE, patch(), result, grad, err );
      break;
    case DIAG:
      valid = of->evaluate_with_Hessian_diagonal( ObjectiveFunction::CALCULATE, patch(), result, grad, diag, err );
      break;
    case HESS:
      hess.initialize( patch(), err );
      ASSERT_NO_ERROR( err );
      valid = of->evaluate_with_Hessian( ObjectiveFunction::CALCULATE, patch(), result, grad, hess, err );
      break;
    default:
      CPPUNIT_ASSERT(false);
  }
  
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT(!valid);
}
  
                                     
void ObjectiveFunctionTests::test_handles_qm_error( 
                                   OFTestMode test_mode,
                                   ObjectiveFunctionTemplate* of )
  
{
  OFTestBadQM metric(true);
  of->set_quality_metric( &metric );
  
  MsqError err;
  vector<Vector3D> grad;
  vector<SymMatrix3D> diag;
  MsqHessian hess;
  double result;
  bool valid;
  
  switch (test_mode) {
    case EVAL:
      valid = of->evaluate( ObjectiveFunction::CALCULATE, patch(), result, OF_FREE_EVALS_ONLY, err );
      break;
    case GRAD:
      valid = of->evaluate_with_gradient( ObjectiveFunction::CALCULATE, patch(), result, grad, err );
      break;
    case DIAG:
      valid = of->evaluate_with_Hessian_diagonal( ObjectiveFunction::CALCULATE, patch(), result, grad, diag, err );
      break;
    case HESS:
      hess.initialize( patch(), err );
      ASSERT_NO_ERROR( err );
      valid = of->evaluate_with_Hessian( ObjectiveFunction::CALCULATE, patch(), result, grad, hess, err );
      break;
    default:
      CPPUNIT_ASSERT(false);
  }
  
  CPPUNIT_ASSERT(err);
}
