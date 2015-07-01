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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TargetMetricTest.hpp
 *  \brief Templatized common code for testing various target metric
 *         implementation types.
 *  \author Jason Kraftcheck 
 */

#include "UnitUtil.hpp"
#include "MsqError.hpp"
#include "MsqMatrix.hpp"

#include "TMetric.hpp"
#include "TMetricBarrier.hpp"
#include "AWMetric.hpp"
#include "AWMetricBarrier.hpp"

// NOTE: Caller must define TARGET_TEST_GROUP to be a quoted string,
//       typically the base file name of the file containing the 
//       calls to TEST_METRIC_*

// Macro arguments:
//  shape_invariant
//  size_invariant
//  orient_invariant
//  barrier

#define REGISTER_BASE_TESTS \
  CPPUNIT_TEST (test_ideal_eval); \
  CPPUNIT_TEST (test_ideal_gradient); \
  CPPUNIT_TEST (test_inverted); \
  CPPUNIT_TEST (test_shape); \
  CPPUNIT_TEST (test_scale); \
  CPPUNIT_TEST (test_orient)

#define REGISTER_GRAD_TESTS \
  CPPUNIT_TEST (compare_eval_and_eval_with_grad); \
  CPPUNIT_TEST (compare_anaytic_and_numeric_grads) 

#define REGISTER_HESS_TESTS \
  CPPUNIT_TEST (compare_eval_with_grad_and_eval_with_hess); \
  CPPUNIT_TEST (compare_anaytic_and_numeric_hess) 

#define BEGIN_TEST_DECL( METRIC, DIM, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ) \
class METRIC ## _ ## DIM ## DTest : public TMetricTest< METRIC, DIM > { public: \
  METRIC ## _ ## DIM ## DTest () : TMetricTest< METRIC, DIM >( (SHAPE_INVAR), (SIZE_INVAR), (ORIENT_INVAR), (BARRIER), IDEAL ) {} \
  CPPUNIT_TEST_SUITE( METRIC ## _ ## DIM ## DTest )

#define END_TEST_DECL(SUITE, DIM, METRIC) \
  CPPUNIT_TEST_SUITE_END(); \
}; \
CPPUNIT_NS::AutoRegisterSuite< METRIC ## _ ## DIM ## DTest > METRIC ## _ ## DIM ## D_UnitRegister ("Unit"); \
CPPUNIT_NS::AutoRegisterSuite< METRIC ## _ ## DIM ## DTest > METRIC ## _ ## DIM ## D_FileRegister (TARGET_TEST_GROUP); \
CPPUNIT_NS::AutoRegisterSuite< METRIC ## _ ## DIM ## DTest > METRIC ## _ ## DIM ## D_BaseRegister ( #SUITE "Test" )


/** Register tests for metric with no derivative implementations */
#define TEST_METRIC_NO_DERIVS_2D( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ) \
  BEGIN_TEST_DECL( METRIC, 2, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ); \
  REGISTER_BASE_TESTS; \
  END_TEST_DECL(METRIC,2,METRIC)

/** Register tests for metric with no derivative implementations */
#define TEST_METRIC_NO_DERIVS_3D( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ) \
  BEGIN_TEST_DECL( METRIC, 3, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIE, IDEALR ); \
  REGISTER_BASE_TESTS; \
  END_TEST_DECL(METRIC,3,METRIC)

/** Register tests for metric with no derivative implementations */
#define TEST_METRIC_NO_DERIVS( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ) \
  TEST_METRIC_NO_DERIVS_2D( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ); \
  TEST_METRIC_NO_DERIVS_3D( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL )

/** Register tests for metric with implementation of analytic gradient */
#define TEST_METRIC_WITH_GRAD_2D( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ) \
  BEGIN_TEST_DECL( METRIC, 2, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ); \
  REGISTER_BASE_TESTS; \
  REGISTER_GRAD_TESTS; \
  END_TEST_DECL(METRIC,2,METRIC)

/** Register tests for metric with implementation of analytic gradient */
#define TEST_METRIC_WITH_GRAD_3D( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ) \
  BEGIN_TEST_DECL( METRIC, 3, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ); \
  REGISTER_BASE_TESTS; \
  REGISTER_GRAD_TESTS; \
  END_TEST_DECL(METRIC,3,METRIC)

/** Register tests for metric with implementation of analytic gradient */
#define TEST_METRIC_WITH_GRAD( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ) \
  TEST_METRIC_WITH_GRAD_2D( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ); \
  TEST_METRIC_WITH_GRAD_3D( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL )

#define TEST_NAMED_METRIC_WITH_HESS_2D( METRIC, NAME, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ) \
  BEGIN_TEST_DECL( METRIC, 2, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ); \
  REGISTER_BASE_TESTS; \
  REGISTER_GRAD_TESTS; \
  REGISTER_HESS_TESTS; \
  END_TEST_DECL(NAME,2,METRIC)

#define TEST_NAMED_METRIC_WITH_HESS_3D( METRIC, NAME, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ) \
  BEGIN_TEST_DECL( METRIC, 3, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ); \
  REGISTER_BASE_TESTS; \
  REGISTER_GRAD_TESTS; \
  REGISTER_HESS_TESTS; \
  END_TEST_DECL(NAME,3,METRIC)

/** Register tests for metric with implementation of analytic gradient and Hessian */
#define TEST_NAMED_METRIC_WITH_HESS( METRIC, NAME, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ) \
  TEST_NAMED_METRIC_WITH_HESS_2D( METRIC, NAME, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ); \
  TEST_NAMED_METRIC_WITH_HESS_3D( METRIC, NAME, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL )

/** Register tests for metric with implementation of analytic gradient and Hessian */
#define TEST_METRIC_WITH_HESS_2D( METRIC,         SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ) \
  TEST_NAMED_METRIC_WITH_HESS_2D( METRIC, METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL )

/** Register tests for metric with implementation of analytic gradient and Hessian */
#define TEST_METRIC_WITH_HESS_3D( METRIC,         SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ) \
  TEST_NAMED_METRIC_WITH_HESS_3D( METRIC, METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL )

/** Register tests for metric with implementation of analytic gradient and Hessian */
#define TEST_METRIC_WITH_HESS( METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ) \
  TEST_NAMED_METRIC_WITH_HESS_2D( METRIC, METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL ); \
  TEST_NAMED_METRIC_WITH_HESS_3D( METRIC, METRIC, SHAPE_INVAR, SIZE_INVAR, ORIENT_INVAR, BARRIER, IDEAL )

#define TEST_NON_QUALITY_METRIC_WITH_HESS_2D( METRIC ) \
  BEGIN_TEST_DECL( METRIC, 2, true, true, true, true, 0.0 ); \
  REGISTER_GRAD_TESTS; \
  REGISTER_HESS_TESTS; \
  END_TEST_DECL(NAME,2,METRIC)

#define TEST_NON_QUALITY_METRIC_WITH_HESS_3D( METRIC ) \
  BEGIN_TEST_DECL( METRIC, 3, true, true, true, true, 0.0 ); \
  REGISTER_GRAD_TESTS; \
  REGISTER_HESS_TESTS; \
  END_TEST_DECL(NAME,3,METRIC)

/** Regsiter tests for a metric that doesn't really measure quality */
#define TEST_NON_QUALITY_METRIC_WITH_HESS( METRIC ) \
  TEST_NON_QUALITY_METRIC_WITH_HESS_2D( METRIC ); \
  TEST_NON_QUALITY_METRIC_WITH_HESS_3D( METRIC ); 

using namespace Mesquite;

const double Avals[][9] = { {0},
                            {2},
                            {2, 1,  // 2x2 values 
                             1, 2}, 
                            {2, 1, 1,  // 3x3 values
                             1, 2, 1, 
                             1, 1, 2} 
                          };
const double Bvals[][9] = { {0},
                            {-0.1},
                            {-0.1,  -0.15,   // 2x2 values
                             -0.25, -0.8 }, 
                            { 1.5, -0.7, -0.8,  // 3x3 values  
                              0.8, -1.3, -0.7, 
                              0.6, -0.9, -2.0} 
                          };
const double Cvals[][9] = { {0},
                            {0.5},
                            {-1.0, 0.5,  // 2x2 values
                              0.0, 1.0 },
                            { 0.5, 0.0,  0.1, // 3x3 values
                              0.5, 1.0,  0.1,
                              0.0, 0.0, -1.5 }
                          };

/**\brief Common tests for all target metric types
 *
 * Commont test framework for implementations of the following types:
 * \c TMetric , \c AWMetric
 */    
template <class Metric,unsigned DIM> class TMetricTest : public CppUnit::TestFixture 
{

  private:

    Metric testMetric;
    const double idealVal;
    const bool shapeInvariant, sizeInvariant, orientInvariant, Barrier;

  public:
    
    typedef MsqMatrix<DIM,DIM> Matrix;
    
    TMetricTest( bool shape_invariant,
                 bool size_invariant,
                 bool orient_invariant,
                 bool barrier,
                 double ideal_val )
      : idealVal(ideal_val),
        shapeInvariant(shape_invariant),
        sizeInvariant(size_invariant),
        orientInvariant(orient_invariant),
        Barrier(barrier),
        Zero(0.0),
        I(1.0),
        A(Avals[DIM]),
        B(Bvals[DIM]),
        C(Cvals[DIM])
      {}
    
      // Some initial matrix values used in many tests
    const Matrix Zero, I, A, B, C;

    void test_ideal_eval();
    void test_ideal_gradient();
    void test_inverted();
    void test_shape();
    void test_scale();
    void test_orient();

    void compare_anaytic_and_numeric_grads();
    void compare_anaytic_and_numeric_hess();
    void compare_eval_and_eval_with_grad();
    void compare_eval_with_grad_and_eval_with_hess();

  private:
    /**\brief Test if metric is or is not sensitive to difference between A and W
     *
     * Given an active matrix A and a target matrix W, test whether or
     * not the metric is sensitive to the difference.  Fail if actual
     * sensitivity to difference is not equal to expected sensitivity
     * passed as the first argument
     */
    void test_non_ideal( bool sensitive, Matrix A, Matrix W );
    

/*************************************************************************
 *               Use overloaded function names to do the stuff
 *               that is different for different base metric types
 *************************************************************************/
 
   // TMetric
  inline bool eval( TMetric& metric, MsqMatrix<DIM,DIM> A, MsqMatrix<DIM,DIM> W, 
                    double& value, MsqError& err )
    {
      return metric.evaluate( A*inverse(W), value, err );
    }
  inline bool grad( TMetric& metric, MsqMatrix<DIM,DIM> A, MsqMatrix<DIM,DIM> W, 
                    double& value, MsqMatrix<DIM,DIM>& dmdA, MsqError& err )
    {
      bool rval = metric.evaluate_with_grad( A*inverse(W), value, dmdA, err );
      dmdA = dmdA * transpose(inverse(W)); 
      return rval; 
    }
  inline bool num_grad( TMetric& metric, MsqMatrix<DIM,DIM> A, MsqMatrix<DIM,DIM> W, 
                        double& value, MsqMatrix<DIM,DIM>& dmdA, MsqError& err )
    {
      bool rval = metric.evaluate_with_grad( A*inverse(W), value, dmdA, err );
      dmdA = dmdA * transpose(inverse(W)); 
      return rval; 
    }
  inline bool hess( TMetric& metric, MsqMatrix<DIM,DIM> A, MsqMatrix<DIM,DIM> W, 
                    double& value, MsqMatrix<DIM,DIM>& dmdA, MsqMatrix<DIM,DIM> d2mdA2[3], MsqError& err )
    {
      bool rval = metric.evaluate_with_hess( A*inverse(W), value, dmdA, d2mdA2, err );
      dmdA = dmdA * transpose(inverse(W)); 
      for (unsigned i = 0; i < DIM*(DIM+1)/2; ++i) d2mdA2[i] = inverse(W) * d2mdA2[i] * transpose(inverse(W));
      return rval; 
    }
  inline bool num_hess( TMetric& metric, MsqMatrix<DIM,DIM> A, MsqMatrix<DIM,DIM> W, 
                        double& value, MsqMatrix<DIM,DIM>& dmdA, MsqMatrix<DIM,DIM> d2mdA2[3], MsqError& err )
    {
      bool rval = metric.evaluate_with_hess( A*inverse(W), value, dmdA, d2mdA2, err );
      dmdA = dmdA * transpose(inverse(W)); 
      for (unsigned i = 0; i < DIM*(DIM+1)/2; ++i) d2mdA2[i] = inverse(W) * d2mdA2[i] * transpose(inverse(W));
      return rval;
    }
 
   // AWMetric
  inline bool eval( AWMetric& metric, MsqMatrix<DIM,DIM> A, MsqMatrix<DIM,DIM> W, 
                    double& value, MsqError& err )
    { 
      bool rval = metric.evaluate( A, W, value, err );  
      return rval;
    }
  inline bool grad( AWMetric& metric, MsqMatrix<DIM,DIM> A, MsqMatrix<DIM,DIM> W, 
                    double& value, MsqMatrix<DIM,DIM>& dmdA, MsqError& err )
    { 
      bool rval = metric.evaluate_with_grad( A, W, value, dmdA, err );   
      return rval;
    }
  inline bool num_grad( AWMetric& metric, MsqMatrix<DIM,DIM> A, MsqMatrix<DIM,DIM> W, 
                    double& value, MsqMatrix<DIM,DIM>& dmdA, MsqError& err )
    { 
      bool rval = metric.evaluate_with_grad( A, W, value, dmdA, err );
      return rval;  
    }
  inline bool hess( AWMetric& metric, MsqMatrix<DIM,DIM> A, MsqMatrix<DIM,DIM> W, 
                    double& value, MsqMatrix<DIM,DIM>& dmdA, MsqMatrix<DIM,DIM> d2mdA2[3], MsqError& err )
    { 
      bool rval = metric.evaluate_with_hess( A, W, value, dmdA, d2mdA2, err );
      return rval;  
    }      
  inline bool num_hess( AWMetric& metric, MsqMatrix<DIM,DIM> A, MsqMatrix<DIM,DIM> W, 
                    double& value, MsqMatrix<DIM,DIM>& dmdA, MsqMatrix<DIM,DIM> d2mdA2[3], MsqError& err )
    { 
      bool  rval = metric.evaluate_with_hess( A, W, value, dmdA, d2mdA2, err ); 
      return rval;  
    }            
 };

#define TMETRIC_FUNC template <class Metric, unsigned DIM> void TMetricTest<Metric,DIM>
#define MAT_TYPE TMetricTest<Metric,DIM>::Matrix

/*************************************************************************
 *          Implement actual (templatized) test code
 *************************************************************************/

TMETRIC_FUNC::test_ideal_eval()
{
  MsqPrintError err(std::cerr);
  double val, eps = 5e-5;
  bool valid;
  
  valid = eval( testMetric, I, I, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( idealVal, val, eps );
  
  valid = eval( testMetric, A, A, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( idealVal, val, eps );
  
  valid = eval( testMetric, B, B, val, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( idealVal, val, eps );
}

TMETRIC_FUNC::test_ideal_gradient()
{
  MsqPrintError err(std::cerr);
  MsqMatrix<DIM,DIM> g;
  double val, eps = 5e-3;
  bool valid;
  
  valid = grad( testMetric, I, I, val, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  ASSERT_MATRICES_EQUAL( Zero, g, eps );
  
  valid = grad( testMetric, A, A, val, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  ASSERT_MATRICES_EQUAL( Zero, g, eps );
  
  valid = grad( testMetric, B, B, val, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  ASSERT_MATRICES_EQUAL( Zero, g, eps );
}

TMETRIC_FUNC::test_inverted() 
{
  MsqPrintError err(std::cerr);
  MsqMatrix<DIM,DIM> V( 1.0 ), W( 1.0 ), g, h[6];
  V(DIM-1,DIM-1) = -1.0;
  double val;
  bool valid;
  
  if (Barrier) {
    valid = eval( testMetric, V, W, val, err );
    if (err.error_code() == err.BARRIER_VIOLATED)
      err.clear();
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(!valid);
    
    valid = grad( testMetric, V, W, val, g, err );
    if (err.error_code() == err.BARRIER_VIOLATED)
      err.clear();
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(!valid);
    
    valid = hess( testMetric, V, W, val, g, h, err );
    if (err.error_code() == err.BARRIER_VIOLATED)
      err.clear();
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(!valid);
  }
  else {
    valid = eval( testMetric, V, W, val, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(valid);
    CPPUNIT_ASSERT( val > idealVal );
    
    valid = grad( testMetric, V, W, val, g, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(valid);
    CPPUNIT_ASSERT( val > idealVal );
    CPPUNIT_ASSERT( sqr_Frobenius(g) > 1e-6 );
  }
}
  
TMETRIC_FUNC::test_non_ideal( bool sensitive,
                              Matrix J,
                              Matrix W )
{
  MsqPrintError err(std::cerr);
  MsqMatrix<DIM,DIM> g;
  double val, eps = 1e-5;
  bool valid;
  if (!sensitive) {
    valid = eval( testMetric, J, W, val, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(valid);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( idealVal, val, eps );

    valid = grad( testMetric, J, W, val, g, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(valid);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( idealVal, val, eps );
    ASSERT_MATRICES_EQUAL( (MsqMatrix<DIM,DIM>(0.0)), g, eps );
  }
  else {
    valid = eval( testMetric, J, W, val, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(valid);
    CPPUNIT_ASSERT( val > idealVal );
  }
}

TMETRIC_FUNC::test_shape()
{
  const double r3 = sqrt(3.0);
  const double U_vals[][9] ={ { 2/r3, 1/r3,
                                1/r3, 2/r3 },
                              { 2/r3, 1/r3, 0, 
                                1/r3, 2/r3, 0, 
                                0,    0,    1 } };
  Matrix U( U_vals[DIM-2] ), W( 1.0 );
  test_non_ideal( !shapeInvariant, U, W );
}
  
TMETRIC_FUNC::test_scale() 
{
  Matrix L( 2.0 ), W( 1.0 );
  test_non_ideal( !sizeInvariant, L, W );
}
  
TMETRIC_FUNC::test_orient() 
{
  const double V_vals[][9] = { { 0, -1, 
                                 1,  0 },
                               { 0, -1,  0, 
                                 1,  0,  0, 
                                 0,  0,  1 } };
  Matrix V( V_vals[DIM-2] ), W( 1.0 );
  test_non_ideal( !orientInvariant, V, W );
}

static double releps( double a ) { return std::max(1e-6,1e-8*fabs(a)); }

TMETRIC_FUNC::compare_eval_and_eval_with_grad()
{
  MsqError err;
  Matrix g;
  bool valid;
  double gv, v;
  
  valid = grad( testMetric, I, A, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = eval( testMetric, I, A, v, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v, gv, releps(v) );
  
  valid = grad( testMetric, A, B, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = eval( testMetric, A, B, v, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v, gv, releps(v) );
  
    // also test inverted for non-barrier metrics
  if (Barrier)
    return;
  
  valid = grad( testMetric, C, I, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = eval( testMetric, C, I, v, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v, gv, releps(v) );
}

TMETRIC_FUNC::compare_eval_with_grad_and_eval_with_hess()
{
  MsqError err;
  Matrix g, h, H[DIM*(DIM+1)/2];
  bool valid;
  double gv, hv;
  
  valid = grad( testMetric, I, A, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = hess( testMetric, I, A, hv, h, H, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( gv, hv, releps(gv) );
  ASSERT_MATRICES_EQUAL( g, h, 1e-5 );
  
  valid = grad( testMetric, A, B, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = hess( testMetric, A, B, hv, h, H, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( gv, hv, releps(gv) );
  ASSERT_MATRICES_EQUAL( g, h, 1e-5 );
  
    // also test inverted for non-barrier metrics
  if (Barrier)
    return;
  
  valid = grad( testMetric, C, I, gv, g, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = hess( testMetric, C, I, hv, h, H, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( gv, hv, releps(gv) );
  ASSERT_MATRICES_EQUAL( g, h, 1e-5 );
}

template <typename M>
double eps_mat( const M& mu ) { return std::max( Frobenius(mu)*1e-2, 1e-4 ); }

TMETRIC_FUNC::compare_anaytic_and_numeric_grads()
{
  const double EPS_VAL = 1e-6;
  
  MsqError err;
  Matrix num, ana;
  bool valid;
  double nval, aval;
  
  Matrix D(I);
  D(0,0) += 1e-5;
  valid = num_grad( testMetric, D, I, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = grad( testMetric, D, I, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps_mat(num) );
  
  valid = num_grad( testMetric, I, A, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = grad( testMetric, I, A, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps_mat(num) );
  
  valid = num_grad( testMetric, A, I, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = grad( testMetric, A, I, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps_mat(num) );
  
  valid = num_grad( testMetric, I, B, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = grad( testMetric, I, B, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps_mat(num) );
  
  valid = num_grad( testMetric, B, I, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = grad( testMetric, B, I, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps_mat(num) );
   
  valid = num_grad( testMetric, A, B, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid  = grad( testMetric, A, B, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps_mat(num) );
  
  valid = num_grad( testMetric, A, I, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = grad( testMetric, A, I, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps_mat(num) );
  
    // also test inverted for non-barrier metrics
  if (Barrier)
    return;
  
  valid = num_grad( testMetric, C, I, nval, num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = grad( testMetric, C, I, aval, ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( nval, aval, EPS_VAL );
  ASSERT_MATRICES_EQUAL( num, ana, eps_mat(num) );
}


TMETRIC_FUNC::compare_anaytic_and_numeric_hess()
{
  const double EPS_VAL = 1e-6;
  
  MsqError err;
  Matrix dmdA_num, dmdA_ana, d2mdA2_num[DIM*(DIM+1)/2], d2mdA2_ana[DIM*(DIM+1)/2];
  bool valid;
  double val_num, val_ana;
  
  valid = num_hess( testMetric, I, I, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = hess( testMetric, I, I, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps_mat(dmdA_num) );
  
  switch (DIM) {
    default:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], eps_mat(d2mdA2_num[3]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], eps_mat(d2mdA2_num[4]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], eps_mat(d2mdA2_num[5]) );
    case 2:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], eps_mat(d2mdA2_num[0]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], eps_mat(d2mdA2_num[1]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], eps_mat(d2mdA2_num[2]) );
  }
  
  valid = num_hess( testMetric, I, A, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = hess( testMetric, I, A, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps_mat(dmdA_num) );
  switch (DIM) {
    default:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], eps_mat(d2mdA2_num[3]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], eps_mat(d2mdA2_num[4]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], eps_mat(d2mdA2_num[5]) );
    case 2:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], eps_mat(d2mdA2_num[0]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], eps_mat(d2mdA2_num[1]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], eps_mat(d2mdA2_num[2]) );
  }
  
  valid = num_hess( testMetric, A, I, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = hess( testMetric, A, I, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps_mat(dmdA_num) );
  switch (DIM) {
    default:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], eps_mat(d2mdA2_num[3]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], eps_mat(d2mdA2_num[4]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], eps_mat(d2mdA2_num[5]) );
    case 2:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], eps_mat(d2mdA2_num[0]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], eps_mat(d2mdA2_num[1]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], eps_mat(d2mdA2_num[2]) );
  }
 
  valid = num_hess( testMetric, B, I, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = hess( testMetric, B, I, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps_mat(dmdA_num) );
  switch (DIM) {
    default:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], eps_mat(d2mdA2_num[3]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], eps_mat(d2mdA2_num[4]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], eps_mat(d2mdA2_num[5]) );
    case 2:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], eps_mat(d2mdA2_num[0]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], eps_mat(d2mdA2_num[1]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], eps_mat(d2mdA2_num[2]) );
  }
  
  valid = num_hess( testMetric, I, B, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = hess( testMetric, I, B, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps_mat(dmdA_num) );
  switch (DIM) {
    default:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], eps_mat(d2mdA2_num[3]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], eps_mat(d2mdA2_num[4]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], eps_mat(d2mdA2_num[5]) );
    case 2:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], eps_mat(d2mdA2_num[0]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], eps_mat(d2mdA2_num[1]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], eps_mat(d2mdA2_num[2]) );
  }
  
  valid = num_hess( testMetric, A, B, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = hess( testMetric, A, B, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps_mat(dmdA_num) );
  switch (DIM) {
    default:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], eps_mat(d2mdA2_num[3]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], eps_mat(d2mdA2_num[4]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], eps_mat(d2mdA2_num[5]) );
    case 2:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], eps_mat(d2mdA2_num[0]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], eps_mat(d2mdA2_num[1]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], eps_mat(d2mdA2_num[2]) );
  }
  
  valid = num_hess( testMetric, B, A, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = hess( testMetric, B, A, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  ASSERT_MATRICES_EQUAL( dmdA_num, dmdA_ana, eps_mat(dmdA_num) );
  switch (DIM) {
    default:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], eps_mat(d2mdA2_num[3]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], eps_mat(d2mdA2_num[4]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], eps_mat(d2mdA2_num[5]) );
    case 2:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], eps_mat(d2mdA2_num[0]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], eps_mat(d2mdA2_num[1]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], eps_mat(d2mdA2_num[2]) );
  }
  
    // also test inverted for non-barrier metrics
  if (Barrier)
    return;
  
  valid = num_hess( testMetric, C, I, val_num, dmdA_num, d2mdA2_num, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  valid = hess( testMetric, C, I, val_ana, dmdA_ana, d2mdA2_ana, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(valid);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val_num, val_ana, EPS_VAL );
  switch (DIM) {
    default:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[3], d2mdA2_ana[3], eps_mat(d2mdA2_num[3]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[4], d2mdA2_ana[4], eps_mat(d2mdA2_num[4]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[5], d2mdA2_ana[5], eps_mat(d2mdA2_num[5]) );
    case 2:
  ASSERT_MATRICES_EQUAL( d2mdA2_num[0], d2mdA2_ana[0], eps_mat(d2mdA2_num[0]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[1], d2mdA2_ana[1], eps_mat(d2mdA2_num[1]) );
  ASSERT_MATRICES_EQUAL( d2mdA2_num[2], d2mdA2_ana[2], eps_mat(d2mdA2_num[2]) );
  }
}

