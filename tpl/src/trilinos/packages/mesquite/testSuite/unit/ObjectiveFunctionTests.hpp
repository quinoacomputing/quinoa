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


/** \file ObjectiveFunctionTests.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_OBJECTIVE_FUNCTION_TESTS_HPP
#define MSQ_OBJECTIVE_FUNCTION_TESTS_HPP

#include "Mesquite.hpp"
#include "PatchDataInstances.hpp"
#include "Mesquite_IdealWeightInverseMeanRatio.hpp"
#include "Mesquite_ObjectiveFunctionTemplate.hpp"
#include "UnitUtil.hpp"
#include "Mesquite_MsqHessian.hpp"

#include <memory>

using namespace Mesquite;
using namespace std;

class ObjectiveFunctionTests {

public:

/** Which function to call */
enum OFTestMode { EVAL, GRAD, DIAG, HESS };

/** Test eval type support for OF templates that provide
 *  block coordinate descent functionality.  
 *  
 *  Note: If OF does not support BCD, this test will fail,
 *        even for the ObjectiveFunction::CALCULATE eval type.
 */
static void test_eval_type( ObjectiveFunction::EvalType,
                            OFTestMode test_mode,
                            ObjectiveFunctionTemplate* of );

/** Verify that if QualityMetric returns invalid but not error
 *  from evaluate, that the ObjectiveFunction does the same */
static void test_handles_invalid_qm( OFTestMode test_mode, ObjectiveFunctionTemplate* of );
                                     
/** Verify that OF correctly handles error in QM::evaluate() */
static void test_handles_qm_error( OFTestMode test_mode, ObjectiveFunctionTemplate* of );

/** Test the behavior of the clone() method */
static void test_clone( ObjectiveFunctionTemplate* of );

/** Test correct handling of QM negate flag */
static void test_negate_flag( OFTestMode test_mode, ObjectiveFunctionTemplate* of );

/** Test OF value */
static void test_value( const double* input_values,
                        unsigned num_input_values,
                        double expected_value,
                        OFTestMode test_mode,
                        ObjectiveFunctionTemplate* of );
                        
/** Compare numerical and analytical gradient values */
static inline void compare_numerical_gradient( ObjectiveFunctionTemplate* of );
static void compare_numerical_gradient( ObjectiveFunction* of );

/** Compare gradients from evaluate_with_gradient with gradient
 *  from evaluate_with_Hessian
 */
static inline void compare_hessian_gradient( ObjectiveFunctionTemplate* of );
static void compare_hessian_gradient( ObjectiveFunction* of );

/** Compare gradients from evaluate_with_gradient with gradient
 *  from evaluate_with_Hessian_diagonal
 */
static inline void compare_diagonal_gradient( ObjectiveFunctionTemplate* of );
static void compare_diagonal_gradient( ObjectiveFunction* of );

/** Compare gradient and diagonal terms from evaluate_with_Hessian
 *  and evaluate_with_Hessian_diagonal
 */
static inline void compare_hessian_diagonal( ObjectiveFunctionTemplate* of );
static void compare_hessian_diagonal( ObjectiveFunction* of );

static void compare_numerical_hessian( ObjectiveFunctionTemplate* of );
static void compare_numerical_hessian_diagonal( ObjectiveFunctionTemplate* of );
static void compare_numerical_hessian( ObjectiveFunction* of, bool diagonal_only );

static PatchData& patch();

private:

static double evaluate_internal( ObjectiveFunction::EvalType type, 
                                 OFTestMode test_mode,
                                 ObjectiveFunction* of );
};

/** The QualityMetric to use for testing purposes
 *
 *  Just pass a specified list of values to the OF 
 */
class OFTestQM : public QualityMetric {
  public:
    
    OFTestQM( ) : negateFlag(1) {}
  
    OFTestQM( const double* values, unsigned num_values )
      : mValues(num_values), negateFlag(1)
      { copy( values, values+num_values, mValues.begin() ); }
    
    void set_values( const double* values, unsigned num_values )
      { 
        mValues.resize(num_values);
        copy( values, values+num_values, mValues.begin() );
      }
      
    void append_values( const double* values, unsigned num_values )
      { copy( values, values+num_values, back_inserter(mValues) ); }
     
    virtual MetricType get_metric_type() const  { return ELEMENT_BASED; }
    
    virtual string get_name() const { return "ObjectiveFunctionTests"; }
    
    virtual int get_negate_flag() const { return negateFlag; }
    
    void set_negate_flag( int value ) { negateFlag = value; }
    
    virtual void get_evaluations( PatchData&, vector<size_t>& h, bool, MsqError& )
      {
        h.resize( mValues.size() );
        for (unsigned i = 0; i < mValues.size(); ++i)
          h[i] = i;
      }
    
    virtual bool evaluate( PatchData&, size_t h, double& v, MsqError& err )
      {
        if (h >= mValues.size()) {
          MSQ_SETERR(err)("handle out of range", MsqError::INVALID_ARG );
          return false;
        }
        v = mValues[h];
        return true;
      }
    
    virtual bool evaluate_with_indices( PatchData& pd, size_t h, double& v, vector<size_t>& i, MsqError& err )
      {
        i.clear();
        for (unsigned j = 0; j < pd.num_free_vertices(); ++j)
          i.push_back(j);
        return evaluate( pd, h, v, err );
      }
    
    virtual bool evaluate_with_gradient( PatchData& pd, size_t h, double& v, vector<size_t>& i, vector<Vector3D>& g, MsqError& err )
      {
        g.clear();
        bool rval = evaluate_with_indices( pd, h, v, i, err );
        // grad values are just used to test negate flag, so just
        // pass back an arbitrary value for each free vertex
        for (unsigned j = 0; j < i.size(); ++j)
          g.push_back( Vector3D(1,0,2) );
        return rval;
      }

    virtual bool evaluate_with_Hessian( PatchData& pd, size_t h, double& v, vector<size_t>& i, vector<Vector3D>& g, vector<Matrix3D>& H, MsqError& err )
      {
        H.clear();
        bool rval = evaluate_with_gradient( pd, h, v, i, g, err );
        // Hessian values are just used to test negate flag, so
        // pass back arbirary values.  
        for (unsigned r = 0; r < i.size(); ++r)
          for (unsigned c = r; c < i.size(); ++c)
            H.push_back( Matrix3D(1.0) );
        return rval;
      }
  
  private:
    
    vector<double> mValues;
    int negateFlag;
};

inline void ObjectiveFunctionTests::compare_numerical_gradient( ObjectiveFunctionTemplate* of )
{
  MsqPrintError err(std::cout);
  IdealWeightInverseMeanRatio metric(err);
  ASSERT_NO_ERROR( err );
  of->set_quality_metric( &metric );
  compare_numerical_gradient( (ObjectiveFunction*)of );
}

inline void ObjectiveFunctionTests::compare_hessian_gradient( ObjectiveFunctionTemplate* of )
{
  MsqPrintError err(std::cout);
  IdealWeightInverseMeanRatio metric(err);
  ASSERT_NO_ERROR( err );
  of->set_quality_metric( &metric );
  compare_hessian_gradient( (ObjectiveFunction*)of );
}

inline void ObjectiveFunctionTests::compare_diagonal_gradient( ObjectiveFunctionTemplate* of )
{
  MsqPrintError err(std::cout);
  IdealWeightInverseMeanRatio metric(err);
  ASSERT_NO_ERROR( err );
  of->set_quality_metric( &metric );
  compare_diagonal_gradient( (ObjectiveFunction*)of );
}

inline void ObjectiveFunctionTests::compare_hessian_diagonal( ObjectiveFunctionTemplate* of )
{
  MsqPrintError err(std::cout);
  IdealWeightInverseMeanRatio metric(err);
  ASSERT_NO_ERROR( err );
  of->set_quality_metric( &metric );
  compare_hessian_diagonal( (ObjectiveFunction*)of );
}

inline void ObjectiveFunctionTests::compare_numerical_hessian( ObjectiveFunctionTemplate* of )
{
  MsqPrintError err(std::cout);
  IdealWeightInverseMeanRatio metric(err);
  ASSERT_NO_ERROR( err );
  of->set_quality_metric( &metric );
  compare_numerical_hessian( (ObjectiveFunction*)of, false );
}

inline void ObjectiveFunctionTests::compare_numerical_hessian_diagonal( ObjectiveFunctionTemplate* of )
{
  MsqPrintError err(std::cout);
  IdealWeightInverseMeanRatio metric(err);
  ASSERT_NO_ERROR( err );
  of->set_quality_metric( &metric );
  compare_numerical_hessian( (ObjectiveFunction*)of, true );
}

#endif
