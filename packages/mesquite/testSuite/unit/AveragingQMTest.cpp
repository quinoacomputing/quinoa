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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */

/*! \file AveragingQMTest.cpp

Unit testing for the AveragingQM class
\author Jasno Kraftcheck
*/
#include "Mesquite.hpp"
#include "Mesquite_AveragingQM.hpp"
#include "Mesquite_IdealElements.hpp"
#include "UnitUtil.hpp"
#include "Mesquite_TopologyInfo.hpp"
#include "Mesquite_PatchData.hpp"
#include "Mesquite_cppunit/extensions/HelperMacros.h"

using namespace Mesquite;

class AveragingQMTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(AveragingQMTest);
  CPPUNIT_TEST (test_average_metrics_liner);
  CPPUNIT_TEST (test_average_metrics_rms);
  CPPUNIT_TEST (test_average_metrics_hms);
  CPPUNIT_TEST (test_average_metrics_minimum);
  CPPUNIT_TEST (test_average_metrics_maximum);
  CPPUNIT_TEST (test_average_metrics_harmonic);
  CPPUNIT_TEST (test_average_metrics_geometric);
  CPPUNIT_TEST (test_average_metrics_sum);
  CPPUNIT_TEST (test_average_metrics_sum_squared);
  CPPUNIT_TEST (test_average_metrics_standard_deviation);
  CPPUNIT_TEST (test_average_metrics_max_over_min);
  CPPUNIT_TEST (test_average_metrics_max_minus_min);
  CPPUNIT_TEST (test_average_metrics_sum_of_ratios_squared);
  CPPUNIT_TEST (test_average_and_weights_linear);
  CPPUNIT_TEST (test_average_and_weights_rms);
  CPPUNIT_TEST (test_average_and_weights_hms);
  CPPUNIT_TEST (test_average_and_weights_minimum);
  CPPUNIT_TEST (test_average_and_weights_maximum);
  CPPUNIT_TEST (test_average_and_weights_harmonic);
  CPPUNIT_TEST (test_average_and_weights_geometric);
  CPPUNIT_TEST (test_average_and_weights_sum);
  CPPUNIT_TEST (test_average_and_weights_sum_squared);
  CPPUNIT_TEST (test_average_and_weights_standard_deviation);
  CPPUNIT_TEST (test_average_and_weights_max_over_min);
  CPPUNIT_TEST (test_average_and_weights_max_minus_min);
  CPPUNIT_TEST (test_average_and_weights_sum_of_ratios_squared);
  CPPUNIT_TEST (test_average_corner_gradients_linear);
  CPPUNIT_TEST (test_average_corner_gradients_rms);
  CPPUNIT_TEST (test_average_corner_gradients_hms);
  CPPUNIT_TEST (test_average_corner_gradients_minimum);
  CPPUNIT_TEST (test_average_corner_gradients_maximum);
  CPPUNIT_TEST (test_average_corner_gradients_harmonic);
  CPPUNIT_TEST (test_average_corner_gradients_geometric);
  CPPUNIT_TEST (test_average_corner_gradients_sum);
  CPPUNIT_TEST (test_average_corner_gradients_sum_squared);
  CPPUNIT_TEST (test_average_corner_gradients_standard_deviation);
  CPPUNIT_TEST (test_average_corner_gradients_max_over_min);
  CPPUNIT_TEST (test_average_corner_gradients_max_minus_min);
  CPPUNIT_TEST (test_average_corner_gradients_sum_of_ratios_squared);
  CPPUNIT_TEST (test_average_corner_hessians_linear);
  CPPUNIT_TEST (test_average_corner_hessians_rms);
  CPPUNIT_TEST (test_average_corner_hessians_hms);
  CPPUNIT_TEST (test_average_corner_hessians_minimum);
  CPPUNIT_TEST (test_average_corner_hessians_maximum);
  CPPUNIT_TEST (test_average_corner_hessians_harmonic);
  CPPUNIT_TEST (test_average_corner_hessians_geometric);
  CPPUNIT_TEST (test_average_corner_hessians_sum);
  CPPUNIT_TEST (test_average_corner_hessians_sum_squared);
  CPPUNIT_TEST (test_average_corner_hessians_standard_deviation);
  CPPUNIT_TEST (test_average_corner_hessians_max_over_min);
  CPPUNIT_TEST (test_average_corner_hessians_max_minus_min);
  CPPUNIT_TEST (test_average_corner_hessians_sum_of_ratios_squared);
  CPPUNIT_TEST (test_average_corner_hessian_diagonals_linear);
  CPPUNIT_TEST (test_average_corner_hessian_diagonals_rms);
  CPPUNIT_TEST (test_average_corner_hessian_diagonals_hms);
  CPPUNIT_TEST (test_average_corner_hessian_diagonals_minimum);
  CPPUNIT_TEST (test_average_corner_hessian_diagonals_maximum);
  CPPUNIT_TEST (test_average_corner_hessian_diagonals_harmonic);
  CPPUNIT_TEST (test_average_corner_hessian_diagonals_geometric);
  CPPUNIT_TEST (test_average_corner_hessian_diagonals_sum);
  CPPUNIT_TEST (test_average_corner_hessian_diagonals_sum_squared);
  CPPUNIT_TEST (test_average_corner_hessian_diagonals_standard_deviation);
  CPPUNIT_TEST (test_average_corner_hessian_diagonals_max_over_min);
  CPPUNIT_TEST (test_average_corner_hessian_diagonals_max_minus_min);
  CPPUNIT_TEST (test_average_corner_hessian_diagonals_sum_of_ratios_squared);
  CPPUNIT_TEST_SUITE_END();
  
  static const double VAL_LIST_1[5];
  static const double VAL_LIST_2[8];
  static const unsigned LEN_LIST_1;
  static const unsigned LEN_LIST_2;
  
  void check_average_and_weights_fails( QualityMetric::AveragingMethod scheme );
  
  void check_average_gradients( QualityMetric::AveragingMethod scheme );
  
  void check_average_gradient_fails( QualityMetric::AveragingMethod scheme );
  
  void check_pmean_hessian_diagonals( QualityMetric::AveragingMethod scheme,
                            double inner_power,
                            double outer_power,
                            bool scale );
                                 
  void check_hessian_diagonal_fails( QualityMetric::AveragingMethod scheme );
  
  void check_pmean_hessian( QualityMetric::AveragingMethod scheme,
                            double inner_power,
                            double outer_power,
                            bool scale );
                                 
  void check_hessian_fails( QualityMetric::AveragingMethod scheme );
  
  void check_average_and_weights( const double* vals, unsigned n,
                                  QualityMetric::AveragingMethod method,
                                  const double* weights );
public:
  
  void test_average_metrics_liner();
  void test_average_metrics_rms();
  void test_average_metrics_hms();
  void test_average_metrics_minimum();
  void test_average_metrics_maximum();
  void test_average_metrics_harmonic();
  void test_average_metrics_geometric();
  void test_average_metrics_sum();
  void test_average_metrics_sum_squared();
  void test_average_metrics_standard_deviation();
  void test_average_metrics_max_over_min();
  void test_average_metrics_max_minus_min();
  void test_average_metrics_sum_of_ratios_squared();
  void test_average_and_weights_linear();
  void test_average_and_weights_rms();
  void test_average_and_weights_hms();
  void test_average_and_weights_minimum();
  void test_average_and_weights_maximum();
  void test_average_and_weights_harmonic();
  void test_average_and_weights_geometric();
  void test_average_and_weights_sum();
  void test_average_and_weights_sum_squared();
  void test_average_and_weights_standard_deviation();
  void test_average_and_weights_max_over_min();
  void test_average_and_weights_max_minus_min();
  void test_average_and_weights_sum_of_ratios_squared();
  void test_average_corner_gradients_linear();
  void test_average_corner_gradients_rms();
  void test_average_corner_gradients_hms();
  void test_average_corner_gradients_minimum();
  void test_average_corner_gradients_maximum();
  void test_average_corner_gradients_harmonic();
  void test_average_corner_gradients_geometric();
  void test_average_corner_gradients_sum();
  void test_average_corner_gradients_sum_squared();
  void test_average_corner_gradients_standard_deviation();
  void test_average_corner_gradients_max_over_min();
  void test_average_corner_gradients_max_minus_min();
  void test_average_corner_gradients_sum_of_ratios_squared();
  void test_average_corner_hessians_linear();
  void test_average_corner_hessians_rms();
  void test_average_corner_hessians_hms();
  void test_average_corner_hessians_minimum();
  void test_average_corner_hessians_maximum();
  void test_average_corner_hessians_harmonic();
  void test_average_corner_hessians_geometric();
  void test_average_corner_hessians_sum();
  void test_average_corner_hessians_sum_squared();
  void test_average_corner_hessians_standard_deviation();
  void test_average_corner_hessians_max_over_min();
  void test_average_corner_hessians_max_minus_min();
  void test_average_corner_hessians_sum_of_ratios_squared();
  void test_average_corner_hessian_diagonals_linear();
  void test_average_corner_hessian_diagonals_rms();
  void test_average_corner_hessian_diagonals_hms();
  void test_average_corner_hessian_diagonals_minimum();
  void test_average_corner_hessian_diagonals_maximum();
  void test_average_corner_hessian_diagonals_harmonic();
  void test_average_corner_hessian_diagonals_geometric();
  void test_average_corner_hessian_diagonals_sum();
  void test_average_corner_hessian_diagonals_sum_squared();
  void test_average_corner_hessian_diagonals_standard_deviation();
  void test_average_corner_hessian_diagonals_max_over_min();
  void test_average_corner_hessian_diagonals_max_minus_min();
  void test_average_corner_hessian_diagonals_sum_of_ratios_squared();
 
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(AveragingQMTest, "AveragingQMTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(AveragingQMTest, "Unit");

const double AveragingQMTest::VAL_LIST_1[5] = { 1, 2, -1, -2, 5 };
const double AveragingQMTest::VAL_LIST_2[8] = { M_PI,
                                                std::exp(1.0),
                                                -20,
                                                8,
                                                M_PI/4,
                                                std::log(2.0),
                                                std::sqrt(2.0),
                                                -1 };
const unsigned AveragingQMTest::LEN_LIST_1 = sizeof(VAL_LIST_1)/sizeof(double);
const unsigned AveragingQMTest::LEN_LIST_2 = sizeof(VAL_LIST_2)/sizeof(double);
                                       
static double pmean( const double* vals, unsigned n, double power )
{
  double result = 0; 
  for(unsigned i = 0; i < n; ++i)
    result += std::pow( vals[i], power );
  return pow( result / n, 1.0/power );
}

static double list_min( const double* vals, unsigned n )
{
  double result = vals[0];
  for (unsigned i = 1; i < n; ++i)
    if (vals[i] < result)
      result = vals[i];
  return result;
}

static double list_max( const double* vals, unsigned n )
{
  double result = vals[0];
  for (unsigned i = 1; i < n; ++i)
    if (vals[i] > result)
      result = vals[i];
  return result;
}

static double sum_sqr( const double* vals, unsigned n )
{
  double result = 0.0;
  for (unsigned i = 0; i < n; ++i)
    result += vals[i]*vals[i];
  return result;
}

static double sum_of_ratios_squared( const double* vals, unsigned n )
{
  double result = 0.0;
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < n; ++j) {
      double ratio = vals[i]/vals[j];
      result += ratio*ratio;
    }
  return result / (n*n);
}

static double geometric_mean( const double* vals, unsigned n )
{
  double result = 1.0;
  for (unsigned i = 0; i < n; ++i)
    result *= vals[i];
  return std::pow( result, 1.0/n );
}

void AveragingQMTest::test_average_metrics_liner()
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( QualityMetric::LINEAR );
  double exp, act;
  
  exp = pmean(VAL_LIST_1, LEN_LIST_1, 1);
  act = aqm.average_metrics( VAL_LIST_1, LEN_LIST_1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
  
  exp = pmean(VAL_LIST_2, LEN_LIST_2, 1);
  act = aqm.average_metrics( VAL_LIST_2, LEN_LIST_2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
}

void AveragingQMTest::test_average_metrics_rms()
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( QualityMetric::RMS );
  double exp, act;
  
  exp = pmean(VAL_LIST_1, LEN_LIST_1, 2);
  act = aqm.average_metrics( VAL_LIST_1, LEN_LIST_1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
  
  exp = pmean(VAL_LIST_2, LEN_LIST_2, 2);
  act = aqm.average_metrics( VAL_LIST_2, LEN_LIST_2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
}

void AveragingQMTest::test_average_metrics_hms()
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( QualityMetric::HMS );
  double exp, act;
  
  exp = pmean(VAL_LIST_1, LEN_LIST_1, -2.0);
  act = aqm.average_metrics( VAL_LIST_1, LEN_LIST_1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
  
  exp = pmean(VAL_LIST_2, LEN_LIST_2, -2.0);
  act = aqm.average_metrics( VAL_LIST_2, LEN_LIST_2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
}

void AveragingQMTest::test_average_metrics_minimum()
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( QualityMetric::MINIMUM );
  double exp, act;
  
  exp = list_min(VAL_LIST_1, LEN_LIST_1);
  act = aqm.average_metrics( VAL_LIST_1, LEN_LIST_1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
  
  exp = list_min(VAL_LIST_2, LEN_LIST_2);
  act = aqm.average_metrics( VAL_LIST_2, LEN_LIST_2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
}

void AveragingQMTest::test_average_metrics_maximum()
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( QualityMetric::MAXIMUM );
  double exp, act;
  
  exp = list_max(VAL_LIST_1, LEN_LIST_1);
  act = aqm.average_metrics( VAL_LIST_1, LEN_LIST_1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
  
  exp = list_max(VAL_LIST_2, LEN_LIST_2);
  act = aqm.average_metrics( VAL_LIST_2, LEN_LIST_2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
}

void AveragingQMTest::test_average_metrics_harmonic()
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( QualityMetric::HARMONIC );
  double exp, act;
  
  exp = pmean(VAL_LIST_1, LEN_LIST_1, -1.0);
  act = aqm.average_metrics( VAL_LIST_1, LEN_LIST_1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
  
  exp = pmean(VAL_LIST_2, LEN_LIST_2, -1.0);
  act = aqm.average_metrics( VAL_LIST_2, LEN_LIST_2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
}

void AveragingQMTest::test_average_metrics_geometric()
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( QualityMetric::GEOMETRIC );
  double exp, act;
  
  exp = geometric_mean( VAL_LIST_1, LEN_LIST_1 );
  act = aqm.average_metrics( VAL_LIST_1, LEN_LIST_1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
  
  exp = geometric_mean( VAL_LIST_2, LEN_LIST_2 );
  act = aqm.average_metrics( VAL_LIST_2, LEN_LIST_2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
}

void AveragingQMTest::test_average_metrics_sum()
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( QualityMetric::SUM );
  double exp, act;
  unsigned i;
  
  for (exp = 0.0, i = 0; i < LEN_LIST_1; exp += VAL_LIST_1[i++]);
  act = aqm.average_metrics( VAL_LIST_1, LEN_LIST_1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
  
  for (exp = 0.0, i = 0; i < LEN_LIST_2; exp += VAL_LIST_2[i++]);
  act = aqm.average_metrics( VAL_LIST_2, LEN_LIST_2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
}

void AveragingQMTest::test_average_metrics_sum_squared()
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( QualityMetric::SUM_SQUARED );
  double exp, act;
  
  exp = sum_sqr( VAL_LIST_1, LEN_LIST_1 );
  act = aqm.average_metrics( VAL_LIST_1, LEN_LIST_1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
  
  exp = sum_sqr( VAL_LIST_2, LEN_LIST_2 );
  act = aqm.average_metrics( VAL_LIST_2, LEN_LIST_2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
}

void AveragingQMTest::test_average_metrics_standard_deviation()
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( QualityMetric::STANDARD_DEVIATION );
  double exp, act, tmp;
  
  tmp = pmean( VAL_LIST_1, LEN_LIST_1, 1.0 );
  exp = sum_sqr( VAL_LIST_1, LEN_LIST_1 );
  exp = exp/LEN_LIST_1 - tmp*tmp;
  act = aqm.average_metrics( VAL_LIST_1, LEN_LIST_1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
  
  tmp = pmean( VAL_LIST_2, LEN_LIST_2, 1.0 );
  exp = sum_sqr( VAL_LIST_2, LEN_LIST_2 );
  exp = exp/LEN_LIST_2 - tmp*tmp;
  act = aqm.average_metrics( VAL_LIST_2, LEN_LIST_2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
}

void AveragingQMTest::test_average_metrics_max_over_min()
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( QualityMetric::MAX_OVER_MIN );
  double exp, act;
  
  exp = list_max( VAL_LIST_1, LEN_LIST_1 )/list_min( VAL_LIST_1, LEN_LIST_1 );
  act = aqm.average_metrics( VAL_LIST_1, LEN_LIST_1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
  
  exp = list_max( VAL_LIST_2, LEN_LIST_2 )/list_min( VAL_LIST_2, LEN_LIST_2 );
  act = aqm.average_metrics( VAL_LIST_2, LEN_LIST_2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
}

void AveragingQMTest::test_average_metrics_max_minus_min()
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( QualityMetric::MAX_MINUS_MIN );
  double exp, act;
  
  exp = list_max( VAL_LIST_1, LEN_LIST_1 )-list_min( VAL_LIST_1, LEN_LIST_1 );
  act = aqm.average_metrics( VAL_LIST_1, LEN_LIST_1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
  
  exp = list_max( VAL_LIST_2, LEN_LIST_2 )-list_min( VAL_LIST_2, LEN_LIST_2 );
  act = aqm.average_metrics( VAL_LIST_2, LEN_LIST_2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
}

void AveragingQMTest::test_average_metrics_sum_of_ratios_squared()
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( QualityMetric::SUM_OF_RATIOS_SQUARED );
  double exp, act;
  
  exp = sum_of_ratios_squared( VAL_LIST_1, LEN_LIST_1 );
  act = aqm.average_metrics( VAL_LIST_1, LEN_LIST_1, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
  
  exp = sum_of_ratios_squared( VAL_LIST_2, LEN_LIST_2 );
  act = aqm.average_metrics( VAL_LIST_2, LEN_LIST_2, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp, act, 1e-8 ); 
}

void AveragingQMTest::check_average_and_weights( const double* vals, unsigned n,
                                       QualityMetric::AveragingMethod method,
                                       const double* weights )
{
  MsqPrintError err(std::cout);
  AveragingQM aqm( method );
  
  CPPUNIT_ASSERT(n > 0);
  std::vector<double> working(n);
  std::copy( vals, vals+n, working.begin() );
  
  double avg1, avg2;
  avg1 = aqm.average_metrics( vals, n, err );
  ASSERT_NO_ERROR(err);
  avg2 = aqm.average_metric_and_weights( arrptr(working), n, err );
  ASSERT_NO_ERROR(err);
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( avg1, avg2, 1e-6 );
  for (unsigned i = 0; i < n; ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL( weights[i], working[i], 1e-6 );
}

void AveragingQMTest::test_average_and_weights_linear()
{
  std::vector<double> weights1( LEN_LIST_1, 1.0/LEN_LIST_1), weights2( LEN_LIST_2, 1.0/LEN_LIST_2 );
  check_average_and_weights( VAL_LIST_1, LEN_LIST_1, QualityMetric::LINEAR, arrptr(weights1) );
  check_average_and_weights( VAL_LIST_2, LEN_LIST_2, QualityMetric::LINEAR, arrptr(weights2) );
}

void AveragingQMTest::test_average_and_weights_rms()
{
  std::vector<double> weights1( LEN_LIST_1 ), weights2( LEN_LIST_2 );
  unsigned i;
  double rms;
  
  rms = pmean( VAL_LIST_1, LEN_LIST_1, 2.0 );
  for (i = 0; i < LEN_LIST_1; ++i)
    weights1[i] = VAL_LIST_1[i]/(LEN_LIST_1*rms);
  check_average_and_weights( VAL_LIST_1, LEN_LIST_1, QualityMetric::RMS, arrptr(weights1) );
  
  rms = pmean( VAL_LIST_2, LEN_LIST_2, 2.0 );
  for (i = 0; i < LEN_LIST_2; ++i)
    weights2[i] = VAL_LIST_2[i]/(LEN_LIST_2*rms);
  check_average_and_weights( VAL_LIST_2, LEN_LIST_2, QualityMetric::RMS, arrptr(weights2) );
}

void AveragingQMTest::test_average_and_weights_hms()
{
  std::vector<double> weights1( LEN_LIST_1 ), weights2( LEN_LIST_2 );
  unsigned i;
  double hms, tmp;
  
  hms = pmean( VAL_LIST_1, LEN_LIST_1, -2.0 );
  for (i = 0; i < LEN_LIST_1; ++i) {
    tmp = hms / VAL_LIST_1[i];
    weights1[i] = tmp*tmp*tmp / LEN_LIST_1;
  }
  check_average_and_weights( VAL_LIST_1, LEN_LIST_1, QualityMetric::HMS, arrptr(weights1) );
  
  hms = pmean( VAL_LIST_2, LEN_LIST_2, -2.0 );
  for (i = 0; i < LEN_LIST_2; ++i) {
    tmp = hms / VAL_LIST_2[i];
    weights2[i] = tmp*tmp*tmp / LEN_LIST_2;
  }
  check_average_and_weights( VAL_LIST_2, LEN_LIST_2, QualityMetric::HMS, arrptr(weights2) );
}

void AveragingQMTest::test_average_and_weights_minimum()
{
  std::vector<double> weights1( LEN_LIST_1, 0.0 ), weights2( LEN_LIST_2, 0.0 );
  const double* ptr;
  
  ptr = std::min_element( VAL_LIST_1, VAL_LIST_1+LEN_LIST_1 );
  weights1[ptr-VAL_LIST_1] = 1.0;
  check_average_and_weights( VAL_LIST_1, LEN_LIST_1, QualityMetric::MINIMUM, arrptr(weights1) );
  
  ptr = std::min_element( VAL_LIST_2, VAL_LIST_2+LEN_LIST_2 );
  weights2[ptr-VAL_LIST_2] = 1.0;
  check_average_and_weights( VAL_LIST_2, LEN_LIST_2, QualityMetric::MINIMUM, arrptr(weights2) );
}  

void AveragingQMTest::test_average_and_weights_maximum()
{
  std::vector<double> weights1( LEN_LIST_1, 0.0 ), weights2( LEN_LIST_2, 0.0 );
  const double* ptr;
  
  ptr = std::max_element( VAL_LIST_1, VAL_LIST_1+LEN_LIST_1 );
  weights1[ptr-VAL_LIST_1] = 1.0;
  check_average_and_weights( VAL_LIST_1, LEN_LIST_1, QualityMetric::MAXIMUM, arrptr(weights1) );
  
  ptr = std::max_element( VAL_LIST_2, VAL_LIST_2+LEN_LIST_2 );
  weights2[ptr-VAL_LIST_2] = 1.0;
  check_average_and_weights( VAL_LIST_2, LEN_LIST_2, QualityMetric::MAXIMUM, arrptr(weights2) );
}  

void AveragingQMTest::test_average_and_weights_harmonic()
{
  std::vector<double> weights1( LEN_LIST_1 ), weights2( LEN_LIST_2 );
  unsigned i;
  double h, tmp;
  
  h = pmean( VAL_LIST_1, LEN_LIST_1, -1.0 );
  for (i = 0; i < LEN_LIST_1; ++i) {
    tmp = h / VAL_LIST_1[i];
    weights1[i] = tmp*tmp / LEN_LIST_1;
  }
  check_average_and_weights( VAL_LIST_1, LEN_LIST_1, QualityMetric::HARMONIC, arrptr(weights1) );
  
  h = pmean( VAL_LIST_2, LEN_LIST_2, -1.0 );
  for (i = 0; i < LEN_LIST_2; ++i) {
    tmp = h / VAL_LIST_2[i];
    weights2[i] = tmp*tmp / LEN_LIST_2;
  }
  check_average_and_weights( VAL_LIST_2, LEN_LIST_2, QualityMetric::HARMONIC, arrptr(weights2) );
}

void AveragingQMTest::test_average_and_weights_geometric()
{
  std::vector<double> weights1( LEN_LIST_1 ), weights2( LEN_LIST_2 );
  unsigned i;
  double g;

  g = geometric_mean( VAL_LIST_1, LEN_LIST_1 );
  for (i = 0; i < LEN_LIST_1; ++i) {
    weights1[i] = g / VAL_LIST_1[i] / LEN_LIST_1;
  }
  check_average_and_weights( VAL_LIST_1, LEN_LIST_1, QualityMetric::GEOMETRIC, arrptr(weights1) );

  g = geometric_mean( VAL_LIST_2, LEN_LIST_2 );
  for (i = 0; i < LEN_LIST_2; ++i) {
    weights2[i] = g / VAL_LIST_2[i] / LEN_LIST_2;
  }
  check_average_and_weights( VAL_LIST_2, LEN_LIST_2, QualityMetric::GEOMETRIC, arrptr(weights2) );
}

void AveragingQMTest::test_average_and_weights_sum()
{
  std::vector<double> weights1( LEN_LIST_1, 1.0), weights2( LEN_LIST_2, 1.0);
  check_average_and_weights( VAL_LIST_1, LEN_LIST_1, QualityMetric::SUM, arrptr(weights1) );
  check_average_and_weights( VAL_LIST_2, LEN_LIST_2, QualityMetric::SUM, arrptr(weights2) );
}

void AveragingQMTest::test_average_and_weights_sum_squared()
{
  std::vector<double> weights1( LEN_LIST_1 ), weights2( LEN_LIST_2 );
  unsigned i;
  
  for (i = 0; i < LEN_LIST_1; ++i) 
    weights1[i] = 2 * VAL_LIST_1[i];
  check_average_and_weights( VAL_LIST_1, LEN_LIST_1, QualityMetric::SUM_SQUARED, arrptr(weights1) );
  
  for (i = 0; i < LEN_LIST_2; ++i) 
    weights2[i] = 2 * VAL_LIST_2[i];
  check_average_and_weights( VAL_LIST_2, LEN_LIST_2, QualityMetric::SUM_SQUARED, arrptr(weights2) );
}

void AveragingQMTest::check_average_and_weights_fails( QualityMetric::AveragingMethod scheme )
{
  MsqError err;
  AveragingQM aqm(scheme);
  double w[] = {1.0, 2.0};
  aqm.average_metric_and_weights( w, 2, err );
  CPPUNIT_ASSERT(err);
}

void AveragingQMTest::test_average_and_weights_standard_deviation()
  { check_average_and_weights_fails( QualityMetric::STANDARD_DEVIATION ); }

void AveragingQMTest::test_average_and_weights_max_over_min()
  { check_average_and_weights_fails( QualityMetric::MAX_OVER_MIN ); }

void AveragingQMTest::test_average_and_weights_max_minus_min()
  { check_average_and_weights_fails( QualityMetric::MAX_MINUS_MIN ); }

void AveragingQMTest::test_average_and_weights_sum_of_ratios_squared()
  { check_average_and_weights_fails( QualityMetric::SUM_OF_RATIOS_SQUARED ); }

void AveragingQMTest::check_average_gradients( QualityMetric::AveragingMethod scheme )
{
  double vals1[] = { 1, -2, -1, 0.5 };
  double vals2[] = { 1, -2, -1, 0.5 };
  Vector3D grads[12], vtxgrads[4];
  for (unsigned i = 0; i < 12; ++i)
    grads[i] = Vector3D( (double)i );

  uint32_t fixed = 0x2; // mark second vertex as fixed
  Vector3D exp;
  MsqError err;
  AveragingQM aqm(scheme);

  double a1 = aqm.average_metric_and_weights( vals1, 4, err );
  ASSERT_NO_ERROR(err);

  double a2 = aqm.average_corner_gradients( QUADRILATERAL,
                                            fixed, 4,
                                            vals2,
                                            grads, vtxgrads,
                                            err );
  ASSERT_NO_ERROR(err);

  CPPUNIT_ASSERT_DOUBLES_EQUAL( a1, a2, 1e-8 );

  exp = vals1[0]*grads[0] + vals1[1]*grads[5] + vals1[3]*grads[10];
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp, vtxgrads[0], 1e-6 );

  exp = vals1[1]*grads[4] + vals1[2]*grads[6] + vals1[3]*grads[11];
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp, vtxgrads[2], 1e-6 );

  exp = vals1[0]*grads[2] + vals1[2]*grads[7] + vals1[3]*grads[9];
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp, vtxgrads[3], 1e-6 );
}


void AveragingQMTest::test_average_corner_gradients_linear()
  { check_average_gradients( QualityMetric::LINEAR ); }

void AveragingQMTest::test_average_corner_gradients_rms()
  { check_average_gradients( QualityMetric::RMS ); }

void AveragingQMTest::test_average_corner_gradients_hms()
  { check_average_gradients( QualityMetric::HMS ); }

void AveragingQMTest::test_average_corner_gradients_minimum()
  { check_average_gradients( QualityMetric::MINIMUM ); }

void AveragingQMTest::test_average_corner_gradients_maximum()
  { check_average_gradients( QualityMetric::MAXIMUM ); }

void AveragingQMTest::test_average_corner_gradients_harmonic()
  { check_average_gradients( QualityMetric::HARMONIC ); }

void AveragingQMTest::test_average_corner_gradients_geometric()
  { check_average_gradients( QualityMetric::GEOMETRIC ); }

void AveragingQMTest::test_average_corner_gradients_sum()
  { check_average_gradients( QualityMetric::SUM ); }

void AveragingQMTest::test_average_corner_gradients_sum_squared()
  { check_average_gradients( QualityMetric::SUM_SQUARED ); }

void AveragingQMTest::check_average_gradient_fails( QualityMetric::AveragingMethod scheme )
{
  double vals[] = { 1, 2, -1, 0.5 };
  Vector3D grads[12], vtxgrads[4];
  for (unsigned i = 0; i < 12; ++i)
    grads[i] = Vector3D( (double)i );
  
  MsqError err;
  AveragingQM aqm(scheme);

  aqm.average_corner_gradients( QUADRILATERAL, 0, 4, vals, grads, vtxgrads, err );
  CPPUNIT_ASSERT(err);
}
  
void AveragingQMTest::test_average_corner_gradients_standard_deviation()
  { check_average_gradient_fails( QualityMetric::STANDARD_DEVIATION ); }

void AveragingQMTest::test_average_corner_gradients_max_over_min()
  { check_average_gradient_fails( QualityMetric::MAX_OVER_MIN ); }

void AveragingQMTest::test_average_corner_gradients_max_minus_min()
  { check_average_gradient_fails( QualityMetric::MAX_MINUS_MIN ); }

void AveragingQMTest::test_average_corner_gradients_sum_of_ratios_squared()
  { check_average_gradient_fails( QualityMetric::SUM_OF_RATIOS_SQUARED ); }



static double pmean_of_triangle_corner_hessians( double inner_power,
                                                 double outer_power,
                                                 const double* v,
                                                 const Vector3D* cg,
                                                 const Matrix3D* ch,
                                                 Vector3D* tg,
                                                 Matrix3D* th,
                                                 bool scale )
{
  Matrix3D op;
  double gf[3], hf[3];
  double nm, m = 0;
  double den = scale ? 3.0: 1.0;
  for (unsigned i = 0; i < 3; ++i) {
	  nm = pow(v[i], inner_power);
	  m += nm;

	  gf[i] = inner_power*nm / v[i] / den;
	  hf[i] = (inner_power-1)*gf[i] / v[i];
  }

  nm = m / den;

  tg[0] = gf[0] * cg[0] + gf[1] * cg[5] + gf[2] * cg[7];
  tg[1] = gf[0] * cg[1] + gf[1] * cg[3] + gf[2] * cg[8];
  tg[2] = gf[0] * cg[2] + gf[1] * cg[4] + gf[2] * cg[6];

  th[0] = hf[0]*op.outer_product( cg[0], cg[0] ) + gf[0]*ch[ 0]
        + hf[1]*op.outer_product( cg[5], cg[5] ) + gf[1]*ch[11]
        + hf[2]*op.outer_product( cg[7], cg[7] ) + gf[2]*ch[15];
  th[3] = hf[0]*op.outer_product( cg[1], cg[1] ) + gf[0]*ch[ 3]
        + hf[1]*op.outer_product( cg[3], cg[3] ) + gf[1]*ch[ 6]
        + hf[2]*op.outer_product( cg[8], cg[8] ) + gf[2]*ch[17];
  th[5] = hf[0]*op.outer_product( cg[2], cg[2] ) + gf[0]*ch[ 5]
        + hf[1]*op.outer_product( cg[4], cg[4] ) + gf[1]*ch[ 9]
        + hf[2]*op.outer_product( cg[6], cg[6] ) + gf[2]*ch[12];
  th[1] = hf[0]*op.outer_product( cg[0], cg[1] ) + gf[0]*ch[ 1]
        + hf[1]*op.outer_product( cg[5], cg[3] ) + gf[1]*transpose(ch[ 8])
        + hf[2]*op.outer_product( cg[7], cg[8] ) + gf[2]*ch[16];
  th[2] = hf[0]*op.outer_product( cg[0], cg[2] ) + gf[0]*ch[ 2]
        + hf[1]*op.outer_product( cg[5], cg[4] ) + gf[1]*transpose(ch[10])
        + hf[2]*op.outer_product( cg[7], cg[6] ) + gf[2]*transpose(ch[13]);
  th[4] = hf[0]*op.outer_product( cg[1], cg[2] ) + gf[0]*ch[ 4]
        + hf[1]*op.outer_product( cg[3], cg[4] ) + gf[1]*ch[ 7]
        + hf[2]*op.outer_product( cg[8], cg[6] ) + gf[2]*transpose(ch[14]);
 

  m = pow(nm, outer_power);
  double g = m * outer_power / nm;
  double h = (outer_power - 1.0) * g / nm;
  for (unsigned r = 0; r < 3; ++r) {
    for (unsigned c = r; c < 3; ++c) {
      *th = g * *th + h * op.outer_product( tg[r], tg[c] );
      ++th;
    }
    tg[r] *= g;
  }

  return m;
}

void AveragingQMTest::check_pmean_hessian_diagonals( 
                                 QualityMetric::AveragingMethod scheme,
                                 double inner_power,
                                 double outer_power,
                                 bool scale )
{
  MsqPrintError err(std::cout);
  AveragingQM aqm(scheme);
  uint32_t fixed = 0; // mark no vertices as fixed

    // define corner values, gradients and Hessians for a triangle
  double vals[] = { 1, 2, 0.5 };
  const Vector3D grads[9] = { Vector3D(0,1,1),
                              Vector3D(0,2,2),
                              Vector3D(0,3,3),
                              Vector3D(0,4,4),
                              Vector3D(0,5,5),
                              Vector3D(0,6,6),
                              Vector3D(0,7,7),
                              Vector3D(0,8,8),
                              Vector3D(0,9,9) };
  Matrix3D hesss[18] = {
    Matrix3D(1.0), Matrix3D(2.0), Matrix3D(-1.),
                   Matrix3D(0.5), Matrix3D(6.0),
                                  Matrix3D(1.1),
    Matrix3D(1.2), Matrix3D(0.3), Matrix3D(0.0),
                   Matrix3D(22.), Matrix3D(-2.),
                                  Matrix3D(4.3),
    Matrix3D(0.8), Matrix3D(1,7), Matrix3D(-.2),
                   Matrix3D(0.5), Matrix3D(7.0),
                                  Matrix3D(1.5)
  };
    // change some values so we catch transpositional errors
  for (unsigned i = 0; i < sizeof(hesss)/sizeof(hesss[0]); ++i)
    hesss[i][2][2] = 0.0;

    // Calculate expected values
  Vector3D aqm_grad[3], exp_grad[3];
  Matrix3D exp_hess[6];
  double aqm_avg, exp_avg;
  exp_avg = pmean_of_triangle_corner_hessians( inner_power, 
                                               outer_power,
                                               vals,
                                               grads, hesss,
                                               exp_grad, exp_hess,
                                               scale );
                                               
    // Test calculation of diagonal from full corner hessian data.
  SymMatrix3D aqm_hess[3] = { 1, 2, 3 };
  aqm_avg = aqm.average_corner_hessian_diagonals( TRIANGLE, fixed, 3,
                                         vals, grads, hesss,
                                         aqm_grad, aqm_hess, 
                                         err );
  ASSERT_NO_ERROR(err);

  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_avg, aqm_avg, 1e-6 );
  
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[0], aqm_grad[0], 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[1], aqm_grad[1], 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[2], aqm_grad[2], 1e-6 );
  
  CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[0], Matrix3D(aqm_hess[0]), 1e-6 );
  CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[3], Matrix3D(aqm_hess[1]), 1e-6 );
  CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[5], Matrix3D(aqm_hess[2]), 1e-6 );
                                               
    // Test calculation of diagonal from diagonal corner hessian data.
  SymMatrix3D diags[9] = { hesss[ 0].upper(), hesss[ 3].upper(), hesss[ 5].upper(),
                           hesss[ 6].upper(), hesss[ 9].upper(), hesss[11].upper(),
                           hesss[12].upper(), hesss[15].upper(), hesss[17].upper() };
  aqm_avg = aqm.average_corner_hessian_diagonals( TRIANGLE, fixed, 3,
                                         vals, grads, diags,
                                         aqm_grad, aqm_hess, 
                                         err );
  ASSERT_NO_ERROR(err);

  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_avg, aqm_avg, 1e-6 );
  
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[0], aqm_grad[0], 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[1], aqm_grad[1], 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[2], aqm_grad[2], 1e-6 );
  
  CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[0], Matrix3D(aqm_hess[0]), 1e-6 );
  CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[3], Matrix3D(aqm_hess[1]), 1e-6 );
  CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[5], Matrix3D(aqm_hess[2]), 1e-6 );
}

void AveragingQMTest::check_hessian_diagonal_fails( QualityMetric::AveragingMethod scheme )
{
    // define corner values, gradients and Hessians for a triangle
  double vals[] = { 1, 2, 0.5 };
  const Vector3D grads[9] = { Vector3D(1.0),
                              Vector3D(2.0),
                              Vector3D(3.0),
                              Vector3D(4.0),
                              Vector3D(5.0),
                              Vector3D(6.0),
                              Vector3D(7.0),
                              Vector3D(8.0),
                              Vector3D(9.0) };
  const Matrix3D hesss[18] = {
    Matrix3D(1.0), Matrix3D(2.0), Matrix3D(-1.),
                   Matrix3D(0.5), Matrix3D(6.0),
                                  Matrix3D(1.1),
    Matrix3D(1.2), Matrix3D(0.3), Matrix3D(0.0),
                   Matrix3D(22.), Matrix3D(-2.),
                                  Matrix3D(4.3),
    Matrix3D(0.8), Matrix3D(1,7), Matrix3D(-.2),
                   Matrix3D(0.5), Matrix3D(7.0),
                                  Matrix3D(1.5)
  };
  
  Vector3D vtx_grad[3];
  SymMatrix3D vtx_hess[6];
  MsqError err;
  AveragingQM aqm(scheme);
  uint32_t fixed = 0x4; // mark third vertex as fixed
  
  aqm.average_corner_hessian_diagonals( TRIANGLE, fixed, 3,
                                        vals, grads, hesss,
                                        vtx_grad, vtx_hess, 
                                        err );
  CPPUNIT_ASSERT(err);
  
  const SymMatrix3D diags[9] = { SymMatrix3D(1.0), 
                                 SymMatrix3D(0.5),
                                 SymMatrix3D(1.0),
                                 SymMatrix3D(1.2),
                                 SymMatrix3D(22.),
                                 SymMatrix3D(4.3),
                                 SymMatrix3D(0.8),
                                 SymMatrix3D(0.5),
                                 SymMatrix3D(1.5) };
  
  aqm.average_corner_hessian_diagonals( TRIANGLE, fixed, 3,
                                        vals, grads, diags,
                                        vtx_grad, vtx_hess, 
                                        err );
  CPPUNIT_ASSERT(err);
}

void AveragingQMTest::test_average_corner_hessian_diagonals_linear()
  { check_pmean_hessian_diagonals( QualityMetric::LINEAR, 1, 1, true ); }

void AveragingQMTest::test_average_corner_hessian_diagonals_rms()
  { check_pmean_hessian_diagonals( QualityMetric::RMS, 2, 0.5, true ); }

void AveragingQMTest::test_average_corner_hessian_diagonals_hms()
  { check_pmean_hessian_diagonals( QualityMetric::HMS, -2, -0.5, true ); }

void AveragingQMTest::test_average_corner_hessian_diagonals_minimum()
  { check_hessian_diagonal_fails( QualityMetric::MINIMUM ); }

void AveragingQMTest::test_average_corner_hessian_diagonals_maximum()
  { check_hessian_diagonal_fails( QualityMetric::MAXIMUM ); }

void AveragingQMTest::test_average_corner_hessian_diagonals_harmonic()
  { check_pmean_hessian_diagonals( QualityMetric::HARMONIC, -1, -1, true ); }

void AveragingQMTest::test_average_corner_hessian_diagonals_geometric()
  { check_hessian_diagonal_fails( QualityMetric::GEOMETRIC ); }

void AveragingQMTest::test_average_corner_hessian_diagonals_sum()
  { check_pmean_hessian_diagonals( QualityMetric::SUM, 1, 1, false ); }

void AveragingQMTest::test_average_corner_hessian_diagonals_sum_squared()
  { check_pmean_hessian_diagonals( QualityMetric::SUM_SQUARED, 2, 1, false ); }

void AveragingQMTest::test_average_corner_hessian_diagonals_standard_deviation()
  { check_hessian_diagonal_fails( QualityMetric::STANDARD_DEVIATION ); }

void AveragingQMTest::test_average_corner_hessian_diagonals_max_over_min()
  { check_hessian_diagonal_fails( QualityMetric::MAX_OVER_MIN ); }

void AveragingQMTest::test_average_corner_hessian_diagonals_max_minus_min()
  { check_hessian_diagonal_fails( QualityMetric::MAX_MINUS_MIN ); }

void AveragingQMTest::test_average_corner_hessian_diagonals_sum_of_ratios_squared()
  { check_hessian_diagonal_fails( QualityMetric::SUM_OF_RATIOS_SQUARED ); }

void AveragingQMTest::check_pmean_hessian( QualityMetric::AveragingMethod scheme,
                                 double inner_power,
                                 double outer_power,
                                 bool scale )
{
    // define corner values, gradients and Hessians for a triangle
  double vals[] = { 1, 2, 0.5 };
  const Vector3D grads[9] = { Vector3D(0,1,1),
                              Vector3D(0,2,2),
                              Vector3D(0,3,3),
                              Vector3D(0,4,4),
                              Vector3D(0,5,5),
                              Vector3D(0,6,6),
                              Vector3D(0,7,7),
                              Vector3D(0,8,8),
                              Vector3D(0,9,9) };
  Matrix3D hesss[18] = {
    Matrix3D(1.0), Matrix3D(2.0), Matrix3D(-1.),
                   Matrix3D(0.5), Matrix3D(6.0),
                                  Matrix3D(1.1),
    Matrix3D(1.2), Matrix3D(0.3), Matrix3D(0.0),
                   Matrix3D(22.), Matrix3D(-2.),
                                  Matrix3D(4.3),
    Matrix3D(0.8), Matrix3D(1,7), Matrix3D(-.2),
                   Matrix3D(0.5), Matrix3D(7.0),
                                  Matrix3D(1.5)
  };
    // change some values so we catch transpositional errors
  for (unsigned i = 0; i < sizeof(hesss)/sizeof(hesss[0]); ++i)
    hesss[i][2][2] = 0.0;

  Vector3D aqm_grad[3], exp_grad[3];
  Matrix3D aqm_hess[6], exp_hess[6];;
  MsqPrintError err(std::cout);
  AveragingQM aqm(scheme);
  uint32_t fixed = 0x4; // mark third vertex as fixed
  double aqm_avg, exp_avg;
  
  exp_avg = pmean_of_triangle_corner_hessians( inner_power, 
                                               outer_power,
                                               vals,
                                               grads, hesss,
                                               exp_grad, exp_hess,
                                               scale );
                                               


  aqm_avg = aqm.average_corner_hessians( TRIANGLE, fixed, 3,
                                         vals, grads, hesss,
                                         aqm_grad, aqm_hess, 
                                         err );
  ASSERT_NO_ERROR(err);

  CPPUNIT_ASSERT_DOUBLES_EQUAL( exp_avg, aqm_avg, 1e-6 );
  
    // only check gradient terms for free vertices
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[0], aqm_grad[0], 1e-6 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( exp_grad[1], aqm_grad[1], 1e-6 );
    
    // only check Hessian terms for free vertices
  CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[0], aqm_hess[0], 1e-6 );
  CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[1], aqm_hess[1], 1e-6 );
  CPPUNIT_ASSERT_MATRICES_EQUAL( exp_hess[3], aqm_hess[3], 1e-6 );
}

void AveragingQMTest::check_hessian_fails( QualityMetric::AveragingMethod scheme )
{
    // define corner values, gradients and Hessians for a triangle
  double vals[] = { 1, 2, 0.5 };
  const Vector3D grads[9] = { Vector3D(1.0),
                              Vector3D(2.0),
                              Vector3D(3.0),
                              Vector3D(4.0),
                              Vector3D(5.0),
                              Vector3D(6.0),
                              Vector3D(7.0),
                              Vector3D(8.0),
                              Vector3D(9.0) };
  const Matrix3D hesss[18] = {
    Matrix3D(1.0), Matrix3D(2.0), Matrix3D(-1.),
                   Matrix3D(0.5), Matrix3D(6.0),
                                  Matrix3D(1.1),
    Matrix3D(1.2), Matrix3D(0.3), Matrix3D(0.0),
                   Matrix3D(22.), Matrix3D(-2.),
                                  Matrix3D(4.3),
    Matrix3D(0.8), Matrix3D(1,7), Matrix3D(-.2),
                   Matrix3D(0.5), Matrix3D(7.0),
                                  Matrix3D(1.5)
  };
  
  Vector3D vtx_grad[3];
  Matrix3D vtx_hess[6];
  MsqError err;
  AveragingQM aqm(scheme);
  uint32_t fixed = 0x4; // mark third vertex as fixed
  
  aqm.average_corner_hessians( TRIANGLE, fixed, 3,
                               vals, grads, hesss,
                               vtx_grad, vtx_hess, 
                               err );
  CPPUNIT_ASSERT(err);
}

void AveragingQMTest::test_average_corner_hessians_linear()
  { check_pmean_hessian( QualityMetric::LINEAR, 1, 1, true ); }

void AveragingQMTest::test_average_corner_hessians_rms()
  { check_pmean_hessian( QualityMetric::RMS, 2, 0.5, true ); }

void AveragingQMTest::test_average_corner_hessians_hms()
  { check_pmean_hessian( QualityMetric::HMS, -2, -0.5, true ); }

void AveragingQMTest::test_average_corner_hessians_minimum()
  { check_hessian_fails( QualityMetric::MINIMUM ); }

void AveragingQMTest::test_average_corner_hessians_maximum()
  { check_hessian_fails( QualityMetric::MAXIMUM ); }

void AveragingQMTest::test_average_corner_hessians_harmonic()
  { check_pmean_hessian( QualityMetric::HARMONIC, -1, -1, true ); }

void AveragingQMTest::test_average_corner_hessians_geometric()
  { check_hessian_fails( QualityMetric::GEOMETRIC ); }

void AveragingQMTest::test_average_corner_hessians_sum()
  { check_pmean_hessian( QualityMetric::SUM, 1, 1, false ); }

void AveragingQMTest::test_average_corner_hessians_sum_squared()
  { check_pmean_hessian( QualityMetric::SUM_SQUARED, 2, 1, false ); }

void AveragingQMTest::test_average_corner_hessians_standard_deviation()
  { check_hessian_fails( QualityMetric::STANDARD_DEVIATION ); }

void AveragingQMTest::test_average_corner_hessians_max_over_min()
  { check_hessian_fails( QualityMetric::MAX_OVER_MIN ); }

void AveragingQMTest::test_average_corner_hessians_max_minus_min()
  { check_hessian_fails( QualityMetric::MAX_MINUS_MIN ); }

void AveragingQMTest::test_average_corner_hessians_sum_of_ratios_squared()
  { check_hessian_fails( QualityMetric::SUM_OF_RATIOS_SQUARED ); }
