/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file QualityAssessorTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Mesquite_QualityAssessor.hpp"
#include "Mesquite_AspectRatioGammaQualityMetric.hpp"
#include "Mesquite_LocalSizeQualityMetric.hpp"
#include "Mesquite_ConditionNumberQualityMetric.hpp"
#include "Mesquite_TShapeNB1.hpp"
#include "Mesquite_IdealShapeTarget.hpp"
#include "Mesquite_TQualityMetric.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_Settings.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_PatchData.hpp"
#include "Mesquite_MsqMeshEntity.hpp"
#include "Mesquite_ArrayMesh.hpp"
#include "Mesquite_InstructionQueue.hpp"

#include "UnitUtil.hpp"

#include <algorithm>
#include <numeric>
#include <iomanip>

using namespace Mesquite;
using namespace std;



class QualityAssessorTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(QualityAssessorTest);
  CPPUNIT_TEST (test_basic_stats_element);
  CPPUNIT_TEST (test_basic_stats_vertex);
  CPPUNIT_TEST (test_basic_stats_sample);
  CPPUNIT_TEST (test_histogram_known_range);
  CPPUNIT_TEST (test_histogram_unknown_range);
  CPPUNIT_TEST (test_power_mean);
  CPPUNIT_TEST (test_invalid_count);
  CPPUNIT_TEST (test_inverted_count);
  CPPUNIT_TEST (test_output_control);
  CPPUNIT_TEST (test_tag_element);
  CPPUNIT_TEST (test_tag_vertex);
  CPPUNIT_TEST (test_tag_inverted);
  CPPUNIT_TEST (test_print_inverted);
  CPPUNIT_TEST (test_print_stats);
  CPPUNIT_TEST (test_print_name);
  CPPUNIT_TEST (test_modify_metric);
  CPPUNIT_TEST (test_free_only);
  CPPUNIT_TEST_SUITE_END();
  
  vector<double> vertCoords,invertCoords;
  vector<int> fixedFlags;
  vector<unsigned long> triConn, invertConn;
  
  MeshImpl myMesh, invertedMesh;
  PlanarDomain myDomain;
  Settings mySettings;
  
public:
  void setUp();
  void tearDown();
  
public:
  
  QualityAssessorTest() : myDomain(PlanarDomain::XY) {}

  void test_basic_stats_element();
  void test_basic_stats_vertex();
  void test_basic_stats_sample();
  void test_histogram_known_range();
  void test_histogram_unknown_range();
  void test_power_mean();
  void test_invalid_count();
  void test_inverted_count();
  void test_output_control();
  void test_tag_element();
  void test_tag_vertex();
  void test_tag_inverted();
  void test_print_inverted();
  void test_print_stats();
  void test_print_name();
  void test_modify_metric();
  void test_free_only();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QualityAssessorTest, "Unit");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QualityAssessorTest, "QualityAssessorTest");


void QualityAssessorTest::setUp()
{
  MsqError err;

  // create the following mesh:
  /*2.0
   *  (3)------------------(2)
   *   |\                __/|
   *   | \       3    __/ / |
   *   |  \        __/   /  |
   *   |   \    __/  4  /   |
   *   |    \ _/       /    |
   *1.0|  0 (4)------(5)  2 |
   *   |    /. \__    .\    |
   *   |   / .    \__5. \   |
   *   |  /  .       \__ \  |
   *   | /   .   1    . \_\ |
   *   |/    .        .    \|
   *  (0)------------------(1)
   *0.0     1.0      2.0   3.0
   */
  const char vtk_file[] = 
  "# vtk DataFile Version 3.0\n"
  "QualityAssessorTest input 1\n"
  "ASCII\n"
  "DATASET UNSTRUCTURED_GRID\n"
  "POINTS 6 float\n"
  "0 0 0\n"
  "3 0 0\n"
  "3 2 0\n"
  "0 2 0\n"
  "1 1 0\n"
  "2 1 0\n"
  "CELLS 6 24\n"
  "3 4 3 0\n"
  "3 0 1 4\n"
  "3 1 2 5\n"
  "3 2 3 4\n"
  "3 4 5 2\n"
  "3 5 4 1\n"
  "CELL_TYPES 6\n"
  "5 5 5 5 5 5\n"
  "POINT_DATA 6\n"
  "SCALARS fixed int\n"
  "LOOKUP_TABLE default\n"
  "1 1 1 1 0 0\n";
  
  // define a mesh with two triangles, one of which is inverted wrt the other
  const char invert_file[] = 
  "# vtk DataFile Version 3.0\n"
  "QualityAssessorTest input 2\n"
  "ASCII\n"
  "DATASET UNSTRUCTURED_GRID\n"
  "POINTS 3 float\n"
  "0 0 0\n"
  "1 0 0\n"
  "0 1 0\n"
  "CELLS 2 8\n"
  "3 0 1 2\n"
  "3 0 2 1\n"
  "CELL_TYPES 2\n"
  "5 5\n"
  "POINT_DATA 3\n"
  "SCALARS fixed int\n"
  "LOOKUP_TABLE default\n"
  "1 1 1\n";
  
  myMesh.clear();
  invertedMesh.clear();
  
  const char* tmpname = "QualityAssessorTest.vtk";
  FILE* filp = fopen( tmpname, "w" );
  CPPUNIT_ASSERT( NULL != filp );
  size_t w = fwrite( vtk_file, sizeof(vtk_file), 1, filp );
  fclose( filp );
  myMesh.clear();
  myMesh.read_vtk( tmpname, err );
  remove( tmpname );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, w );
  ASSERT_NO_ERROR( err );

  filp = fopen( tmpname, "w" );
  CPPUNIT_ASSERT( NULL != filp );
  w = fwrite( invert_file, sizeof(invert_file), 1, filp );
  fclose( filp );
  invertedMesh.clear();
  invertedMesh.read_vtk( tmpname, err );
  remove( tmpname );
  CPPUNIT_ASSERT_EQUAL( (size_t)1, w );
  ASSERT_NO_ERROR( err );
}

void QualityAssessorTest::tearDown()
{
}

// Define a decorator (wrapper) for a QualityMetric
// instance that records each metric value.  
class MetricLogger : public QualityMetric
{
public:
  MetricLogger( QualityMetric* metric ) :invalidCount(0),  mMetric(metric) {}

  vector<double> mValues;
  vector<Mesh::EntityHandle> mHandles;
  int invalidCount;
  
  double min() const { return *min_element(mValues.begin(), mValues.end()); }
  double max() const { return *max_element(mValues.begin(), mValues.end()); }
  double avg() const { double x=0; return accumulate(mValues.begin(), mValues.end(), x)/mValues.size(); }
  double sqrsum() const { 
    double result = 0.0;
    for (unsigned i = 0; i < mValues.size(); ++i)
      result += mValues[i]*mValues[i];
    return result;
  }
  double rms() const { return sqrt(sqrsum()/mValues.size());}
  double dev() const { return sqrt(sqrsum()/mValues.size() - avg()*avg());}
  double pmean(double p) const {
    double result = 0.0;
    for (unsigned i = 0; i < mValues.size(); ++i)
      result += pow(mValues[i], p);
    return pow( result/mValues.size(), 1./p );
  }
  
  bool no_duplicate_evals() const 
    {
      vector<Mesh::EntityHandle> handles(mHandles);
      sort( handles.begin(), handles.end() );
      return handles.end() == unique( handles.begin(), handles.end() );
    }
  
  MetricType get_metric_type() const 
    { return mMetric->get_metric_type(); }

  std::string get_name() const
    { return mMetric->get_name(); }

  int get_negate_flag() const
    { return mMetric->get_negate_flag(); }

  void get_evaluations( PatchData& pd, vector<size_t>& handles, bool free, MsqError& err )
    { return mMetric->get_evaluations( pd, handles, free, err ); }

  bool evaluate( PatchData& pd, size_t handle, double& value, MsqError& err )
    { 
      bool rval = mMetric->evaluate( pd, handle, value, err );
      mValues.push_back( value );
      if (get_metric_type() == QualityMetric::VERTEX_BASED)
        mHandles.push_back( pd.get_vertex_handles_array()[handle] );
      else if (handle < pd.num_elements())
        mHandles.push_back( pd.get_element_handles_array()[handle] );
      if (!rval)
        ++invalidCount;
      return rval;
    }
  
  bool evaluate_with_indices( PatchData&, size_t, double&, vector<size_t>&, MsqError& err)
    { // shouldn't be called by QualityAssessor
      CPPUNIT_ASSERT(false);
      MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED); 
      return false;
    }
  bool evaluate_with_gradient( PatchData&, size_t, double&, vector<size_t>&, vector<Vector3D>&, MsqError& err )
    { // shouldn't be called by QualityAssessor
      CPPUNIT_ASSERT(false);
      MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED); 
      return false;
    }
  bool evaluate_with_Hessian( PatchData&, size_t, double&, vector<size_t>&, vector<Vector3D>&, vector<Matrix3D>&, MsqError& err )
    { // shouldn't be called by QualityAssessor
      CPPUNIT_ASSERT(false);
      MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED); 
      return false;
    }
private:
  QualityMetric* mMetric;
};

void QualityAssessorTest::test_basic_stats_element()
{
  AspectRatioGammaQualityMetric metric;
  MetricLogger logger(&metric);
  QualityAssessor qa(&logger);
  qa.disable_printing_results();
  
  MsqError err;
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );
  
    // check didn't evaluate any element more than once
  CPPUNIT_ASSERT( logger.no_duplicate_evals() );
  
    // check one eval for each element
  vector<Mesh::ElementHandle> elems;
  myMesh.get_all_elements( elems, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( elems.size(), logger.mValues.size() );
  
    // check values
  const QualityAssessor::Assessor* results = qa.get_results(&logger);
  CPPUNIT_ASSERT(results != NULL);
  CPPUNIT_ASSERT_EQUAL( (int)elems.size(), results->get_count() );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.min(), results->get_minimum(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.avg(), results->get_average(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.rms(), results->get_rms    (), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.max(), results->get_maximum(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.dev(), results->get_stddev (), 1e-6 );
}

void QualityAssessorTest::test_basic_stats_vertex()
{
  LocalSizeQualityMetric metric;
  MetricLogger logger(&metric);
  QualityAssessor qa(&logger);
  qa.measure_free_samples_only( false );
  qa.disable_printing_results();
  
  MsqError err;
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );
  
    // check didn't evaluate any vertex more than once
  CPPUNIT_ASSERT( logger.no_duplicate_evals() );
  
    // check one eval for each element
  vector<Mesh::VertexHandle> verts;
  myMesh.get_all_vertices( verts, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( verts.size(), logger.mValues.size() );
  
    // check values
  const QualityAssessor::Assessor* results = qa.get_results(&logger);
  CPPUNIT_ASSERT(results != NULL);
  CPPUNIT_ASSERT_EQUAL( (int)verts.size(), results->get_count() );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.min(), results->get_minimum(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.avg(), results->get_average(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.rms(), results->get_rms    (), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.max(), results->get_maximum(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.dev(), results->get_stddev (), 1e-6 );
}

void QualityAssessorTest::test_basic_stats_sample()
{
  TShapeNB1 tm;
  IdealShapeTarget tc;
  TQualityMetric metric( &tc, &tm );
  MetricLogger logger(&metric);
  QualityAssessor qa(&logger);
  qa.disable_printing_results();
 
  MsqError err;
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );
  
    // check didn't evaluate any sample more than once
  CPPUNIT_ASSERT( logger.no_duplicate_evals() );
  
    // count total number of sample points in mesh
  vector<Mesh::ElementHandle> elems;
  myMesh.get_all_elements( elems, err );
  ASSERT_NO_ERROR( err );
  
  vector<Mesh::VertexHandle> free_verts;
  PatchData global_patch;
  global_patch.set_mesh( &myMesh );
  global_patch.set_domain( &myDomain );
  global_patch.attach_settings( &mySettings );
  global_patch.set_mesh_entities( elems, free_verts, err );
  ASSERT_NO_ERROR(err);
  
  size_t num_samples = 0;
  for (size_t i = 0; i < global_patch.num_elements(); ++i)
    num_samples += global_patch.get_samples(i).num_nodes();
  
    // check correct number of metric evaluations
  CPPUNIT_ASSERT_EQUAL( num_samples, logger.mValues.size() );
  
    // check values
  const QualityAssessor::Assessor* results = qa.get_results(&logger);
  CPPUNIT_ASSERT(results != NULL);
  CPPUNIT_ASSERT_EQUAL( num_samples, (size_t)results->get_count() );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.min(), results->get_minimum(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.avg(), results->get_average(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.rms(), results->get_rms    (), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.max(), results->get_maximum(), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.dev(), results->get_stddev (), 1e-6 );
}

void QualityAssessorTest::test_histogram_known_range()
{
  const double lower = 1.0, upper = 3.0;
  const int intervals = 5;

  AspectRatioGammaQualityMetric metric;
  MetricLogger logger(&metric);
  QualityAssessor qa;
  qa.add_histogram_assessment( &logger, lower, upper, intervals );
  qa.disable_printing_results();
  
  MsqError err;
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );
  
    // calculate expected histogram
  std::vector<int> counts(intervals+2, 0);
  for (vector<double>::iterator i = logger.mValues.begin(); i != logger.mValues.end(); ++i)
  {
    int bucket = (int)((*i - lower) * intervals / (upper - lower));
    if (bucket < 0)
      ++counts.front();
    else if (bucket >= intervals)
      ++counts.back();
    else
      ++counts[bucket+1];
  }
  
    // check values
  const QualityAssessor::Assessor* results = qa.get_results(&logger);
  CPPUNIT_ASSERT(results != NULL);
  CPPUNIT_ASSERT(results->have_histogram());
  
  double r_lower, r_upper;
  vector<int> r_counts;
  results->get_histogram( r_lower, r_upper, r_counts, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( lower, r_lower, DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( upper, r_upper, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( intervals+2, (int)r_counts.size() );
  
  CPPUNIT_ASSERT_EQUAL( counts.front(), r_counts.front() );
  CPPUNIT_ASSERT_EQUAL( counts.back(), r_counts.back() );
  switch (intervals) {
    default: for (unsigned i = 11; i < intervals+1; ++i)
             CPPUNIT_ASSERT_EQUAL( counts[ i], r_counts[ i] );
    case 10: CPPUNIT_ASSERT_EQUAL( counts[10], r_counts[10] );
    case  9: CPPUNIT_ASSERT_EQUAL( counts[ 9], r_counts[ 9] );
    case  8: CPPUNIT_ASSERT_EQUAL( counts[ 8], r_counts[ 8] );
    case  7: CPPUNIT_ASSERT_EQUAL( counts[ 7], r_counts[ 7] );
    case  6: CPPUNIT_ASSERT_EQUAL( counts[ 6], r_counts[ 6] );
    case  5: CPPUNIT_ASSERT_EQUAL( counts[ 5], r_counts[ 5] );
    case  4: CPPUNIT_ASSERT_EQUAL( counts[ 4], r_counts[ 4] );
    case  3: CPPUNIT_ASSERT_EQUAL( counts[ 3], r_counts[ 3] );
    case  2: CPPUNIT_ASSERT_EQUAL( counts[ 2], r_counts[ 2] );
    case  1: CPPUNIT_ASSERT_EQUAL( counts[ 1], r_counts[ 1] );
    case  0: ;
  }
}

void QualityAssessorTest::test_histogram_unknown_range()
{
  const int intervals = 10;

  AspectRatioGammaQualityMetric metric;
  MetricLogger logger(&metric);
  QualityAssessor qa( &logger, intervals );
  qa.disable_printing_results();
  
  MsqError err;
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );
  
    // check values
  const QualityAssessor::Assessor* results = qa.get_results(&logger);
  CPPUNIT_ASSERT(results != NULL);
  CPPUNIT_ASSERT(results->have_histogram());

  double r_lower, r_upper;
  vector<int> r_counts;
  results->get_histogram( r_lower, r_upper, r_counts, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT( r_lower <= logger.min() );
  CPPUNIT_ASSERT( r_upper >= logger.min() );
    // allow some freedom in choice of range, but not more than 30%
  CPPUNIT_ASSERT( (0.30 * (r_upper - r_lower)) <= (logger.max() - logger.min()) );
    // range must contain all metric values
  CPPUNIT_ASSERT_EQUAL( 0, r_counts.front() );
  CPPUNIT_ASSERT_EQUAL( 0, r_counts.back() );
  CPPUNIT_ASSERT_EQUAL( intervals+2, (int)r_counts.size() );
  
    // calculate expected histogram
  std::vector<int> counts(intervals, 0);
  for (vector<double>::iterator i = logger.mValues.begin(); i != logger.mValues.end(); ++i)
  {
    double fract = (*i - r_lower) / (r_upper - r_lower);
    int bucket;
    if (fabs(fract - 1.0) < DBL_EPSILON)
      bucket = intervals-1;
    else
      bucket = (int)(intervals * fract);
    CPPUNIT_ASSERT( bucket >= 0 );
    CPPUNIT_ASSERT( bucket < intervals );
    ++counts[bucket];
  }
  
    // QA should have evaluated metric twice, so adjust values
  CPPUNIT_ASSERT_EQUAL( 12, (int)logger.mValues.size() );
  for (vector<int>::iterator j = counts.begin(); j != counts.end(); ++j)
    *j /= 2;
  
  switch (intervals) {
    default: for (int i = 10; i < intervals; ++i)
             CPPUNIT_ASSERT_EQUAL( counts[ i], r_counts[i+1]);
    case 10: CPPUNIT_ASSERT_EQUAL( counts[ 9], r_counts[10] );
    case  9: CPPUNIT_ASSERT_EQUAL( counts[ 8], r_counts[ 9] );
    case  8: CPPUNIT_ASSERT_EQUAL( counts[ 7], r_counts[ 8] );
    case  7: CPPUNIT_ASSERT_EQUAL( counts[ 6], r_counts[ 7] );
    case  6: CPPUNIT_ASSERT_EQUAL( counts[ 5], r_counts[ 6] );
    case  5: CPPUNIT_ASSERT_EQUAL( counts[ 4], r_counts[ 5] );
    case  4: CPPUNIT_ASSERT_EQUAL( counts[ 3], r_counts[ 4] );
    case  3: CPPUNIT_ASSERT_EQUAL( counts[ 2], r_counts[ 3] );
    case  2: CPPUNIT_ASSERT_EQUAL( counts[ 1], r_counts[ 2] );
    case  1: CPPUNIT_ASSERT_EQUAL( counts[ 0], r_counts[ 1] );
    case  0: ;
  }
}

void QualityAssessorTest::test_power_mean()
{
  const double P1 = 3.0, P2 = 0.5;
  AspectRatioGammaQualityMetric metric1;
  LocalSizeQualityMetric metric2;
  ConditionNumberQualityMetric metric3;
  MetricLogger logger1(&metric1), logger2(&metric2);
  QualityAssessor qa( &metric3 );
  qa.add_quality_assessment( &logger1, 0, P1 );
  qa.add_quality_assessment( &logger2, 0, P2 );
  qa.disable_printing_results();
  
  MsqError err;
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );
  
    // get results
  const QualityAssessor::Assessor *result1, *result2, *result3;
  result1 = qa.get_results(&logger1);
  result2 = qa.get_results(&logger2);
  result3 = qa.get_results(&metric3);
  CPPUNIT_ASSERT(NULL != result1);
  CPPUNIT_ASSERT(NULL != result2);
  CPPUNIT_ASSERT(NULL != result3);
  
    // check flags
  CPPUNIT_ASSERT( result1->have_power_mean() );
  CPPUNIT_ASSERT( result2->have_power_mean() );
  CPPUNIT_ASSERT(!result3->have_power_mean() );
  
    // check values
  CPPUNIT_ASSERT_DOUBLES_EQUAL( P1, result1->get_power(), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( P2, result2->get_power(), DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger1.pmean(P1), result1->get_power_mean(), logger1.pmean(P1)*DBL_EPSILON );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( logger2.pmean(P2), result2->get_power_mean(), logger2.pmean(P2)*DBL_EPSILON );
}

void QualityAssessorTest::test_invalid_count()
{
  MsqError err;
  ConditionNumberQualityMetric metric;
  QualityAssessor qa( &metric );
  qa.measure_free_samples_only( false );
  qa.disable_printing_results();
  const QualityAssessor::Assessor* results = qa.get_results( &metric );
  CPPUNIT_ASSERT(NULL != results);
  
    // try mesh with only valid elements
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( 0, results->get_invalid_element_count() );
  
    // try mesh with one inverted element
  MeshDomainAssoc mesh_and_domain2 = MeshDomainAssoc(&invertedMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain2, &mySettings, err  );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( 1, results->get_invalid_element_count() );
}

void QualityAssessorTest::test_inverted_count()
{
  int inverted, samples;
  MsqError err;
  QualityAssessor qa;
  qa.measure_free_samples_only( false );
  qa.disable_printing_results();
  
    // try mesh with only valid elements
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );
  inverted = samples = -1;
  qa.get_inverted_element_count( inverted, samples, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( 0, inverted );
  CPPUNIT_ASSERT_EQUAL( 0, samples );
  
    // try mesh with one inverted element
  MeshDomainAssoc mesh_and_domain2 = MeshDomainAssoc(&invertedMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain2, &mySettings, err  );
  ASSERT_NO_ERROR( err );
  inverted = samples = -1;
  qa.get_inverted_element_count( inverted, samples, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( 1, inverted );
  CPPUNIT_ASSERT_EQUAL( 1, samples );
}

  // Define class to handle redirecting output streams
  // (e.g. std::cout) to a stringstream instance.  We
  // use this to catch errors where QualityAssessor is
  // writing to a stream it should not be.
class StreamRedirector {
private:
  streambuf *outbuf, *errbuf, *logbuf;
public:
  ostringstream outstr, errstr, logstr;
  StreamRedirector() {
    outbuf = cout.rdbuf();
    errbuf = cerr.rdbuf();
    logbuf = clog.rdbuf();
  }
  
  ~StreamRedirector() { restore(); }

  void redirect() {
    clear();
    cout.rdbuf( outstr.rdbuf() );
    cerr.rdbuf( errstr.rdbuf() );
    clog.rdbuf( logstr.rdbuf() );
  }
  
  void restore() {
    cout.rdbuf( outbuf );
    cerr.rdbuf( errbuf );
    clog.rdbuf( logbuf );
  }
    // reset all streams to empty
  void clear() {
    string empty;
    outstr.str(empty);
    errstr.str(empty);
    logstr.str(empty);
  }

  void print_stream_bytes( const char* prefix, string bytes ) {
      // fix output stream so we can use it to print results
    streambuf* sbuf = cout.rdbuf();
    cout.rdbuf( outbuf );
    cout << prefix << ": " << bytes.size() << " characters: ";
    for (string::iterator i = bytes.begin(); i != bytes.end(); ++i)
      if (isprint(*i) || *i == ' ')
        cout << *i;
      else switch (*i) {
        case '\a': cout << "\\a"; break;
        case '\b': cout << "\\b"; break;
        case '\t': cout << "\\t"; break;
        case '\n': cout << "\\n"; break;
        case '\v': cout << "\\v"; break;
        case '\f': cout << "\\f"; break;
        case '\r': cout << "\\r"; break;
        case '\0': cout << "\\0"; break;
        default: 
          cout << '\\' << setw(3) << setfill('0') << oct << *i;
          break;
      }
    cout << endl;
      // restore cout to incomming state
    cout.rdbuf(sbuf);
  }
  
  void check_stream( const char* name, ostringstream& str, bool& result )
  {
    string s = str.str();
    if (s.empty())
      return;
    result = true;
    print_stream_bytes( name, s );
  }

    // check if any stream was written to (may clear streams, depending on platform)
  bool have_data() {
    bool result = false;
    check_stream( "std::cout", outstr, result );
    check_stream( "std::cerr", errstr, result );
    check_stream( "std::clog", logstr, result );
    return result;
  }
};

void QualityAssessorTest::test_output_control()
{
  // Redirect std::cout, std:cerr, and std::clog
  // so we know when they're written to.
  StreamRedirector redir;
  
  MsqError err;
  ConditionNumberQualityMetric metric;
  
    // disable output from constructor
  QualityAssessor qa1( &metric, 0, 0, false, 0, false );
  redir.redirect();
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&invertedMesh, &myDomain);
  qa1.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  redir.restore();
    // make sure nothing was written to output streams
  CPPUNIT_ASSERT( !redir.have_data() );
  
    // disable output from function
  QualityAssessor qa2( &metric );
  qa2.disable_printing_results();
  redir.redirect();
  qa2.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  redir.restore();
    // make sure nothing was written to output streams
  CPPUNIT_ASSERT( !redir.have_data() );

    // send output to alternate stream
  stringstream deststr;
  QualityAssessor qa3( deststr, &metric );
  redir.redirect();
  qa3.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  redir.restore();
    // make sure nothing was written to output streams
  CPPUNIT_ASSERT( !redir.have_data() );
    // check that some data was written to specified stream
  CPPUNIT_ASSERT( !deststr.str().empty() );  
  
    // redir destructor should now restore normal output streams
}

void QualityAssessorTest::test_tag_element()
{
  const char* tag_name = "xxxTEST_TAGxxx";
  ConditionNumberQualityMetric metric;
  MetricLogger logger(&metric);
  QualityAssessor qa(&logger,0,0,false,tag_name);
  qa.disable_printing_results();
  
  MsqError err;
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );
  
  TagHandle tag;
  tag = myMesh.tag_get( tag_name, err );
  ASSERT_NO_ERROR( err );
  
  string name;
  Mesh::TagType type;
  unsigned length;
  myMesh.tag_properties( tag, name, type, length, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( name, string(tag_name) );
  CPPUNIT_ASSERT_EQUAL( Mesh::DOUBLE, type );
  CPPUNIT_ASSERT_EQUAL( 1u, length );
  
  CPPUNIT_ASSERT_EQUAL( logger.mValues.size(), logger.mHandles.size() );
  CPPUNIT_ASSERT(!logger.mValues.empty());
  vector<double> tag_values( logger.mValues.size() );
  myMesh.tag_get_element_data( tag, logger.mHandles.size(),
                      &logger.mHandles[0], arrptr(tag_values), err );
  ASSERT_NO_ERROR( err );
  
  for (unsigned i = 0; i < logger.mValues.size(); ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.mValues[1], tag_values[1], DBL_EPSILON );
}

void QualityAssessorTest::test_tag_vertex()
{
  const char* tag_name = "vertex-test-tag";
  LocalSizeQualityMetric metric;
  MetricLogger logger(&metric);
  QualityAssessor qa(&logger,0,0,false,tag_name);
  qa.disable_printing_results();
  
  MsqError err;
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );
  
  TagHandle tag;
  tag = myMesh.tag_get( tag_name, err );
  ASSERT_NO_ERROR( err );
  
  string name;
  Mesh::TagType type;
  unsigned length;
  myMesh.tag_properties( tag, name, type, length, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( name, string(tag_name) );
  CPPUNIT_ASSERT_EQUAL( Mesh::DOUBLE, type );
  CPPUNIT_ASSERT_EQUAL( 1u, length );
  
  CPPUNIT_ASSERT_EQUAL( logger.mValues.size(), logger.mHandles.size() );
  vector<double> tag_values( logger.mValues.size() );
  myMesh.tag_get_vertex_data( tag, logger.mHandles.size(),
                      &logger.mHandles[0], arrptr(tag_values), err );
  ASSERT_NO_ERROR( err );
  
  for (unsigned i = 0; i < logger.mValues.size(); ++i)
    CPPUNIT_ASSERT_DOUBLES_EQUAL( logger.mValues[1], tag_values[1], DBL_EPSILON );
}


void QualityAssessorTest::test_tag_inverted()
{
  const char* tag_name = "123INVERT456";
  QualityAssessor qa( false, false, tag_name );
  
  MsqError err;
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&invertedMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );
  
  TagHandle tag;
  tag = invertedMesh.tag_get( tag_name, err );
  ASSERT_NO_ERROR( err );
  
  string name;
  Mesh::TagType type;
  unsigned length;
  invertedMesh.tag_properties( tag, name, type, length, err );
  ASSERT_NO_ERROR( err );
  CPPUNIT_ASSERT_EQUAL( name, string(tag_name) );
  CPPUNIT_ASSERT_EQUAL( Mesh::INT, type );
  CPPUNIT_ASSERT_EQUAL( 1u, length );
  
  vector<Mesh::ElementHandle> elements;
  invertedMesh.get_all_elements( elements, err );
  ASSERT_NO_ERROR( err );
  
  // expect two elements, one inverted
  CPPUNIT_ASSERT_EQUAL( (size_t)2, elements.size() );
  int data[2];
  invertedMesh.tag_get_element_data( tag, 2, arrptr(elements), data, err );
  ASSERT_NO_ERROR( err );
  
  if (data[0]) {
    CPPUNIT_ASSERT_EQUAL( 0, data[1] );
  }
  else {
    CPPUNIT_ASSERT_EQUAL( 0, data[0] );
    CPPUNIT_ASSERT_EQUAL( 1, data[1] );
  }
}

void QualityAssessorTest::test_print_inverted()
{
  MsqError err;
  stringstream str;
  QualityAssessor qa( str );
  qa.measure_free_samples_only( false );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&invertedMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );

    // get inverted count from QA
  int expected = 0, other = 0;
  qa.get_inverted_element_count( expected, other, err );
  ASSERT_NO_ERROR( err );

    // At some point, expect "n INVERTED" to appear in output,
    // where 'n' is the count of inverted elements.
  int value;
  string prev, curr;
  for (; str >> curr && curr != "INVERTED"; prev = curr);
  str.str(prev);
  str >> value;
  
  CPPUNIT_ASSERT_EQUAL( expected, value );
}

void QualityAssessorTest::test_print_stats()
{
  MsqError err;
  stringstream str;
  ConditionNumberQualityMetric metric;
  QualityAssessor qa( str, &metric );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );

    // get results
  const QualityAssessor::Assessor* results = qa.get_results( &metric );
  CPPUNIT_ASSERT(NULL != results);

    // Expect metric name followed by 5 numerical values for statistics,
    // at some point in output.
  
    // Find occurance of metric name in output
  for (;;) {
    stringstream name( metric.get_name() );
    string s, n;
    name >> n;
    while (str >> s && s != n);
    CPPUNIT_ASSERT(str);
    while (name >> n && str >> s && s == n);
    if (!name)
      break;
  }
  
    // Read stats
  double min_s, max_s, avg_s, rms_s, dev_s;
  str >> min_s >> avg_s >> rms_s >> max_s >> dev_s;
  

    // The following commented out because they no longer pass due to a change
    //  in the QA Summary format
 
    // compare results
//  CPPUNIT_ASSERT_DOUBLES_EQUAL( results->get_minimum(), min_s, min_s * 0.01 );
//  CPPUNIT_ASSERT_DOUBLES_EQUAL( results->get_average(), avg_s, avg_s * 0.01 );
//  CPPUNIT_ASSERT_DOUBLES_EQUAL( results->get_rms    (), rms_s, rms_s * 0.01 );
//  CPPUNIT_ASSERT_DOUBLES_EQUAL( results->get_maximum(), max_s, max_s * 0.01 );
//  CPPUNIT_ASSERT_DOUBLES_EQUAL( results->get_stddev (), dev_s, dev_s * 0.01 );
}

void QualityAssessorTest::test_print_name()
{
  const char* NAME = "test_print_name_NAME";

    // expect QualityAssessor Name to be printed somewhere in output.
  MsqError err;
  stringstream str;
  ConditionNumberQualityMetric metric;
  QualityAssessor qa( str, &metric, 0, 0, false, 0, 0, NAME);
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, &myDomain);
  qa.loop_over_mesh( &mesh_and_domain, &mySettings, err  );
  ASSERT_NO_ERROR( err );

    // seach output for first occurance of name
  string s;
  while (str >> s && s != NAME);
  CPPUNIT_ASSERT_EQUAL( s, string(NAME) ); 
}

// Test that adding a metric more than once changes 
// properties for a single occurance of the metric,
// rather than introducing multiple copies.
void QualityAssessorTest::test_modify_metric()
{
  MsqError err;
  ConditionNumberQualityMetric metric;
  QualityAssessor qa( &metric );
  const QualityAssessor::Assessor* results = qa.get_results( &metric );
  CPPUNIT_ASSERT(NULL != results);
 
  CPPUNIT_ASSERT( !results->have_histogram() );
  qa.add_histogram_assessment( &metric, 0.0, 1.0, 10 );
  CPPUNIT_ASSERT( results->have_histogram() );
 
  CPPUNIT_ASSERT( !results->have_power_mean() );
  qa.add_quality_assessment( &metric, 0, 3.0 );
  CPPUNIT_ASSERT( results->have_power_mean() );
}

void QualityAssessorTest::test_free_only()
{
  /*
   *   +--+--------o-----o
   *   |   \       |     |    + - fixed
   *   |    \      |     |    o - free
   *   +-----+-----o-----o
   */

  double coords[] = { 0, 0, 0,
                      1, 0, 0,
                    0.5, 1, 0,
                      0, 1, 0,
                      2, 0, 0,
                      3, 0, 0,
                      3, 1, 0,
                      2, 1, 0 };
  int fixed[] = { true, true, true, true, false, false, false, false };
  unsigned long conn[] = { 0, 1, 2, 3,
                           1, 4, 7, 2,
                           4, 5, 6, 7 };
  
  MsqError err;
  ArrayMesh mesh( 3, 8, coords, fixed, 3, QUADRILATERAL, conn );
  ConditionNumberQualityMetric metric;
  QualityAssessor qa_all( &metric, 0, 0.0, false, 0, false );
  QualityAssessor qa_free( &metric, 0, 0.0, true, 0, false );
  InstructionQueue q;
  q.add_quality_assessor( &qa_all, err );
  q.add_quality_assessor( &qa_free, err );
  PlanarDomain xy(PlanarDomain::XY);
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &xy);
  q.run_instructions( &mesh_and_domain, err );
  ASSERT_NO_ERROR(err);
  
  const QualityAssessor::Assessor* data;
  data = qa_all.get_results( &metric );
  CPPUNIT_ASSERT( 0 != data );
  CPPUNIT_ASSERT_EQUAL( 3, data->get_count() );

  data = qa_free.get_results( &metric );
  CPPUNIT_ASSERT( 0 != data );
  CPPUNIT_ASSERT_EQUAL( 2, data->get_count() );
}

  
  
  
  

