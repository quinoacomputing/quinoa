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


#include "ElementPMeanP.hpp"
#include "VertexPMeanP.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "PatchDataInstances.hpp"
#include "UnitUtil.hpp"
#include "ElemSampleQM.hpp"

using namespace Mesquite;
using std::cout;
using std::cerr;
using std::endl;

class PMeanPMetricTest : public CppUnit::TestFixture
{

private:
  CPPUNIT_TEST_SUITE(PMeanPMetricTest);
  CPPUNIT_TEST(test_get_metric_type);
  CPPUNIT_TEST(test_get_element_evaluations);
  CPPUNIT_TEST(test_get_vertex_evaluations);
  CPPUNIT_TEST(test_element_evaluate);
  CPPUNIT_TEST(test_vertex_evaluate);
  CPPUNIT_TEST(test_indices);
  CPPUNIT_TEST(test_gradient);
  CPPUNIT_TEST(test_hessian);
  CPPUNIT_TEST(test_hessian_diagonal);
  CPPUNIT_TEST_SUITE_END();

  PatchData pd;

public:
  void setUp();

  void test_get_metric_type();
  void test_get_element_evaluations();
  void test_get_vertex_evaluations();
  void test_element_evaluate();
  void test_vertex_evaluate();
  void test_indices();
  void test_gradient();
  void test_hessian();
  void test_hessian_diagonal();
  
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PMeanPMetricTest, "PMeanPMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PMeanPMetricTest, "Unit");

void PMeanPMetricTest::setUp() {
  MsqError err;
  create_four_quads_patch( pd, err );
  CPPUNIT_ASSERT(!err);
}

// bogus metric for testing.  
// returns vertex index value for metric value, {vertex index, 0, 1} 
// for gradint wrt vertex , and { (h1, h2, h1+j2), {h1, 1, h1+h2}, {2*h1, 2*h2, h2} }
// for hessin where h1 and h2 are vertex indices.
class FauxMetric : public ElemSampleQM
{
public:
  MetricType get_metric_type() const { return ELEMENT_BASED; }
  std::string get_name() const { return "FauxMetrix"; }
  int get_negate_flag() const { return 1; }
  void get_evaluations( PatchData& pd, std::vector<size_t>& h, bool free, MsqError& );
  void get_element_evaluations( PatchData&, size_t, std::vector<size_t>&, MsqError& );
  bool evaluate( PatchData& pd, size_t h, double& v, MsqError& err );
  bool evaluate_with_indices( PatchData& pd, size_t h, double& v, 
                              std::vector<size_t>& indices, MsqError& err );
  bool evaluate_with_gradient( PatchData& pd, size_t h, double& v, 
                              std::vector<size_t>& indices, 
                              std::vector<Vector3D>& grads,
                              MsqError& err );
  bool evaluate_with_Hessian( PatchData& pd, size_t h, double& v, 
                              std::vector<size_t>& indices, 
                              std::vector<Vector3D>& grad,
                              std::vector<Matrix3D>& Hess,
                              MsqError& err );
};

void FauxMetric::get_evaluations( PatchData& pd, std::vector<size_t>& h, bool free, MsqError& err )
{
  h.clear();
  for (size_t i = 0; i < pd.num_elements(); ++i)
    get_element_evaluations( pd, i, h, err );
}

void FauxMetric::get_element_evaluations( PatchData& pd, size_t h, std::vector<size_t>& list, MsqError& )
{
  MsqMeshEntity& elem = pd.element_by_index(h);
  for (unsigned i = 0; i  < elem.corner_count(); ++i)
    list.push_back( handle( Sample(0,i), h ) );
}

bool FauxMetric::evaluate( PatchData& pd, size_t h, double& v, MsqError&  )
{
  size_t e = ElemSampleQM::elem( h );
  Sample s = ElemSampleQM::sample( h );
  size_t* verts = pd.element_by_index(e).get_vertex_index_array();
  CPPUNIT_ASSERT_EQUAL( (unsigned short)0, s.dimension );
  v = (double)(verts[s.number]);
  return true;
}

bool FauxMetric::evaluate_with_indices( PatchData& pd, size_t h, double& v, 
                              std::vector<size_t>& indices, MsqError& err )
{
  evaluate( pd, h, v, err );
  indices.resize(3);
  size_t e = ElemSampleQM::elem( h );
  Sample s = ElemSampleQM::sample( h );
  size_t* verts = pd.element_by_index(e).get_vertex_index_array();
  size_t n = pd.element_by_index(e).vertex_count();
  CPPUNIT_ASSERT_EQUAL( (unsigned short)0, s.dimension );
  indices[0] = verts[s.number];
  indices[1] = verts[(s.number+1)%n];
  indices[2] = verts[(s.number+n-1)%n];
  return true;
}

bool FauxMetric::evaluate_with_gradient( PatchData& pd, size_t h, double& v, 
                              std::vector<size_t>& indices, 
                              std::vector<Vector3D>& grads,
                              MsqError& err )
{
  evaluate_with_indices( pd, h, v, indices, err );
  grads.clear();
  for (unsigned i = 0; i < indices.size(); ++i)
    grads.push_back( Vector3D( (double)indices[i], 0, 1 ) );
  return true;
}

bool FauxMetric::evaluate_with_Hessian( PatchData& pd, size_t h, double& v, 
                              std::vector<size_t>& indices, 
                              std::vector<Vector3D>& grad,
                              std::vector<Matrix3D>& hess,
                              MsqError& err )
{
  evaluate_with_gradient( pd, h, v, indices, grad, err );
  hess.clear();
  for (unsigned r = 0; r < indices.size(); ++r)
    for (unsigned c = r; c < indices.size(); ++c)
      hess.push_back( Matrix3D( indices[r],   indices[c], indices[r]+indices[c],
                                indices[r],          1.0, indices[r]+indices[c],
                              2*indices[r], 2*indices[c],            indices[c] ) );
  return true;
}


void PMeanPMetricTest::test_get_metric_type()
{
  FauxMetric m;
  ElementPMeanP e( 1.0, &m );
  VertexPMeanP v( 1.0, &m );
  CPPUNIT_ASSERT( QualityMetric::ELEMENT_BASED == e.get_metric_type() );
  CPPUNIT_ASSERT( QualityMetric::VERTEX_BASED == v.get_metric_type() );
}

void PMeanPMetricTest::test_get_element_evaluations()
{
  MsqError err;
  FauxMetric m;
  ElementPMeanP e( 1.0, &m );
  std::vector<size_t> handles;
    // test that handles array contains all elements
  e.get_evaluations( pd, handles, false, err );
  std::sort( handles.begin(), handles.end() );
  CPPUNIT_ASSERT_EQUAL( pd.num_elements(), handles.size() );
  for (unsigned i = 1; i < handles.size(); ++i)
    CPPUNIT_ASSERT_EQUAL( handles[i-1]+1, handles[i] );
    // test that handles array contains all elements
  e.get_evaluations( pd, handles, true, err );
  std::sort( handles.begin(), handles.end() );
  CPPUNIT_ASSERT_EQUAL( pd.num_elements(), handles.size() );
  for (unsigned i = 1; i < handles.size(); ++i)
    CPPUNIT_ASSERT_EQUAL( handles[i-1]+1, handles[i] );
}
  

void PMeanPMetricTest::test_get_vertex_evaluations()
{
  MsqError err;
  FauxMetric m;
  VertexPMeanP e( 1.0, &m );
  std::vector<size_t> handles;
    // test that handles array contains all vertices
  e.get_evaluations( pd, handles, false, err );
  std::sort( handles.begin(), handles.end() );
  CPPUNIT_ASSERT_EQUAL( pd.num_nodes(), handles.size() );
  for (unsigned i = 1; i < handles.size(); ++i)
    CPPUNIT_ASSERT_EQUAL( handles[i-1]+1, handles[i] );
    // test that handles array contains all vertices
  e.get_evaluations( pd, handles, true, err );
  std::sort( handles.begin(), handles.end() );
  CPPUNIT_ASSERT_EQUAL( pd.num_nodes(), handles.size() );
  for (unsigned i = 1; i < handles.size(); ++i)
    CPPUNIT_ASSERT_EQUAL( handles[i-1]+1, handles[i] );
}

void PMeanPMetricTest::test_vertex_evaluate()
{
  MsqError err;
  FauxMetric m;
  VertexPMeanP m1( 1.0, &m );
  VertexPMeanP m2( 2.0, &m );
  
    // evaluate average around vertex
  double v1, v2;
  m1.evaluate( pd, 0, v1, err );
  CPPUNIT_ASSERT(!err);
  m2.evaluate( pd, 0, v2, err );
  CPPUNIT_ASSERT(!err);

    // get elements adjacent to vertex
  size_t num_elem;
  const size_t* elems = pd.get_vertex_element_adjacencies( 0, num_elem, err );
  CPPUNIT_ASSERT(!err);
  
    // calculate expected values from underyling metric
  double ev1 = 0.0, ev2 = 0.0;
  for (unsigned i = 0; i < num_elem; ++i) {
    const MsqMeshEntity& elem = pd.element_by_index( elems[i] );
    const size_t* verts = elem.get_vertex_index_array();
    const size_t* end = verts + elem.node_count();
    const size_t* p = std::find( verts, end, (size_t)0 );
    CPPUNIT_ASSERT( p < end );
    size_t h = ElemSampleQM::handle( Sample(0,p - verts), elems[i] );
  
    double v;
    m.evaluate( pd, h, v, err );
    CPPUNIT_ASSERT(!err);
    ev1 += v;
    ev2 += v*v;
  }
  
  ev1 /= num_elem;
  ev2 /= num_elem;
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( ev1, v1, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( ev2, v2, 1e-6 );
}

void PMeanPMetricTest::test_element_evaluate()
{
  MsqError err;
  FauxMetric m;
  ElementPMeanP m1( 1.0, &m );
  ElementPMeanP m2( 2.0, &m );
  
    // evaluate average over element
  double v1, v2;
  m1.evaluate( pd, 0, v1, err );
  CPPUNIT_ASSERT(!err);
  m2.evaluate( pd, 0, v2, err );
  CPPUNIT_ASSERT(!err);
  
    // get sample points within element
  std::vector<size_t> handles;
  m.get_element_evaluations( pd, 0, handles, err );
  CPPUNIT_ASSERT(!err);
  
    // calculate expected values from underyling metric
  double ev1 = 0.0, ev2 = 0.0;
  for (unsigned i = 0; i < handles.size(); ++i) {
    double v;
    m.evaluate( pd, handles[i], v, err );
    CPPUNIT_ASSERT(!err);
    ev1 += v;
    ev2 += v*v;
  }
  ev1 /= handles.size();
  ev2 /= handles.size();
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( ev1, v1, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( ev2, v2, 1e-6 );
}


void PMeanPMetricTest::test_indices() 
{
  MsqError err;
  FauxMetric m;
  ElementPMeanP m2( 2.0, &m );
  
  double v1, v2;
  std::vector<size_t> indices;
  m2.evaluate_with_indices( pd, 0, v1, indices, err );
  CPPUNIT_ASSERT(!err);
  m2.evaluate( pd, 0, v2, err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v2, v1, 1e-6 );
  
  std::vector<size_t> vertices;
  pd.element_by_index(0).get_vertex_indices( vertices );
  
  std::sort( vertices.begin(), vertices.end() );
  std::sort( indices.begin(), indices.end() );
  CPPUNIT_ASSERT( vertices == indices );
}

template <typename T>
size_t index_of( const std::vector<T>& v, T a )
{
  return std::find( v.begin(), v.end(), a ) - v.begin();
}

void PMeanPMetricTest::test_gradient()
{
  MsqError err;
  FauxMetric m;
  ElementPMeanP m1( 1.0, &m );
  ElementPMeanP m2( 2.0, &m );
  
    // get vertices for later
  std::vector<size_t> vertices;
  pd.element_by_index(0).get_vertex_indices( vertices );
  
    // evaluate without gradients
  double v1, v2, v3, v4;
  std::vector<size_t> indices1, indices2, indices3, indices4;
  std::vector<Vector3D> grads1, grads2;
  m1.evaluate_with_indices( pd, 0, v1, indices1, err );
  CPPUNIT_ASSERT(!err);
  m2.evaluate_with_indices( pd, 0, v2, indices2, err );
  CPPUNIT_ASSERT(!err);
  
    // evaluate with gradients
  m1.evaluate_with_gradient( pd, 0, v3, indices3, grads1, err );
  CPPUNIT_ASSERT(!err);
  m2.evaluate_with_gradient( pd, 0, v4, indices4, grads2, err );
  CPPUNIT_ASSERT(!err);
  
    // compare value and indices to eval w/out gradient
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v1, v3, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v2, v4, 1e-6 );
  std::sort( indices1.begin(), indices1.end() );
  std::vector<size_t> tm( indices3 );
  std::sort( tm.begin(), tm.end() );
  CPPUNIT_ASSERT( tm == indices1 );
  std::sort( indices2.begin(), indices2.end() );
  tm = indices4 ;
  std::sort( tm.begin(), tm.end() );
  CPPUNIT_ASSERT( tm == indices2 );
  
    // setup evaluation of underying metric
  std::vector<size_t> handles;
  m.get_element_evaluations( pd, 0, handles, err );
  CPPUNIT_ASSERT(!err);
  
    // calculate expected gradients
  std::vector<Vector3D> expected1, expected2, temp;
  expected1.resize( vertices.size(), Vector3D(0,0,0) );
  expected2.resize( vertices.size(), Vector3D(0,0,0) );
  for (unsigned i = 0; i < handles.size(); ++i) {
    double v;
    m.evaluate_with_gradient( pd, handles[i], v, indices3, temp, err );
    CPPUNIT_ASSERT(!err);
    for (unsigned k = 0; k < indices3.size(); ++k) {
      unsigned idx = index_of( vertices, indices3[k] );
      CPPUNIT_ASSERT( idx < vertices.size() );
      expected1[idx] += temp[k];
      expected2[idx] += 2 * v * temp[k];
    }
  }
  for (unsigned i = 0; i < vertices.size(); ++i) {
    expected1[i] /= handles.size();
    expected2[i] /= handles.size();
  }
  
    // compare gradients
  for (unsigned i = 0; i < indices1.size(); ++i) {
    unsigned k = index_of( vertices, indices1[i] );
    CPPUNIT_ASSERT( k < vertices.size() );
    CPPUNIT_ASSERT_VECTORS_EQUAL( expected1[k], grads1[i], 1e-6 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( expected2[k], grads2[i], 1e-6 );
  }
}

void PMeanPMetricTest::test_hessian()
{
  MsqError err;
  FauxMetric m;
  ElementPMeanP m1( 1.0, &m );
  ElementPMeanP m2( 2.0, &m );
  
    // get vertices for later
  std::vector<size_t> vertices;
  pd.element_by_index(0).get_vertex_indices( vertices );
  
    // evaluate gradient
  double v1, v2, v3, v4;
  std::vector<size_t> indices1, indices2, indices3, indices4, tmpi;
  std::vector<Vector3D> grad1, grad2, grad3, grad4;
  std::vector<Matrix3D> hess3, hess4;
  m1.evaluate_with_gradient( pd, 0, v1, indices1, grad1, err );
  CPPUNIT_ASSERT(!err);
  m2.evaluate_with_gradient( pd, 0, v2, indices2, grad2, err );
  CPPUNIT_ASSERT(!err);
  
    // evaluate with Hessian
  m1.evaluate_with_Hessian( pd, 0, v3, indices3, grad3, hess3, err );
  CPPUNIT_ASSERT(!err);
  m2.evaluate_with_Hessian( pd, 0, v4, indices4, grad4, hess4, err );
  CPPUNIT_ASSERT(!err);
  
    // compare value and indices to eval w/out gradient
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v1, v3, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v2, v4, 1e-6 );
    // It isn't a requirement that the index order remain the same
    // for both eval_with_grad and eval_with_Hess, but assuming it 
    // does simplifies a lot of stuff in this test.  Check that the 
    // assumption remains valid.
  CPPUNIT_ASSERT( indices3 == indices1 );
  CPPUNIT_ASSERT( indices4 == indices2 );
    // It isn't a requirement that the index order remain the same
    // for any value of P, but assuming it does simplifies a lot
    // of stuff in this test, so check that the assumption is valid.
  CPPUNIT_ASSERT( indices1 == indices2 ); 
  
    // check that gradient values match
  for (size_t i = 0; i < indices1.size(); ++i) {
    CPPUNIT_ASSERT_VECTORS_EQUAL( grad1[i], grad3[i], 1e-5 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( grad2[i], grad4[i], 1e-5 );
  }
  
    // setup evaluation of underying metric
  std::vector<size_t> handles;
  m.get_element_evaluations( pd, 0, handles, err );
  CPPUNIT_ASSERT(!err);
  
    // calculate expected Hessians
  std::vector<Vector3D> g;
  std::vector<Matrix3D> expected1, expected2, h;
  std::vector<Matrix3D>::iterator h_iter;
  const unsigned N = vertices.size();
  expected1.resize( N*(N+1)/2, Matrix3D(0,0,0,0,0,0,0,0,0) );
  expected2 = expected1;
  Matrix3D outer;
  for (unsigned i = 0; i < handles.size(); ++i) {
    double v;
    m.evaluate_with_Hessian( pd, handles[i], v, tmpi, g, h, err );
    CPPUNIT_ASSERT(!err);
    h_iter = h.begin();
    for (unsigned r = 0; r < tmpi.size(); ++r) {
      unsigned R = index_of( vertices, tmpi[r] );
      CPPUNIT_ASSERT( R < N );
      for (unsigned c = r; c < tmpi.size(); ++c ,++h_iter) {
        unsigned C = index_of( vertices, tmpi[c] );
        CPPUNIT_ASSERT( C < N );
        if (R <= C) {
          unsigned idx = N*R - R*(R+1)/2 + C;
          expected1[idx] += 1.0 / handles.size() * *h_iter;
          expected2[idx] += 2.0 * v / handles.size() * *h_iter;
          outer.outer_product( g[r], g[c] );
          expected2[idx] += 2.0 / handles.size() * outer;
        }
        else {
          unsigned idx = N*C - C*(C+1)/2 + R;
          expected1[idx] += 1.0 / handles.size() * transpose(*h_iter);
          expected2[idx] += 2.0 * v / handles.size() * transpose(*h_iter);
          outer.outer_product( g[c], g[r] );
          expected2[idx] += 2.0 / handles.size() * outer;
        }
      }
    }
  }
  
    // compare Hessians
  unsigned H_idx = 0;
  for (unsigned R = 0; R < vertices.size(); ++R) {
    if (vertices[R] >= pd.num_free_vertices())
      continue;
    unsigned r = index_of( indices3, vertices[R] );
    CPPUNIT_ASSERT(r < indices3.size() );
    for (unsigned C = R; C < vertices.size(); ++C, ++H_idx) {
      if (vertices[C] >= pd.num_free_vertices())
        continue;
      unsigned c = index_of( indices3, vertices[C] );
      CPPUNIT_ASSERT( c < indices3.size() );
      if (r <= c) {
        unsigned idx = indices3.size()*r - r*(r+1)/2 + c;
        CPPUNIT_ASSERT_MATRICES_EQUAL( expected1[H_idx], hess3[idx], 1e-5 );
        CPPUNIT_ASSERT_MATRICES_EQUAL( expected2[H_idx], hess4[idx], 1e-5 );
      }
      else {
        unsigned idx = indices3.size()*c - c*(c+1)/2 + r;
        CPPUNIT_ASSERT_MATRICES_EQUAL( transpose(expected1[H_idx]), hess3[idx], 1e-5 );
        CPPUNIT_ASSERT_MATRICES_EQUAL( transpose(expected2[H_idx]), hess4[idx], 1e-5 );
      }
    }
  }
}

void PMeanPMetricTest::test_hessian_diagonal()
{
  MsqError err;
  FauxMetric m;
  ElementPMeanP m1( 1.0, &m );
  ElementPMeanP m2( 2.0, &m );
  
    // we've already validated the Hessian results in the 
    // previous test, so just check that the diagonal terms
    // match the terms for the full Hessian.
  std::vector<size_t> m1_indices_h, m1_indices_d, m2_indices_h, m2_indices_d;
  std::vector<Vector3D> m1_g_h, m1_g_d, m2_g_h, m2_g_d;
  std::vector<Matrix3D> m1_h_h, m2_h_h;
  std::vector<SymMatrix3D> m1_h_d, m2_h_d;
  double m1_v_h, m1_v_d, m2_v_h, m2_v_d;
  m1.evaluate_with_Hessian( pd, 0, m1_v_h, m1_indices_h, m1_g_h, m1_h_h, err );
  ASSERT_NO_ERROR(err);
  m2.evaluate_with_Hessian( pd, 0, m2_v_h, m2_indices_h, m2_g_h, m2_h_h, err );
  ASSERT_NO_ERROR(err);
  m1.evaluate_with_Hessian_diagonal( pd, 0, m1_v_d, m1_indices_d, m1_g_d, m1_h_d, err );
  ASSERT_NO_ERROR(err);
  m2.evaluate_with_Hessian_diagonal( pd, 0, m2_v_d, m2_indices_d, m2_g_d, m2_h_d, err );
  ASSERT_NO_ERROR(err);
  
    // compare values
  CPPUNIT_ASSERT_DOUBLES_EQUAL( m1_v_h, m1_v_d, 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( m2_v_h, m2_v_d, 1e-6 );
  
    // Assume indices in same order because
    // it simplifiers later code in test
  CPPUNIT_ASSERT( m1_indices_h == m1_indices_d );
  CPPUNIT_ASSERT( m2_indices_h == m1_indices_d );
  
    // compare gradient values
  CPPUNIT_ASSERT_EQUAL(m1_indices_h.size(), m1_g_h.size() );
  CPPUNIT_ASSERT_EQUAL(m2_indices_h.size(), m2_g_h.size() );
  CPPUNIT_ASSERT_EQUAL(m1_indices_d.size(), m1_g_d.size() );
  CPPUNIT_ASSERT_EQUAL(m2_indices_d.size(), m2_g_d.size() );
  for (unsigned i = 0; i < m1_indices_h.size(); ++i) {
    CPPUNIT_ASSERT_VECTORS_EQUAL( m1_g_h[i], m1_g_d[i], 1e-6 );
    CPPUNIT_ASSERT_VECTORS_EQUAL( m2_g_h[i], m2_g_d[i], 1e-6 );
  }
  
    // compare hessian diagonal terms
  CPPUNIT_ASSERT_EQUAL(m1_indices_h.size()*(m1_indices_h.size()+1)/2, m1_h_h.size() );
  CPPUNIT_ASSERT_EQUAL(m2_indices_h.size()*(m2_indices_h.size()+1)/2, m2_h_h.size() );
  CPPUNIT_ASSERT_EQUAL(m1_indices_d.size(), m1_h_d.size() );
  CPPUNIT_ASSERT_EQUAL(m2_indices_d.size(), m2_h_d.size() );
  unsigned h = 0;
  for (unsigned r = 0; r < m1_indices_h.size(); ++r) {
    CPPUNIT_ASSERT_MATRICES_EQUAL( m1_h_h[h], m1_h_d[r], 1e-6 );
    CPPUNIT_ASSERT_MATRICES_EQUAL( m2_h_h[h], m2_h_d[r], 1e-6 );
    h += (m1_indices_h.size() - r);
  }
}
