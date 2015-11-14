/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */

/*! \file TerminationCriterionTest.cpp
    \author kraftche@cae.wisc.edu

Tests for the TerminationCriterion class.. 

*/


#include "Mesquite.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_Vector3D.hpp"
#include "Mesquite_OFEvaluator.hpp"
#include "Mesquite_ObjectiveFunction.hpp"
#include "Mesquite_PatchData.hpp"
#include "Mesquite_TerminationCriterion.hpp"
#include "Mesquite_VertexMover.hpp"
#include "Mesquite_VertexPatches.hpp"

#include "UnitUtil.hpp"
#include "PatchDataInstances.hpp"
#include "Mesquite_ArrayMesh.hpp"
#include "Mesquite_MeshUtil.hpp"
#include "Mesquite_SimpleStats.hpp"
#include "Mesquite_InstructionQueue.hpp"

#include <iostream>
#include <algorithm>
#include <vector>
#include <set>

using namespace Mesquite;


class DummyOF : public ObjectiveFunction
{
public:
  double mValue; //!< Objectve fuction value returned
  Vector3D mGrad; //!< Gradient values for all vertices
  bool mValid;

  DummyOF( double of_value = 0.0,
           Vector3D grad_values = Vector3D(0,0,0) )
    : mValue(0.0), mGrad(grad_values), mValid(true) {}

  bool initialize_block_coordinate_descent( MeshDomainAssoc* domain,
                                            const Settings* settings,
                                            PatchSet* user_set,
                                            MsqError& err );
    
  void initialize_queue( MeshDomainAssoc* , 
                         const Settings* ,
                         MsqError&  ) {}

  bool evaluate( EvalType type,
                 PatchData& pd,
                 double& value_out,
                 bool free,
                 MsqError& err );

  bool evaluate_with_gradient( EvalType type,
                               PatchData& pd,
                               double& value_out,
                               std::vector<Vector3D>& grad_out,
                               MsqError& err );
                               
  ObjectiveFunction* clone() const { return new DummyOF(*this); }
  void clear() {}
  int min_patch_layers() const { return 1; }
};

class TerminationCriterionTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(TerminationCriterionTest);

  CPPUNIT_TEST (test_number_of_iterates_inner);
  CPPUNIT_TEST (test_number_of_iterates_outer);
  CPPUNIT_TEST (test_cpu_time_inner);
  CPPUNIT_TEST (test_cpu_time_outer);

  CPPUNIT_TEST (test_absolute_vertex_movement);
  CPPUNIT_TEST (test_relative_vertex_movement);
  CPPUNIT_TEST (test_absolute_vertex_movement_edge_length);

  CPPUNIT_TEST (test_gradient_L2_norm_absolute);
  CPPUNIT_TEST (test_gradient_Linf_norm_absolute);
  CPPUNIT_TEST (test_gradient_L2_norm_relative);
  CPPUNIT_TEST (test_gradient_Linf_norm_relative);

  CPPUNIT_TEST (test_quality_improvement_absolute);
  CPPUNIT_TEST (test_quality_improvement_relative);
  CPPUNIT_TEST (test_successive_improvements_absolute);
  CPPUNIT_TEST (test_successive_improvements_relative);

  CPPUNIT_TEST (test_vertex_bound);
  CPPUNIT_TEST (test_untangled_mesh);
  CPPUNIT_TEST (test_abs_vtx_movement_culling);

  CPPUNIT_TEST_SUITE_END();
  
  DummyOF objFunc;
  OFEvaluator ofEval;
  
  void test_gradient_common( bool absolute, bool L2 );
  void test_quality_common( bool absolute, bool successive );
  void test_vertex_movement_common( bool absolute );
  void test_cpu_time_common( bool inner );
  
public:

  TerminationCriterionTest()
    : ofEval( &objFunc ) {}
  
    //NUMBER OF ITERATES
  void test_number_of_iterates_inner();
  void test_number_of_iterates_outer();
  
    //CPU TIME
  void test_cpu_time_inner()
    { test_cpu_time_common( true ); }
  void test_cpu_time_outer()
    { test_cpu_time_common( false ); }
  
    // VERTEX MOVEMENT
  void test_absolute_vertex_movement()
    { test_vertex_movement_common( true ); }
  void test_relative_vertex_movement()
    { test_vertex_movement_common( false ); }
  void test_absolute_vertex_movement_edge_length();
  
    //GRADIENT NORM ABSOLUTE
  void test_gradient_L2_norm_absolute()
    {  test_gradient_common( true, true ); }
  void test_gradient_Linf_norm_absolute()
    {  test_gradient_common( true, false ); }
  
    //GRADIENT NORM RELATIVE
  void test_gradient_L2_norm_relative()
    {  test_gradient_common( false, true ); }
  void test_gradient_Linf_norm_relative()
    {  test_gradient_common( false, false ); }
  
    //QUALITY IMPROVEMENT ABSOLUTE
  void test_quality_improvement_absolute()
    { test_quality_common( true, false ); }

    //QUALITY IMPROVEMENT RELATIVE
  void test_quality_improvement_relative()
    { test_quality_common( false, false ); }

    //SUCCESSIVE IMPROVEMENTS ABSOLUTE
  void test_successive_improvements_absolute()
    { test_quality_common( true, true ); }
  
    //SUCCESSIVE IMPROVEMENTS RELATIVE
  void test_successive_improvements_relative()
    { test_quality_common( false, true ); }
  
    //VERTEX BOUND
  void test_vertex_bound();
  
  void test_untangled_mesh();
  
    // test culling on ABSOLUTE_VERTEX_MOVEMENT
  void test_abs_vtx_movement_culling();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TerminationCriterionTest, "TerminationCriterionTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TerminationCriterionTest, "Unit");

  
    //NUMBER OF ITERATES
void TerminationCriterionTest::test_number_of_iterates_inner()
{
  MsqPrintError err(std::cout);
  PatchData pd;
  const int LIMIT = 2;

  TerminationCriterion tc;
  tc.add_iteration_limit(LIMIT);
  tc.reset_inner( pd, ofEval, err );
  ASSERT_NO_ERROR(err);
  tc.reset_patch( pd, err );
  ASSERT_NO_ERROR(err);
  for (int i = 0; i < LIMIT; ++i) {
    CPPUNIT_ASSERT(!tc.terminate());
    CPPUNIT_ASSERT_EQUAL( i, tc.get_iteration_count() );
    tc.accumulate_patch( pd, err );
    ASSERT_NO_ERROR(err);
    tc.accumulate_inner( pd, 0, 0, err );
    ASSERT_NO_ERROR(err);
  }
  CPPUNIT_ASSERT_EQUAL( 2, tc.get_iteration_count() );
  CPPUNIT_ASSERT(tc.terminate());
}

void TerminationCriterionTest::test_number_of_iterates_outer()
{
  MsqPrintError err(std::cout);
  PatchData pd;
  const int LIMIT = 2;

  TerminationCriterion tc;
  tc.add_iteration_limit(LIMIT);
  tc.reset_outer( 0, 0, ofEval, 0, err );
  ASSERT_NO_ERROR(err);
  tc.reset_patch( pd, err );
  ASSERT_NO_ERROR(err);
  for (int i = 0; i < LIMIT; ++i) {
    CPPUNIT_ASSERT(!tc.terminate());
    CPPUNIT_ASSERT_EQUAL( i, tc.get_iteration_count() );
    tc.accumulate_patch( pd, err );
    ASSERT_NO_ERROR(err);
    tc.accumulate_outer( 0, 0, ofEval, 0, err );
    ASSERT_NO_ERROR(err);
  }
  CPPUNIT_ASSERT_EQUAL( 2, tc.get_iteration_count() );
  CPPUNIT_ASSERT(tc.terminate());
}

  
    //CPU TIME
void TerminationCriterionTest::test_cpu_time_common( bool inner )
{
  MsqPrintError err(std::cout);
  PatchData pd;
  double LIMIT = 0.5;

  TerminationCriterion tc;
  tc.add_cpu_time(LIMIT);
  if (inner)
    tc.reset_inner( pd, ofEval, err );
  else
    tc.reset_outer( 0, 0, ofEval, 0, err );
  ASSERT_NO_ERROR(err);
  tc.reset_patch( pd, err );
  ASSERT_NO_ERROR(err);
  
  Timer timer;
  while (timer.since_birth() < 0.5*LIMIT) {
    CPPUNIT_ASSERT(!tc.terminate());
    tc.accumulate_patch( pd, err );
    ASSERT_NO_ERROR(err);
    if (inner)
      tc.accumulate_inner( pd, 0, 0, err );
    else
      tc.accumulate_outer( 0, 0, ofEval, 0, err );
    ASSERT_NO_ERROR(err);
  }
  
  while (timer.since_birth() < 1.1*LIMIT);

  tc.accumulate_patch( pd, err );
  ASSERT_NO_ERROR(err);
  if (inner)
    tc.accumulate_inner( pd, 0, 0, err );
  else
    tc.accumulate_outer( 0, 0, ofEval, 0, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(tc.terminate());
}


void TerminationCriterionTest::test_vertex_movement_common( bool absolute )
{
  MsqPrintError err(std::cout);
  PatchData pd;
  create_twelve_hex_patch( pd, err );
  ASSERT_NO_ERROR(err);
 
  const double LIMIT = 1e-4;
  TerminationCriterion tc;
  if (absolute)
    tc.add_absolute_vertex_movement( LIMIT );
  else
    tc.add_relative_vertex_movement( LIMIT );

  tc.reset_inner( pd, ofEval, err );
  ASSERT_NO_ERROR(err);
  tc.reset_patch( pd, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!tc.terminate());

  const double FIRST_STEP=10.0;
    // move a vertex by 10 units and check that it did not meet criterion
  pd.move_vertex( Vector3D(FIRST_STEP,0,0), 0, err );
  ASSERT_NO_ERROR(err);
  tc.accumulate_inner( pd, 0.0, 0, err );
  ASSERT_NO_ERROR(err);
  tc.accumulate_patch( pd, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!tc.terminate());
  
  double test_limit = LIMIT;
  if (!absolute)
    test_limit *= FIRST_STEP;
    
  int idx = 0;
  for (double step = FIRST_STEP; step > test_limit; step *= 0.09) {
    idx = (idx + 1) % pd.num_free_vertices();
    pd.move_vertex( Vector3D(step,0,0), idx, err );
    ASSERT_NO_ERROR(err);
    
    tc.accumulate_inner( pd, 0.0, 0, err );
    ASSERT_NO_ERROR(err);
    tc.accumulate_patch( pd, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(!tc.terminate());
  }
  
  idx = (idx + 1) % pd.num_free_vertices();
  pd.move_vertex( Vector3D(0.5*test_limit,0,0), idx, err );
  ASSERT_NO_ERROR(err);

  tc.accumulate_inner( pd, 0.0, 0, err );
  ASSERT_NO_ERROR(err);
  tc.accumulate_patch( pd, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(tc.terminate());
}

void TerminationCriterionTest::test_absolute_vertex_movement_edge_length()
{
  MsqPrintError err(std::cout);
  
  // define two-tet mesh where tets share a face
  double coords[] = {  0, -5, 0,
                       0,  5, 0,
                       1,  0, 0,
                       0,  0, 0,
                       0,  0, 1 };
  const unsigned long conn[] = { 4, 3, 2, 0,
                                 2, 3, 4, 1 };
  int fixed[5] = {0};
  ArrayMesh mesh( 3, 5, coords, fixed, 2, TETRAHEDRON, conn );
  
    // calculate beta 
  const double LIMIT = 1e-4; // desired absolute limit
  MeshUtil tool(&mesh);
  SimpleStats len;
  tool.edge_length_distribution( len, err );
  ASSERT_NO_ERROR(err);
  const double beta = LIMIT / (len.average() - len.standard_deviation());
  
    // initialize termination criterion
  TerminationCriterion tc;
  tc.add_absolute_vertex_movement_edge_length( beta );
  MeshDomainAssoc mesh_and_domain2 = MeshDomainAssoc(&mesh, 0);
  tc.initialize_queue( &mesh_and_domain2, 0, err ); ASSERT_NO_ERROR(err);
  
    // get a patch data
  PatchData pd;
  pd.set_mesh( &mesh );
  pd.fill_global_patch( err ); ASSERT_NO_ERROR(err);

    // test termination criteiorn
  tc.reset_inner( pd, ofEval, err );
  ASSERT_NO_ERROR(err);
  tc.reset_patch( pd, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!tc.terminate());

  const double FIRST_STEP=10.0;
    // move a vertex by 10 units and check that it did not meet criterion
  pd.move_vertex( Vector3D(FIRST_STEP,0,0), 0, err );
  ASSERT_NO_ERROR(err);
  tc.accumulate_inner( pd, 0.0, 0, err );
  ASSERT_NO_ERROR(err);
  tc.accumulate_patch( pd, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!tc.terminate());
  
  double test_limit = LIMIT;
    
  int idx = 0;
  for (double step = FIRST_STEP; step > test_limit; step *= 0.09) {
    idx = (idx + 1) % pd.num_free_vertices();
    pd.move_vertex( Vector3D(step,0,0), idx, err );
    ASSERT_NO_ERROR(err);
    
    tc.accumulate_inner( pd, 0.0, 0, err );
    ASSERT_NO_ERROR(err);
    tc.accumulate_patch( pd, err );
    ASSERT_NO_ERROR(err);
    CPPUNIT_ASSERT(!tc.terminate());
  }
  
  idx = (idx + 1) % pd.num_free_vertices();
  pd.move_vertex( Vector3D(0.5*test_limit,0,0), idx, err );
  ASSERT_NO_ERROR(err);

  tc.accumulate_inner( pd, 0.0, 0, err );
  ASSERT_NO_ERROR(err);
  tc.accumulate_patch( pd, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(tc.terminate());
}

static double lenfunc( const Vector3D* vect, int len )
  { return Mesquite::length( vect, len ); }
static double maxfunc( const Vector3D* vect, int len )
  { return Mesquite::Linf( vect, len ); }

void TerminationCriterionTest::test_gradient_common( bool absolute, bool L2 )
{
  MsqPrintError err(std::cout);
  PatchData pd;
  create_twelve_hex_patch( pd, err );
  ASSERT_NO_ERROR(err);
  
  const double LIMIT = 1e-4;
  TerminationCriterion tc;
  if (absolute) {
    if (L2)
      tc.add_absolute_gradient_L2_norm( LIMIT );
    else
      tc.add_absolute_gradient_inf_norm( LIMIT );
  }
  else {
    if (L2)
      tc.add_relative_gradient_L2_norm( LIMIT );
    else
      tc.add_relative_gradient_inf_norm( LIMIT );
  }
  
  double (*func_ptr)(const Vector3D*, int) = L2 ? &lenfunc : &maxfunc;

  double junk, value = 1; 
  objFunc.mGrad = Vector3D(value,value,value);
  tc.reset_inner( pd, ofEval, err );
  ASSERT_NO_ERROR(err);
  tc.reset_patch( pd, err );
  ASSERT_NO_ERROR(err);
  
  std::vector<Vector3D> grad;
  ofEval.evaluate( pd, junk, grad, err );
  ASSERT_NO_ERROR(err);
  
  double limit = LIMIT;
  if (!absolute) 
    limit *= func_ptr(arrptr(grad),pd.num_free_vertices());

  while (func_ptr(arrptr(grad),pd.num_free_vertices()) > limit) {
    CPPUNIT_ASSERT(!tc.terminate());

    value *= 0.1;
    objFunc.mGrad = Vector3D(value,value,value);
    ofEval.evaluate( pd, junk, grad, err );
    ASSERT_NO_ERROR(err);
 
    tc.accumulate_inner( pd, 0.0, arrptr(grad), err );
    ASSERT_NO_ERROR(err);
    tc.accumulate_patch( pd, err );
    ASSERT_NO_ERROR(err);
  }
  
  CPPUNIT_ASSERT(tc.terminate());
}

static bool limit_absolute_quality( double ,
                                    double ,
                                    double curr_value,
                                    double epsilon )
{ return curr_value <= epsilon; }
static bool limit_relative_quality( double init_value,
                                    double ,
                                    double curr_value,
                                    double epsilon )
{ return curr_value <= epsilon*init_value; }
static bool limit_absolute_sucessive( double ,
                                      double prev_value,
                                      double curr_value,
                                      double epsilon )
{ return (prev_value - curr_value) <= epsilon; }
static bool limit_relative_sucessive( double init_value,
                                      double prev_value,
                                      double curr_value,
                                      double epsilon )
{ return (prev_value - curr_value) <= epsilon*(init_value-curr_value); }

void TerminationCriterionTest::test_quality_common( bool absolute, bool successive )
{
  MsqPrintError err(std::cout);
  PatchData pd;
  ASSERT_NO_ERROR(err);
  
  const double LIMIT = 1e-4;
  TerminationCriterion tc;
  bool (*func_ptr)(double, double, double, double);
  if (absolute) {
    if (successive) {
      tc.add_absolute_successive_improvement( LIMIT );
      func_ptr = &limit_absolute_sucessive;
    }
    else {
      tc.add_absolute_quality_improvement( LIMIT );
      func_ptr = &limit_absolute_quality;
    }
  }
  else {
    if (successive) {
      tc.add_relative_successive_improvement( LIMIT );
      func_ptr = &limit_relative_sucessive;
    }
    else {
      tc.add_relative_quality_improvement( LIMIT );
      func_ptr = &limit_relative_quality;
    }
  }

  const double INIT_VALUE = 10.0;
  objFunc.mValue = INIT_VALUE;
  tc.reset_inner( pd, ofEval, err );
  ASSERT_NO_ERROR(err);
  tc.reset_patch( pd, err );
  ASSERT_NO_ERROR(err);
  double prev = HUGE_VAL;

  while (!func_ptr(INIT_VALUE,prev,objFunc.mValue,LIMIT)) {
    CPPUNIT_ASSERT(!tc.terminate());

    prev = objFunc.mValue;
    objFunc.mValue *= 0.1;
    tc.accumulate_inner( pd, objFunc.mValue, 0, err );
    ASSERT_NO_ERROR(err);
    tc.accumulate_patch( pd, err );
    ASSERT_NO_ERROR(err);
  }
  
  CPPUNIT_ASSERT(tc.terminate());
}

//VERTEX BOUND
void TerminationCriterionTest::test_vertex_bound()
{
  MsqPrintError err(std::cout);
  PatchData pd;
  create_twelve_hex_patch( pd, err );
  ASSERT_NO_ERROR(err);
  
    // get bounding dimension for patch
  double maxcoord = 0.0;
  for (size_t i = 0; i < pd.num_nodes(); ++i) 
    for (int d = 0; d < 3; ++d)
      if (fabs(pd.vertex_by_index(i)[d]) > maxcoord)
        maxcoord = fabs(pd.vertex_by_index(i)[d]);
    // add a little bit for rounding error
  maxcoord += 1e-5;
  
  TerminationCriterion tc;
  tc.add_bounded_vertex_movement( maxcoord );
  tc.reset_inner( pd, ofEval, err );
  ASSERT_NO_ERROR(err);
  tc.reset_patch( pd, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!tc.terminate());
  
  int idx = pd.num_free_vertices() - 1;
  Vector3D pos = pd.vertex_by_index(idx);
  pos[0] = 2*maxcoord;
  pd.set_vertex_coordinates( pos, idx, err );
  ASSERT_NO_ERROR(err);
  tc.accumulate_inner( pd, 0.0, 0, err );
  ASSERT_NO_ERROR(err);
  tc.accumulate_patch( pd, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(tc.terminate());
}



//UNTANGLED
void TerminationCriterionTest::test_untangled_mesh()
{
  MsqPrintError err(std::cout);
  PatchData pd;
  create_twelve_hex_patch( pd, err );
  ASSERT_NO_ERROR(err);
  
    // get two opposite vertices in first hexahedral element
  int vtx1 = pd.element_by_index(0).get_vertex_index_array()[0];
  int vtx2 = pd.element_by_index(0).get_vertex_index_array()[7];
  Vector3D saved_coords = pd.vertex_by_index(vtx2);
  Vector3D opposite_coords = pd.vertex_by_index(vtx1);

    // invert the element
  pd.move_vertex( 2*(opposite_coords-saved_coords), vtx2, err );
  ASSERT_NO_ERROR(err);
  int inverted, samples;
  pd.element_by_index(0).check_element_orientation( pd, inverted, samples, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(inverted > 0);
  
    // check initial termination criterion
  TerminationCriterion tc;
  tc.add_untangled_mesh();
  tc.reset_inner( pd, ofEval, err );
  ASSERT_NO_ERROR(err);
  tc.reset_patch( pd, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(!tc.terminate());
  
    // fix the element
  pd.set_vertex_coordinates( saved_coords, vtx2, err );
  ASSERT_NO_ERROR(err);
  pd.element_by_index(0).check_element_orientation( pd, inverted, samples, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(0,inverted);
  
    // check that TC recognized untangled mesh
  tc.accumulate_inner( pd, 0.0, 0, err );
  ASSERT_NO_ERROR(err);
  tc.accumulate_patch( pd, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT(tc.terminate());
}

class TCTFauxOptimizer : public VertexMover
{
public:
  TCTFauxOptimizer( double pertubation_amount ) : mDelta(pertubation_amount) {}
  virtual ~TCTFauxOptimizer() {};
  virtual std::string get_name() const { return "Optimizer for TerminationCriterionTest"; }
  virtual PatchSet* get_patch_set() { return &mPatchSet; }
  virtual void initialize( PatchData& pd, MsqError& err );
  virtual void initialize_mesh_iteration( PatchData& pd, MsqError& err ) {}
  virtual void optimize_vertex_positions( PatchData& pd, MsqError& err ); 
  virtual void terminate_mesh_iteration( PatchData& pd, MsqError& err ) {}
  virtual void cleanup() { all.clear(); }
  int num_passes() const { return numPasses; }
  bool should_have_terminated() const { return perturbFrac > (int)all.size(); }
private:
  std::set<Mesh::VertexHandle> culled, visited;
  std::vector<Mesh::VertexHandle> all;
  VertexPatches mPatchSet;
  int numPasses;   //!< count number of outer iterations
  int perturbFrac; //!< perturb 1/perturbFrac of the free vertices
  double mDelta;   //!< distance to perturb vertices
};

void TCTFauxOptimizer::initialize( PatchData& pd, MsqError& err )
{
  CPPUNIT_ASSERT(all.empty());
  culled.clear();
  visited.clear();
  numPasses = 1;
  
  pd.get_mesh()->get_all_vertices( all, err ); MSQ_ERRRTN(err);
  std::vector<bool> fixed;
  pd.get_mesh()->vertices_get_fixed_flag( &all[0], fixed, all.size(), err );
  size_t w = 0;
  for (size_t r = 0; r < all.size(); ++r)
    if (!fixed[r])
      all[w++] = all[r];
  all.resize(w);
  MSQ_ERRRTN(err);
  
  perturbFrac = 1;
}

void TCTFauxOptimizer::optimize_vertex_positions( PatchData& pd, MsqError& err )
{
  Mesh::VertexHandle free_vtx = pd.get_vertex_handles_array()[0];
  if (visited.insert(free_vtx).second == false) { // already visited this one
    // The inner termination criterion should include an iteration limit of 1.
    // So if we are seeing the same vertex again, this means that we *should*
    // be stating a new pass over the mesh.
  
    // We are presumably starting a new pass over the mesh.
    // Verify that we visisted all of the free, non-culled vertices
    for (size_t i = 0; i < all.size(); ++i) {
      if (culled.find(all[i]) == culled.end()) {
        if (visited.find(all[i]) == visited.end()) {
          std::ostringstream str;
          str << "Did not visit vertex " << i << " (handle " << all[i] << ") in pass " << numPasses << std::endl;
          CPPUNIT_FAIL( str.str() );
        }
      }
    }
    visited.clear();
    visited.insert(free_vtx);
    ++numPasses;
    
    // Check that we terminate when expected
    CPPUNIT_ASSERT(!should_have_terminated());
    
    perturbFrac *= 2; // for each pass, perturb half as many vertices
  }
  
    // check that we are not visiting a culled vertex
  CPPUNIT_ASSERT( culled.find(free_vtx) == culled.end() );
  
    // for each pass, perturb half as many vertices
  size_t idx = std::find( all.begin(), all.end(), free_vtx ) - all.begin();
  CPPUNIT_ASSERT( idx < all.size() ); // not a free vertex????
  if (0 == ((idx+1) % perturbFrac)) {
      // perturb vertex
    double sign = numPasses % 2 == 0 ? 1 : -1;
    Vector3D delta( sign * mDelta, 0, 0 );
    pd.move_vertex( delta, 0, err ); 
    ASSERT_NO_ERROR(err);
      // any adjacent vertices should not be culled
    for (size_t i = 0; i < pd.num_nodes(); ++i)
      culled.erase( pd.get_vertex_handles_array()[i] );
  }
  else {  
      // If we're not moving this vertex, then it should get culled
    culled.insert( free_vtx );
  }
}

void TerminationCriterionTest::test_abs_vtx_movement_culling()
{
  /* define a quad mesh like the following
    16--17--18--19--20--21--22--23
     |   |   |   |   |   |   |   |   Y
     8---9--10--11--12--13--14--15   ^
     |   |   |   |   |   |   |   |   |
     0---1---2---3---4---5---6---7   +-->X
  */

  const int nvtx = 24;
  int fixed[nvtx];
  double coords[3*nvtx];
  for (int i = 0; i < nvtx; ++i) {
    coords[3*i  ] = i/8;
    coords[3*i+1] = i%8;
    coords[3*i+2] = 0;
    fixed[i] = i < 9 || i > 14;
  }
  
  const int nquad = 14;
  unsigned long conn[4*nquad];
  for (int i= 0; i < nquad; ++i) {
    int row = i / 7;
    int idx = i % 7;
    int n0 = 8*row + idx;
    conn[4*i  ] = n0;
    conn[4*i+1] = n0 + 1;
    conn[4*i+2] = n0 + 9;
    conn[4*i+3] = n0 + 8;
  }
  
  const double tol = 0.05;
  PlanarDomain zplane(PlanarDomain::XY, 0.0);
  ArrayMesh mesh( 3, nvtx, coords, fixed, nquad, QUADRILATERAL, conn );
  
    // fill vertex byte with garbage to make sure that quality improver
    // is initializing correctly.
  MsqError err;
  std::vector<unsigned char> junk(nvtx, 255);
  std::vector<Mesh::VertexHandle> vertices;
  mesh.get_all_vertices( vertices, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_EQUAL( junk.size(), vertices.size() );
  mesh.vertices_set_byte( &vertices[0], &junk[0], vertices.size(), err );
  ASSERT_NO_ERROR(err);  
  
    // Define optimizer
  TCTFauxOptimizer smoother( 2*tol );
  TerminationCriterion outer, inner;
  outer.add_absolute_vertex_movement( tol );
  inner.cull_on_absolute_vertex_movement( tol );
  inner.add_iteration_limit( 1 );
  smoother.set_inner_termination_criterion( &inner );
  smoother.set_outer_termination_criterion( &outer );
  
    // No run the "optimizer"
  InstructionQueue q;
  q.set_master_quality_improver( &smoother, err );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &zplane);
  q.run_instructions( &mesh_and_domain, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT( smoother.should_have_terminated() );
  CPPUNIT_ASSERT( smoother.num_passes() > 1 );
}


bool DummyOF::initialize_block_coordinate_descent( MeshDomainAssoc*,
                                                   const Settings*,
                                                   PatchSet*,
                                                   MsqError& err )
{
  MSQ_SETERR(err)(MsqError::NOT_IMPLEMENTED);
  return false;
}

bool DummyOF::evaluate( EvalType, PatchData&, double& value, bool, MsqError& )
{
  value = mValue;
  return mValid;
}

bool DummyOF::evaluate_with_gradient( EvalType,
                                      PatchData& pd,
                                      double& value_out,
                                      std::vector<Vector3D>& grad_out,
                                      MsqError& )
{
  value_out = mValue;
  grad_out.clear();
  grad_out.resize( pd.num_free_vertices(), mGrad );
  return mValid;
}
