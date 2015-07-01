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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 23-Jul-03 at 18:09:51 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
//#include "QuadLagrangeShape.hpp"

// algorithms
#include "IdealShapeTarget.hpp"
#include "TQualityMetric.hpp"
#include "TInverseMeanRatio.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"

#include "PlanarDomain.hpp"
using namespace Mesquite;

/* This is the input mesh topology
     (0)------(16)-----(1)------(17)-----(2)------(18)-----(3)
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
     (19)      0       (20)      1       (21)      2       (22)
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
     (4)------(23)-----(5)------(24)-----(6)------(25)-----(7)
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
     (26)      3       (27)      4       (28)      5       (29)
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
     (8)------(30)-----(9)------(31)-----(10)-----(32)-----(11)
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
     (33)      6       (34)      7       (35)      8       (36)
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
      |                 |                 |                 |
     (12)-----(37)-----(13)-----(38)-----(14)-----(39)-----(15)
*/


const char LINEAR_INPUT_FILE_NAME[]       = SRCDIR "linear_input.vtk";
const char QUADRATIC_INPUT_FILE_NAME[]    = SRCDIR "quadratic_input.vtk";
const char EXPECTED_LINAR_FILE_NAME[]     = SRCDIR "expected_linear_output.vtk";
const char EXPECTED_QUADRATIC_FILE_NAME[] = SRCDIR "expected_quadratic_output.vtk";
const char HOUR_INPUT_FILE_NAME[]         = SRCDIR "hour-quad8.vtk";
const char OUTPUT_FILE_NAME[]             = "smoothed_qudratic_mesh.vtk";
const unsigned NUM_CORNER_VERTICES = 16;
const unsigned NUM_MID_NODES = 24;
const double SPATIAL_COMPARE_TOLERANCE = 4e-6;


void compare_nodes( size_t start_index,
                    size_t end_index,
                    Mesh* mesh1,
                    Mesh* mesh2,
                    MsqError& err )
{
  size_t i, num_verts = end_index - start_index;
  std::vector<MsqVertex> verts1(num_verts), verts2(num_verts);
  std::vector<Mesh::VertexHandle> handles1(num_verts), handles2(num_verts);

/* VertexIterator skips higher-order nodes.
   For now, just assume index == handle

  std::vector<Mesh::VertexHandle>::iterator handle_iter1, handle_iter2;
  
    // Skip start_index vertices
  VertexIterator* iter1 = mesh1->vertex_iterator( err ); MSQ_ERRRTN(err);
  VertexIterator* iter2 = mesh2->vertex_iterator( err ); MSQ_ERRRTN(err);
  for (i = 0; i < start_index; ++i)
  {
    if (iter1->is_at_end())
    {
      MSQ_SETERR(err)("start index out of range for first mesh set", MsqError::INVALID_ARG);
      return;
    }
    if (iter2->is_at_end())
    {
      MSQ_SETERR(err)("start index out of range for second mesh set", MsqError::INVALID_ARG);
      return;
    }
    iter1->operator++();
    iter2->operator++();
  }
  
    // Get handles for vertices
  handle_iter1 = handles1.begin();
  handle_iter2 = handles2.begin();
  for (i = start_index; i < end_index; ++i)
  {
    if (iter1->is_at_end())
    {
      MSQ_SETERR(err)("end index out of range for first mesh set", MsqError::INVALID_ARG);
      return;
    }
    *handle_iter1 = iter1->operator*();
    iter1->operator++();
    ++handle_iter1;
    
    if (iter2->is_at_end())
    {
      MSQ_SETERR(err)("end index out of range for second mesh set", MsqError::INVALID_ARG);
      return;
    }
    *handle_iter2 = iter2->operator*();
    iter2->operator++();
    ++handle_iter2;
  }
*/
  for (i = start_index; i < end_index; ++i)
    handles1[i-start_index] = handles2[i-start_index] = (void*)i;
  
  
    // Get coordinates from handles
  mesh1->vertices_get_coordinates( arrptr(handles1), arrptr(verts1), num_verts, err );
  MSQ_ERRRTN(err);
  mesh2->vertices_get_coordinates( arrptr(handles2), arrptr(verts2), num_verts, err );
  MSQ_ERRRTN(err);
  
    // Compare coordinates
  for (i = 0; i < num_verts; ++i)
  {
    const double diff = (verts1[i] - verts2[i]).length();
    if (diff > SPATIAL_COMPARE_TOLERANCE)
    {
      MSQ_SETERR(err)(MsqError::INTERNAL_ERROR, 
                      "%u%s vertices differ. (%f,%f,%f) vs (%f,%f,%f)",
                      (unsigned)(1+i),
                      i%10 == 0 ? "st" :
                      i%10 == 1 ? "nd" :
                      i%10 == 2 ? "rd" : "th",
                      verts1[i][0],verts1[i][1],verts1[i][2],
                      verts2[i][0],verts2[i][1],verts2[i][2]);
      return;
    }
  }
}
  
  // code copied from testSuite/algorithm_test/main.cpp
InstructionQueue* create_instruction_queue(MsqError& err)
{
  
    // creates an intruction queue
  InstructionQueue* queue1 = new InstructionQueue;

  // creates a mean ratio quality metric ...
  //IdealWeightInverseMeanRatio* mean = new IdealWeightInverseMeanRatio(err); MSQ_ERRZERO(err);
  TargetCalculator* tc = new IdealShapeTarget;
  TQualityMetric* mean = new TQualityMetric( tc, 0, new TInverseMeanRatio );
  
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean, 1, err); MSQ_ERRZERO(err);
  
  // creates the optimization procedures
//   ConjugateGradient* pass1 = new ConjugateGradient( obj_func, err );
  SteepestDescent* pass1 = new SteepestDescent( obj_func );

  //perform optimization globally
  pass1->use_global_patch();
  
  //QualityAssessor* mean_qa = new QualityAssessor(mean);

    //**************Set termination criterion****************

  //perform 1 pass of the outer loop (this line isn't essential as it is
  //the default behavior).
  TerminationCriterion* tc_outer = new TerminationCriterion;
  tc_outer->add_iteration_limit( 1 );
  pass1->set_outer_termination_criterion(tc_outer);
  
  //perform the inner loop until a certain objective function value is
  //reached.  The exact value needs to be determined (about 18095).
  //As a safety, also stop if the time exceeds 10 minutes (600 seconds).
  TerminationCriterion* tc_inner = new TerminationCriterion;
  tc_inner->add_absolute_vertex_movement( 1e-6 ); 
  pass1->set_inner_termination_criterion(tc_inner);
  
  // adds 1 pass of pass1 to mesh_set1
  //queue1->add_quality_assessor(mean_qa,err); MSQ_ERRZERO(err);
  queue1->set_master_quality_improver(pass1, err); MSQ_ERRZERO(err);
  //queue1->add_quality_assessor(mean_qa,err); MSQ_ERRZERO(err);

  return queue1;
}

int do_test( bool slave)
{
  MsqPrintError err(cout);
//  QuadLagrangeShape quad9;
  
    // Create geometry
  Vector3D z(0,0,1), o(0,0,0);
  PlanarDomain geom(z,o);  
  
    // Read in linear input mesh
  cout << "Reading " << LINEAR_INPUT_FILE_NAME << endl;
  MeshImpl* linear_in = new MeshImpl;
  linear_in->read_vtk( LINEAR_INPUT_FILE_NAME, err );
  if (MSQ_CHKERR(err)) return 1;
  
    // Read in expected linear results
  cout << "Reading " << EXPECTED_LINAR_FILE_NAME << endl;
  MeshImpl* linear_ex = new MeshImpl;
  linear_ex->read_vtk( EXPECTED_LINAR_FILE_NAME, err );
  if (MSQ_CHKERR(err)) return 1;
 
    // Read in second copy of quadratic input mesh
  cout << "Reading " << QUADRATIC_INPUT_FILE_NAME << " again" << endl;
  MeshImpl* quadratic_in_2 = new MeshImpl;
  quadratic_in_2->read_vtk( QUADRATIC_INPUT_FILE_NAME, err );
  if (MSQ_CHKERR(err)) return 1;
  
    // Read in expected quadratic results
  cout << "Reading " << EXPECTED_QUADRATIC_FILE_NAME << endl;
  MeshImpl* quadratic_ex = new MeshImpl;
  quadratic_ex->read_vtk( EXPECTED_QUADRATIC_FILE_NAME, err );
  if (MSQ_CHKERR(err)) return 1;
  

    // Smooth linear mesh and check results
  cout << "Smoothing linear elements" << endl;
  InstructionQueue* q1 = create_instruction_queue( err );
  if (MSQ_CHKERR(err)) return 1;
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(linear_in, &geom);
  q1->run_instructions( &mesh_and_domain, err ); 
  if (MSQ_CHKERR(err)) return 1;
  cout << "Checking results" << endl;
  compare_nodes( 0, NUM_CORNER_VERTICES, linear_in, linear_ex, err );
  if (MSQ_CHKERR(err)) {
    MsqError tmperr;
    linear_in->write_vtk("bad_mesh.vtk", tmperr);
    return 1;
  }
  delete q1;
 
    // Smooth corner vertices and adjust mid-side nodes
  cout << "Smoothing quadratic elements" << endl;
  InstructionQueue* q3 = create_instruction_queue( err );
  if (MSQ_CHKERR(err)) return 1;
  if (!slave)
    q3->set_slaved_ho_node_mode(Settings::SLAVE_NONE);
//  q3->set_mapping_function( &quad9 );
  MeshDomainAssoc mesh_and_domain2 = MeshDomainAssoc(quadratic_in_2, &geom);
  q3->run_instructions( &mesh_and_domain2, err ); 
  if (MSQ_CHKERR(err)) return 1;
    // Make sure corner vertices are the same as in the linear case
  cout << "Checking results" << endl;
  compare_nodes( 0, NUM_CORNER_VERTICES, quadratic_in_2, linear_ex, err );
  if (MSQ_CHKERR(err)) return 1;
    // Make sure mid-side vertices are updated correctly
  compare_nodes( NUM_CORNER_VERTICES, NUM_CORNER_VERTICES + NUM_MID_NODES,
                quadratic_in_2, quadratic_ex, err );
  if (MSQ_CHKERR(err)) {
    MsqError tmperr;
    quadratic_in_2->write_vtk("bad_mesh.vtk", tmperr);
    return 1;
  }
  delete q3;
  
  if (MSQ_CHKERR(err)) 
    return 1;
  return 0;
}

int do_smooth_ho()
{
  MsqPrintError err(cout);
//  QuadLagrangeShape quad9;
  
    // Create geometry
  PlanarDomain geom(PlanarDomain::XY);  
  
    // Read in one copy of quadratic input mesh
  cout << "Reading " << HOUR_INPUT_FILE_NAME << endl;
  MeshImpl* quadratic_in = new MeshImpl;
  quadratic_in->read_vtk( HOUR_INPUT_FILE_NAME, err );
  if (MSQ_CHKERR(err)) return 1;
  
    // Read in expected results
  //cout << "Reading " << HOUR_EXPECTED_FILE_NAME << endl;
  //MeshImpl* quadratic_ex = new MeshImpl;
  //quadratic_ex->read_vtk( results, err );
  //if (MSQ_CHKERR(err)) return 1;

    // Smooth linear mesh and check results
  cout << "Smoothing higher-order nodes" << endl;
  InstructionQueue* q1 = create_instruction_queue( err );
  if (MSQ_CHKERR(err)) return 1;
  q1->set_slaved_ho_node_mode(Settings::SLAVE_NONE);
//  q1->set_mapping_function( &quad9 );
    MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(quadratic_in, &geom);
  q1->run_instructions( &mesh_and_domain, err ); 
  if (MSQ_CHKERR(err)) return 1;
  cout << "Checking results" << endl;
  //compare_nodes( 0, NUM_CORNER_VERTICES + NUM_MID_NODES,
  //               quadratic_in, quadratic_ex, err );
  if (MSQ_CHKERR(err)) {
    MsqError tmperr;
    quadratic_in->write_vtk("bad_mesh.vtk", tmperr);
    return 1;
  }
  delete q1;
  quadratic_in->write_vtk("smooth_ho.vtk", err);
  
  if (MSQ_CHKERR(err)) 
    return 1;
  return 0;
}

int main()
{ 
  cout << "Running test with all higher-order nodes slaved." << endl;    
  int result1 = do_test(true);
  cout << "Running test with no higher-order nodes slaved." << endl;    
  int result2 = do_test(false);
  cout << "Running test with only ho-nodes free." << endl;
  int result3 = do_smooth_ho();
  return result1 + result2 + result3;
}
