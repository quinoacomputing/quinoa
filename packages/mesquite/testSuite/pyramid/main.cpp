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
    kraftche@cae.wisc.edu   
   
  ***************************************************************** */


#define TOL 1e-5

#include "meshfiles.h"

#include <iostream>
using std::cout;
using std::endl;
#include <stdlib.h>

#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "LPtoPTemplate.hpp"
#include "LInfTemplate.hpp"

// algorithms
#include "IdealWeightMeanRatio.hpp"
#include "IdealWeightInverseMeanRatio.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
#include "TrustRegion.hpp"
#include "ConditionNumberQualityMetric.hpp"
using namespace Mesquite;

// Use CPPUNIT_ASSERT in code so it's easy to convert to a unit test later.
#define CPPUNIT_ASSERT(A) \
  do { if (!(A)) { \
  std::cout << "Assertion Failed: " << #A << std::endl; \
  std::cout << "  File: " << __FILE__ << std::endl; \
  std::cout << "  Line: " << __LINE__ << std::endl; \
  return true; \
  } } while (false)


// Given a mesh with a single free vertex located at the origin,
// move the vertex to the specified position, smooth the mesh,
// and verify that the vertex was moved back to the origin by
// the smoother.
bool smooth_mesh( Mesh* mesh, Mesh* ref_mesh,
                  Mesh::VertexHandle free_vertex_at_origin, 
                  Vector3D initial_free_vertex_position,
                  QualityMetric* metric );
                  
bool smooth_mixed_mesh( const char* filename );

int main( int argc, char* argv[] )
{
  unsigned i;
  const char* input_file = MESH_FILES_DIR "3D/vtk/mixed/tangled/mixed-hex-pyr-tet.vtk";
  if (argc == 2)
    input_file = argv[1];
  else if (argc != 1)
  {
    std::cerr << "Invalid arguments.\n";
    return 2;
  }
  
  
  Mesquite::MsqPrintError err(cout);
  IdealWeightMeanRatio m1;
  IdealWeightInverseMeanRatio m2(err);
  ConditionNumberQualityMetric m3;
  QualityMetric* metrics[] = { &m1, &m2, &m3, 0 };

    // Read Mesh
  Mesquite::MeshImpl mesh;
  mesh.read_vtk(MESH_FILES_DIR "3D/vtk/pyramids/untangled/12-pyramid-unit-sphere.vtk", err);
  CPPUNIT_ASSERT(!err);
  Mesquite::MeshImpl ideal_mesh;
  ideal_mesh.read_vtk(MESH_FILES_DIR "3D/vtk/pyramids/untangled/12-pyramid-unit-sphere.vtk", err);
  CPPUNIT_ASSERT(!err);

    // Check that the mesh read correctly, and contains what is
    // expected later.

    // Get mesh data
    // Expecting file to contain 12 pyramid elements constructed
    // from 15 vertices.
  std::vector<Mesh::VertexHandle> vert_array;
  std::vector<Mesh::ElementHandle> elem_array;
  std::vector<size_t> conn_offsets;
  mesh.get_all_elements( elem_array, err ); 
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT( elem_array.size() == 12 );
  mesh.elements_get_attached_vertices( arrptr(elem_array),
                                        elem_array.size(),
                                        vert_array,
                                        conn_offsets,
                                        err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT(vert_array.size() == 60);
  CPPUNIT_ASSERT(conn_offsets.size() == 13);
  EntityTopology type_array[12];
  mesh.elements_get_topologies( arrptr(elem_array), type_array, 12, err );
  CPPUNIT_ASSERT(!err);
  
    // Verify element types and number of vertices
  for (i = 0; i < 12; ++i)
  {
    CPPUNIT_ASSERT( type_array[i] == PYRAMID );
    CPPUNIT_ASSERT( conn_offsets[i] == 5*i );
  }
  
    // All pyramids should share a common apex, at the
    // center of the sphere
  Mesh::VertexHandle apex_handle = vert_array[4];
  for (i = 1; i < 12; ++i)
  {
    CPPUNIT_ASSERT( vert_array[5*i+4] == apex_handle );
  }
  
    // Verify that apex is at origin and all other vertices are
    // on unit sphere
  MsqVertex vertices[60];
  mesh.vertices_get_coordinates( arrptr(vert_array), vertices, 60, err );
  CPPUNIT_ASSERT(!err);
  for (i = 0; i < 60; ++i)
  {
    if (vert_array[i] == apex_handle)
      CPPUNIT_ASSERT( vertices[i].within_tolerance_box( Vector3D(0,0,0), 1e-6 ) );
    else
      CPPUNIT_ASSERT( fabs(1.0 - vertices[i].length()) < 1e-6 );
  }
  
    // Try smoothing w/out moving the free vertex and verify that
    // the smoother didn't move the vertex
  Vector3D position(0,0,0);
  for (i = 0; metrics[i] != NULL; ++i)
    CPPUNIT_ASSERT( !smooth_mesh( &mesh, &ideal_mesh, apex_handle, position, metrics[i] ) );
  
    // Now try moving the vertex and see if the smoother moves it back
    // to the origin
  position.set( 0.1, 0.1, 0.1 );
  for (i = 0; metrics[i] != NULL; ++i)
    CPPUNIT_ASSERT( !smooth_mesh( &mesh, &ideal_mesh, apex_handle, position, metrics[i] ) );
  
    // Now try moving the vertex further and see if the smoother moves it back
    // to the origin
  position.set( 0.3, 0.3, 0.3 );
  for (i = 0; metrics[i] != NULL; ++i)
    CPPUNIT_ASSERT( !smooth_mesh( &mesh, &ideal_mesh, apex_handle, position, metrics[i] ) );

    // Now try smoothing a real mixed mesh
  CPPUNIT_ASSERT( !smooth_mixed_mesh( input_file ) );

  return 0;
}
  
  
bool smooth_mesh( Mesh* mesh, Mesh* ref_mesh,
                  Mesh::VertexHandle free_vertex_at_origin, 
                  Vector3D initial_free_vertex_position,
                  QualityMetric* metric )
{
  Mesquite::MsqPrintError err(cout);
  const Vector3D origin( 0, 0, 0 );
  
  // print a little output so we know when we died
  std::cout << 
  "**************************************************************************" 
  << std::endl << 
  "* Smoothing..."
  << std::endl << 
  "* Metric: " << metric->get_name()
  << std::endl << 
  "* Apex position: " << initial_free_vertex_position
  << std::endl //<< 
  //"**************************************************************************" 
  << std::endl;

  
  // Set free vertex to specified position
  mesh->vertex_set_coordinates( free_vertex_at_origin, 
                                initial_free_vertex_position,
                                err );
  CPPUNIT_ASSERT(!err);

  // Create an InstructionQueue
  InstructionQueue Q;

  // Set up objective function
  LPtoPTemplate obj_func(metric, 1, err);
  CPPUNIT_ASSERT(!err);

  // Create solver
  FeasibleNewton solver( &obj_func );
  CPPUNIT_ASSERT(!err);
  solver.use_global_patch();
  CPPUNIT_ASSERT(!err);

  // Set stoping criteria for solver
  TerminationCriterion tc_inner;
  tc_inner.add_absolute_vertex_movement( 1e-6 );
  solver.set_inner_termination_criterion(&tc_inner);
   
  TerminationCriterion tc_outer;
  tc_outer.add_iteration_limit( 1 );
  solver.set_outer_termination_criterion(&tc_outer);
   
  // Add solver to queue
  Q.set_master_quality_improver(&solver, err); 
  CPPUNIT_ASSERT(!err);
 
  // And smooth...
  Q.run_instructions(mesh, err); 
  CPPUNIT_ASSERT(!err);
  
  // Verify that vertex was moved back to origin
  MsqVertex vtx;
  mesh->vertices_get_coordinates( &free_vertex_at_origin, &vtx, 1, err );
  CPPUNIT_ASSERT( !err );
  Vector3D position = vtx;
  
  // print a little output so we know when we died
  std::cout //<< 
  //"**************************************************************************" 
  << std::endl << 
  "* Done Smoothing:"
  << std::endl << 
  "* Metric: " << metric->get_name()
  << std::endl << 
  "* Apex position: " << position
  << std::endl <<  
  "**************************************************************************" 
  << std::endl;
  
  CPPUNIT_ASSERT( position.within_tolerance_box( Vector3D(0,0,0), TOL ) );
  return false;
}



  
bool smooth_mixed_mesh( const char* filename )
{
  Mesquite::MsqPrintError err(cout);
  
  // print a little output so we know when we died
  std::cout << 
  "**************************************************************************" 
  << std::endl << 
  "* Smoothing: " << filename
  << std::endl  <<
  "**************************************************************************" 
  << std::endl;
  
  // The instruction queue to set up
  InstructionQueue Q;
  
  // Use numeric approx of derivitives until analytic solutions
  // are working for pyramids
  IdealWeightInverseMeanRatio mr_metric(err);
  //sRI_DFT dft_metric;
  UntangleBetaQualityMetric un_metric(0);
  CPPUNIT_ASSERT(!err);
  
    // Create Mesh object
  Mesquite::MeshImpl mesh;
  mesh.read_vtk(filename, err);
  CPPUNIT_ASSERT(!err);

  // Set up a preconditioner
  LInfTemplate pre_obj_func( &un_metric );
  ConjugateGradient precond( &pre_obj_func, err ); CPPUNIT_ASSERT(!err);
  precond.use_element_on_vertex_patch();
  TerminationCriterion pre_term, pre_outer;
  //pre_term.add_relative_quality_improvement( 0.1 );
  pre_term .add_iteration_limit( 3 );
  pre_outer.add_iteration_limit( 1 );
  CPPUNIT_ASSERT(!err);
  precond.set_inner_termination_criterion( &pre_term );
  precond.set_outer_termination_criterion( &pre_outer );
  //precond.use_element_on_vertex_patch();

  // Set up objective function
  LPtoPTemplate obj_func(&mr_metric, 1, err);
  CPPUNIT_ASSERT(!err);

  // Create solver
  FeasibleNewton solver( &obj_func );
  CPPUNIT_ASSERT(!err);
  solver.use_global_patch();
  CPPUNIT_ASSERT(!err);

  // Set stoping criteria for solver
  TerminationCriterion tc_inner;
  tc_inner.add_relative_quality_improvement( 0.25 );
  solver.set_inner_termination_criterion(&tc_inner);
   
  TerminationCriterion tc_outer;
  tc_outer.add_iteration_limit( 1 );
  CPPUNIT_ASSERT(!err);
  solver.set_outer_termination_criterion(&tc_outer);

  // Create a QualityAssessor
  Mesquite::QualityAssessor qa;
  qa.add_quality_assessment( &mr_metric );
  qa.add_quality_assessment( &un_metric );
  Q.add_quality_assessor( &qa, err ); 
  CPPUNIT_ASSERT(!err);
 
  // Add untangler to queue
  Q.add_preconditioner( &precond, err ); CPPUNIT_ASSERT(!err);
  Q.add_quality_assessor( &qa, err ); 
  CPPUNIT_ASSERT(!err);
 
  // Add solver to queue
  Q.set_master_quality_improver(&solver, err); 
  CPPUNIT_ASSERT(!err);
  Q.add_quality_assessor( &qa, err ); 
  CPPUNIT_ASSERT(!err);
 
  // And smooth...
  Q.run_instructions(&mesh, err); 
  CPPUNIT_ASSERT(!err);

  return false;
}

