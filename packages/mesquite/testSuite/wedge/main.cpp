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

//#define DO_QUALITY_ASSESSOR

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

#define CPPUNIT_ASSERT_DOUBLES_EQUAL(E,V,T) \
  do { if (fabs(E-V) > T) { \
  std::cout << "Assertion Failed: " << #V << " == " << #E << std::endl; \
  std::cout << "Expected: " << E << "  Got: " << V << std::endl;\
  std::cout << "  File: " << __FILE__ << std::endl; \
  std::cout << "  Line: " << __LINE__ << std::endl; \
  return true; \
  } } while (false)

#define CPPUNIT_ASSERT_EQUAL(E,V) \
  do { if (E != V) { \
  std::cout << "Assertion Failed: " << #V << " == " << #E << std::endl; \
  std::cout << "Expected: " << E << "  Got: " << V << std::endl;\
  std::cout << "  File: " << __FILE__ << std::endl; \
  std::cout << "  Line: " << __LINE__ << std::endl; \
  return true; \
  } } while (false)

#define CPPUNIT_ASSERT_VECTORS_EQUAL(E,V,T) \
  do { if (!E.within_tolerance_box(V,T)) { \
  std::cout << "Assertion Failed: " << #V << " == " << #E << std::endl; \
  std::cout << "Expected: " << E << "  Got: " << V << std::endl;\
  std::cout << "  File: " << __FILE__ << std::endl; \
  std::cout << "  Line: " << __LINE__ << std::endl; \
  return true; \
  } } while (false)


// Given a mesh with a single free vertex located at the origin,
// move the vertex to the specified position, smooth the mesh,
// and verify that the vertex was moved back to the origin by
// the smoother.
bool smooth_mesh( MeshImpl* mesh, Mesh* ref_mesh,
                  Mesh::VertexHandle vertex_8,
                  Mesh::VertexHandle vertex_9, 
                  Vector3D delta,
                  QualityMetric* metric );

const unsigned NUM_ELEM = 6;
const unsigned NUM_VERT = 14;
const unsigned VERT_PER_ELEM = 6;

int main( int argc, char* argv[] )
{
  if (argc != 1)
  {
    std::cerr << "Invalid arguments.\n";
    return 2;
  }
  
  const char* meshfile = MESH_FILES_DIR "3D/vtk/prisms/untangled/6-wedge-prism.vtk";
  unsigned i;

  Mesquite::MsqPrintError err(cout);
  IdealWeightMeanRatio m1;
  IdealWeightInverseMeanRatio m2(err);
  ConditionNumberQualityMetric m3;
  QualityMetric* metrics[] = { &m1, &m2, &m3, 0 };

    // Read Mesh
  Mesquite::MeshImpl mesh;
  mesh.read_vtk( meshfile, err);
  CPPUNIT_ASSERT(!err);
  Mesquite::MeshImpl ideal_mesh;
  ideal_mesh.read_vtk( meshfile, err);
  CPPUNIT_ASSERT(!err);

    // Check that the mesh read correctly, and contains what is
    // expected later.

    // Get mesh data
    // Expecting file to contain 6 wedge elements constructed
    // from 14 vertices.
  std::vector<Mesh::VertexHandle> vert_array;
  std::vector<Mesh::ElementHandle> elem_array;
  std::vector<size_t> conn_offsets;
  mesh.get_all_elements( elem_array, err ); 
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( elem_array.size(), NUM_ELEM );
  mesh.elements_get_attached_vertices( arrptr(elem_array),
                                        elem_array.size(),
                                        vert_array,
                                        conn_offsets,
                                        err );
  CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL(vert_array.size() , VERT_PER_ELEM*NUM_ELEM);
  CPPUNIT_ASSERT_EQUAL(conn_offsets.size() , NUM_ELEM+1);
  EntityTopology type_array[NUM_ELEM];
  mesh.elements_get_topologies( arrptr(elem_array), type_array, NUM_ELEM, err );
  CPPUNIT_ASSERT(!err);
  
    // Verify element types and number of vertices
  for (i = 0; i < NUM_ELEM; ++i)
  {
    CPPUNIT_ASSERT_EQUAL( type_array[i] , PRISM );
    CPPUNIT_ASSERT_EQUAL( conn_offsets[i] , VERT_PER_ELEM*i );
  }
  
    // All wedges should share the 9th and 10th vertices
  const unsigned INDEX_8 = 1, INDEX_9 = 4;
  Mesh::VertexHandle handle8 = vert_array[INDEX_8];
  Mesh::VertexHandle handle9 = vert_array[INDEX_9];
  for (i = 1; i < NUM_ELEM; ++i)
  {
    CPPUNIT_ASSERT_EQUAL( vert_array[VERT_PER_ELEM*i+INDEX_8] , handle8 );
    CPPUNIT_ASSERT_EQUAL( vert_array[VERT_PER_ELEM*i+INDEX_9] , handle9 );
  }
  
    // The input file should be a hexagonal prism decomposed into
    // 6 wedges such that all wedges have one quad face on the 
    // boundary and all wedges share a single edge which is the axis
    // of the prism.
  MsqVertex vertices[NUM_ELEM*VERT_PER_ELEM];
  mesh.vertices_get_coordinates( arrptr(vert_array), vertices, NUM_ELEM*VERT_PER_ELEM, err );
  CPPUNIT_ASSERT(!err);
  for (i = 0; i < NUM_ELEM*VERT_PER_ELEM; ++i)
  {
    if (vert_array[i] == handle8 || vert_array[i] == handle9)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, vertices[i][1], 1e-6 ); 
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, vertices[i][2], 1e-6 ); 
    }
    else
    {
      Vector3D xproj( vertices[i] );
      xproj[0] = 0;
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0, xproj.length(), 1e-6 ); 
    }
  }
  
    // Try smoothing w/out moving the free vertices and verify that
    // the smoother didn't move the vertex
  Vector3D delta(0,0,0);
  for (i = 0; metrics[i] != NULL; ++i)
    CPPUNIT_ASSERT( !smooth_mesh( &mesh, &ideal_mesh, handle8, handle9, delta, metrics[i] ) );
  
    // Now try moving the vertex and see if the smoother moves it back
    // to the origin
  delta.set( 0.1, 0.0, 0.1 );
  for (i = 0; metrics[i] != NULL; ++i)
    CPPUNIT_ASSERT( !smooth_mesh( &mesh, &ideal_mesh, handle8, handle9, delta, metrics[i] ) );
  
    // Now try moving the vertex further and see if the smoother moves it back
    // to the origin
  delta.set( 1.0, 0.0, 1.0 );
  for (i = 0; metrics[i] != NULL; ++i)
    CPPUNIT_ASSERT( !smooth_mesh( &mesh, &ideal_mesh, handle8, handle9, delta, metrics[i] ) );

  return 0;
}
  
  
bool smooth_mesh( MeshImpl* mesh, Mesh* ref_mesh,
                  Mesh::VertexHandle vert1, 
                  Mesh::VertexHandle vert2, 
                  Vector3D delta,
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
  "* Offset: " << delta
  << std::endl //<< 
  //"**************************************************************************" 
  << std::endl;
  
  
  // Set free vertices to specified position
  Mesh::VertexHandle handles[] = { vert1, vert2 };
  MsqVertex coordinates[2];
  mesh->vertices_get_coordinates( handles, coordinates, 2, err );
  CPPUNIT_ASSERT(!err);
  Vector3D coord1 = coordinates[0] + delta;
  Vector3D coord2 = coordinates[1] + delta;
  mesh->vertex_set_coordinates( vert1, coord1, err ); CPPUNIT_ASSERT(!err);
  mesh->vertex_set_coordinates( vert2, coord2, err ); CPPUNIT_ASSERT(!err);

  // Create an InstructionQueue
  InstructionQueue Q;

  // Set up objective function
  LPtoPTemplate obj_func(metric, 1, err);
  CPPUNIT_ASSERT(!err);

  // Create solver
  FeasibleNewton solver( &obj_func );
  CPPUNIT_ASSERT(!err);
  solver.use_global_patch();

  // Set stoping criteria for solver
  TerminationCriterion tc_inner;
  tc_inner.add_absolute_vertex_movement( 1e-7 );
  solver.set_inner_termination_criterion(&tc_inner);
   
  TerminationCriterion tc_outer;
  tc_outer.add_iteration_limit( 1 );
  solver.set_outer_termination_criterion(&tc_outer);

#ifdef DO_QUALITY_ASSESSOR
  QualityAssessor qa( metric, 10 );
  Q.add_quality_assessor( &qa, err ); CPPUNIT_ASSERT(!err);
#endif
   
  // Add solver to queue
  Q.set_master_quality_improver(&solver, err); 
  CPPUNIT_ASSERT(!err);

#ifdef DO_QUALITY_ASSESSOR
  Q.add_quality_assessor( &qa, err ); CPPUNIT_ASSERT(!err);
#endif
 
  // And smooth...
  Q.run_instructions(mesh, err); 
  CPPUNIT_ASSERT(!err);
  
  // Verify that vertices were moved back to origin
  MsqVertex new_coords[2];
  mesh->vertices_get_coordinates( handles, new_coords, 2, err );
  CPPUNIT_ASSERT(!err);
  
  // print a little output so we know when we died
  std::cout //<< 
  //"**************************************************************************" 
  << std::endl << 
  "* Done Smoothing:"
  << std::endl << 
  "* Metric: " << metric->get_name()
  << std::endl << 
  "* Position1: " << new_coords[0][0] << " " << new_coords[0][1] << " " << new_coords[0][2]
  << std::endl << 
  "* Position2: " << new_coords[1][0] << " " << new_coords[1][1] << " " << new_coords[1][2]
  << std::endl <<  
  "**************************************************************************" 
  << std::endl;
  
  CPPUNIT_ASSERT_VECTORS_EQUAL( coordinates[0], new_coords[0], TOL );
  CPPUNIT_ASSERT_VECTORS_EQUAL( coordinates[1], new_coords[1], TOL );
  return false;
}


