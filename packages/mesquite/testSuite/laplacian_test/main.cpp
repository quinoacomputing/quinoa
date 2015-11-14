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
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 23-Jul-03 at 18:04:37 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#include "meshfiles.h"

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

#include "Mesquite.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_Vector3D.hpp"
#include "Mesquite_InstructionQueue.hpp"
#include "Mesquite_PatchData.hpp"
#include "Mesquite_TerminationCriterion.hpp"
#include "Mesquite_QualityAssessor.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_MsqTimer.hpp"

// algorythms
#include "Mesquite_ConditionNumberQualityMetric.hpp"
#include "Mesquite_LInfTemplate.hpp"
#include "Mesquite_SteepestDescent.hpp"
#include "Mesquite_LaplacianSmoother.hpp"
#include "Mesquite_EdgeLengthQualityMetric.hpp"
using namespace Mesquite;

const char DEFAULT_INPUT[] = MESH_FILES_DIR "2D/vtk/quads/untangled/square_quad_2.vtk";

void help(const char* argv0)
{
  std::cerr << "Usage: " << argv0 << " [<input_file>] [<output_file>]" << std::endl
            << "  default input file is: " << DEFAULT_INPUT << std::endl
            << "  defualt is no output file" << std::endl
            << "  Warning: input mesh is assumed to lie in Z=5 plane" << std::endl;
  exit(1);
}

int main(int argc, char* argv[])
{
  const char* input_file = DEFAULT_INPUT;
  const char* output_file = NULL;
  switch (argc) {
    default:
      help(argv[0]);
    case 3:
      if (!strcmp(argv[2],"-h"))
        help(argv[0]);
      output_file = argv[2];
    case 2:
      if (!strcmp(argv[1],"-h"))
        help(argv[0]);
      input_file = argv[1];
    case 1:
      ;
  }

    /* Read a VTK Mesh file */
  MsqPrintError err(cout);
  Mesquite::MeshImpl mesh;
  mesh.read_vtk( input_file, err);
  if (err) return 1;
  
    // creates an intruction queue
  InstructionQueue queue1;
  
    // creates a mean ratio quality metric ...
  ConditionNumberQualityMetric shape_metric;
  EdgeLengthQualityMetric lapl_met;
  lapl_met.set_averaging_method(QualityMetric::RMS);
 
    // creates the laplacian smoother  procedures
  LaplacianSmoother lapl1;
  QualityAssessor stop_qa=QualityAssessor(&shape_metric);
  stop_qa.add_quality_assessment(&lapl_met);
  
    //**************Set stopping criterion****************
  TerminationCriterion sc2;
  sc2.add_iteration_limit( 10 );
  if (err) return 1;
  lapl1.set_outer_termination_criterion(&sc2);
  
    // adds 1 pass of pass1 to mesh_set1
  queue1.add_quality_assessor(&stop_qa,err); 
  if (err) return 1;
  queue1.set_master_quality_improver(&lapl1, err); 
  if (err) return 1;
  queue1.add_quality_assessor(&stop_qa,err); 
  if (err) return 1;
    // adds 1 passes of pass2 to mesh_set1
    //  mesh_set1.add_quality_pass(pass2);
  
    //writeVtkMesh("original_mesh", mesh, err); MSQ_CHKERR(err);
  
  PlanarDomain plane(Vector3D(0,0,1), Vector3D(0,0,5));
  
    // launches optimization on mesh_set1
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &plane);
  Timer t;
  queue1.run_instructions(&mesh_and_domain, err); 
  if (err) return 1;
  double secs = t.since_birth();
  std::cout << "Optimization completed in " << secs << " seconds" << std::endl;
  
  if (output_file) {
    mesh.write_vtk(output_file, err); 
    if (err) return 1;
    std::cout << "Wrote file: " << output_file << std::endl;
  }
  
    // check that smoother is working: 
    // the one free vertex must be at the origin
  if (input_file == DEFAULT_INPUT) {
    std::vector<Mesh::VertexHandle> vertices;
    mesh.get_all_vertices( vertices, err );
    if (err) return 1;

    std::vector<bool> fixed_flags;
    mesh.vertices_get_fixed_flag( arrptr(vertices), fixed_flags, vertices.size(), err );
    if (err) return 1;

    // find one free vertex
    int idx = -1;
    for (unsigned i = 0; i < vertices.size(); ++i) {
      if (fixed_flags[i] == true)
        continue;
      if (idx != -1) {
        std::cerr << "Multiple free vertices in mesh." << std::endl;
        return 1;
      }
      idx = i;
    }

    if (idx == -1) {
      std::cerr << "No free vertex in mesh!!!!!" << std::endl;
      return 1;
    }

    Mesh::VertexHandle vertex = vertices[idx];
    MsqVertex coords;
    mesh.vertices_get_coordinates( &vertex, &coords, 1, err );
    if (err) return 1;

      // calculate distance from origin
    double dist = sqrt( coords[0]*coords[0] + coords[1]*coords[1] );
    if  (dist > 1e-8) {
      std::cerr << "Free vertex not at origin after Laplace smooth." << std::endl
                    << "Expected location: (0,0)" << std::endl
                    << "Actual location: (" << coords[0] << "," << coords[1] << ")" << std::endl;
      return 2;
    }
  }
  
  return 0;
}
