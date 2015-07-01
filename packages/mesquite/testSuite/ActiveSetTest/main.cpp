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
//  LAST-MOD: 23-Jul-03 at 18:10:35 by Thomas Leurent
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
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorithms
#include "ConditionNumberQualityMetric.hpp"
#include "NonSmoothDescent.hpp"

#include "MeshImpl.hpp"
using namespace Mesquite;


int main()
{     
    /* Reads a Mesh file */
  const char *file_name = 
//      MESH_FILES_DIR "2D/vtk/tris/untangled/equil_tri2.vtk";
//      MESH_FILES_DIR "2D/vtk/tris/untangled/tri_20258.vtk";
//      MESH_FILES_DIR "3D/vtk/tets/untangled/tet_1.vtk";
//      MESH_FILES_DIR "3D/vtk/hexes/untangled/cube_tet_2.vtk";
     MESH_FILES_DIR "3D/vtk/tets/untangled//tire.vtk";
  printf("Loading mesh set 1\n");
  MsqPrintError err( cout );
  Mesquite::MeshImpl mesh;
  mesh.read_vtk(file_name, err);
  if (err) return 1;
  
    // Creates an intruction queue
    //  printf("Creating instruction queue\n");
  InstructionQueue queue1;

  // Creates a condition number quality metric 
  //  printf("Creating quality metric\n");
  ConditionNumberQualityMetric cond_no;
  
  // Create the NonSmooth Steepest Descent procedures
  //  printf("creating optimizer\n");
  NonSmoothDescent minmax_method( &cond_no );

  // Set a termination criterion
  TerminationCriterion tc2;
  tc2.add_iteration_limit( 1 );
  minmax_method.set_outer_termination_criterion(&tc2);
  // Set up the quality assessor
  //  printf("Setting up the quality assessor\n");
  QualityAssessor quality_assessor=QualityAssessor(&cond_no);

  // assess the quality of the initial mesh
  queue1.add_quality_assessor(&quality_assessor, err); 
  if (err) return 1;

  // Set the max min method to be the master quality improver
  queue1.set_master_quality_improver(&minmax_method, err); 
  if (err) return 1;

  // assess the quality of the final mesh
  queue1.add_quality_assessor(&quality_assessor, err); 
  if (err) return 1;

  // write out the original mesh
  //  printf("Writing out the original mesh\n");
  mesh.write_vtk("original_mesh.vtk", err); 
  if (err) return 1;

  // launches optimization on mesh_set1
  //  printf("Running the instruction queue\n");
  queue1.run_instructions(&mesh, err); 
  if (err) return 1;

  // write out the smoothed mesh
  //  printf("Writing out the final mesh\n");
  mesh.write_vtk("smoothed_mesh.vtk", err); 
  if (err) return 1;
  
  return 0;
}
