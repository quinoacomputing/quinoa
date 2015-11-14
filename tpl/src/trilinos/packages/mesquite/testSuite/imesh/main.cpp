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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 10-Feb-04 at 22:44:58 by Thomas Leurent
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

#include "Mesquite.hpp"
#include "Mesquite_MsqIMesh.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_InstructionQueue.hpp"
#include "Mesquite_TerminationCriterion.hpp"
#include "Mesquite_QualityAssessor.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_MeshWriter.hpp"

// algorithms
#include "Mesquite_IdealWeightInverseMeanRatio.hpp"
#include "Mesquite_EdgeLengthQualityMetric.hpp"
#include "Mesquite_LPtoPTemplate.hpp"
#include "Mesquite_FeasibleNewton.hpp"
#include "Mesquite_ConjugateGradient.hpp"
#include "Mesquite_SmartLaplacianSmoother.hpp"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include "Mesquite_iBase.h"


using namespace Mesquite;

const char* const default_file_name = MESH_FILES_DIR "3D/VTK/large_box_hex_1000.vtk";

void usage()
{
  cout << "main [-N] [filename]" << endl;
  cout << "  -N : Use native representation instead of TSTT implementation\n";
  cout << "  If no file name is specified, will use \"" 
       << default_file_name << '"' << endl;
  exit (1);
}

  // Construct a MeshTSTT from the file
Mesh* get_imesh_mesh( const char* file_name );

  // Construct a MeshImpl from the file
Mesh* get_native_mesh( const char* file_name );

  // Run FeasibleNewton solver
int run_global_smoother( Mesh* mesh, MsqError& err );

  // Run SmoothLaplacian solver
int run_local_smoother( Mesh* mesh, MsqError& err );

int main(int argc, char* argv[])
{
  Mesquite::MsqPrintError err(cout);
  
  return 0;
    // command line arguments
  const char* file_name = 0;
  bool use_native = false, opts_done = false;
  for (int arg = 1; arg < argc; ++arg)
  {
    if (!opts_done && argv[arg][0] == '-')
    {
      if (!strcmp( argv[arg], "-N"))
        use_native = true;
      else if(!strcmp( argv[arg], "--"))
        opts_done = true;
      else
        usage();
    }
    else if (!file_name)
      file_name = argv[arg];
    else
      usage();
  }
  if (!file_name)
  {
    file_name = default_file_name;
    cout << "No file specified: using default: " << default_file_name << endl;
  }  
  
    // Try running a global smoother on the mesh
  Mesh* mesh = use_native ? 
               get_native_mesh(file_name) : 
               get_imesh_mesh(file_name);
  if (!mesh) {
    std::cerr << "Failed to load input file.  Aborting." << std::endl;
    return 1;
  }
  
  MeshWriter::write_vtk(mesh, "original.vtk", err); 
  if (err) return 1;
  cout << "Wrote \"original.vtk\"" << endl;
  run_global_smoother( mesh, err );
  if (err) return 1;
  
    // Try running a local smoother on the mesh
  mesh = use_native ? 
         get_native_mesh(file_name) : 
         get_imesh_mesh(file_name);
  if (!mesh) {
    std::cerr << "Failed to load input file.  Aborting." << std::endl;
    return 1;
  }
  
  run_local_smoother( mesh, err );
  if (err) return 1;
  
  return 0;
}


int run_global_smoother( Mesh* mesh, MsqError& err )
{
  double OF_value = 0.0001;

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio(err);
  if (err) return 1;
  mean_ratio->set_averaging_method(QualityMetric::SUM, err); 
  if (err) return 1;
  
  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
  if (err) return 1;
  
  // creates the feas newt optimization procedures
  FeasibleNewton* pass1 = new FeasibleNewton( obj_func );
  pass1->use_global_patch();
  if (err) return 1;
  
  QualityAssessor stop_qa( mean_ratio );
  
  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_absolute_vertex_movement( OF_value );
  if (err) return 1;
  TerminationCriterion tc_outer;
  tc_outer.add_iteration_limit( 1 );
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;
   
  // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(pass1, err);
  if (err) return 1;
  
  queue1.add_quality_assessor(&stop_qa, err); 
  if (err) return 1;

  // launches optimization on mesh_set
  queue1.run_instructions(mesh, err);
  if (err) return 1;

  MeshWriter::write_vtk(mesh, "feasible-newton-result.vtk", err); 
  if (err) return 1;
  cout << "Wrote \"feasible-newton-result.vtk\"" << endl;

  //print_timing_diagnostics( cout );
  return 0;
}

int run_local_smoother( Mesh* mesh, MsqError& err )
{
  double OF_value = 0.0001;

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  IdealWeightInverseMeanRatio* mean_ratio = new IdealWeightInverseMeanRatio(err);
  if (err) return 1;
  mean_ratio->set_averaging_method(QualityMetric::SUM, err); 
  if (err) return 1;
  
  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(mean_ratio, 1, err);
  if (err) return 1;
  
  // creates the smart laplacian optimization procedures
  SmartLaplacianSmoother* pass1 = new SmartLaplacianSmoother( obj_func );
  
  QualityAssessor stop_qa( mean_ratio );
  
  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_absolute_vertex_movement( OF_value );
  TerminationCriterion tc_outer;
  tc_outer.add_iteration_limit( 1 );
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

  queue1.add_quality_assessor(&stop_qa, err);
  if (err) return 1;
   
  // adds 1 pass of pass1 to mesh_set
  queue1.set_master_quality_improver(pass1, err);
  if (err) return 1;
  
  queue1.add_quality_assessor(&stop_qa, err); 
  if (err) return 1;

  // launches optimization on mesh_set
  queue1.run_instructions(mesh, err);
  if (err) return 1;

  MeshWriter::write_vtk(mesh, "smart-laplacian-result.vtk", err); 
  if (err) return 1;
  cout << "Wrote \"smart-laplacian-result.vtk\"" << endl;

  //print_timing_diagnostics( cout );
  return 0;
}
 


Mesh* get_imesh_mesh( const char* file_name )
{
  int ierr;
  iMesh_Instance imesh_mesh = 0;
  iMesh_newMesh( NULL, &imesh_mesh, &ierr, 0 );
  if (iBase_SUCCESS != ierr) {
    return 0;
  }
  
  iBase_EntitySetHandle root_set;
  iMesh_getRootSet( imesh_mesh, &root_set, &ierr );
  if (iBase_SUCCESS != ierr) {
    iMesh_dtor( imesh_mesh, &ierr );
    return 0;
  }
  
  iMesh_load( imesh_mesh, root_set, file_name, 0, &ierr, strlen(file_name), 0 );
  if (iBase_SUCCESS != ierr) {
    std::cerr << file_name << ": failed to load file." << std::endl;
    iMesh_dtor( imesh_mesh, &ierr );
    return 0;
  }

  iBase_TagHandle fixed_tag;
  iMesh_getTagHandle( imesh_mesh, "fixed", &fixed_tag, &ierr, strlen("fixed") );
  if (iBase_SUCCESS != ierr) {
    iMesh_dtor( imesh_mesh, &ierr );
    return 0;
  }

  MsqError err;
  Mesh* result = new Mesquite::MsqIMesh( imesh_mesh, root_set, iBase_REGION, err, &fixed_tag );
  if (MSQ_CHKERR(err)) {
    delete result;
    cerr << err << endl;
    return 0;
  }
  
  return result;
}
  


Mesh* get_native_mesh( const char* file_name )
{
  MsqError err;
  MeshImpl* mesh = new MeshImpl;
  mesh->read_vtk( file_name, err );
  if (err)
  {
    cerr << err << endl;
    exit(3);
  }
  
  return mesh;
}


