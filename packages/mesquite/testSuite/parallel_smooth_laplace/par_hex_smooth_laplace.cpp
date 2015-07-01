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
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 24-Jan-12
//  LAST-MOD: 26-Jan-1 by Stephen Kennon
//
//
// DESCRIPTION:
// ============
/*! \file par_hex.cpp

A test of Mesquite's parallel capabilities.  Reads a split vtk file, smooths in parallel using Laplace
smoothing, writes out the result (which can be compared with the "gold" copy of the same name in the
meshFiles VTK directory).

See the Mesquite User's Guide, section "Using Mesquite in Parallel" - this code is very similar to the
example code shown therein.

 */
// DESCRIP-END.
//

#include "meshfiles.h"

#include "MeshImpl.hpp"
#include "MeshUtil.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"
#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "LaplaceWrapper.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

/* Mesquite includes */
#include "ParallelMeshImpl.hpp"
#include "ParallelHelper.hpp"


// algorithms
#include "Randomize.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"
#include "PlanarDomain.hpp"

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

#include <mpi.h>
#include <sstream>

using namespace Mesquite;

#define VTK_3D_DIR MESH_FILES_DIR "3D/vtk/hexes/tangled/"

using namespace std;
  
int main( int argc, char* argv[] )
{
  /* init MPI */
  int rank, nprocs;
  if (MPI_SUCCESS != MPI_Init(&argc, &argv)) {
    cerr << "MPI_Init failed." << endl;
    return 2;
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (nprocs > 2) { cerr << "parallel_laplace_smooth test can only be run with 1 or 2 processors" << std::endl; return 0; }

  const int debug=0;  // 1,2,3 for more debug info
  if (debug)
    {
      MsqDebug::enable(1);
      if (debug > 1) MsqDebug::enable(2);
      if (debug > 2) MsqDebug::enable(3);
    }

  /* create processor-specific file names */
  ostringstream in_name, out_name, gold_name;
  in_name << VTK_3D_DIR << "par_original_hex_mesh." << nprocs << "." << rank << ".vtk";
  gold_name << VTK_3D_DIR << "par_smoothed_hex_mesh." << nprocs << "." << rank << ".vtk";
  out_name << "par_smoothed_hex_mesh." << nprocs << "." << rank << ".vtk";

  /* load different mesh files on each processor */
  Mesquite::MsqError err;
  Mesquite::MeshImpl mesh;
  mesh.read_vtk(in_name.str().c_str(), err);
  if (err) {cerr << err << endl; return 1;}

  /* create parallel mesh instance, specifying tags 
   * containing parallel data */
  Mesquite::ParallelMeshImpl parallel_mesh(&mesh, "GLOBAL_ID", "PROCESSOR_ID");
  Mesquite::ParallelHelperImpl helper;
  helper.set_communicator(MPI_COMM_WORLD);
  helper.set_parallel_mesh(&parallel_mesh);
  parallel_mesh.set_parallel_helper(&helper);

  /* do Laplacian smooth */
  LaplaceWrapper optimizer;
  optimizer.set_vertex_movement_limit_factor(1.e-10);
  optimizer.set_iteration_limit(2000);
  optimizer.enable_culling(false);
  optimizer.run_instructions(&parallel_mesh, err);
  if (err) {cerr << err << endl; return 1; }

  /* write mesh */
  mesh.write_vtk(out_name.str().c_str(),err);
  if (err) {cerr << err << endl; return 1;}

  /* compare mesh with gold copy */
  MeshImpl gold;
  gold.read_vtk(gold_name.str().c_str(),err);
  if (err) {cerr << err << endl; return 1;}

  bool do_print=true;
  double tol = 1.e-4;
  bool diff = MeshUtil::meshes_are_different(mesh, gold, err, tol, do_print);
  if (err) {cerr << err << endl; return 1;}

  if (diff) {cerr << "Error, computed mesh different from gold copy" << std::endl; return 1;}
  
  print_timing_diagnostics(cout);

  MPI_Finalize();
  return 0;
}
