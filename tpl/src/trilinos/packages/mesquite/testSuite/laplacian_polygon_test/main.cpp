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
#include "MsqError.hpp"
#include "MeshImpl.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "PlanarDomain.hpp"
#include "MsqTimer.hpp"

// algorythms
#include "ConditionNumberQualityMetric.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"
#include "LaplacianSmoother.hpp"
#include "EdgeLengthQualityMetric.hpp"
using namespace Mesquite;

const char DEFAULT_INPUT[] = MESH_FILES_DIR "2D/vtk/N-Polygonal/poly1.vtk";

void help(const char* argv0)
{
  std::cerr << "Usage: " << argv0 << " [<input_file>] [<output_file>]" << std::endl
            << "  default input file is: " << DEFAULT_INPUT << std::endl
            << "  defualt is no output file" << std::endl;
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
//  queue1.add_quality_assessor(&stop_qa,err); 
  if (err) return 1;
  queue1.set_master_quality_improver(&lapl1, err); 
  if (err) return 1;
//  queue1.add_quality_assessor(&stop_qa,err); 
  if (err) return 1;
  
  PlanarDomain plane(Vector3D(0,0,1), Vector3D(0,0,0));
  
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
  
  return 0;
}
