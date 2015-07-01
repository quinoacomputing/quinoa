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

// algorythms
#include "IdealWeightInverseMeanRatio.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
using namespace Mesquite;

int main()
{     
  MsqPrintError err(cout);
  Mesquite::MeshImpl mesh;
  mesh.read_vtk(MESH_FILES_DIR "3D/vtk/tets/untangled/tire.vtk", err);
  if (err) return 1;
  
    // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  IdealWeightInverseMeanRatio mean(err);
  if (err) return 1;
  
  LPtoPTemplate obj_func(&mean, 1, err);
  if (err) return 1;
  
  // creates the optimization procedures
//   ConjugateGradient* pass1 = new ConjugateGradient( obj_func, err );
  FeasibleNewton pass1( &obj_func );

  //perform optimization globally
  pass1.use_global_patch();
  if (err) return 1;
  
  QualityAssessor mean_qa=QualityAssessor(&mean);

    //**************Set termination criterion****************

  //perform 1 pass of the outer loop (this line isn't essential as it is
  //the default behavior).
  TerminationCriterion tc_outer;
  tc_outer.add_iteration_limit( 1 ); 
  pass1.set_outer_termination_criterion(&tc_outer);
  
  //perform the inner loop until a certain objective function value is
  //reached.  The exact value needs to be determined (about 18095).
  //As a safety, also stop if the time exceeds 10 minutes (600 seconds).
  TerminationCriterion tc_inner;
  tc_inner.add_absolute_quality_improvement( 13975 ); 
//  tc_inner.add_absolute_quality_improvement( 13964.93818 );
  tc_inner.add_cpu_time( 1800 ); 
  
  pass1.set_inner_termination_criterion(&tc_inner);
  
  //used for cg to get some info
//  pass1->set_debugging_level(2);
  
  // adds 1 pass of pass1 to mesh_set1
  queue1.add_quality_assessor(&mean_qa,err);
  if (err) return 1;
  queue1.set_master_quality_improver(&pass1, err); 
  if (err) return 1;
  queue1.add_quality_assessor(&mean_qa,err);
  if (err) return 1;
  mesh.write_vtk("original_mesh.vtk", err); 
  if (err) return 1;
  
    // launches optimization on mesh_set1
  queue1.run_instructions(&mesh, err); 
  if (err) return 1;
  
  mesh.write_vtk("smoothed_mesh.vtk", err); 
  if (err) return 1;
  print_timing_diagnostics( cout );
  return 0;
}
