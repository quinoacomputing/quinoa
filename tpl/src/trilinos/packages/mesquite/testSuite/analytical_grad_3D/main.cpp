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
//  LAST-MOD: 23-Jul-03 at 18:09:13 by Thomas Leurent
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
#include "SteepestDescent.hpp"

using namespace Mesquite;


int main()
{
  MsqPrintError err(std::cout);
  Mesquite::MeshImpl mesh;
  mesh.read_vtk(MESH_FILES_DIR "3D/vtk/hexes/untangled/hexes_4by2by2.vtk", err);
  
    // creates an intruction queue
  InstructionQueue queue1;
  
    // creates a mean ratio quality metric ...
  IdealWeightInverseMeanRatio mean_ratio(err);
  if (err) return 1;
  ConditionNumberQualityMetric cond_num;
  mean_ratio.set_averaging_method(QualityMetric::LINEAR);
  
    // ... and builds an objective function with it
    //LInfTemplate* obj_func = new LInfTemplate(mean_ratio);
  LPtoPTemplate obj_func(&mean_ratio, 2, err);
  if (err) return 1;
   // creates the steepest descent optimization procedures
  SteepestDescent pass1( &obj_func );
  pass1.use_global_patch();
  //if (err) return 1;
  //pass1.set_maximum_iteration(6);
  
  QualityAssessor stop_qa=QualityAssessor(&mean_ratio);
  if (err) return 1;
  stop_qa.add_quality_assessment(&cond_num);
  if (err) return 1;
  
  
   //**************Set stopping criterion****************
// StoppingCriterion sc1(&stop_qa,1.0,1.8);
    //StoppingCriterion sc2(StoppingCriterion::NUMBER_OF_PASSES,1);
  TerminationCriterion tc2;
  tc2.add_iteration_limit( 1 );
// CompositeAndStoppingCriterion sc(&sc1,&sc2);
  pass1.set_inner_termination_criterion(&tc2);

  // adds 1 pass of pass1 to mesh_set1
//  queue1.add_preconditioner(pass1, err); 
//  if (err) return 1;
  queue1.add_quality_assessor(&stop_qa,err);
  queue1.set_master_quality_improver(&pass1, err); 
  if (err) return 1;
  queue1.add_quality_assessor(&stop_qa,err);
  if (err) return 1;
  // adds 1 passes of pass2 to mesh_set1
//  mesh_set1.add_quality_pass(pass2);

  mesh.write_vtk("original_mesh.vtk", err); 
  if (err) return 1;
  
    // launches optimization on mesh_set1
  queue1.run_instructions(&mesh, err);
  if (err) return 1;
  
  mesh.write_vtk("smoothed_mesh.vtk", err); 
  if (err) return 1;
  
  return 0;
}
