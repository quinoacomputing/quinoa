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
//  LAST-MOD: 23-Jul-03 at 18:08:13 by Thomas Leurent
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

#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_MsqTimer.hpp"
#include "Mesquite.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_Vector3D.hpp"
#include "Mesquite_InstructionQueue.hpp"
#include "Mesquite_PatchData.hpp"
#include "Mesquite_TerminationCriterion.hpp"
#include "Mesquite_QualityAssessor.hpp"

// algorithms
#include "Mesquite_Randomize.hpp"
#include "Mesquite_ConditionNumberQualityMetric.hpp"
#include "Mesquite_UntangleBetaQualityMetric.hpp"
#include "Mesquite_LPtoPTemplate.hpp"
#include "Mesquite_LInfTemplate.hpp"
#include "Mesquite_SteepestDescent.hpp"
#include "Mesquite_ConjugateGradient.hpp"
#include "Mesquite_PlanarDomain.hpp"

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

using namespace Mesquite;

#define VTK_2D_DIR MESH_FILES_DIR "2D/vtk/"

const bool brief_output = true;
const bool write_output = false;

int main( )
{
  Mesquite::MeshImpl mesh;
  MsqPrintError err(cout);
  mesh.read_vtk(VTK_2D_DIR "quads/tangled/tangled_quad.vtk", err);
  if (err) return 1;
  
  // Set Domain Constraint
  Vector3D pnt(0,0,5);
  Vector3D s_norm(0,0,1);
  PlanarDomain msq_geom(s_norm, pnt);
                                                                              
    // creates an intruction queue
  InstructionQueue queue1;
  
    // creates a mean ratio quality metric ...
  ConditionNumberQualityMetric shape_metric;
  UntangleBetaQualityMetric untangle(2);
  Randomize pass0(.05);
    // ... and builds an objective function with it
    //LInfTemplate* obj_func = new LInfTemplate(shape_metric);
  LInfTemplate obj_func(&untangle);
  LPtoPTemplate obj_func2(&shape_metric, 2, err);
  if (err) return 1;
    // creates the steepest descent optimization procedures
  ConjugateGradient pass1( &obj_func, err );
  if (err) return 1;
  
    //SteepestDescent* pass2 = new SteepestDescent( obj_func2 );
  ConjugateGradient pass2( &obj_func2, err );
  if (err) return 1;
  pass2.use_element_on_vertex_patch();
  if (err) return 1;
  pass2.use_global_patch();
  if (err) return 1;
  QualityAssessor stop_qa=QualityAssessor(&shape_metric);
  QualityAssessor stop_qa2=QualityAssessor(&shape_metric);
  if (brief_output) {
    stop_qa.disable_printing_results();
    stop_qa2.disable_printing_results();
  }
  
  stop_qa.add_quality_assessment(&untangle);
    // **************Set stopping criterion**************
    //untangle beta should be 0 when untangled
  TerminationCriterion sc1;
  sc1.add_relative_quality_improvement( 0.000001 );
  TerminationCriterion sc3;
  sc3.add_iteration_limit( 10 );
  TerminationCriterion sc_rand;
  sc_rand.add_iteration_limit( 1 );
  
    //StoppingCriterion sc1(&stop_qa,-1.0,.0000001);
    //StoppingCriterion sc3(&stop_qa2,.9,1.00000001);
    //StoppingCriterion sc2(StoppingCriterion::NUMBER_OF_PASSES,10);
    //StoppingCriterion sc_rand(StoppingCriterion::NUMBER_OF_PASSES,1);
    //either until untangled or 10 iterations
  pass0.set_outer_termination_criterion(&sc_rand);
  pass1.set_outer_termination_criterion(&sc1);
  pass2.set_inner_termination_criterion(&sc3);
  
    // adds 1 pass of pass1 to mesh_set1
  queue1.add_quality_assessor(&stop_qa,err); 
  if (err) return 1;
    //queue1.add_preconditioner(pass0,err);MSQ_CHKERR(err);
    //queue1.add_preconditioner(pass1,err);MSQ_CHKERR(err);
    //queue1.set_master_quality_improver(pass2, err); MSQ_CHKERR(err);
  queue1.set_master_quality_improver(&pass1, err);
  if (err) return 1;
  queue1.add_quality_assessor(&stop_qa2,err);
  if (err) return 1;
  if (write_output)
    mesh.write_vtk("original_mesh.vtk", err);
  if (err) return 1;
  
    // launches optimization on mesh_set1
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &msq_geom);
  queue1.run_instructions(&mesh_and_domain, err);
  if (err) return 1;
  
  if (write_output)
    mesh.write_vtk("smoothed_mesh.vtk", err); 
  if (err) return 1;
  
  if (!brief_output)
    print_timing_diagnostics(cout);
  return 0;
}
