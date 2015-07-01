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

/** \file main.cpp 
 * \brief test NonGradient Solver (barrier and non-barrier).
 * \author Boyd Tidwell
 */
#include "meshfiles.h"
#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "PlanarDomain.hpp"
#include "XYPlanarDomain.hpp"
#include "MeshInterface.hpp"

#include "TShapeNB1.hpp"   
#include "TShapeB1.hpp"
#include "TQualityMetric.hpp"    
#include "IdealShapeTarget.hpp"
#include "MaxTemplate.hpp"
#include "PMeanPTemplate.hpp"
#include "ElementMaxQM.hpp"
#include "ElementAvgQM.hpp"
#include "ElementPMeanP.hpp"
#include "ElemSampleQM.hpp"
#include "TargetCalculator.hpp"   
#include "PMeanPTemplate.hpp"
#include "LPtoPTemplate.hpp"
#include "ElementPMeanP.hpp"
#include "NonGradient.hpp"   

using namespace Mesquite;

int main()
{
  MsqPrintError err(std::cout);
  PlanarDomain xyPlane(PlanarDomain::XY, -5);

  #define FILE_NAME1 "bad_circle_tri.vtk"
  #define FILE_NAME2 "tangled_tri.vtk"
  const char file_name1[] = MESH_FILES_DIR "2D/vtk/tris/untangled/" FILE_NAME1;
  const char file_name2[] = MESH_FILES_DIR "2D/vtk/tris/tangled/" FILE_NAME2;

    // Barrier / Max Objective Function Test

  Mesquite::MeshImpl mesh_max;
  mesh_max.read_vtk(file_name1, err);
  if (err)
  {
    std::cerr << "NonGradient Barrier test: failed to read file." << std::endl;
    return 1;
  }

  IdealShapeTarget target_max;

  TShapeB1 mu;
  TQualityMetric tqMetric_max( &target_max, &mu );
  ElemSampleQM* sampleMetric(&tqMetric_max);
  ElementMaxQM maxMetric( &tqMetric_max );
  ElementPMeanP meanpMetric( 1.0, sampleMetric);

  MaxTemplate maxObjFunction(&maxMetric);  // max(max)
  NonGradient max_opt( &maxObjFunction ); // optimization procedure 

  PMeanPTemplate pmeanpObjFunction(1.0, sampleMetric);
  NonGradient pmeanp_opt (&pmeanpObjFunction);

  LPtoPTemplate PtoPObjMaxfunction(&maxMetric, (short int)1.0, err);  // max(max)

  // Processing for Max Objective Function

  max_opt.setSimplexDiameterScale(0); 
  max_opt.use_element_on_vertex_patch(); // local patch
  max_opt.set_debugging_level(0);

  // Construct and register the Quality Assessor instances
  QualityAssessor max_initial_qa=QualityAssessor(&maxMetric, 10);
  QualityAssessor maxObj_max_optimal_qa=QualityAssessor(&maxMetric, 10);
   
  //**************Set stopping criterion****************
  TerminationCriterion innerTC, outerTC;

  outerTC.add_iteration_limit(40);
  innerTC.add_iteration_limit(20);
  max_opt.set_outer_termination_criterion(&outerTC);
  max_opt.set_inner_termination_criterion(&innerTC);

  // test for barrier violation
  PlanarDomain xyPlane2 (PlanarDomain::XY, 5);
  Mesquite::MeshImpl mesh_bv;
  mesh_bv.read_vtk(file_name2, err);
  if (err)
  {
    std::cerr << "NonGradient Barrier Violation test: failed to read file." << std::endl;
    return 1;
  }

  InstructionQueue queue1; 
  queue1.add_quality_assessor(&max_initial_qa,err);
  if (err) return 1;

  queue1.set_master_quality_improver(&max_opt, err);   // max
  if (err) return 1;

  queue1.add_quality_assessor(&maxObj_max_optimal_qa,err);
  if (err) return 1;

  Mesquite::MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh_max, &xyPlane);
  queue1.run_instructions(&mesh_and_domain, err);
  if (err) return 1;

    // Non-Barrier / Ave Objective Function Test

  Mesquite::MeshImpl mesh_mean;
  mesh_mean.read_vtk(file_name1, err);
  if (err)
  {
    std::cerr << "NonGradient Non-barrier test: failed to read file." << std::endl;
    return 1;
  }

  TShapeNB1 nonBarrier;
  TargetCalculator* target;
  IdealShapeTarget ident_target; 

  target = &ident_target;

  TQualityMetric tqMetric( target, &nonBarrier );
 
  ElementPMeanP mean_metric( 1.0, &tqMetric );
  PMeanPTemplate meanObjfunction( 1.0, &mean_metric );  // ave(ave)
 
  // Processing for Mean Objective Function
  NonGradient mean_opt( &meanObjfunction ); // optimization procedure 
  mean_opt.setExactPenaltyFunction(false); // allow infeasible
  mean_opt.use_element_on_vertex_patch(); // local patch

    // Set Termination Criteria
  TerminationCriterion innerTC_mean, outerTC_mean;
  outerTC_mean.add_iteration_limit(10);
  innerTC_mean.add_iteration_limit(30);
  mean_opt.set_outer_termination_criterion(&outerTC_mean);
  mean_opt.set_inner_termination_criterion(&innerTC_mean);
  mean_opt.set_debugging_level(0);

  // Construct and register the meanObj Quality Assessor instance
  QualityAssessor mean_qa=QualityAssessor(&mean_metric, 10);

  InstructionQueue queue2; 
  queue2.add_quality_assessor(&mean_qa, err);
  if (err) return 1;
  queue2.set_master_quality_improver(&mean_opt, err); 
  if (err) return 1;
  queue2.add_quality_assessor(&mean_qa, err);
  if (err) return 1;

  queue2.run_instructions(&mesh_and_domain, err);
  if (err) return 1;



    // Test barrier target metric using objective function MaxTemplate with inverted mesh
  InstructionQueue queue3; 
  queue3.add_quality_assessor(&max_initial_qa,err);
  if (err) return 1;
  queue3.set_master_quality_improver(&max_opt, err); 
  if (err) return 1;
  Mesquite::MeshDomainAssoc mesh_and_domain2 = MeshDomainAssoc(&mesh_bv, &xyPlane2);
  queue3.run_instructions(&mesh_and_domain2, err);
  if (err.error_code() == err.BARRIER_VIOLATED)
  {
    std::cerr << std::endl << "MaxTemplate OF with inverted mesh test passed" << std::endl;
    err.clear();
  }
  else
    return 1;


    // Test barrier target metric using objective function PMeanPTemplate with inverted mesh
  InstructionQueue queue4; 
  queue4.set_master_quality_improver(&pmeanp_opt, err); 
  if (err) return 1;
  queue4.run_instructions(&mesh_and_domain2, err);
  if (err.error_code() == err.BARRIER_VIOLATED)
  {
    std::cerr << std::endl << "PMeanPTemplate OF with inverted mesh test passed" << std::endl << std::endl;
    err.clear();
  }
  else
    return 1;

  return 0;
}
