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
//  LAST-MOD: 23-Jul-03 at 18:11:05 by Thomas Leurent
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
#include "MeshImpl.hpp"
#include "PlanarDomain.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "MsqError.hpp"
#include "ShapeImprovementWrapper.hpp"
// algorythms
#include "IdealWeightInverseMeanRatio.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

using namespace Mesquite;

int main()
{
  Mesquite::MeshImpl mesh;
  MsqPrintError err(cout);
    //create geometry: plane z=0, normal (0,0,1)
  Vector3D pnt(0,0,0);
  Vector3D s_norm(0,0,1);
  Mesquite::PlanarDomain msq_geom(s_norm, pnt);
     
  mesh.read_vtk(MESH_FILES_DIR "2D/vtk/mixed/untangled/hybrid_3quad_1tri.vtk", err);
  if (err) return 1;
  
    // creates an intruction queue
  InstructionQueue queue1;
  
    // creates a mean ratio quality metric ...
  IdealWeightInverseMeanRatio mean_ratio(err);
  if (err) return 1;
    //   mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
    //   mean_ratio->set_hessian_type(QualityMetric::NUMERICAL_HESSIAN);
    //mean_ratio->set_averaging_method(QualityMetric::SUM, err);
    //MSQ_CHKERR(err);
  
    // ... and builds an objective function with it
  LPtoPTemplate obj_func(&mean_ratio, 2, err);
  if (err) return 1;;
  
    // creates the steepest descent, feas newt optimization procedures
    //ConjugateGradient* pass1 = new ConjugateGradient( &obj_func, err );
  SteepestDescent pass1( &obj_func );
  pass1.use_global_patch();
  if (err) return 1;;
  
  QualityAssessor qa=QualityAssessor(&mean_ratio);
  if (err) return 1;;
  
    // **************Set termination criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_iteration_limit( 1 );
    //_inner.add_absolute_quality_improvement( OF_value );
    //tc_inner.add_absolute_gradient_L2_norm( OF_value );
  TerminationCriterion tc_outer;
    //tc_outer.add_iteration_limit( 1 );
  tc_outer.add_iteration_limit( 1 );
  
  pass1.set_inner_termination_criterion(&tc_inner);
  pass1.set_outer_termination_criterion(&tc_outer);

  queue1.add_quality_assessor(&qa,err); 
  if (err) return 1;
    // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(&pass1, err);
  if (err) return 1;
  queue1.add_quality_assessor(&qa,err); 
  if (err) return 1;
  mesh.write_vtk("original_mesh.vtk",err); 
  if (err) return 1;
  
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &msq_geom);
  queue1.run_instructions(&mesh_and_domain, err); 
  if (err) return 1;
  mesh.write_vtk("smoothed_mesh.vtk",err); 
  if (err) return 1;
    //std::cout<<"\n\nNow running the shape wrapper.\n=n";
    //ShapeImprovementWrapper wrap(100);
    //wrap.run_instructions(mesh_set1, err); MSQ_CHKERR(err);
  print_timing_diagnostics(cout);
  return 0;
}
 
