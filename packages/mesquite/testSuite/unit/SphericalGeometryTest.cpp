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
//    AUTHOR: Michael Brewer
//       ORG: Sandia National Labs
//    E-MAIL: mbrewer@sandia.gov
//
// ORIG-DATE: Jan. 29, 2003
//  LAST-MOD: 25-Feb-04 at 10:49:04 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file SphericalGeometryTest.cpp

Regression testing using the spherical geometry capabilities in
SimplifiedGeometryEngine.
 */
// DESCRIP-END.
//

#include "meshfiles.h"

#include "PatchDataInstances.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include <math.h>

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "PatchData.hpp"
//#include "StoppingCriterion.hpp"
#include "QualityAssessor.hpp"

#include "IdealWeightInverseMeanRatio.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LaplacianSmoother.hpp"
#include "SmartLaplacianSmoother.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"
#include "AspectRatioGammaQualityMetric.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "MultiplyQualityMetric.hpp"
#include "SphericalDomain.hpp"

#include "MeshImpl.hpp"

#include <iostream>
using std::cout;
using std::endl;
using namespace Mesquite;

class SphericalGeometryTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(SphericalGeometryTest);
    //run cg on the quad sphere mesh with condition number L2
  CPPUNIT_TEST (test_cg_mesh_cond_sphere);
    //run cg on the quad sphere mesh with condition number L2
  CPPUNIT_TEST (test_smart_lapl_sphere);
    //run laplacian smoothing on the geo tri mesh
  CPPUNIT_TEST (test_lapl_geo_sphere);
  
  CPPUNIT_TEST_SUITE_END();
  
private:
  double qualTol;//double used for double comparisons
  int pF;//PRINT_FLAG
public:
  void setUp()
  {
      //pF=1;//PRINT_FLAG IS ON
      pF=0;//PRINT_FLAG IS OFF
        //tolerance double
      qualTol=MSQ_MIN;
  }

  void tearDown()
  {
  }
  
public:
  SphericalGeometryTest()
    {}
  
   void test_cg_mesh_cond_sphere()
   {
     Mesquite::MeshImpl mesh;
     Mesquite::MsqPrintError err(cout);
     mesh.read_vtk(MESH_FILES_DIR "2D/vtk/quads/untangled/quads_on_sphere_529.vtk", err);
     CPPUNIT_ASSERT(!err);
     
       //create geometry: sphere, center (2,2,0), radius 3
     Vector3D center(2,2,0);
     SphericalDomain msq_geom(center, 3.0);
     
       // creates an intruction queue
     InstructionQueue queue1;
     
       // creates a mean ratio quality metric ...
     ConditionNumberQualityMetric shape;
     UntangleBetaQualityMetric untan;
     
       // ... and builds an objective function with it
     LPtoPTemplate obj_func(&shape, 2, err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
       // creates the steepest descent optimization procedures
     ConjugateGradient pass1( &obj_func, err );
       //SteepestDescent* pass2 = new SteepestDescent( obj_func );
     pass1.use_global_patch();
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     QualityAssessor qa=QualityAssessor( &shape );
    
       //**********Set stopping criterion  5 iterates ****************
       //StoppingCriterion sc5(StoppingCriterion::NUMBER_OF_PASSES,5);
       //pass1->set_stopping_criterion(&sc5);
     TerminationCriterion sc5;
     sc5.add_iteration_limit( 5 );
     pass1.set_inner_termination_criterion(&sc5);
       //CG's debugging print, increase integer to get more print info
     pass1.set_debugging_level(0);
  
       //  queue1.add_preconditioner(pass2, err); CPPUNIT_ASSERT(!err);
     queue1.set_master_quality_improver(&pass1, err); CPPUNIT_ASSERT(!err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
       // launches optimization on mesh_set1
     MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &msq_geom);
     double orig_qa_val=qa.loop_over_mesh(&mesh_and_domain, 0, err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     queue1.run_instructions(&mesh_and_domain, err); CPPUNIT_ASSERT(!err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     double fin_qa_val=qa.loop_over_mesh(&mesh_and_domain, 0, err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
       //make sure 'quality' improved
     CPPUNIT_ASSERT( fin_qa_val <= orig_qa_val );
   }
   void test_smart_lapl_sphere()
     {
       Mesquite::MeshImpl mesh;
       Mesquite::MsqPrintError err(cout); 
       mesh.read_vtk(MESH_FILES_DIR "2D/vtk/quads/untangled/quads_on_sphere_529.vtk", err);
       
         //create geometry sphere:  ratius 1, centered at (0,0,0)
       Vector3D center(2,2,0);
       SphericalDomain msq_geom(center, 3.0);
  
         // creates an intruction queue
       InstructionQueue queue1;

         // creates an edge length metric ...
       IdealWeightInverseMeanRatio shape_metric(err);
       LInfTemplate shape_func(&shape_metric);
       
         //create the smart laplacian smoother
       SmartLaplacianSmoother s_lapl(&shape_func);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);

         //*******Set stopping criterion 5 iterates  ***********
       TerminationCriterion sc5;
       sc5.add_iteration_limit( 5 );
       s_lapl.set_outer_termination_criterion(&sc5);
         //qa, qi, qa
       queue1.set_master_quality_improver(&s_lapl, err); CPPUNIT_ASSERT(!err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         // launches optimization on mesh_set1
       QualityAssessor qa=QualityAssessor( &shape_metric );
       MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &msq_geom);
       double orig_val=qa.loop_over_mesh(&mesh_and_domain, 0, err);
       
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       queue1.run_instructions(&mesh_and_domain, err); CPPUNIT_ASSERT(!err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);

       double final_val= qa.loop_over_mesh(&mesh_and_domain, 0, err);
  
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         //make sure 'quality' improved
       CPPUNIT_ASSERT( final_val < orig_val );
     }
  
  void test_lapl_geo_sphere()
     {
       Mesquite::MeshImpl mesh;
       Mesquite::MsqPrintError err(cout);
       
       mesh.read_vtk(MESH_FILES_DIR "2D/vtk/tris/untangled/Mesquite_geo_10242.vtk", err);
       
         //create geometry sphere:  ratius 1, centered at (0,0,0)
       Vector3D center(0,0,0);
       Mesquite::SphericalDomain msq_geom(center, 1.0);
  
         // creates an intruction queue
       InstructionQueue queue1;

         // creates an edge length metric ...
       EdgeLengthQualityMetric edg_len;
      
         //create the laplacian smoother
       LaplacianSmoother lapl;
         //Make sure no errors
       CPPUNIT_ASSERT(!err);

         //create a quality assessor
       QualityAssessor qa=QualityAssessor( &edg_len );

         //*******Set stopping criterion 10 iterates  ***********
         //StoppingCriterion sc10(StoppingCriterion::NUMBER_OF_PASSES,10);
         //lapl->set_stopping_criterion(&sc10);
       TerminationCriterion sc10;
       sc10.add_iteration_limit( 10 );
       lapl.set_outer_termination_criterion(&sc10);
         //qa, qi, qa
       queue1.set_master_quality_improver(&lapl, err); 
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         // launches optimization on mesh_set1
       MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &msq_geom);
       double orig_qa_val=qa.loop_over_mesh(&mesh_and_domain, 0, err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       queue1.run_instructions(&mesh_and_domain, err); CPPUNIT_ASSERT(!err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       double fin_qa_val=qa.loop_over_mesh(&mesh_and_domain, 0, err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         //make sure 'quality' improved
       CPPUNIT_ASSERT( fin_qa_val <= orig_qa_val );
     }
  
   
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(SphericalGeometryTest, "SphericalGeometryTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(SphericalGeometryTest, "Regression");
