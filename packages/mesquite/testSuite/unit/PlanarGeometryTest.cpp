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
//    AUTHOR: Michael Brewer
//       ORG: Sandia National Labs
//    E-MAIL: mbrewer@sandia.gov
//
// ORIG-DATE: Jan. 29, 2003
//  LAST-MOD: 25-Feb-04 at 10:49:32 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file PlanarGeometryTest.cpp

Regression testing using the planar geometry capabilities in
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
#include "InstructionQueue.hpp"
#include "PatchData.hpp"
//#include "StoppingCriterion.hpp"
#include "QualityAssessor.hpp"

#include "ConditionNumberQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LaplacianSmoother.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"
#include "AspectRatioGammaQualityMetric.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "MultiplyQualityMetric.hpp"
#include "PlanarDomain.hpp"

#include "MsqFPE.hpp"

#include "UnitUtil.hpp"

#include "MeshImpl.hpp"
#include <iostream>
using std::cout;
using std::endl;

using namespace Mesquite;

class PlanarGeometryTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(PlanarGeometryTest);
    //run steepest descent on the tangled_tri.vtk mesh
  CPPUNIT_TEST (test_plane_tri_tangled);
    //run cg on tangled_quad.vtk mesh
  CPPUNIT_TEST (test_plane_quad_tangled);
    //run cg with asm metric on tri mesh in y=-5 plane
  CPPUNIT_TEST (test_plane_tri_xz);
    // test fit plane
  CPPUNIT_TEST (test_fit_plane);
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
  PlanarGeometryTest()
    {}
  
   void test_plane_tri_tangled()
   {
     Mesquite::MsqPrintError err(cout); 
     Mesquite::MeshImpl mesh;
     
      // This test doesn't use InstructionQueue, so 
      // we need to set up trapping of floating-point
      // exceptions ourself.
     MsqFPE fpe_trap( true );
     
     mesh.read_vtk(MESH_FILES_DIR "2D/vtk/tris/tangled/tangled_tri.vtk", err);
     CPPUNIT_ASSERT(!err);
     
       //create geometry: plane z=5, normal (0,0,1)
     Vector3D pnt(0,0,5);
     Vector3D s_norm(0,0,1);
     Mesquite::PlanarDomain msq_geom(s_norm, pnt);
     
       // creates an intruction queue
     InstructionQueue queue1, queue2;
     
       // creates a mean ratio quality metric ...
     ConditionNumberQualityMetric shape;
     UntangleBetaQualityMetric untan(.1);
     
       // ... and builds an objective function with it (untangle)
     LInfTemplate untan_func(&untan);
     LPtoPTemplate shape_func(&shape,2,err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
       // creates the steepest descent optimization procedures
     SteepestDescent pass1( &untan_func );
     SteepestDescent pass2( &shape_func );
     pass1.use_element_on_vertex_patch();
     pass2.use_global_patch();
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     QualityAssessor stop_qa=QualityAssessor( &untan );
     QualityAssessor qa=QualityAssessor( &shape );
     if(pF==0){
       stop_qa.disable_printing_results();
       qa.disable_printing_results();
     }  
       //**********Set stopping criterion  untangle ver small ********
       //StoppingCriterion sc_qa(&stop_qa,-100,MSQ_MIN);
       //pass1->set_stopping_criterion(&sc_qa);
     TerminationCriterion sc_of, sc_inner;
     sc_of.add_iteration_limit( 10 );
     pass1.set_outer_termination_criterion(&sc_of);
     sc_inner.add_iteration_limit( 6 );
     //sc_inner.add_absolute_gradient_L2_norm( 0.01 );
     pass1.set_inner_termination_criterion(&sc_inner);
    
       //**********Set stopping criterion  5 iterates ****************
       //StoppingCriterion sc5(StoppingCriterion::NUMBER_OF_PASSES,5);
       //pass2->set_stopping_criterion(&sc5);
     TerminationCriterion sc5;
     sc5.add_iteration_limit( 5 );
     pass2.set_inner_termination_criterion(&sc5);
       //TerminationCriterion sc_inner;
       //sc_inner.add_iteration_limit( 5 );
       //pass2->set_inner_termination_criterion(&sc_inner);
       //pass2->set_maximum_iteration(5);
  
     queue1.set_master_quality_improver(&pass1, err); CPPUNIT_ASSERT(!err);
     queue2.set_master_quality_improver(&pass2, err); CPPUNIT_ASSERT(!err);
       //********************UNTANGLE*******************************
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
       // launches optimization on mesh_set1
     MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &msq_geom);
     double orig_qa_val=stop_qa.loop_over_mesh(&mesh_and_domain, 0, err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     queue1.run_instructions(&mesh_and_domain, err); CPPUNIT_ASSERT(!err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     double fin_qa_val=stop_qa.loop_over_mesh(&mesh_and_domain, 0, err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
       //make sure 'quality' improved
     CPPUNIT_ASSERT( (fin_qa_val-orig_qa_val) <= 0.0 );
       //make sure sc_qa really was satisfied
     CPPUNIT_ASSERT( fin_qa_val <= MSQ_MIN );

      //********************SMOOTH*******************************
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
       // launches optimization on mesh_set1
     orig_qa_val=qa.loop_over_mesh(&mesh_and_domain, 0, err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     queue2.run_instructions(&mesh_and_domain, err); CPPUNIT_ASSERT(!err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
     fin_qa_val=qa.loop_over_mesh(&mesh_and_domain, 0, err);
       //Make sure no errors
     CPPUNIT_ASSERT(!err);
       //make sure 'quality' improved
     CPPUNIT_ASSERT( (fin_qa_val-orig_qa_val) <= 0.0 );
     print_timing_diagnostics(cout);
   }
  
  void test_plane_quad_tangled()
     {
       Mesquite::MeshImpl mesh;
       MsqPrintError err(cout); 
       mesh.read_vtk(MESH_FILES_DIR "2D/vtk/quads/tangled/tangled_quad.vtk", err);
       CPPUNIT_ASSERT(!err);

         //create geometry: plane z=5, normal (0,0,1)
       Vector3D pnt(0,0,5);
       Vector3D s_norm(0,0,1);
       Mesquite::PlanarDomain msq_geom(s_norm, pnt);
       
         // creates an intruction queue
       InstructionQueue queue1, queue2;

         //creates a mean ratio quality metric ...
       ConditionNumberQualityMetric shape;
       UntangleBetaQualityMetric untan(.1);
  
         // ... and builds an objective function with it (untangle)
       LInfTemplate untan_func(&untan);
       LPtoPTemplate shape_func(&shape,2,err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         // creates the cg optimization procedures
       ConjugateGradient pass1( &untan_func, err );
       ConjugateGradient pass2( &shape_func, err );
       pass1.use_element_on_vertex_patch();
       pass2.use_global_patch();
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       QualityAssessor stop_qa=QualityAssessor( &untan );
       QualityAssessor qa=QualityAssessor( &shape );
         //turn off printing if print flag not set.
       if(pF==0){
         stop_qa.disable_printing_results();
         qa.disable_printing_results();
       }
       //**********Set stopping criterion  untangle ver small ********
       //StoppingCriterion sc_qa(&stop_qa,-100,MSQ_MIN);
       //pass1->set_stopping_criterion(&sc_qa);
       TerminationCriterion sc_of;
       sc_of.add_iteration_limit( 10 );
       pass1.set_outer_termination_criterion(&sc_of);
       
         //**********Set stopping criterion  5 iterates ****************
         //StoppingCriterion sc5(StoppingCriterion::NUMBER_OF_PASSES,5);
         //pass2->set_stopping_criterion(&sc5);
       TerminationCriterion sc5;
       sc5.add_iteration_limit( 5 );
       pass2.set_inner_termination_criterion(&sc5);
         //pass2->set_maximum_iteration(5);
         //TerminationCriterion sc_inner;
         //sc_inner.add_iteration_limit( 5 );
         //pass2->set_inner_termination_criterion(&sc_inner);
       queue1.set_master_quality_improver(&pass1, err); CPPUNIT_ASSERT(!err);
       queue2.set_master_quality_improver(&pass2, err); CPPUNIT_ASSERT(!err);
         //********************UNTANGLE*******************************
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         // launches optimization on mesh_set1
       MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &msq_geom);
       double orig_qa_val=stop_qa.loop_over_mesh(&mesh_and_domain, 0, err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       queue1.run_instructions(&mesh_and_domain, err); CPPUNIT_ASSERT(!err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       double fin_qa_val=stop_qa.loop_over_mesh(&mesh_and_domain, 0, err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         //make sure 'quality' improved
       CPPUNIT_ASSERT( (fin_qa_val-orig_qa_val) <= 0.0 );
         //make sure sc_qa really was satisfied
       CPPUNIT_ASSERT( fin_qa_val <= MSQ_MIN );
       
         //********************SMOOTH*******************************
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         // launches optimization on mesh_set1
       orig_qa_val=qa.loop_over_mesh(&mesh_and_domain, 0, err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       queue2.run_instructions(&mesh_and_domain, err); CPPUNIT_ASSERT(!err);
       //Make sure no errors
       CPPUNIT_ASSERT(!err);
       fin_qa_val=qa.loop_over_mesh(&mesh_and_domain, 0, err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         //make sure 'quality' improved
       CPPUNIT_ASSERT( (fin_qa_val-orig_qa_val) <= 0.0 );
       print_timing_diagnostics(cout);
     }
  
  void test_plane_tri_xz()
     {
       MsqPrintError err(cout); 
       Mesquite::MeshImpl mesh;
       mesh.read_vtk(MESH_FILES_DIR "2D/vtk/tris/untangled/tri_5_xz.vtk", err);
       CPPUNIT_ASSERT(!err);

         //create geometry: plane y=5, normal (0,1,0)
       Vector3D pnt(0,-5,0);
       Vector3D s_norm(0,-1,0);
       Mesquite::PlanarDomain msq_geom(s_norm, pnt);
       
         // creates an intruction queue
       InstructionQueue queue1;
       
         //creates a asm quality metric ...
       ConditionNumberQualityMetric smooth;
       
         // ... and builds an objective function with it (untangle)
       LPtoPTemplate smooth_func(&smooth,1,err);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
         // creates the cg optimization procedures
       ConjugateGradient pass1( &smooth_func, err );
         //pass1->set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH, err,1 ,1);
       pass1.use_global_patch();
       pass1.set_debugging_level(1);
         //Make sure no errors
       CPPUNIT_ASSERT(!err);
       QualityAssessor qa=QualityAssessor( &smooth );
       
         //**********Set stopping criterion  5 iterates ****************
       TerminationCriterion sc5;
       sc5.add_iteration_limit( 5 );
       pass1.set_inner_termination_criterion(&sc5);
         //StoppingCriterion sc5(StoppingCriterion::NUMBER_OF_PASSES,5);
         //pass1->set_stopping_criterion(&sc5);
       TerminationCriterion sc_inner;
       sc_inner.add_iteration_limit( 5 );
       pass1.set_inner_termination_criterion(&sc_inner);
         //pass1->set_maximum_iteration(5);
       
       queue1.set_master_quality_improver(&pass1, err); CPPUNIT_ASSERT(!err);
         //********************UNTANGLE*******************************
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
       CPPUNIT_ASSERT( (fin_qa_val-orig_qa_val) <= 0.0 );
       print_timing_diagnostics(cout);
     }
  
   void test_fit_plane();
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PlanarGeometryTest, "PlanarGeometryTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PlanarGeometryTest, "Regression");

void PlanarGeometryTest::test_fit_plane()
{
  MsqPrintError err(std::cerr);
  PlanarDomain plane(PlanarDomain::XY,-1);
  const double epsilon = 1e-8;
  
  MeshImpl mesh1;
  mesh1.read_vtk(MESH_FILES_DIR "2D/vtk/tris/untangled/bad_circle_tri.vtk", err);
  ASSERT_NO_ERROR(err);
  plane.fit_vertices( &mesh1, err, epsilon );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,0,1), plane.get_normal(), epsilon );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 5, plane.get_coeff(), epsilon );
  
  MeshImpl mesh2;
  mesh2.read_vtk(MESH_FILES_DIR "2D/vtk/tris/untangled/equil_tri.vtk", err);
  ASSERT_NO_ERROR(err);
  plane.fit_vertices( &mesh2, err, epsilon );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,0,1), plane.get_normal(), epsilon );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, plane.get_coeff(), epsilon );
  
  MeshImpl mesh3;
  mesh3.read_vtk(MESH_FILES_DIR "2D/vtk/quads/untangled/quads_4by2.vtk", err);
  ASSERT_NO_ERROR(err);
  plane.fit_vertices( &mesh3, err, epsilon );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,0,1), plane.get_normal(), epsilon );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( -2, plane.get_coeff(), epsilon );
  
  MeshImpl mesh4;
  mesh4.read_vtk(MESH_FILES_DIR "2D/vtk/tris/untangled/tri_5_xz.vtk", err);
  ASSERT_NO_ERROR(err);
  plane.fit_vertices( &mesh4, err, epsilon );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(0,-1,0), plane.get_normal(), epsilon );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( -5, plane.get_coeff(), epsilon );
}
