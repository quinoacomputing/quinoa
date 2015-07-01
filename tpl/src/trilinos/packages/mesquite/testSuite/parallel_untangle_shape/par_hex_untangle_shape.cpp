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
#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "LaplaceWrapper.hpp"
#include "ShapeImprovementWrapper.hpp"
#include "UntangleWrapper.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

/* Mesquite includes */
#include "ParallelMeshImpl.hpp"
#include "ParallelHelper.hpp"
#include "MsqDebug.hpp"
#include "Settings.hpp"
//#include "ShapeImprovementWrapper.hpp"
//#include "UntangleWrapper.hpp"

#include "IdealWeightInverseMeanRatio.hpp" 

#include "UntangleBetaQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "ConjugateGradient.hpp"
#include "SteepestDescent.hpp"
//#include "FeasibleNewton.hpp"

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
#define VTK_2D_DIR MESH_FILES_DIR "2D/vtk/quads/tangled/"

using namespace std;
  
const double DEF_UNT_BETA = 1e-8;
const double DEF_SUC_EPS = 1e-4;

class ParShapeImprover
{
  int innerIter;
  double gradNorm;
public:

  ParShapeImprover(int inner_iterations=100, double grad_norm=1.e-8) : innerIter(inner_iterations),gradNorm(grad_norm) {}

  class ParShapeImprovementWrapper : public Wrapper {
     
  public:  
        
    //Constructor sets the instructions in the queue.
    ParShapeImprovementWrapper(int inner_iterations = 100,
                               double cpu_time = 0.0, 
                               double grad_norm =1.e-8,
                               int parallel_iterations = 10)
      : innerIter(inner_iterations),
        maxTime(cpu_time), 
        gradNorm(grad_norm),
        untBeta(DEF_UNT_BETA),
        successiveEps(DEF_SUC_EPS),
        parallelIterations(parallel_iterations),
        m_do_untangle_only(false)
    {}


  protected:

	    void run_wrapper( MeshDomainAssoc* mesh_and_domain,
			      ParallelMesh* pmesh,
			      Settings* settings,
			      QualityAssessor* qa,
			      MsqError& err );
	      
	  private:

	    int innerIter;
	    double maxTime, gradNorm;
	    // constants
	    const double untBeta;
	    const double successiveEps;
	    int parallelIterations;
	  public:
	    bool m_do_untangle_only;
	  };

	  static int count_invalid_elements(Mesh &mesh, MeshDomain *domain=0);

	  void run(Mesh &mesh, MeshDomain *domain, MsqError& err, bool always_smooth=true, int debug=0);

	};

	void ParShapeImprover::ParShapeImprovementWrapper::run_wrapper( MeshDomainAssoc* mesh_and_domain,
									ParallelMesh* pmesh,
									Settings* settings,
									QualityAssessor* qa,
									MsqError& err )
	{
	  int rank, nprocs;
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	  // Define an untangler
	  //UntangleBetaQualityMetric untangle_metric( untBeta );
	  UntangleBetaQualityMetric untangle_metric( 1.e-6 );

	  bool check_untangle = true;
	  if (check_untangle)
	    {
	      std::cout << "\nP[" << rank << "]  ParShapeImprover.... running QA with untangle_metric before... " << std::endl;
	      InstructionQueue q1;
	      QualityAssessor qa_untangle(&untangle_metric);
	      q1.add_quality_assessor(&qa_untangle, err); MSQ_ERRRTN(err);
	      q1.run_common( mesh_and_domain, pmesh, settings, err ); 
	      std::cout << "\nP[" << rank << "]  ParShapeImprover.... running QA with untangle_metric... before... done " << std::endl;
	    }

	  LPtoPTemplate untangle_func( 2, &untangle_metric );
	  ConjugateGradient untangle_solver( &untangle_func );
	  //untangle_solver.set_debugging_level(3);

	  //SteepestDescent untangle_solver( &untangle_func );
	  TerminationCriterion untangle_inner("<type:untangle_inner>"), untangle_outer("<type:untangle_outer>");
	  untangle_solver.use_global_patch();

	  //untangle_inner.add_absolute_gradient_L2_norm( gradNorm );
	  //untangle_inner.add_absolute_successive_improvement( successiveEps );
	  //untangle_inner.add_relative_successive_improvement( 1.e-6 );
	  //untangle_inner.add_untangled_mesh();

	  untangle_inner.write_iterations("untangle.gpt", err);

	  // For parallel runs, we generally need to have the inner and outer TerminationCriterion
	  //   have the same criteria else we can get an infinite loop (see VertexMover::loop_over_mesh)
	  untangle_inner.add_absolute_quality_improvement( 0.0 );
	  untangle_inner.add_iteration_limit( 20 );

	  untangle_outer.add_absolute_quality_improvement( 0.0 );
	  untangle_outer.add_iteration_limit( pmesh ? parallelIterations : 1 );

	  untangle_solver.set_inner_termination_criterion( &untangle_inner );
	  untangle_solver.set_outer_termination_criterion( &untangle_outer );

	  // define shape improver
	  IdealWeightInverseMeanRatio inverse_mean_ratio;
	  inverse_mean_ratio.set_averaging_method( QualityMetric::LINEAR );
	  LPtoPTemplate obj_func( 2, &inverse_mean_ratio );

	  ConjugateGradient shape_solver( &obj_func );
	  TerminationCriterion term_inner("<type:shape_inner>"), term_outer("<type:shape_outer>");
	  term_inner.write_iterations("shape.gpt", err);

	  shape_solver.use_global_patch();
	  qa->add_quality_assessment( &inverse_mean_ratio );

	  // For parallel runs, we generally need to have the inner and outer TerminationCriterion
	  //   have the same criteria else we can get an infinite loop (see VertexMover::loop_over_mesh)
	  term_inner.add_absolute_gradient_L2_norm( gradNorm );
	  term_inner.add_absolute_vertex_movement(0.0);
	  term_inner.add_iteration_limit( innerIter );

	  term_outer.add_absolute_gradient_L2_norm( gradNorm );
	  term_outer.add_absolute_vertex_movement(0.0);
	  term_outer.add_iteration_limit( pmesh ? parallelIterations : 1 );

	  //term_outer.add_absolute_quality_improvement( 1.e-6 );
	  //!term_outer.add_relative_successive_improvement( successiveEps );

	  shape_solver.set_inner_termination_criterion( &term_inner );
	  shape_solver.set_outer_termination_criterion( &term_outer );

	  // Apply CPU time limit to untangler
	  if (maxTime > 0.0)
	    untangle_inner.add_cpu_time( maxTime );
	  
	  Timer totalTimer;

	  // Run untangler
	  std::cout << "\nP[" << rank << "] " << " ParShapeImprovementWrapper: running untangler...\n " << std::endl;
	  bool use_untangle_wrapper = false;
	  if (use_untangle_wrapper)
	    {
	      UntangleWrapper uw;
	      //uw.set_untangle_metric(UntangleWrapper::BETA);
	      uw.run_instructions(mesh_and_domain, err);
	    }
	  else
	    {
	      InstructionQueue q1;
	      QualityAssessor qa_untangle(&untangle_metric);
	      q1.add_quality_assessor(&qa_untangle, err); MSQ_ERRRTN(err);
	      q1.set_master_quality_improver( &untangle_solver, err ); MSQ_ERRRTN(err);
	      q1.add_quality_assessor(&qa_untangle, err); MSQ_ERRRTN(err);
	      q1.run_common( mesh_and_domain, pmesh, settings, err ); 
	    }
	  std::cout << "\nP[" << rank << "] " << " ParShapeImprovementWrapper: running untangler... done\n " << std::endl;
	  std::cout << "\nP[" << rank << "] " << " ParShapeImprovementWrapper: MsqError after untangler: " << err << std::endl;

	  bool check_quality_after_untangler = true;
	  if (check_quality_after_untangler)
	    {
              Mesh* mesh = mesh_and_domain->get_mesh();
              MeshDomain* domain = mesh_and_domain->get_domain();
	      int num_invalid = count_invalid_elements(*mesh, domain);
	      std::cout << "\nP[" << rank << "] " << " ParShapeImprover num_invalid after untangler= " << num_invalid << " " 
			<< (num_invalid ? " ERROR still have invalid elements after Mesquite untangle" : 
			    " SUCCESS: untangled invalid elements ")
			<< std::endl;

	      if (check_untangle)
		{
		  std::cout << "\nP[" << rank << "]  ParShapeImprover.... running QA with untangle_metric " << std::endl;
		  InstructionQueue q1;
		  QualityAssessor qa_untangle(&untangle_metric);
		  q1.add_quality_assessor(&qa_untangle, err); MSQ_ERRRTN(err);
		  q1.run_common( mesh_and_domain, pmesh, settings, err ); 
		  std::cout << "\nP[" << rank << "]  ParShapeImprover.... running QA with untangle_metric... done " << std::endl;
		}

	      if (num_invalid) return;
	    }
	  if (m_do_untangle_only) return;
	  MSQ_ERRRTN(err);

	  // If limited by CPU time, limit next step to remaning time
	  if (maxTime > 0.0) {
	    double remaining = maxTime - totalTimer.since_birth();
	    if (remaining <= 0.0 ){
	      MSQ_DBGOUT(2) << "Optimization is terminating without perfoming shape improvement." << std::endl;
	      remaining = 0.0;
	    }
	    term_inner.add_cpu_time( remaining );
	  }
	  
	  // Run shape improver
	  InstructionQueue q2;
	  std::cout << "\nP[" << rank << "] " << " ParShapeImprovementWrapper: running shape improver... \n" << std::endl;
	  q2.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
	  q2.set_master_quality_improver( &shape_solver, err ); MSQ_ERRRTN(err);
	  q2.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
	  q2.run_common( mesh_and_domain, pmesh, settings, err ); 
	  std::cout << "\nP[" << rank << "] " << " ParShapeImprovementWrapper: running shape improver... done \n" << std::endl;
	  MSQ_ERRRTN(err);
	}


	int ParShapeImprover::count_invalid_elements(Mesh &mesh, MeshDomain *domain)
	{
	  MsqError err;
	  InstructionQueue q;
	  
	  IdealWeightInverseMeanRatio metric;
	  metric.set_averaging_method( QualityMetric::LINEAR );

	  // Check for inverted elements in the mesh
	  QualityAssessor inv_check( &metric );
	  //inv_check.disable_printing_results();
	  q.add_quality_assessor( &inv_check, err );  MSQ_ERRZERO(err);
	  Settings settings;
	  // bug?  should we pass in pmesh?
          Mesh* mesh_ptr = &mesh; 
          MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(mesh_ptr, domain);
	  q.run_common( &mesh_and_domain, 0, &settings, err ); MSQ_ERRZERO(err);
	  const QualityAssessor::Assessor* inv_b = inv_check.get_results( &metric );
	  int num_invalid = inv_b->get_invalid_element_count();
	  return num_invalid;
	}

	void ParShapeImprover::run(Mesh &mesh, MeshDomain *domain, MsqError& err, bool always_smooth, int debug)
	{
	  int rank, nprocs;
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	  if (debug)
	    {
	      MsqDebug::enable(1);
	      if (debug > 1) MsqDebug::enable(2);
	      if (debug > 2) MsqDebug::enable(3);
	    }

	  ParallelMesh *pmesh = dynamic_cast<ParallelMesh *>(&mesh);
	  std::cout << "P[" << rank << "] " << " ParShapeImprover::run: pmesh= " << pmesh << std::endl;

	  MsqError mErr;
	  int num_invalid = 0;
	  bool check_quality=true;
	  if (check_quality)
	    {
	      num_invalid = count_invalid_elements(mesh, domain);
	      std::cout << "\nP[" << rank << "] " << " ParShapeImprover num_invalid before= " << num_invalid 
			<< (num_invalid ? " WARNING: invalid elements exist before Mesquite smoothing" : 
			    (!always_smooth ? "WARNING: no smoothing requested since always_smooth=false" : " "))
			<< std::endl;
	    }

	  if (num_invalid || always_smooth)
	    {
	      bool use_canned_wrapper = false;
	      if (use_canned_wrapper)
		{
		  ShapeImprovementWrapper siw(mErr);
		  if (pmesh)
		    siw.run_instructions(pmesh, domain, mErr);
		  else
                  {
                    MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, domain);
		    siw.run_instructions(&mesh_and_domain, mErr);
                  }
		}
	      else
		{
		  int  msq_debug             = debug; // 1,2,3 for more debug info
		  bool always_smooth_local   = false;

		  bool do_untangle_only = false;
		  ParShapeImprover::ParShapeImprovementWrapper siw(innerIter,0.0,gradNorm,100);
		  siw.m_do_untangle_only = do_untangle_only;
		  if (pmesh)
		    siw.run_instructions(pmesh, domain, mErr);
		  else
                  {
                    MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, domain);
                    siw.run_instructions(&mesh_and_domain, mErr);
                  }

		}

	      std::cout << "\nP[" << rank << "] " << " ParShapeImprover: MsqError after ShapeImprovementWrapper: " << mErr << std::endl;

	      if (check_quality)
		{
		  num_invalid = count_invalid_elements(mesh, domain);
		  std::cout << "\nP[" << rank << "] " << " ParShapeImprover num_invalid after= " << num_invalid << " " 
			    << (num_invalid ? " ERROR still have invalid elements after Mesquite smoothing" : 
				" SUCCESS: smoothed and removed invalid elements ")
			    << std::endl;
		}

	      MSQ_ERRRTN(mErr);

	    }
	}

	static int test(std::string filename_prefix, std::string mesh_topology_name, MeshDomain *domain=0)
	{
	  int rank, nprocs;
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	  if (nprocs > 2) { cerr << "parallel_untangle_shape::test(" << mesh_topology_name << " can only be run with 1 or 2 processors" << std::endl; return 0; }

	  /* create processor-specific file names */
  ostringstream in_name, out_name, gold_name;
  in_name << filename_prefix << "par_untangle_original_" << mesh_topology_name << "_mesh." << nprocs << "." << rank << ".vtk";
  gold_name << filename_prefix << "par_untangle_smoothed_" << mesh_topology_name << "_mesh." << nprocs << "." << rank << ".vtk";
  out_name << "par_untangle_smoothed_" << mesh_topology_name << "_mesh." << nprocs << "." << rank << ".vtk";

  cout << "in_name= " << in_name.str() << " gold_name= " << gold_name.str() << " out_name= " << out_name.str() << std::endl;

  /* load different mesh files on each processor */
  MsqError err;
  MeshImpl mesh;
  mesh.read_vtk(in_name.str().c_str(), err);
  if (err) {cerr << err << endl; return 1;}

  /* create parallel mesh instance, specifying tags 
   * containing parallel data */
  ParallelMeshImpl parallel_mesh(&mesh, "GLOBAL_ID", "PROCESSOR_ID");
  ParallelHelperImpl helper;
  helper.set_communicator(MPI_COMM_WORLD);
  helper.set_parallel_mesh(&parallel_mesh);
  parallel_mesh.set_parallel_helper(&helper);

  /* do Laplacian smooth */
  //LaplaceWrapper optimizer;
  //optimizer.run_instructions(&parallel_mesh, err);

  int  msq_debug             = 0; // 1,2,3 for more debug info
  bool always_smooth         = true;
  int innerIter = 100;
  double gradNorm = 1.e-9;

  ParShapeImprover si(innerIter, gradNorm);
  //Mesh *pmesh = &parallel_mesh;
  si.run(parallel_mesh, domain, err, always_smooth, msq_debug);
  if (err) {cerr << err << endl; return 1; }

  /* write mesh */
  mesh.write_vtk(out_name.str().c_str(),err);
  if (err) {cerr << err << endl; return 1;}

  //std::cout << "P[ " << rank <<"] reading gold..." << std::endl;

  /* compare mesh with gold copy */
  MeshImpl gold;
  gold.read_vtk(gold_name.str().c_str(),err);
  if (err) {cerr << err << endl; return 1;}

  //std::cout << "P[ " << rank <<"] read gold, checking mesh diff..." << std::endl;
  bool do_print=true;
  double tol = 1.e-4;
  bool diff = MeshUtil::meshes_are_different(mesh, gold, err, tol, do_print);
  if (err) {cerr << err << endl; return 1;}
  //std::cout << "P[ " << rank <<"] read gold, checking mesh diff...done" << std::endl;

  if (diff) {cerr << "Error, computed mesh different from gold copy" << std::endl; return 1;}
  
  print_timing_diagnostics(cout);

  return 0;
}

int main( int argc, char* argv[] )
{
  /* init MPI */
  if (MPI_SUCCESS != MPI_Init(&argc, &argv)) {
    cerr << "MPI_Init failed." << endl;
    return 2;
  }

  Vector3D pnt(0,0,0);
  Vector3D s_norm(0,0,1);
  PlanarDomain msq_geom(s_norm, pnt);

  int t1 = test(VTK_2D_DIR, "quad", &msq_geom); 
  if (t1) return t1;
  t1 = test(VTK_3D_DIR, "hex");
  if (t1) return t1;

  MPI_Finalize();
  return 0;
}
