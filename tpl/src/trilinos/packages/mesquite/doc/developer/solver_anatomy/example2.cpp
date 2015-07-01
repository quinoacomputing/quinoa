#include "TerminationCriterion.hpp"
#include "FeasibleNewton.hpp"
#include "PMeanPTemplate.hpp"
#include "JacobianMetric.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "ReferenceMesh.hpp"
#include "MeshImpl.hpp"
#include "UnitWeight.hpp"
#include "Target2DShape.hpp"
#include "SamplePoints.hpp"
#include "LinearFunctionSet.hpp"

using namespace Mesquite;

int main()
{
  SamplePoints samples(true);
  Target2DShape metric_2d;
  UnitWeight weight;
  
  MsqError err;
  MeshImpl ref_mesh_data;
  ref_mesh_data.read_vtk( "refmesh.vtk", err );
  if (err) {
    std::cerr << "Error reading file: \"refmesh.vtk\""
              << std::endl << err << std::endl;
    exit( 1 );
  }
  
  ReferenceMesh ref_mesh( &ref_mesh_data );
  RefMeshTargetCalculator target( &ref_mesh );
  
  JacobianMetric metric( &samples,
                         &target,
                         &weight,
                         &metric_2d,
                         NULL );
  
  PMeanPTemplate obj_func( 2.0, &metric );
  
  TerminationCriterion inner_term, outer_term;
  inner_term.add_criterion_type_with_double( 
TerminationCriterion::VERTEX_MOVEMENT_ABSOLUTE, 0.01, err );
  inner_term.add_criterion_type_with_int( 
TerminationCriterion::CPU_TIME, 60, err );
  outer_term.add_criterion_type_with_int( 
TerminationCriterion::NUMBER_OF_ITERATES, 1, err );

  FeasibleNewton solver( &obj_func );
  solver.set_inner_termination_criterion( &inner_term );
  solver.set_outer_termination_criterion( &outer_term );
  
  solver.use_global_patch();
  
  LinearFunctionSet mapping_functions;
  
  return 0;
}
