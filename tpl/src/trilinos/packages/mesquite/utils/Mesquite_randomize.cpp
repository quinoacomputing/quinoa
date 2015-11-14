#include "Mesquite_Randomize.hpp"
#include "Mesquite_InstructionQueue.hpp"
#include "Mesquite_QualityAssessor.hpp"
#include "Mesquite_MeshImpl.hpp"

#include "Mesquite_PatchData.hpp"
#include "Mesquite_MsqVertex.hpp"
#include "Mesquite_IdealWeightInverseMeanRatio.hpp"
#include "Mesquite_PMeanPTemplate.hpp"
#include "Mesquite_TerminationCriterion.hpp"
#include <assert.h>

#include "Mesquite_domain.hpp"

#include <iostream>
#include <iomanip>
#include <memory>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

using namespace Mesquite;


const char INVALID_FLAG = 'i';
const char PERCENT_FLAG = 'p';
const char UNOPTIMIZE_FLAG = 'u';



class UnOptimizer : public VertexMover
{
public:
  
  UnOptimizer( ObjectiveFunction* of ) : objectiveFunction(of) {}
  
  virtual ~UnOptimizer() {}
  
  virtual std::string get_name() const;
  
  virtual PatchSet* get_patch_set();
  
protected:

    virtual void initialize(PatchData &pd, MsqError &err);
    
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);
                                         
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    
    virtual void cleanup();

private:
  
    std::vector<size_t> adjVtxList;
    VertexPatches patchSet;
    ObjectiveFunction* objectiveFunction;
};

std::string UnOptimizer::get_name() const { return "UnOptimize"; }
PatchSet* UnOptimizer::get_patch_set() { return &patchSet; }
void UnOptimizer::initialize( PatchData&, MsqError& ) {}
void UnOptimizer::initialize_mesh_iteration(PatchData& , MsqError& ) {}
void UnOptimizer::terminate_mesh_iteration(PatchData& , MsqError& ) {}
void UnOptimizer::cleanup() {}

void UnOptimizer::optimize_vertex_positions( PatchData &pd, MsqError &err) {
  assert( pd.num_free_vertices() == 1 && pd.vertex_by_index(0).is_free_vertex() );
  std::vector<Vector3D> grad(1);
  double val, junk, coeff;
  bool state;
  
  state = objectiveFunction->evaluate_with_gradient( ObjectiveFunction::CALCULATE,
                                                     pd, val, grad, err );
  MSQ_ERRRTN(err);
  if (!state) {
    MSQ_SETERR(err)(MsqError::INVALID_MESH);
    return;
  }
  grad[0] /= grad[0].length();
  
  PatchDataVerticesMemento* memento = pd.create_vertices_memento( err ); MSQ_ERRRTN(err);
  std::auto_ptr<PatchDataVerticesMemento> deleter( memento );
  pd.get_minmax_edge_length( junk, coeff );
  
  for (int i = 0; i < 100; ++i) {
    pd.set_free_vertices_constrained( memento, arrptr(grad), 1, coeff, err ); MSQ_ERRRTN(err);
    state = objectiveFunction->evaluate( ObjectiveFunction::CALCULATE, pd, val, true, err );
    MSQ_ERRRTN(err);
    if (state)
      break;
    coeff *= 0.5;
  }
  if (!state) {
    pd.set_to_vertices_memento( memento, err );
  }
}


int main( int argc, char* argv[] )
{
  const double default_fraction = 0.05;
  const double zero = 0.0;
  int one = 1;
  CLArgs::ToggleArg allow_invalid( false );
  CLArgs::DoubleRangeArg rand_percent( default_fraction, &zero, 0 );
  CLArgs::IntRangeArg unoptimize( 0, &one, 0 );
  
  CLArgs args( "vtkrandom",
               "Randomize mesh vertex locations.",
               "Read VTK file, randomize locations of containded vertices, and re-write file." );
  args.toggle_flag( INVALID_FLAG, "Allow inverted elements in output", &allow_invalid );
  args.double_flag( PERCENT_FLAG, "fract", "Randomize fraction", &rand_percent );
  args.int_flag( UNOPTIMIZE_FLAG, "N", "Use UnOptimizer with N passes rather than Randomize", &unoptimize );
  add_domain_args( args );
  args.add_required_arg( "input_file" );
  args.add_required_arg( "output_file" );

  std::vector<std::string> files;
  if (!args.parse_options( argc, argv, files, std::cerr )) {
    args.print_usage( std::cerr );
    exit(1);
  }
  std::string input_file = files[0];
  std::string output_file = files[1];
  
  MsqError err;
  MeshImpl mesh;
  mesh.read_vtk( input_file.c_str(), err );
  if (err) {
    std::cerr << "ERROR READING FILE: " << input_file << std::endl
                    << err << std::endl;
    return 2;
  }
  MeshDomain* domain = process_domain_args( &mesh );

  TerminationCriterion tc;
  QualityAssessor qa( false );
  InstructionQueue q;
  Randomize op( rand_percent.value() );
  IdealWeightInverseMeanRatio metric;
  PMeanPTemplate of( 1, &metric );
  UnOptimizer op2( &of );
  if (unoptimize.seen()) {
    tc.add_iteration_limit( unoptimize.value() );
    op2.set_outer_termination_criterion( &tc );
    q.add_preconditioner( &op, err );
    q.set_master_quality_improver( &op2, err );
  }
  else {
    q.set_master_quality_improver( &op, err );
  }
  q.add_quality_assessor( &qa, err );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, domain);
  q.run_instructions( &mesh_and_domain, err );
  if (err) {
    std::cerr << err << std::endl;
    return 3;
  }

  int inverted, junk;
  if (qa.get_inverted_element_count( inverted, junk, err ) && inverted ) {
    if (allow_invalid.value())
      std::cerr << "Warning: output mesh contains " << inverted << " inverted elements" << std::endl;
    else {
      std::cerr << "Error: output mesh contains " << inverted << " inverted elements" << std::endl;
      return 4;
    }
  }
  
  mesh.write_vtk( output_file.c_str(), err );
  if (err) {
    std::cerr << "ERROR WRITING FILE: " << output_file << std::endl
                    << err << std::endl;
    return 2;
  }
  
  return 0;
}

