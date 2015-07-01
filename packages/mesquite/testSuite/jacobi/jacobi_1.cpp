/**\file jacobi_1.cpp
  *\brief simple convergence test for Jacobi optimization
  *\author Jason Kraftcheck <kraftche@cae.wisc.edu>
  *
  *Create simple quad mesh and test for convergence.
  */

#include "ArrayMesh.hpp"
#include "PlanarDomain.hpp"
#include "IdealShapeTarget.hpp"
#include "TShapeB1.hpp"
#include "TQualityMetric.hpp"
#include "PMeanPTemplate.hpp"
#include "SteepestDescent.hpp"
#include "TerminationCriterion.hpp"
#include "InstructionQueue.hpp"
#include "MsqError.hpp"
#include "QualityAssessor.hpp"
#include "MeshWriter.hpp"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

using namespace Mesquite;

/* Mesh:
 
  12---13---14---15
  |    |    |    | 
  |    |    |    |
  8----9----10---11
  |    |    |    | 
  |    |    |    |
  4----5----6----7  
  |    |    |    | 
  |    |    |    |
  0----1----2----3

*/

double dist( double* p1, double* p2 ) {
  double d1 = p1[0] - p2[0];
  double d2 = p1[1] - p2[1];
  double d3 = p1[2] - p2[2];
  return sqrt(d1*d2*d3);
}

void usage( const char* argv0 )
{
  fprintf(stderr,"Usage: %s [-i <init_vtk_file>] [-f <final_vtk_file>]\n", argv0);
  exit(1);
}

int main(int argc, char* argv[])
{
  const char* initial_mesh_file = 0;
  const char* final_mesh_file = 0;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i],"-f")) {
      ++i;
      if (i == argc)
        usage(argv[0]);
      final_mesh_file = argv[i];
    }
    else if (!strcmp(argv[i],"-i")) {
      ++i;
      if (i == argc)
        usage(argv[0]);
      initial_mesh_file = argv[i];
    }
    else
      usage( argv[i] );
  }


    // create mesh
  const int intervals = 3;
  const double perturb = 0.3;
  const int nvtx = (intervals+1) * (intervals+1);
  double coords[nvtx*3];
  double exp_coords[nvtx*3];
  int fixed[nvtx];
  for (int i = 0; i < nvtx; ++i) {
    double* c = coords + 3*i;
    double* e = exp_coords + 3*i;
    int row = i / (intervals+1);
    int col = i % (intervals+1);
    double xdelta, ydelta;
    if (row > 0 && row < intervals && col > 0 && col < intervals) {
      fixed[i] = 0;
      xdelta = row % 2 ? -perturb : 0;
      ydelta = col % 2 ? perturb : -perturb;
    }
    else {
      fixed[i] = 1;
      xdelta = ydelta = 0.0;
    }
    c[0] = col + xdelta;
    c[1] = row + ydelta;
    c[2] = 0.0;
    e[0] = col;
    e[1] = row;
    e[2] = 0.0;
  }
  const int nquad = intervals * intervals;
  unsigned long conn[nquad * 4];
  for (int i = 0; i < nquad; ++i) {
    unsigned long* c = conn + 4*i;
    int row = i / intervals;
    int col = i % intervals;
    int n0 = (intervals+1)*row + col;
    c[0] = n0;
    c[1] = n0+1;
    c[2] = n0+intervals+2;
    c[3] = n0+intervals+1;
  }
  MsqPrintError err(std::cerr);
  ArrayMesh mesh( 3, nvtx, coords, fixed, nquad, QUADRILATERAL, conn );
  PlanarDomain zplane( PlanarDomain::XY );
  if (initial_mesh_file) {
    MeshWriter::write_vtk( &mesh, initial_mesh_file, err );
    if (err) {
      fprintf(stderr,"%s: failed to write file\n", initial_mesh_file );
      return 1;
    }
  }
  
    // do optimization
  const double eps = 0.01;
  IdealShapeTarget w;
  TShapeB1 mu;
  TQualityMetric metric( &w, &mu );
  PMeanPTemplate func( 1.0, &metric );
  SteepestDescent solver( &func );
  solver.use_element_on_vertex_patch();
  solver.do_jacobi_optimization();
  TerminationCriterion inner, outer;
  inner.add_absolute_vertex_movement( 0.5*eps );
  outer.add_absolute_vertex_movement( 0.5*eps );
  QualityAssessor qa( &metric );
  InstructionQueue q;
  q.add_quality_assessor( &qa, err );
  q.set_master_quality_improver( &solver, err );
  q.add_quality_assessor( &qa, err );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &zplane);
  q.run_instructions( &mesh_and_domain, err );
  if (err)
    return 2;
  if (final_mesh_file) {
    MeshWriter::write_vtk( &mesh, final_mesh_file, err );
    if (err) {
      fprintf(stderr,"%s: failed to write file\n", final_mesh_file );
      return 1;
    }
  }
  
    // check final vertex positions
  int invalid = 0;
  for (int i = 0; i < nvtx; ++i) {
    if (dist( coords + 3*i, exp_coords + 3*i ) > eps) {
      ++invalid;
      printf("Vertex %d at (%f,%f,%f), expected at (%f,%f,%f)\n", i,
        coords[3*i], coords[3*i+1], coords[3*i+2],
        exp_coords[3*i], exp_coords[3*i+1], exp_coords[3*i+2] );
    }
  }
  
  return invalid ? 2 : 0;
}
