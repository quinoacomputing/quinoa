#include "mesquite_version.h"

#include "TMPQualityMetric.hpp"
#include "ReferenceMesh.hpp"
#include "RefMeshTargetCalculator.hpp"

#include "Target2DShapeSizeOrient.hpp"
typedef Mesquite::Target2DShapeSizeOrient TMetric;
#define TMetric_NAME "ShapeSizeOrient"

#include "QualityAssessor.hpp"

#include "PlanarDomain.hpp"
#include "SphericalDomain.hpp"
#include "CylinderDomain.hpp"
#include "ConicDomain.hpp"
#include "InstructionQueue.hpp"

#include "SteepestDescent.hpp"
typedef Mesquite::SteepestDescent SolverType;

#include "TerminationCriterion.hpp"
#include "PMeanPTemplate.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "IdealWeightInverseMeanRatio.hpp"
#include "ElementPMeanP.hpp"
#include "MsqTimer.hpp"

#include <iostream>

using namespace Mesquite;

const double EPS = 1e-3;

bool find_sphere( Mesh* mesh, SphericalDomain* dom, MsqError& err );
bool find_cyl( Mesh* mesh, CylinderDomain* dom, MsqError& err );
bool find_cone( Mesh* mesh, ConicDomain* dom, MsqError& err );

void usage( const char* argv0 )
{
  std::cerr << "Usage: " << argv0 << " <input_file> <reference_mesh> [-o <output_file>] [-i <init_file>] [-p <plot_file>] [-t <timestep_base>] " << std::endl;
  exit(1);
}

int main( int argc, char* argv[] )
{
  const char *input_file = 0, *reference_file = 0,*output_file = 0, *init_file = 0, 
             *plot_file = 0, *timestep_base = 0;
  const char** arg = 0;
  for (int i = 1; i < argc; ++i) {
    if (arg) {
      if (*arg) {
        std::cerr << "File specified multiple times" <<std::endl;
        usage(argv[0]);
      }
      *arg = argv[i];
      arg = 0;
    } 
    else if (argv[i][0] == '-') {
      if (!argv[i][1] || argv[i][2])
        usage(argv[0]);
      switch (argv[i][1]) {
        case 'o': arg = &output_file;   break;
        case 'i': arg = &init_file;     break;
        case 'p': arg = &plot_file;     break;
        case 't': arg = &timestep_base; break;
        default: usage(argv[0]);
      }
    }
    else if (!input_file)
      input_file = argv[i];
    else if (!reference_file)
      reference_file = argv[i];
    else 
      usage(argv[0]);
  }
  if (arg || !input_file || !reference_file) usage(argv[0]);

  MsqPrintError err(std::cerr);
  MeshImpl mesh;
  mesh.read_vtk( input_file, err );
  if (err) {
    std::cerr << input_file << " : failed to read input file" << std::endl;
    return 2;
  }

  MeshImpl rmesh;
  rmesh.read_vtk( reference_file, err );
  if (err) {
    std::cerr << reference_file << " : failed to read reference mesh file" << std::endl;
    return 2;
  }
  
  SphericalDomain sphere( Vector3D(0,0,0), 0.0 );
  CylinderDomain cylinder( 0, Vector3D(1,0,0), Vector3D(0,0,0) );
  ConicDomain cone;
  MeshDomain* geom = 0;
  if (find_sphere(&mesh, &sphere, err)) {
    geom = &sphere;
    std::cout << "Using spherical domain" << std::endl;
  }
  else if (err) 
    return 3;
  else if (find_cyl(&mesh, &cylinder, err)) {
    geom = &cylinder;
    std::cout << "Using cylinderical domain" << std::endl;
    mesh.mark_skin_fixed( err );
    if (err) {
      std::cerr << "Skinning failed" << std::endl;
      return 3;
    }
  }
  else if (err)
    return 3;
  else if (find_cone(&mesh, &cone, err)) {
    geom = &cone;
    std::cout << "Using conic domain" << std::endl;
    mesh.mark_skin_fixed( err );
    if (err) {
      std::cerr << "Skinning failed" << std::endl;
      return 3;
    }
  }
  else if (err)
    return 3;
  else {
    std::cerr << "Failed to identify mesh domain" << std::endl;
    return 2;
  }

  ReferenceMesh ref_mesh(&rmesh);
  RefMeshTargetCalculator tc(&ref_mesh);;

  TMetric mu;
  TMPQualityMetric qm( &tc, &mu, 0 );
  ElementPMeanP avg_qm( 1.0, &qm );

  PMeanPTemplate of(1.0, &qm);
  SolverType qi(&of);
  TerminationCriterion inner;
  inner.add_absolute_vertex_movement( 1e-4 );
  if (plot_file) 
    inner.write_iterations( plot_file, err );
  if (timestep_base)
    inner.write_mesh_steps( timestep_base );
  
  qi.set_inner_termination_criterion(&inner);
  IdealWeightInverseMeanRatio imr;
  QualityAssessor qa(&qm, 0, 0.0, false, TMetric_NAME "_init", true, "Inverted_init");
  QualityAssessor qa2(&qm, 0, 0.0, false, TMetric_NAME, true, "Inverted");
  qa.add_quality_assessment(&avg_qm, 0, 0, "Average " TMetric_NAME "_init");
  qa2.add_quality_assessment(&avg_qm, 0, 0, "Average " TMetric_NAME );
  qa.add_quality_assessment(&imr, 0, 0, "InverseMeanRatio_init");
  qa2.add_quality_assessment(&imr, 0, 0, "InverseMeanRatio");
  InstructionQueue q1;
  q1.add_quality_assessor(&qa,err);
  //q1.run_instructions( &mesh, geom, err );
  if (init_file) {
    q1.run_instructions( &mesh, geom, err );
    if (err) {
      std::cerr << "Assessment failed!!!" << std::endl;
      return 3;
    }
    mesh.write_vtk( init_file, err );
    if (err) {
      std::cerr << init_file << " : failed to write output file" << std::endl;
      return 2;
    }
    std::cout << "Wrote initial mesh to file: " << init_file << std::endl;
  }

  InstructionQueue q2;
  q2.set_master_quality_improver(&qi,err);
  Timer timer;
  q2.run_instructions( &mesh, geom, err );
  if (err) {
    std::cerr << "Optimization failed!!!" << std::endl;
    return 3;
  }
  std::cout << "InstructionQueue returned after " << timer.since_birth() << " seconds" << std::endl;
  
  InstructionQueue q3;
  q3.add_quality_assessor(&qa2,err);
  q3.run_instructions( &mesh, geom, err );
  
  if (output_file) {
    mesh.write_vtk( output_file, err );
    if (err) {
      std::cerr << output_file << " : failed to write output file" << std::endl;
      return 2;
    }
    std::cout << "Wrote final mesh to file: " << output_file << std::endl;
  }
  
  return 0;
}

bool find_sphere( Mesh* mesh, SphericalDomain* dom, MsqError& err )
{
  std::vector<Mesh::VertexHandle> verts;
  mesh->get_all_vertices( verts, err );
  if (err) return false;
  
  std::vector<MsqVertex> coords(verts.size());
  mesh->vertices_get_coordinates( &verts[0], &coords[0], verts.size(), err );
  if (err) return false;
  
  Vector3D min(coords[0]), max(coords[0]);
  for (size_t i = 1; i < coords.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      if (coords[i][j] < min[j])
        min[j] = coords[i][j];
      if (coords[i][j] > max[j])
        max[j] = coords[i][j];
    }
  }
  
  max *= 0.5;
  min *= 0.5;
  for (int j = 0; j < 3; ++j) {
    max[j] = round(20*max[j])/20;
    min[j] = round(20*min[j])/20;
  }
  Vector3D diff = max - min;
  Vector3D cent = max + min;
  double rad = diff[0];
  bool ok = true;
  for (size_t i = 0; i < coords.size(); ++i) {
    diff = coords[i]-cent;
    if (fabs(diff.length() - rad) > EPS) {
      //std::cerr << "Vertex " << i << " is not on mean sphere" << std::endl;
      //std::cerr << "  cent: " << cent << ", rad: " << rad << ", coords: " << coords[i] << std::endl;
      ok = false;
    }
  }

 dom->set_sphere( cent, rad );
 return ok;
}
  
bool find_cyl( Mesh* mesh, CylinderDomain* dom, MsqError& err )
{
  std::vector<Mesh::VertexHandle> verts;
  mesh->get_all_vertices( verts, err );
  if (err) return false;
  
  std::vector<MsqVertex> coords(verts.size());
  mesh->vertices_get_coordinates( &verts[0], &coords[0], verts.size(), err );
  if (err) return false;
  
  Vector3D min(coords[0]), max(coords[0]);
  for (size_t i = 1; i < coords.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      if (coords[i][j] < min[j])
        min[j] = coords[i][j];
      if (coords[i][j] > max[j])
        max[j] = coords[i][j];
    }
  }
  
  Vector3D mid = 0.5*(min+max);
  min.set(HUGE_VAL,HUGE_VAL,HUGE_VAL);
  max.set(0,0,0);
  for (size_t i = 0; i < coords.size(); ++i) {
    Vector3D diff = mid - coords[i];
    for (int j = 0; j < 3; ++j) {
      double x1 = diff[(j+1)%3];
      double x2 = diff[(j+2)%3];
      double dist_sqr = x1*x1 + x2*x2;
      if (min[j] > dist_sqr)
        min[j] = dist_sqr;
      if (max[j] < dist_sqr)
        max[j] = dist_sqr;
    }
  }
  
  max -= min;
  int axis = (max[0] < max[1]) && (max[0] < max[2]) ? 0 : (max[1] < max[2]) ? 1 : 2;
  double rad = round(20*(sqrt(min[axis]+0.5*max[axis])))/20;
  for (int j = 0; j < 3; ++j) 
    mid[j] = round(20*mid[j])/20;
  
  bool ok = true;
  const int a = (axis+1)%3, b = (axis+2)%3;
  for (size_t i = 0; i < coords.size(); ++i) {
    Vector3D diff = mid - coords[i];
    for (int j = 0; j < 3; ++j) {
      double dist = sqrt(diff[a]*diff[a]+diff[b]*diff[b]);
      if (fabs(dist - rad) > EPS)
        ok = false;
    }
  }
  
  Vector3D axis_vec(0,0,0);
  axis_vec[axis] = 1.0;
  *dom = CylinderDomain( rad, axis_vec, mid );
  return ok;
}
  
bool find_cone( Mesh* mesh, ConicDomain* dom, MsqError& err )
{
  ConicDomain cone( 2, 10 );

  std::vector<Mesh::VertexHandle> verts;
  mesh->get_all_vertices( verts, err );
  if (err) return false;
  
  std::vector<MsqVertex> coords(verts.size());
  mesh->vertices_get_coordinates( &verts[0], &coords[0], verts.size(), err );
  if (err) return false;
  
  for (size_t i = 0; i < coords.size(); ++i) {
    Vector3D p(coords[i]);
    cone.snap_to( 0, p );
    if ((p - coords[i]).length_squared() > 1e-5)
      return false;
  }
  
  *dom = cone;
  return true;
}
