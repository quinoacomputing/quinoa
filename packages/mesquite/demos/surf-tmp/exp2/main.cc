#include "mesquite_version.h"

#include "TMPQualityMetric.hpp"
#include "TargetCalculator.hpp"
#include "ReferenceMesh.hpp"

#include "Target2DShapeSizeOrientBarrier.hpp"
typedef Mesquite::Target2DShapeSizeOrientBarrier TMetric;
#define TMetric_NAME "|T-I|^2"

#include "QualityAssessor.hpp"

#include "XYRectangle.hpp"
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
#include "TagVertexMesh.hpp"

#include "PatchData.hpp"

#include <iostream>

using namespace Mesquite;

const int WIDTH = 10;
const int HEIGHT = 7;

class AlignTC : public TargetCalculator {
private:
  ReferenceMeshInterface* initMesh;
  int normalAxis;
public:
  AlignTC( ReferenceMeshInterface* init_mesh, int normal ) 
    : initMesh(init_mesh), normalAxis(normal) {}
  bool get_3D_target( PatchData&, size_t, Sample, MsqMatrix<3,3>&, MsqError& err )
    { MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT); return false; }
  bool get_surface_target( PatchData&, size_t, Sample, MsqMatrix<3,2>&, MsqError& );
  bool get_2D_target( PatchData&, size_t, Sample, MsqMatrix<2,2>&, MsqError& err )
    { MSQ_SETERR(err)(MsqError::INVALID_STATE); return false; }
  bool have_surface_orient() const 
    { return true; }
    
  static MsqVector<3> calc_b( int normal_idx, const Vector3D& coord )
  {
    int width_idx = (normal_idx+1)%3;
    int height_idx = (normal_idx+2)%3;
    MsqVector<3> b;
    b[width_idx] = 20*(coord[width_idx] - 5);
    b[height_idx] = 14*(coord[height_idx] + 5);
    b[normal_idx] = 0;
    b /= length(b);
    return b;
  }
};

void usage( const char* argv0 )
{
  std::cerr << "Usage: " << argv0 << " <input_file> [-x|-y|-z] [-o <output_file>] [-i <init_file>] [-p <plot_file>] [-t <timestep_base>] " << std::endl;
  std::cerr << "Input mesh must be a mesh of a rectangle in the XY plane where X is in [0,10] and Y is in [0,7]" << std::endl;
  exit(1);
}

void rotate_mesh( int normal_axis, Mesh& mesh );
void tag_vertices_with_field_vects( Mesh& mesh, int normal_axis );

int main( int argc, char* argv[] )
{
  XYRectangle::Plane normal_axis = XYRectangle::XY; // default to Z
  const char *input_file = 0, *output_file = 0, *init_file = 0, 
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
        case 'x': normal_axis = XYRectangle::YZ; break;
        case 'y': normal_axis = XYRectangle::ZX; break;
        case 'z': normal_axis = XYRectangle::XY; break;
        default: usage(argv[0]);
      }
    }
    else if (!input_file)
      input_file = argv[i];
    else 
      usage(argv[0]);
  }
  if (arg || !input_file) usage(argv[0]);

  MsqPrintError err(std::cerr);
  MeshImpl mesh;
  mesh.read_vtk( input_file, err );
  if (err) {
    std::cerr << input_file << " : failed to read input file" << std::endl;
    return 2;
  }

  if (normal_axis != XYRectangle::XY)
    rotate_mesh( normal_axis, mesh );

  XYRectangle plane( WIDTH, HEIGHT, 0, 0, 0, normal_axis );
  MeshDomain* geom = &plane;
  plane.setup( &mesh, err );
  if (err) {
    std::cerr << "Domain setup failed" << std::endl;
    return 3;
  }

  tag_vertices_with_field_vects( mesh, normal_axis );

  TagVertexMesh init_mesh( err, &mesh );
  ReferenceMesh ref_mesh(&init_mesh);
  AlignTC tc(&ref_mesh, normal_axis);
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
  q1.add_tag_vertex_mesh( &init_mesh, err );
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


bool AlignTC::get_surface_target( PatchData& pd, size_t element, Sample sample, 
                                  MsqMatrix<3,2>& W, MsqError& err )
{
  MsqMatrix<3,2> J;
  get_refmesh_Jacobian_2D( initMesh, pd, element, sample, J, err );
  MSQ_ERRZERO(err);
  
  MsqMeshEntity& elem = pd.element_by_index( element );
  const EntityTopology type = elem.get_element_type();
    
  assert(sample.dimension == 0);
  Mesh::VertexHandle vtx_h = pd.get_vertex_handles_array()[elem.get_vertex_index_array()[sample.number]];
  Vector3D coord;
  initMesh->get_reference_vertex_coordinates( &vtx_h, 1, &coord, err );
  MSQ_ERRZERO(err);
  
  MsqVector<3> b, u(0.0);
  u[normalAxis] = 1.0;
  b = calc_b( normalAxis, coord );
  
  MsqVector<3> b_s = (b - u * (b % u))/length(b * u);
  MsqMatrix<3,2> V;
  V.set_column( 0, b_s );
  V.set_column( 1, u * b_s );
  
  W = size(J) * V * shape(J);
  return true;  
}

void tag_vertices_with_field_vects( Mesh& mesh, int normal_axis )
{
  MsqPrintError err(std::cerr);
  std::vector<Mesh::VertexHandle> handles;
  mesh.get_all_vertices( handles, err );
  std::vector<MsqVertex> coords(handles.size());
  std::vector< MsqVector<3> > b(handles.size());
  mesh.vertices_get_coordinates( &handles[0], &coords[0], handles.size(), err );
  for (size_t i = 0; i < coords.size(); ++i) 
    b[i] = AlignTC::calc_b( normal_axis, coords[i] );
  
  TagHandle t = mesh.tag_create( "Vector Field", Mesh::DOUBLE, 3, 0, err );
  mesh.tag_set_vertex_data( t, handles.size(), &handles[0], &b[0], err );
}

void rotate_mesh( int normal_axis, Mesh& mesh )
{
  MsqPrintError err(std::cerr);
  std::vector<Mesh::VertexHandle> handles;
  mesh.get_all_vertices( handles, err );
  std::vector<MsqVertex> coords(handles.size());
  mesh.vertices_get_coordinates( &handles[0], &coords[0], handles.size(), err );
  for (size_t i = 0; i < handles.size(); ++i) {
    Vector3D c;
    
    if (normal_axis == XYRectangle::YZ) {
      c[0] = coords[i][2];
      c[1] = coords[i][0];
      c[2] = coords[i][1];
    }
    else if (normal_axis == XYRectangle::ZX) {
      c[0] = coords[i][1];
      c[1] = coords[i][2];
      c[2] = coords[i][0];
    }
    else {
      c = coords[i];
    }
    mesh.vertex_set_coordinates( handles[i], c, err );
  }
}
