#include "mesquite_version.h"

#include "TMPQualityMetric.hpp"
#include "ReferenceMesh.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "TagVertexMesh.hpp"

#include "Target2DShapeSizeOrient.hpp"
#include "Target2DShape.hpp"
typedef Mesquite::Target2DShape TMetric;
#define TMetric_NAME "Shape"

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
#include <set>
#include <iterator>
#include <algorithm>

using namespace Mesquite;

const double EPS = 1e-3;

// figure out which vertices in the input mesh correspond
// to which in the refernce mesh and store the coordinates of
// the refernece mesh vertex in a tag on each input mesh vertex.
//
// NOTE: this makes some assumptions about the input.  It probably
//       won't work w/ any other input mesh.
void match_reference_mesh( Mesh& input_mesh,
                           Mesh& ref_mesh,
                           const char* vertex_tag_name );

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
  mesh.mark_skin_fixed( err );
  assert(!err);

  MeshImpl rmesh;
  rmesh.read_vtk( reference_file, err );
  if (err) {
    std::cerr << reference_file << " : failed to read reference mesh file" << std::endl;
    return 2;
  }
  
  PlanarDomain plane( PlanarDomain::XZ, 10 );
  MeshDomain* geom = &plane;
  
  const char* REF_COORD_TAG = "reference coords";
  match_reference_mesh( mesh, rmesh, REF_COORD_TAG );
  //mesh.write_vtk( "debug.vtk", err );
  //return 0;
  
  TagVertexMesh tmesh( err, &mesh, false, REF_COORD_TAG );
  ReferenceMesh ref_mesh(&tmesh);
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
  q1.run_instructions( &mesh, geom, err );
  if (err) {
    std::cerr << "Assessment failed!!!" << std::endl;
    return 3;
  }
  if (init_file) {
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

class XSort {
  public:
    XSort( Mesh& mesh ) : mMesh(mesh) {}
    bool operator()( Mesh::VertexHandle v1, Mesh::VertexHandle v2 ) {
      MsqError err;
      Mesh::VertexHandle handles[] = {v1, v2};
      MsqVertex coords[2];
      mMesh.vertices_get_coordinates( handles, coords, 2, err );
      assert(!err);
      return coords[0][0] < coords[1][0];
    }
    Mesh& mMesh;
};

// figure out which vertices in the input mesh correspond
// to which in the refernce mesh and store the coordinates of
// the refernece mesh vertex in a tag on each input mesh vertex.
void match_reference_mesh( Mesh& mesh,
                           Mesh& ref_mesh,
                           const char* vertex_tag_name )
{
  double EPSILON = 1e-3;

  MsqPrintError err(std::cerr);
  std::vector<Mesh::VertexHandle> verts, ref_verts;
  mesh.get_all_vertices( verts, err );
  ref_mesh.get_all_vertices( ref_verts, err );
  if (verts.size() != ref_verts.size()) {
    std::cerr << "Input and reference mesh incompatible." << std::endl;
    exit(1);
  }
  std::vector<MsqVertex> coords(verts.size()), ref_coords(ref_verts.size());
  mesh.vertices_get_coordinates( &verts[0], &coords[0], verts.size(), err );
  ref_mesh.vertices_get_coordinates( &ref_verts[0], &ref_coords[0], ref_verts.size(), err );
  
  Mesh::VertexHandle zero = 0;
  TagHandle match_tag = mesh.tag_create( "match", Mesh::HANDLE, 1, &zero, err );
  
    // Assume that there exist two lines of nodes, one at the maximum
    // z value of the target surface and one at the maximum y value of
    // the source surface, for which the nodes map by x coordinate.
  
  double inp_max_z = -HUGE_VAL, ref_max_y = -HUGE_VAL;
  for (size_t i = 0; i < coords.size(); ++i) {
    if (coords[i][2] > inp_max_z)
      inp_max_z = coords[i][2];
    if (ref_coords[i][1] > ref_max_y)
      ref_max_y = ref_coords[i][1];
  }
  
  std::vector<Mesh::VertexHandle> inp_line, ref_line;
  for (size_t i = 0; i < coords.size(); ++i) {
    if (coords[i][2] > inp_max_z - EPSILON)
      inp_line.push_back( verts[i] );
    if (ref_coords[i][1] > ref_max_y - EPSILON)
      ref_line.push_back( ref_verts[i] );
  }
  if (inp_line.size() != ref_line.size()) {
    std::cerr << "Input and reference mesh incompatible." << std::endl;
    exit(1);
  }
  
  std::sort( inp_line.begin(), inp_line.end(), XSort(mesh) );
  std::sort( ref_line.begin(), ref_line.end(), XSort(ref_mesh) );
  
  mesh.tag_set_vertex_data( match_tag, inp_line.size(), &inp_line[0], &ref_line[0], err );
  assert(!err);
  
  std::vector<Mesh::VertexHandle> next_line, conn, matchv, matchp, ref_conn;
  std::vector<Mesh::ElementHandle> elems, ref_elems, tmp;
  std::vector<size_t> offsets;
  std::set<Mesh::ElementHandle> ref_matched, inp_matched;
  while (!inp_line.empty()) {
    for (size_t i = 0; i < inp_line.size(); ++i) {
      bool all_matched = true;
      offsets.clear();
      elems.clear();
      mesh.vertices_get_attached_elements( &inp_line[i], 1, elems, offsets, err );
      assert(!err);
      
      for (size_t j = 0; j < elems.size(); ++j) {
        if (inp_matched.find( elems[j] ) != inp_matched.end())
          continue;
        
        conn.clear();
        offsets.clear();
        mesh.elements_get_attached_vertices( &elems[j], 1, conn, offsets, err );
        assert(!err);
        
        assert(conn.size() == 4); // assume linear quads
        matchv.resize( conn.size() );
        mesh.tag_get_vertex_data( match_tag, conn.size(), &conn[0], &matchv[0], err );
        assert(!err);
        
        matchp.clear();
        for (size_t k = 0; k < matchv.size(); ++k)
          if (matchv[k])
            matchp.push_back( matchv[k] );
        if (matchp.size() < 2) {
          all_matched = false;
          continue;
        }
        
        tmp.clear();
        offsets.clear();
        ref_mesh.vertices_get_attached_elements( &matchp[0], matchp.size(), tmp, offsets, err );
        ref_elems.clear();
        std::sort( tmp.begin(), tmp.begin() + offsets[1] );
        std::sort( tmp.begin() + offsets[1], tmp.begin() + offsets[2] );
        std::set_intersection( tmp.begin(), tmp.begin() + offsets[1],
                               tmp.begin() + offsets[1],
                               tmp.begin() + offsets[2], 
                               std::back_inserter( ref_elems ) );
        for (size_t k = 2; k < matchp.size(); ++k) {
          offsets[0] = ref_elems.size();
          std::copy( ref_elems.begin(), ref_elems.end(), tmp.begin() );
          ref_elems.clear();
          std::sort( tmp.begin() + offsets[k], tmp.begin() + offsets[k+1] );
          std::set_intersection( tmp.begin(), tmp.begin() + offsets[0],
                                 tmp.begin() + offsets[k],
                                 tmp.begin() + offsets[k+1], 
                                 std::back_inserter( ref_elems ) );
        }
        
        tmp.swap( ref_elems );
        ref_elems.clear();
        std::set_difference( tmp.begin(), tmp.end(), 
                             ref_matched.begin(), ref_matched.end(), 
                             std::back_inserter( ref_elems ) );
        
        bool found_match = false;
        for (size_t k = 0; k < ref_elems.size() && !found_match; ++k) {
          offsets.clear();
          ref_conn.clear();
          ref_mesh.elements_get_attached_vertices( &ref_elems[k], 1, ref_conn, offsets, err );
          bool match_forward = true;
          matchp.clear();
          assert(ref_conn.size() == 4);
          int off;
          for (off = 0; !matchv[off]; ++off);
          int ref = std::find( ref_conn.begin(), ref_conn.end(), matchv[off] ) - ref_conn.begin();
          assert(ref < 4);
          found_match = true;
          for (int v = 0; v < 4; ++v) {
            if (matchv[(v+off)%4] &&
                matchv[(v+off)%4] != ref_conn[(v+ref)%4]) {
              match_forward = false;
              break;
            }
          }
          if (!match_forward) {
            for (int v = 0; v < 4; ++v) {
              if (matchv[(v+off)%4] &&
                  matchv[(v+off)%4] != ref_conn[(4-v+ref)%4]) {
                found_match = false;
                break;
              }
            }
          }
          
          if (!found_match)
            continue;
          
          int step = match_forward ? 1 : -1;
          for (int v = 0; v < 4; ++v) {
            if (!matchv[(v+off)%4]) {
              next_line.push_back( conn[(v+off)%4] );
              matchv[(v+off)%4] = ref_conn[(4+v*step+ref)%4];
            }
          }
          
          Mesh::ElementHandle h = ref_elems[k];
          ref_elems.clear();
          ref_elems.push_back(h);
        }
        if (!found_match) {
          std::cerr << "Input and reference mesh incompatible." << std::endl;
          exit(1);
        }
        
        ref_matched.insert( ref_elems.front() );
        inp_matched.insert( elems[j] );
        mesh.tag_set_vertex_data( match_tag, conn.size(), &conn[0], &matchv[0], err );        
      }
      
      if (!all_matched)
        next_line.push_back( inp_line[i] );
    }
    
    inp_line.swap( next_line );
    next_line.clear();
  }
  
  mesh.tag_get_vertex_data( match_tag, verts.size(), &verts[0], &ref_verts[0], err );
  assert(!err);
  mesh.tag_delete( match_tag, err );
  ref_mesh.vertices_get_coordinates( &ref_verts[0], &ref_coords[0], ref_verts.size(), err );
  assert(!err);
  std::vector<Vector3D> tag_coords(ref_coords.size());
  std::copy( ref_coords.begin(), ref_coords.end(), tag_coords.begin() );
  
  
  TagHandle coord_tag = mesh.tag_create( vertex_tag_name, Mesh::DOUBLE, 3, 0, err );
  assert(!err);
  mesh.tag_set_vertex_data( coord_tag, verts.size(), &verts[0], &tag_coords[0], err );
  assert(!err);
  
  TagHandle debug_tag = mesh.tag_create( "debug_vect", Mesh::DOUBLE, 3, 0, err );
  assert(!err);
  for (size_t i = 0; i < verts.size(); ++i)
    tag_coords[i] -= coords[i];
  mesh.tag_set_vertex_data( debug_tag, verts.size(), &verts[0], &tag_coords[0], err );
  assert(!err);
}
