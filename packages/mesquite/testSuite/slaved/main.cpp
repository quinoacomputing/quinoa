/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file main.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */


#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_TriLagrangeShape.hpp"
#include "Mesquite_QuadLagrangeShape.hpp"
#include "Mesquite_IdealShapeTarget.hpp"
#include "Mesquite_TShapeNB1.hpp"
#include "Mesquite_TQualityMetric.hpp"
#include "Mesquite_SteepestDescent.hpp"
#include "Mesquite_SlaveBoundaryVertices.hpp"
#include "Mesquite_PMeanPTemplate.hpp"
#include "Mesquite_InstructionQueue.hpp"
#include "Mesquite_MsqVertex.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_IdealShapeTarget.hpp"
#include "Mesquite_QualityAssessor.hpp"
#include "Mesquite_PatchData.hpp"

#include <iostream>
#include <algorithm>
#include <assert.h>

using namespace Mesquite;

const char DEFAULT_INPUT_FILE[] = SRCDIR "input.vtk";

void usage( const char* argv0 )
{
  std::cout << "Usage: " << argv0 << " [<input_file>] [-o <output_file_base>]" << std::endl;
  exit(1);
}

int check_slaved_coords( Mesh& mesh, Settings::HigherOrderSlaveMode mode,
                         std::string name, bool have_slaved_flag,
                         MsqError& err );

// compare vertex coodinates between two topologically equivalent meshes
int compare_node_coords( Mesh& mesh1, Mesh& mesh2, MsqError& err );

// verify that no element corner node is marked as slaved
int check_no_slaved_corners( Mesh& mesh, MsqError& err );

// compare slaved flags in mesh to slaved state in global patch data
int check_global_patch_slaved( Mesh& mesh, MsqError& err );

// tag vertices marked as slaved in the PatchData
void tag_patch_slaved( Mesh& mesh, 
                       Settings::HigherOrderSlaveMode mode,
                       MsqError& err );

int main( int argc, char* argv[] )
{
  const char* input_file = DEFAULT_INPUT_FILE;
  const char* output_file_base = 0;
  bool expect_output_base = false;
  for (int i = 1; i < argc; ++i) {
    if (expect_output_base) {
      output_file_base = argv[i];
      expect_output_base = false;
    }
    else if (!strcmp(argv[i],"-o"))
      expect_output_base = true;
    else if (input_file != DEFAULT_INPUT_FILE)
      usage(argv[0]);
    else
      input_file = argv[i];
  }
  if (expect_output_base)
    usage(argv[0]);
  
  MsqPrintError err(std::cerr);
  SlaveBoundaryVertices slaver(1);
  TShapeNB1 tmetric;
  IdealShapeTarget target;
  TQualityMetric metric( &target, &tmetric );
  PMeanPTemplate of( 1.0, &metric );
  SteepestDescent improver( &of );
  TerminationCriterion inner;
  inner.add_absolute_vertex_movement( 1e-3 );
  improver.set_inner_termination_criterion( &inner );
  QualityAssessor assess( &metric );
  InstructionQueue q;
  q.set_master_quality_improver( &improver, err );
  q.add_quality_assessor( &assess, err );
  
  TriLagrangeShape trishape;
  QuadLagrangeShape quadshape;
  q.set_mapping_function( &trishape );
  q.set_mapping_function( &quadshape );
  
  const int NUM_MODES = 4;
  
  Settings::HigherOrderSlaveMode modes[NUM_MODES] = 
    { Settings::SLAVE_NONE, 
      Settings::SLAVE_ALL,
      Settings::SLAVE_CALCULATED,
      Settings::SLAVE_FLAG
    };
  
  std::string names[NUM_MODES] = 
    { "NONE",
      "ALL",
      "CALCULATED",
      "FLAG" };
  
  MeshImpl meshes[NUM_MODES];
  std::vector<MeshDomainAssoc> meshes_and_domains;

  bool have_slaved_flag = true;
  std::vector<bool> flag(1);
  for (int i = 0; i < NUM_MODES; ++i) {
    std::cout << std::endl
              << "-----------------------------------------------" << std::endl
              << "     Mode: " << names[i] << std::endl
              << "-----------------------------------------------" << std::endl;
    
    meshes[i].read_vtk( input_file, err );
    if (err) return 1;
  
    if (modes[i] == Settings::SLAVE_CALCULATED) {
        q.add_vertex_slaver( &slaver, err );
    }
    else if (modes[i] == Settings::SLAVE_FLAG) {
      std::vector<Mesh::VertexHandle> verts;
      meshes[i].get_all_vertices( verts, err );
      if (err) return 1;
      meshes[i].vertices_get_slaved_flag( arrptr(verts), flag, 1, err );
      if (err) {
        have_slaved_flag = false;
        std::cout << "Skipped because input file does not contain slaved attribute" << std::endl;
        err.clear();
        continue;
      }
    }

  
    if (have_slaved_flag && modes[i] == Settings::SLAVE_FLAG) {
      if (!check_no_slaved_corners( meshes[i], err ) || err)
        return 1;
      if (!check_global_patch_slaved( meshes[i], err ) || err)
        return 1;
    }

  
    PlanarDomain plane;
    plane.fit_vertices( &meshes[i], err );
    if (err) return 1;

    q.set_slaved_ho_node_mode( modes[i] );
    meshes_and_domains.push_back(MeshDomainAssoc(&meshes[i], &plane));
    q.run_instructions( &meshes_and_domains[i], err );
    if (err) return 1;

    if (modes[i] == Settings::SLAVE_CALCULATED) {
        q.remove_vertex_slaver( &slaver, err );
    }
    
    if (output_file_base) {
      tag_patch_slaved( meshes[i], modes[i], err );
      std::string name(output_file_base);
      name += "-";
      name += names[i];
      name += ".vtk";
      meshes[i].write_vtk( name.c_str(), err );
      if (err) return 1;
    }
  }
  
  int exit_code = 0;
  if (input_file == DEFAULT_INPUT_FILE) {
    for (int i = 0; i < NUM_MODES; ++i) {
      std::cout << std::endl
                << "-----------------------------------------------" << std::endl
                << "     Mode: " << names[i] << std::endl
                << "-----------------------------------------------" << std::endl;
    
      exit_code += check_slaved_coords( meshes[i], modes[i], names[i], have_slaved_flag, err );
      if (err) return 1;
    }
    
      // flags should correspond to same slaved nodes as calculated,
      // so resulting meshes should be identical.
    if (have_slaved_flag) {
      int flag_idx = std::find( modes, modes+NUM_MODES, Settings::SLAVE_FLAG ) - modes;
      int calc_idx = std::find( modes, modes+NUM_MODES, Settings::SLAVE_CALCULATED ) - modes;
      exit_code += compare_node_coords( meshes[flag_idx], meshes[calc_idx], err );
      if (err) return 1;    
    }
  }
  
  return exit_code;
}

Vector3D get_slaved_coords( Mesh& mesh, Mesh::VertexHandle vertex, MsqError& err )
{
  std::vector<Mesh::ElementHandle> elem;
  std::vector<size_t> off;
  mesh.vertices_get_attached_elements( &vertex, 1, elem, off, err );
  if (MSQ_CHKERR(err)) return Vector3D(0.0);
  
  std::vector<Mesh::VertexHandle> verts;
  mesh.elements_get_attached_vertices( arrptr(elem), 1, verts, off, err );
  if (MSQ_CHKERR(err)) return Vector3D(0.0);
  
  EntityTopology type;
  mesh.elements_get_topologies( arrptr(elem), &type, 1, err );
  if (MSQ_CHKERR(err)) return Vector3D(0.0);
  
  size_t idx = std::find( verts.begin(), verts.end(), vertex ) - verts.begin();
  unsigned side_dim, side_num;
  TopologyInfo::side_from_higher_order( type, verts.size(), idx, side_dim, side_num, err);
  if (MSQ_CHKERR(err)) return Vector3D(0.0);
  
    // just return the mean of the corner vertices defining the side.
    // this should always be correct for mid-edge nodes but is a bit
    // dubious for mid-face or mid-region nodes.  But our default input
    // file (the only file for which we should be in this function) will
    // have only mid-edge nodes.
  unsigned n;
  const unsigned* side_vtx = TopologyInfo::side_vertices( type, side_dim, side_num, n, err );
  if (MSQ_CHKERR(err)) 
    return Vector3D(0.0);

  Vector3D sum(0.0);
  for (unsigned i = 0; i < n; ++i) {
    MsqVertex coords;
    Mesh::VertexHandle vtx = verts[side_vtx[i]];
    mesh.vertices_get_coordinates( &vtx, &coords, 1, err );
    if (MSQ_CHKERR(err)) return Vector3D(0.0);
    sum += coords;
  }
  return sum / n;
}

int check_slaved_coords( Mesh& mesh, 
                         Settings::HigherOrderSlaveMode mode,
                         std::string name, 
                         bool have_slaved_flag,
                         MsqError& err )
{
  const double EPSILON = 1e-4;

  // We can distinguish between nodes near the boundary and those
  // further inside based on the slaved attribute flag in the
  // the default input file.  The default input file should always
  // have it defined
  if (!have_slaved_flag) {
    std::cerr << "slaved flag not specified in input file. Cannot validate results." << std::endl;
    return 1;
  }
  
  // Get list of all higher-order nodes
  std::vector<Mesh::ElementHandle> elems;
  std::vector<Mesh::VertexHandle> verts;
  std::vector<size_t> offsets;
  mesh.get_all_elements( elems, err ); MSQ_ERRZERO(err);
  std::vector<EntityTopology> types(elems.size());
  mesh.elements_get_topologies( arrptr(elems), arrptr(types), elems.size(), err );
  MSQ_ERRZERO(err);
  mesh.elements_get_attached_vertices( arrptr(elems), elems.size(),
                                       verts, offsets, err );
  MSQ_ERRZERO(err);
  std::vector<Mesh::VertexHandle>::iterator r, e, w = verts.begin();
  for (size_t i = 0; i < elems.size(); ++i) {
    r = verts.begin() + offsets[i] + TopologyInfo::corners( types[i] );
    e = verts.begin() + offsets[i+1];
    w = std::copy( r, e, w );
  }
  std::sort( verts.begin(), w );
  verts.erase( std::unique( verts.begin(), w ), verts.end() );
  
  // Get lists of slaved and non-slaved free vertices, where 'slaved'
  // are those vertices far from the boundary that would be slaved if
  // mode == SLAVE_FLAG, not those that were actually slaved during
  // the optimization
  std::vector<bool> fixed, slaved;
  mesh.vertices_get_fixed_flag( arrptr(verts), fixed, verts.size(), err );
  if (MSQ_CHKERR(err)) {
    return 1;
  }
  mesh.vertices_get_slaved_flag( arrptr(verts), slaved, verts.size(), err );
  if (MSQ_CHKERR(err)) {
    return 1;
  }
  std::vector<Mesh::VertexHandle> free, slave;
  for (size_t i = 0; i < verts.size(); ++i) {
    if (!fixed[i] && slaved[i])
      slave.push_back( verts[i] );
    else if (!fixed[i] && !slaved[i])
      free.push_back( verts[i] );
  }
  
    // get all coordinates
  std::vector<MsqVertex> free_coords(free.size()), slave_coords(slave.size());
  mesh.vertices_get_coordinates( arrptr(free), arrptr(free_coords), free.size(), err );
  MSQ_ERRZERO(err);
  mesh.vertices_get_coordinates( arrptr(slave), arrptr(slave_coords), slave.size(), err );
  MSQ_ERRZERO(err);
  
  int error_count = 0;
  
    // Expect vertices near boundary to be at slave vertex positions
    // if mode was SLAVE_ALL
  if (mode == Settings::SLAVE_ALL) {
    for (size_t i = 0; i < free.size(); ++i) {
      Vector3D exp = get_slaved_coords( mesh, free[i], err );
      MSQ_ERRZERO(err);
      exp -= free_coords[i];
      if (exp.length() > EPSILON) {
        std::cerr << "Slaved vertex " << (size_t)free[i] << " not at expected slaved location" << std::endl;
        ++error_count;
      }
    }
  }
    // Otherwise, given the default input mesh, at least some of the 
    // vertices should be somewhere other than the slaved position
  else {
    int not_at_slaved_count = 0;
  
    for (size_t i = 0; i < free.size(); ++i) {
      Vector3D exp = get_slaved_coords( mesh, free[i], err );
      MSQ_ERRZERO(err);
      exp -= free_coords[i];
      if (exp.length() > EPSILON) {
        ++not_at_slaved_count;
      }
    }
    
    if (0 == not_at_slaved_count) {
      std::cerr << "All non-slaved vertices at slaved vertex locations" << std::endl;
      error_count += free.size();
    }
  }
  
    // expect all interior vertices to be at roughly the slaved location
  if (mode != Settings::SLAVE_NONE) {
    for (size_t i = 0; i < slave.size(); ++i) {
      Vector3D exp = get_slaved_coords( mesh, slave[i], err );
      MSQ_ERRZERO(err);
      exp -= slave_coords[i];
      if (exp.length() > EPSILON) {
        std::cerr << "Interior vertex " << (size_t)slave[i] << " not at expected slaved location" << std::endl;
        ++error_count;
      }
    }
  }
  
  return error_count;
}

int compare_node_coords( Mesh& mesh1, Mesh& mesh2, MsqError& err )
{
  const double EPSILON = 1e-4;

  std::vector<Mesh::VertexHandle> verts1, verts2;
  mesh1.get_all_vertices( verts1, err ); MSQ_ERRZERO(err);
  mesh2.get_all_vertices( verts2, err ); MSQ_ERRZERO(err);
  std::vector<MsqVertex> coords1(verts1.size()), coords2(verts2.size());
  mesh1.vertices_get_coordinates( arrptr(verts1), arrptr(coords1), verts1.size(), err );
  MSQ_ERRZERO(err);
  mesh2.vertices_get_coordinates( arrptr(verts2), arrptr(coords2), verts2.size(), err );
  MSQ_ERRZERO(err);
  
  int error_count = 0;
  assert(verts1.size() == verts2.size());
  for (size_t i = 0; i < verts1.size(); ++i) {
    assert( verts1[i] == verts2[i] );
    Vector3D diff = coords1[i] - coords2[i];
    if (diff.length() > EPSILON) {
      std::cerr << "Vertex coordinates differ between calculated and flagged "
                   "meshes for vertex " << (size_t)verts1[i] << std::endl;
      ++error_count;
    }
  }
  
  return error_count;
}

int check_no_slaved_corners( Mesh& mesh, MsqError& err )
{
  std::vector<Mesh::ElementHandle> elems;
  std::vector<Mesh::VertexHandle> verts;
  std::vector<EntityTopology> types;
  std::vector<size_t> offsets;
  mesh.get_all_elements( elems, err );  MSQ_ERRZERO(err);
  types.resize( elems.size() );
  mesh.elements_get_topologies( arrptr(elems), arrptr(types), elems.size(), err );  MSQ_ERRZERO(err);
  mesh.elements_get_attached_vertices( arrptr(elems), elems.size(), verts, offsets, err );
  MSQ_ERRZERO(err);
  
  std::vector<bool> slaved;
  mesh.vertices_get_slaved_flag( arrptr(verts), slaved, verts.size(), err );
  MSQ_ERRZERO(err);
  
  int error_count = 0;
  for (size_t i = 0; i < elems.size(); ++i) {
    unsigned n = TopologyInfo::corners( types[i] );
    for (unsigned j = 0; j < n; ++j) {
      if (slaved[offsets[i]+j]) {
        std::cerr << "Element " << (size_t)elems[i] << " corner " << j << " is slaved" <<std::endl;
        ++error_count;
      }
    }
  }
  
  return error_count == 0;
}

int check_global_patch_slaved( Mesh& mesh, MsqError& err )
{
  Settings s;
  s.set_slaved_ho_node_mode( Settings::SLAVE_FLAG );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, 0);
  Instruction::initialize_vertex_byte( &mesh_and_domain, &s, err ); MSQ_ERRZERO(err);
  
  PatchData pd;
  pd.attach_settings( &s );
  pd.set_mesh( &mesh );
  pd.fill_global_patch( err ); MSQ_ERRZERO(err);

  std::vector<bool> fixed, slaved;
  mesh.vertices_get_fixed_flag( pd.get_vertex_handles_array(), 
                                fixed, pd.num_nodes(), err );
  MSQ_ERRZERO(err);
  mesh.vertices_get_slaved_flag( pd.get_vertex_handles_array(), 
                                 slaved, pd.num_nodes(), err );
  MSQ_ERRZERO(err);
  
  const size_t first_free = 0;
  const size_t first_slaved = pd.num_free_vertices();
  const size_t first_fixed = pd.num_free_vertices() + pd.num_slave_vertices();
  int error_count = 0;
  for (size_t i = first_free; i < first_slaved; ++i) {
    if (fixed[i]) {
      std::cerr << "Vertex " << (size_t)pd.get_vertex_handles_array()[i]
                << " is fixed in mesh but free in PatchData" << std::endl;
      ++error_count;
    }
    if (slaved[i]) {
      std::cerr << "Vertex " << (size_t)pd.get_vertex_handles_array()[i]
                << " is slaved in mesh but free in PatchData" << std::endl;
      ++error_count;
    }
  }
  for (size_t i = first_slaved; i < first_fixed; ++i) {
    if (fixed[i]) {
      std::cerr << "Vertex " << (size_t)pd.get_vertex_handles_array()[i]
                << " is fixed in mesh but slaved in PatchData" << std::endl;
      ++error_count;
    }
    else if (!slaved[i]) {
      std::cerr << "Vertex " << (size_t)pd.get_vertex_handles_array()[i]
                << " is free in Mesh but slaved in PatchData" << std::endl;
      ++error_count;
    }
  }
  for (size_t i = first_fixed; i < pd.num_nodes(); ++i) {
    if (!fixed[i]) {
      std::cerr << "Vertex " << (size_t)pd.get_vertex_handles_array()[i]
                << " is not fixed in mesh but is in PatchData" << std::endl;
      ++error_count;
    }
  }
  return 0 == error_count;
}

void tag_patch_slaved( Mesh& mesh, 
                       Settings::HigherOrderSlaveMode mode,
                       MsqError& err )
{
  int zero = 0;
  TagHandle tag = mesh.tag_create( "pd_slaved", Mesh::INT, 1, &zero, err );
  MSQ_ERRRTN(err);
  
  Settings s;
  s.set_slaved_ho_node_mode( mode );
  PatchData pd;
  pd.attach_settings( &s );
  pd.set_mesh( &mesh );
  pd.fill_global_patch( err ); 
  MSQ_ERRRTN(err);

  const Mesh::VertexHandle* verts = pd.get_vertex_handles_array() + pd.num_free_vertices();
  std::vector<int> ones( pd.num_slave_vertices(), 1 );
  mesh.tag_set_vertex_data( tag, pd.num_slave_vertices(), verts, arrptr(ones), err );
  MSQ_ERRRTN(err);
}
