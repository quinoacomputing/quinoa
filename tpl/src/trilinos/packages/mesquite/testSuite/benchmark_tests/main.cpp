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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 13-08-12 by Boyd Tidwell
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

  Benchmark timing tests for various Mesquite wrappers.

 */
// DESCRIP-END.
//

#include "meshfiles.h"

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

#include "Mesquite.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_Vector3D.hpp"
#include "Mesquite_InstructionQueue.hpp"
#include "Mesquite_PatchData.hpp"
#include "Mesquite_TerminationCriterion.hpp"
#include "Mesquite_QualityAssessor.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_MsqTimer.hpp"

// algorythms
#include "Mesquite_LInfTemplate.hpp"
#include "Mesquite_EdgeLengthQualityMetric.hpp"
#include "Mesquite_LaplaceWrapper.hpp"
#include "Mesquite_UntangleWrapper.hpp"

#include "Mesquite_ShapeImprover.hpp"
#include "Mesquite_MsqTimer.hpp"
#include "Mesquite_SizeAdaptShapeWrapper.hpp"
#include "Mesquite_SphericalDomain.hpp"
#include "Mesquite_PaverMinEdgeLengthWrapper.hpp"
#include "Mesquite_DeformingDomainWrapper.hpp"
#include "Mesquite_MeshDomain1D.hpp"

using namespace Mesquite;

const char shape_improv_file_name_1[] = MESH_FILES_DIR "3D/vtk/hexes/untangled/1000hex-block-internal-bias.vtk";
const char shape_improv_file_name_2[] = MESH_FILES_DIR "3D/vtk/tets/untangled/tire.vtk";
const char laplacian_file_name_1[] = MESH_FILES_DIR "2D/vtk/quads/untangled/square_quad_10_rand.vtk";
const char laplacian_file_name_2[] = MESH_FILES_DIR "2D/vtk//quads/untangled/shashkov_quad.vtk";
const char untangle_file_name_1[] = MESH_FILES_DIR "2D/vtk/quads/tangled/tangled_horse1.vtk";
const char untangle_file_name_2[] = MESH_FILES_DIR "2D/vtk/quads/untangled//shest_grid32.vtk";
const char size_adapt_shape_file_name_1[] = MESH_FILES_DIR "2D/vtk/quads/untangled/bias-sphere-quads.vtk";
const char min_edge_length_file_name_1[] = MESH_FILES_DIR "2D/vtk/quads/untangled/quads_4by2_bad.vtk";
const char min_edge_length_file_name_2[] = MESH_FILES_DIR "2D/vtk/quads/untangled/shashkov_quad.vtk";
const char deforming_domain_file_name_1[] = MESH_FILES_DIR "3D/vtk/hexes/untangled/sph-10-zsquare.vtk";

void help(const char* argv0)
{
  std::cerr << "Parameters are not supported for this test" << std::endl;
            exit(1);
}
typedef std::vector<Mesh::ElementHandle> elem_vec_t;

// Assume mesh is on a sphere of radius 10 with an axis equal to Z.
// Return elements near poles and elements near equator.
void find_z10_extreme_elements( Mesh& mesh,
                                elem_vec_t& polar_elems,
                                elem_vec_t& equatorial_elems,
                                MsqError& err );

// Get areas of quads and tris
void elem_areas( Mesh& mesh, const elem_vec_t& elems, 
                 double& min, double& mean, double& max,
                 MsqError& err );

void classify_boundary( Mesh* mesh,
                        Mesh::VertexHandle corners_out[4],
                        std::vector<Mesh::VertexHandle> curves_out[4],
                        MsqError& err );
const double Z = 7.0;
// size of new domain
const double HD = 4.5;

int main(int argc, char* argv[])
{

  std::cout << std::endl << "********* Wrappers Timing Tests **********" 
            << std::endl << "Version "  << version_string(true) 
            << std::endl << std::endl;

  Mesquite::MsqPrintError err(cout);
  Mesquite::MeshImpl mesh;

// #################### Begin ShapeImprover tests ###################

  ShapeImprover si_wrapper;
  mesh.read_vtk(shape_improv_file_name_1, err);

  Timer t;  
  si_wrapper.run_instructions(&mesh, err); 
  if (err) return 1;
  double si_s_secs = t.since_birth();
  std::cout << std::endl << "ShapeImprover small file optimization completed in " 
            << si_s_secs << " seconds" << std::endl;

  mesh.clear();
  mesh.read_vtk(shape_improv_file_name_2, err); 
  
  t.reset();
  si_wrapper.run_instructions(&mesh, err); 
  if (err) return 1;
  double si_l_secs = t.since_birth();
  std::cout << std::endl << "ShapeImprover large file optimization completed in " 
            << si_l_secs << " seconds" << std::endl;

// #################### Begin LaplacianWrapper tests ###################
  
  Vector3D pnt1(0,0,5);
  Vector3D s_norm(0,0,1);
  Mesquite::PlanarDomain msq_geom(s_norm, pnt1);
  
  LaplaceWrapper lp_wrapper;

  mesh.clear();
  mesh.read_vtk(laplacian_file_name_1, err);
  if (err) return 1;

  MeshDomainAssoc mesh_and_domain4 = MeshDomainAssoc(&mesh, &msq_geom);
  t.reset();
  lp_wrapper.run_instructions(&mesh_and_domain4, err); 
  if (err) return 1;
  double lp_s_secs = t.since_birth();
  std::cout << std::endl << "LaplacianWrapper small file optimization completed in " 
            << lp_s_secs << " seconds" << std::endl;

  Vector3D pnt2(0,0,0);
  Mesquite::PlanarDomain msq_geom2(s_norm, pnt2);

  mesh.clear();
  mesh.read_vtk(laplacian_file_name_2, err);
  if (err) return 1;  

  MeshDomainAssoc mesh_and_domain5 = MeshDomainAssoc(&mesh, &msq_geom2);
  t.reset();
  lp_wrapper.run_instructions(&mesh_and_domain5, err); 
  if (err) return 1;
  double lp_l1_secs = t.since_birth();
  std::cout << std::endl << "LaplacianWrapper large file (term crit=0.001) completed in " 
            << lp_l1_secs << " seconds" << std::endl;

  mesh.clear();
  mesh.read_vtk(laplacian_file_name_2, err);
  if (err) return 1;  

  lp_wrapper.set_vertex_movement_limit_factor(0.1);
  t.reset();
  lp_wrapper.run_instructions(&mesh_and_domain5, err); 
  if (err) return 1;
  double lp_l2_secs = t.since_birth();
  std::cout << std::endl << "LaplacianWrapper large file (term crit=0.1) completed in " 
            << lp_l2_secs << " seconds" << std::endl;


// #################### Begin UntangleWrapper::BETA tests ###################

  mesh.clear();
  mesh.read_vtk(untangle_file_name_1, err);
  if (err) return 1;

  std::vector<Mesh::VertexHandle> verts;
  mesh.get_all_vertices( verts, err );
  if (err || verts.empty()) return 1;
  MsqVertex coords;
  mesh.vertices_get_coordinates( arrptr(verts), &coords, 1, err );
  if (err) return 1;
  Vector3D norm(0,0,1);
  PlanarDomain u_domain( norm, coords );

  UntangleWrapper::UntangleMetric metric = UntangleWrapper::BETA;
  UntangleWrapper un_wrapper (metric);
  un_wrapper.set_vertex_movement_limit_factor( 0.005 );

  MeshDomainAssoc mesh_and_domain3 = MeshDomainAssoc(&mesh, &u_domain);
  t.reset();
  un_wrapper.run_instructions( &mesh_and_domain3, err );
  if (err) return 1;

  double unb_s_secs = t.since_birth();
  std::cout << std::endl << "UntangleWrapper::BETA small file optimization completed in " 
            << unb_s_secs << " seconds" << std::endl;

  mesh.clear();
  mesh.read_vtk(untangle_file_name_2, err);
  if (err) return 1;

   // get domain
  verts.clear();
  mesh.get_all_vertices( verts, err );
  if (err || verts.empty()) return 1;
  MsqVertex coords2;
  mesh.vertices_get_coordinates( arrptr(verts), &coords2, 1, err );
  if (err) return 1;

  PlanarDomain un_domain2( norm, coords2 );
 
    MeshDomainAssoc mesh_and_domain6 = MeshDomainAssoc(&mesh, &un_domain2);
  t.reset();
  un_wrapper.run_instructions( &mesh_and_domain6, err );
  if (err) return 1;

  double unb_l1_secs = t.since_birth();
  std::cout << std::endl << "UntangleWrapper::BETA large file (term crit=0.005) completed in " 
            << unb_l1_secs << " seconds" << std::endl;

  mesh.clear();
  mesh.read_vtk(untangle_file_name_2, err);
  if (err) return 1;

   // get domain
  verts.clear();
  mesh.get_all_vertices( verts, err );
  if (err || verts.empty()) return 1;
  MsqVertex coords3;
  mesh.vertices_get_coordinates( arrptr(verts), &coords3, 1, err );
  if (err) return 1;

  PlanarDomain un_domain3( norm, coords3 );
 
  un_wrapper.set_vertex_movement_limit_factor( 0.1 );
  MeshDomainAssoc mesh_and_domain7 = MeshDomainAssoc(&mesh, &un_domain3);
  t.reset();
  un_wrapper.run_instructions( &mesh_and_domain7, err );
  if (err) return 1;

  double unb_l2_secs = t.since_birth();
  std::cout << std::endl << "UntangleWrapper::BETA large file (term crit=0.1) completed in " 
            << unb_l2_secs << " seconds" << std::endl;


// #################### Begin UntangleWrapper::SIZE tests ###################

  mesh.clear();
  mesh.read_vtk(untangle_file_name_1, err);
  if (err) return 1;

  verts.clear();
  mesh.get_all_vertices( verts, err );
  if (err || verts.empty()) return 1;
  MsqVertex coords2a;
  mesh.vertices_get_coordinates( arrptr(verts), &coords2a, 1, err );
  if (err) return 1;
  PlanarDomain u_domain3( norm, coords2a );

  UntangleWrapper::UntangleMetric metric2 = UntangleWrapper::SIZE;
  UntangleWrapper un_wrapper2s(metric2);
  UntangleWrapper un_wrapper2l(metric2);

  MeshDomainAssoc mesh_and_domain8 = MeshDomainAssoc(&mesh, &u_domain3);
  t.reset();
  un_wrapper2s.run_instructions( &mesh_and_domain8, err );
  if (err) return 1;
  double uns_s_secs = t.since_birth();
  std::cout << std::endl << "UntangleWrapper::SIZE small file optimization completed in " 
            << uns_s_secs << " seconds" << std::endl;

  mesh.clear();
  mesh.read_vtk(untangle_file_name_2, err);
  if (err) return 1;

   // get domain
  verts.clear();
  mesh.get_all_vertices( verts, err );
  if (err || verts.empty()) return 1;
  MsqVertex coords4;
  mesh.vertices_get_coordinates( arrptr(verts), &coords4, 1, err );
  if (err) return 1;

  PlanarDomain un_domain4( norm, coords4 );
  un_wrapper2s.set_vertex_movement_limit_factor( 0.005 );

  MeshDomainAssoc mesh_and_domain9 = MeshDomainAssoc(&mesh, &un_domain4);
  t.reset();
  un_wrapper2s.run_instructions( &mesh_and_domain9, err );
  if (err) return 1;

  double uns_l1_secs = t.since_birth();
  std::cout << std::endl << "UntangleWrappe::SIZE large file (term crit=0.005) completed in " 
            << uns_l1_secs << " seconds" << std::endl;

  mesh.clear();
  mesh.read_vtk(untangle_file_name_2, err);
  if (err) return 1;

  mesh.get_all_vertices( verts, err );
  if (err || verts.empty()) return 1;

  mesh.vertices_get_coordinates( arrptr(verts), &coords4, 1, err );
  if (err) return 1;

  un_wrapper2l.set_vertex_movement_limit_factor( 0.1 );
  t.reset();
  un_wrapper2l.run_instructions( &mesh_and_domain9, err );
  if (err) return 1;

  double uns_l2_secs = t.since_birth();
  std::cout << std::endl << "UntangleWrappe::SIZE large file (term crit=0.1) completed in " 
            << uns_l2_secs << " seconds" << std::endl;


// #################### Begin UntangleWrapper::SHAPESIZE tests ###################

  mesh.clear();
  mesh.read_vtk(untangle_file_name_1, err);
  if (err) return 1;

  verts.clear();
  mesh.get_all_vertices( verts, err );
  if (err || verts.empty()) return 1;
  MsqVertex coords5;
  mesh.vertices_get_coordinates( arrptr(verts), &coords5, 1, err );
  if (err) return 1;
  PlanarDomain u_domain5( norm, coords3 );

  UntangleWrapper::UntangleMetric metric3 = UntangleWrapper::SHAPESIZE;
  UntangleWrapper un_wrapper3(metric3);

  MeshDomainAssoc mesh_and_domain10 = MeshDomainAssoc(&mesh, &u_domain5);
  t.reset();
  un_wrapper3.run_instructions( &mesh_and_domain10, err );
  if (err) return 1;

  double unss_s_secs = t.since_birth();
  std::cout << std::endl << "UntangleWrapper::SHAPESIZE small file optimization completed in " 
            << unss_s_secs << " seconds" << std::endl;

  mesh.clear();
  mesh.read_vtk(untangle_file_name_2, err);
  if (err) return 1;

   // get domain
  verts.clear();
  mesh.get_all_vertices( verts, err );
  if (err || verts.empty()) return 1;
  MsqVertex coords6;
  mesh.vertices_get_coordinates( arrptr(verts), &coords6, 1, err );
  if (err) return 1;

  PlanarDomain un_domain6( norm, coords6 );
 
  MeshDomainAssoc mesh_and_domain11 = MeshDomainAssoc(&mesh, &un_domain6);
  t.reset();
  un_wrapper3.run_instructions( &mesh_and_domain11, err );
  if (err) return 1;

  double unss_l_secs = t.since_birth();
  std::cout << std::endl << "UntangleWrapper::SHAPESIZE large file optimization completed in " 
            << unss_l_secs << " seconds" << std::endl;

  // #################### Begin SizeAdaptShapeWrapper tests ###################

  mesh.clear();
  mesh.read_vtk(size_adapt_shape_file_name_1, err);
  if (err) return 1;
 
  elem_vec_t polar, equatorial;
  find_z10_extreme_elements( mesh, polar, equatorial, err ); 
  if (err) return 1;
  
  double eq_min, eq_max, eq_mean, pol_min, pol_max, pol_mean;
  elem_areas( mesh, polar, pol_min, pol_mean, pol_max, err ); 
  if (err) return 1;
  elem_areas( mesh, equatorial, eq_min, eq_mean, eq_max, err ); 
  if (err) return 1;
  
  SphericalDomain geom( Vector3D(0,0,0), 10.0 );
  SizeAdaptShapeWrapper sas_wrapper1(1e-2, 50);
  SizeAdaptShapeWrapper sas_wrapper2(1e-1, 50);

  MeshDomainAssoc mesh_and_domain12 = MeshDomainAssoc(&mesh, &geom);
  t.reset();
  sas_wrapper1.run_instructions( &mesh_and_domain12, err);
  if (err) return 1;
  double sas1_secs = t.since_birth();
  std::cout << std::endl << "SizeAdaptShapeWrapper (term crit=0.01) completed in " 
            << sas1_secs << " seconds" << std::endl;

  mesh.clear();
  mesh.read_vtk(size_adapt_shape_file_name_1, err);
  if (err) return 1;

  t.reset();
  sas_wrapper2.run_instructions( &mesh_and_domain12, err);
  if (err) return 1;
  double sas2_secs = t.since_birth();
  std::cout << std::endl << "SizeAdaptShapeWrapper (term crit=0.1) completed in " 
            << sas2_secs << " seconds" << std::endl;

  // #################### Begin PaverMinEdgeLengthWrapper tests ###################

  PaverMinEdgeLengthWrapper mel_wrapper1(.005, 50);
  PaverMinEdgeLengthWrapper mel_wrapper2(.1, 50);

  mesh.clear();
  mesh.read_vtk(min_edge_length_file_name_1, err); 
  
  t.reset();
  mel_wrapper1.run_instructions(&mesh, err); 
  if (err) return 1;
  double mel_s_secs = t.since_birth();
  std::cout << std::endl << "PaverMinEdgeLengthWrapper small file optimization completed in " 
            << mel_s_secs << " seconds" << std::endl;

  mesh.clear();
  mesh.read_vtk(min_edge_length_file_name_2, err); 
  
  t.reset();
  mel_wrapper1.run_instructions(&mesh, err); 
  if (err) return 1;
  double mel1_l_secs = t.since_birth();
  std::cout << std::endl << "PaverMinEdgeLengthWrapper large file (term crit=0.005) completed in " 
            << mel1_l_secs << " seconds" << std::endl;


  mesh.clear();
  mesh.read_vtk(min_edge_length_file_name_2, err); 
  t.reset();
  mel_wrapper2.run_instructions(&mesh, err); 
  if (err) return 1;
  double mel2_l_secs = t.since_birth();
  std::cout << std::endl << "PaverMinEdgeLengthWrapper large file (term crit=0.1) completed in " 
            << mel2_l_secs << " seconds" << std::endl;

  // #################### Begin DeformingDomainWrapper tests ###################

  
    // load mesh
  mesh.clear();
  mesh.read_vtk( deforming_domain_file_name_1, err ); 
  if (MSQ_CHKERR(err)) return 1;

  std::vector<Mesh::VertexHandle> curves[4];
  Mesh::VertexHandle corners[4];
  classify_boundary( &mesh, corners, curves, err );
  if (MSQ_CHKERR(err)) return 1;
  
    // new, "deformed" domain will be an 2HDx2HD planar square
  const double corner_coords[][3] = { {-HD,-HD, Z},
                                      { HD,-HD, Z},
                                      { HD, HD, Z},
                                      {-HD, HD, Z} };
  LineDomain lines[4] = { 
    LineDomain( Vector3D(corner_coords[0]), Vector3D( 1, 0, 0) ),
    LineDomain( Vector3D(corner_coords[1]), Vector3D( 0, 1, 0) ),
    LineDomain( Vector3D(corner_coords[2]), Vector3D(-1, 0, 0) ),
    LineDomain( Vector3D(corner_coords[3]), Vector3D( 0,-1, 0) ) };
  PlanarDomain surface( PlanarDomain::XY, Z );
  
    // save initial mesh state
  DeformingCurveSmoother curve_tool;
  for (int i = 0; i < 4; ++i) {
    curve_tool.store_initial_mesh( &mesh, &curves[i][0], curves[i].size(), &lines[i], err );
    if (MSQ_CHKERR(err)) return 1;
  }
  DeformingDomainWrapper dd_wrapper;
  dd_wrapper.store_initial_mesh( &mesh, err );
  if (MSQ_CHKERR(err)) return 1;
  
    // move corner vertices to new location
  for (int i = 0; i < 4; ++i) {
    Vector3D vect(corner_coords[i]);
    mesh.vertex_set_coordinates( corners[i], vect, err );
    if (MSQ_CHKERR(err)) return 1;
  }
  std::vector<bool> fixed(4,true);
  mesh.vertices_set_fixed_flag( corners, fixed, 4, err );
  if (MSQ_CHKERR(err)) return 1;
  
    // smooth curves
  for (int i = 0; i < 4; ++i) {
    curve_tool.smooth_curve( &mesh, &curves[i][0], curves[i].size(), &lines[i],
                             DeformingCurveSmoother::PROPORTIONAL, err );
    if (MSQ_CHKERR(err)) return 1;
    fixed.resize(curves[i].size(),true);
    mesh.vertices_set_fixed_flag( &curves[i][0], fixed, curves[i].size(), err );
    if (MSQ_CHKERR(err)) return 1;
  }
  
  MeshDomainAssoc mesh_and_domain1 = MeshDomainAssoc(&mesh, &surface);
  t.reset();
  dd_wrapper.run_instructions( &mesh_and_domain1, err );
  if (MSQ_CHKERR(err)) return 1;
  double dd_secs = t.since_birth();
  std::cout << std::endl << "DeformingDomainWrapper file (term crit=0.01) completed in " 
            << dd_secs << " seconds" << std::endl;

    // Do it all again for the next test
  mesh.clear();
  mesh.read_vtk( deforming_domain_file_name_1, err ); 
  if (MSQ_CHKERR(err)) return 1;

  std::vector<Mesh::VertexHandle> curves2[4];
  Mesh::VertexHandle corners2[4];
  classify_boundary( &mesh, corners2, curves2, err );
  if (MSQ_CHKERR(err)) return 1;
  
    // new, "deformed" domain will be an 2HDx2HD planar square
  const double corner_coords2[][3] = { {-HD,-HD, Z},
                                      { HD,-HD, Z},
                                      { HD, HD, Z},
                                      {-HD, HD, Z} };
  LineDomain lines2[4] = { 
    LineDomain( Vector3D(corner_coords2[0]), Vector3D( 1, 0, 0) ),
    LineDomain( Vector3D(corner_coords2[1]), Vector3D( 0, 1, 0) ),
    LineDomain( Vector3D(corner_coords2[2]), Vector3D(-1, 0, 0) ),
    LineDomain( Vector3D(corner_coords2[3]), Vector3D( 0,-1, 0) ) };
  PlanarDomain surface2( PlanarDomain::XY, Z );
  
    // save initial mesh state
  DeformingCurveSmoother curve_tool2;
  for (int i = 0; i < 4; ++i) {
    curve_tool2.store_initial_mesh( &mesh, &curves2[i][0], curves2[i].size(), &lines2[i], err );
    if (MSQ_CHKERR(err)) return 1;
  }
  DeformingDomainWrapper dd_wrapper2;
  dd_wrapper2.store_initial_mesh( &mesh, err );
  if (MSQ_CHKERR(err)) return 1;
  
    // move corner vertices to new location
  for (int i = 0; i < 4; ++i) {
    Vector3D vect(corner_coords2[i]);
    mesh.vertex_set_coordinates( corners2[i], vect, err );
    if (MSQ_CHKERR(err)) return 1;
  }
  std::vector<bool> fixed2(4,true);
  mesh.vertices_set_fixed_flag( corners2, fixed2, 4, err );
  if (MSQ_CHKERR(err)) return 1;
  
    // smooth curves
  for (int i = 0; i < 4; ++i) {
    curve_tool2.smooth_curve( &mesh, &curves2[i][0], curves2[i].size(), &lines2[i],
                              DeformingCurveSmoother::PROPORTIONAL, err );
    if (MSQ_CHKERR(err)) return 1;
    fixed2.resize(curves2[i].size(),true);
    mesh.vertices_set_fixed_flag( &curves2[i][0], fixed2, curves2[i].size(), err );
    if (MSQ_CHKERR(err)) return 1;
  }
  
  dd_wrapper2.set_vertex_movement_limit_factor(0.1);
  MeshDomainAssoc mesh_and_domain2 = MeshDomainAssoc(&mesh, &surface2);
  t.reset();
  dd_wrapper2.run_instructions( &mesh_and_domain2, err );
  if (MSQ_CHKERR(err)) return 1;
  double dd_secs2 = t.since_birth();
  std::cout << std::endl << "DeformingDomainWrapper file (term crit=0.1) completed in " 
            << dd_secs2 << " seconds" << std::endl;

  // Timing Summary
  std::cout << std::endl << "********* Wrappers Timing Summary **********" 
            << std::endl << "Version "  << version_string(true) 
            << std::endl << std::endl;
  std::cout << "ShapeImprover small file optimization completed in " 
            << si_s_secs << " seconds" << std::endl;
  std::cout << "ShapeImprover large file optimization completed in " 
            << si_l_secs << " seconds" << std::endl;
  std::cout << "LaplacianWrapper small file optimization completed in " 
            << lp_s_secs << " seconds" << std::endl;
  std::cout << "LaplacianWrapper large file optimization (term crit=0.001) in " 
            << lp_l1_secs << " seconds" << std::endl;
  std::cout << "LaplacianWrapper large file optimization (term crit=0.1) in " 
            << lp_l2_secs << " seconds" << std::endl;
  std::cout << "UntangleWrapper::BETA small file optimization completed in " 
            << unb_s_secs << " seconds" << std::endl;
  std::cout << "UntangleWrapper::BETA large file (term crit=0.005) completed in " 
            << unb_l1_secs << " seconds" << std::endl;
  std::cout << "UntangleWrapper::BETA large file (term crit=0.1) completed in " 
            << unb_l2_secs << " seconds" << std::endl;
  std::cout << "UntangleWrapper::SIZE small file optimization completed in " 
            << uns_s_secs << " seconds" << std::endl;
  std::cout << "UntangleWrapper::SIZE large file (term crit=0.005) completed in " 
            << uns_l1_secs << " seconds" << std::endl;
  std::cout << "UntangleWrapper::SIZE large file (term crit=0.1) completed in " 
            << uns_l2_secs << " seconds" << std::endl;
  std::cout << "UntangleWrapper::SHAPESIZE small file optimization completed in " 
            << unss_s_secs << " seconds" << std::endl;
  std::cout << "UntangleWrapper::SHAPESIZE large file optimization completed in " 
            << unss_l_secs << " seconds" << std::endl;
  std::cout << "SizeAdaptShapeWrapper (term crit=0.01) completed in " 
            << sas1_secs << " seconds" << std::endl;
  std::cout << "SizeAdaptShapeWrapper (term crit=0.1) completed in " 
            << sas2_secs << " seconds" << std::endl;
  std::cout << "PaverMinEdgeLengthWrapper small file optimization completed in " 
            << mel_s_secs << " seconds" << std::endl;
  std::cout << "PaverMinEdgeLengthWrapper large file (term crit=0.005) completed in " 
            << mel1_l_secs << " seconds" << std::endl;
  std::cout << "PaverMinEdgeLengthWrapper large file (term crit=0.1) completed in " 
            << mel2_l_secs << " seconds" << std::endl;
  std::cout << "DeformingDomainWrapper file (term crit=0.01) completed in " 
            << dd_secs << " seconds" << std::endl;
  std::cout << "DeformingDomainWrapper file (term crit=0.1) completed in " 
            << dd_secs2 << " seconds" << std::endl;
  
  return 0;
}

void find_z10_extreme_elements( Mesh& mesh,
                                elem_vec_t& polar_elems,
                                elem_vec_t& equatorial_elems,
                                MsqError& err )
{
  elem_vec_t elems;
  mesh.get_all_elements( elems, err ); MSQ_ERRRTN(err);
  
  std::vector<Mesh::VertexHandle> verts;
  std::vector<MsqVertex> coords;
  std::vector<size_t> junk;
  for (elem_vec_t::iterator i = elems.begin(); i != elems.end(); ++i) {
    verts.clear(); junk.clear();
    mesh.elements_get_attached_vertices( &*i, 1, verts, junk, err ); MSQ_ERRRTN(err);
    coords.resize(verts.size());
    mesh.vertices_get_coordinates( arrptr(verts), arrptr(coords), verts.size(), err ); MSQ_ERRRTN(err);
    
    for (std::vector<MsqVertex>::iterator j = coords.begin(); j != coords.end(); ++j) {
      double z = (*j)[2];
      if (fabs(z) < 1e-6) {
        equatorial_elems.push_back(*i);
        break;
      }
      else if (fabs(z) - 10 < 1e-6) {
        polar_elems.push_back(*i);
        break;
      }
    }
  }
}

void elem_areas( Mesh& mesh, const elem_vec_t& elems, 
                 double& min, double& mean, double& max,
                 MsqError& err )
{
  min = HUGE_VAL;
  max = -1;
  mean = 0.0;
  
  std::vector<EntityTopology> types(elems.size());
  mesh.elements_get_topologies( arrptr(elems), arrptr(types), elems.size(), err ); MSQ_ERRRTN(err);
  
  std::vector<Mesh::VertexHandle> verts;
  std::vector<MsqVertex> coords;
  std::vector<size_t> junk;
  for (size_t i = 0; i < elems.size(); ++i) {
    verts.clear(); junk.clear();
    mesh.elements_get_attached_vertices( &elems[i], 1, verts, junk, err ); MSQ_ERRRTN(err);
    coords.resize(verts.size());
    mesh.vertices_get_coordinates( arrptr(verts), arrptr(coords), verts.size(), err ); MSQ_ERRRTN(err);
    
    Vector3D v1, v2;
    double area;
    if (types[i] == TRIANGLE) {
      assert(coords.size() == 3);
      v1 = coords[1] - coords[0];
      v2 = coords[2] - coords[0];
      area = 0.5 * (v1 * v2).length();
    }
    else if (types[i] == QUADRILATERAL) {
      assert(coords.size() == 4);
      v1 = coords[0] + coords[1] - coords[2] - coords[3];
      v2 = coords[0] + coords[3] - coords[1] - coords[2];
      area = 0.25 * (v1 * v2).length();
    }
    else {
      MSQ_SETERR(err)("Input file contains volume elements", MsqError::UNSUPPORTED_ELEMENT);
      return;
    }
    
    if (min > area)
      min = area;
    if (max < area)
      max = area;
    mean += area;
  }
  
  mean /= elems.size();
}

#define CHKMESH( COND, ERR ) do { \
  if (!(COND)) { \
    MSQ_SETERR(err)("Unexpected mesh topology", MsqError::INVALID_MESH); \
    return; \
  } } while (false)

void classify_boundary( Mesh* mesh,
                        Mesh::VertexHandle corners_out[4],
                        std::vector<Mesh::VertexHandle> curves_out[4],
                        MsqError& err )
{
  std::vector<Mesh::VertexHandle> verts;
  mesh->get_all_vertices( verts, err ); MSQ_ERRRTN(err);
  
    // Find the corner vertex that has negative X and Y coordinates
  Mesh::VertexHandle start;
  bool have_start = false;
  std::vector<Mesh::ElementHandle> elems;
  std::vector<size_t> offsets;
  for (size_t i = 0; i < verts.size(); ++i) {
    elems.clear();
    offsets.clear();
    mesh->vertices_get_attached_elements( &verts[i], 1, elems, offsets, err );
    MSQ_ERRRTN(err);
    if (elems.size() == 1) {
      MsqVertex coords;
      mesh->vertices_get_coordinates( &verts[i], &coords, 1, err ); MSQ_ERRRTN(err);
      if (coords[0] < 0.0 && coords[1] < 0.0) {
        CHKMESH( !have_start, err );
        have_start = true;
        start = verts[i];
      }
    }
  }
  CHKMESH( have_start, err );
  
    // starting at a a corner vertex, find skin vertices
  std::vector<Mesh::VertexHandle> boundary;
  boundary.push_back(start);
  elems.clear();
  offsets.clear();
  mesh->vertices_get_attached_elements( &start, 1, elems, offsets, err ); MSQ_ERRRTN(err);
  Mesh::ElementHandle prev = elems.front();
  corners_out[0] = start;
  int ncorner = 1;
  do {
    verts.clear();
    offsets.clear();
    mesh->elements_get_attached_vertices( &prev, 1, verts, offsets, err );
    size_t idx = std::find(verts.begin(), verts.end(), boundary.back()) - verts.begin();
    CHKMESH( idx < verts.size(), err );
    
    Mesh::VertexHandle next = verts[(idx+1)%verts.size()];
    elems.clear();
    offsets.clear();
    mesh->vertices_get_attached_elements( &next, 1, elems, offsets, err ); MSQ_ERRRTN(err);
    CHKMESH( elems.size() == 1 || elems.size() == 2, err );
    if (elems.size() == 2) {
      idx = std::find(elems.begin(), elems.end(), prev) - elems.begin();
      CHKMESH( idx < 2, err );
      prev = elems[1-idx];
    }

    if (elems.size() == 1) {
      CHKMESH( ncorner < 5, err );
      if (ncorner < 4) // don't add start vertx twice
        corners_out[ncorner] = next;
      ++ncorner;
    }
    
    boundary.push_back( next );
  } while (boundary.front() != boundary.back());
  CHKMESH( ncorner == 5, err );
  
    // find curve vertices
  size_t idx = 0;
  for (int c = 0; c < 4; ++c) {
    Mesh::VertexHandle s, e;
    s = corners_out[c];
    e = corners_out[(c+1)%4];
    CHKMESH(idx < boundary.size(), err );
    CHKMESH(boundary[idx] == s, err);
    for (; boundary[idx] != e; ++idx)
      curves_out[c].push_back( boundary[idx] );
    curves_out[c].push_back( e );
  }
}
