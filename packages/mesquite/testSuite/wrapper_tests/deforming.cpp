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


/** \file deform.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Mesquite_DeformingDomainWrapper.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_MeshDomain1D.hpp"
#include "Mesquite_PlanarDomain.hpp"
#include "Mesquite_QualityAssessor.hpp"
#include "Mesquite_MsqVertex.hpp"

#include <iostream>
#include <cstdlib>
#include <algorithm>

using namespace Mesquite;

const char INPUT_FILE[] = SRCDIR "sph-10-zsquare.vtk";
//const char INPUT_FILE[] = "test.vtk";
const double Z = 7.0;
// size of new domain
const double HD = 4.5;

void usage( const char* argv0 )
{
  std::cerr << "Usage: " << argv0 << " [-t <file>] [-d <file>] [-f <file>]" << std::endl
            << "\t-t : write initial mesh with tags containing reference mesh" << std::endl
            << "\t-d : write file containing deformed mesh with smoothed curves " << std::endl
            << "\t-f : write file containing final mesh" << std::endl
            << std::endl;
  std::exit(1);
}

void classify_boundary( Mesh* mesh,
                        Mesh::VertexHandle corners_out[4],
                        std::vector<Mesh::VertexHandle> curves_out[4],
                        MsqError& err );

void cond_write_file( MeshImpl& mesh, const char* filename ) 
{
  if (filename) {
    MsqPrintError err(std::cerr);
    mesh.write_vtk( filename, err );
    if (MSQ_CHKERR(err)) {
      std::cerr << filename << ": failed to write file" << std::endl;
      exit(1);
    }
    std::cout << "Wrote file: " << filename << std::endl;
  }
}

int main( int argc, char* argv[] )
{
  const char* deformed_file = 0;
  const char* smoothed_file = 0;
  const char* tag_file = 0;
  struct { const char* flag; const char** name_ptr; } flags[] = 
    { { "-d", &deformed_file },
      { "-t", &tag_file      },
      { "-f", &smoothed_file },
      { 0, 0 } };
  
  for (int i = 1; i < argc; ++i) {
    int j;
    for (j = 0; flags[j].flag && strcmp( flags[j].flag, argv[i] ); ++j);
    if (!flags[j].flag) {
      std::cerr << "Invalid argument: \"" << argv[i] << '"' << std::endl;
      usage(argv[0]);
    }
    else if (++i == argc) {
      std::cerr << "Expected argument following \"" << argv[i-1] << '"' << std::endl;
      usage(argv[0]);
    }
    *(flags[j].name_ptr) = argv[i];
  }
  
    // load mesh
  MsqPrintError err(std::cerr);
  MeshImpl mesh;
  mesh.read_vtk( INPUT_FILE, err ); 
  if (MSQ_CHKERR(err)) return 1;
  
    // find boundary vertices
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
  DeformingDomainWrapper wrapper;
  wrapper.store_initial_mesh( &mesh, err );
  if (MSQ_CHKERR(err)) return 1;
  
  cond_write_file( mesh, tag_file );
  
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
  
  cond_write_file( mesh, deformed_file );
  
    // smooth surface mesh
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &surface);
  wrapper.run_instructions( &mesh_and_domain, err );
  if (MSQ_CHKERR(err)) return 1;
  
  cond_write_file( mesh, smoothed_file );
  return wrapper.quality_assessor().invalid_elements();
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
                             
