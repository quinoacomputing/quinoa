/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file main.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */



#include "meshfiles.h"

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <algorithm>

#include "Mesquite.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_SizeAdaptShapeWrapper.hpp"
#include "Mesquite_SphericalDomain.hpp"
#include "Mesquite_MsqVertex.hpp"

using namespace Mesquite;

const char DEFAULT_INPUT[] = MESH_FILES_DIR "2D/vtk/quads/untangled/bias-sphere-quads.vtk";

void help(const char* argv0)
{
  std::cerr << "Usage: " << argv0 << " [<output_file>]" << std::endl
            << "  Input file is always: " << DEFAULT_INPUT << std::endl
            << "  defualt is no output file" << std::endl;
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

int main(int argc, char* argv[])
{
  const char* input_file = DEFAULT_INPUT;
  const char* output_file = NULL;
  switch (argc) {
    default:
      help(argv[0]);
    case 2:
      if (!strcmp(argv[1],"-h"))
        help(argv[0]);
      output_file = argv[1];
    case 1:
      ;
  }

    /* Read a VTK Mesh file */
  MsqPrintError err(cout);
  Mesquite::MeshImpl mesh;
  mesh.read_vtk( input_file, err);
  if (err) {
    std::cerr << input_file << ": file read failed" << endl;
    return 1;
  }
  
  elem_vec_t polar, equatorial;
  find_z10_extreme_elements( mesh, polar, equatorial, err ); 
  if (err) return 1;
  
  double eq_min, eq_max, eq_mean, pol_min, pol_max, pol_mean;
  elem_areas( mesh, polar, pol_min, pol_mean, pol_max, err ); 
  if (err) return 1;
  elem_areas( mesh, equatorial, eq_min, eq_mean, eq_max, err ); 
  if (err) return 1;
  
    /* Run optimizer */
  SphericalDomain geom( Vector3D(0,0,0), 10.0 );
  SizeAdaptShapeWrapper smoother(1e-2);
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &geom);
  smoother.run_instructions( &mesh_and_domain, err);
  if (err) return 1;
  
  if (output_file) {
    mesh.write_vtk( output_file, err );
    if (err) {
      std::cerr << output_file << ": file write failed" << endl;
      return 1;
    }
    else {
      std::cout << "Wrote file: " << output_file << endl;
    }
  }
  
  double eq2_min, eq2_max, eq2_mean, pol2_min, pol2_max, pol2_mean;
  elem_areas( mesh, polar, pol2_min, pol2_mean, pol2_max, err ); 
  if (err) return 1;
  elem_areas( mesh, equatorial, eq2_min, eq2_mean, eq2_max, err ); 
  if (err) return 1;
  
  double eq_min_pct = 100*fabs(eq_min - eq2_min)/eq_min;
  double eq_max_pct = 100*fabs(eq_max - eq2_max)/eq_max;
  double eq_mean_pct = 100*fabs(eq_mean - eq2_mean)/eq_mean;
  double pol_min_pct = 100*fabs(pol_min - pol2_min)/pol_min;
  double pol_max_pct = 100*fabs(pol_max - pol2_max)/pol_max;
  double pol_mean_pct = 100*fabs(pol_mean - pol2_mean)/pol_mean;
  
  std::printf("\n\n");
  std::printf( "AREAS:      Initial     Final       Difference  Change\n" );
  std::printf( "Polar:\n");
  std::printf( "  Minimum:  %10f  %10f  %10f  %10f%%\n", pol_min, pol2_min, fabs(pol_min - pol2_min), pol_min_pct );
  std::printf( "  Average:  %10f  %10f  %10f  %10f%%\n", pol_mean, pol2_mean, fabs(pol_mean - pol2_mean), pol_mean_pct );
  std::printf( "  Maximum:  %10f  %10f  %10f  %10f%%\n", pol_max, pol2_max, fabs(pol_max - pol2_max), pol_max_pct );
  std::printf( "Equatorial:\n");
  std::printf( "  Minimum:  %10f  %10f  %10f  %10f%%\n", eq_min, eq2_min, fabs(eq_min - eq2_min), eq_min_pct );
  std::printf( "  Average:  %10f  %10f  %10f  %10f%%\n", eq_mean, eq2_mean, fabs(eq_mean - eq2_mean), eq_mean_pct );
  std::printf( "  Maximum:  %10f  %10f  %10f  %10f%%\n", eq_max, eq2_max, fabs(eq_max - eq2_max), eq_max_pct );
  
  bool success = pol_min_pct < 6.0 && pol_max_pct < 6.5 && eq_min_pct < 25.0 && eq_max_pct < 6.5;
  return !success;
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
