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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 23-Jul-03 at 18:04:37 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#include "meshfiles.h"

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>
#include <vector>
#include <algorithm>

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "MeshImpl.hpp"
#include "ViscousCFDTetShapeWrapper.hpp"
#include "MsqVertex.hpp"

using namespace Mesquite;

const char DEFAULT_INPUT[] = MESH_FILES_DIR "3D/vtk/tets/untangled/flat-tet-sphere.vtk";

void help(const char* argv0)
{
  std::cerr << "Usage: " << argv0 << " [<input_file>] [<output_file>]" << std::endl
            << "  default input file is: " << DEFAULT_INPUT << std::endl
            << "  defualt is no output file" << std::endl
            << "  Warning: input mesh is assumed to lie in Z=5 plane" << std::endl;
  exit(1);
}

/* For each tet, calculate ratio of min/max dihedral angle. 
 * Return stats for all tets */
void tet_dihedral_angle_ratios( Mesh& mesh,
                                double& ratio_min,
                                double& ratio_avg,
                                double& ratio_max,
                                MsqError& err );

int main(int argc, char* argv[])
{
  const char* input_file = DEFAULT_INPUT;
  const char* output_file = NULL;
  switch (argc) {
    default:
      help(argv[0]);
    case 3:
      if (!strcmp(argv[2],"-h"))
        help(argv[0]);
      output_file = argv[2];
    case 2:
      if (!strcmp(argv[1],"-h"))
        help(argv[0]);
      input_file = argv[1];
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
  
  double min_ini, avg_ini, max_ini;
  tet_dihedral_angle_ratios( mesh, min_ini, avg_ini, max_ini, err );
  
    /* Run optimizer */
  ViscousCFDTetShapeWrapper smoother(1e-2);
  smoother.run_instructions( &mesh, err);
  if (err) return 1;
  
  double min_fin, avg_fin, max_fin;
  tet_dihedral_angle_ratios( mesh, min_fin, avg_fin, max_fin, err );
  
  cout << "\tMinimum \tAverage \tMaximum "<<endl
       << "Init\t"<<min_ini<<'\t'<<avg_ini<<'\t'<<max_ini<<endl
       << "Fini\t"<<min_fin<<'\t'<<avg_fin<<'\t'<<max_fin<<endl;
  
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
      
  if (min_ini > min_fin || avg_ini > avg_fin || max_ini > max_fin) {
    std::cerr << "Dihedral handle ratio decreased" << endl;
    return 1;
  }
  
  return 0;
}

static inline double
da( double dot )
{ return 180 - (180/M_PI)*acos( dot ); }

void tet_dihedral_angle_ratios( Mesh& mesh,
                                double& ratio_min,
                                double& ratio_avg,
                                double& ratio_max,
                                MsqError& err )
{
  std::vector<Mesh::VertexHandle> verts;
  std::vector<Mesh::ElementHandle> tets;
  std::vector<MsqVertex> coords(4);
  std::vector<size_t> junk;
  mesh.get_all_elements( tets, err );
  
  ratio_min = HUGE_VAL;
  ratio_max = -HUGE_VAL;
  ratio_avg = 0;
  size_t count = 0;
  
  for (std::vector<Mesh::ElementHandle>::iterator i = tets.begin();
       i != tets.end(); ++i) {
       
    Mesh::ElementHandle e = *i;
    EntityTopology type;
    mesh.elements_get_topologies( &e, &type, 1, err );
    assert(!err);
    if (type != TETRAHEDRON)
      continue;
       
    verts.clear();
    mesh.elements_get_attached_vertices( &e, 1, verts, junk, err );
    assert(!err);
    assert(verts.size() == 4);
    mesh.vertices_get_coordinates( arrptr(verts), arrptr(coords), 4, err );
    assert(!err);
    
    Vector3D v01 = coords[1] - coords[0];
    Vector3D v02 = coords[2] - coords[0];
    Vector3D v31 = coords[1] - coords[3];
    Vector3D v32 = coords[2] - coords[3];

    Vector3D n012 = ~(v02 * v01);
    Vector3D n013 = ~(v31 * v01);
    Vector3D n023 = ~(v02 * v32);
    Vector3D n123 = ~(v31 * v32);
  
    double ds[] = { da(n012 % n013),
                    da(n012 % n123),
                    da(n012 % n023),
                    da(n013 % n023),
                    da(n013 % n123),
                    da(n023 % n123) };
    
    double ratio = *std::min_element( ds, ds+6 ) / *std::max_element( ds, ds+6 );
    if (ratio < ratio_min)
      ratio_min = ratio;
    if (ratio > ratio_max)
      ratio_max = ratio;
    ratio_avg += ratio;
    count++;
  }
  
  if (count)
    ratio_avg /= count;
}
