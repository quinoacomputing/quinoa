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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
  
#include "Mesquite_all_headers.hpp"
#include <ostream>
using namespace Mesquite;
int main(int argc, char* argv[])
{
  MsqError err;
  
  if (argc != 2) {
    std::cerr << "Expected mesh file names as single argument." << std::endl;
    exit (EXIT_FAILURE);
  }

  // new code starts here
  //... 
  Mesquite::MeshImpl my_mesh;
  my_mesh.read_vtk(argv[1], err);
  if (err)
  {

    std::cout << err << std::endl;
    return 1;
  }

  my_mesh.write_vtk("original_mesh.vtk",err);

  Vector3D normal(0,0,1);
  Vector3D point(0,0,5);
  PlanarDomain my_mesh_plane(normal, point);

  Mesquite::ShapeImprover mesh_quality_algorithm;
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&my_mesh, &my_mesh_plane);
  mesh_quality_algorithm.run_instructions(&mesh_and_domain, err);
    //Should check the error object after the instruction is ran
    // to see whether the instructions were all successful.
  if (err)
  {
    std::cout << err << std::endl;
    return 1;
  }

  my_mesh.write_vtk("smoothed_mesh.vtk",err);

  return 0;
}
