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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 14-Nov-02 at 17:33:19 by Thomas Leurent
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

#include "Mesquite.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_MsqError.hpp"
#include "Mesquite_PlanarDomain.hpp"
// algorythms
#include "Mesquite_LaplaceWrapper.hpp"


#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

using namespace Mesquite;


int main()
{
  Mesquite::MeshImpl mesh;
  MsqPrintError err(cout);
  mesh.read_vtk(MESH_FILES_DIR "2D/vtk/quads/untangled/square_quad_2.vtk", err);
  if (err) return 1;
  
     //create geometry: plane z=0, normal (0,0,1)
  Vector3D pnt(0,0,5);
  Vector3D s_norm(0,0,1);
  Mesquite::PlanarDomain msq_geom(s_norm, pnt);
  
    // creates an intruction queue
  LaplaceWrapper laplacian_smoother;
  
  mesh.write_vtk("original_mesh.vtk", err); 
  if (err) return 1;
  
    // launches optimization on mesh_set1
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &msq_geom);
  laplacian_smoother.run_instructions(&mesh_and_domain, err); 
  if (err) return 1;
 
  mesh.write_vtk("smoothed_mesh.vtk", err); 
  if (err) return 1;
}
