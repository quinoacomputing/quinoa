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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Sep-03 at 10:57:52
//  LAST-MOD: 21-Sep-03 at 15:03:52 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

Calls the Mesquite wrappers. First command line argument is the mesh file.

 */
// DESCRIP-END.
//
#include "meshfiles.h"

#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

// algorythms
#include "ShapeImprovementWrapper.hpp"
#include "MsqTimer.hpp"

using namespace Mesquite;

int main(int argc, char* argv[])
{
  Mesquite::MsqPrintError err(cout);
  char file_name[256];
  
  // command line arguments
  if (argc==1) {
	
     cout << "No command line argument given for mesh file.\n"
          << "Defaulting to " MESH_FILES_DIR "3D/vtk/tets/untangled/tire.vtk\n";
     strcpy(file_name,MESH_FILES_DIR "3D/vtk/tets/untangled/tire.vtk");
  } 
  else if (argc>2) 
    cout << "Too many command line arguments.\n" << endl;
  else if (argc==2) {
    cout << " given 1 command line arguments.\n";
    strcpy(file_name, argv[1]);
  }
  
  Mesquite::MeshImpl mesh;
  mesh.read_vtk(file_name, err);
  
  // creates a wrapper
  ShapeImprovementWrapper wrapper;

//  mesh->write_vtk("original_mesh",err); MSQ_CHKERR(err);
  
  // launches optimization on mesh_set1
  wrapper.run_instructions(&mesh, err); 
  if (err) return 1;
  
//  mesh->write_vtk("smoothed_mesh", err); MSQ_CHKERR(err);
  print_timing_diagnostics(cout);
  return 0;
}
