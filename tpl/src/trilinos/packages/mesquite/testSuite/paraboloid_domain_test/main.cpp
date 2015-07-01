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
//   SUMMARY: Tests a paraboloid mesh for compatibility with the ParaboloidDomain.
//     USAGE:
//
// ORIG-DATE: 25-Jan-2013
//
//    AUTHOR: Boyd Tidwell
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

#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "MeshInterface.hpp"

using namespace Mesquite;


class ParaboloidDomain : public MeshDomain
{
  public:

    virtual void closest_point( Mesh::VertexHandle handle,
                                const Vector3D& position,
                                Vector3D& closest,
                                Vector3D& normal,
                                MsqError& err ) const 
    {
      closest = Vector3D(position[0], position[1], position[0]*position[0] + position[1]*position[1]);
    };

    virtual void snap_to(Mesh::VertexHandle entity_handle,
                         Vector3D &coordinate) const {};
    
    virtual void vertex_normal_at(Mesh::VertexHandle entity_handle,
                                  Vector3D &coordinate) const {};
    virtual void element_normal_at(Mesh::ElementHandle entity_handle,
                                   Vector3D &coordinate) const {};
                          
    virtual void vertex_normal_at( const Mesh::VertexHandle* handles,
                                   Vector3D coordinates[],
                                   unsigned count,
                                   MsqError& err ) const {};
                                
     virtual void domain_DoF( const Mesh::EntityHandle* handle_array,
                             unsigned short* dof_array,
                             size_t num_handles,
                             MsqError& err ) const {};
                             
       

};



int main()
{     
  MsqPrintError err(cout);
  Mesquite::MeshImpl mesh;
  mesh.read_vtk(MESH_FILES_DIR "2D/vtk/quads/untangled/paraboloid.vtk", err); 
//  mesh.read_vtk("/home/bktidwell/tmp/paraboloid.vtk", err);
  if (err) return 1;
  
  ParaboloidDomain domain;

  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &domain, true, true, false);

  std::cout << "Paraboloid Domain Test Passes" << std::endl;


  return 0;
}
