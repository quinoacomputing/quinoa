/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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
 
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */

/*!
  \file   MsqIRel.cpp
  \brief  


  \author Jason Kraftcheck
  \date   2007-08-14
*/

#include "MsqIRel.hpp"
#include "MsqDebug.hpp"
#include "MsqVertex.hpp"
#include "MsqError.hpp"
#include "MsqIBase.hpp"

namespace MESQUITE_NS
{


/***************** DomainTSTT class methods *********************/

MsqIRel::MsqIRel( iGeom_Instance geom,
                  iRel_Instance relate_iface,
                  iRel_PairHandle relate_instance ) 
  : MsqCommonIGeom( geom ), 
    relateIface( relate_iface ),
    relateInstance( relate_instance )
{
}

MsqIRel::~MsqIRel() {}


void MsqIRel::snap_to( Mesh::VertexHandle handle,
                           Vector3D& coordinate ) const
{
  int ierr, dim;
  iBase_EntityHandle geom;
  
  ierr = geom_from_mesh( handle, geom, dim );
  if (iBase_SUCCESS != ierr) {
    process_itaps_error( ierr );
    return;
  }
  
  if (dim < 3) {  
    ierr = move_to( geom, coordinate );
    if (iBase_SUCCESS != ierr) {
      process_itaps_error( ierr );
      return;
    }
  }
}

void MsqIRel::vertex_normal_at( Mesh::VertexHandle handle,
                                    Vector3D& coordinate ) const
{
  int ierr, dim;
  iBase_EntityHandle geom;
  
  ierr = geom_from_mesh( handle, geom, dim );
  if (iBase_SUCCESS != ierr) {
    process_itaps_error( ierr );
    return;
  }
  
  if (dim == 2) {
    ierr = normal( geom, coordinate );
    if (iBase_SUCCESS != ierr) {
      process_itaps_error( ierr );
      return;
    }
  }
  else {
    assert(0);
  }
}

void MsqIRel::element_normal_at( Mesh::ElementHandle handle,
                                     Vector3D& coordinate ) const
{
  MsqIRel::vertex_normal_at( handle, coordinate );
}

void MsqIRel::vertex_normal_at( const Mesh::VertexHandle* handle,
                                    Vector3D coordinates[],
                                    unsigned count,
                                    MsqError& err ) const
{
  int ierr, dim;
  iBase_EntityHandle geom;
  for (unsigned i = 0; i < count; ++i) {
    ierr = geom_from_mesh( handle[i], geom, dim );
    if (iBase_SUCCESS != ierr) {
      process_itaps_error( ierr );
      return;
    }

    if (dim != 2) {
      MSQ_SETERR(err)("Cannot get normal for non-surface geometry", MsqError::INVALID_ARG );
      return;
    }

    ierr = normal( geom, coordinates[i] );
    if (iBase_SUCCESS != ierr) {
      process_itaps_error( ierr );
      return;
    }
  }
}

void MsqIRel::domain_DoF( const Mesh::VertexHandle* handle_array,
                              unsigned short* dof_array,
                              size_t count,
                              MsqError& err ) const
{
  int ierr;
  
  geomHandles.resize( count );
  ierr = geom_from_mesh( handle_array, arrptr(geomHandles), dof_array, count );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)(process_itaps_error( ierr ), MsqError::INTERNAL_ERROR);
    return;
  }
}
    


void MsqIRel::closest_point( Mesh::VertexHandle handle,
                                 const Vector3D& position,
                                 Vector3D& closest,
                                 Vector3D& normal,
                                 MsqError& err ) const
{
  int ierr, dim;
  iBase_EntityHandle geom;
  
  ierr = geom_from_mesh( handle, geom, dim );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)(process_itaps_error( ierr ), MsqError::INTERNAL_ERROR);
    return;
  }

  if (dim != 2) {
    MSQ_SETERR(err)("Cannot get normal for non-surface geometry", MsqError::INVALID_ARG );
    return;
  }

  ierr = closest_and_normal( geom, position, closest, normal );
  if (iBase_SUCCESS != ierr) {
    MSQ_SETERR(err)(process_itaps_error( ierr ), MsqError::INTERNAL_ERROR);
    return;
  }
}


int MsqIRel::geom_from_mesh( Mesh::EntityHandle mesh_ent_handle,
                             iBase_EntityHandle& geom_handle,
                             int& geom_dim ) const
{
    // get geometric entity
  int ierr;
  iRel_getEntEntRelation( relateIface,
                             relateInstance,
                             (iBase_EntityHandle)mesh_ent_handle,
                             true,
                             &geom_handle,
                             &ierr );
  if (iBase_SUCCESS != ierr)
    return ierr;
  
    // get dimension of geometric entities
  int one = 1, one_too = 1, *type_ptr = &geom_dim;
  iGeom_getArrType( geomIFace, &geom_handle, 1, &type_ptr, &one, &one_too, &ierr );
  if (iBase_SUCCESS != ierr)
    return ierr;
  
  return iBase_SUCCESS;
}


int MsqIRel::geom_from_mesh( const Mesh::EntityHandle* handles,
                             iBase_EntityHandle* geom_handles,
                             unsigned short* dims,
                             size_t count ) const
{
  int ierr, dim;
  for (size_t i = 0; i < count; ++i) {
    ierr = geom_from_mesh( handles[i], geom_handles[i], dim );
    if (iBase_SUCCESS != ierr)
      return ierr;
    dims[i] = dim;
  }
  
  return iBase_SUCCESS;
}                    


} // namespace Mesquite

