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
  \file   MsqIGeom.cpp
  \brief  


  \author Jason Kraftcheck
  \date   2007-08-14
*/

#include "MsqIGeom.hpp"
#include "MsqDebug.hpp"
#include "MsqVertex.hpp"
#include "MsqError.hpp"
#include "MsqIBase.hpp"

namespace MESQUITE_NS
{




/***************** MsqIGeom class methods *********************/

MsqIGeom::MsqIGeom( iGeom_Instance geom, iBase_EntityHandle geom_ent_handle ) 
  : MsqCommonIGeom( geom ), 
    geomEntHandle( geom_ent_handle )
{
}

MsqIGeom::~MsqIGeom() {}


void MsqIGeom::snap_to( Mesh::VertexHandle,
                        Vector3D& coordinate ) const
{
  int ierr = move_to( geomEntHandle, coordinate );
  if (iBase_SUCCESS != ierr)
    process_itaps_error(ierr);
}

void MsqIGeom::vertex_normal_at( Mesh::VertexHandle,
                                 Vector3D& coordinate ) const
{
  int ierr = normal( geomEntHandle, coordinate );
  if (iBase_SUCCESS != ierr)
    process_itaps_error(ierr);
}

void MsqIGeom::element_normal_at( Mesh::ElementHandle,
                                  Vector3D& coordinate ) const
{
  int ierr = normal( geomEntHandle, coordinate );
  if (iBase_SUCCESS != ierr)
    process_itaps_error(ierr);
}


void MsqIGeom::vertex_normal_at( const Mesh::VertexHandle*,
                                 Vector3D coordinates[],
                                 unsigned count,
                                 MsqError& err ) const
{
  int ierr = normal( geomEntHandle, coordinates, count );
  if (iBase_SUCCESS != ierr)
    MSQ_SETERR(err)(process_itaps_error(ierr), MsqError::INTERNAL_ERROR);
}

void MsqIGeom::closest_point( Mesh::VertexHandle handle,
                              const Vector3D& position,
                              Vector3D& closest,
                              Vector3D& normal,
                              MsqError& err ) const
{
  int ierr = closest_and_normal( geomEntHandle, position, closest, normal );
  if (iBase_SUCCESS != ierr)
    MSQ_SETERR(err)(process_itaps_error(ierr), MsqError::INTERNAL_ERROR);
}

void MsqIGeom::domain_DoF( const Mesh::VertexHandle* ,
                           unsigned short* dof_array,
                           size_t num_vertices,
                           MsqError& err ) const
{
  unsigned short dim;
  int ierr = get_dimension( &geomEntHandle, &dim, 1 );
  if (iBase_SUCCESS != ierr)
    MSQ_SETERR(err)(process_itaps_error(ierr), MsqError::INTERNAL_ERROR);
  std::fill( dof_array, dof_array + num_vertices, dim );
}




/***************** GeomTSTTCommon class methods *********************/

MsqCommonIGeom::MsqCommonIGeom( iGeom_Instance geom )
  : geomIFace( geom )
{
}

MsqCommonIGeom::~MsqCommonIGeom() {}



int MsqCommonIGeom::move_to( iBase_EntityHandle geom, Vector3D& coord ) const
{
  double x, y, z;
  int ierr;
  iGeom_getEntClosestPt( geomIFace, geom, coord[0], coord[1], coord[2], &x, &y, &z, &ierr );
  coord.set( x, y, z );
  return ierr;
}

 
 
int MsqCommonIGeom::normal( iBase_EntityHandle geom, Vector3D& coord ) const
{
  double i, j, k;
  int ierr;
  iGeom_getEntNrmlXYZ( geomIFace, geom, coord[0], coord[1], coord[2], &i, &j, &k, &ierr );
  coord.set( i, j, k );
  return ierr;
}
 
int MsqCommonIGeom::normal( iBase_EntityHandle geom, Vector3D coords[], unsigned count ) const
{
  geomHandles.resize( count, geom );
  return normal( arrptr(geomHandles), coords, count );
}
 
int MsqCommonIGeom::normal( const iBase_EntityHandle* geom_handles, 
                         Vector3D coords[], 
                         unsigned count ) const
{
    // going to assume this in the following reinterpret_cast, so
    // check to make sure it is true
  assert( sizeof(Vector3D) == 3*sizeof(double) );
  
    // copy input coordinates into array
  coordArray.clear();
  coordArray.insert( coordArray.begin(), &coords[0][0], &coords[0][0] + 3*count );
   
    // define junk variables required for ITAPS "consistancy"
  int junk_1 = count*3, junk_2 = count*3;
  double* norm_ptr = &coords[0][0];
  
    // get the normals
  int ierr;
  iGeom_getArrNrmlXYZ( geomIFace, 
                       geom_handles,
                       count,
                       iBase_INTERLEAVED,
                       arrptr(coordArray),
                       count*3,
                       &norm_ptr,
                       &junk_1,
                       &junk_2,
                       &ierr ); 
  
  return ierr;
}

int MsqCommonIGeom::closest_and_normal( iBase_EntityHandle geom, 
                                     const Vector3D& position,
                                     Vector3D& closest,
                                     Vector3D& normal ) const
{
  int ierr;
  iGeom_getEntNrmlPlXYZ( geomIFace, geom, 
                         position[0], position[1], position[2], 
                         &closest[0], &closest[1], &closest[2],
                         &normal[0],  &normal[1],  &normal[2],
                         &ierr );
  return ierr;
}

                         
int MsqCommonIGeom::get_dimension( const iBase_EntityHandle* geom_handle, 
                                unsigned short* dof_out,
                                size_t count ) const
{
  int ierr;
  typeArray.resize( count );
  
    // define junk variables required for ITAPS "consistancy"
  int junk_1 = count, junk_2 = count;
  int* type_ptr = arrptr(typeArray);
  
    // get the types
  iGeom_getArrType( geomIFace, geom_handle, count, &type_ptr, &junk_1, &junk_2, &ierr );
  
    // convert from int to unsigned short
  std::copy( typeArray.begin(), typeArray.end(), dof_out );
  return ierr;
}    

} // namespace Mesquite

