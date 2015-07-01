
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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */
/*!
  \file   MsqIGeom.hpp
  \brief  Mesquite::MeshDomain implemented on ITAPS iGeom API
  \author Jason Kraftcheck
  \date   2007-08-14
*/

#ifndef MSQ_IGEOM_HPP
#define MSQ_IGEOM_HPP

#include "Mesquite.hpp"
#include "MeshInterface.hpp"
#include "iGeom.h"
#include <vector>

namespace MESQUITE_NS
{

/**\brief Common code for specific implementations of MeshDomain on ITAPS interfaces.
 *
 * This class contains the common functionality used by concrete implementations
 * of MeshDomain on the ITAPS geometry interface.
 */
class MsqCommonIGeom : public MeshDomain
{
public:

    /**\param geom The ITAPS geometry interface implementation to query */
  MsqCommonIGeom( iGeom_Instance geom );
  
  virtual ~MsqCommonIGeom();

    /** Evaluate the closest point to the input position on the specified
     *  geometric entity and return the result in the passed position 
     *  argument (move the passed position onto the geometry.)
     */
  int move_to( iBase_EntityHandle geom_handle, Vector3D& coord ) const;
  
    /** Given a geometric entity and a position, evaluate the normal 
     *  on the geometric entity at the closest point on that entity
     *  to the input position, and pass back the result in the input
     *  coord vector.
     */
  int normal ( iBase_EntityHandle geom_handle, Vector3D& coord ) const ;
  
    /** Given a geometric entity and a position, evaluate the normal 
     *  on the geometric entity at the closest point on that entity
     *  to the input position, and pass back the result in the input
     *  coord vector.
     */
  int normal ( iBase_EntityHandle geom_handle, Vector3D coords[], unsigned count ) const;
  
    /** Given a geometric entity and a position, evaluate the normal 
     *  on the geometric entity at the closest point on that entity
     *  to the input position, and pass back the result in the input
     *  coord vector.
     */
  int normal ( const iBase_EntityHandle geom_handles[], Vector3D coords[], unsigned count ) const;
    
    /** Given a geometric entity and a position, get point on 
     *  the geometric entity closest to the input position, and
     *  the surface normal at that position.
     */
  int closest_and_normal( iBase_EntityHandle geom_handle,
                           const Vector3D& position,
                           Vector3D& closest, 
                           Vector3D& normal ) const;
                         
  int get_dimension( iBase_EntityHandle const* geom_handle, 
                      unsigned short* dof_out,
                      size_t count ) const ;
                      
  iGeom_Instance geomIFace;

private:  
  mutable std::vector<iBase_EntityHandle> geomHandles;
  mutable std::vector<double> coordArray;
  mutable std::vector<int> typeArray;
};


/**\brief A Mesquite::MeshDomain implemented on top of the ITAPS iGeom API.
 *
 * Simple MeshDomain class implementatation that queries a single iGeom
 * entity for all geometric queries.  Suitable for use when the entire
 * mesh to be smoothed lies on a single geometric surface.
 */
class MsqIGeom : public MsqCommonIGeom
{
public:

  MsqIGeom( iGeom_Instance geom,
            iBase_EntityHandle geom_ent_handle );

  virtual ~MsqIGeom();

  void snap_to( Mesh::VertexHandle entity_handle,
                Vector3D& coordinat ) const;

  void vertex_normal_at( Mesh::VertexHandle entity_handle,
                         Vector3D& coordinate ) const;

  void element_normal_at( Mesh::ElementHandle entity_handle,
                          Vector3D& coordinate ) const;
  
  void vertex_normal_at( const Mesh::VertexHandle* handles,
                         Vector3D coordinates[],
                         unsigned count,
                         MsqError& err ) const;

  void closest_point( Mesh::VertexHandle handle,
                      const Vector3D& position,
                      Vector3D& closest,
                      Vector3D& normal,
                      MsqError& err ) const;

  void domain_DoF( const Mesh::VertexHandle* handle_array,
                   unsigned short* dof_array,
                   size_t num_vertices,
                   MsqError& err ) const;
private:
  
    /** A handle for the geometry entity to evaluate */
  iBase_EntityHandle geomEntHandle;
};

}

#endif
