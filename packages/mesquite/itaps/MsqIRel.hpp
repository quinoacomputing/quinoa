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


/** \file MsqIRel.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_MSQ_IREL_HPP
#define MSQ_MSQ_IREL_HPP

#include "Mesquite.hpp"
#include "MsqIGeom.hpp"
#include "iRel.h"

namespace MESQUITE_NS {


/* General MeshDomain on iGeom & iRel implementation */

class MsqIRel : public MsqCommonIGeom
{
public:

  MsqIRel( iGeom_Instance geom,
           iRel_Instance irel_iface,
           iRel_PairHandle irel_instance );

  virtual ~MsqIRel();

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
                      
protected:

    /** Get geometric entity owning a mesh entity */
  int geom_from_mesh( Mesh::EntityHandle  mesh_handle_in,
                      iBase_EntityHandle& geom_handle_out,
                      int& geom_dimension_out ) const;
  
  int geom_from_mesh( Mesh::EntityHandle const* mesh_handles_in,
                      iBase_EntityHandle      * geom_handles_out,
                      unsigned short* geom_dimensions_out,
                      size_t count ) const;

private:

    /** ITAPS interface implementation for mesh->geometry association */
  iRel_Instance  relateIface;
  iRel_PairHandle relateInstance;
  
    /** temporary storage of geometry entity handles */
  mutable std::vector<iBase_EntityHandle> geomHandles;
};
} // namespace Mesquite

#endif
