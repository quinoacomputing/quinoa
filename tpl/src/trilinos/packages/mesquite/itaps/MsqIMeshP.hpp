/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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

    (2009) kraftche@cae.wisc.edu   isenburg@cs.unc.edu 

  ***************************************************************** */

/** \file MsqIMeshP.hpp
 *  \brief Adaptor for ITAPS iMesh interface
 *  \author Jason Kraftcheck 
 *  \author Martin Isenburg
 */

#ifndef MSQ_I_MESH_P_HPP
#define MSQ_I_MESH_P_HPP

#include "iMeshP.h"
#include "MsqIMesh.hpp"
#include "ParallelMeshInterface.hpp"

namespace MESQUITE_NS {

/**\class MsqIMeshP
 *\brief Mesquite iMesh Adapter
 *
 * Adpater for interfacing parallel Mesquite with an application
 * that provides the ITAPS iMeshP interface for interacting with
 * mesh data.
 */
  class MsqIMeshP : public MsqIMesh, virtual public ParallelMesh
{
public:
//********* Functions that are NOT inherited ************

  MsqIMeshP();
  virtual ~MsqIMeshP();

  MsqIMeshP(iMesh_Instance imesh, iMeshP_PartitionHandle partition,
	   iBase_EntitySetHandle meshset, iBase_EntityType element_dimension,
           MsqError& err,
	   const iBase_TagHandle* fixed_tag = 0,
	   const iBase_TagHandle* slaved_tag= 0 );
  
  MsqIMeshP(iMesh_Instance imesh, iMeshP_PartitionHandle partition,
            iBase_EntityType element_dimension, 
            MsqError& err,
	    const iBase_TagHandle* fixed_tag = 0,
	    const iBase_TagHandle* slaved_tag= 0 );
  
//********** Inherited Functions from ParallelMesh ****************

  /*! Get global ids for given vertices.
   */
  virtual void vertices_get_global_id(const VertexHandle vert_array[],
				      size_t gid[],
				      size_t num_vtx,
				      MsqError& err);
         
  /*! Get processor ids for given vertices.
   */
  virtual void vertices_get_processor_id(const VertexHandle vert_array[],
					 int pid[],
					 size_t num_vtx,
					 MsqError& err);
     
protected:

  iMeshP_PartitionHandle partitionInstance;

};

} // namespace Mesquite

#endif
