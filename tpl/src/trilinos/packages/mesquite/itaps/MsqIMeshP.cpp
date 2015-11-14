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
 
    (c) 2009 isenburg@cs.unc.edu   kraftche@cae.wisc.edu
   
  ***************************************************************** */
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 26-March-08 at 10:26:21
//  LAST-MOD: 15-Nov-04 by kraftche@cae.wisc.edu
//
/*! \file ParallelMeshImpl.cpp

\brief This files contains a parallel mesh implementation that can be used
to run mesquite by default.

    \author Jason Kraftcheck
    \author Martin Isenburg
    \date 2008-3-26
 */

#include "MsqIMeshP.hpp"
#include "MsqError.hpp"
#include "MsqDebug.hpp"
#include <stdlib.h>

namespace MESQUITE_NS
{

  MsqIMeshP::MsqIMeshP( iMesh_Instance mesh, iMeshP_PartitionHandle partition,
			iBase_EntitySetHandle meshset, iBase_EntityType type,
                        MsqError& err,
			const iBase_TagHandle* fixed_tag,
			const iBase_TagHandle* slaved_tag )
    : MsqIMesh(mesh, meshset, type, err, fixed_tag, slaved_tag), 
      partitionInstance(partition)
  {
  }

  MsqIMeshP::MsqIMeshP( iMesh_Instance mesh, iMeshP_PartitionHandle partition,
                        iBase_EntityType element_dimension,
                        MsqError& err,
			const iBase_TagHandle* fixed_tag,
			const iBase_TagHandle* slaved_tag  )
    : MsqIMesh(mesh, element_dimension, err, fixed_tag, slaved_tag), 
      partitionInstance(partition)
  {
  }

  MsqIMeshP::~MsqIMeshP()
  {
  }

//**************** Parallel Methods ******************************

void MsqIMeshP::vertices_get_global_id(const VertexHandle vert_array[],
					      size_t gid[],
					      size_t num_vtx,
					      MsqError& err)
{
  int itaps_err;
  // get a local part id
  iMeshP_PartHandle *parts = 0;
  int parts_allocated = 0;
  int parts_size = 0;
  iMeshP_getLocalParts(meshInstance, partitionInstance, &parts, &parts_allocated, &parts_size, &itaps_err);
  iMeshP_Part part_id;
  iMeshP_getPartIdFromPartHandle(meshInstance, partitionInstance, parts[0], &part_id, &itaps_err);
  if (parts_allocated) free(parts);
  // get rank of local part
  int rank;
  iMeshP_getRankOfPart(meshInstance, partitionInstance, part_id, &rank, &itaps_err);
  // get global ids for all vertex handles
  for (unsigned i = 0; i < num_vtx; i++)
  {
    iMeshP_getEntOwnerPart(meshInstance, partitionInstance, (iBase_EntityHandle)(vert_array[i]), &part_id, &itaps_err);
    int part_rank;
    iMeshP_getRankOfPart(meshInstance, partitionInstance, part_id, &part_rank, &itaps_err);
    if (part_rank == rank) {
      gid[i] = (size_t)(vert_array[i]);
    }
    else {
      iBase_EntityHandle handle;
      iMeshP_getOwnerCopy(meshInstance, partitionInstance, (iBase_EntityHandle)(vert_array[i]), &part_id, &handle, &itaps_err);
      gid[i] = (size_t)handle;
    }
  }
}

void MsqIMeshP::vertices_get_processor_id(const VertexHandle vert_array[],
						 int pid[],
						 size_t num_vtx,
						 MsqError& err)
{
  int itaps_err;
  for (unsigned i = 0; i < num_vtx; i++)
  {
    iMeshP_Part part_id;
    iMeshP_getEntOwnerPart(meshInstance, partitionInstance, (iBase_EntityHandle)(vert_array[i]), &part_id, &itaps_err);
    iMeshP_getRankOfPart(meshInstance, partitionInstance, part_id, &(pid[i]), &itaps_err);
  }
}

} // namespace Mesquite
