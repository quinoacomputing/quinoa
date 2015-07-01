/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef SUNDANCE_INCREMENTALLYCREATABLEMESH_H
#define SUNDANCE_INCREMENTALLYCREATABLEMESH_H



#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
/**
 * IncrementallyCreatableMesh is an interface for the creation of meshes
 * by adding vertices, then elements.
 */
class IncrementallyCreatableMesh : public MeshBase
{
public:
  /** Construct an empty mesh of the given dimension, distributed over
   * processors in the MPI communicator*/
  IncrementallyCreatableMesh(int dim, const MPIComm& comm, 
    const MeshEntityOrder& meshOrder)
    : MeshBase(dim, comm, meshOrder) {;}

  /** virtual dtor */
  virtual ~IncrementallyCreatableMesh(){;}

  /** Optional preallocation of space for an estimated number of vertices */
  virtual void estimateNumVertices(int nPts) {;}

  /** Optional preallocation of space for an estimated number of elements */
  virtual void estimateNumElements(int nElems) {;}

      

  /** 
   * Add new new vertex to the mesh.
   * \param globalIndex the GID of the new vertex
   * \param x the spatial position of the new vertex
   * \param ownerProcID the processor that "owns" this vertex 
   * \param label a label for this vertex (optionally used in setting loads, boundary
   * conditions, etc)
   * \return the LID of the vertex.
   */
  virtual int addVertex(int globalIndex, const Point& x,
    int ownerProcID, int label) = 0 ;

  /** 
   * Add a new element to the mesh.
   * \param globalIndex the GID of the new element
   * \param vertexGIDs tuple of GIDs for the vertices defining this element. 
   * \param ownerProcID the processor that "owns" this element
   * \param label a label for this element (optionally used in setting loads, 
   * material properties, etc)
   * \return the LID of the element
   */
  virtual int addElement(int globalIndex, const Array<int>& vertexGIDs,
    int ownerProcID, int label) = 0 ;


  /** 
   * 
   */
  virtual void freezeTopology() {}
};

}



#endif
