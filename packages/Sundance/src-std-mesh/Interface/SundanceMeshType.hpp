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

#ifndef SUNDANCE_MESHTYPE_H
#define SUNDANCE_MESHTYPE_H

#include "SundanceDefs.hpp"
#include "SundanceMeshTypeBase.hpp"
#include "SundanceMesh.hpp"
#include "PlayaHandle.hpp"

namespace Sundance
{
/**
 * Class MeshType is a user-level object for specification of which
 * internal mesh representation is to be used when building or reading
 * a mesh. An example of using a MeshType to control the creation 
 * of a mesh with a TriangleMeshReader is as follows: 
 * \code
 * MeshType meshType = new BasicSimplicialMeshType();
 * MeshSource meshSrc = new TriangleMeshReader("meshFile", meshType, MPIComm::world());
 * \endcode
 * The internal representation of the mesh will be as a BasicSimplicialMesh
 * object. 
 */
class MeshType : public Playa::Handle<MeshTypeBase>
{
public:
  /** Construct an empty mesh type object */
  MeshType();

  /** Construct from a raw pointer to a mesh type subtype */
  MeshType(Playa::Handleable<MeshTypeBase>* rawPtr);

  /** Construct from a smart pointer to a mesh type subtype */
  MeshType(const RCP<MeshTypeBase>& smartPtr);

  /** Create a mesh of the given dimension */
  Mesh createEmptyMesh(int dim, const MPIComm& comm) const ;
    
};
}

#endif
