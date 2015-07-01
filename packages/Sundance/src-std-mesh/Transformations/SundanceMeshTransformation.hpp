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

#ifndef SUNDANCE_MESHFILTER_H
#define SUNDANCE_MESHFILTER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshTransformationBase.hpp"
#include "PlayaHandle.hpp"

namespace Sundance
{
  /**
   * MeshTransformation is the user-level interface for mesh filters, i.e.,
   * objects that take an input mesh and produce a new mesh. Examples
   * of filter operations are refinement, load balancing,
   * and extrusion from 2D to 3D. 
   *
   * <h4> Example: </h4> extrude a 2D mesh into 2D
   * \code
   * // create a 2D mesh 
   * MeshType meshType = new BasicSimplicialMeshType();
   * MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, 10, 1,
   *                                                    0.0, 1.0, 10, 1,
   *                                                    meshType);
   * Mesh mesh2D = mesher.getMesh();
   * // create a filter for extruding 2 levels between z=0.0 and z=0.2
   * MeshTransformation extruder = new ExtrusionMeshTransformation(0.0, 0.2, 2);
   * // perform the extrusion
   * Mesh mesh3D = extruder.apply(mesh2D);
   * \endcode
   */
  class MeshTransformation : public Playa::Handle<MeshTransformationBase>
  {
  public:
    /** Construct an empty mesh filter object */
    MeshTransformation();

    /** Construct from a raw pointer to a mesh filter subtype */
    MeshTransformation(Playa::Handleable<MeshTransformationBase>* rawPtr);

    /** Construct from a smart pointer to a mesh filter subtype */
    MeshTransformation(const RCP<MeshTransformationBase>& smartPtr);

    /** apply the filter to create a new mesh */
    Mesh apply(const Mesh& inputMesh) const ;

    const bool& serializeLocal() const {return serializeLocal_;}

    bool& serializeLocal() {return serializeLocal_;}
  private:
    bool serializeLocal_;
    
  };
}

#endif
