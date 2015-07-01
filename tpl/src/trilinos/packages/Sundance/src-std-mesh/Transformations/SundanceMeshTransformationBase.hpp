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

#ifndef SUNDANCE_MESHFILTERBASE_H
#define SUNDANCE_MESHFILTERBASE_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshType.hpp"
#include "PlayaHandleable.hpp"
#include "Teuchos_Describable.hpp"
#include "PlayaPrintable.hpp"
#include "SundanceNoncopyable.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceIncrementallyCreatableMesh.hpp"

namespace Sundance
{
/**
 * MeshSourceBase provides the internal interface for mesh filters, i.e.,
 * objects that take an input mesh and produce a new mesh. Examples
 * of filter operations are refinement, load balancing,
 * and extrusion from 2D to 3D. 
 * The action of a mesh filter should be independent
 * of the internal mesh representation used. To allow user-level
 * specification of the type of internal mesh representation to be
 * used, a MeshTransformationBase is constructed with a MeshType object
 * which acts as a factory to produce empty output meshes.
 *
 * If the
 * communicator has more than one processor, the mesh created will
 * be distributed.
 *
 * <h4> Writing your own MeshTransformationBase subtype </h4>
 *
 * The only method you will need to override is
 * <ul>
 * <li> <tt>virtual Mesh apply(const Mesh& inputMesh) const = 0 </tt> 
 * </ul>
 * which is where you do the filter action and return an output
 * mesh. This method
 * should usually physically create the mesh with a call to createMesh(),
 * ensuring that the correct mesh representation type is created
 * using the MeshType factory with which the filter was constructed.
 *
 * See the ExtrustionMeshTransformation source code for a very simple
 * example of how to write a mesh filter subtype. 
 *
 * Optionally, you can override the description() method to 
 * provide more informative descriptive output than the std::string
 * <tt>"MeshTransformationBase[unknown subtype]".</tt>
 */
class MeshTransformationBase : public Playa::Handleable<MeshTransformationBase>,
                               public Playa::Printable,
                               public Teuchos::Describable,
                               public Noncopyable,
                               public ObjectWithClassVerbosity<MeshTransformationBase>
{
public:
  /** Construct with a mesh type, which specifies
   *  the type of mesh to be built when the filter is applied. */
  MeshTransformationBase(const MeshType& meshType)
    : meshType_(meshType) {;}

  /** virtual dtor */
  virtual ~MeshTransformationBase(){;}

      
  /** Apply the filter to the given input mesh, 
   *  producing an output mesh */
  virtual Mesh apply(const Mesh& inputMesh) const = 0 ;

  /** \name Printable interface */
  //@{
  /** Print to a stream */
  virtual void print(std::ostream& os) const {os << description();}
  //@}

  /** \name Describable interface */
  //@{
  /** Print to a stream */
  virtual std::string description() const 
    {return "MeshTransformationBase[unknown subtype]";}
  //@}
      
protected:

  /** createMesh() allocates the mesh object with a call to 
   * meshType's createMesh() method. */
  Mesh createMesh(int dim, const MPIComm& comm) const ;

private:
  /** */
  MeshType meshType_;


};

}


#endif
