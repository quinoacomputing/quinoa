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

#ifndef SUNDANCE_MESHSOURCE_H
#define SUNDANCE_MESHSOURCE_H

#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"

#include "PlayaHandle.hpp"

namespace Sundance
{
  /**
   * MeshSource is the user-level interface for objects such as
   * mesh generators and mesh file readers. A MeshSource can
   * create a mesh object with the <tt> getMesh() </tt> method,
   * and if node and element attributes are available, it can
   * access them with the getAttributes() method.
   *
   * Example: read input from a file celled "meshFile" 
   * in Shewchuk's Triangle format, and
   * create a mesh of type BasicSimplicialMesh distributed over the
   * MPI communicator MPI_COMM_WORLD.
   * \code
   * MeshType meshType = new BasicSimplicialMeshType();
   * MeshSource meshSrc = new TriangleMeshReader("meshFile", meshType, MPIComm::world());
   * \endcode
   */
  class MeshSource : public Playa::Handle<MeshSourceBase>
  {
  public:
    /** Construct an empty mesh source object */
    MeshSource();

    /** Construct from a raw pointer to a mesh source subtype */
    MeshSource(Playa::Handleable<MeshSourceBase>* rawPtr);

    /** Construct from a smart pointer to a mesh source subtype */
    MeshSource(const RCP<MeshSourceBase>& smartPtr);

    /** Create and return a mesh */
    Mesh getMesh() const ;

    /** Get any attributes associated with the nodes and elements in the
     * mesh. If no attributes exist, the arrays are empty. If the mesh
     * does not exist, it will be created with a cell to getMesh(). */
    void getAttributes(RCP<Array<Array<double> > >& nodeAttributes,
                       RCP<Array<Array<double> > >& elemAttributes) const ;

    /** Return the mesh type to be used by default if no MeshType
     * is given in a MeshSource subtype ctor. The default mesh type
     * can be set by including a specifer such as
     * <pre>
     * <DefaultMesh type="BasicSimplicial"/>
     * </pre>
     * as a child in the XML configuration file. 
     */
    static MeshType& defaultMeshType() ;

    /** access to the MPI communicator */
    const MPIComm& comm() const ;

    static bool& staggerOutput() {static bool rtn=false; return rtn;}

  private:
  };
}

#endif
