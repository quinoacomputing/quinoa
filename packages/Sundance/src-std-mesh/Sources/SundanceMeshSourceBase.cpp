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


#include "SundanceMeshSourceBase.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;


MeshSourceBase::MeshSourceBase(const MeshType& meshType,
  int verbosity,
  const MPIComm& comm)
  : ObjectWithVerbosity(verbosity),
    cachedMesh_(),
    hasCachedMesh_(),
    meshType_(meshType),
    comm_(comm),
    elemAttributes_(rcp(new Array<Array<double> >())),
    nodeAttributes_(rcp(new Array<Array<double> >()))
{
}

MeshSourceBase::MeshSourceBase(const ParameterList& params)
  : cachedMesh_(),
    hasCachedMesh_(),
    meshType_(new BasicSimplicialMeshType()),
    comm_(MPIComm::world()),
    elemAttributes_(rcp(new Array<Array<double> >())),
    nodeAttributes_(rcp(new Array<Array<double> >()))
{
  
}

Mesh MeshSourceBase::getMesh() const
{
  Tabs tabs;
  
  /* if we don't have a cached mesh, build one */
  if (!hasCachedMesh_)
    {
      Mesh rtn =  fillMesh();
      if (verb() > 0)
        {
          std::cerr << tabs << "got a mesh with " << rtn.numCells(0)
               << " nodes and " << rtn.numCells(rtn.spatialDim())
               << " maximal cells" << std::endl;
        }
      return rtn;
    }
  return cachedMesh_;
}

void MeshSourceBase
::getAttributes(RCP<Array<Array<double> > >& nodeAttributes,
                RCP<Array<Array<double> > >& elemAttributes) const
{
  nodeAttributes = nodeAttributes_;
  elemAttributes = elemAttributes_;
}

Mesh MeshSourceBase::createMesh(int dim) const
{
  cachedMesh_ = meshType_.createEmptyMesh(dim, comm_);
  cachedMesh_.ptr()->setVerb(verb());
  hasCachedMesh_ = true;
  
  return cachedMesh_;
}
