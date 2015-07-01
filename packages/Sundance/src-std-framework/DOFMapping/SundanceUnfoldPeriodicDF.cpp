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

#include "SundanceMap.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceUnfoldPeriodicDF.hpp"
#include "SundancePeriodicMesh1D.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceMeshSource.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceDOFMapBase.hpp"
#include "PlayaMPIContainerComm.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "PlayaLoadableVector.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorSpaceImpl.hpp"
#include "PlayaVectorImpl.hpp"
#endif


namespace Sundance
{
using namespace Teuchos;
using namespace Playa;

Mesh unfoldPeriodicMesh(const Mesh& mesh)
{
  const MeshBase* mb = mesh.ptr().get();
  const PeriodicMesh1D* pm = dynamic_cast<const PeriodicMesh1D*>(mb);

  TEUCHOS_TEST_FOR_EXCEPT(pm==0);

  int numElems = mesh.numCells(1);
  double a = mesh.nodePosition(0)[0];
  double b = mesh.nodePosition(numElems)[0];

  MeshType meshType = new BasicSimplicialMeshType();
  int verb=0;
  MeshSource mesher = new PartitionedLineMesher(a, b, numElems, meshType, verb,
    MPIComm::self());

  Mesh rtn = mesher.getMesh();

  return rtn;
}


DiscreteSpace unfoldPeriodicDiscreteSpace(const DiscreteSpace& space)
{
  VectorType<double> vecType = space.vecType();
  Mesh mesh = unfoldPeriodicMesh(space.mesh());
  BasisArray basis = space.basis();

  return DiscreteSpace(mesh, basis, vecType);
}

Expr unfoldPeriodicDiscreteFunction(const Expr& f, const string& name)
{
  const DiscreteFunction* df = DiscreteFunction::discFunc(f);
  TEUCHOS_TEST_FOR_EXCEPT(df==0);


  DiscreteSpace perSpace = df->discreteSpace();
  DiscreteSpace space = unfoldPeriodicDiscreteSpace(perSpace);
  
  Mesh oldMesh = perSpace.mesh();
  Mesh newMesh = space.mesh();

  df->updateGhosts();
  RCP<GhostView<double> > ghostView = df->ghostView();

  Vector<double> newVec = space.createVector();


  const RCP<DOFMapBase>& oldMap = perSpace.map();
  const RCP<DOFMapBase>& newMap = space.map();

  /* Copy the element DOFs to the new vector. There are no duplicated
   * elements so this map is one-to-one */
  Array<int> oldCellLID(oldMesh.numCells(1));
  for (int i=0; i<oldMesh.numCells(1); i++)
  {
    oldCellLID[i] = i;
  }
  RCP<const Set<int> > funcs = oldMap->allowedFuncsOnCellBatch(1, oldCellLID);

  Array<Array<int> > oldElemDofs;
  Array<Array<int> > newElemDofs;
  Array<int> oldNNodes;
  RCP<const MapStructure> oldMapStruct ;
  RCP<const MapStructure> newMapStruct ;

  if (funcs->size())
  {
    oldMapStruct 
      = oldMap->getDOFsForCellBatch(1, oldCellLID, *funcs, oldElemDofs, 
        oldNNodes, 0);
    
    newMapStruct 
      = newMap->getDOFsForCellBatch(1, oldCellLID, *funcs, newElemDofs, 
        oldNNodes, 0);
    
    for (int chunk=0; chunk<oldMapStruct->numBasisChunks(); chunk++)
    {
      int nf = oldMapStruct->numFuncs(chunk);
      int nDofsPerCell = oldNNodes[chunk] * nf;
      for (int c=0; c<oldCellLID.size(); c++)
      {
        for (int d=0; d<nDofsPerCell; d++)
        {
          int oldDof = oldElemDofs[chunk][nDofsPerCell*c + d];
          int newDof = newElemDofs[chunk][nDofsPerCell*c + d];
          loadable(newVec)->setElement(newDof, ghostView->getElement(oldDof));
        }
      }
    }
  }

  /* Copy the vertex dofs to the new vector. The data at v=0 are duplicated */
  Array<int> oldVertLID(oldMesh.numCells(0));
  Array<int> newVertLID(newMesh.numCells(0));
  for (int i=0; i<oldMesh.numCells(0); i++)
  {
    oldVertLID[i] = i;
  }
  for (int i=0; i<newMesh.numCells(0); i++)
  {
    newVertLID[i] = i;
  }
  if (funcs->size())
  {
    funcs = oldMap->allowedFuncsOnCellBatch(0, oldCellLID);
    
    oldMapStruct = oldMap->getDOFsForCellBatch(0, oldVertLID, *funcs, 
      oldElemDofs, 
      oldNNodes, 0);
    newMapStruct 
      = newMap->getDOFsForCellBatch(0, newVertLID, *funcs, newElemDofs, 
        oldNNodes, 0);
    
    for (int chunk=0; chunk<oldMapStruct->numBasisChunks(); chunk++)
    {
      int nf = oldMapStruct->numFuncs(chunk);
      int nDofsPerCell = oldNNodes[chunk] * nf;
      for (int c=0; c<newVertLID.size(); c++)
      {
        if (c < oldVertLID.size())
        {
          for (int d=0; d<nDofsPerCell; d++)
          {
            int oldDof = oldElemDofs[chunk][nDofsPerCell*c + d];
            int newDof = newElemDofs[chunk][nDofsPerCell*c + d];
            loadable(newVec)->setElement(newDof, ghostView->getElement(oldDof));
          }
        }
        else
        {
          for (int d=0; d<nDofsPerCell; d++)
          {
            int oldDof = oldElemDofs[chunk][d];
            int newDof = newElemDofs[chunk][nDofsPerCell*c + d];
            loadable(newVec)->setElement(newDof, ghostView->getElement(oldDof));
          }
        }
      }
    }
  }

  return new DiscreteFunction(space, newVec, name);
}


}


