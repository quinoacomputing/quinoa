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
#include "SundanceOrderedTuple.hpp"
#include "SundanceNodalDOFMap.hpp"
#include "SundanceLagrange.hpp"
#include "PlayaMPIContainerComm.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Teuchos;
using Playa::MPIComm;
using Playa::MPIContainerComm;

NodalDOFMap::NodalDOFMap(const Mesh& mesh, 
  int nFuncs, 
  const CellFilter& maxCellFilter,
  int setupVerb)
  : SpatiallyHomogeneousDOFMapBase(mesh, nFuncs, setupVerb),
    maxCellFilter_(maxCellFilter),
    dim_(mesh.spatialDim()),
    nFuncs_(nFuncs),
    nElems_(mesh.numCells(mesh.spatialDim())),
    nNodes_(mesh.numCells(0)),
    nNodesPerElem_(mesh.numFacets(mesh.spatialDim(),0,0)),
    elemDofs_(nElems_ * nFuncs * nNodesPerElem_),
    nodeDofs_(mesh.numCells(0)*nFuncs_, -1),
    structure_(rcp(new MapStructure(nFuncs_, rcp(new Lagrange(1)))))
{
  init();
}


RCP<const MapStructure> 
NodalDOFMap::getDOFsForCellBatch(int cellDim,
  const Array<int>& cellLID,
  const Set<int>& requestedFuncSet,
  Array<Array<int> >& dofs,
  Array<int>& nNodes,
  int verbosity) const
{
  TimeMonitor timer(batchedDofLookupTimer());

  Tabs tab;
  SUNDANCE_MSG3(verbosity, 
               tab << "NodalDOFMap::getDOFsForCellBatch(): cellDim=" << cellDim
               << " cellLID=" << cellLID);


  dofs.resize(1);
  nNodes.resize(1);

  int nCells = cellLID.size();


  if (cellDim == dim_)
    {
      int dofsPerElem = nFuncs_ * nNodesPerElem_;
      nNodes[0] = nNodesPerElem_;
      dofs[0].resize(nCells * dofsPerElem);
      Array<int>& dof0 = dofs[0];
      
      for (int c=0; c<nCells; c++)
        {
          for (int i=0; i<dofsPerElem; i++)
            {
              dof0[c*dofsPerElem + i] = elemDofs_[cellLID[c]*dofsPerElem+i];
            }
        }
    }
  else if (cellDim == 0)
    {
      nNodes[0] = 1;
      dofs[0].resize(nCells * nFuncs_);
      Array<int>& dof0 = dofs[0];
      
      for (int c=0; c<nCells; c++)
        {
          for (int i=0; i<nFuncs_; i++)
            {
              dof0[c*nFuncs_ + i] = nodeDofs_[cellLID[c]*nFuncs_+i];
            }
        }
    }
  else
    {
      int nFacets = mesh().numFacets(cellDim, cellLID[0], 0);
      nNodes[0] = nFacets;
      int dofsPerCell = nFuncs_ * nNodes[0];
      dofs[0].resize(nCells * dofsPerCell);
      Array<int>& dof0 = dofs[0];
      Array<int> facetLIDs(nCells * nNodes[0]);
      Array<int> dummy(nCells * nNodes[0]);
      mesh().getFacetLIDs(cellDim, cellLID, 0, facetLIDs, dummy);

      for (int c=0; c<nCells; c++)
        {
          for (int f=0; f<nFacets; f++)
            {
              int facetCellLID = facetLIDs[c*nFacets+f];
              for (int i=0; i<nFuncs_; i++)
                {
                  dof0[(c*nFuncs_+i)*nFacets + f] 
                    = nodeDofs_[facetCellLID*nFuncs_+i];
                }
            }
        }
    }

  return structure_;
}



void NodalDOFMap::init() 
{ 
  Tabs tab;

  SUNDANCE_MSG1(setupVerb(), tab << "initializing nodal DOF map");

  Array<Array<int> > remoteNodes(mesh().comm().getNProc());
  Array<int> hasProcessedCell(nNodes_, 0);

  /* start the DOF count at zero. */
  int nextDOF = 0;
  int owner;

  int nFacets = mesh().numFacets(dim_, 0, 0);
  Array<int> elemLID(nElems_);
  Array<int> facetLID(nFacets*nElems_);
  Array<int> facetOrientations(nFacets*nElems_);

  /* Assign node DOFs for locally owned 0-cells */
  CellSet maxCells = maxCellFilter_.getCells(mesh());

  int cc = 0;
  for (CellIterator iter=maxCells.begin(); iter != maxCells.end(); iter++, cc++)
    {
      int c = *iter;
      elemLID[cc] = c;
    }
  /* look up the LIDs of the facets */
  mesh().getFacetLIDs(dim_, elemLID, 0, facetLID, facetOrientations);
  
  for (int c=0; c<nElems_; c++)
    {
      /* for each facet, process its DOFs */
      for (int f=0; f<nFacets; f++)
        {
          /* if the facet's DOFs have been assigned already,
           * we're done */
          int fLID = facetLID[c*nFacets+f];
          if (hasProcessedCell[fLID] == 0)
            {
              /* the facet may be owned by another processor */
              if (isRemote(0, fLID, owner))
                {
                  int facetGID 
                    = mesh().mapLIDToGID(0, fLID);
                  remoteNodes[owner].append(facetGID);
                  
                }
              else /* we can assign a DOF locally */
                {
                  /* assign DOFs */
                  for (int i=0; i<nFuncs_; i++)
                    {
                      nodeDofs_[fLID*nFuncs_ + i] = nextDOF;
                      nextDOF++;
                    }
                }
              hasProcessedCell[fLID] = 1;
            }
        }
    }

  
  /* Compute offsets for each processor */
  int localCount = nextDOF;
  computeOffsets(localCount);
  
  /* Resolve remote DOF numbers */
  shareRemoteDOFs(remoteNodes);


  /* Assign DOFs for elements */
  for (int c=0; c<nElems_; c++)
    {
      /* set the element DOFs given the dofs of the facets */
      for (int f=0; f<nFacets; f++)
        {
          int fLID = facetLID[c*nFacets+f];
          for (int i=0; i<nFuncs_; i++)
            {
              elemDofs_[(c*nFuncs_+i)*nFacets + f] = nodeDofs_[fLID*nFuncs_ + i]; 
            }
        }
    }

}


void NodalDOFMap::computeOffsets(int localCount)
{
  Array<int> dofOffsets;
  int totalDOFCount;
  int myOffset = 0;
  int np = mesh().comm().getNProc();
  if (np > 1)
    {
      MPIContainerComm<int>::accumulate(localCount, dofOffsets, totalDOFCount,
                                        mesh().comm());
      myOffset = dofOffsets[mesh().comm().getRank()];

      int nDofs = nNodes_ * nFuncs_;
      for (int i=0; i<nDofs; i++)
        {
          if (nodeDofs_[i] >= 0) nodeDofs_[i] += myOffset;
        }
    }
  else
    {
      totalDOFCount = localCount;
    }
  
  setLowestLocalDOF(myOffset);
  setNumLocalDOFs(localCount);
  setTotalNumDOFs(totalDOFCount);
}


void NodalDOFMap::shareRemoteDOFs(const Array<Array<int> >& outgoingCellRequests)
{
  int np = mesh().comm().getNProc();
  if (np==1) return;
  int rank = mesh().comm().getRank();

  Array<Array<int> > incomingCellRequests;
  Array<Array<int> > outgoingDOFs(np);
  Array<Array<int> > incomingDOFs;

  SUNDANCE_MSG2(setupVerb(),  
               "p=" << mesh().comm().getRank()
               << "synchronizing DOFs for cells of dimension 0");
  SUNDANCE_MSG2(setupVerb(),  
               "p=" << mesh().comm().getRank()
               << " sending cell reqs d=0, GID=" 
               << outgoingCellRequests);

  /* share the cell requests */
  MPIContainerComm<int>::allToAll(outgoingCellRequests, 
                                  incomingCellRequests,
                                  mesh().comm());
  
  /* get DOF numbers for the zeroth function index on every node that's been 
   * requested by someone else */
  for (int p=0; p<np; p++)
    {
      if (p==rank) continue;
      const Array<int>& requestsFromProc = incomingCellRequests[p];
      int nReq = requestsFromProc.size();

      SUNDANCE_MSG4(setupVerb(), "p=" << mesh().comm().getRank() 
                            << " recv'd from proc=" << p
                            << " reqs for DOFs for cells " 
                            << requestsFromProc);

      outgoingDOFs[p].resize(nReq);

      for (int c=0; c<nReq; c++)
        {
          int GID = requestsFromProc[c];
          SUNDANCE_MSG3(setupVerb(),  
                       "p=" << rank
                       << " processing zero-cell with GID=" << GID); 
          int LID = mesh().mapGIDToLID(0, GID);
          SUNDANCE_MSG3(setupVerb(),  
                       "p=" << rank
                       << " LID=" << LID << " dofs=" 
                       << nodeDofs_[LID*nFuncs_]);
          outgoingDOFs[p][c] = nodeDofs_[LID*nFuncs_];
          SUNDANCE_MSG3(setupVerb(),  
                       "p=" << rank
                       << " done processing cell with GID=" << GID);
        }
    }

  SUNDANCE_MSG2(setupVerb(),  
               "p=" << mesh().comm().getRank()
               << "answering DOF requests for cells of dimension 0");

  /* share the DOF numbers */
  MPIContainerComm<int>::allToAll(outgoingDOFs,
                                  incomingDOFs,
                                  mesh().comm());

  SUNDANCE_MSG2(setupVerb(),  
               "p=" << mesh().comm().getRank()
               << "communicated DOF answers for cells of dimension 0" );

  
  /* now assign the DOFs from the other procs */

  for (int p=0; p<mesh().comm().getNProc(); p++)
    {
      if (p==mesh().comm().getRank()) continue;
      const Array<int>& dofsFromProc = incomingDOFs[p];
      int numCells = dofsFromProc.size();
      for (int c=0; c<numCells; c++)
        {
          int cellGID = outgoingCellRequests[p][c];
          int cellLID = mesh().mapGIDToLID(0, cellGID);
          int dof = dofsFromProc[c];
          for (int i=0; i<nFuncs_; i++)
            {
              nodeDofs_[cellLID*nFuncs_ + i] = dof+i;
              addGhostIndex(dof+i);
            }
        }
    }
}
