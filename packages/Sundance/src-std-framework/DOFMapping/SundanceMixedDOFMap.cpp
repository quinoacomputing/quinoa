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

#include "SundanceMixedDOFMap.hpp"
#include "SundanceBasisDOFTopologyBase.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "PlayaMPIContainerComm.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"


using namespace Sundance;
using namespace Teuchos;
using Playa::MPIComm;
using Playa::MPIContainerComm;
using Playa::MPIDataType;
using Playa::MPIOp;


static Time& mixedDOFCtorTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("mixed DOF map init"); 
  return *rtn;
}
static Time& maxDOFBuildTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("max-cell dof table init"); 
  return *rtn;
}

MixedDOFMap::MixedDOFMap(const Mesh& mesh, 
  const Array<RCP<BasisDOFTopologyBase> >& basis,
  const CellFilter& maxCells,
  int setupVerb)
  : SpatiallyHomogeneousDOFMapBase(mesh, basis.size(), setupVerb), 
    maxCells_(maxCells),
    dim_(mesh.spatialDim()),
    dofs_(mesh.spatialDim()+1),
    maximalDofs_(),
    haveMaximalDofs_(false),
    localNodePtrs_(),
    nNodesPerCell_(),
    totalNNodesPerCell_(),
    cellHasAnyDOFs_(dim_+1),
    numFacets_(mesh.spatialDim()+1),
    originalFacetOrientation_(2),
    hasBeenAssigned_(mesh.spatialDim()+1),
    structure_(),
    nFuncs_()
{
  TimeMonitor timer(mixedDOFCtorTimer());
  Tabs tab;
  SUNDANCE_MSG1(setupVerb, tab << "building mixed DOF map");

  Sundance::Map<OrderedHandle<BasisDOFTopologyBase>, int> basisToChunkMap;
  Array<RCP<BasisDOFTopologyBase> > chunkBases;
  Array<Array<int> > chunkFuncs;
  
  int chunk = 0;
  int nBasis = basis.size();
  for (int i=0; i<nBasis; i++)
  {
    OrderedHandle<BasisDOFTopologyBase> bh = basis[i];
    if (!basisToChunkMap.containsKey(bh))
    {
      chunkBases.append(basis[i]);
      basisToChunkMap.put(bh, chunk);
      chunkFuncs.append(tuple(i));
      chunk++;
    }
    else
    {
      int b = basisToChunkMap.get(bh);
      chunkFuncs[b].append(i);
    }
  }

  Tabs tab1;
  SUNDANCE_MSG2(setupVerb, tab1 << "basis to chunk map = " << basisToChunkMap);


  structure_ = rcp(new MapStructure(basis.size(), chunkBases, chunkFuncs));

  nFuncs_.resize(chunkBases.size());
  for (int i=0; i<nFuncs_.size(); i++) 
    nFuncs_[i] = chunkFuncs[i].size();

  allocate(mesh);

  initMap();

  buildMaximalDofTable();

  /* do a sanity check */
  checkTable();

  SUNDANCE_MSG1(setupVerb, tab << "done building mixed DOF map");
}


void MixedDOFMap::allocate(const Mesh& mesh)
{
  Tabs tab;

  /* gather functions into chunks sharing identical basis functions */
  SUNDANCE_MSG1(setupVerb(), tab << "grouping like basis functions");
  

  /* now that we know the number of basis chunks, allocate arrays */
  localNodePtrs_.resize(nBasisChunks());
  nNodesPerCell_.resize(nBasisChunks());
  totalNNodesPerCell_.resize(nBasisChunks());
  nDofsPerCell_.resize(nBasisChunks());
  totalNDofsPerCell_.resize(nBasisChunks());
  maximalDofs_.resize(nBasisChunks());

  for (int b=0; b<nBasisChunks(); b++)
  {
    localNodePtrs_[b].resize(mesh.spatialDim()+1);
    nNodesPerCell_[b].resize(mesh.spatialDim()+1);
    totalNNodesPerCell_[b].resize(mesh.spatialDim()+1);
    nDofsPerCell_[b].resize(mesh.spatialDim()+1);
    totalNDofsPerCell_[b].resize(mesh.spatialDim()+1);
  }
  

  /* compute node counts for each cell dimension and each basis type */
  SUNDANCE_MSG1(setupVerb(), tab << "working out DOF map node counts");
  
  numFacets_.resize(dim_+1);

  for (int d=0; d<=dim_; d++)
  {
    Tabs tab1;
    SUNDANCE_MSG2(setupVerb(), tab << "allocating maps for cell dim=" << d);
    /* record the number of facets for each cell type so we're
     * not making a bunch of mesh calls */
    numFacets_[d].resize(d);
    for (int fd=0; fd<d; fd++) numFacets_[d][fd]=mesh.numFacets(d, 0, fd);
    SUNDANCE_MSG3(setupVerb(), tab1 << "num facets for dimension " << d << " is " 
      << numFacets_[d]);

    cellHasAnyDOFs_[d] = false;
    dofs_[d].resize(nBasisChunks());

    int numCells = mesh.numCells(d);
    hasBeenAssigned_[d].resize(numCells);

    for (int b=0; b<nBasisChunks(); b++)
    {
      Tabs tab2;
      SUNDANCE_MSG2(setupVerb(), tab1 << "basis chunk=" << b);
      SUNDANCE_MSG2(setupVerb(), tab2 << "basis=" << basis(b)->description());
      int nNodes = basis(b).ptr()->nReferenceDOFsWithFacets(mesh.cellType(dim_), mesh.cellType(d));
      SUNDANCE_MSG2(setupVerb(), tab2 << "nNodes per cell=" << nNodes);
      if (nNodes == 0)
      {
        nNodesPerCell_[b][d] = 0;
        nDofsPerCell_[b][d] = 0;
      }
      else
      {
        /* look up the node pointer for this cell and for all of its
         * facets */
        basis(b).ptr()->getReferenceDOFs(mesh.cellType(dim_),
          mesh.cellType(d), 
          localNodePtrs_[b][d]);
              
              
        SUNDANCE_MSG3(setupVerb(), tab2 << "node ptrs are " 
          << localNodePtrs_[b][d]);
              
        /* with the node pointers in hand, we can work out the number
         * of nodes per cell in this dimension for this basis */
        if (localNodePtrs_[b][d][d].size() > 0) 
        {
          nNodesPerCell_[b][d] = localNodePtrs_[b][d][d][0].size();
          if (nNodesPerCell_[b][d] > 0) cellHasAnyDOFs_[d] = true;
        }
        else
        {
          nNodesPerCell_[b][d] = 0;
        }
        nDofsPerCell_[b][d] = nNodesPerCell_[b][d] * nFuncs(b);
      }

      /* if the cell is of intermediate dimension and it has DOFs, we
       * need to assign GIDs for the cells of this dimension. Otherwise,
       * we can skip intermediate GID assignment, saving some parallel
       * communication */
      if (nDofsPerCell_[b][d] > 0 && d > 0 && d < mesh.spatialDim())
      {
        Mesh& tmpMesh = const_cast<Mesh&>(mesh);
        tmpMesh.assignIntermediateCellGIDs(d);
      }
          
      SUNDANCE_MSG3(setupVerb(), tab2 << 
        "num nodes is " 
        << nNodesPerCell_[b][d]);

      totalNNodesPerCell_[b][d] = nNodesPerCell_[b][d];
      for (int dd=0; dd<d; dd++) 
      {
        totalNNodesPerCell_[b][d] 
          += numFacets_[d][dd]*nNodesPerCell_[b][dd];
      }
      totalNDofsPerCell_[b][d] = totalNNodesPerCell_[b][d] * nFuncs(b);
      SUNDANCE_MSG3(setupVerb(), tab2 << "totalNDofsPerCell " << totalNDofsPerCell_[b][d]);

      if (nNodes > 0)
      {
        /* allocate the DOFs array for this dimension */
        dofs_[d][b].resize(nDofsPerCell_[b][d] * numCells);
              
        /* set all DOFs to a marker value */
        int nDof = dofs_[d][b].size();
        Array<int>& dofs = dofs_[d][b];
        int marker = uninitializedVal();
        for (int i=0; i<nDof; i++) 
        {
          dofs[i] = marker;
        }
      }
      /* allocate the maximal dof array */
      if (d==dim_)
      {
        maximalDofs_[b].resize(totalNDofsPerCell_[b][d]*numCells);
      }
    }

    /* allocate the array of original facet orientations */
    if (d > 0 && d < dim_) 
    {
      originalFacetOrientation_[d-1].resize(numCells);
    }
      
  }
  SUNDANCE_MSG1(setupVerb(), tab << "done allocating DOF map");
}

void MixedDOFMap::initMap()
{
  Tabs tab;
  SUNDANCE_MSG1(setupVerb(), tab << "initializing DOF map");
  /* start the DOF count at zero. */
  int nextDOF = 0;

  /* Space in which to keep a list of remote cells needed by each processor
   * for each dimension. The first index is dimension, the second proc, the
   * third cell number. */
  Array<Array<Array<int> > > remoteCells(mesh().spatialDim()+1);

  for (int d=0; d<remoteCells.size(); d++) 
    remoteCells[d].resize(mesh().comm().getNProc());
  
  /* Loop over maximal cells in the order specified by the cell iterator.
   * This might be reordered relative to the mesh. 
   *
   * At each maximal cell, we'll run through the facets and 
   * assign DOFs. That will take somewhat more work, but gives much 
   * better cache locality for the matrix because all the DOFs for
   * each maximal element and its facets are grouped together. */
  
  CellSet cells = maxCells_.getCells(mesh());
  CellIterator iter;
  for (iter=cells.begin(); iter != cells.end(); iter++)
  {
    /* first assign any DOFs associated with the maximal cell */
    int cellLID = *iter;
    int owner;
      
    if (cellHasAnyDOFs_[dim_])
    {
      /* if the maximal cell is owned by another processor,
       * put it in the list of cells for which we need to request 
       * DOF information from another processor */
      if (isRemote(dim_, cellLID, owner))
      {
        int cellGID = mesh().mapLIDToGID(dim_, cellLID);
        SUNDANCE_MSG4(setupVerb(), "proc=" << comm().getRank() 
          << " thinks d-" << dim_ 
          << " cell GID=" << cellGID
          << " is owned remotely by p=" 
          << owner);
        remoteCells[dim_][owner].append(cellGID); 
      }
      else /* the cell is locally owned, so we can 
            * set its DOF numbers now */
      {
        for (int b=0; b<nBasisChunks(); b++)
        {
          setDOFs(b, dim_, cellLID, nextDOF);
        }
      }
    }

    /* Now assign any DOFs associated with the facets. */
    for (int d=0; d<dim_; d++)
    {
      if (cellHasAnyDOFs_[d])
      {
        int nf = numFacets_[dim_][d];
        Array<int> facetLID(nf);
        Array<int> facetOrientations(nf);
        /* look up the LIDs of the facets */
        mesh().getFacetArray(dim_, cellLID, d, 
          facetLID, facetOrientations);
        /* for each facet, process its DOFs */
        for (int f=0; f<nf; f++)
        {
          /* if the facet's DOFs have been assigned already,
           * we're done */
          if (!hasBeenAssigned(d, facetLID[f]))
          {
            markAsAssigned(d, facetLID[f]);
            /* the facet may be owned by another processor */
            if (isRemote(d, facetLID[f], owner))
            {
              int facetGID 
                = mesh().mapLIDToGID(d, facetLID[f]);
              SUNDANCE_MSG4(setupVerb(), "proc=" << comm().getRank() 
                << " thinks d-" << d 
                << " cell GID=" << facetGID
                << " is owned remotely by p=" << owner);
              remoteCells[d][owner].append(facetGID);
            }
            else /* we can assign a DOF locally */
            {
              /* assign DOF */
              for (int b=0; b<nBasisChunks(); b++)
              {
                setDOFs(b, d, facetLID[f], nextDOF);
              }
              /* record the orientation wrt the maximal cell */
              if (d > 0) 
              {
                originalFacetOrientation_[d-1][facetLID[f]] 
                  = facetOrientations[f];
              }
            }
          }
        }
      }
    }
  }
    

  /* Done with first pass, in which we have assigned DOFs for all
   * local processors. We now have to share DOF information between
   * processors */

  int numLocalDOFs = nextDOF;
  if (mesh().comm().getNProc() > 1)
  {
    for (int d=0; d<=dim_; d++)
    {
      if (cellHasAnyDOFs_[d])
      {
        computeOffsets(d, numLocalDOFs);
        shareDOFs(d, remoteCells[d]);
      }
    }
  }
  else
  {
    setLowestLocalDOF(0);
    setNumLocalDOFs(numLocalDOFs);
    setTotalNumDOFs(numLocalDOFs);
  }
  SUNDANCE_MSG1(setupVerb(), tab << "done initializing DOF map");
}

void MixedDOFMap::shareDOFs(int cellDim,
  const Array<Array<int> >& outgoingCellRequests)
{
  int np = mesh().comm().getNProc();
  int rank = mesh().comm().getRank();

  Array<Array<int> > incomingCellRequests;
  Array<Array<int> > outgoingDOFs(np);
  Array<Array<int> > incomingDOFs;

  SUNDANCE_MSG2(setupVerb(),  
    "p=" << mesh().comm().getRank()
    << "synchronizing DOFs for cells of dimension " << cellDim);
  SUNDANCE_MSG2(setupVerb(),  
    "p=" << mesh().comm().getRank()
    << " sending cell reqs d=" << cellDim << " GID=" << outgoingCellRequests);

  /* share the cell requests */
  MPIContainerComm<int>::allToAll(outgoingCellRequests, 
    incomingCellRequests,
    mesh().comm());
  
  /* we send the following information in response:
   * (1) The first DOF for each chunk for the requested cell
   * (2) The orientation if the cell is an edge or face 
   */
  int blockSize = 0;
  bool sendOrientation = false;
  for (int b=0; b<nBasisChunks(); b++)
  {
    int nDofs = nDofsPerCell_[b][cellDim];
    if (nDofs > 0) blockSize++;
    if (nDofs > 1 && cellDim > 0 && cellDim < dim_) sendOrientation = true;
  }
  blockSize += sendOrientation;

  SUNDANCE_MSG2(setupVerb(),  
    "p=" << rank
    << "recvd DOF requests for cells of dimension " << cellDim
    << " GID=" << incomingCellRequests);

  /* get orientations and DOF numbers for the first node of every cell that's been 
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

    outgoingDOFs[p].resize(nReq * blockSize);

    for (int c=0; c<nReq; c++)
    {
      int GID = requestsFromProc[c];
      SUNDANCE_MSG3(setupVerb(),  
        "p=" << rank
        << " processing cell with d=" << cellDim 
        << " GID=" << GID);
      int LID = mesh().mapGIDToLID(cellDim, GID);
      SUNDANCE_MSG3(setupVerb(),  
        "p=" << rank
        << " LID=" << LID << " dofs=" << dofs_[cellDim]);
      int blockOffset = 0;
      for (int b=0; b<nBasisChunks(); b++)
      {
        if (nDofsPerCell_[b][cellDim] == 0) continue;
        outgoingDOFs[p][blockSize*c+blockOffset] 
          = getInitialDOFForCell(cellDim, LID, b);
        blockOffset++;
      }
      if (sendOrientation)
      {
        outgoingDOFs[p][blockSize*(c+1) - 1] 
          = originalFacetOrientation_[cellDim-1][LID];
      }
      SUNDANCE_MSG3(setupVerb(),  
        "p=" << rank
        << " done processing cell with GID=" << GID);
    }
  }
 

  SUNDANCE_MSG2(setupVerb(),  
    "p=" << mesh().comm().getRank()
    << "answering DOF requests for cells of dimension " << cellDim);

  /* share the DOF numbers */
  MPIContainerComm<int>::allToAll(outgoingDOFs,
    incomingDOFs,
    mesh().comm());

  SUNDANCE_MSG2(setupVerb(),  
    "p=" << mesh().comm().getRank()
    << "communicated DOF answers for cells of dimension " << cellDim);

  
  /* now assign the DOFs from the other procs */

  for (int p=0; p<mesh().comm().getNProc(); p++)
  {
    if (p==mesh().comm().getRank()) continue;
    const Array<int>& dofsFromProc = incomingDOFs[p];
    int numCells = dofsFromProc.size()/blockSize;
    for (int c=0; c<numCells; c++)
    {
      int cellGID = outgoingCellRequests[p][c];
      int cellLID = mesh().mapGIDToLID(cellDim, cellGID);
      int blockOffset = 0;
      for (int b=0; b<nBasisChunks(); b++)
      {
        if (nDofsPerCell_[b][cellDim] == 0) continue;
        int dof = dofsFromProc[blockSize*c+blockOffset];
        setDOFs(b, cellDim, cellLID, dof, true);
        blockOffset++;
      }
      if (sendOrientation) 
      {
        originalFacetOrientation_[cellDim-1][cellLID] 
          = dofsFromProc[blockSize*(c+1)-1];
      }
    }
  }
  
}



void MixedDOFMap::setDOFs(int basisChunk, int cellDim, int cellLID, 
  int& nextDOF, bool isRemote)
{
  Tabs tab;
  SUNDANCE_MSG3(setupVerb(), tab << "setting DOFs for " << cellDim 
    << "-cell " << cellLID);
  int nDofs = nDofsPerCell_[basisChunk][cellDim];
  if (nDofs==0) return;

  int* ptr = getInitialDOFPtrForCell(cellDim, cellLID, basisChunk);
  
  if (isRemote)
  {
    for (int i=0; i<nDofs; i++, nextDOF++) 
    {
      ptr[i] = nextDOF;;
      addGhostIndex(nextDOF);
    }
  }
  else
  {
    for (int i=0; i<nDofs; i++,nextDOF++) 
    {
      ptr[i] = nextDOF;
    }
  }
}



RCP<const MapStructure> MixedDOFMap
::getDOFsForCellBatch(int cellDim,
  const Array<int>& cellLID,
  const Set<int>& requestedFuncSet,
  Array<Array<int> >& dofs,
  Array<int>& nNodes,
  int verbosity) const 
{
  TimeMonitor timer(batchedDofLookupTimer());

  Tabs tab;
  SUNDANCE_MSG3(verbosity, 
    tab << "getDOFsForCellBatch(): cellDim=" << cellDim
    << " cellLID=" << cellLID);

  dofs.resize(nBasisChunks());
  nNodes.resize(nBasisChunks());

  int nCells = cellLID.size();

  if (cellDim == dim_)
  {
    Tabs tab1;

    if (!haveMaximalDofs_) 
    {
      buildMaximalDofTable();
    }

    SUNDANCE_MSG4(verbosity, tab1 << "getting dofs for maximal cells");

    for (int b=0; b<nBasisChunks(); b++)
    {
      nNodes[b] = totalNNodesPerCell_[b][cellDim];
      dofs[b].resize(nNodes[b]*nFuncs(b)*nCells);
      int dofsPerCell = nFuncs(b)*nNodes[b];
      Array<int>& chunkDofs = dofs[b];
      for (int c=0; c<nCells; c++)
      {
        for (int i=0; i<dofsPerCell; i++)
        {
          chunkDofs[c*dofsPerCell + i] 
            = maximalDofs_[b][cellLID[c]*dofsPerCell+i];
        }
      }
    }
  }
  else
  {
    Tabs tab1;
    SUNDANCE_MSG4(verbosity, tab1 << "getting dofs for non-maximal cells");
  
    static Array<Array<int> > facetLID(3);
    static Array<Array<int> > facetOrientations(3);
    static Array<int> numFacets(3);

    for (int d=0; d<cellDim; d++) 
    {
      numFacets[d] = mesh().numFacets(cellDim, cellLID[0], d);
      mesh().getFacetLIDs(cellDim, cellLID, d, facetLID[d], 
        facetOrientations[d]);
    }

    for (int b=0; b<nBasisChunks(); b++)
    {
      nNodes[b] = totalNNodesPerCell_[b][cellDim];
      dofs[b].resize(nNodes[b]*nFuncs(b)*nCells);
      int dofsPerCell = nFuncs(b)*nNodes[b];
          
      Array<int>& toPtr = dofs[b];
      int nf = nFuncs(b);

      for (int c=0; c<nCells; c++)
      {
        Tabs tab2;
        SUNDANCE_MSG4(verbosity, tab2 << "cell=" << c);
        int offset = dofsPerCell*c;

        /* first get the DOFs for the nodes associated with 
         * the cell's interior */
        SUNDANCE_MSG4(verbosity, tab2 << "doing interior nodes");
        int nInteriorNodes = nNodesPerCell_[b][cellDim];
        //              int nInteriorNodes = localNodePtrs_[b][cellDim][cellDim][0].size();
        if (nInteriorNodes > 0)
        {
          if (cellDim==0 || nInteriorNodes <= 1) /* orientation-independent */
          {
            // const int* fromPtr 
            //                        = getInitialDOFPtrForCell(cellDim, cellLID[c], b);

            for (int func=0; func<nf; func++)
            {
              for (int n=0; n<nInteriorNodes; n++)
              {
                int ptr = localNodePtrs_[b][cellDim][cellDim][0][n];
                toPtr[offset + func*nNodes[b] + ptr] 
                  = dofs_[cellDim][b][cellLID[c]*nDofsPerCell_[b][cellDim]+func*nNodesPerCell_[b][cellDim]+n];
                //                                = fromPtr[func*nNodes[b] + n];
              }
            }
          }
          else
          {
            int sign = originalFacetOrientation_[cellDim-1][cellLID[c]];
            int nInteriorNodes = localNodePtrs_[b][cellDim][cellDim][0].size();
            //     const int* fromPtr 
            //                        = getInitialDOFPtrForCell(cellDim, cellLID[c], b);
                 
            for (int func=0; func<nf; func++)
            {
              for (int m=0; m<nInteriorNodes; m++)
              {
                int n = m;
                if (sign<0) n = nInteriorNodes-1-m;
                int ptr = localNodePtrs_[b][cellDim][cellDim][0][m];
                toPtr[offset + func*nNodes[b] + ptr]  = dofs_[cellDim][b][cellLID[c]*nDofsPerCell_[b][cellDim]+func*nNodesPerCell_[b][cellDim]+n];
                //    = fromPtr[func*nNodes[b] + n];
              }
            }
          }
        }

        /* now do the facets */
        for (int d=0; d<cellDim; d++)
        {
          Tabs tab2;
          SUNDANCE_MSG4(verbosity, tab2 << "facet dim=" << d);
          if (nNodesPerCell_[b][d] == 0) continue;
          for (int f=0; f<numFacets[d]; f++)
          {
            Tabs tab3;
            int facetID = facetLID[d][c*numFacets[d]+f];
            SUNDANCE_MSG4(verbosity, tab2 << "f=" << f << " facetLID=" << facetID);
            if (localNodePtrs_[b][cellDim].size()==0) continue;
            int nFacetNodes = localNodePtrs_[b][cellDim][d][f].size();
            //const int* fromPtr = getInitialDOFPtrForCell(d, facetID, b);
            int* toPtr1 = &(dofs[b][dofsPerCell*c]);
            const int* nodePtr = &(localNodePtrs_[b][cellDim][d][f][0]);
            for (int func=0; func<nf; func++)
            {
              if (d == 0 || nFacetNodes <= 1) /* orientation-independent */
              {
                for (int n=0; n<nFacetNodes; n++)
                {
                  int ptr = nodePtr[n];
                  toPtr1[func*nNodes[b] + ptr] //= fromPtr[func*nNodes[b] + n];
                    = dofs_[d][b][facetID*nDofsPerCell_[b][d]+func*nNodesPerCell_[b][d]+n];
                }
              }
              else /* orientation-dependent */
              {
                int facetOrientation = facetOrientations[d][c*numFacets[d]+f]
                  * originalFacetOrientation_[d-1][facetID];
                for (int m=0; m<nFacetNodes; m++)
                {
                  int n = m;
                  if (facetOrientation<0) n = nFacetNodes-1-m;
                  int ptr = nodePtr[n];
                  toPtr1[func*nNodes[b]+ptr] 
                    = dofs_[d][b][facetID*nDofsPerCell_[b][d]+func*nNodesPerCell_[b][d]+n];
                  //                                    = fromPtr[func*nNodes[b]+n];
                }
              }
            }
          }
        }
      }
    }
  }
  return structure_;
}    

void MixedDOFMap::buildMaximalDofTable() const
{
  TimeMonitor timer(maxDOFBuildTimer());
  Tabs tab;
  int cellDim = dim_;
  int nCells = mesh().numCells(dim_);

  SUNDANCE_MSG2(setupVerb(), tab << "building dofs for maximal cells");

  Array<Array<int> > facetLID(3);
  Array<Array<int> > facetOrientations(3);
  Array<int> numFacets(3);

  Array<int> cellLID(nCells);

  for (int c=0; c<nCells; c++) cellLID[c]=c;
  
  for (int d=0; d<cellDim; d++) 
  {
    numFacets[d] = mesh().numFacets(cellDim, cellLID[0], d);
    mesh().getFacetLIDs(cellDim, cellLID, d, 
      facetLID[d], facetOrientations[d]);
  }

  Array<int> nInteriorNodes(nBasisChunks());
  Array<int> nNodes(nBasisChunks());
  for (int b = 0; b<nBasisChunks(); b++)
  {
    if (localNodePtrs_[b][cellDim].size() != 0)
    {
      nInteriorNodes[b] = localNodePtrs_[b][cellDim][cellDim][0].size();
    }
    else 
    {
      nInteriorNodes[b] = 0;
    }
    nNodes[b] = totalNNodesPerCell_[b][cellDim];
  }

  for (int c=0; c<nCells; c++)
  {
    Tabs tab1;
    SUNDANCE_MSG4(setupVerb(), tab1 << "working on cell=" << c 
      << " LID=" << cellLID[c]);
    /* first get the DOFs for the nodes associated with 
     * the cell's interior */
    SUNDANCE_MSG4(setupVerb(), tab1 << "doing interior nodes");
    for (int b=0; b<nBasisChunks(); b++)
    {
      SUNDANCE_MSG4(setupVerb(), tab1 << "basis chunk=" << b); 
      if (nInteriorNodes[b]>0)
      {
        SUNDANCE_MSG4(setupVerb(), tab1<< "nInteriorNodes = " <<nInteriorNodes[b]);
        //const int* fromPtr = getInitialDOFPtrForCell(dim_, cellLID[c], b);
        int* toPtr = &(maximalDofs_[b][nNodes[b]*nFuncs(b)*cellLID[c]]);
        int nf = nFuncs(b);
        for (int func=0; func<nf; func++)
        {
          for (int n=0; n<nInteriorNodes[b]; n++)
          {
                      
            int ptr = localNodePtrs_[b][cellDim][cellDim][0][n];
            toPtr[func*nNodes[b] + ptr] = //fromPtr[func*nNodes[b] + n];
              dofs_[cellDim][b][cellLID[c]*nDofsPerCell_[b][cellDim]+func*nNodesPerCell_[b][cellDim]+n];
          }
        }
      }
    }
      
    SUNDANCE_MSG4(setupVerb(), tab1 << "doing facet nodes");
    /* now get the DOFs for the nodes on the facets */
    for (int d=0; d<cellDim; d++)
    {
      Tabs tab2;
      SUNDANCE_MSG4(setupVerb(), tab2 << "facet dim=" << d);

      for (int f=0; f<numFacets[d]; f++)
      {
        Tabs tab3;
        int facetID = facetLID[d][c*numFacets[d]+f];
        SUNDANCE_MSG4(setupVerb(), tab2 << "f=" << f << " facetLID=" << facetID);

        for (int b=0; b<nBasisChunks(); b++)
        {
          Tabs tab4;
          SUNDANCE_MSG4(setupVerb(), tab4 << "basis chunk=" << b); 
          SUNDANCE_MSG4(setupVerb(), tab4 << "num nodes per cell=" << nNodesPerCell_[b][d]); 
          int nf = nFuncs(b);
          if (nDofsPerCell_[b][d]==0) continue;
          int nFacetNodes = 0;
          if (localNodePtrs_[b][cellDim].size()!=0)
            nFacetNodes = localNodePtrs_[b][cellDim][d][f].size();
          if (nFacetNodes == 0) continue;
          //  const int* fromPtr = getInitialDOFPtrForCell(d, facetID, b);
          SUNDANCE_MSG4(setupVerb(), tab4 << "dof table size=" << maximalDofs_[b].size()); 
          int* toPtr = &(maximalDofs_[b][nNodes[b]*nFuncs(b)*cellLID[c]]);
          const int* nodePtr = &(localNodePtrs_[b][cellDim][d][f][0]);
          for (int func=0; func<nf; func++)
          {
            Tabs tab5;
            SUNDANCE_MSG4(setupVerb(), tab5 << "func=" << func); 
            if (d == 0 || nFacetNodes <= 1) /* orientation-independent */
            {
              Tabs tab6;
              for (int n=0; n<nFacetNodes; n++)
              {
                SUNDANCE_MSG4(setupVerb(), tab5 << "n=" << n); 
                int ptr = nodePtr[n];
                SUNDANCE_MSG4(setupVerb(), tab6 << "ptr=" << ptr); 

                toPtr[func*nNodes[b] + ptr] 
                  = dofs_[d][b][facetID*nDofsPerCell_[b][d]+func*nNodesPerCell_[b][d]+n];
                //= fromPtr[func*nNodes[b] + n];
              }
            }
            else /* orientation-dependent */
            {
              int facetOrientation = facetOrientations[d][c*numFacets[d]+f]
                * originalFacetOrientation_[d-1][facetID];
              for (int m=0; m<nFacetNodes; m++)
              {
                int n = m;
                if (facetOrientation<0) n = nFacetNodes-1-m;
                int ptr = nodePtr[m];
                toPtr[func*nNodes[b]+ptr] 
                  = dofs_[d][b][facetID*nDofsPerCell_[b][d]+func*nNodesPerCell_[b][d]+n];
                //= fromPtr[func*nNodes[b]+n];
              }
            }
          }
        }
      }
    }
  }

  haveMaximalDofs_ = true;

  SUNDANCE_MSG2(setupVerb(), tab << "done building dofs for maximal cells");
}





void MixedDOFMap::computeOffsets(int dim, int localCount)
{
  if (setupVerb() > 2)
  {
    comm().synchronize();
    comm().synchronize();
    comm().synchronize();
    comm().synchronize();
  }
  SUNDANCE_MSG2(setupVerb(), 
    "p=" << mesh().comm().getRank()
    << " sharing offsets for DOF numbering for dim=" << dim);

  SUNDANCE_MSG2(setupVerb(), 
    "p=" << mesh().comm().getRank()
    << " I have " << localCount << " cells");

  Array<int> dofOffsets;
  int totalDOFCount;
  MPIContainerComm<int>::accumulate(localCount, dofOffsets, totalDOFCount,
    mesh().comm());
  int myOffset = dofOffsets[mesh().comm().getRank()];

  SUNDANCE_MSG2(setupVerb(), 
    "p=" << mesh().comm().getRank()
    << " back from MPI accumulate");

  if (setupVerb() > 2)
  {
    comm().synchronize();
    comm().synchronize();
    comm().synchronize();
    comm().synchronize();
  }

  for (int chunk=0; chunk<dofs_[dim].size(); chunk++)
  {
    for (int n=0; n<dofs_[dim][chunk].size(); n++)
    {
      if (dofs_[dim][chunk][n] >= 0) dofs_[dim][chunk][n] += myOffset;
    }
  }

  setLowestLocalDOF(myOffset);
  setNumLocalDOFs(localCount);
  setTotalNumDOFs(totalDOFCount);

  SUNDANCE_MSG2(setupVerb(), 
    "p=" << mesh().comm().getRank() 
    << " done sharing offsets for DOF numbering for dim=" << dim);
  if (setupVerb() > 2)
  {
    comm().synchronize();
    comm().synchronize();
    comm().synchronize();
    comm().synchronize();
  }

}                           



void MixedDOFMap::checkTable() const 
{
  int bad = 0;
  for (int d=0; d<dofs_.size(); d++)
  {
    for (int chunk=0; chunk<dofs_[d].size(); chunk++)
    {
      const Array<int>& dofs = dofs_[d][chunk];
      for (int n=0; n<dofs.size(); n++)
      {
        if (dofs[n] < 0) bad = 1;
      }
    }
  }
  
  int anyBad = bad;
  comm().allReduce((void*) &bad, (void*) &anyBad, 1, 
    MPIDataType::intType(), MPIOp::sumOp());
  TEUCHOS_TEST_FOR_EXCEPTION(anyBad > 0, std::runtime_error, "invalid DOF map");
}
