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

#include "SundanceMixedDOFMapHN.hpp"
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

static Time& mixedHNDOFCtorTimer()
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("mixed hanging-node DOF map init");
  return *rtn;
}
static Time& maxDOFBuildTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("max-cell dof table init"); 
  return *rtn;
}

MixedDOFMapHN::MixedDOFMapHN(const Mesh& mesh,
  const Array<RCP<BasisDOFTopologyBase> >& basis,
  const CellFilter& maxCells,
  int setupVerb)
 : HNDoFMapBaseHomogeneous(mesh, basis.size(), setupVerb),
    maxCells_(maxCells),
    dim_(mesh.spatialDim()),
    dofs_(mesh.spatialDim()+1),
    maximalDofs_(),
    haveMaximalDofs_(false),
    localNodePtrs_(),
    hasCellHanging_(),
    isElementHanging_(mesh.spatialDim()+1),
    HN_To_globalFacetsLID_(),
    HN_To_globalFacetsDim_(),
    HN_To_coeffs_(),
    maxCellLIDwithHN_to_TrafoMatrix_(),
    matrixStore_(),
    nNodesPerCell_(),
    totalNNodesPerCell_(),
    cellHasAnyDOFs_(dim_+1),
    numFacets_(mesh.spatialDim()+1),
    originalFacetOrientation_(2),
    hasBeenAssigned_(mesh.spatialDim()+1),
    structure_(),
    nFuncs_()
{
  TimeMonitor timer(mixedHNDOFCtorTimer());
  Tabs tab;

  SUNDANCE_MSG1(setupVerb, tab << "building mixed hanging node DOF map");

  Sundance::Map<OrderedHandle<BasisDOFTopologyBase>, int> basisToChunkMap;
  Array<RCP<BasisDOFTopologyBase> > chunkBases;
  Array<Array<int> > chunkFuncs;
  
  int chunk = 0;
  int nBasis = basis.size();
  basis_.resize(nBasis);

  nPoints_ = mesh.numCells(0);

  for (int i=0; i<nBasis; i++)
  {
    OrderedHandle<BasisDOFTopologyBase> bh = basis[i];
    if (!basisToChunkMap.containsKey(bh))
    {
      chunkBases.append(basis[i]);
      basisToChunkMap.put(bh, chunk);
      chunkFuncs.append(tuple(i));
      basis_[chunk] = rcp_dynamic_cast<BasisFamilyBase>(basis[i]);
      chunk++;
    }
    else
    {
      int b = basisToChunkMap.get(bh);
      chunkFuncs[b].append(i);
    }
  }
  basis_.resize(chunk);

  // initialize the matrix store
  matrixStore_.init(chunk);

  Tabs tab1;
  SUNDANCE_MSG1(setupVerb, tab1 << "basis to chunk map = " << basisToChunkMap);


  structure_ = rcp(new MapStructure(basis.size(), chunkBases, chunkFuncs));

  nFuncs_.resize(chunkBases.size());
  for (int i=0; i<nFuncs_.size(); i++) 
    nFuncs_[i] = chunkFuncs[i].size();

  allocate(mesh);

  initMap();

  buildMaximalDofTable();

  /* do a sanity check */
  //checkTable(); // todo: not implemented yet (would be useful in parallel case)

  SUNDANCE_MSG1(setupVerb, tab << "done building mixed HN DOF map");

}


void MixedDOFMapHN::allocate(const Mesh& mesh)
{
  Tabs tab;

  /* gather functions into chunks sharing identical basis functions */
  SUNDANCE_MSG1(setupVerb(),tab << "grouping like basis functions");

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
  SUNDANCE_MSG1(setupVerb(),tab << "working out DOF map node counts");
  
  numFacets_.resize(dim_+1);
  isElementHanging_.resize(dim_+1);
  hasCellHanging_.resize(mesh.numCells(dim_),false);

  for (int d=0; d<=dim_; d++)
  {
    Tabs tab1;
    SUNDANCE_MSG1(setupVerb(),tab << "allocating maps for cell dim=" << d);
    /* record the number of facets for each cell type so we're
     * not making a bunch of mesh calls */
    numFacets_[d].resize(d);

    for (int fd=0; fd<d; fd++) numFacets_[d][fd]=mesh.numFacets(d, 0, fd);
    SUNDANCE_MSG1(setupVerb(),tab1 << "num facets for dimension " << d << " is "
      << numFacets_[d]);

    cellHasAnyDOFs_[d] = false;
    dofs_[d].resize(nBasisChunks());

    int numCells = mesh.numCells(d);
    hasBeenAssigned_[d].resize(numCells);

    // set the arrays which show which element is hanging or cell has hanging DoFs
    isElementHanging_[d].resize(numCells,false);
    for (int c=0; c<numCells; c++) { isElementHanging_[d][c] = false; }
    if (d == dim_){
    	for (int c=0; c<numCells; c++) { hasCellHanging_[c] = false; }
    }

    for (int b=0; b<nBasisChunks(); b++)
    {
      Tabs tab2;
      SUNDANCE_MSG2(setupVerb(),tab1 << "basis chunk=" << b);
      SUNDANCE_MSG2(setupVerb(),tab2 << "basis=" << basis(b)->description());
      int nNodes = basis(b).ptr()->nReferenceDOFsWithFacets(mesh.cellType(dim_), mesh.cellType(d));
      SUNDANCE_MSG2(setupVerb(),tab2 << "nNodes per cell=" << nNodes);
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
              
              
        SUNDANCE_MSG3(setupVerb(),tab2 << "node ptrs are "
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
          
      SUNDANCE_MSG3(setupVerb(),tab2 <<
        "num nodes is " 
        << nNodesPerCell_[b][d]);

      totalNNodesPerCell_[b][d] = nNodesPerCell_[b][d];
      for (int dd=0; dd<d; dd++) 
      {
        totalNNodesPerCell_[b][d] 
          += numFacets_[d][dd]*nNodesPerCell_[b][dd];
      }
      totalNDofsPerCell_[b][d] = totalNNodesPerCell_[b][d] * nFuncs(b);
      SUNDANCE_MSG3(setupVerb(),tab2 << "totalNDofsPerCell " << totalNDofsPerCell_[b][d]);

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
  SUNDANCE_MSG1(setupVerb(),tab << "done allocating DOF map");
}

void MixedDOFMapHN::initMap()
{
  Tabs tab;
  SUNDANCE_MSG1(setupVerb(),tab << "initializing DOF map");
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
    	// cells can not be hanging
        int cellGID = mesh().mapLIDToGID(dim_, cellLID);
        SUNDANCE_MSG4(setupVerb(),"proc=" << comm().getRank()
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
        // look up the LIDs of the facets

        mesh().getFacetArray(dim_, cellLID, d, 
          facetLID, facetOrientations);
        // for each facet, process its DOFs
        for (int f=0; f<nf; f++)
        {
          // if the facet's DOFs have been assigned already,
          // we're done
          if (!hasBeenAssigned(d, facetLID[f]))
          {
            markAsAssigned(d, facetLID[f]);
            // the facet may be owned by another processor, if the facet is hanging element then not
            if ( isRemote(d, facetLID[f], owner) && (!mesh().isElementHangingNode(d,facetLID[f])))
            {
              int facetGID 
                = mesh().mapLIDToGID(d, facetLID[f]);
              SUNDANCE_MSG4(setupVerb(),"proc=" << comm().getRank()
                << " thinks d-" << d 
                << " cell GID=" << facetGID
                << " is owned remotely by p=" << owner);
              remoteCells[d][owner].append(facetGID);
            }
            else // we can assign a DOF locally
            {
              // assign DOF
              for (int b=0; b<nBasisChunks(); b++)
              {
                setDOFs(b, d, facetLID[f], nextDOF);
              }
              // record the orientation wrt the maximal cell
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
  SUNDANCE_MSG1(setupVerb(),tab << "done initializing DOF map , numLocalDOFs:" << numLocalDOFs);
}


void MixedDOFMapHN::setDOFs(int basisChunk, int cellDim, int cellLID,
  int& nextDOF, bool isRemote)
{
  int nDofs = nDofsPerCell_[basisChunk][cellDim];
  if (nDofs==0) return;
  Tabs tab;
  SUNDANCE_MSG4(setupVerb(),tab << "setting DOFs for " << cellDim
    << "-cell " << cellLID <<  " , basisChunk:" << basisChunk <<" , nDofs:" <<nDofs << " , nextDOF:" << nextDOF );

  int* ptr = getInitialDOFPtrForCell(cellDim, cellLID, basisChunk);
  //int setupVerb = 6;

  //dofs_[cellDim][basisChunk][cellLID*nDofsPerCell_[basisChunk][cellDim]]
  //SUNDANCE_MSG1(setupVerb," dofs_[cellDim][basisChunk].size():" << dofs_[cellDim][basisChunk].size() );

  if (isRemote)
  {
	if (mesh().isElementHangingNode(cellDim, cellLID)){
		isElementHanging_[cellDim][cellLID] = true;
		for (int i=0; i<nDofs; i++) {//hanging cell, no global DoF
			ptr[i] = -1;
			//addGhostIndex(nextDOF);
		}
	}
	else
	{  // not hanging cell
		for (int i=0; i<nDofs; i++, nextDOF++){
			ptr[i] = nextDOF;
			addGhostIndex(nextDOF);
		}
	}
  }
  else
  {
	if (mesh().isElementHangingNode(cellDim, cellLID)){
		isElementHanging_[cellDim][cellLID] = true;
		for (int i=0; i<nDofs; i++) {//hanging cell, no global DoF
			//SUNDANCE_MSG1(setupVerb," cellDim:" << cellDim << ",cellLID:" << cellLID << ",basisChunk:" << basisChunk << "  ELEM HANGING");
			//SUNDANCE_MSG1(setupVerb," cellDim:" << cellDim << ",cellLID:" << cellLID << ",basisChunk:" << basisChunk << "  dof:" << nextDOF);
			//SUNDANCE_MSG1(setupVerb," writing at index:" << (cellLID*nDofsPerCell_[basisChunk][cellDim]+ i) );
			//SUNDANCE_MSG1(setupVerb," ptr[i]:" << ptr[i] );
			ptr[i] = -1;
		}
	}
	else
	{  // not hanging DoF
		for (int i=0; i<nDofs; i++,nextDOF++) {
			//SUNDANCE_MSG1(setupVerb," cellDim:" << cellDim << ",cellLID:" << cellLID << ",basisChunk:" << basisChunk << "  dof:" << nextDOF);
			//SUNDANCE_MSG1(setupVerb," writing at index:" << (cellLID*nDofsPerCell_[basisChunk][cellDim]+ i) );
			//SUNDANCE_MSG1(setupVerb," ptr[i]:" << ptr[i] );
			ptr[i] = nextDOF;
		}
	}
  }
}





void MixedDOFMapHN::shareDOFs(int cellDim,
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

    SUNDANCE_MSG4(setupVerb(),"p=" << mesh().comm().getRank()
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


RCP<const MapStructure> MixedDOFMapHN
::getDOFsForCellBatch(int cellDim,
  const Array<int>& cellLID,
  const Set<int>& requestedFuncSet,
  Array<Array<int> >& dofs,
  Array<int>& nNodes,
  int verbosity) const 
{
  TimeMonitor timer(batchedDofLookupTimer());

  Tabs tab;
  //verbosity = 6; // hard code, eliminate this
  SUNDANCE_MSG2(verbosity,
    tab << "getDOFsForCellBatch(): cellDim=" << cellDim
    << " cellLID=" << cellLID);

  dofs.resize(nBasisChunks());
  nNodes.resize(nBasisChunks());

  int nCells = cellLID.size();

  if (cellDim == dim_)
  {
    Tabs tab1;

    SUNDANCE_MSG4(verbosity,tab1 << "getting dofs for maximal cells : " << cellLID );

    for (int b=0; b<nBasisChunks(); b++)
    {
      nNodes[b] = totalNNodesPerCell_[b][cellDim];
      dofs[b].resize(nNodes[b]*nFuncs(b)*nCells);
      int dofsPerCell = nFuncs(b)*nNodes[b];
      Array<int>& chunkDofs = dofs[b];
      // Here nothing to do for hanging DoFs
      // -> this should be the only place which is called with hanging nodes
      // -> since no facet integral with hanging DoFs should occur
      for (int c=0; c<nCells; c++)
      {
        for (int i=0; i<dofsPerCell; i++)
        {
          chunkDofs[c*dofsPerCell + i] 
            = maximalDofs_[b][cellLID[c]*dofsPerCell+i];
        }
      }
      SUNDANCE_MSG3(verbosity,tab1 << "dofs for chunk-" << b << " :" << chunkDofs);
    }
  }
  else
  {
    Tabs tab1;
    SUNDANCE_MSG3(verbosity,tab1 << "getting dofs for non-maximal cells");
  
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
        SUNDANCE_MSG3(verbosity,tab2 << "cell=" << c << "    ,cellLID[c]:" << cellLID[c]);
        int offset = dofsPerCell*c;

        /* first get the DOFs for the nodes associated with 
         * the cell's interior */
        // nothing to do here for hanging DoFs, since no facet integration should occure with hanging DoF
        SUNDANCE_MSG3(verbosity,tab2 << "doing interior nodes");
        int nInteriorNodes = nNodesPerCell_[b][cellDim];
        //              int nInteriorNodes = localNodePtrs_[b][cellDim][cellDim][0].size();
        if (nInteriorNodes > 0)
        {
          if (cellDim==0 || nInteriorNodes <= 1) /* orientation-independent */
          {

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
          SUNDANCE_MSG3(verbosity,tab2 << "facet dim=" << d);
          if (nNodesPerCell_[b][d] == 0) continue;
          for (int f=0; f<numFacets[d]; f++)
          {
            Tabs tab3;
            int facetID = facetLID[d][c*numFacets[d]+f];
            SUNDANCE_MSG3(verbosity,tab2 << "f=" << f << " facetLID=" << facetID);
            if (localNodePtrs_[b][cellDim].size()==0) continue;
            int nFacetNodes = localNodePtrs_[b][cellDim][d][f].size();
            //const int* fromPtr = getInitialDOFPtrForCell(d, facetID, b);
            int* toPtr1 = &(dofs[b][dofsPerCell*c]);
            const int* nodePtr = &(localNodePtrs_[b][cellDim][d][f][0]);

            // nothing to do here for hanging DoFs, since no facet integration should occur with hanging DoF

            for (int func=0; func<nf; func++)
            {
              if (d == 0 || nFacetNodes <= 1) /* orientation-independent */
              {
                for (int n=0; n<nFacetNodes; n++)
                {
                  int ptr = nodePtr[n];
                  toPtr1[func*nNodes[b] + ptr]
                    = dofs_[d][b][facetID*nDofsPerCell_[b][d]+func*nNodesPerCell_[b][d]+n];
                  //SUNDANCE_MSG3(verbosity,tab2 << "ptr=" << ptr << " , dof= " << toPtr1[func*nNodes[b] + ptr]);
                  //SUNDANCE_MSG1(verbosity," read from index:" << (facetID*nDofsPerCell_[b][d]) << " , add func:" << func*nNodesPerCell_[b][d]+n );
                  // if this is negative , then we have to correct it
                  if (toPtr1[func*nNodes[b] + ptr] < 0 )
                     {
                	  int fIndex = 1 , tmp1 = 1;
                	  // get any cell which has this element
                	  //int maxCell = mesh().maxCofacetLID( cellDim , cellLID[c] , 0 , fIndex);
                	  // we replaced the line above with this below, due to possible bug in the 3D mesh
                	  int maxCell = mesh().maxCofacetLID( d , facetID , 0 , fIndex);
                	  //SUNDANCE_MSG3(verbosity,tab2 << "mesh().isElementHangingNode("<<cellDim <<","<<cellLID[c]<<")=" <<
                		//	  mesh().isElementHangingNode( cellDim , cellLID[c]); );
                	  //SUNDANCE_MSG3(verbosity,tab2 << "mesh().isElementHangingNode("<<d <<","<<facetID<<")=" <<
                		//	  mesh().isElementHangingNode( d , facetID); );
                	  //SUNDANCE_MSG3(verbosity,tab2 << " Point[" << facetID << "]=" <<  mesh().nodePosition(facetID) );
                      Array<int> facetsLIDs;
                      mesh().returnParentFacets( maxCell , d , facetsLIDs , tmp1 );
                      // return the parent facet at the same index
                      //SUNDANCE_MSG3(verbosity,tab2 << "maxCell=" << maxCell << " fIndex=" << fIndex );
                      //SUNDANCE_MSG3(verbosity,tab2 << "facetsLIDs[0]=" << facetsLIDs[0] << " , facetsLIDs[1]=" << facetsLIDs[1] <<
                    	//	  " , facetsLIDs[2]=" << facetsLIDs[2] << " , facetsLIDs[3]=" << facetsLIDs[3]);
                      //SUNDANCE_MSG3(verbosity,tab2 << "f=" << f << " facetLID=" << facetID << " replaced by facetsLIDs[fIndex]:" << facetsLIDs[fIndex]);
                      toPtr1[func*nNodes[b] + ptr] //= fromPtr[func*nNodes[b] + n];
                                          = dofs_[d][b][facetsLIDs[fIndex]*nDofsPerCell_[b][d]+func*nNodesPerCell_[b][d]+n];
                     }
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

                  // if this is negative , then we have to correct it
                  if (toPtr1[func*nNodes[b]+ptr] < 0)
                     {
                	  int fIndex = 1 , tmp1 = 1;
                	  // get any cell which has this element
                	  //int maxCell = mesh().maxCofacetLID( cellDim , cellLID[c] , 0 , fIndex);
                	  // we replaced the line above with this below, due to possible bug in the 3D mesh
                	  int maxCell = mesh().maxCofacetLID( d , facetID , 0 , fIndex);
                      Array<int> facetsLIDs;
                      mesh().returnParentFacets( maxCell , d , facetsLIDs , tmp1 );
                      // return the parent facet at the same index
                      // return the parent facet at the same index
                      //SUNDANCE_MSG3(verbosity,tab2 << "f=" << f << " facetLID=" << facetID << " replaced by facetsLIDs[fIndex]:" << facetsLIDs[fIndex]);
                      toPtr1[func*nNodes[b] + ptr] //= fromPtr[func*nNodes[b] + n];
                                          = dofs_[d][b][facetsLIDs[fIndex]*nDofsPerCell_[b][d]+func*nNodesPerCell_[b][d]+n];
                     }
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

void MixedDOFMapHN::buildMaximalDofTable()
{
  TimeMonitor timer(maxDOFBuildTimer());
  Tabs tab;
  int cellDim = dim_;
  int nCells = mesh().numCells(dim_);

  SUNDANCE_MSG1(setupVerb(),tab << "building dofs for maximal cells");

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
	/* localNodePtrs_[basisChunk][cellDim][facetDim][facetNumber][nodeNumber] */
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

//--------------------- Loop over all cells -------------------------

  // List of global DoFs for each basis chunk and function nr
  Array< Array<Array<int> > >globalDoFs_Cell(nBasisChunks());
  // trafo matrix for each basis chunk
  Array<Array<double> > trafoMatrix(nBasisChunks());
  // the row in the traf matrix
  Array<int>    localDoFs;
  Array<int>    parentFacetDim;
  Array<int>    parentFacetIndex;
  Array<int>    parentFacetNode;
  Array<int>    retFacets;
  Array<double> coefs;


  for (int c=0; c<nCells; c++)
  {
	// first for each cell we should find out is it has hanging elements
	// it is enough if we check the points if they are hanging
	for (int jj = 0 ; jj < numFacets[0] ; jj++){
		if (mesh().isElementHangingNode(0,facetLID[0][ c*numFacets[0] + jj])){
			  //    if order > 0, but only for that basis ... (not as a cell)
			  //       but even if we mark for this as a hanging it is OK
			  hasCellHanging_[cellLID[c]] = true;
			  break;
	    }
	}

    Tabs tab1;
    SUNDANCE_MSG3(setupVerb(),tab1 << "working on cell=" << c
      << " LID=" << cellLID[c]);
    /* first get the DOFs for the nodes associated with 
     * the cell's interior */
    SUNDANCE_MSG3(setupVerb(),tab1 << "doing interior nodes");
    for (int b=0; b<nBasisChunks(); b++)
    {
      // Initialization for cells with hanging nodes we need to create
      if (hasCellHanging_[cellLID[c]]){
    	// ----------- NH -------------
    	globalDoFs_Cell[b].resize( nFuncs(b));
    	trafoMatrix[b].resize(nNodes[b]*nNodes[b], 0.0);
    	for (int jj = 0; jj < nNodes[b]*nNodes[b] ; jj++) trafoMatrix[b][jj] = 0.0;
        for (int func=0; func<nFuncs(b); func++)
        {
        	globalDoFs_Cell[b][func].resize(nNodes[b]);
        	for (int kk = 0 ; kk < nNodes[b] ; kk++ ){
        	   globalDoFs_Cell[b][func][kk] = -1;
        	}
        }
      }

      SUNDANCE_MSG3(setupVerb(),tab1 << "basis chunk=" << b);
      if (nInteriorNodes[b]>0)
      {
    	SUNDANCE_MSG3(setupVerb(),tab1<< "nInteriorNodes = " <<nInteriorNodes[b]);
        int* toPtr = &(maximalDofs_[b][nNodes[b]*nFuncs(b)*cellLID[c]]);
        int nf = nFuncs(b);
        for (int func=0; func<nf; func++)
        {
          for (int n=0; n<nInteriorNodes[b]; n++)
          {
            /* dof(cellLID, chunk, func, node)
             * = maximalDofs_[chunk][(cellLID*nFunc + func)*nNode + node]; */
        	/* dof(cellDim, cellLID, chunk, func, node)
        	 * = dofs_[cellDim][chunk][(cellLID*nFunc + func)*nNode + node] */
        	/* localNodePtrs_[basisChunk][cellDim][facetDim][facetNumber][nodeNumber] */

            int ptr = localNodePtrs_[b][cellDim][cellDim][0][n]; // ptr is the local DoF form that function (there is only one facet)

            if (hasCellHanging_[cellLID[c]])
            {
            	SUNDANCE_MSG1(setupVerb(), tab1 << "Cell has HDoF cellLID[c]=" << cellLID[c] );
            	int dof_tmp = dofs_[cellDim][b][cellLID[c]*nDofsPerCell_[b][cellDim]+func*nNodesPerCell_[b][cellDim]+n];
                // dof_tmp can not be negative since it is inside the cell and can not be hanging
            	SUNDANCE_MSG1(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , globalDoFs_Cell[b][func]:"<<globalDoFs_Cell[b][func]<<" dof=" << dof_tmp);
            	globalDoFs_Cell[b][func][ptr] = dof_tmp;
            	SUNDANCE_MSG1(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , ptr=" << ptr <<" dof=" << dof_tmp);
                if (func == 0){
                  trafoMatrix[b][nNodes[b]*ptr + ptr] = 1.0;
                }
            }
            // internal nodes can not be hanging
            toPtr[func*nNodes[b] + ptr] = //fromPtr[func*nNodes[b] + n];
               dofs_[cellDim][b][cellLID[c]*nDofsPerCell_[b][cellDim]+func*nNodesPerCell_[b][cellDim]+n];
            SUNDANCE_MSG1(setupVerb(), tab1 << " dof=" << dofs_[cellDim][b][cellLID[c]*nDofsPerCell_[b][cellDim]+func*nNodesPerCell_[b][cellDim]+n]);
          }
        }
      }
    }
      
    SUNDANCE_MSG3(setupVerb(),tab1 << "doing facet nodes");
    /* now get the DOFs for the nodes on the facets */
    for (int d=0; d<cellDim; d++)
    {
      Tabs tab2;
      SUNDANCE_MSG3(setupVerb(),tab2 << "facet dim=" << d);

      for (int f=0; f<numFacets[d]; f++)
      {
        Tabs tab3;
        int facetID = facetLID[d][c*numFacets[d]+f];
        SUNDANCE_MSG3(setupVerb(),tab2 << "f=" << f << " facetLID=" << facetID);

        for (int b=0; b<nBasisChunks(); b++)
        {
          Tabs tab4;
          SUNDANCE_MSG3(setupVerb(),tab4 << "basis chunk=" << b);
          SUNDANCE_MSG3(setupVerb(),tab4 << "num nodes per cell=" << nNodesPerCell_[b][d]);
          int nf = nFuncs(b);
          if (nDofsPerCell_[b][d]==0) continue;
          int nFacetNodes = 0;
          if (localNodePtrs_[b][cellDim].size()!=0)
            nFacetNodes = localNodePtrs_[b][cellDim][d][f].size();
          if (nFacetNodes == 0) continue;
          //  const int* fromPtr = getInitialDOFPtrForCell(d, facetID, b);
          SUNDANCE_MSG3(setupVerb(),tab4 << "dof table size=" << maximalDofs_[b].size());
          /* = maximalDofs_[chunk][(cellLID*nFunc + func)*nNode + node]; */
          int* toPtr = &(maximalDofs_[b][nNodes[b]*nFuncs(b)*cellLID[c]]);
          /* localNodePtrs_[basisChunk][cellDim][facetDim][facetNumber][nodeNumber] */
          const int* nodePtr = &(localNodePtrs_[b][cellDim][d][f][0]);

          for (int func=0; func<nf; func++)
          {
            Tabs tab5;
            SUNDANCE_MSG3(setupVerb(),tab5 << "func=" << func);
            if (d == 0 || nFacetNodes <= 1) /* orientation-independent */
            {
              Tabs tab6;
              for (int n=0; n<nFacetNodes; n++)
              {
            	SUNDANCE_MSG3(setupVerb(),tab5 << "n=" << n);
                int ptr = nodePtr[n];
                SUNDANCE_MSG3(setupVerb(),tab6 << "ptr=" << ptr);
                /* dof(cellLID, chunk, func, node)
                 * = maximalDofs_[chunk][(cellLID*nFunc + func)*nNode + node]; */
            	/* dof(cellDim, cellLID, chunk, func, node)
            	   * = dofs_[cellDim][chunk][(cellLID*nFunc + func)*nNode + node] */
                if (hasCellHanging_[cellLID[c]]){
                	int parentCellLID = 0;
                	int dof_tmp = dofs_[d][b][facetID*nDofsPerCell_[b][d]+func*nNodesPerCell_[b][d]+n];
                	SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , d="<<d<<", f="<<f<<", n="<<n<<", dof=" << dof_tmp);
                    if (dof_tmp < 0){
                    	Array<int> constFacetLID;
                    	Array<int> constFacetDim;
                       // get the collaborating global DoFs, add the coefficients to the row of the trafo matrix (if func == 0)
                        basis_[b]->getConstrainsForHNDoF(
              				    mesh().indexInParent(cellLID[c]), cellDim ,
              				    mesh().maxChildren(), d , f , n,
              				    localDoFs, parentFacetDim,
              			    	parentFacetIndex, parentFacetNode, coefs );
                        SUNDANCE_MSG3(setupVerb(), tab1 << " After  basis_[b]->getConstrainsForHNDoF:");
                        constFacetLID.resize(localDoFs.size());
                        constFacetDim.resize(localDoFs.size());
                        for (int indexi = 0 ; indexi < localDoFs.size() ; indexi++)
                        {
                            // localNodePtrs_[basisChunk][cellDim][facetDim][facetNumber][nodeNumber]
                        	int ptr1 = localNodePtrs_[b][cellDim][parentFacetDim[indexi]]
                        	            [parentFacetIndex[indexi]][parentFacetNode[indexi]];
                        	// get the dof belonging to the parent cell
                        	mesh().returnParentFacets(cellLID[c] , parentFacetDim[indexi], retFacets , parentCellLID );
                        	int parFacetLID = retFacets[parentFacetIndex[indexi]];
                        	int parFacetNode = parentFacetNode[indexi];
                        	SUNDANCE_MSG3(setupVerb(), tab1 << "returnParentFacets , parFacetLID=" << parFacetLID <<" parFacetNode=" << parFacetNode );
                        	SUNDANCE_MSG3(setupVerb(), tab1 << "returnParentFacets , retFacets=" << retFacets );
                        	// store the facet LID and its index
                        	constFacetDim[indexi] = parentFacetDim[indexi];
                        	constFacetLID[indexi] = parFacetLID;
                        	// set the DoF in the array, and add line to the
                       	    dof_tmp = dofs_[parentFacetDim[indexi]][b]
                       	     [parFacetLID * nDofsPerCell_[b][parentFacetDim[indexi]]+func*nNodesPerCell_[b][parentFacetDim[indexi]]+parFacetNode];

                        	globalDoFs_Cell[b][func][ptr1] = dof_tmp;

                        	SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , ptr=" << ptr << "ptr1=" << ptr1 <<" dof=" << dof_tmp);
							if (func == 0){
								trafoMatrix[b][nNodes[b]*ptr + ptr1] = coefs[indexi];
							}
							SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , DONE Fill");
                        }
                        //only for hanging points store the components
                        if ( (d == 0) && (func == 0)) {
                        	HN_To_globalFacetsLID_.put( nPoints_*b + facetID , constFacetLID);
                        	HN_To_globalFacetsDim_.put( nPoints_*b + facetID , constFacetDim);
                        	HN_To_coeffs_.put( nPoints_*b + facetID , coefs);
                        }
                    }else {
                       // just put the DoF in the list , and in the matrix add one row (unit) (if func == 0)
                    	SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , ptr=" << ptr <<" dof=" << dof_tmp);
                        //SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , trafoMatrix[b].size()" <<
                        	//	trafoMatrix[b].size() << ", nNodes[b]:" << nNodes[b]);
                        //SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , globalDoFs_Cell[b].size():" <<
                        	//	globalDoFs_Cell[b].size() << ", globalDoFs_Cell[b][func].size():" << globalDoFs_Cell[b][func].size());
                    	globalDoFs_Cell[b][func][ptr] = dof_tmp;
                        if (func == 0){
                          trafoMatrix[b][nNodes[b]*ptr + ptr] = 1.0;
                        }
                    }
                } else {
                	toPtr[func*nNodes[b] + ptr]
                	      = dofs_[d][b][facetID*nDofsPerCell_[b][d]+func*nNodesPerCell_[b][d]+n];
                	SUNDANCE_MSG3(setupVerb(), tab1 << " dof=" << dofs_[d][b][facetID*nDofsPerCell_[b][d]+func*nNodesPerCell_[b][d]+n]);
                }
              }
            }
            else /* orientation-dependent */
            {
              int facetOrientation = facetOrientations[d][c*numFacets[d]+f]
                * originalFacetOrientation_[d-1][facetID];
              for (int m=0; m<nFacetNodes; m++)
              {
                int n = m;
                if (facetOrientation<0) n = nFacetNodes-1-m; // this means we have to take the DoFs in reverse order
                int ptr = nodePtr[m];
                /* dof(cellLID, chunk, func, node)
                 * = maximalDofs_[chunk][(cellLID*nFunc + func)*nNode + node]; */
            	/* dof(cellDim, cellLID, chunk, func, node)
            	   * = dofs_[cellDim][chunk][(cellLID*nFunc + func)*nNode + node] */

                if (hasCellHanging_[cellLID[c]])
                {
                	int parentCellLID = 0;
                	int dof_tmp = dofs_[d][b][facetID*nDofsPerCell_[b][d]+func*nNodesPerCell_[b][d]+n];
                	SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , d="<<d<<", f="<<f<<", n="<<n<<", dof=" << dof_tmp);
                    if (dof_tmp < 0){
                    	Array<int> constFacetLID;
                    	Array<int> constFacetDim;
                       // get the collaborating global DoFs, add the coefficients to the row of the trafo matrix (if func == 0)
                        basis_[b]->getConstrainsForHNDoF(
              				    mesh().indexInParent(cellLID[c]), cellDim ,
              				    mesh().maxChildren(), d , f , n,
              				    localDoFs, parentFacetDim,
              			    	parentFacetIndex, parentFacetNode, coefs );
                        SUNDANCE_MSG3(setupVerb(), tab1 << " After  basis_[b]->getConstrainsForHNDoF:");
                        constFacetLID.resize(localDoFs.size());
                        constFacetDim.resize(localDoFs.size());
                        for (int indexi = 0 ; indexi < localDoFs.size() ; indexi++)
                        {
                            // localNodePtrs_[basisChunk][cellDim][facetDim][facetNumber][nodeNumber]
                        	int ptr1 = localNodePtrs_[b][cellDim][parentFacetDim[indexi]]
                        	            [parentFacetIndex[indexi]][parentFacetNode[indexi]];
                        	// get the dof belonging to the parent cell
                        	mesh().returnParentFacets(cellLID[c] , parentFacetDim[indexi], retFacets , parentCellLID );
                        	int parFacetLID = retFacets[parentFacetIndex[indexi]];
                        	int parFacetNode = parentFacetNode[indexi];
                        	SUNDANCE_MSG3(setupVerb(), tab1 << "returnParentFacets , parFacetLID=" << parFacetLID <<" parFacetNode=" << parFacetNode );
                        	SUNDANCE_MSG3(setupVerb(), tab1 << "returnParentFacets , retFacets=" << retFacets );
                        	// store the facet LID and its index
                        	constFacetDim[indexi] = parentFacetDim[indexi];
                        	constFacetLID[indexi] = parFacetLID;
                        	// set the DoF in the array, and add line to the
                       	   dof_tmp = dofs_[parentFacetDim[indexi]][b]
                       	     [parFacetLID * nDofsPerCell_[b][parentFacetDim[indexi]]+func*nNodesPerCell_[b][parentFacetDim[indexi]]+parFacetNode];

                        	globalDoFs_Cell[b][func][ptr1] = dof_tmp;

                        	SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , ptr=" << ptr1 <<" dof=" << dof_tmp);
							if (func == 0){
								trafoMatrix[b][nNodes[b]*ptr + ptr1] = coefs[indexi];
							}
                        }
                        //only for hanging points store the components
                        if ( (d == 0) && (func == 0)) {
                        	HN_To_globalFacetsLID_.put( nPoints_*b + facetID , constFacetLID);
                        	HN_To_globalFacetsDim_.put( nPoints_*b + facetID , constFacetDim);
                        	HN_To_coeffs_.put( nPoints_*b + facetID , coefs);
                        }
                    }else{
                       // just put the DoF in the list , and in the matrix add one row (unit) (if func == 0)
                        //returnInsertionPosition( globalDoFs_Cell[b][func] , dof_tmp , insertPos);
                    	SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , ptr=" << ptr <<" dof=" << dof_tmp);
                    	globalDoFs_Cell[b][func][ptr] = dof_tmp;
                        if (func == 0){
                          trafoMatrix[b][nNodes[b]*ptr + ptr] = 1.0;
                        }
                    }
                } else {
                	toPtr[func*nNodes[b]+ptr]
                	      = dofs_[d][b][facetID*nDofsPerCell_[b][d]+func*nNodesPerCell_[b][d]+n];
                	SUNDANCE_MSG3(setupVerb(), tab1 << " dof=" << dofs_[d][b][facetID*nDofsPerCell_[b][d]+func*nNodesPerCell_[b][d]+n]);
                }
              }
            }
          } //func
        } //basischunk
      } // facetInd
    } // facetDim

    //in case of cells with hanging node write the collected information to the maximalDofs_ and store the trafo matrix
    if (hasCellHanging_[cellLID[c]]){
		// store tranfo Matrix
    	Array<int> matrixIndexes(nBasisChunks());
    	for (int b=0; b<nBasisChunks(); b++){
    		matrixIndexes[b] = matrixStore_.addMatrix( b , trafoMatrix[b] );
    	}
		maxCellLIDwithHN_to_TrafoMatrix_.put( cellLID[c] , matrixIndexes );

    	for (int b=0; b<nBasisChunks(); b++)
		{
    		// display trafo matrix per basis chunk
    		SUNDANCE_MSG3(setupVerb(), "trafoMatrix[b]=" << trafoMatrix[b]);
			for (int func=0; func<nFuncs(b); func++)
			{
				// display global DoFs for this basis chunk and this function inside the chunk
				SUNDANCE_MSG1(setupVerb(), "globalDoFs_Cell[b][func]=" << globalDoFs_Cell[b][func]);
				for (int jj=0 ; jj < nNodes[b] ; jj++){
				  // store global DoFs for this cell
				  maximalDofs_[b][(cellLID[c]*nFuncs(b) + func)*nNodes[b] + jj] = globalDoFs_Cell[b][func][jj];
				}
			}
		}
    }
  } // -------------- cell iteration -------------

  haveMaximalDofs_ = true;

  SUNDANCE_MSG1(setupVerb(),tab << "done building dofs for maximal cells");
}


void MixedDOFMapHN::getTrafoMatrixForCell(
	    int cellLID,
	    int funcID,
	    int& trafoMatrixSize,
	    bool& doTransform,
	    Array<double>& transfMatrix ) const {

	trafoMatrixSize = 0;

	if (maxCellLIDwithHN_to_TrafoMatrix_.containsKey(cellLID))
	{
		// return the transformation matrix from the Map
		const Array<int> &matrixIndexes = maxCellLIDwithHN_to_TrafoMatrix_.get( cellLID );
		matrixStore_.getMatrix( chunkForFuncID(funcID) , matrixIndexes[chunkForFuncID(funcID)] , transfMatrix );
		//transfMatrix = (maxCellLIDwithHN_to_TrafoMatrix_.get( cellLID ))[chunkForFuncID(funcID)]; // this should return a valid array

    // KL added cast to double to avoid compilation problems on windows
		trafoMatrixSize = sqrt((double) transfMatrix.size());
		SUNDANCE_MSG1(setupVerb(), "getTrafoMatrixForCell() cellLID:" << cellLID << ",funcID:" <<
				funcID << ",chunkForFuncID(funcID):" << chunkForFuncID(funcID) << ", trafoMatrixSize:" << trafoMatrixSize);
		//SUNDANCE_MSG1(setupVerb(), "getTrafoMatrixForCell() Matrix:" << std::endl << transfMatrix );
		doTransform = true;
	}
	else  // no transformation needed, return false
	{
		doTransform = false;
	}
}

void MixedDOFMapHN::getTrafoMatrixForFacet(
	  int cellDim,
	  int cellLID,
	  int facetIndex,
	  int funcID,
	  int& trafoMatrixSize,
	  bool& doTransform,
	  Array<double>& transfMatrix ) const {

	int fIndex;
	int maxCellLID;
	// here we ask for cofacet 0 , assuming that these are anyhow boundary facets
	SUNDANCE_MSG2(setupVerb() , "NodalDOFMapHN::getTrafoMatrixForFacet() cellDim :" << cellDim << ", cellLID:" << cellLID);
	maxCellLID = mesh().maxCofacetLID( cellDim, cellLID, 0 , fIndex);
	SUNDANCE_MSG2(setupVerb() , "NodalDOFMapHN::getTrafoMatrixForFacet() testing :" << maxCellLID);

	// todo: we might pre filter cases when this is not necessary

	if (maxCellLIDwithHN_to_TrafoMatrix_.containsKey(maxCellLID))
	{
		const Array<int> &matrixIndexes = maxCellLIDwithHN_to_TrafoMatrix_.get( maxCellLID );
		matrixStore_.getMatrix( chunkForFuncID(funcID) , matrixIndexes[chunkForFuncID(funcID)] , transfMatrix );
		doTransform = true;
	}
	else  // no transformation needed, return false
	{
		doTransform = false;
	}
}

/** Function for nodal plotting */
void MixedDOFMapHN::getDOFsForHNCell(
	  int cellDim,
	  int cellLID,
      int funcID,
      Array<int>& dofs ,
      Array<double>& coefs ) const {

	// here we can make the asumption that there is only one node at this element
	// and that the element is a point (since this is called only for at the plotting)

	// treat the cells only of dimension zero
	if (  (cellDim == 0) && ( HN_To_globalFacetsLID_.containsKey(nPoints_*chunkForFuncID(funcID) + cellLID)) )
	{
		Array<int> facetLIDs;
		Array<int> facetDim;
		SUNDANCE_MSG1(setupVerb(), "getDOFsForHNCell() cellDim:"<<cellDim<<" cellLID:"<<
				cellLID<<"funcID:" << funcID <<",chunkForFuncID(funcID)" << chunkForFuncID(funcID));
		facetLIDs = HN_To_globalFacetsLID_.get( nPoints_*chunkForFuncID(funcID) + cellLID );
		facetDim = HN_To_globalFacetsDim_.get( nPoints_*chunkForFuncID(funcID) + cellLID );
		dofs.resize(facetLIDs.size());
		int b = chunkForFuncID(funcID);
		int func = indexForFuncID(funcID);
		// return the DoFs belonging to the points and function which collaborate to the hanging node
		for (int ind = 0 ; ind < facetLIDs.size() ; ind++){
			// dofs_[cellDim][chunk][(cellLID*nFunc + func)*nNode + node]
			dofs[ind] = dofs_[facetDim[ind]][b]
			     [facetLIDs[ind]*nDofsPerCell_[b][facetDim[ind]]+func*nNodesPerCell_[b][facetDim[ind]] + 0 ]; // node = 0
		}
		// get the coefficients
		coefs = HN_To_coeffs_.get( nPoints_*chunkForFuncID(funcID) + cellLID );
		SUNDANCE_MSG1(setupVerb(),  "b=" << b);
		SUNDANCE_MSG1(setupVerb(),  "func=" << func);
		SUNDANCE_MSG1(setupVerb(),  "facetLIDs=" << facetLIDs);
		SUNDANCE_MSG1(setupVerb(),  "facetDim = " << facetDim);
		SUNDANCE_MSG1(setupVerb(),  "dofs=" << dofs);
		SUNDANCE_MSG1(setupVerb(),  "coefs = " << coefs);
	}
}


void MixedDOFMapHN::computeOffsets(int dim, int localCount)
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
      // increment only when this is not marked as hanging
      if (dofs_[dim][chunk][n] >= 0) {
    	   dofs_[dim][chunk][n] += myOffset;
      }
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


/* not used right now */
void MixedDOFMapHN::checkTable() const
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
