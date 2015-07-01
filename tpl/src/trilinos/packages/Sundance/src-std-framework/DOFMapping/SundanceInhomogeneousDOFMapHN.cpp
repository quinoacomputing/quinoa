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

#include "SundanceInhomogeneousDOFMapHN.hpp"
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

InhomogeneousDOFMapHN::InhomogeneousDOFMapHN(const Mesh& mesh,
  const Array<RCP<BasisDOFTopologyBase> >& basis,
  const Array<Map<Set<int>, CellFilter> >& funcSetToDomainMap,
  int setupVerb)
 : HNDoFMapBase(mesh, basis.size(), setupVerb),
   DOFMapBase(mesh, setupVerb),
    dim_(mesh.spatialDim()),
    dofs_(mesh.spatialDim()+1),
    funcDefined_(mesh.spatialDim()+1),
    maximalDofs_(),
    haveMaximalDofs_(false),
    localNodePtrs_(),
    maxSubdomains_(),
    allFuncs_(),
    elemFuncSets_(),
    hasCellHanging_(),
    isElementHanging_(mesh.spatialDim()+1),
    HN_To_globalFacetsLID_(),
    HN_To_globalFacetsDim_(),
    HN_To_coeffs_(),
    maxCellLIDwithHN_to_TrafoMatrix_(),
    matrixStore_(),
    nNodesPerCell_(),
    cellHasAnyDOFs_(dim_+1),
    numFacets_(mesh.spatialDim()+1),
    originalFacetOrientation_(2),
    hasBeenAssigned_(mesh.spatialDim()+1),
    structure_(),
    nFuncs_()
{
  TimeMonitor timer(mixedHNDOFCtorTimer());
  Tabs tab;

  SUNDANCE_MSG1(setupVerb, tab << "building mixed inhomogeneous hanging node DOF map");

  int nMaxCell = mesh.numCells(dim_);

// ---------- DEAL WITH INHOMOGENITY ---------------
// we only consider the maxcell filters, if considering only this results in error then it is the users fault ;-)

  for (Map<Set<int>, CellFilter>::const_iterator
         i=funcSetToDomainMap[dim_].begin(); i!=funcSetToDomainMap[dim_].end(); i++)
  {
	  allFuncs_.merge(i->first);
	  elemFuncSetsDomain_.append(i->first);
      maxSubdomains_.append(i->second);
  }

  nrAllFuncs_ = allFuncs_.size();
  Array<Set<int> > maxCellToFuncSetMap(mesh.numCells(dim_));
  Array<Array<int> > cellLIDs(maxSubdomains_.size());
  funcDomains_.resize(nrAllFuncs_);

  SUNDANCE_MSG1(setupVerb, tab << " getting the cell set for each cell ");
  for (int r=0; r<maxSubdomains_.size(); r++)
  {
    int d = maxSubdomains_[r].dimension(mesh);
    CellSet cells = maxSubdomains_[r].getCells(mesh);
    SUNDANCE_MSG2(setupVerb, " dim:" << d << " , domain " << maxSubdomains_[r] << " has functions "
      << elemFuncSetsDomain_[r]);

    for (CellIterator c=cells.begin(); c!=cells.end(); c++)
    {
      int cellLID = *c;
      cellLIDs[r].append(cellLID);
    }
    if (cellLIDs[r].size()==0) continue;

    // for each filter add to the cell the function (merge)
    for (int c=0; c<cellLIDs[r].size(); c++)
    {
    	maxCellToFuncSetMap[cellLIDs[r][c]].merge(elemFuncSetsDomain_[r]);
    }
    // add for each function the cell filter where it is defined
    Array<Array<int> > elemFuncs = tuple(elemFuncSetsDomain_[r].elements());
    for (int i=0; i<elemFuncs[0].size(); i++)
    {
      funcDomains_[elemFuncs[0][i]] = funcDomains_[elemFuncs[0][i]] + maxSubdomains_[r];
    }
  }

  elemFuncSets_.resize(0);
  maxCellFuncSetsIndex_.resize(nMaxCell);
  Map<Set<int>, int> fsToFSIndexMap;

  SUNDANCE_MSG1(setupVerb, tab << " For each cell get the correct function combination ");
  /** for each cell store the index in the set array */
  for (int n=0; n<nMaxCell; n++)
  {
    const Set<int>& f = maxCellToFuncSetMap[n];
    if (f.size()==0) continue;
    int funcComboIndex;
    if (!fsToFSIndexMap.containsKey(f))
    {
      funcComboIndex = elemFuncSets_.size();
      fsToFSIndexMap.put(f, funcComboIndex);
      SUNDANCE_MSG2(setupVerb, " New combination found : " << f );
      elemFuncSets_.append(f);
      funcSetCellCount_.append(1);
    }
    else
    {
      funcComboIndex = fsToFSIndexMap.get(f);
      funcSetCellCount_[funcComboIndex]++;
    }
    SUNDANCE_MSG2(setupVerb, " CellLID: "<< n << " combiIndex:" << funcComboIndex << " , function set: "<< f );
    maxCellFuncSetsIndex_[n] = funcComboIndex;
  }
  SUNDANCE_MSG1(setupVerb, tab << " END deal with inhomogeneity variables ");
// ---------- END DEAL WITH INHOMOGENITY ---------------

  Sundance::Map<OrderedHandle<BasisDOFTopologyBase>, int> basisToChunkMap;
  Array<RCP<BasisDOFTopologyBase> > chunkBases;
  Array<Array<int> > chunkFuncs;
  
  int chunk = 0;
  int nBasis = basis.size();
  basis_.resize(nBasis);

  nPoints_ = mesh.numCells(0);
  // this must be satisfied
  //TEUCHOS_TEST_FOR_EXCEPTION( nrAllFuncs_ != nBasis , std::runtime_error,
  //   " nrAllFuncs_ != nBasis , nrAllFuncs_=" <<  nrAllFuncs_ << " , nBasis = " << nBasis );
  if ( nrAllFuncs_ < nBasis) nrAllFuncs_ = nBasis;

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

  SUNDANCE_MSG1(setupVerb, tab << "done building mixed HN DOF map");

}


void InhomogeneousDOFMapHN::allocate(const Mesh& mesh)
{
  Tabs tab;

  /* gather functions into chunks sharing identical basis functions */
  SUNDANCE_MSG1(setupVerb(),tab << "InhomogeneousDOFMapHN::allocate() grouping like basis functions");

  /* now that we know the number of basis chunks, allocate arrays */
  localNodePtrs_.resize(nBasisChunks());
  nNodesPerCell_.resize(nBasisChunks());
  nDofsPerCell_.resize(nBasisChunks());
  maximalDofs_.resize(mesh.numCells(dim_));
  totalNNodesPerCell_.resize(nBasisChunks());

  for (int b=0; b<nBasisChunks(); b++)
  {
    localNodePtrs_[b].resize(mesh.spatialDim()+1);
    nNodesPerCell_[b].resize(mesh.spatialDim()+1);
    nDofsPerCell_[b].resize(mesh.spatialDim()+1);
    totalNNodesPerCell_[b].resize(mesh.spatialDim()+1);
  }

  /* compute node counts for each cell dimension and each basis type */
  SUNDANCE_MSG1(setupVerb(),tab << "working out DOF map node counts");
  
  numFacets_.resize(dim_+1);
  isElementHanging_.resize(dim_+1);
  hasCellHanging_.resize(mesh.numCells(dim_),false);

  /* allocate the maximal dof array */
  for (int i = 0 ; i < mesh.numCells(dim_) ; i++){
    maximalDofs_[i].resize(nrAllFuncs_);
  }


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

    int numCells = mesh.numCells(d);
    hasBeenAssigned_[d].resize(numCells);
    dofs_[d].resize(numCells);
    funcDefined_[d].resize(nrAllFuncs_);

    // set the arrays which show which element is hanging or cell has hanging DoFs
    isElementHanging_[d].resize(numCells,false);
    for (int c=0; c<numCells; c++) { isElementHanging_[d][c] = false; }
    if (d == dim_){
    	for (int c=0; c<numCells; c++) { hasCellHanging_[c] = false; }
    }

    /* allocate the DOFs array for this dimension */
    for (int i = 0 ; i < numCells ; i++){
    	dofs_[d][i].resize(nrAllFuncs_);
    }
    for (int f = 0 ; f < nrAllFuncs_ ; f++){
		funcDefined_[d][f].resize( numCells , 0 );
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
          mesh.cellType(d),  localNodePtrs_[b][d]);
              
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
        nDofsPerCell_[b][d] = nNodesPerCell_[b][d];
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
      SUNDANCE_MSG3(setupVerb(),tab2 << "num nodes is "  << nNodesPerCell_[b][d]);

      totalNNodesPerCell_[b][d] = nNodesPerCell_[b][d];
      for (int dd=0; dd<d; dd++)
      {
        totalNNodesPerCell_[b][d]
          += numFacets_[d][dd]*nNodesPerCell_[b][dd];
      }
      SUNDANCE_MSG3(setupVerb(),tab2 << "totalNNodesPerCell_[b][d] " << totalNNodesPerCell_[b][d]);
    }

    /* allocate the array of original facet orientations */
    if (d > 0 && d < dim_) 
    {
      originalFacetOrientation_[d-1].resize(numCells);
    }

  } // end from dimension loop

  SUNDANCE_MSG1(setupVerb(),tab << "done allocating DOF map");
}

void InhomogeneousDOFMapHN::initMap()
{
  Tabs tab;
  SUNDANCE_MSG1(setupVerb(),tab << "InhomogeneousDOFMapHN::initMap() initializing DOF map");
  /* start the DOF count at zero. */
  int nextDOF = 0;

  /* Space in which to keep a list of remote cells needed by each processor
   * for each dimension. The first index is dimension, the second proc, the
   * third cell number. */
  Array<Array<Array<int> > > remoteCells(mesh().spatialDim()+1);
  Array<int> functionValid(nrAllFuncs_);

  for (int d=0; d<remoteCells.size(); d++) 
    remoteCells[d].resize(mesh().comm().getNProc());

  int nrMaxCell = mesh().numCells(dim_);

  // first we set the "funcDefined_" variable which shows at each element which function
  // is defined and which is not
  SUNDANCE_MSG1(setupVerb(),tab << " starting detetecting which cell has which functions");
  for (int cellLID = 0; cellLID < nrMaxCell ; cellLID++)
  {
	  Tabs tab2;
	  for (int nrf=0; nrf < nrAllFuncs_ ; nrf++) functionValid[nrf] = 0;

	  const Set<int>& functs = elemFuncSets_[maxCellFuncSetsIndex_[cellLID]];

	  SUNDANCE_MSG1(setupVerb() , tab2 << " maxCell LID:" << cellLID << " has functions:" << functs);
	  // here we store on this maxCell which functions are defined
	  for (int nrf=0; nrf < nrAllFuncs_ ; nrf++) {
	       if (functs.contains(nrf)){
	    	   functionValid[nrf] = 1;
	    	   funcDefined_[dim_][nrf][cellLID] = 1;
	    	   SUNDANCE_MSG1(setupVerb() , tab2 << " functionID:" << nrf << " added to function set" );
	       } else {
	    	   //funcDefined_[dim_][nrf][cellLID] = 0;
	    	   SUNDANCE_MSG1(setupVerb() , tab2 << " functionID:" << nrf << " is not defined on this cell" );
	       }
	  }

	  // mark all the facets that they have also the functions which the maxcell has
	  for (int d=0; d<dim_; d++)
	  {
        int nf = numFacets_[dim_][d];
        Array<int> facetLID(nf);
        Array<int> facetOrientations(nf);
        // look up the LIDs of the facets
        mesh().getFacetArray(dim_, cellLID, d,  facetLID, facetOrientations);

        // for each facet, process its DOFs
        for (int f=0; f<nf; f++)
        {
        	for (int nrf=0; nrf < nrAllFuncs_ ; nrf++){
 	    	   SUNDANCE_MSG1(setupVerb() , tab2 << " dim:" << d << " facetLID:" << facetLID[f]
 	    	                     << " functionID:" << nrf << " is defined:" << functionValid[nrf]);
        		if (functionValid[nrf] == 1) {
        			funcDefined_[d][nrf][facetLID[f]] = 1;
        		} else {
        			//funcDefined_[d][nrf][facetLID[f]] = 0;
        		}
        	}
        }
	  }
  } // end of cell iteration
  

  /* we run trough each cell and for each element we test which functions have DOFs */
  for (int cellLID = 0; cellLID < nrMaxCell ; cellLID++)
  {
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
        for (int nrf=0; nrf < nrAllFuncs_ ; nrf++)
        {
          // assign DoFs in the interior of the cell
          // we try to assign DoFs to all functions i nthis cell
          setDOFs(nrf, dim_, cellLID, nextDOF);
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

        mesh().getFacetArray(dim_, cellLID, d,  facetLID, facetOrientations);

        // for each facet, process its DOFs
        for (int f=0; f<nf; f++)
        {
          // if the facet's DOFs have been assigned already, we're done
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
              for (int nrf=0; nrf < nrAllFuncs_ ; nrf++)
              {
                // we try assign DoFs for all functions
                setDOFs( nrf , d, facetLID[f], nextDOF);
              }
              // record the orientation wrt the maximal cell
              if (d > 0) 
              {
                originalFacetOrientation_[d-1][facetLID[f]] 
                  = facetOrientations[f];
              }
            }
          } // end from the hasbeedAssigned
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


void InhomogeneousDOFMapHN::setDOFs(int funcID, int cellDim, int cellLID,
  int& nextDOF, bool isRemote)
{
  int basisChunk = chunkForFuncID(funcID);
  int nDofs = nDofsPerCell_[basisChunk][cellDim];
  if (nDofs==0) return;

  Tabs tab;
  SUNDANCE_MSG4(setupVerb(),tab << "setting DOFs for " << cellDim
    << "-cell " << cellLID <<  " , basisChunk:" << basisChunk <<" , nDofs:" <<nDofs << " , nextDOF:" << nextDOF
    << " , def:" << funcDefined_[cellDim][funcID][cellLID]);

  // this means the function is not defined on this element
  // so just resize the vector to zero
  if (funcDefined_[cellDim][funcID][cellLID] <= 0){
	  SUNDANCE_MSG4(setupVerb(),tab << "No DOFs defined on dim:" << cellDim << "-cell " << cellLID << " funcID:" << funcID);
	  dofs_[cellDim][cellLID][funcID].resize(0);
	  return;
  }

  // resize the dofs_ to the number of DoFs of this function on this element
  dofs_[cellDim][cellLID][funcID].resize(nDofs,-1);
  int* ptr = getInitialDOFPtrForCell(cellDim, cellLID, funcID );

  // display the size of that part of the dofs_ array
  //SUNDANCE_MSG1(setupVerb(), " dofs_[cellDim][cellLID][funcID].size(): " << dofs_[cellDim][cellLID][funcID].size() );

  if (isRemote)
  {
	// for remote nodes where one function does not exist we'll get -1, in this case do nothing and exit
	if ( nextDOF < 0 ) return;

	if (mesh().isElementHangingNode(cellDim, cellLID)){
		isElementHanging_[cellDim][cellLID] = true;
		for (int i=0; i<nDofs; i++) {//hanging cell, no global DoF
			ptr[i] = -1;
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
			//SUNDANCE_MSG1(setupVerb()," cellDim:" << cellDim << ",cellLID:" << cellLID << ",basisChunk:" << basisChunk << "  ELEM HANGING");
			ptr[i] = -1;
		}
	}
	else
	{  // not hanging DoF
		for (int i=0; i<nDofs; i++,nextDOF++) {
			//SUNDANCE_MSG1(setupVerb()," cellDim:" << cellDim << ",cellLID:" << cellLID << ",basisChunk:" << basisChunk << "  dof:" << nextDOF);
			//SUNDANCE_MSG1(setupVerb()," ptr[i]:" << ptr[i] );
			ptr[i] = nextDOF;
		}
	}
  }
}



void InhomogeneousDOFMapHN::shareDOFs(int cellDim,
  const Array<Array<int> >& outgoingCellRequests)
{

  int np = mesh().comm().getNProc();
  int rank = mesh().comm().getRank();

  Array<Array<int> > incomingCellRequests;
  Array<Array<int> > outgoingDOFs(np);
  Array<Array<int> > incomingDOFs;

  SUNDANCE_MSG2(setupVerb(),
    "InhomogeneousDOFMapHN::shareDOFs , p=" << mesh().comm().getRank()
    << "synchronizing DOFs for cells of dimension " << cellDim);
  SUNDANCE_MSG2(setupVerb(),
    "p=" << mesh().comm().getRank()
    << " sending cell reqs d=" << cellDim << " GID=" << outgoingCellRequests);

  // share the cell requests
  MPIContainerComm<int>::allToAll(outgoingCellRequests, 
    incomingCellRequests,
    mesh().comm());
  
  // we send the following information in response:
   // (1) The first DOF for each chunk for the requested cell
   // (2) The orientation if the cell is an edge or face

  int blockSize = 0;
  bool sendOrientation = false;
  for (int b=0; b<nBasisChunks(); b++)
  {
	// here we have to sum up all the functions which are defined
	int nDofs = nDofsPerCell_[b][cellDim];
	// per chunck we have to transmit all the functions DoFs (some might not exist)
	if (nDofs > 0) blockSize = blockSize + nFuncs(b);
	if (nDofs > 1 && cellDim > 0 && cellDim < dim_) sendOrientation = true;
  }
  blockSize += sendOrientation;

  SUNDANCE_MSG2(setupVerb(),
    "p=" << rank
    << "recvd DOF requests for cells of dimension " << cellDim
    << " GID=" << incomingCellRequests);

  // get orientations and DOF numbers for the first node of every cell that's been
  //  requested by someone else
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
        << " GID=" << GID << " nReq=" << nReq << " , blockSize=" << blockSize);
      int LID = mesh().mapGIDToLID(cellDim, GID);
      SUNDANCE_MSG3(setupVerb(),
        "p=" << rank
        << " LID=" << LID << " dofs=" << dofs_[cellDim][LID]);
      int blockOffset = 0;
      for (int b=0; b<nBasisChunks(); b++)
      {
        if (nDofsPerCell_[b][cellDim] <= 0) continue;
        for (int funcI = 0 ; funcI < nFuncs(b) ; funcI++){
        	// get the first function ID from this chunk
        	int funcID = getfuncID( b , funcI );
			//assign value only if the function is defined, -1 otherwise
        	//SUNDANCE_MSG3(setupVerb(), "Getting DoF for b="<<b<<", funcI="<< funcI <<
        	//		", funcID=" << funcID << " is defined:" << funcDefined_[cellDim][funcID][LID]);
			if (funcDefined_[cellDim][funcID][LID] > 0){
				outgoingDOFs[p][blockSize*c+blockOffset]
				                = getInitialDOFForCell(cellDim, LID, funcID );
			} else {
				outgoingDOFs[p][blockSize*c+blockOffset] = -1;
			}
            // increase the blockoffset
			blockOffset++;
        }
      }
      // add the orientation
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

  // share the DOF numbers
  MPIContainerComm<int>::allToAll(outgoingDOFs,
    incomingDOFs,
    mesh().comm());

  SUNDANCE_MSG2(setupVerb(),
    "p=" << mesh().comm().getRank()
    << "communicated DOF answers for cells of dimension " << cellDim );

  SUNDANCE_MSG2(setupVerb(), "My rank:" << mesh().comm().getRank() << ",  nr. Proc:" << mesh().comm().getNProc());
  
  // now assign the DOFs from the other procs
  for (int p=0; p<mesh().comm().getNProc(); p++)
  {
	SUNDANCE_MSG3(setupVerb(), "Test p="<<p);
    if (p==mesh().comm().getRank()) continue;

    const Array<int>& dofsFromProc = incomingDOFs[p];
    int numCells = dofsFromProc.size()/blockSize;

    //
	//SUNDANCE_MSG3(setupVerb(), "Setting from p="<<p<<", GID="<< outgoingCellRequests[p]);
	//SUNDANCE_MSG3(setupVerb(), "received from p="<<p<<", dofsFromProc="<< dofsFromProc);

    for (int c=0; c<numCells; c++)
    {
      int cellGID = outgoingCellRequests[p][c];
      int cellLID = mesh().mapGIDToLID(cellDim, cellGID);
      int blockOffset = 0;
      for (int b=0; b<nBasisChunks(); b++)
      {
        if (nDofsPerCell_[b][cellDim] == 0) continue;
        for (int funcI = 0 ; funcI < nFuncs(b) ; funcI++){
        	// get the first function ID from this chunk
        	int funcID = getfuncID( b , funcI );
        	int dof = dofsFromProc[blockSize*c+blockOffset];
        	// set the recieved DoF
        	SUNDANCE_MSG3(setupVerb(), "Set DoF b="<<b<<", funcI="<< funcI <<
        			", funcID=" << funcID << " is defined:" << funcDefined_[cellDim][funcID][cellLID] <<
        			 " dof=" << dof );
        	setDOFs(funcID, cellDim, cellLID, dof, true);
        	blockOffset++;
        }
      }
      if (sendOrientation) 
      {
        originalFacetOrientation_[cellDim-1][cellLID] 
          = dofsFromProc[blockSize*(c+1)-1];
      }
    }
  }
}


RCP<const MapStructure> InhomogeneousDOFMapHN
::getDOFsForCellBatch(int cellDim,
  const Array<int>& cellLID,
  const Set<int>& requestedFuncSet,
  Array<Array<int> >& dofs,
  Array<int>& nNodes,
  int verbosity) const 
{
  TimeMonitor timer(batchedDofLookupTimer());

  Tabs tab;
  //verbosity = 6;
  SUNDANCE_MSG2(verbosity,
    tab << "getDOFsForCellBatch(): cellDim=" << cellDim
    << " cellLID=" << cellLID);


  int nCells = cellLID.size();
  // the requested function elements
  Array<int> funcs = requestedFuncSet.elements();
  int nrReqFuncs = funcs.size();
  SUNDANCE_MSG2(verbosity, tab << "nrReqFuncs=" << nrReqFuncs << " functions requested:" << requestedFuncSet);

  dofs.resize(nBasisChunks());
  nNodes.resize(nBasisChunks());

  // here we specify which DoF to return at function places where
  // these values might be input for "ghostView_" object
  // as a dummy value this will be used, this dummy value should always exist
  // now we take first DoF which is local on the processor
  int DoFs_Non =  lowestLocalDOF()+1; //-1;

  if (cellDim == dim_)
  {
    Tabs tab1;

    SUNDANCE_MSG4(verbosity,tab1 << "getting dofs for maximal cells : " << cellLID );

    for (int b = 0 ; b < nBasisChunks() ; b++)
    {
    	// set the value and resize the dof array
    	nNodes[b] = totalNNodesPerCell_[ b ][cellDim];
    	// resize the array with, and put dummy DoFs in the array
    	dofs[b].resize( nNodes[b]* nFuncs(b) * nCells );
    	for (int tmpI = 0; tmpI < nNodes[b]* nFuncs(b) * nCells ; tmpI++)
    		dofs[b][tmpI] = DoFs_Non;
    }

    for (int b = 0; b < nBasisChunks() ; b++)
    {
     for (int funcI = 0; funcI < nFuncs(b) ; funcI++)
     {
      //int b = chunkForFuncID(nfc);
      //int func = indexForFuncID(nfc);
      int funcID = getfuncID(b , funcI );

      int nrFunc = nFuncs(b);
      int dofsPerCell = nNodes[b];
      Array<int>& chunkDofs = dofs[b];
      // for each cell just copy the DoFs of this function
      for (int c=0; c<nCells; c++)
      {
    	// if the function is not defined then skip this cell
        if (funcDefined_[cellDim][funcID][cellLID[c]] <= 0) {
            SUNDANCE_MSG3(verbosity,tab << " Skip this cell since function is not defined on this cell:" << cellLID[c]);
        	continue;
        }

        for (int i=0; i<dofsPerCell; i++)
        {
          // copy the values from the "maximalDofs_" into the target array
          //SUNDANCE_MSG3(verbosity,tab1 << "i:" << i << " , fi:" << fi << ", funcs[fi]:" << funcs[fi]);
          chunkDofs[c*dofsPerCell*nrFunc + dofsPerCell*funcI + i]
            = maximalDofs_[cellLID[c]][funcID][i];
        }
      }
      SUNDANCE_MSG3(verbosity,tab1 << "dofs for chunk-" << b << " , funcI:" << funcI << " :" << chunkDofs);
     }// end function index in chunk
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

	// set the value and resize the dof array
    for (int b = 0 ; b < nBasisChunks() ; b++)
    {
    	nNodes[b] = totalNNodesPerCell_[b][cellDim];
    	// resize the array with, and put dummy DoFs in the array
    	dofs[b].resize( nNodes[b]* nFuncs(b) * nCells , -1);
    	for (int tmpI = 0; tmpI < nNodes[b]* nFuncs(b) * nCells ; tmpI++)
    		dofs[b][tmpI] = DoFs_Non;
    }

    // for each function that has been requested
    for (int b = 0; b < nBasisChunks() ; b++)
    {
     for (int funcI = 0; funcI < nFuncs(b) ; funcI++)
     {
      //int b = chunkForFuncID(nfc);
      //int func = indexForFuncID(nfc);
      int funcID = getfuncID(b , funcI );
      int nrFunc = nFuncs(b); //mapFunctions[b_r].size();
      int dofsPerCell = nrFunc*nNodes[b];
          
      Array<int>& toPtr = dofs[b];

      for (int c=0; c<nCells; c++)
      {
        Tabs tab2;
        SUNDANCE_MSG3(verbosity,tab2 << "cell=" << c << "    ,cellLID[c]:" << cellLID[c]);
        int offset = dofsPerCell*c;

    	// if the function is not defined then skip this cell
        if (funcDefined_[cellDim][funcID][cellLID[c]] <= 0) {
            SUNDANCE_MSG3(verbosity,tab2 << " Skip this cell since function is not defined on this cell:" << cellLID[c]);
        	continue;
        }

        // first get the DOFs for the nodes associated with
        // the cell's interior
        // nothing to do here for hanging DoFs, since no facet integration should occur with hanging DoF
        SUNDANCE_MSG3(verbosity,tab2 << "doing interior nodes");
        int nInteriorNodes = nNodesPerCell_[b][cellDim];
        if (nInteriorNodes > 0)
        {
          if (cellDim==0 || nInteriorNodes <= 1) // orientation-independent
          {
            for (int n=0; n<nInteriorNodes; n++)
            {
                int ptr = localNodePtrs_[b][cellDim][cellDim][0][n];
                toPtr[offset + funcI*nNodes[b] + ptr]
                  = dofs_[cellDim][cellLID[c]][funcID][n];
            }
          }
          else
          {
            int sign = originalFacetOrientation_[cellDim-1][cellLID[c]];
            int nInteriorNodes = localNodePtrs_[b][cellDim][cellDim][0].size();

            for (int m=0; m<nInteriorNodes; m++)
            {
               int n = m;
               if (sign<0) n = nInteriorNodes-1-m;
               int ptr = localNodePtrs_[b][cellDim][cellDim][0][m];
               toPtr[offset + funcI*nNodes[b] + ptr]
                 = dofs_[cellDim][cellLID[c]][funcID][n];
            }
          }
        } // end if interior nodes > 0

        // now do the facetsnInteriorNodes
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
            int* toPtr1 = &( dofs[b][dofsPerCell*c]);
            const int* nodePtr = &(localNodePtrs_[b][cellDim][d][f][0]);

            // hanging DoF can occure here only in 3D for 2D surf)
            if (d == 0 || nFacetNodes <= 1) // orientation-independent
            {
                for (int n=0; n<nFacetNodes; n++)
                {
                  int ptr = nodePtr[n];
                  toPtr1[funcI*nNodes[b] + ptr] = dofs_[d][facetID][funcID][n];

                  // if this is negative , then we have to correct it
                  if (toPtr1[funcI*nNodes[b] + ptr] < 0 )
                  {
                	  int fIndex = 1 , tmp1 = 1;
                	  // get any cell which has this element
                	  int maxCell = mesh().maxCofacetLID( cellDim , cellLID[c] , 0 , fIndex);
                      Array<int> facetsLIDs;
                      mesh().returnParentFacets( maxCell , d , facetsLIDs , tmp1 );
                      // return the parent facet at the same index
                      //SUNDANCE_MSG3(verbosity,tab2 << "f=" << f << " facetLID=" << facetID << " replaced by facetsLIDs[fIndex]:" << facetsLIDs[fIndex]);
                      toPtr1[funcI*nNodes[b] + ptr] = dofs_[d][facetsLIDs[fIndex]][funcID][n];
                   }
                }
              }
              else // orientation-dependent
              {
                int facetOrientation = facetOrientations[d][c*numFacets[d]+f]
                  * originalFacetOrientation_[d-1][facetID];
                for (int m=0; m<nFacetNodes; m++)
                {
                  int n = m;
                  if (facetOrientation<0) n = nFacetNodes-1-m;
                  int ptr = nodePtr[n];
                  toPtr1[funcI*nNodes[b] + ptr] = dofs_[d][facetID][funcID][n];

                  // if this is negative , then we have to correct it
                  if (toPtr1[funcI*nNodes[b]+ptr] < 0)
                  {
                	  int fIndex = 1 , tmp1 = 1;
                	  // get any cell which has this element
                	  int maxCell = mesh().maxCofacetLID( cellDim , cellLID[c] , 0 , fIndex);
                      Array<int> facetsLIDs;
                      mesh().returnParentFacets( maxCell , d , facetsLIDs , tmp1 );
                      // return the parent facet at the same index
                      //SUNDANCE_MSG3(verbosity,tab2 << "f=" << f << " facetLID=" << facetID << " replaced by facetsLIDs[fIndex]:" << facetsLIDs[fIndex]);
                      toPtr1[funcI*nNodes[b] + ptr] = dofs_[d][facetsLIDs[fIndex]][funcID][n];
                   }
                }
              } // orientation dependent branch
          } // end loop over facets
        } // end loop of the facets of the facet
      } // end loop over cells
      SUNDANCE_MSG3(verbosity,tab1 << "dofs for chunk-" << b << " , dofs[b]:" << dofs[b]);
     } // end loop over functions in chunk
    } // end loop over chunks
  } // else branch from the facet case
  //*/

  return structure_;
  // return the structure which was created at the beginning
  //return rtn;
}    

void InhomogeneousDOFMapHN::buildMaximalDofTable()
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

  // set up the array of the cell LIDs,
  for (int c=0; c<nCells; c++) cellLID[c]=c;
  
  // get the facets and the number of facets
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
	// localNodePtrs_[basisChunk][cellDim][facetDim][facetNumber][nodeNumber]
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
  Array<Array<int> > globalDoFs_Cell(nrAllFuncs_);
  // trafo matrix for each basis chunk
  Array<Array<double> > trafoMatrix(nBasisChunks());
  // the row in the traf matrix
  Array<int>    localDoFs;
  Array<int>    parentFacetDim;
  Array<int>    parentFacetIndex;
  Array<int>    parentFacetNode;
  Array<int>    retFacets;
  Array<double> coefs;
  int           maxCellLID;


  for (int c=0; c<nCells; c++)
  {
	maxCellLID = cellLID[c];

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
    // first get the DOFs for the nodes associated with
    // the cell's interior


    SUNDANCE_MSG3(setupVerb(),tab1 << "doing interior nodes");
    for (int b = 0; b < nBasisChunks() ; b++)
    {

     SUNDANCE_MSG3(setupVerb(),tab1 << "basis chunk=" << b);

     // loop over each function in this chunk
     for (int func = 0; func < nFuncs(b) ; func++)
     {
      //int b = chunkForFuncID(nfc);
      //int func = indexForFuncID(nfc);
      int nfc = getfuncID(b , func );

      // Initialization for cells with hanging nodes we need to create
      if (hasCellHanging_[cellLID[c]]){
       	// ----------- NH -------------
    	if (func == 0){
    		// do this only once for a chunk
    		trafoMatrix[b].resize(nNodes[b]*nNodes[b], 0.0);
    		for (int jj = 0; jj < nNodes[b]*nNodes[b] ; jj++) trafoMatrix[b][jj] = 0.0;
    	}

        globalDoFs_Cell[nfc].resize(nNodes[b]);
        for (int kk = 0 ; kk < nNodes[b] ; kk++ ){
           	globalDoFs_Cell[nfc][kk] = -1;
        }
      }


      // if on the maxcell the this function is not defined then skip this
      if (funcDefined_[dim_][nfc][maxCellLID] <= 0 ){
      	  SUNDANCE_MSG3(setupVerb(),tab1 << " function " << nfc << " , on this maxCell not defined: " << maxCellLID );
      	  maximalDofs_[maxCellLID][nfc].resize(0);
    	  continue;
      }
      else
      {   // in this case we resize the
      	  SUNDANCE_MSG3(setupVerb(),tab1 << " function " << nfc << " , is defined on maxCell : " << maxCellLID << " resize to:" << nNodes[b]);
    	  maximalDofs_[maxCellLID][nfc].resize(nNodes[b]);
      }

      if (nInteriorNodes[b]>0)
      {
    	SUNDANCE_MSG3(setupVerb(),tab1<< "nInteriorNodes = " <<nInteriorNodes[b]);
        int* toPtr = &(maximalDofs_[maxCellLID][nfc][0]);
          for (int n=0; n<nInteriorNodes[b]; n++)
          {
        	// copy all the interior nodes to the maximal Cell DOF map
        	int ptr = localNodePtrs_[b][cellDim][cellDim][0][n];

            if (hasCellHanging_[cellLID[c]])
            {
            	SUNDANCE_MSG1(setupVerb(), tab1 << "Cell has HDoF cellLID[c]=" << cellLID[c] );
            	int dof_tmp = dofs_[cellDim][maxCellLID][nfc][n];
                // dof_tmp can not be negative since it is inside the cell and can not be hanging
            	//SUNDANCE_MSG1(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , globalDoFs_Cell[nfc][n]:"<<globalDoFs_Cell[nfc][n]<<" dof=" << dof_tmp);
            	globalDoFs_Cell[nfc][ptr] = dof_tmp;
            	//SUNDANCE_MSG1(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , ptr=" << ptr <<" dof=" << dof_tmp);
                if (func == 0){
                  trafoMatrix[b][nNodes[b]*ptr + ptr] = 1.0;
                }
            }
            // internal nodes can not be hanging
            toPtr[ptr] = dofs_[cellDim][maxCellLID][nfc][n];
            //SUNDANCE_MSG1(setupVerb(), tab1 << " dof=" << dofs_[cellDim][maxCellLID][nfc][n];);
          }
      }
     } // end function in chunk
    }// end chunk
      
    SUNDANCE_MSG3(setupVerb(),tab1 << "doing facet nodes");
    // now get the DOFs for the nodes on the facets
    for (int d=0; d<cellDim; d++)
    {
      Tabs tab2;
      SUNDANCE_MSG3(setupVerb(),tab2 << "facet dim=" << d);

      for (int f=0; f<numFacets[d]; f++)
      {
        Tabs tab3;
        int facetID = facetLID[d][c*numFacets[d]+f];
        SUNDANCE_MSG3(setupVerb(),tab2 << "f=" << f << " facetLID=" << facetID);

        for (int b = 0; b < nBasisChunks() ; b++)
        {
         for (int func = 0; func < nFuncs(b) ; func++)
         {
          //int b = chunkForFuncID(nfc);
          //int func = indexForFuncID(nfc);
  		  int nfc = getfuncID(b , func );

          // if on the maxcell the this function is not defined then skip this
          if (funcDefined_[dim_][nfc][maxCellLID] <= 0 ){
          	  SUNDANCE_MSG3(setupVerb(),tab2 << " function= " << nfc << " , on this maxCell not defined: " << maxCellLID );
          	  continue;
          }
          else
          {
          	  SUNDANCE_MSG3(setupVerb(),tab2 << " function= " << nfc << " , is defined on maxcell: " << maxCellLID );
          }
          Tabs tab4;
          SUNDANCE_MSG3(setupVerb(),tab4 << "basis chunk=" << b);
          SUNDANCE_MSG3(setupVerb(),tab4 << "function index in  chunk=" << func);
          SUNDANCE_MSG3(setupVerb(),tab4 << "num nodes per cell=" << nNodesPerCell_[b][d]);
          if (nDofsPerCell_[b][d]==0) continue;
          int nFacetNodes = 0;
          if (localNodePtrs_[b][cellDim].size()!=0)
            nFacetNodes = localNodePtrs_[b][cellDim][d][f].size();
          if (nFacetNodes == 0) continue;

          SUNDANCE_MSG3(setupVerb(),tab4 << "dof table size=" << maximalDofs_[maxCellLID][nfc].size());
          // = maximalDofs_[LID][funcID][DoFIndex];
          int* toPtr = &(maximalDofs_[maxCellLID][nfc][0]);
          // localNodePtrs_[basisChunk][cellDim][facetDim][facetNumber][nodeNumber]
          const int* nodePtr = &(localNodePtrs_[b][cellDim][d][f][0]);

            Tabs tab5;
            SUNDANCE_MSG3(setupVerb(),tab5 << "func=" << func);
            if (d == 0 || nFacetNodes <= 1) // orientation-independent
            {
              Tabs tab6;
              for (int n=0; n<nFacetNodes; n++)
              {
            	//SUNDANCE_MSG3(setupVerb(),tab5 << "n=" << n);
                int ptr = nodePtr[n];
                //SUNDANCE_MSG3(setupVerb(),tab6 << "ptr=" << ptr);
                // = maximalDofs_[LID][funcID][DoFIndex];
            	// = dofs_[cellDim][cellLID][funcID][node]
                if (hasCellHanging_[cellLID[c]]){
                	int parentCellLID = 0;
                	int dof_tmp = dofs_[d][facetID][nfc][n];
                	//SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , d="<<d<<", f="<<f<<", n="<<n<<", dof=" << dof_tmp);
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
                        	//SUNDANCE_MSG3(setupVerb(), tab1 << "returnParentFacets , parFacetLID=" << parFacetLID <<" parFacetNode=" << parFacetNode );
                        	//SUNDANCE_MSG3(setupVerb(), tab1 << "returnParentFacets , retFacets=" << retFacets );
                        	// store the facet LID and its index
                        	constFacetDim[indexi] = parentFacetDim[indexi];
                        	constFacetLID[indexi] = parFacetLID;
                        	// set the DoF in the array, and add line to the
                        	// = dofs_[cellDim][cellLID][funcID][node]
                        	dof_tmp = dofs_[parentFacetDim[indexi]][parFacetLID][nfc][parFacetNode];

                        	globalDoFs_Cell[nfc][ptr1] = dof_tmp;

                        	//SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , ptr=" << ptr << "ptr1=" << ptr1 <<" dof=" << dof_tmp);
							if (func == 0){
								trafoMatrix[b][nNodes[b]*ptr + ptr1] = coefs[indexi];
							}
							//SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , DONE Fill");
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
                        //SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , globalDoFs_Cell[nfc].size():" <<
                        	//	globalDoFs_Cell[nfc].size() << ", globalDoFs_Cell[nfc].size():" << globalDoFs_Cell[nfc].size());
                    	globalDoFs_Cell[nfc][ptr] = dof_tmp;
                        if (func == 0){
                          trafoMatrix[b][nNodes[b]*ptr + ptr] = 1.0;
                        }
                    }
                } else {
                	// = dofs_[cellDim][cellLID][funcID][node]
                	toPtr[ptr] = dofs_[d][facetID][nfc][n];
                	//SUNDANCE_MSG3(setupVerb(), tab1 << " dof=" << dofs_[d][facetID][nfc][n] );
                }
              }
            }
            else // orientation-dependent
            {
              int facetOrientation = facetOrientations[d][c*numFacets[d]+f]
                * originalFacetOrientation_[d-1][facetID];
              for (int m=0; m<nFacetNodes; m++)
              {
                int n = m;
                if (facetOrientation<0) n = nFacetNodes-1-m; // this means we have to take the DoFs in reverse order
                int ptr = nodePtr[m];
                // = maximalDofs_[LID][funcID][DoFIndex];
             	// = dofs_[cellDim][cellLID][funcID][node]

                if (hasCellHanging_[cellLID[c]])
                {
                	int parentCellLID = 0;
                	int dof_tmp = dofs_[d][facetID][nfc][n];
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
                        	//SUNDANCE_MSG3(setupVerb(), tab1 << "returnParentFacets , parFacetLID=" << parFacetLID <<" parFacetNode=" << parFacetNode );
                        	//SUNDANCE_MSG3(setupVerb(), tab1 << "returnParentFacets , retFacets=" << retFacets );
                        	// store the facet LID and its index
                        	constFacetDim[indexi] = parentFacetDim[indexi];
                        	constFacetLID[indexi] = parFacetLID;
                        	// set the DoF in the array, and add line to the
                           	// = dofs_[cellDim][cellLID][funcID][node]
                            dof_tmp = dofs_[parentFacetDim[indexi]][parFacetLID][nfc][parFacetNode];

                        	globalDoFs_Cell[nfc][ptr1] = dof_tmp;

                        	//SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , ptr=" << ptr1 <<" dof=" << dof_tmp);
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
                        //returnInsertionPosition( globalDoFs_Cell[nfc] , dof_tmp , insertPos);
                    	//SUNDANCE_MSG3(setupVerb(), tab1 << "Cell has HDoF cellLID[c] , ptr=" << ptr <<" dof=" << dof_tmp);
                    	globalDoFs_Cell[nfc][ptr] = dof_tmp;
                        if (func == 0){
                          trafoMatrix[b][nNodes[b]*ptr + ptr] = 1.0;
                        }
                    }
                } else {
                   	// = dofs_[cellDim][cellLID][funcID][node]
                    toPtr[ptr] = dofs_[d][facetID][nfc][n];
                	//SUNDANCE_MSG3(setupVerb(), tab1 << " dof=" << dofs_[d][facetID][nfc][n] );
                }
              }
            } // end orientation dependent
         } //end functions in chunk
        } //end chunks
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
				int funcID = getfuncID(b , func );
				// if the function is not defined then just skip to the next one
				if (funcDefined_[dim_][funcID][maxCellLID] <= 0 ) continue;
				// display global DoFs for this basis chunk and this function inside the chunk
				SUNDANCE_MSG1(setupVerb(), "funcID=" << funcID << "globalDoFs_Cell[funcID]=" << globalDoFs_Cell[funcID]);
				for (int jj=0 ; jj < nNodes[b] ; jj++){
				  // store global DoFs for this cell
				  maximalDofs_[maxCellLID][funcID][jj] = globalDoFs_Cell[funcID][jj];
				}
			}
		}
    }
  } // -------------- cell iteration -------------
  haveMaximalDofs_ = true;

  SUNDANCE_MSG1(setupVerb(),tab << "done building dofs for maximal cells");
}


void InhomogeneousDOFMapHN::getTrafoMatrixForCell(
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
				funcID << ",chunkForFuncID(funcID):" << chunkForFuncID(funcID) <<
				" , indexForFuncID(funcID)" << indexForFuncID(funcID) << ", trafoMatrixSize:" << trafoMatrixSize);
		SUNDANCE_MSG1(setupVerb(), "getTrafoMatrixForCell() Matrix:" << std::endl << transfMatrix );
		doTransform = true;
	}
	else  // no transformation needed, return false
	{
		doTransform = false;
	}
}

void InhomogeneousDOFMapHN::getTrafoMatrixForFacet(
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
	SUNDANCE_MSG2(setupVerb() , "NodalDOFMapHN::getTrafoMatrixForFacet() cellDim :" << cellDim << ", cellLID:" << cellLID
			 << ",chunkForFuncID(funcID):" << chunkForFuncID(funcID) << " , indexForFuncID(funcID)" << indexForFuncID(funcID) );
	maxCellLID = mesh().maxCofacetLID( cellDim, cellLID, 0 , fIndex);
	SUNDANCE_MSG2(setupVerb() , "NodalDOFMapHN::getTrafoMatrixForFacet() testing :" << maxCellLID);

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
void InhomogeneousDOFMapHN::getDOFsForHNCell(
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
			dofs[ind] = dofs_[facetDim[ind]][facetLIDs[ind]][funcID][0]; // node = 0
		}
		// get the coefficients
		coefs = HN_To_coeffs_.get( nPoints_*chunkForFuncID(funcID) + cellLID );
		SUNDANCE_MSG1(setupVerb(),  "b=" << b);
		SUNDANCE_MSG1(setupVerb(),  "func=" << func);
		SUNDANCE_MSG1(setupVerb(),  "funcID=" << funcID);
		SUNDANCE_MSG1(setupVerb(),  "facetLIDs=" << facetLIDs);
		SUNDANCE_MSG1(setupVerb(),  "facetDim = " << facetDim);
		SUNDANCE_MSG1(setupVerb(),  "dofs=" << dofs);
		SUNDANCE_MSG1(setupVerb(),  "coefs = " << coefs);
	}
}

RCP<const Set<int> >
InhomogeneousDOFMapHN::allowedFuncsOnCellBatch(int cellDim,
  const Array<int>& cellLID) const
{
  Set<int> rtn;

  if (cellDim != dim_)
  {
	SUNDANCE_MSG3(setupVerb(),  "InhomogeneousDOFMapHN::allowedFuncsOnCellBatch cellDim = " << cellDim);
    for (int c=0; c<cellLID.size(); c++)
    {
      // for each maxcofacet make a merge of the functions
      int numMaxCoFs = mesh().numMaxCofacets( cellDim, cellLID[c] ) , tmp , maxCoFacet;
	  SUNDANCE_MSG3(setupVerb(),  "Get all the maxCoFacets , numMaxCoFs = " << numMaxCoFs);
	  Set<int> rtn_tmp;
      for (int i = 0 ; i < numMaxCoFs ; i++){
    	  maxCoFacet = mesh().maxCofacetLID( cellDim, cellLID[c] , i, tmp );
    	  SUNDANCE_MSG3(setupVerb(),  " maxCoFacet = " << maxCoFacet << " ComboIndex:" << maxCellFuncSetsIndex_[maxCoFacet]);
    	  SUNDANCE_MSG3(setupVerb(),  " function set = " << elemFuncSets_[maxCellFuncSetsIndex_[maxCoFacet]] );
    	  rtn_tmp.merge( elemFuncSets_[maxCellFuncSetsIndex_[maxCoFacet]] );
      }
      // make intersect
      if (c == 0){
    	  rtn = rtn_tmp;
      }
      else {
    	  rtn = rtn.intersection(rtn_tmp);
      }
    } // end loop over cells
  }
  else
  {
    rtn = elemFuncSets_[maxCellFuncSetsIndex_[cellLID[0]]];
	SUNDANCE_MSG3(setupVerb(),  "InhomogeneousDOFMapHN::allowedFuncsOnCellBatch cellDim = " << cellDim);
    for (int c=1; c<cellLID.size(); c++)
    {
  	    SUNDANCE_MSG3(setupVerb(),  " maxCoFacet merge set of cell:" << cellLID[c] << " , " << elemFuncSets_[maxCellFuncSetsIndex_[cellLID[c]]]);
    	rtn = rtn.intersection(elemFuncSets_[maxCellFuncSetsIndex_[cellLID[c]]]);
    } // end loop over cells
  }

  return rcp(new Set<int>(rtn));
}

void InhomogeneousDOFMapHN::computeOffsets(int dim, int localCount)
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
    << " back from MPI accumulate , myOffset=" << myOffset);

  if (setupVerb() > 2)
  {
    comm().synchronize();
    comm().synchronize();
    comm().synchronize();
    comm().synchronize();
  }

  for (int cellNr=0 ; cellNr < dofs_[dim].size() ; cellNr++)
  {
    for (int funcID=0 ; funcID < dofs_[dim][cellNr].size() ;  funcID++)
    {
    	if (funcDefined_[dim][funcID][cellNr] > 0){
         SUNDANCE_MSG2(setupVerb(), " cellNr = " << cellNr << " , funcID=" << funcID );
    	 for (int nrDoF=0 ; nrDoF < dofs_[dim][cellNr][funcID].size() ;  nrDoF++)
         {
        	// increment only when this is not marked as hanging
        	if (dofs_[dim][cellNr][funcID][nrDoF] >= 0) {
        		 dofs_[dim][cellNr][funcID][nrDoF] += myOffset;
        	}
         }
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


void InhomogeneousDOFMapHN::print(std::ostream& os) const
{
	// todo: implement this
	/*
  Tabs tab0;
  std::cout << tab0 << "dof map: " << std::endl;
  for (int d=0; d<=dim_; d++)
  {
    Tabs tab1;
    std::cout << tab1 << d << "-cells: " << std::endl;
    for (int n=0; n<mesh().numCells(d); n++)
    {
      Tabs tab2;
      std::cout << tab2 << "d=" << d << " cell=" << n << std::endl;
      RCP<const Set<int> > funcs = allowedFuncsOnCellBatch(d, tuple(n));
      for (Set<int>::const_iterator f=funcs->begin(); f!=funcs->end(); f++)
      {
        Tabs tab3;
        Array<int> dofs;
        getDOFsForCell(d, n, *f, dofs);
        std::cout << tab3 << " f=" << *f << " dofs=" << dofs << std::endl;
      }
    }
  }*/
}
