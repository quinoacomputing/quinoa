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

#include "SundanceHomogeneousDOFMap.hpp"
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


static Time& dofLookupTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("unbatched dof lookup"); 
  return *rtn;
}

static Time& dofBatchLookupTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("batched dof lookup"); 
  return *rtn;
}

HomogeneousDOFMap::HomogeneousDOFMap(const Mesh& mesh, 
  const BasisFamily& basis,
  int numFuncs, 
  int setupVerb)
  : DOFMapBase(mesh, setupVerb), 
    dim_(mesh.spatialDim()),
    dofs_(mesh.spatialDim()+1),
    maximalDofs_(),
    haveMaximalDofs_(false),
    localNodePtrs_(mesh.spatialDim()+1),
    nNodesPerCell_(mesh.spatialDim()+1),
    totalNNodesPerCell_(mesh.spatialDim()+1, 0),
    numFacets_(mesh.spatialDim()+1),
    originalFacetOrientation_(2),
    basisIsContinuous_(false)
{
  verbosity() = DOFMapBase::classVerbosity();
  
  CellFilter maximalCells = new MaximalCellFilter();
  cellSets().append(maximalCells.getCells(mesh));
  cellDimOnCellSets().append(mesh.spatialDim());

  allocate(mesh, basis, numFuncs);
  initMap();
}


HomogeneousDOFMap::HomogeneousDOFMap(const Mesh& mesh, 
                                     const BasisFamily& basis,
                                     const Array<CellFilter>& subregions,
                                     int numFuncs)
  : DOFMapBase(mesh), 
    dim_(mesh.spatialDim()),
    dofs_(mesh.spatialDim()+1),
    maximalDofs_(),
    haveMaximalDofs_(false),
    localNodePtrs_(mesh.spatialDim()+1),
    nNodesPerCell_(mesh.spatialDim()+1),
    totalNNodesPerCell_(mesh.spatialDim()+1, 0),
    numFacets_(mesh.spatialDim()+1),
    originalFacetOrientation_(2),
    basisIsContinuous_(false)
{
  verbosity() = DOFMapBase::classVerbosity();
  
  for (int r=0; r<subregions.size(); r++)
    {
      cellSets().append(subregions[r].getCells(mesh));
      cellDimOnCellSets().append(subregions[r].dimension(mesh));
    }

  allocate(mesh, basis, numFuncs);
  initMap();
}


void HomogeneousDOFMap::allocate(const Mesh& mesh, 
                                 const BasisFamily& basis,
                                 int numFuncs)
{
  Tabs tab;
  SUNDANCE_MSG1(setupVerb(), tab << "allocating DOF map for nFuncs=" << numFuncs);
  Array<int> fid(numFuncs);
  for (int f=0; f<numFuncs; f++) fid[f] = f;
  funcIDOnCellSets().append(fid);


  
  for (int d=0; d<=dim_; d++)
    {
      Tabs tab1;
      SUNDANCE_MSG2(setupVerb(), tab1 << "allocating d=" << d);
      /* record the number of facets for each cell type so we're
       * not making a bunch of mesh calls */
      numFacets_[d].resize(d);
      for (int fd=0; fd<d; fd++) numFacets_[d][fd]=mesh.numFacets(d, 0, fd);
      SUNDANCE_MSG3(setupVerb(), tab1 << "num facets for dimension " << d << " is " 
                         << numFacets_[d]);
          
      /* look up the node pointer for this cell and for all of its
       * facets */
      basis.ptr()->getLocalDOFs(mesh.cellType(d), localNodePtrs_[d]);


      SUNDANCE_MSG3(setupVerb(), tab1 << "node ptrs for dimension " << d << " are " 
                         << localNodePtrs_[d]);

      /* with the node pointers in hand, we can work out the number
       * of nodes per cell in this dimension */
      if (localNodePtrs_[d][d].size() > 0) 
        {
          nNodesPerCell_[d] = localNodePtrs_[d][d][0].size();
        }
      else
        {
          nNodesPerCell_[d] = 0;
        }
      SUNDANCE_MSG3(setupVerb(), tab1 << 
                         "num nodes for dimension " << d << " is " 
                         << nNodesPerCell_[d]);

      totalNNodesPerCell_[d] = nNodesPerCell_[d];
      for (int dd=0; dd<d; dd++) 
        {
          totalNNodesPerCell_[d] += numFacets_[d][dd]*nNodesPerCell_[dd];
        }

      /* we know from the mesh the number of cells in this dimension */
      if (nNodesPerCell_[d] > 0)
        {
          dofs_[d].resize(mesh.numCells(d));
        }
      else
        {
          dofs_[d].resize(0);
        }

      if (d > 0 && d < dim_) originalFacetOrientation_[d-1].resize(mesh.numCells(d));

      /* If any nodes are associated with the facets, then we know we have
       * a continuous basis function */
      if (d < dim_ && nNodesPerCell_[d] > 0) basisIsContinuous_ = true;


      /* now that we know the number of nodes per cell for this dimension,
       * we can allocate space for the DOFs in this dimension */
      int numCells = dofs_[d].size();
      for (int c=0; c<numCells; c++)
        {
          dofs_[d][c].resize(funcIDList().size() * nNodesPerCell_[d]);
          /* set everything to uninitializedVal() */
          for (int i=0; i<dofs_[d][c].size(); i++) 
            {
              dofs_[d][c][i] = uninitializedVal();
            }
        }
    }
  SUNDANCE_MSG1(setupVerb(), tab << "done allocating DOF map");
}

void HomogeneousDOFMap::initMap()
{
  Tabs tab;
  SUNDANCE_MSG1(setupVerb(), tab << "initializing DOF map");
  /* start the DOF count at zero. */
  int nextDOF = 0;

  /* Space in which to keep a list of remote cells needed by each processor
   * for each dimension. The first index is dimension, the second proc, the
   * third cell number. */
  Array<Array<Array<int> > > remoteCells(mesh().spatialDim()+1,
                                         mesh().comm().getNProc());
  
  for (int r=0; r<numCellSets(); r++)
    {
      /* Loop over maximal cells in the order specified by the cell iterator.
       * This might be reordered relative to the mesh. 
       *
       * At each maximal cell, we'll run through the facets and 
       * assign DOFs. That will take somewhat more work, but gives much 
       * better cache locality for the matrix because all the DOFs for
       * each maximal element and its facets are grouped together. */

      CellSet cells = cellSet(r);
      CellIterator iter;
      for (iter=cells.begin(); iter != cells.end(); iter++)
        {
          /* first assign any DOFs associated with the maximal cell */
          int cellLID = *iter;
          int owner;
      
          if (nNodesPerCell_[dim_] > 0)
            {
              /* if the maximal cell is owned by another processor,
               * put it in the list of cells for which we need to request 
               * DOF information from another processor */
              if (isRemote(dim_, cellLID, owner))
                {
                  int dummy=0;
                  int cellGID = mesh().mapLIDToGID(dim_, cellLID);
                  remoteCells[dim_][owner].append(cellGID); 
                  setDOFs(dim_, cellLID, dummy);
                }
              else /* the cell is locally owned, so we can 
                    * set its DOF numbers now */
                {
                  setDOFs(dim_, cellLID, nextDOF);
                }
            }

          /* Now assign any DOFs associated with the facets. */
          /* We can skip this step if the basis is discontinuous at element
           * boundaries, because the facets will own no nodes */
          if (basisIsContinuous_)
            {
              for (int d=0; d<dim_; d++)
                {
                  if (nNodesPerCell_[d] > 0)
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
                              /* the facet may be owned by another processor */
                              if (isRemote(d, facetLID[f], owner))
                                {
                                  int dummy=0;
                                  int facetGID = mesh().mapLIDToGID(d, facetLID[f]);
                                  remoteCells[d][owner].append(facetGID);
                                  setDOFs(d, facetLID[f], dummy);
                                }
                              else /* we can assign a DOF locally */
                                {
                                  /* assign DOF */
                                  setDOFs(d, facetLID[f], nextDOF);
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
        }
    }
  /* Done with first pass, in which we have assigned DOFs for all
   * local processors. We now have to share DOF information between
   * processors */

  if (mesh().comm().getNProc() > 1)
    {
      for (int d=0; d<=dim_; d++)
        {
          if (nNodesPerCell_[d] > 0)
            {
              computeOffsets(d, nextDOF);
              shareDOFs(d, remoteCells[d]);
            }
        }
    }
  else
    {
      setLowestLocalDOF(0);
      setNumLocalDOFs(nextDOF);
      setTotalNumDOFs(nextDOF);
    }
  SUNDANCE_MSG1(setupVerb(), tab << "done initializing DOF map");
}

void HomogeneousDOFMap::shareDOFs(int cellDim,
                                  const Array<Array<int> >& outgoingCellRequests)
{
  int np = mesh().comm().getNProc();
  int rank = mesh().comm().getRank();

  if (np == 1) return;

  TEUCHOS_TEST_FOR_EXCEPTION(np != outgoingCellRequests.size(), std::runtime_error,
                     "incorrect size of outgoingCellRequests array. Size is "
                     << outgoingCellRequests.size() << ", should be np=" << np);



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
      int blockSize = 1;
      if (cellDim > 0 && cellDim < dim_)
        {
          outgoingDOFs[p].resize(2 * nReq);
          blockSize = 2;
        }
      else
        {
          outgoingDOFs[p].resize(nReq);
        }
      SUNDANCE_MSG3(setupVerb(),  
                   "p=" << mesh().comm().getRank() << " recv'd from proc=" << p
                   << " reqs for DOFs for cells " << requestsFromProc);
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
          outgoingDOFs[p][blockSize*c] = dofs_[cellDim][LID][0];
          if (cellDim > 0 && cellDim < dim_) 
            {
              outgoingDOFs[p][blockSize*c+1] = originalFacetOrientation_[cellDim-1][LID];
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
      int blockSize = 1;
      if (cellDim > 0 && cellDim < dim_)
        {
          blockSize = 2;
        }

      for (int c=0; c<dofsFromProc.size()/blockSize; c++)
        {
          int cellGID = outgoingCellRequests[p][c];
          int cellLID = mesh().mapGIDToLID(cellDim, cellGID);
          int dof = dofsFromProc[blockSize*c];
          setDOFs(cellDim, cellLID, dof, true);
          if (cellDim > 0 && cellDim < dim_) 
            {
              originalFacetOrientation_[cellDim-1][cellLID] = dofsFromProc[blockSize*c+1];
            }
        }
    }
  
}



void HomogeneousDOFMap::setDOFs(int cellDim, int cellLID, 
                                int& nextDOF,
                                bool isRemote)
{
  Tabs tab;
  SUNDANCE_MSG3(setupVerb(), tab << "setting DOFs for " << cellDim << "-cell " << cellLID);
  Array<int>& cellDOFs = dofs_[cellDim][cellLID];
  
  int nn = nNodesPerCell_[cellDim];
  int nf = funcIDList().size();

  for (int i=0; i<nn; i++)
    {
      for (int f=0; f<nf; f++)
        {
          int k = nextDOF++;
          cellDOFs[funcIDList()[f] + nf*i] = k;
          if (isRemote) addGhostIndex(k);
        }
    }
}



void HomogeneousDOFMap::getDOFsForCellBatch(int cellDim, 
                                            const Array<int>& cellLID,
                                            Array<int>& dofs,
                                            int& nNodes) const 
{
  TimeMonitor timer(dofBatchLookupTimer());

  Tabs tab;
  SUNDANCE_MSG3(setupVerb(), 
               tab << "getDOFsForCellBatch(): cellDim=" << cellDim
               << " cellLID=" << cellLID);

  if (cellLID.size()==0) return;

  if (cellDim == dim_)
    {
      if (!haveMaximalDofs_) 
        {
          buildMaximalDofTable();
        }

      SUNDANCE_MSG4(setupVerb(), tab << "getting dofs for maximal cells");

      nNodes = totalNNodesPerCell_[cellDim];
      int nCells = cellLID.size();
      int nTotalCells = mesh().numCells(cellDim);
      int nf = funcIDList().size();
      dofs.resize(totalNNodesPerCell_[cellDim] * cellLID.size() * nf);
      
      for (int c=0; c<cellLID.size(); c++)
        {
          Tabs tab1;
          SUNDANCE_MSG4(setupVerb(), tab1 << "cell=" << c);
          for (int fid=0; fid<nf; fid++)
            {
              Tabs tab2;
              SUNDANCE_MSG4(setupVerb(), tab2 << "f= " << fid);
              for (int n=0; n<nNodes; n++) 
                {
                  Tabs tab3;
                  dofs[(fid*nCells + c)*nNodes + n] 
                    = maximalDofs_[(fid*nTotalCells + cellLID[c])*nNodes + n];
                  SUNDANCE_MSG4(setupVerb(), tab3 << "n=" << n << " dof=" 
                                        << dofs[(fid*nCells + c)*nNodes + n]);
                }
              
            }
        }
    }
  else
    {
      int nf = funcIDList().size();
      int nCells = cellLID.size();
      dofs.resize(totalNNodesPerCell_[cellDim] * cellLID.size() * nf);

      Tabs tab1;
      SUNDANCE_MSG4(setupVerb(), tab1 << "getting dofs for non-maximal cells");
  
      static Array<Array<int> > facetLID(3);
      static Array<Array<int> > facetOrientations(3);
      static Array<int> numFacets(3);

      for (int d=0; d<cellDim; d++) 
        {
          numFacets[d] = mesh().numFacets(cellDim, cellLID[0], d);
          mesh().getFacetLIDs(cellDim, cellLID, d, facetLID[d], 
                              facetOrientations[d]);
        }

  
      int nInteriorNodes = localNodePtrs_[cellDim][cellDim][0].size();
  
      nNodes = totalNNodesPerCell_[cellDim];
      
      for (int c=0; c<cellLID.size(); c++)
        {
          Tabs tab2;
          SUNDANCE_MSG4(setupVerb(), tab2 << "cell=" << c);
          SUNDANCE_MSG4(setupVerb(), tab2 << "doing interior nodes");

          /* first get the DOFs for the nodes associated with 
           * the cell's interior */
          if (cellDim==0 || nInteriorNodes <= 1) /* orientation-independent */
            {
              for (int m=0; m<nInteriorNodes; m++)
                {
                  Tabs tab3;
                  int ptr = localNodePtrs_[cellDim][cellDim][0][m];
                  for (int f=0; f<nf; f++)
                    {
                      dofs[(f*nCells + c)*nNodes+ptr] 
                        = dofs_[cellDim][cellLID[c]][f + nf*m];
                      SUNDANCE_MSG4(setupVerb(), tab3 << "n=" << m << " f=" << f 
                                            << " dof=" 
                                            <<  dofs[(f*nCells + c)*nNodes+ptr]);
                    }
                }
            }
          else
            {
              int sign = originalFacetOrientation_[cellDim-1][cellLID[c]];
              for (int m=0; m<nInteriorNodes; m++)
                {
                  int n = m;
                  if (sign < 0) n = nInteriorNodes - 1 - m;
                  Tabs tab3;
                  int ptr = localNodePtrs_[cellDim][cellDim][0][m];
                  for (int f=0; f<nf; f++)
                    {
                      dofs[(f*nCells + c)*nNodes+ptr] 
                        = dofs_[cellDim][cellLID[c]][f + nf*n];
                      SUNDANCE_MSG4(setupVerb(), tab3 << "n=" << m << " f=" << f 
                                            << " dof=" 
                                            <<  dofs[(f*nCells + c)*nNodes+ptr]);
                    }
                }
            }

          SUNDANCE_MSG4(setupVerb(), tab1 << "doing facet nodes");
          /* now get the DOFs for the nodes on the facets */
          for (int d=0; d<cellDim; d++)
            {
              Tabs tab3;
              SUNDANCE_MSG4(setupVerb(), tab3 << "d= " << d);
              for (int f=0; f<numFacets[d]; f++)
                {
                  Tabs tab4;
                  int facetID = facetLID[d][c*numFacets[d]+f];
                  
                  SUNDANCE_MSG4(setupVerb(), tab4 << "f= " << f << " facetLID = " << facetID);
                  int nFacetNodes = localNodePtrs_[cellDim][d][f].size();
                  if (d==0 || nFacetNodes <= 1) /* orientation-independent */
                    {
                      Tabs tab5;
                      for (int n=0; n<nFacetNodes; n++)
                        {
                          SUNDANCE_MSG3(setupVerb(), 
                                       tab5 << "n=" << n);
                          int ptr = localNodePtrs_[cellDim][d][f][n];
                          SUNDANCE_MSG3(setupVerb(), 
                                       tab5 << "local ptr=" << ptr);
                          for (int funcID=0; funcID<nf; funcID++)
                            {
                              SUNDANCE_MSG3(setupVerb(), 
                                           tab5 << "found dof=" 
                                           << dofs_[d][facetID][funcID + nf*n]);
                              dofs[(funcID*nCells + c)*nNodes+ptr] 
                                = dofs_[d][facetID][funcID + nf*n];
                            }
                        }
                    }
                  else /* orientation-dependent */
                    {
                      Tabs tab5;
                      int facetOrientation 
                        = facetOrientations[d][c*numFacets[d]+f] 
                        * originalFacetOrientation_[d-1][facetID];
                      for (int m=0; m<nFacetNodes; m++)
                        {
                          int n = m;
                          if (facetOrientation<0) n = nFacetNodes-1-m;
                          SUNDANCE_MSG3(setupVerb(), 
                                       tab5 << "n=" << n);
                          int ptr = localNodePtrs_[cellDim][d][f][m];
                          SUNDANCE_MSG3(setupVerb(), 
                                       tab5 << "local ptr=" << ptr);
                          for (int funcID=0; funcID<nf; funcID++)
                            {
                              SUNDANCE_MSG3(setupVerb(), 
                                           tab5 << "found dof=" 
                                           << dofs_[d][facetID][funcID + nf*n]);
                              dofs[(funcID*nCells + c)*nNodes+ptr] 
                                = dofs_[d][facetID][funcID + nf*n];
                            }
                        }
                    }
                }
            }
        }
    }
}    

void HomogeneousDOFMap::buildMaximalDofTable() const
{
  Tabs tab;
  int cellDim = dim_;
  int nf = funcIDList().size();
  int nCells = mesh().numCells(dim_);

  SUNDANCE_MSG2(setupVerb(), tab << "building dofs for maximal cells");

  SUNDANCE_MSG3(setupVerb(), tab << "nf=" << nf
               << " total nNodes=" << totalNNodesPerCell_[cellDim]);
  
  static Array<Array<int> > facetLID(3);
  static Array<Array<int> > facetOrientations(3);
  static Array<int> numFacets(3);

  Array<int> cellLID(nCells);

  for (int c=0; c<cellLID.size(); c++) cellLID[c]=c;
  
  for (int d=0; d<cellDim; d++) 
    {
      numFacets[d] = mesh().numFacets(cellDim, cellLID[0], d);
      mesh().getFacetLIDs(cellDim, cellLID, d, facetLID[d], facetOrientations[d]);
    }

  
  int nInteriorNodes = localNodePtrs_[cellDim][cellDim][0].size();
  

  int nNodes = totalNNodesPerCell_[cellDim];

  maximalDofs_.resize(nCells*nf*nNodes);

  for (int c=0; c<nCells; c++)
    {
      Tabs tab1;
      SUNDANCE_MSG4(setupVerb(), tab1 << "working on cell=" << c << " LID=" << cellLID[c]);
      /* first get the DOFs for the nodes associated with 
       * the cell's interior */
      SUNDANCE_MSG4(setupVerb(), tab1 << "doing interior nodes");
      for (int n=0; n<nInteriorNodes; n++)
        {
          int ptr = localNodePtrs_[cellDim][cellDim][0][n];
          for (int f=0; f<nf; f++)
            {
              maximalDofs_[(f*nCells + c)*nNodes+ptr] 
                = dofs_[cellDim][c][f + nf*n];
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
              int nFacetNodes = localNodePtrs_[cellDim][d][f].size();
              if (d == 0 || nFacetNodes <= 1) /* orientation-independent */
                {
                  for (int n=0; n<nFacetNodes; n++)
                    {
                      Tabs tab4;
                      SUNDANCE_MSG3(setupVerb(), 
                                   tab4 << "n=" << n);
                      int ptr = localNodePtrs_[cellDim][d][f][n];
                      for (int funcID=0; funcID<nf; funcID++)
                        {
                          SUNDANCE_MSG3(setupVerb(), 
                                       tab4 << "found dof=" 
                                       << dofs_[d][facetID][funcID + nf*n])
                            maximalDofs_[(funcID*nCells + c)*nNodes+ptr] 
                            = dofs_[d][facetID][funcID + nf*n];
                        }
                    }
                }
              else /* orientation-dependent */
                {
                  Tabs tab4;
                  int facetOrientation = facetOrientations[d][c*numFacets[d]+f]
                    * originalFacetOrientation_[d-1][facetID];
                  for (int m=0; m<nFacetNodes; m++)
                    {
                      Tabs tab4;
                      int n = m;
                      if (facetOrientation<0) n = nFacetNodes-1-m;
                      SUNDANCE_MSG3(setupVerb(), 
                                   tab4 << "n=" << n);
                      int ptr = localNodePtrs_[cellDim][d][f][m];
                      for (int funcID=0; funcID<nf; funcID++)
                        {
                          SUNDANCE_MSG3(setupVerb(), 
                                       tab4 << "found dof=" 
                                       << dofs_[d][facetID][funcID + nf*n])
                            maximalDofs_[(funcID*nCells + c)*nNodes+ptr] 
                            = dofs_[d][facetID][funcID + nf*n];
                        }
                    }
                }
            }
        }
    }
  haveMaximalDofs_ = true;
}





void HomogeneousDOFMap::computeOffsets(int dim, int localCount)
{
  if (verbosity() > 2)
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

  if (verbosity() > 2)
    {
      comm().synchronize();
      comm().synchronize();
      comm().synchronize();
      comm().synchronize();
    }

  for (int c=0; c<dofs_[dim].size(); c++)
    {
      if (hasBeenAssigned(dim, c))
        {
          for (int n=0; n<dofs_[dim][c].size(); n++) 
            {
              dofs_[dim][c][n] += myOffset;
            }
        }
    }

  setLowestLocalDOF(myOffset);
  setNumLocalDOFs(localCount);
  setTotalNumDOFs(totalDOFCount);

  SUNDANCE_MSG2(setupVerb(), 
               "p=" << mesh().comm().getRank() 
               << " done sharing offsets for DOF numbering for dim=" << dim);
  if (verbosity() > 2)
    {
      comm().synchronize();
      comm().synchronize();
      comm().synchronize();
      comm().synchronize();
    }

}                           



void HomogeneousDOFMap::print(std::ostream& os) const
{
  int myRank = mesh().comm().getRank();

  Tabs tabs;

  os << "DOFS = " << dofs_ << std::endl;

  for (int p=0; p<mesh().comm().getNProc(); p++)
    {
      mesh().comm().synchronize();
      mesh().comm().synchronize();
      if (p == myRank)
        {
          os << tabs << 
            "========= DOFMap on proc p=" << p << " =============" << std::endl;
          for (int d=dim_; d>=0; d--)
            {
              Tabs tabs1;
              os << tabs1 << "dimension = " << d << std::endl;
              for (int c=0; c<mesh().numCells(d); c++)
                {
                  Tabs tabs2;
                  os << tabs2 << "Cell LID=" << c << " GID=" 
                     << mesh().mapLIDToGID(d, c) << std::endl;
                  for (int f=0; f<funcIDList().size(); f++)
                    {
                      Tabs tabs3;
                      Array<int> dofs;
                      getDOFsForCell(d, c, funcIDList()[f], dofs);
                      os << tabs3 << "f=" << funcIDList()[f] << " " 
                         << dofs << std::endl;
                      if (false)
                        {
                          os << tabs3 << "{";
                          for (int i=0; i<dofs.size(); i++)
                            {
                              if (i != 0) os << ", ";
                              if (isLocalDOF(dofs[i])) os << "L";
                              else os << "R";
                            }
                          os << "}" << std::endl;
                        }
                    }
                }
            }
        }
      mesh().comm().synchronize();
      mesh().comm().synchronize();
    }
}


