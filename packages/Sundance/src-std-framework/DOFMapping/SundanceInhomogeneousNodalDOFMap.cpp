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
#include "SundanceInhomogeneousNodalDOFMap.hpp"
#include "SundanceLagrange.hpp"
#include "PlayaMPIContainerComm.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Teuchos;
using Playa::MPIComm;
using Playa::MPIContainerComm;


InhomogeneousNodalDOFMap
::InhomogeneousNodalDOFMap(const Mesh& mesh, 
  const Array<Map<Set<int>, CellFilter> >& funcSetToDomainMap,
  int setupVerb)
  : DOFMapBase(mesh, setupVerb),
    dim_(mesh.spatialDim()),
    basis_(rcp(new Lagrange(1))),
    nTotalFuncs_(),
    funcDomains_(),
    nodeDofs_(),
    elemDofs_(),
    nodeToFuncSetIndexMap_(mesh.numCells(0)),
    elemToFuncSetIndexMap_(mesh.numCells(dim_)),
    elemFuncSets_(),
    nodalFuncSets_(),
    nodeToOffsetMap_(mesh.numCells(0)),
    elemToOffsetMap_(mesh.numCells(0)),
    funcIndexWithinNodeFuncSet_(),
    elemStructure_(),
    nodeStructure_()
{
  SUNDANCE_MSG1(setupVerb, "in InhomogeneousNodalDOFMap ctor");
  SUNDANCE_MSG4(setupVerb, "func set to domain map " << funcSetToDomainMap);

  /* count the total number of functions across all subdomains */
  Set<int> allFuncs;
  for (int d=dim_; d>=0; d--)
  {
    for (Map<Set<int>, CellFilter>::const_iterator 
           i=funcSetToDomainMap[d].begin(); i!=funcSetToDomainMap[d].end(); i++)
    {
      allFuncs.merge(i->first);
    }
  }

  
  nTotalFuncs_ = allFuncs.size();
  SUNDANCE_MSG2(setupVerb, "found " << nTotalFuncs_ << " functions");
  

  /* get flat arrays of subdomains and function arrays */
  Array<CellFilter> subdomains;
  Array<CellFilter> maxSubdomains;
  Array<Array<int> > funcArrays;
  Array<Set<int> > funcSets;

  for (int d=dim_; d>=0; d--)
  {
    for (Map<Set<int>, CellFilter>::const_iterator 
           i=funcSetToDomainMap[d].begin(); i!=funcSetToDomainMap[d].end(); i++)
    {
      subdomains.append(i->second);
      if (d==dim_) 
      {
        elemFuncSets_.append(i->first);
        maxSubdomains.append(i->second);
      }

      Array<int> myFuncs = i->first.elements();
      funcArrays.append(myFuncs);
      funcSets.append(i->first);
    }
  }


  /* determine the functions appearing at each node */



  Array<Array<int> > cellLIDs(subdomains.size());
  Array<Array<int> > facetLIDs(subdomains.size());
  Array<Array<int> > facetOrientations(subdomains.size());
  Array<Set<int> > nodeToFuncSetMap(mesh.numCells(0));

  for (int r=0; r<subdomains.size(); r++)
  {
    int d = subdomains[r].dimension(mesh);
    CellSet cells = subdomains[r].getCells(mesh);
    SUNDANCE_MSG4(setupVerb, "domain " << subdomains[r] << " has functions "
      << funcSets[r]);
      
    for (CellIterator c=cells.begin(); c!=cells.end(); c++)
    {
      int cellLID = *c;
      cellLIDs[r].append(cellLID);
    }
    if (cellLIDs[r].size()==0) continue; 
    int nFacets = mesh.numFacets(d, cellLIDs[r][0], 0);
    if (d==0)
    {
      for (int c=0; c<cellLIDs[r].size(); c++)
      {
        int fLID = cellLIDs[r][c];
        nodeToFuncSetMap[fLID].merge(funcSets[r]);
      }
    }
    else
    {
      Array<int>& facetLID = facetLIDs[r];
      facetLID.resize(nFacets*cellLIDs[r].size());
      facetOrientations[r].resize(nFacets*cellLIDs[r].size());
      mesh.getFacetLIDs(d, cellLIDs[r], 0, facetLID, facetOrientations[r]);
      for (int c=0; c<cellLIDs[r].size(); c++)
      {
        for (int f=0; f<nFacets; f++)
        {
          int fLID = facetLID[c*nFacets+f];
          nodeToFuncSetMap[fLID].merge(funcSets[r]);
        }
      }
    }
  }

  /* find the distinct combinations of functions at the nodes */

  Map<Set<int>, int> fsToFSIndexMap;
  Array<int> funcSetNodeCounts;
  int nNodes = mesh.numCells(0);

  for (int n=0; n<nNodes; n++)
  {
    const Set<int>& f = nodeToFuncSetMap[n];
    int funcComboIndex;
    if (!fsToFSIndexMap.containsKey(f)) 
    {
      funcComboIndex = nodalFuncSets_.size();
      fsToFSIndexMap.put(f, funcComboIndex);
      nodalFuncSets_.append(f);
      funcSetNodeCounts.append(1);
      nodeToOffsetMap_[n] = 0;
    }
    else
    {
      funcComboIndex = fsToFSIndexMap.get(f);
      nodeToOffsetMap_[n] = f.size() * funcSetNodeCounts[funcComboIndex];
      funcSetNodeCounts[funcComboIndex]++;
    }
    nodeToFuncSetIndexMap_[n] = funcComboIndex;
  }

  SUNDANCE_MSG2(setupVerb, "nodal func sets = " << nodalFuncSets_);

  nodeDofs_.resize(nodalFuncSets_.size());
  nodeStructure_.resize(nodalFuncSets_.size());
  funcIndexWithinNodeFuncSet_.resize(nodalFuncSets_.size());

  for (int f=0; f<nodalFuncSets_.size(); f++)
  {
    Array<int> funcs = nodalFuncSets_[f].elements();
    int nFuncs = funcs.size();
    funcIndexWithinNodeFuncSet_[f].resize(nTotalFuncs_, -1);
    for (int i=0; i<nFuncs; i++)
    {
      funcIndexWithinNodeFuncSet_[f][funcs[i]] = i;
    }
    nodeDofs_[f].resize(nFuncs * funcSetNodeCounts[f], -1);
    Array<RCP<BasisDOFTopologyBase> > nodeBases = tuple(basis_);
    Array<Array<int> > nodeFuncs = tuple(nodalFuncSets_[f].elements());
    nodeStructure_[f] = rcp(new MapStructure(nTotalFuncs_, 
        nodeBases, nodeFuncs));
  }



  /* second pass to assign node DOFs */
  int nextDOF = 0;
  Array<int> hasProcessedNode(nNodes);
  Array<Array<int> > remoteNodes(mesh.comm().getNProc());

  for (int r=0; r<subdomains.size(); r++)
  {
    int d = subdomains[r].dimension(mesh);

    if (cellLIDs[r].size()==0) continue; 

    if (d==0)
    {
      for (int c=0; c<cellLIDs[r].size(); c++)
      {
        int fLID = cellLIDs[r][c];
        int funcSetIndex = nodeToFuncSetIndexMap_[fLID];
        int nFuncs = nodalFuncSets_[funcSetIndex].size();
        int dofOffset = nodeToOffsetMap_[fLID];
        if (!hasProcessedNode[fLID])
        {
          assignNode(fLID, funcSetIndex, dofOffset, nFuncs,
            remoteNodes, hasProcessedNode, nextDOF);
        }
      }
    }
    else
    {
      Array<int>& facetLID = facetLIDs[r];
      int nFacets = mesh.numFacets(d, cellLIDs[r][0], 0);      
      for (int c=0; c<cellLIDs[r].size(); c++)
      {
        for (int f=0; f<nFacets; f++)
        {
          int fLID = facetLID[c*nFacets+f];
          int funcSetIndex = nodeToFuncSetIndexMap_[fLID];
          int dofOffset = nodeToOffsetMap_[fLID];
          int nFuncs = nodalFuncSets_[funcSetIndex].size();
          if (!hasProcessedNode[fLID])
          {
            assignNode(fLID, funcSetIndex, dofOffset, nFuncs,
              remoteNodes, hasProcessedNode, nextDOF);
          }
        }
      }
    }
  }


  /* Compute offsets for each processor */
  int localCount = nextDOF;
  computeOffsets(localCount);
  
  /* Resolve remote DOF numbers */
  shareRemoteDOFs(remoteNodes);

  /* Build dof tables for elements */
  elemDofs_.resize(maxSubdomains.size());
  int count = 0;
  Array<Array<int> > elemFuncArrays;
  
  funcDomains_.resize(nTotalFuncs_);

  for (int r=0; r<subdomains.size(); r++)
  {
    int d = subdomains[r].dimension(mesh);
    if (d != dim_) continue;
    if (cellLIDs[r].size()==0) continue; 
    int nFacets = mesh.numFacets(d, cellLIDs[r][0], 0);      
    Array<int>& facetLID = facetLIDs[r];
    Array<Array<int> > dofs(1);
    dofs[0].resize(funcArrays[r].size() * cellLIDs[r].size() * nFacets);
    getFunctionDofs(d, cellLIDs[r], facetLID, funcArrays[r], dofs);
    elemDofs_[count] = dofs[0];
    elemFuncArrays.append(funcArrays[r]);
    Array<RCP<BasisDOFTopologyBase> > elemBases = tuple(basis_);
    Array<Array<int> > elemFuncs = tuple(elemFuncSets_[count].elements());
    elemStructure_.append(rcp(new MapStructure(nTotalFuncs_, 
          elemBases, elemFuncs)));
    for (int c=0; c<cellLIDs[r].size(); c++)
    {
      elemToFuncSetIndexMap_[cellLIDs[r][c]] = count;
    }
    for (int i=0; i<elemFuncs[0].size(); i++)
    {
      funcDomains_[elemFuncs[0][i]] = funcDomains_[elemFuncs[0][i]] + subdomains[r];
    }
    count++;
  }
  
}


void InhomogeneousNodalDOFMap::getFunctionDofs(int cellDim,
  const Array<int>& cellLID,
  const Array<int>& facetLID,
  const Array<int>& funcs,
  Array<Array<int> >& dofs) const
{
  Array<int>& dofChunk = dofs[0];
  dofChunk.clear();

  if (cellDim != 0)
  {
    int nFacets = mesh().numFacets(cellDim, cellLID[0], 0);
    dofChunk.resize(funcs.size() * cellLID.size() * nFacets, -1);
      
      
    for (int c=0; c<cellLID.size(); c++)
    {
      for (int f=0; f<nFacets; f++)
      {
        int fLID = facetLID[c*nFacets + f];
        int fci = nodeToFuncSetIndexMap_[fLID];
        int nodeOffset = nodeToOffsetMap_[fLID];
        for (int i=0; i<funcs.size(); i++)
        {
          int funcIndex = funcIndexWithinNodeFuncSet_[fci][funcs[i]];
          if (funcIndex >= 0)
          {
            dofChunk[(c*funcs.size()+i)*nFacets + f] 
            = nodeDofs_[fci][nodeOffset+funcIndex];
        
          }
        }
      }
    }
  }
  else
  {
    dofChunk.resize(funcs.size() * cellLID.size(), -1);

    for (int c=0; c<cellLID.size(); c++)
    {
      int fci = nodeToFuncSetIndexMap_[cellLID[c]];
      int nodeOffset = nodeToOffsetMap_[cellLID[c]];
      for (int i=0; i<funcs.size(); i++)
      {
        int funcIndex = funcIndexWithinNodeFuncSet_[fci][funcs[i]];
        if (funcIndex >= 0)
        {
          dofChunk[c*funcs.size()+ i] 
            = nodeDofs_[fci][nodeOffset+funcIndex];
        }
      }
    }
  }
}

void InhomogeneousNodalDOFMap::assignNode(int fLID,
  int funcSetIndex,
  int dofOffset,
  int nFuncs,
  Array<Array<int> >& remoteNodes,
  Array<int>& hasProcessedCell,
  int& nextDOF)
{
  /* the facet may be owned by another processor */
  int owner;
  if (isRemote(0, fLID, owner))
  {
    int facetGID 
      = mesh().mapLIDToGID(0, fLID);
    remoteNodes[owner].append(facetGID);
      
  }
  else /* we can assign a DOF locally */
  {
    /* assign DOFs */
    Array<int>& myDofs = nodeDofs_[funcSetIndex];
    for (int i=0; i<nFuncs; i++)
    {
      myDofs[dofOffset + i] = nextDOF;
      nextDOF++;
    }
  }
  hasProcessedCell[fLID] = 1;
}

void InhomogeneousNodalDOFMap::computeOffsets(int localCount)
{
  TEUCHOS_TEST_FOR_EXCEPTION(mesh().comm().getNProc() != 1,
    std::runtime_error,
    "parallel inhomogeneous DOF maps not yet supported");
  
  int totalDOFCount = localCount;
  int myOffset = 0;
  setLowestLocalDOF(myOffset);
  setNumLocalDOFs(localCount);
  setTotalNumDOFs(totalDOFCount);
}

void InhomogeneousNodalDOFMap::shareRemoteDOFs(const Array<Array<int> >& remoteNodes)
{
  TEUCHOS_TEST_FOR_EXCEPTION(mesh().comm().getNProc() != 1,
    std::runtime_error,
    "parallel inhomogeneous DOF maps not yet supported");
}

RCP<const Set<int> >
InhomogeneousNodalDOFMap::allowedFuncsOnCellBatch(int cellDim,
  const Array<int>& cellLID) const 
{
  Set<int> rtn;

  if (cellDim==0)
  {
    rtn = nodalFuncSets_[nodeToFuncSetIndexMap_[cellLID[0]]];
    for (int c=0; c<cellLID.size(); c++) 
    {
      const Set<int>& s = nodalFuncSets_[nodeToFuncSetIndexMap_[cellLID[c]]];
      rtn = rtn.intersection(s);
    }
  }
  else if (cellDim==dim_)
  {
    rtn = elemFuncSets_[elemToFuncSetIndexMap_[cellLID[0]]];
    for (int c=0; c<cellLID.size(); c++) 
    {
      const Set<int>& s = elemFuncSets_[elemToFuncSetIndexMap_[cellLID[c]]];
      rtn = rtn.intersection(s);
    }
  }
  else
  {
    Array<int> facetLID;
    Array<int> facetOrientations;
    mesh().getFacetLIDs(cellDim, cellLID, 0, facetLID, facetOrientations);
    rtn = nodalFuncSets_[nodeToFuncSetIndexMap_[facetLID[0]]];
    for (int f=0; f<facetLID.size(); f++)
    {
      const Set<int>& s = nodalFuncSets_[nodeToFuncSetIndexMap_[facetLID[f]]];
      rtn = rtn.intersection(s);
    }
  }
  return rcp(new Set<int>(rtn));
}


RCP<const MapStructure> 
InhomogeneousNodalDOFMap::getDOFsForCellBatch(int cellDim,
  const Array<int>& cellLID,
  const Set<int>& requestedFuncSet,
  Array<Array<int> >& dofs,
  Array<int>& nNodes,
  int verb) const 
{
  TimeMonitor timer(batchedDofLookupTimer());
  Tabs tab0;

  SUNDANCE_MSG2(verb, tab0 << "in InhomNodalDOFMap::getDOFsForCellBatch()");

  if (cellDim==0)
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab1 << "cell dim = " << cellDim);
    bool isHomogeneous = true;
    int firstFuncSet = nodeToFuncSetIndexMap_[cellLID[0]];
    for (int c=0; c<cellLID.size(); c++)
    {
      if (nodeToFuncSetIndexMap_[cellLID[c]] != firstFuncSet) 
      {
        isHomogeneous = false;
        break;
      }
    }

    dofs.resize(1);
    nNodes.resize(1);
    nNodes[0] = 1;
    Array<int> dummyFacets;

    if (isHomogeneous)
    {
      const Set<int>& funcSet = nodalFuncSets_[firstFuncSet];
      TEUCHOS_TEST_FOR_EXCEPT(requestedFuncSet.setDifference(funcSet).size() != 0);
      Array<int> funcs = funcSet.elements();
      getFunctionDofs(cellDim, cellLID, dummyFacets, funcs, dofs);
      return nodeStructure_[firstFuncSet];
    }
    else
    {
      Array<int> funcs = requestedFuncSet.elements();
      getFunctionDofs(cellDim, cellLID, dummyFacets, funcs, dofs);
      RCP<const MapStructure> rtn 
        = rcp(new MapStructure(nTotalFuncs_, basis_, tuple(funcs)));
      return rtn;
    }
  }
  else if (cellDim==dim_)
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab1 << "cell dim = " << cellDim);
    bool isHomogeneous = true;
    int firstFuncSet = elemToFuncSetIndexMap_[cellLID[0]];
    SUNDANCE_MSG2(verb, tab1 << "first func set = " << firstFuncSet);
    for (int c=0; c<cellLID.size(); c++)
    {
      if (elemToFuncSetIndexMap_[cellLID[c]] != firstFuncSet) 
      {
        isHomogeneous = false;
        break;
      }
    }

    dofs.resize(1);
    nNodes.resize(1);
    nNodes[0] = mesh().numFacets(cellDim, cellLID[0], 0);

    if (isHomogeneous)
    {
      Tabs tab2;
      Array<int> facetLID;
      Array<int> facetOrientations;
      mesh().getFacetLIDs(cellDim, cellLID, 0, facetLID, facetOrientations);
      if (verb >= 2)
      {
        Out::os() << tab2 << "cellLID = " << cellLID << std::endl;
        Out::os() << tab2 << "facetLID = " << facetLID << std::endl;
        Out::os() << tab2 << "elem func sets = " << elemFuncSets_ << std::endl;
      }
      const Set<int>& funcSet = elemFuncSets_[firstFuncSet];
      TEUCHOS_TEST_FOR_EXCEPT(requestedFuncSet.setDifference(funcSet).size() != 0);
      Array<int> funcs = funcSet.elements();
      getFunctionDofs(cellDim, cellLID, facetLID, funcs, dofs);
      return elemStructure_[firstFuncSet];
    }
    else
    {
      Array<int> facetLID;
      Array<int> facetOrientations;
      mesh().getFacetLIDs(cellDim, cellLID, 0, facetLID, facetOrientations);
      Array<int> funcs = requestedFuncSet.elements();
      getFunctionDofs(cellDim, cellLID, facetLID, funcs, dofs);
      RCP<const MapStructure> rtn 
        = rcp(new MapStructure(nTotalFuncs_, basis_, tuple(funcs)));
      return rtn;
    }
  }
  else
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab1 << "cell dim = " << cellDim);
    Array<int> facetLID;
    Array<int> facetOrientations;
    mesh().getFacetLIDs(cellDim, cellLID, 0, facetLID, facetOrientations);
    dofs.resize(1);
    nNodes.resize(1);
    nNodes[0] = mesh().numFacets(cellDim, cellLID[0], 0);

    Array<int> funcs = requestedFuncSet.elements();

    getFunctionDofs(cellDim, cellLID, facetLID, funcs, dofs);
    RCP<const MapStructure> rtn 
      = rcp(new MapStructure(nTotalFuncs_, basis_, tuple(funcs)));
    return rtn;
  }
}



void InhomogeneousNodalDOFMap::print(std::ostream& os) const
{
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
  }
}
