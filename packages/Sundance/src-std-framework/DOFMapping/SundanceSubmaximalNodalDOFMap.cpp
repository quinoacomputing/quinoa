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
#include "SundanceSubmaximalNodalDOFMap.hpp"
#include "SundanceLagrange.hpp"
#include "PlayaMPIContainerComm.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Teuchos;
using Playa::MPIComm;
using Playa::MPIContainerComm;


SubmaximalNodalDOFMap
::SubmaximalNodalDOFMap(const Mesh& mesh, 
  const CellFilter& cf,
  int nFuncs,
  int setupVerb)
  : DOFMapBase(mesh, setupVerb),
    dim_(0),
    nTotalFuncs_(nFuncs),
    domain_(cf),
    domains_(tuple(cf)),
    nodeLIDs_(),
    nodeDOFs_(),
    lidToPtrMap_(),
    mapStructure_()
{
  Tabs tab0(0);
  SUNDANCE_MSG1(setupVerb, tab0 << "in SubmaximalNodalDOFMap ctor");
  Tabs tab1;
  SUNDANCE_MSG2(setupVerb, tab1 << "domain " << domain_);
  SUNDANCE_MSG2(setupVerb, tab1 << "N funcs " << nFuncs);

  const MPIComm& comm = mesh.comm();
  int rank = comm.getRank();
  int nProc = comm.getNProc();
  
  dim_ = cf.dimension(mesh);  
  TEUCHOS_TEST_FOR_EXCEPT(dim_ != 0);

  CellSet nodes = cf.getCells(mesh);
  int nc = nodes.numCells();
  nodeLIDs_.reserve(nc);
  nodeDOFs_.reserve(nc);

  Array<Array<int> > remoteNodes(nProc);
  
  int nextDOF = 0;
  int k=0; 
  for (CellIterator c=nodes.begin(); c!=nodes.end(); c++, k++)
  {
    int nodeLID = *c;
    lidToPtrMap_.put(nodeLID, k);
    nodeLIDs_.append(nodeLID);
    int remoteOwner = rank;
    if (isRemote(0, nodeLID, remoteOwner))
    {
      int GID = mesh.mapLIDToGID(0, nodeLID);
      remoteNodes[remoteOwner].append(GID);
      for (int f=0; f<nFuncs; f++) nodeDOFs_.append(-1);
    }
    else
    {
      for (int f=0; f<nFuncs; f++) nodeDOFs_.append(nextDOF++);
    }
  }

  /* Compute offsets for each processor */
  int localCount = nextDOF;
  computeOffsets(localCount);
  
  /* Resolve remote DOF numbers */
  shareRemoteDOFs(remoteNodes);

  BasisFamily basis = new Lagrange(1);
  mapStructure_ = rcp(new MapStructure(nTotalFuncs_, basis.ptr()));
}

void SubmaximalNodalDOFMap::computeOffsets(int localCount)
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
    
    int nDofs = nodeDOFs_.size();
    for (int i=0; i<nDofs; i++)
    {
      if (nodeDOFs_[i] >= 0) nodeDOFs_[i] += myOffset;
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


void SubmaximalNodalDOFMap::shareRemoteDOFs(
  const Array<Array<int> >& outgoingNodeRequests)
{
  int np = mesh().comm().getNProc();
  if (np==1) return;
  int rank = mesh().comm().getRank();

  int nFuncs = nTotalFuncs_;

  Array<Array<int> > incomingNodeRequests;
  Array<Array<int> > outgoingDOFs(np);
  Array<Array<int> > incomingDOFs;

  SUNDANCE_MSG2(setupVerb(),  
    "p=" << mesh().comm().getRank()
    << "synchronizing DOFs for cells of dimension 0");
  SUNDANCE_MSG2(setupVerb(),  
    "p=" << mesh().comm().getRank()
    << " sending cell reqs d=0, GID=" 
    << outgoingNodeRequests);

  /* share the cell requests */
  MPIContainerComm<int>::allToAll(outgoingNodeRequests, 
    incomingNodeRequests,
    mesh().comm());
  
  /* get DOF numbers for the zeroth function index on every node that's been 
   * requested by someone else */
  for (int p=0; p<np; p++)
  {
    if (p==rank) continue;
    const Array<int>& requestsFromProc = incomingNodeRequests[p];
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
        << nodeDOFs_[LID*nFuncs]);
      outgoingDOFs[p][c] = nodeDOFs_[LID*nFuncs];
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
      int cellGID = outgoingNodeRequests[p][c];
      int cellLID = mesh().mapGIDToLID(0, cellGID);
      int dof = dofsFromProc[c];
      for (int i=0; i<nFuncs; i++)
      {
        nodeDOFs_[cellLID*nFuncs + i] = dof+i;
        addGhostIndex(dof+i);
      }
    }
  }
}



RCP<const Set<int> >
SubmaximalNodalDOFMap::allowedFuncsOnCellBatch(int cellDim,
  const Array<int>& cellLID) const 
{
  Set<int> rtn;
  for (int f=0; f<nTotalFuncs_; f++) rtn.put(f);
  return rcp(new Set<int>(rtn));
}


RCP<const MapStructure> 
SubmaximalNodalDOFMap::getDOFsForCellBatch(int cellDim,
  const Array<int>& cellLID,
  const Set<int>& requestedFuncSet,
  Array<Array<int> >& dofs,
  Array<int>& nNodes,
  int verb) const 
{
  TimeMonitor timer(batchedDofLookupTimer());
  Tabs tab0;

  SUNDANCE_MSG2(verb, tab0 << "in SubmaximalNodalDOFMap::getDOFsForCellBatch()");

  TEUCHOS_TEST_FOR_EXCEPT(cellDim != 0);

  dofs.resize(1);
  nNodes.resize(1);
  nNodes[0] = 1;

  dofs[0].reserve(cellLID.size()*nTotalFuncs_);

  for (int c=0; c<cellLID.size(); c++)
  {
    int lid = cellLID[c];
    int ptr = lidToPtrMap_.get(lid);
    for (int f=0; f<nTotalFuncs_; f++)
    {
      int q = ptr*nTotalFuncs_;
      dofs[0].append(nodeDOFs_[q + f]);
    }
  }

  return mapStructure_;
}



void SubmaximalNodalDOFMap::print(std::ostream& os) const
{
  Tabs tab0(0);
  Out::os() << tab0 << "submaximal nodal dof map: " << std::endl;
  Tabs tab1;
  Out::os() << tab1 << "0-cells: " << std::endl;
  for (int c=0; c<nodeLIDs_.size(); c++)
  {
    Tabs tab2;
    int gid = mesh().mapLIDToGID(0, nodeLIDs_[c]);
    Out::os() << tab2 << "G=" << gid << " dofs={";
    for (int f=0; f<nTotalFuncs_; f++)
    {
      if (f != 0) Out::os() << ", ";
      Out::os() << nodeDOFs_[c*nTotalFuncs_ + f];
    }
    Out::os() << "}" << std::endl;
  }
}
