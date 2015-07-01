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
#include "SundanceDOFMapBase.hpp"
#include "PlayaMPIContainerComm.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Teuchos;
using Playa::MPIComm;
using Playa::MPIContainerComm;




DOFMapBase::DOFMapBase(const Mesh& mesh, int setupVerb)
  : setupVerb_(setupVerb),
    localProcID_(mesh.comm().getRank()),
    mesh_(mesh),
    lowestLocalDOF_(),
    numLocalDOFs_(),
    numDOFs_(),
    ghostIndices_(rcp(new Array<int>()))
{}


void DOFMapBase::getDOFsForCell(int cellDim, int cellLID,
  int funcID,
  Array<int>& dofs) const
{
  TimeMonitor timer(dofLookupTimer());
  
  Array<Array<int> > allDofs;
  Array<int> nNodes;
  RCP<const MapStructure> s 
    = getDOFsForCellBatch(cellDim, tuple(cellLID), makeSet(funcID), allDofs, nNodes,0);


  int chunkNumber = s->chunkForFuncID(funcID);
  int funcIndex = s->indexForFuncID(funcID);
  dofs.resize(nNodes[chunkNumber]);
  for (int i=0; i<nNodes[chunkNumber]; i++)
  {
    dofs[i] = allDofs[chunkNumber][nNodes[chunkNumber]*funcIndex + i];
  }
}

Time& DOFMapBase::dofLookupTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("unbatched dof lookup"); 
  return *rtn;
}

Time& DOFMapBase::batchedDofLookupTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("batched dof lookup"); 
  return *rtn;
}

