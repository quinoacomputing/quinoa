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

#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceLocalDOFMap.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


using std::endl;
using std::setw;


LocalDOFMap::LocalDOFMap(int numBlocks, int verb)
  : verb_(verb),
    isUsed_(numBlocks),
    hasCells_(false),
    nLocalNodesPerChunk_(rcp(new Array<Array<int> >(numBlocks))),
    mapStruct_(rcp(new Array<RCP<const MapStructure> >(numBlocks))),
    localDOFs_(rcp(new Array<Array<Array<int> > >(numBlocks))),
    cellLID_(),
    activeCellDim_(-1),
    maxCellDim_(-1)
{}

int LocalDOFMap::nCells() const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!hasCells(), std::runtime_error,
    "cells not valid when LocalDOFMap::nCells() called");
  return cellLID_->size();
}

void  LocalDOFMap::markAsUnused() 
{
  for (int b=0; b<numBlocks(); b++) 
  {
    isUsed_[b] = 0;
  }
  hasCells_ = false;
}


bool  LocalDOFMap::isUnused() const 
{
  for (int b=0; b<numBlocks(); b++) 
  {
    if (isUsed_[b]) return false;
  }
  return true;
}

void LocalDOFMap::verifyValidBlock(int b) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(b < 0 || b>=numBlocks(), std::runtime_error,
    "block index " << b << " out of range [0," << numBlocks() << ")");
}

void LocalDOFMap::setCells(
  int activeCellDim, 
  int maxCellDim,
  const RCP<const Array<int> >& cellLID)
{
  activeCellDim_ = activeCellDim;
  maxCellDim_ = maxCellDim_;
  cellLID_ = cellLID;
  hasCells_=true;
}


void LocalDOFMap::fillBlock(int b, const RCP<DOFMapBase>& globalMap,
  const Array<Set<int> >& requiredFuncs)
{
  Tabs tab;
  verifyValidBlock(b);

  SUNDANCE_MSG2(verb_, tab << "getting local DOFs for block=" << b 
    << ", active cell dim="
    << activeCellDim_);

  mapStruct(b) = globalMap->getDOFsForCellBatch(
    activeCellDim_,
    *cellLID_,
    requiredFuncs[b],
    localDOFs(b),
    nLocalNodesPerChunk(b),
    verb_);
}

std::ostream& LocalDOFMap::print(std::ostream& os) const
{
  Tabs tab0;
  bool noneAreUsed=true;
  for (int b=0;  b<numBlocks(); b++) 
  {
    if (isUsed(b)) {noneAreUsed = false; break;}
  }
  if (noneAreUsed)
  {
    os << std::endl << tab0 << "[empty local DOF map]";
  }
  else
  {
    os << std::endl << tab0 << "LocalDOFMap[" << std::endl;
    Tabs tab;
    os << tab << "num cells=" << cellLID_->size() << ", active cellDim="
       << activeCellDim_ << ", maxCellDim=" << maxCellDim_ << std::endl;
    os << tab << "num blocks=" << numBlocks() << std::endl << std::endl;

    for (int b=0; b<numBlocks(); b++)
    {
      Tabs tab1;
      os << tab1 << "block " << b << " of " << numBlocks()
         << *mapStruct(b) << std::endl;
      const Array<Array<int> >& dofs = localDOFs(b);
      int nChunks = mapStruct(b)->numBasisChunks();
      
      for (int c=0; c<cellLID_->size(); c++)
      {
        Tabs tab2;
        os << tab2 << "cell LID=" << (*cellLID_)[c] << std::endl;
        for (int chunk=0; chunk<nChunks; chunk++)
        {
          Tabs tab3;
          int nFuncs = mapStruct(b)->numFuncs(chunk);
          const Array<int>& funcs = mapStruct(b)->funcs(chunk);
          int nNodes = nLocalNodesPerChunk(b)[chunk];
          for (int f=0; f<nFuncs; f++)
          {
            os << tab3 << "fid=" << funcs[f] << " dofs=";
            int funcOffset = mapStruct(b)->indexForFuncID(funcs[f]);
            for (int n=0; n<nNodes; n++)
            {
              int dof = dofs[chunk][(c*nFuncs+funcOffset)*nNodes+n];
              os << setw(6) << dof;
              if (n < nNodes-1) os << ", ";
            }
            os << std::endl;
          }
        }
      }
      
    }
    os << tab0 << "] ### End LocalDOFMap" << std::endl;
  }
  return os;
}
