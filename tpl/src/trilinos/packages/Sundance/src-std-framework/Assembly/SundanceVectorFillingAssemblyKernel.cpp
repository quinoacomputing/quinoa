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
#include "SundanceVectorFillingAssemblyKernel.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#endif

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace Playa;
using std::setw;
using std::endl;
      
static Time& vecInsertTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("vector insertion"); 
  return *rtn;
}

VectorFillingAssemblyKernel::VectorFillingAssemblyKernel(
  const Array<RCP<DOFMapBase> >& dofMap,
  const Array<RCP<Array<int> > >& isBCIndex,
  const Array<int>& lowestLocalIndex,
  Array<Vector<double> >& b,
  bool partitionBCs,
  int verbosity
  )
  : AssemblyKernelBase(verbosity),
    b_(b),
    vec_(b.size()),
    mapBundle_(dofMap, isBCIndex, lowestLocalIndex, partitionBCs, verbosity)
{
  Tabs tab0;

  SUNDANCE_MSG1(verb(), tab0 << "VectorFillingAssemblyKernel ctor");
  
  int numBlocks = dofMap.size();

  for (int i=0; i<b_.size(); i++)
  {
    vec_[i].resize(numBlocks);
    for (int block=0; block<numBlocks; block++)
    {
      Tabs tab1;
      SUNDANCE_MSG1(verb(), tab1 << "getting vector for block b=" 
        << block << " of " << numBlocks);
      Vector<double> vecBlock = b[i].getBlock(block);
      vec_[i][block] = rcp_dynamic_cast<LoadableVector<double> >(vecBlock.ptr());

      TEUCHOS_TEST_FOR_EXCEPTION(vec_[i][block].get()==0, std::runtime_error,
        "vector block " << block << " is not loadable");
      vecBlock.zero();
    }
  }
  SUNDANCE_MSG1(verb(), tab0 << "done VectorFillingAssemblyKernel ctor");
}


void VectorFillingAssemblyKernel::buildLocalDOFMaps(
  const RCP<StdFwkEvalMediator>& mediator,
  IntegrationCellSpecifier intCellSpec,
  const Array<Set<int> >& requiredFuncs) 
{
  mapBundle_.buildLocalDOFMaps(mediator, intCellSpec, requiredFuncs,
    verb());
}


void VectorFillingAssemblyKernel::insertLocalVectorBatch(
  bool isBCRqc,
  bool useCofacetCells,
  const Array<int>& funcID,  
  const Array<int>& funcBlock, 
  const Array<int>& mvIndices, 
  const Array<double>& localValues) const
{
  TimeMonitor timer(vecInsertTimer());
  Tabs tab0;

  SUNDANCE_MSG1(verb(), tab0 << "inserting local vector batch");
  SUNDANCE_MSG4(verb(), tab0 << "vector values are " << localValues);

  const MapBundle& mb = mapBundle_;
  int nCells = mb.nCells();

  for (int i=0; i<funcID.size(); i++)
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb(), tab1 << "function ID = "<< funcID[i] 
      << " of " << funcID.size());
    SUNDANCE_MSG2(verb(), tab1 << "is BC eqn = " << isBCRqc);
    SUNDANCE_MSG2(verb(), tab1 << "num cells = " << nCells);
    SUNDANCE_MSG2(verb(), tab1 << "using cofacet cells = " << useCofacetCells);
    SUNDANCE_MSG2(verb(), tab1 << "multivector index = " 
      << mvIndices[i]);

    /* First, find the block associated with the current function
     * so that we can find the appropriate DOF information */
    int block = funcBlock[i];

    const RCP<DOFMapBase>& dofMap = mb.dofMap(block);
    int lowestLocalRow = mb.lowestLocalIndex(block);

    int chunk = mb.mapStruct(block, useCofacetCells)->chunkForFuncID(funcID[i]);
    SUNDANCE_MSG2(verb(), tab1 << "chunk = " << chunk);

    int funcIndex = mb.mapStruct(block, useCofacetCells)->indexForFuncID(funcID[i]);
    SUNDANCE_MSG2(verb(), tab1 << "func offset into local DOF map = " 
      << funcIndex);

    const Array<int>& dofs = mb.localDOFs(block, useCofacetCells, chunk);
    SUNDANCE_MSG4(verb(), tab1 << "local dofs = " << dofs);    

    int nFuncs = mb.mapStruct(block, useCofacetCells)->numFuncs(chunk);
    SUNDANCE_MSG2(verb(), tab1 << "num funcs in chunk = " << nFuncs);

    int nNodes = mb.nNodesInChunk(block, useCofacetCells, chunk);
    SUNDANCE_MSG2(verb(), tab1 << "num nodes in chunk = " << nNodes);

    const Array<int>& isBCIndex = *(mb.isBCIndex(block));

    /* At this point, we can start to load the elements */
    int r=0;
    RCP<Playa::LoadableVector<double> > vecBlock 
      = vec_[mvIndices[i]][block];

    FancyOStream& os = Out::os();

    for (int c=0; c<nCells; c++)
    {
      Tabs tab2;
      SUNDANCE_MSG2(verb(), tab2 << "cell = " << c << " of " << nCells);
      for (int n=0; n<nNodes; n++, r++)
      {
        Tabs tab3;
        int rowIndex = dofs[(c*nFuncs + funcIndex)*nNodes + n];
        int localRowIndex = rowIndex - lowestLocalRow;
        if (verb() >= 2) os << tab3 << "n=" << setw(4) << n 
                            << " G=" << setw(8) << rowIndex 
                            << " L=" << setw(8) << localRowIndex;
        if (!(dofMap->isLocalDOF(rowIndex))
          || isBCRqc!=isBCIndex[localRowIndex]) 
        {
          if (verb() >= 2)
          {
            if (!dofMap->isLocalDOF(rowIndex)) 
            {
              os << " --- skipping (is non-local) ---" << std::endl;
            }
            else if (!isBCRqc && isBCIndex[localRowIndex])
            {
              os << " --- skipping (is BC row) ---" << std::endl;
            }
            else
            {
              os << " --- skipping (is non-BC row) ---" << std::endl;
            }
          }
        }
        else
        {
          if (verb() >= 2) os << setw(15) << localValues[r] << std::endl;
          vecBlock->addToElement(rowIndex, localValues[r]);
        }
      }
    }
  }
  SUNDANCE_MSG1(verb(), tab0 << "...done vector insertion");
}


