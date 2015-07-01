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
#include "SundanceMapBundle.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceStdFwkEvalMediator.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

      

MapBundle::MapBundle(
  const Array<RCP<DOFMapBase> >& dofMap,
  const Array<RCP<Array<int> > >& isBCIndex,
  const Array<int>& lowestLocalIndex,
  bool partitionBCs,
  int verb
  )
  : verb_(verb),
    dofMap_(dofMap),
    isBCIndex_(isBCIndex),
    lowestLocalIndex_(lowestLocalIndex),
    localDOFMap_(rcp(new LocalDOFMap(dofMap_.size(), verb))),
    cofacetLocalDOFMap_(rcp(new LocalDOFMap(dofMap_.size(), verb)))
{}

int MapBundle::nCells() const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(localDOFMap_->isUnused() 
    && cofacetLocalDOFMap_->isUnused(), std::runtime_error,
    "no local DOF maps defined in call to MapBundle::nCells()");

  
  if (cofacetLocalDOFMap_->isUnused())
  {
    return localDOFMap_->nCells();
  }
  else if (localDOFMap_->isUnused())
  {
    return cofacetLocalDOFMap_->nCells();
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(localDOFMap_->nCells() != cofacetLocalDOFMap_->nCells(),
      std::runtime_error,
      "mismatched cell counts in MapBundle::nCells()");
    return cofacetLocalDOFMap_->nCells();
  }
}

RCP<const Array<int> > MapBundle::workSet(int block,
  bool useCofacets) const
{
  return chooseMap(block, useCofacets)->cellLIDs();
}


const RCP<LocalDOFMap>& MapBundle::chooseMap(
  int block, bool useCofacets) const
{
  if (useCofacets)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(cofacetLocalDOFMap_->isUnused(block),
      std::runtime_error,
      "request for unavailable cofacet-based local map for block = " << block);
    return cofacetLocalDOFMap_;
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(localDOFMap_->isUnused(block),
      std::runtime_error,
      "request for unavailable local map for block = " << block);
    return localDOFMap_;
  }
}




void MapBundle::buildLocalDOFMaps(
  const RCP<StdFwkEvalMediator>& mediator,
  IntegrationCellSpecifier intCellSpec,
  const Array<Set<int> >& requiredFuncs,
  int verbosity)
{
  Tabs tab;

  int numBlocks = dofMap_.size();

  localDOFMap_->markAsUnused();
  cofacetLocalDOFMap_->markAsUnused();
  localDOFMap_->setVerb(verbosity);
  cofacetLocalDOFMap_->setVerb(verbosity);

  int maxCellDim = mediator->maxCellDim();
  int cellDim = mediator->cellDim();

  SUNDANCE_MSG3(verbosity, tab << "cell dim=" << cellDim);
  SUNDANCE_MSG3(verbosity, tab << "max cell dim=" << maxCellDim);

  for (int b=0; b<numBlocks; b++)
  {   
    Tabs tab2;
    SUNDANCE_MSG3(verbosity, tab2 << "getting dofs for block " 
      << b << " of " << numBlocks);
        
    if (intCellSpec != AllTermsNeedCofacets)
    {
      Tabs tab3;
      SUNDANCE_MSG3(verbosity, tab3 << "getting ordinary dofs");

      if (!localDOFMap_->hasCells()) 
      {
        SUNDANCE_MSG3(verbosity, tab3 << "setting cells of dim " 
          << cellDim);
        localDOFMap_->setCells(cellDim, maxCellDim, mediator->cellLID());
      }     
      localDOFMap_->fillBlock(b, dofMap_[b], requiredFuncs);
      localDOFMap_->markAsUsed(b);
    }
    else
    {
      Tabs tab3;
      SUNDANCE_MSG3(verbosity, tab3 << "ordinary dofs not needed for block " << b);
    }

        
    if (intCellSpec != NoTermsNeedCofacets)
    {
      Tabs tab3;
      SUNDANCE_MSG3(verbosity, tab3 << "getting cofacet dofs");
      SUNDANCE_MSG3(verbosity, tab3 << "cofacet cells="
        << *(mediator->cofacetCellLID()));

      if (!cofacetLocalDOFMap_->hasCells()) 
        cofacetLocalDOFMap_->setCells(maxCellDim, 
          maxCellDim, mediator->cofacetCellLID());

      cofacetLocalDOFMap_->fillBlock(b, dofMap_[b], requiredFuncs);
      cofacetLocalDOFMap_->markAsUsed(b);
    }
    else
    {
      Tabs tab3;
      SUNDANCE_MSG3(verbosity, tab3 << "cofacet dofs not needed for block " << b);
    }
        
  }

  if (intCellSpec != AllTermsNeedCofacets)
  {
    SUNDANCE_MSG4(verbosity, tab << "local DOF values " << *localDOFMap_);
  }

  if (intCellSpec != NoTermsNeedCofacets)
  {
    SUNDANCE_MSG4(verbosity, tab << "local cofacet DOF values " 
      << *cofacetLocalDOFMap_);
  }
}
