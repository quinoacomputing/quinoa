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

#include "SundanceSubsetCellFilter.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

SubsetCellFilter::SubsetCellFilter(const CellFilter& superset,
                                   const CellPredicate& predicate)
  : CellFilterBase(), superset_(superset), predicate_(predicate)
{
  int verb=0;
  SUNDANCE_MSG3(verb, "creating subset cell filter: [" 
    << predicate.description()
    << "]");
  setName("Subset[sup=" + superset.description() + ", pred="
    + predicate.description()+"]");
}

SubsetCellFilter::~SubsetCellFilter()
{
//  Out::os() << "~SubsetCellFilter()" << std::endl;
}

XMLObject SubsetCellFilter::toXML() const 
{
  XMLObject rtn("SubsetCellFilter");
  rtn.addAttribute("id", Teuchos::toString(id()));
  rtn.addChild(predicate_.toXML());
  return rtn;
}


#ifdef OLD_CELL_FILTER
bool SubsetCellFilter::lessThan(const CellFilterStub* other) const
{
  const SubsetCellFilter* S 
    = dynamic_cast<const SubsetCellFilter*>(other);

  TEUCHOS_TEST_FOR_EXCEPTION(S==0,
                     std::logic_error,
                     "argument " << other->toXML() 
                     << " to SubsetCellFilter::lessThan() should be "
                     "a SubsetCellFilter pointer.");

  return OrderedPair<CellFilter, CellPredicate>(superset_, predicate_)
    < OrderedPair<CellFilter, CellPredicate>(S->superset_, S->predicate_);
}
#endif

CellSet SubsetCellFilter::internalGetCells(const Mesh& mesh) const
{
  SUNDANCE_OUT(this->verb() > 1,
                   "SubsetCellFilter::internalGetCells()");
  CellSet super = superset_.getCells(mesh);

  int dim = superset_.dimension(mesh);

  CellType cellType = mesh.cellType(dim);

  predicate_.setMesh(mesh, dim);

  ExplicitCellSet* rtn = new ExplicitCellSet(mesh, dim, cellType);

  Set<int>& cells = rtn->cells();

  const CellPredicateBase* pred = predicate_.ptr().get();


  Array<int> cellLID;

  cellLID.reserve(mesh.numCells(dim));

  for (CellIterator i=super.begin(); i != super.end(); i++)
    {
      cellLID.append(*i);
    }

  Array<int> testResults(cellLID.size());
  pred->testBatch(cellLID, testResults);

  for (int i=0; i<cellLID.size(); i++)
    {
      SUNDANCE_OUT(this->verb() > 2,
                   "SubsetCellFilter is testing " << cellLID[i]);
      if (testResults[i]) 
        {
          SUNDANCE_OUT(this->verb() > 2,
                       "accepted " << cellLID[i]);
          cells.insert(cellLID[i]);
        }
      else
        {
          SUNDANCE_OUT(this->verb() > 2,
                       "rejected " << cellLID[i]);
        }
    }

  return rtn;
}
