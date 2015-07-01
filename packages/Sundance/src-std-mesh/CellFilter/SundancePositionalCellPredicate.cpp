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

#include "SundancePositionalCellPredicate.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

bool PositionalCellPredicate::lessThan(const CellPredicateBase* other) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(dynamic_cast<const PositionalCellPredicate*>(other) == 0,
                     std::logic_error,
                     "argument " << other->toXML() 
                     << " to PositionalCellPredicate::lessThan() should be "
                     "a PositionalCellPredicate pointer.");

  return func_.get() < dynamic_cast<const PositionalCellPredicate*>(other)->func_.get();
}

void PositionalCellPredicate::testBatch(const Array<int>& cellLID,
                                        Array<int>& results) const
{
  results.resize(cellLID.size());

  if (cellDim()==0)
    {
      for (int i=0; i<cellLID.size(); i++)
        {
          results[i] = (*func_)(mesh().nodePosition(cellLID[i]));
        }
    }
  else
    {
      Array<int> facetLIDs;
      Array<int> facetSigns;
      int nf = mesh().numFacets(cellDim(), cellLID[0], 0);
      mesh().getFacetLIDs(cellDim(), cellLID, 0, facetLIDs, facetSigns);
      for (int c=0; c<cellLID.size(); c++)
        {
          results[c] = true;
          for (int f=0; f<nf; f++)
            {
              int fLID = facetLIDs[c*nf + f];
              if ((*func_)(mesh().nodePosition(fLID)) == false)
                {
                  results[c] = false;
                  continue;
                }
            }
        }
    }
}

XMLObject PositionalCellPredicate::toXML() const 
{
  XMLObject rtn("PositionalCellPredicate");
  return rtn;
}

bool PointCellPredicateFunctor::operator()(const Point& x) const
{
  for (int i=0; i<x_.dim(); i++)
  {
    if (fabs(x[i]-x_[i]) > tol_) return false;
  }
  return true;
}


bool CoordinateValueCellPredicateFunctor::operator()(const Point& x) const
{
  return (fabs(x[direction_]-value_) < tol_);
}
