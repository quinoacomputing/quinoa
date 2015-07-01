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

#include "SundanceLabelCellPredicate.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

bool LabelCellPredicate::lessThan(const CellPredicateBase* other) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(dynamic_cast<const LabelCellPredicate*>(other) == 0,
                     std::logic_error,
                     "argument " << other->toXML() 
                     << " to LabelCellPredicate::lessThan() should be "
                     "a LabelCellPredicate pointer.");

  const LabelCellPredicate* op 
    = dynamic_cast<const LabelCellPredicate*>(other);
  
  return (labelIndices_ < op->labelIndices_);
}

void LabelCellPredicate::testBatch(const Array<int>& cellLID,
                                   Array<int>& results) const
{
  mesh().getLabels(cellDim(), cellLID, results);
  for (int i=0; i<cellLID.size(); i++)
    {
      int cellLabel = results[i];
      results[i] = labelIndices_.contains(cellLabel);
    }
}

XMLObject LabelCellPredicate::toXML() const 
{
  XMLObject rtn("LabelCellPredicate");
  rtn.addAttribute("label", labelIndices_.toString());
  return rtn;
}

