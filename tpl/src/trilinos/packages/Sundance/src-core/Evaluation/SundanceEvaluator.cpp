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

#include "SundanceEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceSet.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;




Evaluator::Evaluator()
  : numClients_(0),
    numCalls_(0),
    vectorResultCache_(),
    constantResultCache_(),
    constantIndexMap_(),
    vectorIndexMap_(),
    vectorIndices_(),
    constantIndices_()
{}



void Evaluator::eval(const EvalManager& mgr,
                     Array<double>& constantResults,
                     Array<RCP<EvalVector> >& vectorResults) const
{
  Tabs tabs;

  if (numCalls_ == 0)
    {
      internalEval(mgr, constantResultCache_, vectorResultCache_);
    }
  
  numCalls_++;

  /* Go ahead and copy the constant results every time, 
   * since this is cheap */
  if (constantResultCache_.size() > 0) constantResults = constantResultCache_;

  /* If all clients have called, we can return the original data
   * which can then be changed by the client. */
  if (numCalls_ == numClients_)
    {
      SUNDANCE_VERB_MEDIUM(tabs << "surrendering cached results");
      vectorResults = vectorResultCache_;
    }
  else /* Otherwise, make a copy of the data */
    {
      SUNDANCE_VERB_MEDIUM(tabs << "cloning cached results");
      vectorResults.resize(vectorResultCache_.size());
      for (int i=0; i < vectorResults.size(); i++)
        {
          vectorResults[i] = vectorResultCache_[i]->clone();
        }
    }
}


void Evaluator::addConstantIndex(int index, int constantIndex)
{
  TEUCHOS_TEST_FOR_EXCEPTION(constantIndexMap_.containsKey(index), std::logic_error,
                     "duplicate index " << index 
                     << " found in Evaluator::addConstantIndex");
  constantIndexMap_.put(index, constantIndex);
  constantIndices_.append(index);
}

void Evaluator::addVectorIndex(int index, int vectorIndex)
{
  TEUCHOS_TEST_FOR_EXCEPTION(vectorIndexMap_.containsKey(index), std::logic_error,
                     "duplicate index " << index 
                     << " found in Evaluator::addVectorIndex");
  vectorIndexMap_.put(index, vectorIndex);
  vectorIndices_.append(index);
}

