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

#include "SundanceDiscreteFuncEvaluator.hpp"
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



DiscreteFuncElementEvaluator
::DiscreteFuncElementEvaluator(const DiscreteFuncElement* expr, 
                               const EvalContext& context)
  : SubtypeEvaluator<DiscreteFuncElement>(expr, context), 
    mi_(this->sparsity()->numDerivs()),
    miToIndexMap_(),
    stringReps_()
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "initializing discrete func evaluator for " 
                    << expr->toString());

  SUNDANCE_VERB_MEDIUM(tabs << "return sparsity " << std::endl << *(this->sparsity()));

  static Array<string> coordNames;
  if (coordNames.size() != 3)
    {
      coordNames.resize(3);
      coordNames[0] = "x";
      coordNames[1] = "y";
      coordNames[2] = "z";
    }
  std::string funcName = expr->name();

  for (int i=0; i<this->sparsity()->numDerivs(); i++)
    {
      /* Make sure that every derivative we're to evaluate is either
      * zero-order or a spatial derivative */
      if (this->sparsity()->deriv(i).order()==0) 
        {
          mi_[i] = MultiIndex();
        }
      else 
        {
      
          TEUCHOS_TEST_FOR_EXCEPTION(!this->sparsity()->isSpatialDeriv(i), std::logic_error,
                             "DiscreteFuncElementEvaluator ctor found "
                             "an entry in the sparsity superset that is not "
                             "a spatial derivative. "
                             "The bad entry is " << this->sparsity()->deriv(i) 
                             << ". The superset is " 
                             << *(this->sparsity)());

          mi_[i] = this->sparsity()->multiIndex(i);
        }
      addVectorIndex(i,i);
      TEUCHOS_TEST_FOR_EXCEPTION(miToIndexMap_.containsKey(mi_[i]), std::logic_error,
                         "DiscreteFuncElementEvaluator ctor detected a "
                         "duplicate multiindex");

      miToIndexMap_.put(mi_[i], i);

      if (mi_[i].order()==0)
        {
          stringReps_.append(funcName);
        }
      else
        {
          int dir = mi_[i].firstOrderDirection();
          std::string deriv = "D[" + funcName + ", " + coordNames[dir] + "]";
          stringReps_.append(deriv);
        }
    }

  
}

bool DiscreteFuncElementEvaluator::hasMultiIndex(const MultiIndex& mi) const
{
  Tabs tabs;
  bool rtn = miToIndexMap_.containsKey(mi);
  SUNDANCE_VERB_MEDIUM(tabs << "checking for mi=" << mi << " for " 
                       << expr()->toString()
                       << std::endl << tabs 
                       << " sparsity " << std::endl << *(this->sparsity()));
  
  return rtn;
}

int DiscreteFuncElementEvaluator::miIndex(const MultiIndex& mi) const
{
  return miToIndexMap_.get(mi);
}


void DiscreteFuncElementEvaluator
::internalEval(const EvalManager& mgr,
               Array<double>& constantResults,
               Array<RCP<EvalVector> >& vectorResults) const 
{
  Tabs tabs(0);
  SUNDANCE_MSG1(mgr.verb(),
    tabs << "DiscreteFuncElementEvaluator::eval: expr=" 
    << expr()->toString());

  vectorResults.resize(mi_.size());
  for (int i=0; i<mi_.size(); i++)
    {
      Tabs tab2;
      vectorResults[i] = mgr.popVector();
      TEUCHOS_TEST_FOR_EXCEPTION(!vectorResults[i]->isValid(), 
                         std::logic_error,
                         "invalid evaluation vector allocated in "
                         "DiscreteFuncElementEvaluator::internalEval()");
      SUNDANCE_MSG2(mgr.verb(),tab2<< "setting string rep " << stringReps_[i]);
      vectorResults[i]->setString(stringReps_[i]);
    }
  mgr.evalDiscreteFuncElement(expr(), mi_, vectorResults);
  mgr.stack().setVecSize(vectorResults[0]->length());
  
  if (mgr.verb() > 1)
    {
      Out::os() << tabs << "results " << std::endl;
      mgr.showResults(Out::os(), sparsity(), vectorResults,
                            constantResults);
    }
  SUNDANCE_MSG1(mgr.verb(), tabs << "DiscreteFuncEvaluator::eval() done"); 
}

