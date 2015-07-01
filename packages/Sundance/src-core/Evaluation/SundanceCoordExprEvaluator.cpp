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

#include "SundanceSubtypeEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceCoordExpr.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceSet.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


CoordExprEvaluator::CoordExprEvaluator(const CoordExpr* expr, 
                                       const EvalContext& context)
  : SubtypeEvaluator<CoordExpr>(expr, context), 
    doValue_(false),
    doDeriv_(false),
    stringRep_(expr->toString())
{
  int verb = context.setupVerbosity();
  Tabs tabs;
  SUNDANCE_MSG1(verb, tabs << "initializing coord expr evaluator for " 
                    << expr->toString());
  SUNDANCE_MSG2(verb, tabs << "return sparsity " << std::endl << tabs << *(this->sparsity)());

  TEUCHOS_TEST_FOR_EXCEPTION(this->sparsity()->numDerivs() > 2, std::logic_error,
                     "CoordExprEvaluator ctor found a sparsity table "
                     "with more than two entries. The bad sparsity table is "
                     << *(this->sparsity)());

  /* 
   * There are only two possible entries in the nozeros table for a
   * coordinate expression: a zeroth derivative, and a first-order
   * spatial derivative in the same direction as the expr's coordinate. 
   */
  
  for (int i=0; i<this->sparsity()->numDerivs(); i++)
    {
      const MultipleDeriv& d = this->sparsity()->deriv(i);

      /* for a zeroth-order derivative, evaluate the coord expr */
      if (d.order()==0) 
        {
          doValue_ = true;
          addVectorIndex(i, 0);
        }
      else /* for a first-order deriv, make sure it's in the proper direction,
            * then evaluate the spatial derivative. */
        {
          TEUCHOS_TEST_FOR_EXCEPTION(!this->sparsity()->isSpatialDeriv(i), std::logic_error,
                             "CoordExprEvaluator ctor found an entry in the "
                             "sparsity superset that is not a spatial derivative. "
                             "The bad entry is " << this->sparsity()->deriv(i) 
                             << ". The superset is " 
                             << *(this->sparsity)());

          const MultiIndex& mi = this->sparsity()->multiIndex(i);
          
          TEUCHOS_TEST_FOR_EXCEPTION(mi.order() != 1, std::logic_error,
                             "CoordExprEvaluator ctor found a multiindex of "
                             "order != 1. Bad multiindex is " << mi.toString());
          
          TEUCHOS_TEST_FOR_EXCEPTION(mi[expr->dir()]!=1, std::logic_error,
                             "CoordExprEvaluator sparsity pattern has an "
                             "element corresponding to differentiation wrt "
                             "a coordinate direction other than that of the "
                             "coord expr's direction");
          doDeriv_ = true;
          addConstantIndex(i, 0);
        }
    }
  
}



void CoordExprEvaluator::internalEval(const EvalManager& mgr,
  Array<double>& constantResults,
  Array<RCP<EvalVector> >& vectorResults) const 
{
  Tabs tabs(0);

  SUNDANCE_MSG1(mgr.verb(), tabs << "CoordExprEvaluator::eval() expr=" << expr()->toString());

  SUNDANCE_MSG2(mgr.verb(), tabs << "sparsity = " << std::endl 
    << *(this->sparsity)())

  if (doValue_)
    {
      Tabs tab2;
      SUNDANCE_MSG3(mgr.verb(), tab2 << "computing value");
      vectorResults.resize(1);
      vectorResults[0] = mgr.popVector();
      mgr.evalCoordExpr(expr(), vectorResults[0]);
      mgr.stack().setVecSize(vectorResults[0]->length());
      vectorResults[0]->setString(stringRep_);
    }
  
  if (doDeriv_)
    {
      Tabs tab2;
      SUNDANCE_MSG3(mgr.verb(), tab2 << "computing derivative");
      constantResults.resize(1);
      constantResults[0] = 1.0;
    }

  if (mgr.verb() > 1)
    {
      Tabs tab1;
      Out::os() << tab1 << "results " << std::endl;
      mgr.showResults(Out::os(), this->sparsity(), vectorResults,
                            constantResults);
    }

}

