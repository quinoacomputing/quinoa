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
#include "SundanceCellDiameterExpr.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceSet.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


CellDiameterExprEvaluator::CellDiameterExprEvaluator(
  const CellDiameterExpr* expr, 
  const EvalContext& context)
  : SubtypeEvaluator<CellDiameterExpr>(expr, context), 
    stringRep_(expr->toString())
{

  Tabs tabs;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG1(verb, tabs << "initializing cell diameter expr evaluator for " 
                    << expr->toString());
  SUNDANCE_MSG2(verb, tabs << "return sparsity " << std::endl << *(this->sparsity)());

  TEUCHOS_TEST_FOR_EXCEPTION(this->sparsity()->numDerivs() > 1, std::logic_error,
                     "CellDiameterExprEvaluator ctor found a sparsity table "
                     "with more than one entry. The bad sparsity table is "
                     << *(this->sparsity)());

  /* 
   * There is only one possible entry in the nozeros table for a
   * cell diameter expression: a zeroth derivative.
   */
  
  for (int i=0; i<this->sparsity()->numDerivs(); i++)
    {
      const MultipleDeriv& d = this->sparsity()->deriv(i);

      TEUCHOS_TEST_FOR_EXCEPTION(d.order()!=0, std::logic_error,
                         "CellDiameterExprEvaluator ctor found an entry in the "
                         "sparsity superset that is not a zeroth-order derivative. "
                         "The bad entry is " << this->sparsity()->deriv(i) 
                         << ". The superset is " 
                         << *(this->sparsity)());
      addVectorIndex(i, 0);
    }
}



void CellDiameterExprEvaluator::internalEval(const EvalManager& mgr,
  Array<double>& constantResults,
  Array<RCP<EvalVector> >& vectorResults) const 
{
  Tabs tabs(0);

  SUNDANCE_MSG2(mgr.verb(), tabs << "CellDiameterExprEvaluator::eval() expr=" 
    << expr()->toString());

  if (mgr.verb() > 2)
    {
      Out::os() << tabs << "sparsity = " 
                << std::endl << tabs << *(this->sparsity)() << std::endl;
    }

  if (this->sparsity()->numDerivs() > 0)
    {
      vectorResults.resize(1);
      vectorResults[0] = mgr.popVector();
      SUNDANCE_MSG3(mgr.verb(), tabs << "forwarding to evaluation manager");
      mgr.evalCellDiameterExpr(expr(), vectorResults[0]);
      mgr.stack().setVecSize(vectorResults[0]->length());
      if (EvalVector::shadowOps()) vectorResults[0]->setString(stringRep_);
    }
  else
  {
    SUNDANCE_MSG4(mgr.verb(), tabs << "no results requested");
  }

  if (mgr.verb() > 1)
    {
      Out::os() << tabs << "results " << std::endl;
      mgr.showResults(Out::os(), this->sparsity(), vectorResults,
		      constantResults);
    }
}

