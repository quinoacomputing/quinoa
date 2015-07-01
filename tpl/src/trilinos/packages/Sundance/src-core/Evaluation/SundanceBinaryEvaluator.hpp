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

#ifndef SUNDANCE_BINARYEVALUATOR_H
#define SUNDANCE_BINARYEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceSubtypeEvaluator.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceEvalManager.hpp"

namespace Sundance 
{
class EvalContext;
  


/**
 * 
 */
template <class ExprType> class BinaryEvaluator 
  : public SubtypeEvaluator<ExprType>
{
public:
  /** */
  BinaryEvaluator(const ExprType* expr,
    const EvalContext& context)
    : SubtypeEvaluator<ExprType>(expr, context),
      leftExpr_(expr->leftEvaluatable()),
      rightExpr_(expr->rightEvaluatable()),
      leftSparsity_(leftExpr_->sparsitySuperset(context)),
      rightSparsity_(rightExpr_->sparsitySuperset(context)),
      leftEval_(leftExpr_->evaluator(context)),
      rightEval_(rightExpr_->evaluator(context))
    {
      Tabs tab;
      leftEval_->addClient();
      rightEval_->addClient();
    }

  /** */
  virtual ~BinaryEvaluator(){;}

  /** */
  virtual void resetNumCalls() const 
    {
      leftEval_->resetNumCalls();
      rightEval_->resetNumCalls();
      Evaluator::resetNumCalls();
    }

protected:
      
  /** */
  const RCP<SparsitySuperset>& leftSparsity() const 
    {return leftSparsity_;}
      
  /** */
  const RCP<SparsitySuperset>& rightSparsity() const 
    {return rightSparsity_;}

  /** */
  const EvaluatableExpr* leftExpr() const {return leftExpr_;}

  /** */
  const EvaluatableExpr* rightExpr() const {return rightExpr_;}

  /** */
  const RCP<Evaluator>& leftEval() const {return leftEval_;}

  /** */
  const RCP<Evaluator>& rightEval() const {return rightEval_;}

  /** */
  void evalChildren(const EvalManager& mgr,
    Array<double>& leftConstResults,
    Array<RCP<EvalVector> >& leftVecResults,
    Array<double>& rightConstResults,
    Array<RCP<EvalVector> >& rightVecResults) const 
    {
      Tabs tabs;
      SUNDANCE_MSG2(mgr.verb(), 
        tabs << "Evaluating left and right children: "
        << std::endl << tabs << "left=" << leftExpr_->toString()
        << std::endl << tabs << "right=" << rightExpr_->toString());
      SUNDANCE_MSG2(mgr.verb(),
        tabs << "Evaluating left=" << leftExpr_->toString());
      leftEval()->eval(mgr, leftConstResults, leftVecResults);
        
      SUNDANCE_MSG2(mgr.verb(),
        tabs << "Evaluating right=" << rightExpr_->toString());
      rightEval()->eval(mgr, rightConstResults, rightVecResults);
    }

private:
  const EvaluatableExpr* leftExpr_;

  const EvaluatableExpr* rightExpr_;

  RCP<SparsitySuperset> leftSparsity_;

  RCP<SparsitySuperset> rightSparsity_;

  RCP<Evaluator> leftEval_;

  RCP<Evaluator> rightEval_;
};
}


#endif
