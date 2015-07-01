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

#ifndef SUNDANCE_UNARYEVALUATOR_H
#define SUNDANCE_UNARYEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceSubtypeEvaluator.hpp"
#include "SundanceEvaluatableExpr.hpp"


namespace Sundance 
{
class EvalContext;


/**
 * 
 */
template <class ExprType> class UnaryEvaluator 
  : public SubtypeEvaluator<ExprType>
{
public:
  /** */
  UnaryEvaluator(const ExprType* expr,
    const EvalContext& context)
    : SubtypeEvaluator<ExprType>(expr, context),
      argExpr_(expr->evaluatableArg()),
      argSparsitySuperset_(argExpr_->sparsitySuperset(context)),
      argEval_(argExpr_->evaluator(context))
    {
      try
      {
        Tabs tab;
        int verb = context.evalSetupVerbosity();

        SUNDANCE_MSG2(verb, tab << "UnaryEvaluator ctor: expr = " << expr->toString());

        SUNDANCE_MSG2(verb, tab << "arg sparsity superset maxOrder: " 
          << argSparsitySuperset_->maxOrder());
            
        argEval_->addClient();
            
        SUNDANCE_MSG2(verb, tab << "done unary evalulator ctor");
      }
      catch(std::exception& e)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, 
          "exception detected in UnaryEvaluator: expr="
          << expr->toString() << std::endl
          << "arg=" << expr->evaluatableArg()->toString() << std::endl
          << "exception=" << e.what());
      }
    }

  /** */
  virtual ~UnaryEvaluator(){;}

  /** */
  virtual void resetNumCalls() const 
    {
      argEval_->resetNumCalls();
      Evaluator::resetNumCalls();
    }

protected:

  /** */
  const RCP<SparsitySuperset>& argSparsitySuperset() const 
    {return argSparsitySuperset_;}
      
  /** */
  const EvaluatableExpr* argExpr() const {return argExpr_;}

  /** */
  const RCP<Evaluator>& argEval() const {return argEval_;}
      

  /** */
  void evalOperand(const EvalManager& mgr,
    Array<double>& argConstantResults,
    Array<RCP<EvalVector> >& argVectorResults) const 
    {
      Tabs tabs;
      SUNDANCE_MSG1(this->verb(),  tabs << "UnaryEvaluator: evaluating operand: ");
      argEval()->eval(mgr, argConstantResults, argVectorResults);
      SUNDANCE_MSG1(this->verb(),  tabs << "UnaryEvaluator: done eval operand ");
    }
private:
  const EvaluatableExpr* argExpr_;

  RCP<SparsitySuperset> argSparsitySuperset_;

  RCP<Evaluator> argEval_;
};
}

#endif
