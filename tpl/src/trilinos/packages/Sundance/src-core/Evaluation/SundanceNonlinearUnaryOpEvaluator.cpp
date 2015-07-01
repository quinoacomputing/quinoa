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

#include "SundanceNonlinearUnaryOpEvaluator.hpp"
#include "SundanceEvalManager.hpp"

#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceNonlinearUnaryOp.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;





NonlinearUnaryOpEvaluator
::NonlinearUnaryOpEvaluator(const NonlinearUnaryOp* expr,
                            const EvalContext& context)
  : ChainRuleEvaluator(expr, context),
    op_(expr->op()),
    maxOrder_(-1),
    argIsConstant_(true),
    argValueIndex_(-1)
{
  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "NonlinearUnaryOp evaluator ctor for " 
                    << expr->toString());

  
  const SparsitySuperset* sArg = childSparsity(0);
  SUNDANCE_VERB_HIGH(tabs << "arg sparsity " << *sArg);
  maxOrder_ = expr->sparsitySuperset(context)->maxOrder();

  /* See if the argument is non-constant */
  for (int i=0; i<sArg->numDerivs(); i++)
    {
      if (sArg->state(i) == VectorDeriv) 
        {
          argIsConstant_ = false;
          break;
        }
    }

  /* We will need derivatives of the operator wrt the argument value for
   * orders up to maxOrder. Define their locations in the argDeriv array. */
  for (int order=0; order<=maxOrder_; order++)
    {
      MultiSet<int> df;
      /* The index set for the N-th order arg deriv is the arg index (0)
       * replicated N times. */
      for (int i=0; i<order; i++) df.put(0);

      if (argIsConstant_) addConstArgDeriv(df, order);
      else addVarArgDeriv(df, order);
    }

  
  /* Find the index of the argument value (zeroth-order deriv) in the 
   * vector of derivatives of the argument */
  MultipleDeriv d0;
  TEUCHOS_TEST_FOR_EXCEPTION(!sArg->containsDeriv(d0), std::logic_error,
                     "NonlinearUnaryOpEvaluator::ctor did not find zeroth-order "
                     "derivative of argument");
  int d0Index = sArg->getIndex(d0);
  const Evaluator* argEv = childEvaluator(0);
  if (argIsConstant_) argValueIndex_ = argEv->constantIndexMap().get(d0Index);
  else argValueIndex_ = argEv->vectorIndexMap().get(d0Index);


  SUNDANCE_VERB_HIGH(tabs << "arg is constant: " << argIsConstant_);
  SUNDANCE_VERB_HIGH(tabs << "arg value index: " << argValueIndex_);
  /* Call init() at the base class to set up chain rule evaluation */
  init(expr, context);
}


void NonlinearUnaryOpEvaluator
::evalArgDerivs(const EvalManager& mgr,
                const Array<RCP<Array<double> > >& constArgRes,
                const Array<RCP<Array<RCP<EvalVector> > > >& vArgResults,
                Array<double>& constArgDerivs,
                Array<RCP<EvalVector> >& varArgDerivs) const
{
  Tabs tabs;
  SUNDANCE_MSG1(mgr.verb(), tabs 
    << "NonlinearUnaryOpEvaluator::evalArgDerivs() for " 
    << expr()->toString());
  
  if (argIsConstant_)
    {
      double argValue = (*(constArgRes[0]))[argValueIndex_];
      constArgDerivs.resize(maxOrder_+1);
      switch(maxOrder_)
        {
        case 0:
          op_->eval0(&argValue, 1, &(constArgDerivs[0]));
          break;
        case 1:
          op_->eval1(&argValue, 1, &(constArgDerivs[0]), &(constArgDerivs[1]));
          break;
        case 2:
          op_->eval2(&argValue, 1, &(constArgDerivs[0]), &(constArgDerivs[1]),
                     &(constArgDerivs[2]));
          break;
        case 3:
          op_->eval3(&argValue, 1, &(constArgDerivs[0]), &(constArgDerivs[1]),
                     &(constArgDerivs[2]), &(constArgDerivs[3]));
          break;
        default:
          TEUCHOS_TEST_FOR_EXCEPT(true);
        }
    }
  else
    {
      Tabs tab1;
      SUNDANCE_MSG2(mgr.verb(), 
        tab1 << "argument is a vector: argValIndex=" << argValueIndex_);
      SUNDANCE_MSG2(mgr.verb(), 
        tab1 << "num vector results =" << vArgResults.size());

      const Array<RCP<EvalVector> >& av = *(vArgResults[0]);

      SUNDANCE_MSG2(mgr.verb(), 
        tab1 << "av size=" << av.size());


      const double* argValue = av[argValueIndex_]->start();
      varArgDerivs.resize(maxOrder_+1);
      varArgDerivs[0] = mgr.stack().popVector();
      int nx = varArgDerivs[0]->length();

      const std::string& argStr = (*(vArgResults[0]))[argValueIndex_]->str();
      if (EvalVector::shadowOps())
        {
          varArgDerivs[0]->setString(op_->name() + "(" + argStr + ")");
        }

      double* f = varArgDerivs[0]->start();
      double* df_dArg = 0 ;
      if (maxOrder_ >= 1) 
        {
          varArgDerivs[1] = mgr.stack().popVector();
          if (EvalVector::shadowOps())
            {
              varArgDerivs[1]->setString(op_->name() + "'(" + argStr + ")");
            }
          df_dArg = varArgDerivs[1]->start();
        }
      double* d2f_dArg2 = 0 ;
      if (maxOrder_ >= 2) 
        {
          varArgDerivs[2] =mgr.stack().popVector();
          if (EvalVector::shadowOps())
            {
              varArgDerivs[2]->setString(op_->name() + "''(" + argStr + ")");
            }
          d2f_dArg2 = varArgDerivs[2]->start();
        }
      double* d3f_dArg3 = 0 ;
      if (maxOrder_ >= 3) 
        {
          UnaryFunctor::fdStep()=1.0e-3;
          varArgDerivs[3] = mgr.stack().popVector();
          if (EvalVector::shadowOps())
            {
              varArgDerivs[3]->setString(op_->name() + "'''(" + argStr + ")");
            }
          d3f_dArg3 = varArgDerivs[3]->start();
        }

      switch(maxOrder_)
        {
        case 0:
          op_->eval0(argValue, nx, f);
          break;
        case 1:
          TEUCHOS_TEST_FOR_EXCEPT(df_dArg==0);
          op_->eval1(argValue, nx, f, df_dArg);
          break;
        case 2:
          TEUCHOS_TEST_FOR_EXCEPT(df_dArg==0);
          TEUCHOS_TEST_FOR_EXCEPT(d2f_dArg2==0);
          op_->eval2(argValue, nx, f, df_dArg, d2f_dArg2);
          break;
        case 3:
          TEUCHOS_TEST_FOR_EXCEPT(df_dArg==0);
          TEUCHOS_TEST_FOR_EXCEPT(d2f_dArg2==0);
          TEUCHOS_TEST_FOR_EXCEPT(d3f_dArg3==0);
          op_->eval3(argValue, nx, f, df_dArg, d2f_dArg2, d3f_dArg3);
          break;
        default:
          TEUCHOS_TEST_FOR_EXCEPT(true);
        }
    }

  SUNDANCE_MSG1(mgr.verb(),
    tabs << "done NonlinearUnaryOpEvaluator::evalArgDerivs() for " 
    << expr()->toString());
}



