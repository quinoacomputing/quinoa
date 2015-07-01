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

#include "SundanceUserDefOpEvaluator.hpp"
#include "SundanceUserDefOpCommonEvaluator.hpp"
#include "SundanceUserDefOpElement.hpp"
#include "SundanceEvalManager.hpp"

#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceUserDefOp.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;





UserDefOpEvaluator
::UserDefOpEvaluator(const UserDefOpElement* expr,
                     const RCP<const UserDefOpCommonEvaluator>& commonEval,
                     const EvalContext& context)
  : ChainRuleEvaluator(expr, context),
    argValueIndex_(expr->numChildren()),
    argValueIsConstant_(expr->numChildren()),
    functor_(expr->functorElement()),
    commonEval_(commonEval),
    maxOrder_(0),
    numVarArgDerivs_(0),
    numConstArgDerivs_(0),
    allArgsAreConstant_(true)
{
  Tabs tab1;
  SUNDANCE_VERB_LOW(tab1 << "initializing user defined op evaluator for " 
                    << expr->toString());
  Array<int> orders = findRequiredOrders(expr, context);

  for (int i=0; i<orders.size(); i++) 
    {
      if (orders[i] > maxOrder_) maxOrder_ = orders[i];
    }
  commonEval->updateMaxOrder(maxOrder_);

  SUNDANCE_VERB_HIGH(tab1 << "setting arg deriv indices");
  
  
  /* Find the mapping from argument derivatives to indices in the 
   * functor's vector of return values */
  Map<MultiSet<int>, int> varArgDerivs;
  Map<MultiSet<int>, int> constArgDerivs;
  expr->getArgDerivIndices(orders, varArgDerivs, constArgDerivs);
  numVarArgDerivs_ = varArgDerivs.size();
  numConstArgDerivs_ = constArgDerivs.size();
  typedef Map<MultiSet<int>, int>::const_iterator iter;
  for (iter i=varArgDerivs.begin(); i!=varArgDerivs.end(); i++)
    {
      Tabs tab2;
      SUNDANCE_VERB_EXTREME(tab2 << "variable arg deriv " << i->first 
                            << " will be at index "
                            << i->second);
      addVarArgDeriv(i->first, i->second);
    }
  
  for (iter i=constArgDerivs.begin(); i!=constArgDerivs.end(); i++)
    {
      Tabs tab2;
      SUNDANCE_VERB_EXTREME(tab2 << "constant arg deriv " << i->first 
                            << " will be at index "
                            << i->second);
      addConstArgDeriv(i->first, i->second);
    }

  /* Find the indices to the zeroth derivative of each argument */
  
  for (int i=0; i<expr->numChildren(); i++)
    {
      const SparsitySuperset* sArg = childSparsity(i);
      int numConst=0;
      int numVec=0;
      for (int j=0; j<sArg->numDerivs(); j++)
        {
          if (sArg->deriv(j).order() == 0) 
            {
              if (sArg->state(j)==VectorDeriv)
                {
                  argValueIndex_[i] = numVec;              
                  allArgsAreConstant_ = false;
                }
              else
                {
                  argValueIndex_[i] = numConst;              
                }
              break;
            }
          if (sArg->state(j) == VectorDeriv) 
            {
              numVec++;
            }
          else
            {
              numConst++;
            }
        }
    }
  
  /* Call init() at the base class to set up chain rule evaluation */
  init(expr, context);
}




void UserDefOpEvaluator::resetNumCalls() const
{
  commonEval()->markCacheAsInvalid();
  ChainRuleEvaluator::resetNumCalls();
}




Array<int> UserDefOpEvaluator::findRequiredOrders(const ExprWithChildren* expr, 
                                                  const EvalContext& context)
{
  Tabs tab0;
  SUNDANCE_VERB_HIGH(tab0 << "finding required arg deriv orders");

  Set<int> orders;
  
  const Set<MultipleDeriv>& R = expr->findR(context);
  typedef Set<MultipleDeriv>::const_iterator iter;

  for (iter md=R.begin(); md!=R.end(); md++)
    {
      Tabs tab1;
      
      int N = md->order();
      if (N > maxOrder_) maxOrder_ = N;
      if (N==0) orders.put(N);
      for (int n=1; n<=N; n++)
        {
          const Set<MultiSet<int> >& QW = expr->findQ_W(n, context);
          for (Set<MultiSet<int> >::const_iterator q=QW.begin(); q!=QW.end(); q++)
            {
              orders.put(q->size());
            }
        }
    }
  SUNDANCE_VERB_HIGH(tab0 << "arg deriv orders=" << orders);
  return orders.elements();
}




void UserDefOpEvaluator
::evalArgDerivs(const EvalManager& mgr,
                const Array<RCP<Array<double> > >& constArgVals,
                const Array<RCP<Array<RCP<EvalVector> > > >& varArgVals,
                Array<double>& constArgDerivs,
                Array<RCP<EvalVector> >& varArgDerivs) const
{
  if (!commonEval()->cacheIsValid())
    {
      commonEval()->evalAllComponents(mgr, constArgVals, varArgVals);
    }
  if (allArgsAreConstant_)
    {
      constArgDerivs = commonEval()->constArgDerivCache(myIndex());
    }
  else
    {
      varArgDerivs = commonEval()->varArgDerivCache(myIndex());
    }
}

