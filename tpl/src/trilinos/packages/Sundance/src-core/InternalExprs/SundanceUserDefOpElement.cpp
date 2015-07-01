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

#include "SundanceUserDefOpElement.hpp"
#include "SundanceUserDefOp.hpp"
#include "SundanceUserDefOpEvaluator.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Teuchos;


UserDefOpElement::UserDefOpElement(const Array<RCP<ScalarExpr> >& args,
  const RCP<Sundance::Map<EvalContext, RCP<const UserDefOpCommonEvaluator> > >& evalMap,
  const RCP<const UserDefFunctorElement>& functorElement)
  : ExprWithChildren(args), 
    commonEvaluatorsMap_(evalMap),
    functorElement_(functorElement)
{}



void UserDefOpElement::getArgDerivIndices(const Array<int>& orders,
  Map<MultiSet<int>, int>& varArgDerivs,
  Map<MultiSet<int>, int>& constArgDerivs) const
{
  int n = functorElement()->numArgs();

  int counter = 0;

  for (int o=0; o<orders.size(); o++)
  {
    int order = orders[o];
      
    if (order==0)
    {
      varArgDerivs.put(MultiSet<int>(), counter++);
    }
    else if (order==1)
    {
      for (int i=0; i<n; i++)
      {
        varArgDerivs.put(makeMultiSet<int>(i), counter++);
      }
    }
    else if (order==2)
    {
      for (int i=0; i<n; i++)
      {
        for (int j=0; j<=i; j++)
        {
          varArgDerivs.put(makeMultiSet<int>(i,j), counter++);
        }
      }
    }
    else
    {
      TEUCHOS_TEST_FOR_EXCEPTION(order > functorElement()->maxOrder() || order < 0, std::runtime_error,
        "order " << order << " not supported by functor " 
        << functorElement()->masterName());
    }
  }
}

std::ostream& UserDefOpElement::toText(std::ostream& os, bool paren) const 
{
  os << functorElement()->name() << "(";
  for (int i=0; i<numChildren(); i++)
  {
    os << child(i).toString();
    if (i < numChildren()-1) os << ",";
  }
  os << ")";
  return os;
}

XMLObject UserDefOpElement::toXML() const
{
  XMLObject rtn("UserDefOpElement");
  XMLObject args("Arguments");
  for (int i=0; i<numChildren(); i++)
  {
    args.addChild(child(i).toXML());
  }
  rtn.addChild(args);
  rtn.addAttribute("op", functorElement()->name());
  return rtn;
}


Evaluator* UserDefOpElement::createEvaluator(const EvaluatableExpr* expr,
  const EvalContext& context) const
{
  const UserDefOpElement* me = dynamic_cast<const UserDefOpElement*>(expr);
  TEUCHOS_TEST_FOR_EXCEPTION(me == 0, std::logic_error,
    "cast failed in UserDefOpElement::createEvaluator()");
  return new UserDefOpEvaluator(me, getCommonEvaluator(context),
    context);
}



RCP<const UserDefOpCommonEvaluator> 
UserDefOpElement::getCommonEvaluator(const EvalContext& context) const 
{
  if (!commonEvaluatorsMap_->containsKey(context))
  {
    commonEvaluatorsMap_->put(context, 
      rcp(new UserDefOpCommonEvaluator(functorElement()->master(), this, context)));
  }
  return commonEvaluatorsMap_->get(context);
}
