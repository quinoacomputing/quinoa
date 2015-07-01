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

#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceSymbolicFunc.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceZeroExpr.hpp"

#include "SundanceDerivSet.hpp"
#include "PlayaTabs.hpp"


using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;



SymbolicFuncElement::SymbolicFuncElement(const std::string& name,
  const std::string& suffix,
  const FunctionIdentifier& fid,
  const RCP<const CommonFuncDataStub>& data)
	: EvaluatableExpr(), FuncElementBase(name, suffix, fid),
    commonData_(data),
    evalPt_(),
    evalPtDerivSetIndices_()
{}

void SymbolicFuncElement::registerSpatialDerivs(const EvalContext& context, 
  const Set<MultiIndex>& miSet) const
{
  evalPt()->registerSpatialDerivs(context, miSet);
  EvaluatableExpr::registerSpatialDerivs(context, miSet);
}

void SymbolicFuncElement::accumulateFuncSet(Set<int>& funcDofIDs, 
  const Set<int>& activeSet) const
{
  if (activeSet.contains(fid().dofID())) 
  {
    funcDofIDs.put(fid().dofID());
  }
}


Set<MultipleDeriv> 
SymbolicFuncElement::internalFindW(int order, const EvalContext& context) const
{
  Tabs tab(0);
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tab << "SFE::internalFindW(order=" << order << ") for "
    << toString());

  Set<MultipleDeriv> rtn;

  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab1 << "findW() for eval point");
    evalPt()->findW(order, context);
  }

  if (order==0) 
  {
    Tabs tab1;
    if (!evalPtIsZero()) 
    {
      SUNDANCE_MSG5(verb, tab1 << "value of " << toString() << " is nonzero" );
      rtn.put(MultipleDeriv());
    }
    else
    {
      SUNDANCE_MSG5(verb,  tab1 << "value of " << toString() << " is zero" );
    }
  }
  else if (order==1)
  {
    Deriv d = funcDeriv(this);
    MultipleDeriv md;
    md.put(d);
    rtn.put(md);
  }
  
  SUNDANCE_MSG2(verb,  tab << "SFE: W[" << order << "] = " << rtn );

  return rtn;
}


Set<MultipleDeriv> 
SymbolicFuncElement::internalFindV(int order, const EvalContext& context) const
{
  Tabs tab(0);
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tab << "SFE::internalFindV(order=" << order << ") for "
    << toString());
  Set<MultipleDeriv> rtn;

  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab1 << "findV() for eval point");
    evalPt()->findV(order, context);
  }

  if (order==0) 
  {
    if (!evalPtIsZero()) rtn.put(MultipleDeriv());
  }

  SUNDANCE_MSG5(verb,  tab << "SFE: V = " << rtn );
  SUNDANCE_MSG5(verb,  tab << "SFE: R = " << findR(order, context) );
  rtn = rtn.intersection(findR(order, context));

  SUNDANCE_MSG2(verb,  tab << "SFE: V[" << order << "] = " << rtn );
  return rtn;
}


Set<MultipleDeriv> 
SymbolicFuncElement::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab(0);
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tab << "SFE::internalFindC(order=" << order << ") for "
    << toString());
  Set<MultipleDeriv> rtn;

  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab1 << "findC() for eval point");
    evalPt()->findC(order, context);
  }

  if (order==1)
  {
    Deriv d(funcDeriv(this));
    MultipleDeriv md;
    md.put(d);
    rtn.put(md);
  }

  rtn = rtn.intersection(findR(order, context));

  SUNDANCE_MSG2(verb,  tab << "SFE: C[" << order << "] = " << rtn );
  return rtn;
}


RCP<Array<Set<MultipleDeriv> > > SymbolicFuncElement
::internalDetermineR(const EvalContext& context,
  const Array<Set<MultipleDeriv> >& RInput) const
{
  Tabs tab(0);
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tab << "SFE::internalDetermineR() for "
    << toString());
  {
    Tabs tab1;
    SUNDANCE_MSG3(verb, tab1 << "determineR() for eval point");
    evalPt()->determineR(context, RInput);

    SUNDANCE_MSG3(verb, tab1 << "SFE::internalDetermineR() for "
      << toString() << " delegating to EE");
    return EvaluatableExpr::internalDetermineR(context, RInput);
  }
}

bool SymbolicFuncElement::evalPtIsZero() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(evalPt_.get()==0, std::logic_error,
    "null evaluation point in SymbolicFuncElement::evalPtIsZero()");
  bool isZero = 0 != dynamic_cast<const ZeroExpr*>(evalPt());
  bool isTest = 0 != dynamic_cast<const TestFuncElement*>(this);
  return isZero || isTest;
}

void SymbolicFuncElement::substituteZero() const 
{
  evalPt_ = rcp(new ZeroExpr());
}

void SymbolicFuncElement
::substituteFunction(const RCP<DiscreteFuncElement>& u0) const
{
  evalPt_ = u0;
}



bool SymbolicFuncElement::isIndependentOf(const Expr& u) const 
{
  Expr uf = u.flatten();
  for (int i=0; i<uf.size(); i++)
  {
    const ExprBase* p = uf[i].ptr().get();
    const SymbolicFuncElement* f = dynamic_cast<const SymbolicFuncElement*>(p);
    TEUCHOS_TEST_FOR_EXCEPTION(f==0, std::logic_error, "expected a list of functions, "
      " got " << u);
    if (fid().dofID() == f->fid().dofID()) return false;
  }
  return true;
}


bool SymbolicFuncElement::isLinearForm(const Expr& u) const
{
  Expr uf = u.flatten();
  for (int i=0; i<uf.size(); i++)
  {
    const ExprBase* p = uf[i].ptr().get();
    const SymbolicFuncElement* f = dynamic_cast<const SymbolicFuncElement*>(p);
    TEUCHOS_TEST_FOR_EXCEPTION(f==0, std::logic_error, "expected a list of functions, "
      " got " << u);
    if (fid().dofID() == f->fid().dofID()) return true;
  }
  return false;
}


