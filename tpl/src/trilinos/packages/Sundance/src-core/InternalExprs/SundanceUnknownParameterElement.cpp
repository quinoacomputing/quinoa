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

#include "SundanceUnknownParameterElement.hpp"


#include "SundanceDerivSet.hpp"
#include "PlayaTabs.hpp"


using namespace Sundance;
using namespace Teuchos;

UnknownParameterElement
::UnknownParameterElement(const std::string& name,
  const std::string& suffix,
  const FunctionIdentifier& fid)
	: UnknownFuncElement(rcp(new UnknownFuncDataStub()), name, 
  suffix, fid),
  SpatiallyConstantExpr()
{}



Set<MultipleDeriv> 
UnknownParameterElement::internalFindW(int order, const EvalContext& context) const
{
  Tabs tab;
  SUNDANCE_MSG2(context.setupVerbosity(), 
    tab << "in UPE::internalFindW, order=" << order);
  Set<MultipleDeriv> rtn;

  if (order==0) 
  {
    if (!evalPtIsZero()) rtn.put(MultipleDeriv());
  }
  else if (order==1)
  {
    MultipleDeriv md = makeMultiDeriv(funcDeriv(this));
    rtn.put(md);
  }

  return rtn;
}



Set<MultipleDeriv> 
UnknownParameterElement::internalFindV(int order, const EvalContext& context) const
{
  Tabs tab;
  SUNDANCE_MSG2(context.setupVerbosity(), 
    tab << "UPE::internalFindV is a no-op");
  Set<MultipleDeriv> rtn;

  return rtn;
}


Set<MultipleDeriv> 
UnknownParameterElement::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab;
  SUNDANCE_MSG2(context.setupVerbosity(), 
    tab << "in UPE::internalFindC, order=" << order);
  Set<MultipleDeriv> rtn;

  if (order==0)
  {
    MultipleDeriv md;
    if (!evalPtIsZero()) rtn.put(md);
  }

  if (order==1)
  {
    MultipleDeriv md = makeMultiDeriv(funcDeriv(this));
    rtn.put(md);
  }
  return rtn.intersection(UnknownFuncElement::findR(order, context));
}



Evaluator* UnknownParameterElement
::createEvaluator(const EvaluatableExpr* expr,
  const EvalContext& context) const 
{
  return SymbolicFuncElement::createEvaluator(expr, context);
}


const Parameter* UnknownParameterElement::parameterValue() const 
{
  const Parameter* p = dynamic_cast<const Parameter*>(evalPt());
  TEUCHOS_TEST_FOR_EXCEPTION(p==0, std::logic_error, 
    "UnknownParameter evalPt() is not a Parameter");
  return p;
}

Parameter* UnknownParameterElement::parameterValue()  
{
  Parameter* p = dynamic_cast<Parameter*>(evalPt());
  TEUCHOS_TEST_FOR_EXCEPTION(p==0, std::logic_error, 
    "UnknownParameter evalPt() is not a Parameter");
  return p;
}


bool UnknownParameterElement::lessThan(const ScalarExpr* other) const
{
  const UnknownParameterElement* p 
    = dynamic_cast<const UnknownParameterElement*>(other);
  TEUCHOS_TEST_FOR_EXCEPT(p==0);

  if (name() < p->name()) return true;

  TEUCHOS_TEST_FOR_EXCEPTION(name()==p->name() && this == p, std::runtime_error,
    "detected two different parameters with the same name");
  return false;
}




XMLObject UnknownParameterElement::toXML() const 
{
  XMLObject rtn("UnknownParameterElement");
  rtn.addAttribute("name", name());
  return rtn;
}

