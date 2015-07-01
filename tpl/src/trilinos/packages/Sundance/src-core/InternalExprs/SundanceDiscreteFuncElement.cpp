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

#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"

#include "SundanceDeriv.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


DiscreteFuncElement
::DiscreteFuncElement(const RCP<DiscreteFuncDataStub>& data,
  const std::string& name,
  const std::string& suffix,
  const FunctionIdentifier& fid, int myIndex)
	: EvaluatableExpr(), 
    FuncElementBase(name, suffix, fid),
    commonData_(data),
    miSet_(),
    myIndex_(myIndex)
{}


RCP<Array<Set<MultipleDeriv> > > DiscreteFuncElement
::internalDetermineR(const EvalContext& context,
  const Array<Set<MultipleDeriv> >& RInput) const
{
  Tabs tab(0);
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tab << "DFE::internalDetermineR() for "
    << toString());
  Array<Set<MultipleDeriv> > RIn = RInput;
  Set<MultiIndex> miSet = activeSpatialDerivs(context);

  for (Set<MultiIndex>::const_iterator i=miSet.begin(); i!=miSet.end(); i++)
  {
    const MultiIndex& mi = *i;
    int order = mi.order();
    if (order==0) RIn[0].put(MultipleDeriv());
    if (order==1) RIn[1].put(MultipleDeriv(coordDeriv(mi)));
  }

  return EvaluatableExpr::internalDetermineR(context, RIn);
}


Set<MultipleDeriv> 
DiscreteFuncElement::internalFindW(int order, const EvalContext& context) const
{
  Tabs tab(0);
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tab << "DFE::internalFindW(order=" << order << ") for "
    << toString());
  Set<MultipleDeriv> rtn;

  Set<MultiIndex> miSet = activeSpatialDerivs(context);

  if (order==0) 
  {
    if (miSet.contains(MultiIndex())) rtn.put(MultipleDeriv());
  }
  if (order==1)
  {
    for (Set<MultiIndex>::const_iterator i=miSet.begin(); i!=miSet.end(); i++)
    {
      const MultiIndex& mi = *i;
      int diffOrder = mi.order();
      if (diffOrder==1) 
        rtn.put(MultipleDeriv(coordDeriv(mi)));
    }
  }

  SUNDANCE_MSG3(verb, tab << "W[" << order << "]=" << rtn);
  SUNDANCE_MSG3(verb, tab << "done with DFE::internalFindW(" << order << ") for "
    << toString());
  return rtn;
}

Set<MultipleDeriv> 
DiscreteFuncElement::internalFindV(int order, const EvalContext& context) const
{
  Tabs tab(0);
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tab << "DFE::internalFindV(order=" << order << ") for "
    << toString());
  Set<MultipleDeriv> rtn;
  Set<MultiIndex> miSet = activeSpatialDerivs(context);

  if (order==0) 
  {
    if (miSet.contains(MultiIndex())) rtn.put(MultipleDeriv());
  }
  if (order==1)
  {
    for (Set<MultiIndex>::const_iterator i=miSet.begin(); i!=miSet.end(); i++)
    {
      const MultiIndex& mi = *i;
      int diffOrder = mi.order();
      if (diffOrder==1) 
        rtn.put(MultipleDeriv(coordDeriv(mi)));
    }
  }
  
  rtn = rtn.intersection(findR(order, context));

  SUNDANCE_MSG2(verb, tab << "V[" << order << "]=" << rtn);
  SUNDANCE_MSG2(verb, tab << "done with DFE::internalFindV(" << order << ") for "
    << toString());
  return rtn;
}

Set<MultipleDeriv> 
DiscreteFuncElement::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab(0);
  SUNDANCE_MSG2(context.setupVerbosity(), 
    tab << "DFE::internalFindC is a no-op");
  Set<MultipleDeriv> rtn;
  return rtn;
}

void DiscreteFuncElement::addMultiIndex(const MultiIndex& newMi) const
{
  miSet_.put(newMi);
}

XMLObject DiscreteFuncElement::toXML() const 
{
	XMLObject rtn("DiscreteFuncElement");
	rtn.addAttribute("name", name());
	return rtn;
}


bool DiscreteFuncElement::lessThan(const ScalarExpr* other) const
{
  const DiscreteFuncElement* p 
    = dynamic_cast<const DiscreteFuncElement*>(other);
  TEUCHOS_TEST_FOR_EXCEPT(p==0);

  return fid() < p->fid();
}
