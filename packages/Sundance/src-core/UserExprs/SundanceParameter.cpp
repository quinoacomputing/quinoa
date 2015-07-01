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

#include "SundanceParameter.hpp"
#include "SundanceSymbolicFunc.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;

Parameter::Parameter(const double& value, const std::string& name)
	: DiscreteFuncElement(rcp(new ParameterData(value)), 
    name, "", makeFuncID(0), 0),
    SpatiallyConstantExpr()
{;}

void Parameter::setValue(const double& value)
{
  data()->setValue(value);
}

const double& Parameter::value() const {return data()->value();}


RCP<Array<Set<MultipleDeriv> > > Parameter
::internalDetermineR(const EvalContext& context,
                     const Array<Set<MultipleDeriv> >& RInput) const
{
  Tabs tab;
  SUNDANCE_VERB_HIGH(tab << "Parameter::internalDetermineR() for "
                     << toString());
  return EvaluatableExpr::internalDetermineR(context, RInput);
}

Set<MultipleDeriv> 
Parameter::internalFindW(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;

  if (order==0) rtn.put(MultipleDeriv());

  return rtn;
}

Set<MultipleDeriv> 
Parameter::internalFindV(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;

  return rtn;
}


Set<MultipleDeriv> 
Parameter::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG3(verb, tab 
    << "Parameter::internalFindC(order=" << order << ") for "
    << toString());
  Set<MultipleDeriv> rtn;

  if (order==0)
    {
      rtn.put(MultipleDeriv());
    }

  rtn = rtn.intersection(findR(order, context));

  SUNDANCE_MSG3(verb,  tab << "Parameter: C[" << order << "] = " << rtn );
  return rtn;
}


Evaluator* Parameter::createEvaluator(const EvaluatableExpr* expr,
                                      const EvalContext& context) const 
{
  return new ConstantEvaluator(this, context);
}


XMLObject Parameter::toXML() const 
{
	XMLObject rtn("Parameter");
	rtn.addAttribute("name", name());
	rtn.addAttribute("value", Teuchos::toString(value()));
	return rtn;
}


const ParameterData* Parameter::data() const
{
  const ParameterData* pd = dynamic_cast<const ParameterData*>(commonData().get());

  TEUCHOS_TEST_FOR_EXCEPTION(pd==0, std::logic_error, "cast failed");

  return pd;
}



ParameterData* Parameter::data()
{
  ParameterData* pd = dynamic_cast<ParameterData*>(commonData());

  TEUCHOS_TEST_FOR_EXCEPTION(pd==0, std::logic_error, "cast failed");

  return pd;
}
