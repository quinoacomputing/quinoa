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

#include "SundanceConstantExpr.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;

ConstantExpr::ConstantExpr(const double& value)
	: SpatiallyConstantExpr(), value_(value)
{}



Set<MultipleDeriv> 
ConstantExpr::internalFindW(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;

  if (order==0) rtn.put(MultipleDeriv());
  Tabs tab;
  SUNDANCE_MSG2(context.setupVerbosity(), 
    tab << "ConstantExpr::internalFindW found" << rtn << " for order="
    << order);

  return rtn;
}

Set<MultipleDeriv> 
ConstantExpr::internalFindV(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;
  Tabs tab;
  SUNDANCE_MSG2(context.setupVerbosity(), 
    tab << "ConstantExpr::internalFindV is a no-op");
  return rtn;
}


Set<MultipleDeriv> 
ConstantExpr::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab;
  SUNDANCE_MSG2(context.setupVerbosity(), 
    tab << "ConstantExpr::internalFindC is forwarding to findR()");
  return findR(order, context);
}


bool ConstantExpr::lessThan(const ScalarExpr* other) const
{
  const ConstantExpr* c = dynamic_cast<const ConstantExpr*>(other);
  TEUCHOS_TEST_FOR_EXCEPTION(c==0, std::logic_error, "cast should never fail at this point");
  return value() < c->value();
}


std::ostream& ConstantExpr::toText(std::ostream& os, bool /* paren */) const 
{
	os << value();
	return os;
}


XMLObject ConstantExpr::toXML() const 
{
	XMLObject rtn("Constant");
	rtn.addAttribute("value", Teuchos::toString(value()));
	return rtn;
}


