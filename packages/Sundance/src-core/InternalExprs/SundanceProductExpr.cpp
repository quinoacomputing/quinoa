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


#include "SundanceProductExpr.hpp"
#include "SundanceProductEvaluator.hpp"
#include "SundanceDeriv.hpp"

#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;



ProductExpr::ProductExpr(const RCP<ScalarExpr>& left,
  const RCP<ScalarExpr>& right)
	: BinaryExpr(left, right, 1)
{}


Evaluator* ProductExpr::createEvaluator(const EvaluatableExpr* expr,
  const EvalContext& context) const
{
  return new ProductEvaluator(dynamic_cast<const ProductExpr*>(expr), context);
}

bool ProductExpr::isHungryDiffOp() const
{
  return rightScalar()->isHungryDiffOp();
}



const std::string& ProductExpr::xmlTag() const 
{
	static std::string timesStr = "Times";
	static std::string divideStr = "Divide";
	if (sign() < 0) return divideStr;
	return timesStr;
}

const std::string& ProductExpr::opChar() const 
{
	static std::string timesStr = "*";
	static std::string divideStr = "/";
	if (sign() < 0) return divideStr;
	return timesStr;
}



Set<MultiSet<int> > ProductExpr::internalFindQ_W(int order, const EvalContext& context) const
{
  Tabs tab0(0);
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tab0 << "ProdExpr::internalFindQ_W(" << order << ")");

  Set<MultiSet<int> > rtn;
  if (order > 2) return rtn;

  if (order==2)
  {
    rtn.put(makeMultiSet<int>(0,1));
    return rtn;
  }

  Tabs tab1;
  SUNDANCE_MSG3(verb, tab1 << "calling findW(0) for left");
  const Set<MultipleDeriv>& wLeft 
    = leftEvaluatable()->findW(0, context);
  SUNDANCE_MSG3(verb, tab1 << "found wLeft(0)=" << wLeft);

  SUNDANCE_MSG3(verb, tab1 << "calling findW(0) for right");
  const Set<MultipleDeriv>& wRight
    = rightEvaluatable()->findW(0, context);
  SUNDANCE_MSG3(verb, tab1 << "found wRight(0)=" << wRight);
  
  if (order==0)
  {
    if (wLeft.size() > 0)
    {
      rtn.put(makeMultiSet<int>(0));
    }
    if (wRight.size() > 0)
    {
      rtn.put(makeMultiSet<int>(1));
    }
  }
  
  if (order==1)
  {
    if (wLeft.size() > 0) rtn.put(makeMultiSet<int>(1));
    if (wRight.size() > 0) rtn.put(makeMultiSet<int>(0));
  }
  
  SUNDANCE_MSG2(verb, tab0 << "Q_W[" << order << "]=" << rtn);
  return rtn;
}


Set<MultiSet<int> > ProductExpr::internalFindQ_V(int order, const EvalContext& context) const
{
  Tabs tab0(0);
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tab0 << "ProdExpr::internalFindQ_V(" << order << ")");

  Set<MultiSet<int> > rtn;
  if (order > 1) return rtn;

  const Set<MultipleDeriv>& vLeft 
    = leftEvaluatable()->findV(0, context);
  const Set<MultipleDeriv>& vRight
    = rightEvaluatable()->findV(0, context);

  if (order==0)
  {
    if (vLeft.size() > 0)
    {
      rtn.put(makeMultiSet<int>(0));
    }
    if (vRight.size() > 0)
    {
      rtn.put(makeMultiSet<int>(1));
    }
  }

  if (order==1)
  {
    if (vLeft.size() > 0) rtn.put(makeMultiSet<int>(1));
    if (vRight.size() > 0) rtn.put(makeMultiSet<int>(0));
  }

  SUNDANCE_MSG2(verb, tab0 << "Q_V[" << order << "]=" << rtn);  
  return rtn;
}

