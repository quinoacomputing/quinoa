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

#include "SundanceSumExpr.hpp"
#include "SundanceExpr.hpp"
#include "PlayaTabs.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceOut.hpp"



using namespace Sundance;
using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;


SumExpr::SumExpr(const RCP<ScalarExpr>& left,
  const RCP<ScalarExpr>& right, int sign)
	: BinaryExpr(left, right, sign), sumTree_()
{
  /*
    Expr L = Expr::handle(left);
    Expr R = Expr::handle(right);

    sumTree_ = L.getSumTree();
    Map<Expr, int> rightTree = R.getSumTree();

    for (Map<Expr, int>::const_iterator i=rightTree.begin(); i!=rightTree.end(); i++)
    {
    int leftCount = 0;
    if (sumTree_.containsKey(i->first))
    {
    leftCount = sumTree_[i->first];
    }
    int rightCount = sign * i->second;
    sumTree_.put(i->first, leftCount + rightCount);
    }
  */
}

bool SumExpr::isHungryDiffOp() const
{
  return leftScalar()->isHungryDiffOp() || rightScalar()->isHungryDiffOp();
}


const std::string& SumExpr::xmlTag() const 
{
	static std::string plusStr = "Plus";
	static std::string minusStr = "Minus";
	if (sign() < 0) return minusStr;
	return plusStr;
}

const std::string& SumExpr::opChar() const 
{
	static std::string plusStr = "+";
	static std::string minusStr = "-";
	if (sign() < 0) return minusStr;
	return plusStr;
}


bool SumExpr::everyTermHasTestFunctions() const
{
  return leftEvaluatable()->everyTermHasTestFunctions()
    && rightEvaluatable()->everyTermHasTestFunctions();
}

bool SumExpr::isLinearInTests() const
{
  bool leftHasTests = leftScalar()->hasTestFunctions();
  bool rightHasTests = rightScalar()->hasTestFunctions();

  bool leftIsLinear = leftScalar()->isLinearInTests();
  bool rightIsLinear = rightScalar()->isLinearInTests();

  return (!leftHasTests || leftIsLinear) && (!rightHasTests || rightIsLinear);
}


bool SumExpr::isLinearForm(const Expr& u) const 
{
  bool LL = leftScalar()->isLinearForm(u);
  bool RL = rightScalar()->isLinearForm(u);
  bool LI = leftScalar()->isIndependentOf(u);
  bool RI = rightScalar()->isIndependentOf(u);

  return ( (LL && (RL || RI)) || (RL && (LL || LI)) );
}

bool SumExpr::isQuadraticForm(const Expr& u) const
{
  bool LQ = leftScalar()->isQuadraticForm(u);
  bool RQ = rightScalar()->isQuadraticForm(u);
  bool LL = leftScalar()->isLinearForm(u);
  bool RL = rightScalar()->isLinearForm(u);
  bool LI = leftScalar()->isIndependentOf(u);
  bool RI = rightScalar()->isIndependentOf(u);

  return ( (LQ && (RQ || RL || RI)) || (RQ && (LQ || LL || LI))); 
}
