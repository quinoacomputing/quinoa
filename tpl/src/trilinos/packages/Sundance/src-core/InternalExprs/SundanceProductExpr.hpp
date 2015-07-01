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



#ifndef SUNDANCE_PRODUCTEXPR_H
#define SUNDANCE_PRODUCTEXPR_H

#include "SundanceBinaryExpr.hpp"


namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;





/** 
 * ProductExpr represents a product of two scalar-valued expressions
 */
class ProductExpr : public BinaryExpr
{
public:
  /** */
  ProductExpr(const RCP<ScalarExpr>& a, 
    const RCP<ScalarExpr>& b);

  /** virtual dtor */
  virtual ~ProductExpr() {;}

  /** Indicate whether this expression is a "hungry"
   * differential operator that is awaiting an argument. */
  virtual bool isHungryDiffOp() const ;

  /** */
  virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
    const EvalContext& context) const ;
      
  /** */
  virtual Set<MultiSet<int> > internalFindQ_W(int order, 
    const EvalContext& context) const ;
      
  /** */
  virtual Set<MultiSet<int> > internalFindQ_V(int order, 
    const EvalContext& context) const ;

      

  /** */
  virtual bool isProduct() const {return true;}

  /** 
   * Indicate whether the expression is linear 
   * with respect to test functions */
  virtual bool isLinearInTests() const 
    {
      bool leftIsLinear = leftScalar()->isLinearInTests() ;
      bool rightIsLinear = rightScalar()->isLinearInTests() ;

      bool leftHasTests = leftScalar()->hasTestFunctions() ;
      bool rightHasTests = rightScalar()->hasTestFunctions() ;
          
      return (leftIsLinear && !rightHasTests) 
        || (!leftHasTests && rightIsLinear);
    }
      
  /** Indicate whether the expression is linear in the given 
   * functions */
  virtual bool isLinearForm(const Expr& u) const 
    {
      bool L = leftScalar()->isLinearForm(u);
      bool R = rightScalar()->isLinearForm(u);
      return ( (L || R) && !(L==R) );
    }
      
  /** Indicate whether the expression is a quadratic form in the given 
   * functions */
  virtual bool isQuadraticForm(const Expr& u) const
    {
      bool LQ = leftScalar()->isQuadraticForm(u);
      bool RQ = rightScalar()->isQuadraticForm(u);
      bool LL = leftScalar()->isLinearForm(u);
      bool RL = leftScalar()->isLinearForm(u);
      bool LI = leftScalar()->isIndependentOf(u);
      bool RI = rightScalar()->isIndependentOf(u);

      if (LI && RQ) return true;
      if (RI && LQ) return true;
      if (LL && RL) return true;
      return false;
    }

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}
protected:
  /** */
  virtual bool parenthesizeSelf() const {return true;}
  /** */
  virtual bool parenthesizeOperands() const {return true;}
  /** */
  virtual const std::string& xmlTag() const ;
  /** */
  virtual const std::string& opChar() const ;

private:

};
}

#endif
