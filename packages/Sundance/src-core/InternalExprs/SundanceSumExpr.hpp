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

#ifndef SUNDANCE_SUMEXPR_H
#define SUNDANCE_SUMEXPR_H

#include "SundanceBinaryExpr.hpp"
#include "SundanceSumEvaluator.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;



/**
 * SumExpr is the internal representation of an addition or subtraction
 * node in the expression tree.
 */
class SumExpr : public BinaryExpr,
                public GenericEvaluatorFactory<SumExpr, SumEvaluator>
{
public:
  /** */
  SumExpr(const RCP<ScalarExpr>& a, 
    const RCP<ScalarExpr>& b, int sign);

  /** virtual dtor */
  virtual ~SumExpr() {;}

  /** */
  virtual bool isHungryDiffOp() const ;

  /** */
  virtual bool isLinear() const {return true;}

          

  /** 
   * Indicate whether the expression is linear 
   * with respect to test functions */
  virtual bool isLinearInTests() const ;

  /** 
   * Indicate whether every term in the expression contains test functions */
  virtual bool everyTermHasTestFunctions() const ;
          

  /** Indicate whether the expression is linear in the given 
   * functions */
  virtual bool isLinearForm(const Expr& u) const ;

  /** Indicate whether the expression is a 
   * quadratic form in the given functions */
  virtual bool isQuadraticForm(const Expr& u) const ;

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}
          
  /** */
  const Map<Expr, int>& getSumTree() const {return sumTree_;}

protected:
  /** */
  virtual bool parenthesizeSelf() const {return true;}
  /** */
  virtual bool parenthesizeOperands() const {return false;}
  /** */
  virtual const std::string& xmlTag() const ;
  /** */
  virtual const std::string& opChar() const ;

private:
  Map<Expr, int> sumTree_;


};
}
#endif
