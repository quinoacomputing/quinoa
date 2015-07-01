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

#ifndef SUNDANCE_DIFFOP_H
#define SUNDANCE_DIFFOP_H

#include "SundanceDefs.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceDiffOpEvaluator.hpp"


namespace Sundance
{
class Expr;

using namespace Sundance;
using namespace Teuchos;




/**
 *
 */
class DiffOp : public UnaryExpr
{
public:
  /** ctor */
  DiffOp(const MultiIndex& op, const RCP<ScalarExpr>& arg);

  /** virtual destructor */
  virtual ~DiffOp() {;}


  /** 
   * Indicate whether the expression is linear 
   * with respect to test functions 
   */
  virtual bool isLinearInTests() const 
    {return evaluatableArg()->isLinearInTests();}
      
      
  /** Indicate whether the expression is linear in the given 
   * functions */
  virtual bool isLinearForm(const Expr& u) const 
    {
      return evaluatableArg()->isLinearForm(u);
    }
      
  /** Indicate whether the expression is at most 
   * quadratic in the given functions */
  virtual bool isQuadraticForm(const Expr& u) const
    {
      return evaluatableArg()->isQuadraticForm(u);
    }

  /**
   * Find the maximum differentiation order acting on discrete
   * functions in this expression. 
   */
  virtual int maxDiffOrderOnDiscreteFunctions() const 
    {
      int rtn = evaluatableArg()->maxDiffOrderOnDiscreteFunctions();
      if (evaluatableArg()->hasDiscreteFunctions()) 
      {
        rtn += mi_.order();
      }
      return rtn;
    }


  /** Write a simple text description suitable
   * for output to a terminal */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;

  /** Write in XML */
  virtual XMLObject toXML() const ;


  /** */
  virtual Set<MultipleDeriv> 
  internalFindW(int order, const EvalContext& context) const ;


  /** */
  virtual Set<MultipleDeriv> 
  internalFindV(int order, const EvalContext& context) const ;


  /** */
  virtual Set<MultipleDeriv> 
  internalFindC(int order, const EvalContext& context) const ;


  /** */
  virtual RCP<Array<Set<MultipleDeriv> > > 
  internalDetermineR(const EvalContext& context,
    const Array<Set<MultipleDeriv> >& RInput) const ;

  /** */
  void requestMultiIndexAtEvalPoint(const MultiIndex& mi,
    const MultipleDeriv& u,
    const EvalContext& context) const ;
      
      
      
  /** */
  const Deriv& myCoordDeriv() const {return myCoordDeriv_;}

  /** */
  const MultiIndex& mi() const {return mi_;}

  /** Get the functions that are required in the evaluation
   * of the multiple deriv d */
  const Sundance::Set<Deriv>& requiredFunctions(const MultipleDeriv& d) const 
    {return requiredFunctions_[d];}

  /** */
  bool requiresFunctionsToEval(const MultipleDeriv& d) const 
    {return requiredFunctions_.containsKey(d);}

    
      
     

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}


  /** */
  virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
    const EvalContext& context) const ;

     

  /** */
  virtual void registerSpatialDerivs(const EvalContext& context, 
    const Set<MultiIndex>& miSet) const ;

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;

private:

      
      

  MultiIndex mi_;


  Deriv myCoordDeriv_;

  mutable Map<MultipleDeriv, Sundance::Set<Deriv>, 
              increasingOrder<MultipleDeriv> > requiredFunctions_;

  mutable bool ignoreFuncTerms_;
};
}


#endif
