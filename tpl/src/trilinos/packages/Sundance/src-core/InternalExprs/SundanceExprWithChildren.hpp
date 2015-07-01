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

#ifndef SUNDANCE_EXPRWITHCHILDREN_H
#define SUNDANCE_EXPRWITHCHILDREN_H

#include "SundanceDefs.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceCombinatorialUtils.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceExpr.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;


/** 
 * ExprWithChildren is a base class for any evaluatable expression
 * that has child nodes, for example, sums and unary operators.
 * ExprWithChildren adds nothing new to the expr interface, but 
 * provides some common utilities for getting children
 * and recursing to children.
 */
class ExprWithChildren : public virtual EvaluatableExpr
{
public:
  /** construct with a list of child operands */
  ExprWithChildren(const Array<RCP<ScalarExpr> >& children);

  /** virtual dtor */
  virtual ~ExprWithChildren() {;}

  /**
   * Do preprocessing to set up sparse evaluation in the given region 
   */
  virtual void setupEval(const EvalContext& context) const ;

  /** Determine whether this expression is constant. It will
   * be constant if all children are constant. */
  virtual bool isConstant() const ;

  /** Return the number of children */
  int numChildren() const {return children_.size();}
          
  /** downcast the i-th to an evaluatable expr */
  const EvaluatableExpr* evaluatableChild(int i) const ;

  /** downcast the i-th to a scalar expr */
  const ScalarExpr* scalarChild(int i) const 
    {return children_[i].get();}

  /** Get a handle to the i-th child */
  Expr child(int i) const {return Expr::handle(children_[i]);}


  /**
   * Find the maximum differentiation order acting on discrete
   * functions in this expression. 
   */
  virtual int maxDiffOrderOnDiscreteFunctions() const ;

  /**
   * Indicate whether this expression contains discrete functions.
   */
  virtual bool hasDiscreteFunctions() const ;

  /** */
  virtual Set<MultipleDeriv> 
  internalFindW(int order, const EvalContext& context) const ;

  /** */
  virtual Set<MultipleDeriv> 
  internalFindC(int order, const EvalContext& context) const ;
          
  /** */
  virtual Set<MultipleDeriv> 
  internalFindV(int order, const EvalContext& context) const ;

  /** */
  virtual void displayNonzeros(std::ostream& os, 
    const EvalContext& context) const ;

  /** */
  const Set<MultiSet<int> >& findQ_W(int order, 
    const EvalContext& context) const ;

  /** */
  const Set<MultiSet<int> >& findQ_C(int order, 
    const EvalContext& context) const ;

  /** */
  const Set<MultiSet<int> >& findQ_V(int order, 
    const EvalContext& context) const ;

  /** */
  virtual Set<MultiSet<int> > 
  internalFindQ_W(int order, 
    const EvalContext& context) const ;

  /** */
  virtual Set<MultiSet<int> > 
  internalFindQ_V(int order, 
    const EvalContext& context) const ;

  /** */
  virtual Set<MultiSet<int> > 
  internalFindQ_C(int order, 
    const EvalContext& context) const ;

  /** */
  const Set<MultiSet<int> >& getI_N() const ;

  /** */
  Set<MultiSet<int> > indexSetProduct(const Set<MultiSet<int> >& a,
    const Set<MultiSet<int> >& b) const ;
          
  /** Return true if any child returns true. The sum expression
   * will override this requiring all children to return true */
  virtual bool everyTermHasTestFunctions() const ;

  /** Test whether this expr contains a test function. 
   * If any child contains a test, return true. */
  virtual bool hasTestFunctions() const ;

  /** Test whether this expr contains an unknown function. 
   * If any child contains an unknown, return true. */
  virtual bool hasUnkFunctions() const ;

  /** */
  virtual void showSparsity(std::ostream& os, 
    const EvalContext& context) const ;

  /** */
  virtual void getUnknowns(Set<int>& unkID, Array<Expr>& unks) const ;

  /** */
  virtual void getTests(Set<int>& varID, Array<Expr>& vars) const ;

          
  /** */
  virtual int countNodes() const ;

  /** */
  virtual bool isLinear() const {return false;}

  /** */
  virtual bool isProduct() const {return false;}

  /** Indicate whether the expression is independent of the given 
   * functions */
  virtual bool isIndependentOf(const Expr& u) const ;


  /** */
  Set<MultiSet<int> > subsetContainingIndex(const Set<MultiSet<int> >& s,
    int index) const ;

  /** */
  virtual RCP<Array<Set<MultipleDeriv> > > 
  internalDetermineR(const EvalContext& context,
    const Array<Set<MultipleDeriv> >& RInput) const ;

  /** Determine whether the given child is needed to compute derivatives
   * of the given order */
  bool childIsRequired(int childIndex, int diffOrder,
    const EvalContext& context) const ;



  /** */
  Set<MultipleDeriv> product(const Array<int>& J, const Array<int>& K,
    DerivSubsetSpecifier dss,
    const EvalContext& context) const ;

  /** Append to the set of func IDs present in this expression. */
  virtual void accumulateFuncSet(Set<int>& funcIDs, 
    const Set<int>& activeSet) const ;

  /** */
  virtual void registerSpatialDerivs(const EvalContext& context, 
    const Set<MultiIndex>& miSet) const ;

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;

private:
  Array<RCP<ScalarExpr> > children_;

  static Map<int, Set<MultiSet<int> > >& cachedI_N()
    {static Map<int, Set<MultiSet<int> > > rtn; return rtn;}

  mutable Array<Map<EvalContext, Set<MultiSet<int> > > > contextToQWMap_;
  mutable Array<Map<EvalContext, Set<MultiSet<int> > > > contextToQVMap_;
  mutable Array<Map<EvalContext, Set<MultiSet<int> > > > contextToQCMap_;
};          

      

/** \relates ExprWithChildren */
Array<Array<std::pair<int, Array<MultipleDeriv> > > > chainRuleDerivsOfArgs(int nArgs,
  const MultiSet<int>& bSet,
  const MultipleDeriv& c);

/**  \relates ExprWithChildren */
Array<Array<Array<int> > > bStructure(const Array<int>& b,
  const Array<Array<int> >& tmp);

/**  \relates ExprWithChildren 
 * Return the set of (k,l) tuples appearing in the Constantine
 * and Savits formulation of the multivariable, multiargument
 * chain rule. 
 * \param s 
 * \param lambda
 * \param nu
 */
Array<OrderedPair<Array<MultiSet<int> >, Array<MultipleDeriv> > >
chainRuleTerms(int s, 
  const MultiSet<int>& lambda,
  const MultipleDeriv& nu) ;
      

/** Return all subsets of a multiset. */
Set<MultipleDeriv> multisetSubsets(const MultipleDeriv& x);

/** Return the multiplicity of a chain rule term */
int chainRuleMultiplicity(const MultipleDeriv& nu,
  const Array<MultiSet<int> >& K,
  const Array<MultipleDeriv>& L);
}

#endif
