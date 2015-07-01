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

#ifndef SUNDANCE_EVALUATABLEEXPR_H
#define SUNDANCE_EVALUATABLEEXPR_H



#include "SundanceDefs.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceFuncSetAccumulator.hpp"
#include "SundanceMap.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceEvalContext.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;
using namespace Sundance;



class Expr;
class MultipleDeriv;

class EvalManager;
class Evaluator;
class EvaluatorFactory;

/**
 *
 */
enum DerivSubsetSpecifier {AllNonzeros, 
                           RequiredNonzeros,
                           VariableNonzeros,
                           ConstantNonzeros};

/**
 *
 */

class EvaluatableExpr : public virtual ScalarExpr,
                        public virtual EvaluatorFactory,
                        public virtual FuncSetAccumulator,
                        public ObjectWithClassVerbosity<EvaluatableExpr>
{
  typedef OrderedQuartet<EvalContext, 
                         Set<MultiIndex>,
                         Set<MultiSet<int> >,
                         bool> NonzeroSpecifier ;


public:
  /** Ctor is empty, but has some internal initialization to do
   * and so must be called by all subclasses */
  EvaluatableExpr();

  /** virtual dtor */
  virtual ~EvaluatableExpr(){;}


  /** \name Evaluation */
  //@{
  /**
   * Evaluate this expression in the given region, putting the results
   * of the evaluation in the results argument. 
   */
  void evaluate(const EvalManager& mgr,
    Array<double>& constantResults,
    Array<RCP<EvalVector> >& vectorResults) const ;
  //@}

  /** \name Preprocessing */
  //@{
  /**
   * 
   */
  virtual void setupEval(const EvalContext& context) const ;

      


  /** Return the set of all nonzero derivatives
   * required in the given context */
  RCP<SparsitySuperset> sparsitySuperset(const EvalContext& context) const ;

  //@}
      



  /** Utility to downcast an expression to an evaluatable expr. Throws
   * an exception if the cast fails. */
  static const EvaluatableExpr* getEvalExpr(const Expr& expr);

  /** Return the evaluator to be used for the given context */
  const RCP<Evaluator>& evaluator(const EvalContext& context) const; 

  /** */
  virtual void showSparsity(std::ostream& os, 
    const EvalContext& context) const ;


  /** */
  virtual int countNodes() const ;

  /**
   * Find the maximum differentiation order acting on discrete
   * functions in this expression. This is needed to identify
   * expressions where cofacet Jacobians are needed to transform  
   * discrete function derivatives.
   *
   * The base class implementation is a no-op.
   */
  virtual int maxDiffOrderOnDiscreteFunctions() const {return 0;}

  /**
   * Indicate whether this expression contains discrete functions.
   */
  virtual bool hasDiscreteFunctions() const {return false;}

  /** */
  virtual bool nodesHaveBeenCounted() const {return nodesHaveBeenCounted_;}

  /** */
  static int maxFuncDiffOrder() {static int rtn=3; return rtn;}


  /** */
  virtual Set<MultipleDeriv> 
  internalFindW(int order, const EvalContext& context) const = 0 ;

  /** Find spatially-variable functional derivatives. Default implementation
   * returns R */
  virtual Set<MultipleDeriv> 
  internalFindV(int order, const EvalContext& context) const
    {return findR(order, context);}

  /** Find spatially-constant functional derivatives. Default implementation
   * returns the empty set */
  virtual Set<MultipleDeriv> 
  internalFindC(int order, const EvalContext& context) const
    {Set<MultipleDeriv> rtn; return rtn;}

      

  /** */
  const Set<MultipleDeriv>& 
  findDerivSubset(int order,
    const DerivSubsetSpecifier& dss,
    const EvalContext& context) const ;
  /** */
  const Set<MultipleDeriv>& 
  findDerivSubset(const DerivSubsetSpecifier& dss,
    const EvalContext& context) const ;

  /** */
  const Set<MultipleDeriv>& findW(int order, 
    const EvalContext& context) const ;

  /** */
  const Set<MultipleDeriv>& findR(int order, 
    const EvalContext& context) const ;
  /** */
  const Set<MultipleDeriv>& findV(int order, 
    const EvalContext& context) const ;
  /** */
  const Set<MultipleDeriv>& findC(int order, 
    const EvalContext& context) const ;


  /** */
  const Set<MultipleDeriv>& findW(const EvalContext& context) const ;

  /** */
  const Set<MultipleDeriv>& findR(const EvalContext& context) const ;
  /** */
  const Set<MultipleDeriv>& findV(const EvalContext& context) const ;
  /** */
  const Set<MultipleDeriv>& findC(const EvalContext& context) const ;

  /** */
  virtual void displayNonzeros(std::ostream& os, 
    const EvalContext& context) const ;


  /** */
  Set<MultipleDeriv> setProduct(const Set<MultipleDeriv>& a,
    const Set<MultipleDeriv>& b) const ;
  /** */
  Set<MultipleDeriv> setDivision(const Set<MultipleDeriv>& a,
    const Set<MultipleDeriv>& b) const ;
  /** */
  Set<MultiSet<int> > setDivision(const Set<MultiSet<int> >& s,
    int index) const ;
      
  /** */
  void determineR(const EvalContext& context,
    const Array<Set<MultipleDeriv> >& RInput) const ;

  /** */
  virtual RCP<Array<Set<MultipleDeriv> > > 
  internalDetermineR(const EvalContext& context,
    const Array<Set<MultipleDeriv> >& RInput) const ;

  /** */
  const Set<MultipleDeriv>& getR(int order, const EvalContext& context) const ;

      
  /** */
  Array<Set<MultipleDeriv> > 
  computeInputR(const EvalContext& context,
    const Array<Set<MultiSet<int> > >& funcIDCombinations,
    const Array<Set<MultiIndex> >& spatialDerivs) const ;


  /** */
  virtual void registerSpatialDerivs(const EvalContext& context, 
    const Set<MultiIndex>& miSet) const ;
      
  /** */
  static Time& evaluationTimer() ;
      
protected:

  /** Record the evaluator to be used for the given context */
  void registerEvaluator(const EvalContext& context,
    const RCP<Evaluator>& evaluator) const 
    {return evaluators_.put(context, evaluator);}

  /** */
  static bool isEvaluatable(const ExprBase* expr);

  /** */
  Map<EvalContext, RCP<Evaluator> >& evaluators() const 
    {return evaluators_;}


  /** */
  int maxOrder(const Set<MultiIndex>& m) const ;

  /** */
  const Set<MultiIndex>& activeSpatialDerivs(const EvalContext& context) const ;

  /** */
  std::string derivType(const DerivSubsetSpecifier& dss) const;

  /** */
  void printR(int verb, const RCP<Array<Set<MultipleDeriv> > >& R) const ;

private:
      
  /** */
  static int numDerivSubsetTypes() {static int rtn=4; return rtn;}

  /** 
   * evaluators, indexed by context 
   */
  mutable Map<EvalContext, RCP<Evaluator> > evaluators_;

  /** 
   * supersets of nonzero derivatives to be computed, index by
   * context
   */
  mutable Map<EvalContext, RCP<SparsitySuperset> > sparsity_;

  /** Polynomial order of the dependency upon each coordinate direction */
  Array<int> orderOfDependency_;

  /** Set of function combinations appearing in nonzero mixed partials */ 
  Set<MultiSet<int> > funcIDSet_;

  /** */
  Set<int> funcDependencies_;

  mutable Set<NonzeroSpecifier> knownNonzeros_;

  mutable bool nodesHaveBeenCounted_; 

  typedef Array<Map<EvalContext, Set<MultipleDeriv> > > contextToDSSMap_ele_t;
  mutable Array<contextToDSSMap_ele_t> contextToDSSMap_;

  mutable bool rIsReady_;

  mutable Array<Map<EvalContext, Set<MultipleDeriv> > > allOrdersMap_;

  mutable Map<EvalContext, Set<MultiIndex> > activeSpatialDerivMap_;
};
}

#endif
