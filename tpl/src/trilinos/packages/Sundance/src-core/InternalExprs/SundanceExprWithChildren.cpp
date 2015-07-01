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

#include "SundanceExprWithChildren.hpp"
#include "SundanceCombinatorialUtils.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceExpr.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceNullEvaluator.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceUnaryExpr.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;




ExprWithChildren::ExprWithChildren(const Array<RCP<ScalarExpr> >& children)
	: EvaluatableExpr(), 
    children_(children),
    contextToQWMap_(4),
    contextToQVMap_(4),
    contextToQCMap_(4)
{}

bool ExprWithChildren::lessThan(const ScalarExpr* other) const
{
  const ExprWithChildren* c = dynamic_cast<const ExprWithChildren*>(other);
  TEUCHOS_TEST_FOR_EXCEPTION(c==0, std::logic_error, "cast should never fail at this point");

  if (children_.size() < c->children_.size()) return true;
  if (children_.size() > c->children_.size()) return false;
  
  for (int i=0; i<children_.size(); i++)
  {
    Expr me = Expr::handle(children_[i]);
    Expr you = Expr::handle(c->children_[i]);
    if (me.lessThan(you)) return true;
    if (you.lessThan(me)) return false;
  }
  return false;
}



bool ExprWithChildren::isConstant() const
{
  for (int i=0; i<children_.size(); i++) 
  {
    if (!children_[i]->isConstant()) return false;
  }
  return true;
}

bool ExprWithChildren::isIndependentOf(const Expr& u) const
{
  for (int i=0; i<children_.size(); i++) 
  {
    if (!children_[i]->isIndependentOf(u)) return false;
  }
  return true;
}

const EvaluatableExpr* ExprWithChildren::evaluatableChild(int i) const
{
  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(children_[i].get());

  TEUCHOS_TEST_FOR_EXCEPTION(e==0, std::logic_error, 
    "ExprWithChildren: cast of child [" 
    << children_[i]->toString()
    << " to evaluatable expr failed");

  return e;
}

int ExprWithChildren::maxDiffOrderOnDiscreteFunctions() const
{
  int biggest = -1;
  for (int i=0; i<numChildren(); i++)
  {
    int x = evaluatableChild(i)->maxDiffOrderOnDiscreteFunctions();
    if (x > biggest) biggest = x;
  }
  return biggest;
}

bool ExprWithChildren::hasDiscreteFunctions() const
{
  for (int i=0; i<numChildren(); i++)
  {
    if (evaluatableChild(i)->hasDiscreteFunctions()) return true;
  }
  return false;
}


void ExprWithChildren::accumulateFuncSet(Set<int>& funcIDs, 
  const Set<int>& activeSet) const
{
  for (int i=0; i<children_.size(); i++) 
  {
    children_[i]->accumulateFuncSet(funcIDs, activeSet);
  }
}

void ExprWithChildren::registerSpatialDerivs(const EvalContext& context, 
  const Set<MultiIndex>& miSet) const
{
  for (int i=0; i<children_.size(); i++) 
  {
    evaluatableChild(i)->registerSpatialDerivs(context, miSet);
  }
  EvaluatableExpr::registerSpatialDerivs(context, miSet);
}

void ExprWithChildren::setupEval(const EvalContext& context) const
{
  Tabs tabs;
  int verb = context.evalSetupVerbosity();
  SUNDANCE_MSG1(verb, tabs << "expr " + toString() 
    << ": creating evaluators for children");
  RCP<SparsitySuperset> ss = sparsitySuperset(context);
  SUNDANCE_MSG2(verb, tabs << "my sparsity superset = " 
    << std::endl << *ss);


  /* Create the evaluators for the children first so that I can refer to
   * them when I create my own evaluator */
  if (ss->numDerivs() > 0)
  {
    for (int i=0; i<children_.size(); i++)
    {
      Tabs tabs1;
      SUNDANCE_MSG1(verb, tabs1 << "creating evaluator for child " 
        << evaluatableChild(i)->toString());
      evaluatableChild(i)->setupEval(context);
    }
  }

  if (!evaluators().containsKey(context))
  {
    RCP<Evaluator> eval;
    if (ss->numDerivs()>0)
    {
      eval = rcp(createEvaluator(this, context));
    }
    else
    {
      SUNDANCE_MSG2(verb, 
        tabs << "EWC: no results needed... creating null evaluator");
      eval = rcp(new NullEvaluator());
    }
    evaluators().put(context, eval);
  }
}

void ExprWithChildren::showSparsity(std::ostream& os, 
  const EvalContext& context) const
{
  Tabs tab0;
  os << tab0 << "Node: " << toString() << std::endl;
  sparsitySuperset(context)->print(os);
  for (int i=0; i<children_.size(); i++)
  {
    Tabs tab1;
    os << tab1 << "Child " << i << std::endl;
    evaluatableChild(i)->showSparsity(os, context);
  }
}


bool ExprWithChildren::everyTermHasTestFunctions() const
{
  for (int i=0; i<children_.size(); i++)
  {
    if (evaluatableChild(i)->everyTermHasTestFunctions()) return true;
  }
  return false;
}


bool ExprWithChildren::hasTestFunctions() const
{
  for (int i=0; i<children_.size(); i++)
  {
    if (evaluatableChild(i)->hasTestFunctions()) return true;
  }
  return false;
}


bool ExprWithChildren::hasUnkFunctions() const
{
  for (int i=0; i<children_.size(); i++)
  {
    if (evaluatableChild(i)->hasUnkFunctions()) return true;
  }
  return false;
}


void ExprWithChildren::getUnknowns(Set<int>& unkID, Array<Expr>& unks) const
{
  for (int i=0; i<children_.size(); i++)
  {
    const RCP<ExprBase>& e = children_[i];
    const UnknownFuncElement* u 
      = dynamic_cast<const UnknownFuncElement*>(e.get());
    if (u != 0)
    {
      Expr expr(e);
      if (!unkID.contains(u->fid().dofID())) 
      {
        unks.append(expr);
        unkID.put(u->fid().dofID());
      }
    }
    else
    {
      evaluatableChild(i)->getUnknowns(unkID, unks);
    }
  }
}

void ExprWithChildren::getTests(Set<int>& varID, Array<Expr>& vars) const
{
  for (int i=0; i<children_.size(); i++)
  {
    const RCP<ExprBase>& e = children_[i];
    const TestFuncElement* u 
      = dynamic_cast<const TestFuncElement*>(e.get());
    if (u != 0)
    {
      Expr expr(e);
      if (!varID.contains(u->fid().dofID())) 
      {
        vars.append(expr);
        varID.put(u->fid().dofID());
      }

    }
    else
    {
      evaluatableChild(i)->getTests(varID, vars);
    }
  }
}


int ExprWithChildren::countNodes() const
{
  if (nodesHaveBeenCounted()) 
  {
    return 0;
  }

  /* count self */
  int count = EvaluatableExpr::countNodes();

  /* count children */
  for (int i=0; i<children_.size(); i++)
  {
    if (!evaluatableChild(i)->nodesHaveBeenCounted())
    {
      count += evaluatableChild(i)->countNodes();
    }
  }
  return count;
}

Set<MultipleDeriv> 
ExprWithChildren::product(const Array<int>& J, const Array<int>& K,
  DerivSubsetSpecifier dss,
  const EvalContext& context) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(J.size() != K.size(), std::logic_error,
    "mismatched index set sizes");
  TEUCHOS_TEST_FOR_EXCEPTION(J.size() == 0, std::logic_error,
    "empty index set");

  Set<MultipleDeriv> rtn 
    = evaluatableChild(J[0])->findDerivSubset(K[0], dss, context);
  
  for (int i=1; i<J.size(); i++)
  {
    const Set<MultipleDeriv>& S 
      = evaluatableChild(J[i])->findDerivSubset(K[i], dss, context);
    rtn = setProduct(rtn, S);
  }

  return rtn;
}

Set<MultipleDeriv> 
ExprWithChildren::internalFindV(int order, const EvalContext& context) const
{
  Tabs tab0;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tab0 << "EWC::internalFindV(" << order 
    << ") for " << toString());

  Set<MultipleDeriv> rtn;



  /* we'll dealt with zero order derivatives specially */
  if (order==0) 
  {
    Tabs tab1;
    SUNDANCE_MSG3(verb, tab1 << "case: order=0"); 
    for (int i=0; i<numChildren(); i++)
    {
      if (!childIsRequired(i, order, context)) continue;
      const Set<MultipleDeriv>& childV0 
        = evaluatableChild(i)->findV(0, context);        
      rtn.merge(childV0);
    }
    const Set<MultipleDeriv>& R0 = findR(0, context);
    rtn = rtn.intersection(R0);

    SUNDANCE_MSG3(verb, tab0 << "V[" << order << "]=" << rtn);
    SUNDANCE_MSG3(verb, tab0 << "done with EWC::internalFindV(" << order
      << ") for " << toString());
    return rtn;
  }

  /* now do arbitrary order derivatives with the multivariable chain rule*/
  SUNDANCE_MSG3(verb, tab0 << "getting required set"); 
  const Set<MultipleDeriv>& RM = findR(order, context);
  SUNDANCE_MSG3(verb, tab0 << "RM=" << RM); 
  Array<Array<Array<int> > > comp = compositions(order);
  for (int i=1; i<=order; i++) 
  {
    Tabs tab1;
    SUNDANCE_MSG3(verb, tab1 << "i=" << i << " out of order=" << order);
    const Set<MultiSet<int> >& QC = findQ_C(i, context);
    SUNDANCE_MSG3(verb, tab1 << "QC[" << i << "] = " << QC);
    for (Set<MultiSet<int> >::const_iterator j=QC.begin(); j!=QC.end(); j++)
    {
      Tabs tab2;
      Array<int> J = j->elements();
      const Array<Array<int> >& K = comp[J.size()-1];
      SUNDANCE_MSG3(verb, tab2 << "J=" << J);
      SUNDANCE_MSG3(verb, tab2 << "K=" << K);

      for (int k=0; k<K.size(); k++)
      {
        Tabs tab3;
        Set<MultipleDeriv> RProd = product(J, K[k], 
          RequiredNonzeros, context);
        Set<MultipleDeriv> CProd = product(J, K[k], 
          ConstantNonzeros, context);
        SUNDANCE_MSG3(verb, tab3 << "CProd = " << CProd);
        SUNDANCE_MSG3(verb, tab3 << "RProd = " << RProd);
        rtn.merge(RProd.setDifference(CProd));
      }
    }

    const Set<MultiSet<int> >& QV = findQ_V(i, context);
    SUNDANCE_MSG3(verb, tab1 << "QV[" << i << "] = " << QV);
    for (Set<MultiSet<int> >::const_iterator j=QV.begin(); j!=QV.end(); j++)
    {
      Tabs tab2;
      Array<int> J = j->elements();
      const Array<Array<int> >& K = comp[J.size()-1];
      SUNDANCE_MSG3(verb, tab2 << "J=" << J);
      SUNDANCE_MSG3(verb, tab2 << "K=" << K);

      for (int k=0; k<K.size(); k++)
      {
        Tabs tab3;
        Set<MultipleDeriv> RProd = product(J, K[k], 
          RequiredNonzeros, context);
        SUNDANCE_MSG3(verb, tab3 << "RProd = " << RProd);
        rtn.merge(RProd);
      }
    }
  }

  SUNDANCE_MSG3(verb, tab0 << "rtn before intersection = " << rtn);
  rtn = rtn.intersection(RM);
  SUNDANCE_MSG2(verb, tab0 << "V[" << order << "]=" << rtn);
  SUNDANCE_MSG2(verb, tab0 << "done with EWC::internalFindV(" << order
    << ") for " << toString());

  return rtn;
}


Set<MultipleDeriv> 
ExprWithChildren::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab0;
  Set<MultipleDeriv> rtn;
  bool started = false;
  int verb = context.setupVerbosity();

  SUNDANCE_MSG2(verb, tab0 << "EWC::internalFindC(" << order 
    << ") for " << toString());

  /* we'll dealt with zero order derivatives specially */
  if (order==0) 
  {
    for (int i=0; i<numChildren(); i++)
    {
      if (!childIsRequired(i, 0, context)) continue;
      const Set<MultipleDeriv>& childV0 
        = evaluatableChild(i)->findV(0, context);        
      rtn.merge(childV0);
    }
    const Set<MultipleDeriv>& R0 = findR(0, context);
    rtn = R0.setDifference(rtn);
    SUNDANCE_MSG3(verb, tab0 << "C[" << order << "]=" << rtn);
    SUNDANCE_MSG3(verb, tab0 << "done with EWC::internalFindC(" << order
      << ") for " << toString());
    return rtn;
  }


  /* now do arbitrary order derivatives with the multivariable chain rule*/
  Array<Array<Array<int> > > comp = compositions(order);
  const Set<MultipleDeriv>& RM = findR(order, context);
  for (int i=1; i<=order; i++) 
  {
    Tabs tab1;
    SUNDANCE_MSG3(verb, tab1 << "i=" << i << " up to order=" << order);


    const Set<MultiSet<int> >& QC = findQ_C(i, context);
    SUNDANCE_MSG5(verb, tab1 << "finding CProd union (R\\RProd) over QC");      
    SUNDANCE_MSG5(verb, tab1 << "QC = " << QC);


    for (Set<MultiSet<int> >::const_iterator j=QC.begin(); j!=QC.end(); j++)
    {
      Tabs tab2;
      Array<int> J = j->elements();
      const Array<Array<int> >& K = comp[J.size()-1];
      SUNDANCE_MSG5(verb, tab2 << "J=" << J);
      SUNDANCE_MSG5(verb, tab2 << "K=" << K);

      for (int k=0; k<K.size(); k++)
      {
        Tabs tab3;
        Set<MultipleDeriv> RProd = product(J, K[k], 
          RequiredNonzeros, context);
        Set<MultipleDeriv> CProd = product(J, K[k], 
          ConstantNonzeros, context);
        SUNDANCE_MSG5(verb, tab3 << "CProd = " << CProd);
        SUNDANCE_MSG5(verb, tab3 << "RProd = " << RProd);
        Set<MultipleDeriv> X = CProd.setUnion(RM.setDifference(RProd));
        if (!started) 
        {
          rtn = X;
          started = true;
        }
        else 
        {
          rtn = rtn.intersection(X);
        }
      }
    }

    const Set<MultiSet<int> >& QV = findQ_V(i, context);
    SUNDANCE_MSG5(verb, tab1 << "finding R\\RProd over QV");      
    SUNDANCE_MSG5(verb, tab1 << "QV = " << QV);

    for (Set<MultiSet<int> >::const_iterator j=QV.begin(); j!=QV.end(); j++)
    {
      Tabs tab2;
      Array<int> J = j->elements();
      const Array<Array<int> >& K = comp[J.size()-1];
      SUNDANCE_MSG5(verb, tab2 << "J=" << J);
      SUNDANCE_MSG5(verb, tab2 << "K=" << K);

      for (int k=0; k<K.size(); k++)
      {
        Tabs tab3;
        Set<MultipleDeriv> RProd = product(J, K[k], 
          RequiredNonzeros, context);
        Set<MultipleDeriv> X = RM.setDifference(RProd);
        if (!started) 
        {
          rtn = X;
          started = true;
        }
        else 
        {
          rtn = rtn.intersection(X);
        }
        SUNDANCE_MSG5(verb, tab3 << "RProd = " << RProd);
      }
    }
  }

  rtn = rtn.intersection(RM);
  SUNDANCE_MSG2(verb, tab0 << "C[" << order << "]=" << rtn);
  SUNDANCE_MSG2(verb, tab0 << "done with EWC::internalFindC(" << order
    << ") for " << toString());
  return rtn;
}



Set<MultipleDeriv> 
ExprWithChildren::internalFindW(int order, const EvalContext& context) const
{
  Tabs tab0;
  int verb = context.setupVerbosity();
  Set<MultipleDeriv> rtn;
  SUNDANCE_MSG2(verb, tab0 << "EWC::internalFindW(order="
    << order << ") for " << toString());  
  Tabs tabz;
  /* we'll deal with zero order derivatives specially */
  if (order==0) 
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab1 << "** CASE: order=0");  
    /* If I am an arbitrary nonlinear expression, I cannot be known to
     * be zero regardless of the state of my arguments. Return the
     * zeroth-order derivative. */
    if (!(isLinear() || isProduct())) 
    {
      Tabs tab2;
      SUNDANCE_MSG3(verb, tab2 << "I am neither product nor linear");  
      rtn.put(MultipleDeriv());
      SUNDANCE_MSG3(verb, tab2 << "W[" << order << "]=" << rtn);
      SUNDANCE_MSG3(verb, tab2 << "done with EWC::internalFindW for "
        << toString());
      return rtn;
    }

    /* At this point, I've dealt with arbitrary nonlinear exprs so 
     * I know I'm either a product or a linear expr */

    SUNDANCE_MSG3(verb, tab1 << "getting Q");  
    const Set<MultiSet<int> >& Q = findQ_W(0, context);

    /* If there are no nonzero terms, a linear combination or a product
     * will be zero. Return the empty set. */
    if (Q.size()==0)
    {
      Tabs tab2;
      SUNDANCE_MSG3(verb, tab2 << "W[" << order << "]=" << rtn);
      SUNDANCE_MSG3(verb, tab2 << "done with EWC::internalFindW for "
        << toString());
      return rtn;
    }
      
    /* if I'm a linear combination and any term is nonzero, I am nonzero */
    if (isLinear())
    {
      Tabs tab2;
      SUNDANCE_MSG3(verb, tab2 << "I am a linear combination");  
      rtn.put(MultipleDeriv());      
      SUNDANCE_MSG3(verb, tab2 << "W[" << order << "]=" << rtn);
      SUNDANCE_MSG3(verb, tab2 << "done with EWC::internalFindW for "
        << toString());    
      return rtn;
    }

    /* The only possibility remaining is that I'm a product.
     * If any term is zero, I am zero. If a term is 
     * known to be zero, it's index will not appear in Q, so that 
     * comparing the size of Q and the number of children tells me
     * if I'm zero.
     */
    if (((int) Q.size()) == numChildren())
    {
      rtn.put(MultipleDeriv());
    }
    SUNDANCE_MSG2(verb, tab1 << "I am a product");  
    SUNDANCE_MSG2(verb, tab1 << "W[" << order << "]=" << rtn);
    SUNDANCE_MSG2(verb, tab0 << "done with EWC::internalFindW for "
      << toString());
    return rtn;
  }


  /* now do arbitrary order derivatives with the multivariable chain rule*/
  Array<Array<Array<int> > > comp = compositions(order);
  for (int i=1; i<=order; i++) 
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab1 << "** CASE order=" << i);  

    const Set<MultiSet<int> >& QW = findQ_W(i, context);

    SUNDANCE_MSG3(verb, tab1 << "QW=" << QW);  
      
    for (Set<MultiSet<int> >::const_iterator j=QW.begin(); j!=QW.end(); j++)
    {
      Array<int> J = j->elements();
      const Array<Array<int> >& K = comp[J.size()-1];

      for (int k=0; k<K.size(); k++)
      {
        Set<MultipleDeriv> WProd = product(J, K[k], AllNonzeros, context);
        rtn.merge(WProd);
      }
    }
  }

  SUNDANCE_MSG2(verb, tabz << "W[" << order << "]=" << rtn);
  SUNDANCE_MSG2(verb, tab0 << "done with EWC::internalFindW for "
    << toString());
  return rtn;
}

const Set<MultiSet<int> >& 
ExprWithChildren::findQ_W(int order, 
  const EvalContext& context) const
{
  Tabs tab1;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tab1 << "finding Q_W");  
  if (!contextToQWMap_[order].containsKey(context))
  {
    contextToQWMap_[order].put(context, internalFindQ_W(order, context));
  }
  else
  {
    SUNDANCE_MSG5(verb, tab1 << "using previously computed Q_W");  
  }
  return contextToQWMap_[order].get(context);

}

const Set<MultiSet<int> >& 
ExprWithChildren::findQ_V(int order, 
  const EvalContext& context) const
{
  if (!contextToQVMap_[order].containsKey(context))
  {
    contextToQVMap_[order].put(context, internalFindQ_V(order, context));
  }
  return contextToQVMap_[order].get(context);
}

const Set<MultiSet<int> >& 
ExprWithChildren::findQ_C(int order, 
  const EvalContext& context) const
{
  if (!contextToQCMap_[order].containsKey(context))
  {
    contextToQCMap_[order].put(context, internalFindQ_C(order, context));
  }
  return contextToQCMap_[order].get(context);
}


const Set<MultiSet<int> >& ExprWithChildren::getI_N() const
{
  if (!cachedI_N().containsKey(numChildren()))
  {
    Set<MultiSet<int> > x;
    for (int i=0; i<numChildren(); i++)
    {
      x.put(makeMultiSet<int>(i));
    }
    cachedI_N().put(numChildren(), x);
  }
  return cachedI_N().get(numChildren());
}

Set<MultiSet<int> > ExprWithChildren::indexSetProduct(const Set<MultiSet<int> >& a,
  const Set<MultiSet<int> >& b) const
{
  Set<MultiSet<int> > rtn;
  for (Set<MultiSet<int> >::const_iterator i=a.begin(); i!=a.end(); i++)
  {
    for (Set<MultiSet<int> >::const_iterator j=b.begin(); j!=b.end(); j++)
    {
      MultiSet<int> ab = (*i).merge(*j);
      rtn.put(ab);
    }
  }
  return rtn;
}


Set<MultiSet<int> > ExprWithChildren::internalFindQ_V(int order, 
  const EvalContext& context) const
{
  Tabs tab0;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tab0 << "EWC::internalFindQ_V(order=" << order <<")");
  Set<MultiSet<int> > rtn;

  if (!isLinear())
  {
    bool isVar = false;
    for (int i=0; i<numChildren(); i++)
    {
      if (childIsRequired(i,order,context) && evaluatableChild(i)->findV(0, context).size() != 0) 
      {
        isVar=true;
        break;
      }
    }
    if (isVar) rtn = findQ_W(order, context); 
  }

  SUNDANCE_MSG2(verb, tab0 << "Q_V = " << rtn);
  return rtn;
}

Set<MultiSet<int> > ExprWithChildren::internalFindQ_C(int order, 
  const EvalContext& context) const
{
  if (isLinear()) return findQ_W(order,context);

  return findQ_W(order,context).setDifference(findQ_V(order, context));
}


Set<MultiSet<int> > ExprWithChildren
::internalFindQ_W(int order, 
  const EvalContext& context) const
{
  Tabs tab0;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tab0 << "in internalFindQ_W");
  Set<MultiSet<int> > rtn;
  const Set<MultiSet<int> >& I_N = getI_N();
  SUNDANCE_MSG2(verb, tab0 << "I_N=" << I_N);

  if (isLinear())
  {
    Tabs tab1;
    /* first derivatives of the sum wrt the arguments are 
     * always nonzero */
    if (order==1) 
    {
      SUNDANCE_MSG3(verb, tab1 << "is linear, order=1, so Q_W = I_N");
      return I_N;
    }
    /* zeroth derivatives are nonzero if terms are nonzero */
    if (order==0)
    {
      SUNDANCE_MSG3(verb, tab1 << "is linear, order=0");
      for (int i=0; i<numChildren(); i++)
      {
        const Set<MultipleDeriv>& W_i = evaluatableChild(i)->findW(0,context);
        SUNDANCE_MSG3(verb, tab1 << "is linear, order=0");
        if (W_i.size() > 0) 
        {
          Tabs tab2;
          SUNDANCE_MSG3(verb, tab2 << "child " << i << " is nonzero");
          rtn.put(makeMultiSet(i));
        }
        else
        {
          Tabs tab2;
          SUNDANCE_MSG3(verb, tab2 << "child " << i << " is zero");
        }
      }
    }
  }
  else
  {
    Tabs tab1;
    SUNDANCE_MSG3(verb, tab1 << "is nonlinear, using all index set products");
    rtn = I_N;
      
    for (int i=1; i<order; i++)
    {
      rtn = indexSetProduct(rtn, I_N);
    }
  }
  SUNDANCE_MSG2(verb, tab0 << "rtn=" << rtn);
  return rtn;
}

bool ExprWithChildren::childIsRequired(int index, int diffOrder,
  const EvalContext& context) const
{
  const Set<MultiSet<int> >& Q = findQ_W(diffOrder, context);
  for (Set<MultiSet<int> >::const_iterator it=Q.begin(); it != Q.end(); it++)
  {
    if (it->contains(index)) return true;
  }
  return true;
}

RCP<Array<Set<MultipleDeriv> > > ExprWithChildren
::internalDetermineR(const EvalContext& context,
  const Array<Set<MultipleDeriv> >& RInput) const
{
  Tabs tab0;
  int verb = context.setupVerbosity();
  RCP<Array<Set<MultipleDeriv> > > rtn 
    = rcp(new Array<Set<MultipleDeriv> >(RInput.size()));

  SUNDANCE_MSG2(verb, tab0 << "in internalDetermineR() for " << toString());
  SUNDANCE_MSG2(verb, tab0 << "RInput = " << RInput);



  for (int i=0; i<RInput.size(); i++)
  {
    if (RInput[i].size()==0) continue;
    const Set<MultipleDeriv>& Wi = findW(i, context);
    (*rtn)[i] = RInput[i].intersection(Wi);
  }

  int maxOrder = rtn->size()-1;

  const Set<MultiSet<int> >& Q1 = findQ_W(1, context);
  const Set<MultiSet<int> >& Q2 = findQ_W(2, context);
  const Set<MultiSet<int> >& Q3 = findQ_W(3, context);

  SUNDANCE_MSG3(verb, tab0 << "Q1 = " << Q1);
  SUNDANCE_MSG3(verb, tab0 << "Q2 = " << Q2);
  SUNDANCE_MSG3(verb, tab0 << "Q3 = " << Q3);

  for (int i=0; i<numChildren(); i++)
  {
    Tabs tab1;
    MultiSet<int> mi = makeMultiSet(i);
    SUNDANCE_MSG4(verb, tab1 << "i=" << i << ", Q1_i = " << mi );
    TEUCHOS_TEST_FOR_EXCEPTION(mi.size() != 1, std::logic_error, "unexpected multiset size");
    //      int i = *(mi.begin());
    Set<MultipleDeriv> R11;
    Set<MultipleDeriv> R12;
    Set<MultipleDeriv> R13;
    Set<MultipleDeriv> R22;
    Set<MultipleDeriv> R23;
    Set<MultipleDeriv> R33;
    if (maxOrder >= 1) 
    {
      R11 = (*rtn)[1];
      if (maxOrder >=2) R22 = (*rtn)[2];
      if (maxOrder >=3) R33 = (*rtn)[3];
    }
    if (maxOrder >= 2)
    {
      Tabs tab2;
      Set<MultiSet<int> > jSet = setDivision(Q2, i);
      SUNDANCE_MSG5(verb, tab2 << "Q2/i = " << jSet);
      for (Set<MultiSet<int> >::const_iterator 
             j=jSet.begin(); j!=jSet.end(); j++)
      {
        Tabs tab3;
        TEUCHOS_TEST_FOR_EXCEPTION(j->size()!=1, std::logic_error, 
          "unexpected set size");
        int jIndex = *(j->begin());
        SUNDANCE_MSG5(verb,  tab3 << "j=" << jIndex );
        const Set<MultipleDeriv>& W1j = evaluatableChild(jIndex)->findW(1, context);
        Set<MultipleDeriv> ROverW = setDivision((*rtn)[2], W1j);
        R12.merge(ROverW);
        SUNDANCE_MSG5(verb,  tab3 << "R2=" << (*rtn)[2] );
        SUNDANCE_MSG5(verb,  tab3 << "W1(j)=" << W1j );
        SUNDANCE_MSG5(verb,  tab3 << "R2/W1(j)=" << ROverW );

        if (maxOrder>=3)
        {
          const Set<MultipleDeriv>& W2j 
            = evaluatableChild(jIndex)->findW(2, context);             
          R13.merge(setDivision((*rtn)[3], W2j));
          R23.merge(setDivision((*rtn)[3], W1j));
        }
      }
    }
    if (maxOrder >= 3)
    {
      Set<MultiSet<int> > jkSet = setDivision(Q3, i);
      for (Set<MultiSet<int> >::const_iterator 
             jk=jkSet.begin(); jk!=jkSet.end(); jk++)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(jk->size()!=2, std::logic_error, 
          "unexpected set size");
        Array<int> jka = jk->elements();
        int j = jka[0];
        int k = jka[1];
        const Set<MultipleDeriv>& W1j = evaluatableChild(j)->findW(1, context);
        const Set<MultipleDeriv>& W1k = evaluatableChild(k)->findW(1, context);
        R13.merge(setDivision((*rtn)[3], setProduct(W1j, W1k)));
      }
    }
    SUNDANCE_MSG5(verb,  tab1 << "R11 = " << R11 );
    SUNDANCE_MSG5(verb,  tab1 << "R12 = " << R12 );
    SUNDANCE_MSG5(verb,  tab1 << "R13 = " << R13 );
    SUNDANCE_MSG5(verb,  tab1 << "R22 = " << R22 );
    SUNDANCE_MSG5(verb,  tab1 << "R23 = " << R23 );
    SUNDANCE_MSG5(verb,  tab1 << "R33 = " << R33 );
    Set<MultipleDeriv> R1 = R11;
    R1.merge(R12);
    R1.merge(R13);
    Set<MultipleDeriv> R2 = R22;
    R2.merge(R23);
    Set<MultipleDeriv> R3 = R33;

    Set<MultipleDeriv> R0;
    bool childIsNeeded = (*rtn)[0].size() > 0;
    if (!childIsNeeded && R1.size() > 0) childIsNeeded = childIsRequired(i, 2, context);
    if (!childIsNeeded && R2.size() > 0) childIsNeeded = childIsRequired(i, 3, context);
    if (!childIsNeeded && R3.size() > 0) childIsNeeded = childIsRequired(i, 4, context);
    if (childIsNeeded) R0.put(MultipleDeriv());
      
    Array<Set<MultipleDeriv> > RChild;
      
    RChild.append(R0);
    if (maxOrder >= 1) RChild.append(R1);
    if (maxOrder >= 2) RChild.append(R2);
    if (maxOrder >= 3) RChild.append(R3);
    SUNDANCE_MSG4(verb,  tab1 << "RChild = " << RChild );
    evaluatableChild(i)->determineR(context, RChild);
  }

  printR(verb, rtn);
  SUNDANCE_MSG2(verb, tab0 << "done with EWC::internalDetermineR for "
    << toString());
  
  return rtn;
}




void ExprWithChildren::displayNonzeros(std::ostream& os, const EvalContext& context) const 
{
  Tabs tabs0(0);
  os << tabs0 << "Nonzeros of " << toString() << std::endl;
  if (context.setupVerbosity() > 4)
  {
    os << tabs0 << "Diving into children " << std::endl;
  }

  for (int i=0; i<numChildren(); i++)
  {
    Tabs tab1;
    if (context.setupVerbosity() > 4) os << tab1 << "Child " << i << std::endl;
    evaluatableChild(i)->displayNonzeros(os, context);
  }

  if (context.setupVerbosity() > 4)
    os << tabs0 << "Printing nonzeros for parent " << toString() << std::endl;
  const Set<MultipleDeriv>& W = findW(context);
  const Set<MultipleDeriv>& R = findR(context);
  const Set<MultipleDeriv>& C = findC(context);

  
  for (Set<MultipleDeriv>::const_iterator i=W.begin(); i != W.end(); i++)
  {
    Tabs tab1;
    std::string state = "Variable";
    if (C.contains(*i)) state = "Constant";
    if (!R.contains(*i)) state = "Not Required";
    os << tab1 << std::setw(25) << std::left << i->toString() << ": " << state << std::endl;
  }
}


namespace Sundance {
 
Array<Array<std::pair<int, Array<MultipleDeriv> > > >  
chainRuleDerivsOfArgs(int nArgs,
  const MultiSet<int>& bSet,
  const MultipleDeriv& c)
{
  Array<Array<std::pair<int, Array<MultipleDeriv> > > > rtn;

  /* convert to tuple representation of b */
  Array<int> b(nArgs, 0);
  int J = 0;
  for (MultiSet<int>::const_iterator i=bSet.begin(); i!=bSet.end(); i++)
  {
    b[*i]++;
    J++;
  }
      
  /* count orders of each functional deriv in c */
  Sundance::Map<Deriv, int> counts;
  Array<Deriv> d;
  typedef Sundance::Map<Deriv, int>::const_iterator iter;
  for (MultipleDeriv::const_iterator i=c.begin(); i!=c.end(); i++)
  {
    if (!counts.containsKey(*i)) counts.put(*i, 1);
    else counts[*i]++;
    d.append(*i);
  }

  Array<Array<Array<Array<int> > > > a;
  Array<int> s;
  for (iter i=counts.begin(); i!=counts.end(); i++)
  {
    Array<Array<int> > tmp = nonNegCompositions(i->second, J);
    Array<Array<Array<int> > > ai = bStructure(b, tmp);
    a.append(ai);
    s.append(ai.size());
  }

  Array<Array<int> > ic = indexCombinations(s);



  Array<Array<Array<Array<int> > > > all;
  for (int i=0; i<ic.size(); i++)
  {
    bool good = true;
    int numFuncs = ic[i].size();
    Array<Array<Array<int> > > tmp;

    Array<Array<int> > aTot(nArgs);
    for (int j=0; j<nArgs; j++)
    {
      if (b[j] > 0) aTot[j].resize(b[j]);
      for (int k=0; k<b[j]; k++) aTot[j][k] = 0;
    }
    for (int f=0; f<numFuncs; f++)
    {
      bool skip = false;
      const Array<Array<int> >& e = a[f][ic[i][f]];
      for (int j=0; j<nArgs; j++)
      {
        for (int k=0; k<b[j]; k++)
        {
          aTot[j][k] += e[j][k];
        }
      }
      if (!skip) tmp.append(e);
    }

    for (int j=0; j<nArgs; j++)
    {
      for (int k=0; k<b[j]; k++)
      {
        if (aTot[j][k] == 0) good = false;
      }
    }
    if (good) all.append(tmp);
  }

      
  for (int p=0; p<all.size(); p++)
  {
    Array<std::pair<int, Array<MultipleDeriv> > > terms;
    for (int j=0; j<nArgs; j++)
    {
      std::pair<int, Array<MultipleDeriv> > factors;
      factors.first = j;
      for (int k=0; k<b[j]; k++)
      {
        MultipleDeriv md;
        for (int i=0; i<all[p].size(); i++)
        {
          int order = all[p][i][j][k];
          if (order > 0)
          {
            md.put(d[i]);
          }
        }
        if (md.order() > 0) factors.second.append(md);
      }
      if (factors.second.size() > 0) terms.append(factors);
    }
    rtn.append(terms);
  }

  return rtn;
}


Array<Array<Array<int> > > bStructure(const Array<int>& b,
  const Array<Array<int> >& tmp)
{
  Array<Array<Array<int> > > rtn(tmp.size());
      
  for (int p=0; p<tmp.size(); p++)
  {
    int count=0;
	  rtn[p].resize(b.size());
    for (int j=0; j<b.size(); j++)
    {
      rtn[p][j].resize(b[j]);
      for (int k=0; k<b[j]; k++, count++)
      {
        rtn[p][j][k] = tmp[p][count];
      }
    }
  }
  return rtn;
}
    
Array<OrderedPair<Array<MultiSet<int> >, Array<MultipleDeriv> > >
chainRuleTerms(int s, 
  const MultiSet<int>& lambda,
  const MultipleDeriv& nu)
{
  Array<Array<MultiSet<int> > > allK = multisetCompositions(s, lambda);

  Array<MultipleDeriv> allL = multisetSubsets(nu).elements();

  Array<Array<int> > indexTuples = distinctIndexTuples(s, allL.size());

  Array<Array<MultipleDeriv> > allLTuples;
  for (int i=0; i<indexTuples.size(); i++)
  {
    Array<MultipleDeriv> t(s);
    for (int p=0; p<s; p++) t[p] = allL[indexTuples[i][p]];
    allLTuples.append(t);
  }

  Array<OrderedPair<Array<MultiSet<int> >, Array<MultipleDeriv> > > rtn;
  for (int i=0; i<allLTuples.size(); i++)
  {
    for (int j=0; j<allK.size(); j++)
    {
      MultipleDeriv result;
      for (int p=0; p<s; p++)
      {
        for (unsigned int q=0; q<allK[j][p].size(); q++)
        {
          result = result.product(allLTuples[i][p]);
        }
      }
      if (result==nu) 
      {
        OrderedPair<Array<MultiSet<int> >, Array<MultipleDeriv> >
          kl(allK[j], allLTuples[i]);
        rtn.append(kl);
      }
    }
  }
  return rtn;
}

Set<MultipleDeriv> multisetSubsets(const MultipleDeriv& nu)
{
  /* We'll generate the subsets by traversing them in bitwise order.
   * For a multiset having N elements, there are up to 2^N subsets each
   * of which can be described by a N-bit number with the i-th
   * bit indicating whether the i-th element is in the subset. Note that
   * with a multiset, repetitions can occur so we need to record the
   * results in a Set object to eliminate duplicates. 
   */

  /* Make an indexable array of the elements. This will be convenient
   * because we'll need to access the i-th element after reading
   * the i-th bit. */
  Array<Deriv> elements = nu.elements();

  /* Compute the maximum number of subsets. This number will be reached
   * only in the case of no repetitions. */
  int n = elements.size();
  int maxNumSubsets = pow2(n);
    
  Set<MultipleDeriv> rtn;
    
    
  /* Loop over subsets in bitwise order. We start the count at 1 
     to avoid including the empty subset */
  for (int i=1; i<maxNumSubsets; i++)
  {
    Array<int> bits = bitsOfAnInteger(i, n);
    MultipleDeriv md;
    for (int j=0; j<n; j++) 
    {
      if (bits[j] == 1) md.put(elements[j]);
    }
    rtn.put(md);
  }
  return rtn;
}


  
int chainRuleMultiplicity(const MultipleDeriv& nu,
  const Array<MultiSet<int> >& K,
  const Array<MultipleDeriv>& L)
{
  int rtn = factorial(nu);
  for (int i=0; i<K.size(); i++)
  {
    rtn = rtn/factorial(K[i]);
    int lFact = factorial(L[i]);
    int lPow = 1;
    int kNorm = K[i].size();
    for (int j=0; j<kNorm; j++) lPow *= lFact;
    TEUCHOS_TEST_FOR_EXCEPT(rtn % lPow != 0);
    rtn = rtn/lPow;
  }
  return rtn;
}


}



