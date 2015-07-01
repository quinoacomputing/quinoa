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

#include "SundanceEvaluatableExpr.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceNullEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceExpr.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Utils.hpp"

using namespace Sundance;
using namespace Teuchos;
using std::endl;

TEUCHOS_TIMER(evalTimer, "Symbolic Evaluation")

  EvaluatableExpr::EvaluatableExpr()
    : ScalarExpr(), 
      evaluators_(),
      sparsity_(),
      nodesHaveBeenCounted_(false),
      contextToDSSMap_(numDerivSubsetTypes(), contextToDSSMap_ele_t(maxFuncDiffOrder()+1)),
      rIsReady_(false),
      allOrdersMap_(numDerivSubsetTypes())
{}


Time& EvaluatableExpr::evaluationTimer()
{
  return evalTimer();
}

string EvaluatableExpr::derivType(const DerivSubsetSpecifier& dss) const
{
  switch(dss)
  {
    case RequiredNonzeros:
      return "R";
    case ConstantNonzeros:
      return "C";
    case VariableNonzeros:
      return "V";
    default:
      return "W";
  }
}

void EvaluatableExpr::registerSpatialDerivs(const EvalContext& context, 
  const Set<MultiIndex>& miSet) const
{
  Tabs tab;
  if (activeSpatialDerivMap_.containsKey(context))
  {
    Set<MultiIndex> s = activeSpatialDerivMap_[context];
    s.merge(miSet);
    activeSpatialDerivMap_.put(context, s);
  }
  else
  {
    activeSpatialDerivMap_.put(context, miSet);
  }
}

const Set<MultiIndex>& EvaluatableExpr
::activeSpatialDerivs(const EvalContext& context) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!activeSpatialDerivMap_.containsKey(context),
    std::logic_error,
    "Unknown context " << context);
  return activeSpatialDerivMap_[context];
}

RCP<SparsitySuperset> 
EvaluatableExpr::sparsitySuperset(const EvalContext& context) const 
{
  Tabs tab;
  
  SUNDANCE_MSG2(context.setupVerbosity(), 
    tab << "getting sparsity superset for " << toString());

  RCP<SparsitySuperset> rtn;

  if (sparsity_.containsKey(context))
  {
    Tabs tab1;
    SUNDANCE_MSG2(context.setupVerbosity(), 
      tab1 << "reusing previously computed data...");
    rtn = sparsity_.get(context);
  }
  else
  {
    Tabs tab1;
    SUNDANCE_MSG2(context.setupVerbosity(), 
      tab1 << "computing from scratch...");
    const Set<MultipleDeriv>& R = findR(context);
    const Set<MultipleDeriv>& C = findC(context);
    const Set<MultipleDeriv>& V = findV(context);
    if (context.setupVerbosity() > 4)
    {
      Out::os() << tab1 << "R=" << R << endl;
      Out::os() << tab1 << "C=" << C << endl;
      Out::os() << tab1 << "V=" << V << endl;
    }
    rtn = rcp(new SparsitySuperset(C.intersection(R), V.intersection(R)));
    sparsity_.put(context, rtn);
  }
  return rtn;
}



const RCP<Evaluator>&
EvaluatableExpr::evaluator(const EvalContext& context) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!evaluators_.containsKey(context), std::runtime_error, 
    "Evaluator not found for context " << context);
  return evaluators_.get(context);
}

void EvaluatableExpr::evaluate(const EvalManager& mgr,
  Array<double>& constantResults,
  Array<RCP<EvalVector> >& vectorResults) const
{
  TimeMonitor timer(evalTimer());
  evaluator(mgr.getRegion())->eval(mgr, constantResults, vectorResults);
}

void EvaluatableExpr::setupEval(const EvalContext& context) const
{
  Tabs tabs0(0);
  int verb = context.evalSetupVerbosity();
  SUNDANCE_MSG1(verb, tabs0 << "setupEval() for " << this->toString());
  if (!evaluators_.containsKey(context))
  {
    Tabs tabs;
    SUNDANCE_MSG2(verb, tabs << "creating new evaluator...");
    SUNDANCE_MSG2(verb, tabs << "my sparsity superset = " << std::endl
      << *sparsitySuperset(context));
    RCP<Evaluator> eval;
    if (sparsitySuperset(context)->numDerivs()>0)
    {
      SUNDANCE_MSG2(verb, tabs << "calling createEvaluator()");
      eval = rcp(createEvaluator(this, context));
    }
    else
    {
      SUNDANCE_MSG2(verb, 
        tabs << "EE: no results needed... creating null evaluator");
      eval = rcp(new NullEvaluator());
    }
    evaluators_.put(context, eval);
  }
  else
  {
    Tabs tabs;
    SUNDANCE_MSG2(verb, tabs << "reusing existing evaluator...");
  }
}

void EvaluatableExpr::showSparsity(std::ostream& os, 
  const EvalContext& context) const
{
  Tabs tab0;
  os << tab0 << "Node: " << toString() << std::endl;
  sparsitySuperset(context)->print(os);
}



int EvaluatableExpr::maxOrder(const Set<MultiIndex>& m) const 
{
  int rtn = 0;
  for (Set<MultiIndex>::const_iterator i=m.begin(); i != m.end(); i++)
  {
    rtn = std::max(i->order(), rtn);
  }
  return rtn;
}

const EvaluatableExpr* EvaluatableExpr::getEvalExpr(const Expr& expr)
{
  const EvaluatableExpr* rtn 
    = dynamic_cast<const EvaluatableExpr*>(expr[0].ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(rtn==0, std::logic_error,
    "cast of " << expr 
    << " failed in EvaluatableExpr::getEvalExpr()");
  TEUCHOS_TEST_FOR_EXCEPTION(expr.size() != 1, std::logic_error,
    "non-scalar expression " << expr
    << " in EvaluatableExpr::getEvalExpr()");

  return rtn;
}


bool EvaluatableExpr::isEvaluatable(const ExprBase* expr) 
{
  return dynamic_cast<const EvaluatableExpr*>(expr) != 0;
}


int EvaluatableExpr::countNodes() const
{
  nodesHaveBeenCounted_ = true;
  return 1;
}



const Set<MultipleDeriv>& 
EvaluatableExpr::findW(int order,
  const EvalContext& context) const
{
  return findDerivSubset(order, AllNonzeros, context);
}

const Set<MultipleDeriv>& 
EvaluatableExpr::findV(int order,
  const EvalContext& context) const
{
  return findDerivSubset(order, VariableNonzeros, context);
}

const Set<MultipleDeriv>& 
EvaluatableExpr::findC(int order,
  const EvalContext& context) const
{
  return findDerivSubset(order, ConstantNonzeros, context);
}

const Set<MultipleDeriv>& 
EvaluatableExpr::findR(int order,
  const EvalContext& context) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!rIsReady_, std::logic_error,
    "findR() cannot be used for initial computation of the "
    "R subset. Calling object is " << toString());
  return findDerivSubset(order, RequiredNonzeros, context);
}


const Set<MultipleDeriv>& 
EvaluatableExpr::findDerivSubset(int order,
  const DerivSubsetSpecifier& dss,
  const EvalContext& context) const
{
  Tabs tabs;
  int verb = context.setupVerbosity();
  SUNDANCE_MSG2(verb, tabs << "finding " << derivType(dss) 
    << "[" << order << "] for " << toString());
  if (contextToDSSMap_[dss][order].containsKey(context))
  {
    Tabs tab1;
    SUNDANCE_MSG3(verb, tab1 << "...reusing previously computed data");
  }
  else
  {
    Tabs tab1;
    SUNDANCE_MSG3(verb, tab1 << "...need to call find<Set>");
    switch(dss)
    {
      case AllNonzeros:
        contextToDSSMap_[dss][order].put(context, internalFindW(order, context));
        break;
      case RequiredNonzeros:
        break;
      case VariableNonzeros:
        contextToDSSMap_[dss][order].put(context, internalFindV(order, context));
        break;
      case ConstantNonzeros:
        contextToDSSMap_[dss][order].put(context, internalFindC(order, context));
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "this should never happen");
    }
  }


  if (!contextToDSSMap_[dss][order].containsKey(context))
  {
    contextToDSSMap_[dss][order].put(context, Set<MultipleDeriv>());
  }
  const Set<MultipleDeriv>& rtn = contextToDSSMap_[dss][order].get(context);
  SUNDANCE_MSG2(verb, tabs << "found " << rtn);
  return rtn;
}


Set<MultipleDeriv> 
EvaluatableExpr::setProduct(const Set<MultipleDeriv>& a,
  const Set<MultipleDeriv>& b) const
{
  Set<MultipleDeriv> rtn;
  for (Set<MultipleDeriv>::const_iterator i=a.begin(); i!=a.end(); i++)
  {
    for (Set<MultipleDeriv>::const_iterator j=b.begin(); j!=b.end(); j++)
    {
      rtn.put(i->product(*j));
    }
  }
  return rtn;
}

Set<MultiSet<int> > EvaluatableExpr::setDivision(const Set<MultiSet<int> >& s,
  int index) const 
{
  Set<MultiSet<int> > rtn;
  for (Set<MultiSet<int> >::const_iterator 
         i=s.begin(); i!=s.end(); i++)
  {
    if (i->contains(index))
    {
      MultiSet<int> e = *i;
      e.erase(e.find(index));
      rtn.put(e);
    }
  }
  return rtn;
}


Set<MultipleDeriv>  
EvaluatableExpr::setDivision(const Set<MultipleDeriv>& a,
  const Set<MultipleDeriv>& b) const
{
  Set<MultipleDeriv> rtn;
  for (Set<MultipleDeriv>::const_iterator i=a.begin(); i!=a.end(); i++)
  {
    for (Set<MultipleDeriv>::const_iterator j=b.begin(); j!=b.end(); j++)
    {
      MultipleDeriv c = i->factorOutDeriv(*j);
      if (c.size() != 0) rtn.put(c);
      if (*i == *j) rtn.put(MultipleDeriv());
    }
  }
  return rtn;
}

void EvaluatableExpr::determineR(const EvalContext& context,
  const Array<Set<MultipleDeriv> >& RInput) const
{
  Tabs tabs;
  RCP<Array<Set<MultipleDeriv> > > additionToR 
    =  internalDetermineR(context, RInput);
  
  for (int i=0; i<RInput.size(); i++)
  {
    if ((*additionToR)[i].size()==0) continue;
    if (!contextToDSSMap_[RequiredNonzeros][i].containsKey(context))
    {
      contextToDSSMap_[RequiredNonzeros][i].put(context, (*additionToR)[i]); 
    }
    else
    {
      contextToDSSMap_[RequiredNonzeros][i][context].merge((*additionToR)[i]);
    }
  }
  rIsReady_ = true;

}

RCP<Array<Set<MultipleDeriv> > > EvaluatableExpr
::internalDetermineR(const EvalContext& context,
  const Array<Set<MultipleDeriv> >& RInput) const
{
  Tabs tab0(0);
  int verb = context.setupVerbosity();
  RCP<Array<Set<MultipleDeriv> > > rtn 
    = rcp(new Array<Set<MultipleDeriv> >(RInput.size()));

  SUNDANCE_MSG2(verb, tab0 << "EE::internalDetermineR() for " << toString() );
  SUNDANCE_MSG2(verb, tab0 << "RInput = " << RInput );

  for (int i=0; i<RInput.size(); i++)
  {
    Tabs tab1;
    if (RInput[i].size()==0) continue;
    const Set<MultipleDeriv>& Wi = findW(i, context);
    SUNDANCE_MSG3(verb, tab1 << "W[" << i << "] = " << Wi );
    (*rtn)[i] = RInput[i].intersection(Wi);
  }

  printR(verb, rtn);
  return rtn;
}


Array<Set<MultipleDeriv> > EvaluatableExpr
::computeInputR(const EvalContext& context,
  const Array<Set<MultiSet<int> > >& funcIDCombinations,
  const Array<Set<MultiIndex> >& spatialDerivs) const 
{
  int verb = context.setupVerbosity();
  Tabs tabs(0);
  SUNDANCE_MSG2(verb, tabs << "in EE::computeInputR()");
  Array<Set<MultipleDeriv> > rtn(funcIDCombinations.size());
  for (int order=0; order<funcIDCombinations.size(); order++)
  {
    Tabs tabs1;
    SUNDANCE_MSG2(verb, tabs << "order=" << order);
    SUNDANCE_MSG2(verb, tabs1 << "funcID combs=" 
      << funcIDCombinations[order]);
    const Set<MultipleDeriv>& W = findW(order, context);
    SUNDANCE_MSG2(verb, tabs1 << "W[order=" << order << "]=" << W);

    for (Set<MultipleDeriv>::const_iterator i=W.begin(); i!=W.end(); i++)
    {
      Tabs tab2;
      const MultipleDeriv& nonzeroDeriv = *i;
      if (nonzeroDeriv.isInRequiredSet(funcIDCombinations[order],
          spatialDerivs[order]))
      {
        SUNDANCE_MSG2(verb, tab2 << "md=" << nonzeroDeriv 
          << " added to inputR");
        rtn[order].put(nonzeroDeriv);
      }
      else
      {
        SUNDANCE_MSG2(verb, tab2 << "md=" << nonzeroDeriv << " not needed");
      }
    }
  }
  SUNDANCE_MSG2(verb, tabs << "inputR=" << rtn);
  return rtn;
}


const Set<MultipleDeriv>& EvaluatableExpr::findW(const EvalContext& context) const
{
  return findDerivSubset(AllNonzeros, context);
}

const Set<MultipleDeriv>& EvaluatableExpr::findR(const EvalContext& context) const
{
  return findDerivSubset(RequiredNonzeros, context);
}

const Set<MultipleDeriv>& EvaluatableExpr::findC(const EvalContext& context) const
{
  return findDerivSubset(ConstantNonzeros, context);
}

const Set<MultipleDeriv>& EvaluatableExpr::findV(const EvalContext& context) const
{
  return findDerivSubset(VariableNonzeros, context);
}

const Set<MultipleDeriv>& 
EvaluatableExpr::findDerivSubset(const DerivSubsetSpecifier& dss,
  const EvalContext& context) const
{
  if (!allOrdersMap_[dss].containsKey(context))
  {
    Set<MultipleDeriv> tmp;
    for (int i=0; i<=3; i++)
    {
      tmp.merge(findDerivSubset(i, dss, context));
    }
    allOrdersMap_[dss].put(context, tmp);
  }

  return allOrdersMap_[dss].get(context);
}

void EvaluatableExpr::displayNonzeros(std::ostream& os, const EvalContext& context) const 
{
  Tabs tabs0;
  os << tabs0 << "Nonzeros of " << toString() << std::endl;

  const Set<MultipleDeriv>& W = findW(context);
  const Set<MultipleDeriv>& R = findR(context);
  const Set<MultipleDeriv>& C = findC(context);
  const Set<MultipleDeriv>& V = findV(context);

  TEUCHOS_TEST_FOR_EXCEPT(C.intersection(V).size() != 0);

  for (Set<MultipleDeriv>::const_iterator i=W.begin(); i != W.end(); i++)
  {
    Tabs tab1;
    std::string state = "Variable";
    if (C.contains(*i)) state = "Constant";
    if (!R.contains(*i)) state = "Not Required";
    os << tab1 << std::setw(25) << std::left << i->toString() << ": " << state << std::endl;
  }
}


void EvaluatableExpr::printR(int verb,
  const RCP<Array<Set<MultipleDeriv> > >& R) const
{
  Tabs tab(0);
  const Array<Set<MultipleDeriv> >& r = *R;
  for (int order=0; order<r.size(); order++)
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb, tab << "order=" << order);
    SUNDANCE_MSG2(verb, tab1 << "R[" << order << "]=" << r[order]);
  }
}
