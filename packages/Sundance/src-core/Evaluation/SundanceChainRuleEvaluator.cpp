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

#include "SundanceChainRuleEvaluator.hpp"
#include "SundanceCombinatorialUtils.hpp"

#include "SundanceUnknownFuncElement.hpp"
#include "SundanceEvalManager.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceSet.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"


using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


ChainRuleEvaluator::ChainRuleEvaluator(const ExprWithChildren* expr, 
  const EvalContext& context)
  : SubtypeEvaluator<ExprWithChildren>(expr, context), 
    expansions_(),
    childEvaluators_(expr->numChildren()),
    childSparsity_(expr->numChildren()),
    constArgDerivMap_(),
    varArgDerivMap_(),
    zerothDerivResultIndex_(-1),
    zerothDerivIsConstant_(false)
{
  Tabs tabs;
  SUNDANCE_MSG1(context.setupVerbosity(),
    tabs << "ChainRuleEvaluator base class ctor for " 
    << expr->toString());
  for (int i=0; i<numChildren(); i++)
  {
    childEvaluators_[i] = expr->evaluatableChild(i)->evaluator(context);
    childEvaluators_[i]->addClient();
    childSparsity_[i] = expr->evaluatableChild(i)->sparsitySuperset(context);
  }
}

Sundance::Map<OrderedPair<int, int>, Array<Array<int> > >& ChainRuleEvaluator::compMap()
{
  static Map<OrderedPair<int, int>, Array<Array<int> > > rtn;
  return rtn;
}

void ChainRuleEvaluator::resetNumCalls() const
{
  for (int i=0; i<numChildren(); i++)
  {
    childEvaluators_[i]->resetNumCalls();
  }
  Evaluator::resetNumCalls();
}


void ChainRuleEvaluator::addConstArgDeriv(const MultiSet<int>& df, int index)
{
  constArgDerivMap_.put(df, index);
}

void ChainRuleEvaluator::addVarArgDeriv(const MultiSet<int>& df, int index)
{
  varArgDerivMap_.put(df, index);
}

int ChainRuleEvaluator::constArgDerivIndex(const MultiSet<int>& df) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!constArgDerivMap_.containsKey(df), std::logic_error,
    "argument derivative " << df << " not found in constant "
    "argument derivative map");

  return constArgDerivMap_.get(df);
}

int ChainRuleEvaluator::varArgDerivIndex(const MultiSet<int>& df) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!varArgDerivMap_.containsKey(df), std::logic_error,
    "argument derivative " << df << " not found in variable "
    "argument derivative map");

  return varArgDerivMap_.get(df);
}


const Array<Array<int> >& ChainRuleEvaluator::nComps(int N, int n) const
{
  OrderedPair<int,int> key(n,N);
  if (!compMap().containsKey(key))
  {
    compMap().put(key, compositions(N)[n-1]);
  }
  return compMap().get(key);
}


double ChainRuleEvaluator::fact(int n) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(n<0, std::logic_error, "negative argument " << n << " to factorial");
  if (n==0 || n==1) return 1.0;
  return n*fact(n-1);
}

double ChainRuleEvaluator::choose(int N, int n) const
{
  return fact(N)/fact(n)/fact(N-n);
}

double ChainRuleEvaluator::stirling2(int n, int k) const
{
  if (n < k) return 0;
  if (n == k) return 1;
  if (k<=0) return 0;
  if (k==1) return 1;
  if (n-1 == k) return choose(n, 2);
  return k*stirling2(n-1, k) + stirling2(n-1, k-1);
}


MultipleDeriv ChainRuleEvaluator::makeMD(const Array<Deriv>& d) 
{
  MultipleDeriv rtn;
  for (int i=0; i<d.size(); i++)
  {
    rtn.put(d[i]);
  }
  return rtn;
}


Set<MultiSet<MultipleDeriv> > ChainRuleEvaluator::chainRuleBins(const MultipleDeriv& d,
  const MultiSet<int>& q)
{
  int n = q.size();
  Array<Array<Array<Deriv> > > bins = binnings(d, n);

  Set<MultiSet<MultipleDeriv> > rtn;

  for (int i=0; i<bins.size(); i++)
  {
    MultiSet<MultipleDeriv> b;
    for (int j=0; j<bins[i].size(); j++)
    {
      b.put(makeMD(bins[i][j]));
    }
    rtn.put(b);
  }


  return rtn;
}


int ChainRuleEvaluator::derivComboMultiplicity(const MultiSet<MultipleDeriv>& b) const
{
  /* ugly brute force code -- there has got to be a better way to compute multiplicities */

  MultipleDeriv dTot;
  Array<MultiSet<Deriv> > derivSets(b.size());
  Array<Array<Deriv> > derivArrays(b.size());
  Set<Deriv> allDerivs;
  int k=0;
  bool allDerivsAreDistinct = true;
  bool allDerivsAreIdentical = true;
  for (MultiSet<MultipleDeriv>::const_iterator i=b.begin(); i!=b.end(); i++, k++)
  {
    for (MultipleDeriv::const_iterator j=i->begin(); j!=i->end(); j++)
    {
      derivSets[k].put(*j);
      derivArrays[k].append(*j);
      dTot.put(*j);
      if (allDerivs.contains(*j)) allDerivsAreDistinct = false;
      if (allDerivs.size()>0 && !allDerivs.contains(*j)) allDerivsAreIdentical = false;
      allDerivs.put(*j);
    }
  }
  int totOrder = dTot.order();

  /* eliminate 4th order or higher */
  TEUCHOS_TEST_FOR_EXCEPTION(totOrder > 3, std::logic_error,
    "deriv order " << totOrder << " not supported");

  if (b.size()==1) return 1;  /* handles case with a single multiple deriv */
  if (totOrder == (int) b.size()) return 1; /* handles case with N first derivs */

  /* The only remaining cases are total order = 3, grouped in 2 bins. 
   * Order=3 in 1 bin and order=3 in 3 bins have been delat with already. 
   * The ways of grouping 3 derivatives in 2 bins are:
   * {(a,bc), (b,ac), (c, ab)}.
   * If all three first derivs are the same, we have (a,aa) with multiplicity 3.
   * If all are distinct, we have (a,bc) with multiplicity 1. 
   * If two coincide, we have either
   * (a,ba) with multiplicity 2, or
   * (b,aa) with multiplicity 1.
   */
  TEUCHOS_TEST_FOR_EXCEPTION(derivArrays.size() != 2, std::logic_error,
    "unexpected size=" << derivArrays.size());

  if (allDerivsAreIdentical) return 3;
  if (allDerivsAreDistinct) return 1;

  if (derivArrays[0].size()==1) 
  {
    if (derivSets[1].contains(derivArrays[0][0])) return 2;
    return 1;
  }
  else
  {
    if (derivSets[0].contains(derivArrays[1][0])) return 2;
    return 1;
  }
}


void ChainRuleEvaluator::init(const ExprWithChildren* expr, 
  const EvalContext& context)
{
  int verb = context.setupVerbosity();

  typedef Array<OrderedPair<Array<MultiSet<int> >, Array<MultipleDeriv> > > CR;
  Tabs tabs;
  SUNDANCE_MSG1(verb, tabs << "ChainRuleEvaluator::init() for " 
    << expr->toString());

  const Set<MultipleDeriv>& C = expr->findC(context);
  const Set<MultipleDeriv>& R = expr->findR(context);

  Array<Set<MultipleDeriv> > argV(expr->numChildren());
  Array<Set<MultipleDeriv> > argC(expr->numChildren());
  Array<Set<MultipleDeriv> > argR(expr->numChildren());

  for (int i=0; i<numChildren(); i++)
  {
    argV[i] = expr->evaluatableChild(i)->findV(context);
    argC[i] = expr->evaluatableChild(i)->findC(context);
    argR[i] = expr->evaluatableChild(i)->findR(context);
  }
  SUNDANCE_MSG3(verb, tabs << "sparsity = " << *(this->sparsity()));
  typedef Set<MultipleDeriv>::const_iterator iter;

  int count=0;
  int vecResultIndex = 0;
  int constResultIndex = 0;
  for (iter md=R.begin(); md!=R.end(); md++, count++)
  {
    Tabs tab1;
    SUNDANCE_MSG3(verb, tab1 << "working out evaluator for " << *md);
    int N = md->order();
    bool resultIsConstant = C.contains(*md);
    int resultIndex;
    if (resultIsConstant)
    {
      Tabs tab2;
      SUNDANCE_MSG3(verb, tab2 << "result is constant, const index=" << constResultIndex);
      addConstantIndex(count, constResultIndex);
      resultIndex = constResultIndex;
      constResultIndex++;
    }
    else
    {
      Tabs tab2;
      SUNDANCE_MSG3(verb, tab2 << "result is variable, vec index=" << vecResultIndex);
      addVectorIndex(count, vecResultIndex);
      resultIndex = vecResultIndex;
      vecResultIndex++;
    }

    SUNDANCE_MSG3(verb, tab1 << "order=" << N);
      
    if (N==0)
    {
      Tabs tab2;
      SUNDANCE_MSG3(verb, tab2 << "zeroth deriv index=" << resultIndex);
      zerothDerivIsConstant_ = resultIsConstant;
      zerothDerivResultIndex_ = resultIndex;
      continue;
    }


      
    RCP<ChainRuleSum> sum 
      = rcp(new ChainRuleSum(*md, resultIndex, resultIsConstant));

    const MultipleDeriv& nu = *md;

    for (int n=1; n<=N; n++)
    {
      Tabs tab2;
      SUNDANCE_MSG3(verb, tab2 << "n=" << n);
      const Set<MultiSet<int> >& QW = expr->findQ_W(n, context);
      const Set<MultiSet<int> >& QC = expr->findQ_C(n, context);
      SUNDANCE_MSG3(verb, tab2 << "Q_W=" << QW);
      SUNDANCE_MSG3(verb, tab2 << "Q_C=" << QC);
      for (Set<MultiSet<int> >::const_iterator 
             j=QW.begin(); j!=QW.end(); j++)
      {
        Tabs tab3;
        const MultiSet<int>& lambda = *j;
        SUNDANCE_MSG3(verb, tab3 << "arg index set =" << lambda);
        bool argDerivIsConstant = QC.contains(lambda);
        int argDerivIndex = -1;
        if (argDerivIsConstant) 
        {
          argDerivIndex = constArgDerivIndex(lambda);
        }
        else 
        {
          argDerivIndex = varArgDerivIndex(lambda);
        }
        Array<DerivProduct> pSum;
        for (int s=1; s<=N; s++)
        {
          Tabs tab4;
          SUNDANCE_MSG3(verb, tab4 << "preparing chain rule terms for "
            "s=" << s << ", lambda=" << lambda << ", nu=" << nu);
          CR p = chainRuleTerms(s, lambda, nu);
          for (CR::const_iterator j=p.begin(); j!=p.end(); j++)
          {
            Tabs tab5;
            Array<MultiSet<int> > K = j->first();
            Array<MultipleDeriv> L = j->second();
            SUNDANCE_MSG3(verb, tab5 << "K=" << K << std::endl << tab5 << "L=" << L);
            double weight = chainRuleMultiplicity(nu, K, L);
            SUNDANCE_MSG3(verb, tab5 << "weight=" << weight);
            DerivProduct prod(weight);
            bool termIsZero = false;
            for (int j=0; j<K.size(); j++)
            {
              for (MultiSet<int>::const_iterator 
                     k=K[j].begin(); k!=K[j].end(); k++)
              {
                int argIndex = *k;
                const MultipleDeriv& derivOfArg = L[j];
                const RCP<SparsitySuperset>& argSp 
                  = childSparsity_[argIndex];
                const RCP<Evaluator>& argEv
                  = childEvaluators_[argIndex];
                               
                int rawValIndex = argSp->getIndex(derivOfArg);
                SUNDANCE_MSG3(verb, tab5 << "argR=" 
                  << argR[argIndex]);
                if (argV[argIndex].contains(derivOfArg))
                {
                  SUNDANCE_MSG3(verb, tab5 << "mdArg is variable");
                  int valIndex 
                    = argEv->vectorIndexMap().get(rawValIndex);
                  prod.addVariableFactor(IndexPair(argIndex, valIndex));
                }
                else if (argC[argIndex].contains(derivOfArg))
                {
                  SUNDANCE_MSG3(verb, tab5 << "mdArg is constant");
                  int valIndex 
                    = argEv->constantIndexMap().get(rawValIndex);
                  prod.addConstantFactor(IndexPair(argIndex, valIndex));
                }
                else
                {
                  SUNDANCE_MSG3(verb, tab5 << "mdArg is zero");
                  termIsZero = true;
                  break;
                }
              }
              if (termIsZero) break;
            }
            if (!termIsZero) pSum.append(prod);
          }
        }
        sum->addTerm(argDerivIndex, argDerivIsConstant, pSum);
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(sum->numTerms()==0, std::logic_error,
      "Empty sum in chain rule expansion for derivative "
      << *md);
    expansions_.append(sum);
  }

  SUNDANCE_MSG3(verb, tabs << "num constant results: " 
    << this->sparsity()->numConstantDerivs());

  SUNDANCE_MSG3(verb, tabs << "num var results: " 
    << this->sparsity()->numVectorDerivs());

  
}



void ChainRuleEvaluator::internalEval(const EvalManager& mgr,
  Array<double>& constantResults,
  Array<RCP<EvalVector> >& vectorResults) const 
{
  TimeMonitor timer(chainRuleEvalTimer());
  Tabs tabs(0);

  SUNDANCE_MSG1(mgr.verb(), tabs << "ChainRuleEvaluator::eval() expr=" 
    << expr()->toString());

  
  SUNDANCE_MSG2(mgr.verb(), tabs << "max diff order = " << mgr.getRegion().topLevelDiffOrder());
  SUNDANCE_MSG2(mgr.verb(), tabs << "return sparsity " << std::endl << tabs << *(this->sparsity()));
  
  constantResults.resize(this->sparsity()->numConstantDerivs());
  vectorResults.resize(this->sparsity()->numVectorDerivs());

  SUNDANCE_MSG3(mgr.verb(),tabs << "num constant results: " 
    << this->sparsity()->numConstantDerivs());

  SUNDANCE_MSG3(mgr.verb(),tabs << "num var results: " 
    << this->sparsity()->numVectorDerivs());

  Array<RCP<Array<double> > > constantArgResults(numChildren());
  Array<RCP<Array<RCP<EvalVector> > > > varArgResults(numChildren());

  Array<double> constantArgDerivs;
  Array<RCP<EvalVector> > varArgDerivs;

  for (int i=0; i<numChildren(); i++)
  {
    Tabs tab1;
    SUNDANCE_MSG3(mgr.verb(), tab1 << "computing results for child #"
      << i);
                         
    constantArgResults[i] = rcp(new Array<double>());
    varArgResults[i] = rcp(new Array<RCP<EvalVector> >());
    childEvaluators_[i]->eval(mgr, *(constantArgResults[i]), 
      *(varArgResults[i]));
    if (mgr.verb() > 3)
    {
      Out::os() << tabs << "constant arg #" << i << 
        " results:" << *(constantArgResults[i]) << std::endl;
      Out::os() << tabs << "variable arg #" << i << " derivs:" << std::endl;
      for (int j=0; j<varArgResults[i]->size(); j++)
      {
        Tabs tab1;
        Out::os() << tab1 << j << " ";
        (*(varArgResults[i]))[j]->print(Out::os());
        Out::os() << std::endl;
      }
    }
  }

  evalArgDerivs(mgr, constantArgResults, varArgResults,
    constantArgDerivs, varArgDerivs);

  
  if (mgr.verb() > 2)
  {
    Out::os() << tabs << "constant arg derivs:" << constantArgDerivs << std::endl;
    Out::os() << tabs << "variable arg derivs:" << std::endl;
    for (int i=0; i<varArgDerivs.size(); i++)
    {
      Tabs tab1;
      Out::os() << tab1 << i << " ";
      varArgDerivs[i]->print(Out::os());
      Out::os() << std::endl;
    }
  }
  

  for (int i=0; i<expansions_.size(); i++)
  {
    Tabs tab1;
    int resultIndex = expansions_[i]->resultIndex();
    bool isConstant = expansions_[i]->resultIsConstant();
    SUNDANCE_MSG3(mgr.verb(), tab1 << "doing expansion for deriv #" << i
      << ", result index="
      << resultIndex << std::endl << tab1
      << "deriv=" << expansions_[i]->deriv());
    if (isConstant)
    {
      expansions_[i]->evalConstant(mgr, constantArgResults, constantArgDerivs,
        constantResults[resultIndex]);
    }
    else
    {
      expansions_[i]->evalVar(mgr, constantArgResults, varArgResults,
        constantArgDerivs, varArgDerivs,
        vectorResults[resultIndex]);
    }
  }

  if (zerothDerivResultIndex_ >= 0)
  {
    SUNDANCE_MSG3(mgr.verb(), tabs << "processing zeroth-order deriv");
    Tabs tab1;
    SUNDANCE_MSG3(mgr.verb(), tab1 << "result index = " << zerothDerivResultIndex_);
    if (zerothDerivIsConstant_)
    {
      Tabs tab2;
      SUNDANCE_MSG3(mgr.verb(), tab2 << "zeroth-order deriv is constant");
      constantResults[zerothDerivResultIndex_] = constantArgDerivs[0];
    }
    else
    {
      Tabs tab2;
      SUNDANCE_MSG3(mgr.verb(), tab2 << "zeroth-order deriv is variable");
      vectorResults[zerothDerivResultIndex_] = varArgDerivs[0];
    }
  }


  if (mgr.verb() > 1)
  {
    Tabs tab1;
    Out::os() << tab1 << "chain rule results " << std::endl;
    mgr.showResults(Out::os(), this->sparsity(), vectorResults,
      constantResults);
  }

  SUNDANCE_MSG1(mgr.verb(), tabs << "ChainRuleEvaluator::eval() done"); 
}




namespace Sundance {

MultipleDeriv makeDeriv(const Expr& a)
{
  const UnknownFuncElement* aPtr
    = dynamic_cast<const UnknownFuncElement*>(a[0].ptr().get());

  TEUCHOS_TEST_FOR_EXCEPT(aPtr==0);

  Deriv d = funcDeriv(aPtr);
  MultipleDeriv rtn;
  rtn.put(d);
  return rtn;
}


MultipleDeriv makeDeriv(const Expr& a, const Expr& b)
{
  const UnknownFuncElement* aPtr
    = dynamic_cast<const UnknownFuncElement*>(a[0].ptr().get());

  TEUCHOS_TEST_FOR_EXCEPT(aPtr==0);

  const UnknownFuncElement* bPtr
    = dynamic_cast<const UnknownFuncElement*>(b[0].ptr().get());

  TEUCHOS_TEST_FOR_EXCEPT(bPtr==0);

  Deriv da = funcDeriv(aPtr);
  Deriv db = funcDeriv(bPtr);
  MultipleDeriv rtn;
  rtn.put(da);
  rtn.put(db);
  return rtn;
}



MultipleDeriv makeDeriv(const Expr& a, const Expr& b, const Expr& c)
{
  const UnknownFuncElement* aPtr
    = dynamic_cast<const UnknownFuncElement*>(a[0].ptr().get());

  TEUCHOS_TEST_FOR_EXCEPT(aPtr==0);

  const UnknownFuncElement* bPtr
    = dynamic_cast<const UnknownFuncElement*>(b[0].ptr().get());

  TEUCHOS_TEST_FOR_EXCEPT(bPtr==0);

  const UnknownFuncElement* cPtr
    = dynamic_cast<const UnknownFuncElement*>(c[0].ptr().get());

  TEUCHOS_TEST_FOR_EXCEPT(cPtr==0);

  Deriv da = funcDeriv(aPtr);
  Deriv db = funcDeriv(bPtr);
  Deriv dc = funcDeriv(cPtr);
  MultipleDeriv rtn;
  rtn.put(da);
  rtn.put(db);
  rtn.put(dc);
  return rtn;
}

} // namespace Sundance
