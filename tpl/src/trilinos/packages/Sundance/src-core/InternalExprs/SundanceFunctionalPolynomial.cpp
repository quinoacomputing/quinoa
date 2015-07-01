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

#include "SundanceFunctionalPolynomial.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceExpr.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceDerivOfSymbFunc.hpp"
#include "SundanceUnaryExpr.hpp"


using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;



FunctionalPolynomial::FunctionalPolynomial(const RCP<ScalarExpr>& expr)
  : funcs_(),
    funcMultiIndices_(),
    coeffs_(),
    keys_()
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isConvertibleToPoly(expr.get()), std::logic_error, 
                     "FunctionalPolynomial ctor called with an argument that "
                     "cannot be converted to a polynomial");
  int funcID;
  MultiIndex mi;
  Deriv deriv;

  const SymbolicFuncElement* s 
    = dynamic_cast<const SymbolicFuncElement*>(expr.get());
  if (s != 0)
    {
      mi = MultiIndex();
      deriv = funcDeriv(s);
    }
  
  const DerivOfSymbFunc* d
    = dynamic_cast<const DerivOfSymbFunc*>(expr.get());
  if (d != 0)
    {
      deriv = d->representMeAsFunctionalDeriv();
    }

  MultipleDeriv mu;
  mu.put(deriv);


  if (d != 0 || s != 0)
    {
      funcs_.put(funcID, expr);
      Set<MultiIndex> miSet;
      miSet.put(mi);
      funcMultiIndices_.put(funcID, miSet);

      Expr coeff = 1.0;
      RCP<ScalarExpr> cPtr = rcp_dynamic_cast<ScalarExpr>(coeff.ptr());

      coeffs_.resize(2);
      keys_.resize(2);
      coeffs_[1].put(mu, cPtr);
      keys_[1].put(mu);
    }
  else
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
                         "impossible case in FunctionalPolynomial ctor");
    }
}


FunctionalPolynomial::FunctionalPolynomial(const Map<int, RCP<ScalarExpr> >& funcs,
                                           const Map<int, Set<MultiIndex> >& funcMultiIndices,
                                           const Array<Map<MultipleDeriv, RCP<ScalarExpr> > > & coeffs)
  : funcs_(funcs),
    funcMultiIndices_(funcMultiIndices),
    coeffs_(coeffs),
    keys_(coeffs.size())
{
  typedef Map<MultipleDeriv, RCP<ScalarExpr> > termMap;

  for (int i=0; i < coeffs_.size(); i++)
    {
      for (termMap::const_iterator 
             j = coeffs_[i].begin(); j != coeffs_[i].end(); j++)
        {
          const MultipleDeriv& key = j->first;
          keys_[i].put(key);
        }
    }
}


Set<MultipleDeriv> 
FunctionalPolynomial::internalFindW(int order, const EvalContext& context) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
                     "FunctionalPolynomial not implemented");
  return Set<MultipleDeriv> ();
}


RCP<FunctionalPolynomial> FunctionalPolynomial::
addPoly(const FunctionalPolynomial* other, int sign) const 
{
  typedef Map<MultipleDeriv, RCP<ScalarExpr> > termMap;
  Map<int, RCP<ScalarExpr> > newFuncs = funcs_;
  Map<int, Set<MultiIndex> > newFuncMultiIndices = funcMultiIndices_;
  Array<Map<MultipleDeriv, RCP<ScalarExpr> > > newCoeffs = coeffs_;

  if (other->coeffs_.size() > coeffs_.size()) 
    newCoeffs.resize(other->coeffs_.size());


  for (int i=0; i < other->coeffs_.size(); i++)
    {
      
      for (termMap::const_iterator 
             j = other->coeffs_[i].begin(); j != other->coeffs_[i].end(); j++)
        {
          const MultipleDeriv& key = j->first;
          Expr right = Expr::handle(j->second);
          Expr sum;
      
          if (newCoeffs[i].containsKey(key))
            {
              Expr left = Expr::handle(newCoeffs[i].get(key));

              if (sign > 0) sum = left + right;
              else sum = left - right;
            }
          else
            {
              if (sign > 0) sum = right;
              else sum = -right;
            }
      
          const ScalarExpr* se = dynamic_cast<const ScalarExpr*>(sum.ptr().get());
          TEUCHOS_TEST_FOR_EXCEPTION(se==0, std::logic_error,
                             "Sum could not be cast to scalar expression");
          newCoeffs[i].put(key, rcp_dynamic_cast<ScalarExpr>(sum.ptr()));
        }
    }
  
  for (Map<int, RCP<ScalarExpr> >::const_iterator 
         i = other->funcs_.begin(); i != other->funcs_.end(); i++)
    {
      newFuncs.put(i->first, i->second);
    }

  for (Map<int, Set<MultiIndex> >::const_iterator 
         i = other->funcMultiIndices_.begin(); 
       i != other->funcMultiIndices_.end(); i++)
    {
      newFuncMultiIndices.put(i->first, i->second);
    }

  

  return rcp(new FunctionalPolynomial(newFuncs, newFuncMultiIndices, newCoeffs));
}

RCP<FunctionalPolynomial> FunctionalPolynomial::
multiplyPoly(const FunctionalPolynomial* other) const 
{
  typedef Map<MultipleDeriv, RCP<ScalarExpr> > termMap;
  Map<int, RCP<ScalarExpr> > newFuncs;
  Map<int, Set<MultiIndex> > newFuncMultiIndices;
  Array<Map<MultipleDeriv, RCP<ScalarExpr> > > newCoeffs;

  newCoeffs.resize(coeffs_.size() + other->coeffs_.size() - 1);

  for (int i=0; i < coeffs_.size(); i++)
    {
      for (int j = 0; j<other->coeffs_.size(); j++)
        {
          for (termMap::const_iterator 
                 me = coeffs_[i].begin(); me != coeffs_[i].end(); me++)
            {
              const MultipleDeriv& myKey = me->first;
              Expr myExpr = Expr::handle(me->second);
              for (termMap::const_iterator 
                     you = other->coeffs_[j].begin(); 
                   you != other->coeffs_[j].end(); you++)
                {
                  const MultipleDeriv& yourKey = you->first;
                  Expr yourExpr = Expr::handle(you->second);
                  MultipleDeriv newKey = myKey.product(yourKey);
                  int order = i+j;
                  Expr newTerm = myExpr*yourExpr;
                  if (newCoeffs[order].containsKey(newKey))
                    {
                      Expr old = Expr::handle(newCoeffs[i].get(newKey));
                      newTerm = newTerm + old;
                    }
                  const ScalarExpr* se 
                    = dynamic_cast<const ScalarExpr*>(newTerm.ptr().get());
                  TEUCHOS_TEST_FOR_EXCEPTION(se==0, std::logic_error,
                                     "New coeff could not be cast to scalar expression");
                  newCoeffs[order].put(newKey, 
                                       rcp_dynamic_cast<ScalarExpr>(newTerm.ptr()));
                }
            }
        }
    }
  
  for (Map<int, RCP<ScalarExpr> >::const_iterator 
         i = funcs_.begin(); i != funcs_.end(); i++)
    {
      newFuncs.put(i->first, i->second);
    }
  for (Map<int, RCP<ScalarExpr> >::const_iterator 
         i = other->funcs_.begin(); i != other->funcs_.end(); i++)
    {
      newFuncs.put(i->first, i->second);
    }

  for (Map<int, Set<MultiIndex> >::const_iterator 
         i = funcMultiIndices_.begin(); 
       i != funcMultiIndices_.end(); i++)
    {
      newFuncMultiIndices.put(i->first, i->second);
    }
  for (Map<int, Set<MultiIndex> >::const_iterator 
         i = other->funcMultiIndices_.begin(); 
       i != other->funcMultiIndices_.end(); i++)
    {
      Set<MultiIndex> miSet = i->second;
      if (newFuncMultiIndices.containsKey(i->first))
        {
          miSet.merge(newFuncMultiIndices.get(i->first));
        }
      newFuncMultiIndices.put(i->first, miSet);
    }

  return rcp(new FunctionalPolynomial(newFuncs, newFuncMultiIndices, newCoeffs));
}

RCP<FunctionalPolynomial> FunctionalPolynomial::
addFunction(const RCP<ScalarExpr>& u, int sign) const 
{
  RCP<FunctionalPolynomial> other = rcp(new FunctionalPolynomial(u));
  return addPoly(other.get(), sign);
}

RCP<FunctionalPolynomial> FunctionalPolynomial::
multiplyScalar(const RCP<ScalarExpr>& alpha) const 
{
  typedef Map<MultipleDeriv, RCP<ScalarExpr> > termMap;

  Array<Map<MultipleDeriv, RCP<ScalarExpr> > > newCoeffs = coeffs_;

  Expr alphaExpr = Expr::handle(alpha);

  for (int i=0; i < coeffs_.size(); i++)
    {
      for (termMap::const_iterator 
             j = coeffs_[i].begin(); j != coeffs_[i].end(); j++)
        {
          const MultipleDeriv& key = j->first;
          Expr oldCoeff = Expr::handle(j->second);
          Expr newCoeff = alphaExpr*oldCoeff;

          const ScalarExpr* se 
            = dynamic_cast<const ScalarExpr*>(newCoeff.ptr().get());
          TEUCHOS_TEST_FOR_EXCEPTION(se==0, std::logic_error,
                             "Coefficient could not be cast to "
                             "scalar expression");
      
          newCoeffs[i].put(key, rcp_dynamic_cast<ScalarExpr>(newCoeff.ptr()));
        }
    }
  
  return rcp(new FunctionalPolynomial(funcs_, funcMultiIndices_, newCoeffs));
}


Evaluator* FunctionalPolynomial::createEvaluator(const EvaluatableExpr* expr,
                                                 const EvalContext& context) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "poly eval not ready");
  return (Evaluator*) 0;
}



bool FunctionalPolynomial::isConvertibleToPoly(const ScalarExpr* expr)
{
  const SymbolicFuncElement* s 
    = dynamic_cast<const SymbolicFuncElement*>(expr);
  
  const DerivOfSymbFunc* d
    = dynamic_cast<const DerivOfSymbFunc*>(expr);

  const FunctionalPolynomial* p
    = dynamic_cast<const FunctionalPolynomial*>(expr);

  return (s != 0 || d != 0 || p != 0);
}


RCP<FunctionalPolynomial> FunctionalPolynomial::toPoly(const RCP<ScalarExpr>& expr)
{
  const FunctionalPolynomial* p
    = dynamic_cast<const FunctionalPolynomial*>(expr.get());

  if (p != 0) return rcp_dynamic_cast<FunctionalPolynomial>(expr);

  Expr rtn = new FunctionalPolynomial(expr);
  return rcp_dynamic_cast<FunctionalPolynomial>(rtn.ptr());
}


std::ostream& FunctionalPolynomial::toText(std::ostream& os, bool paren) const
{
  os << evalString();
  return os;
}


XMLObject FunctionalPolynomial::toXML() const
{
  XMLObject rtn("Polynomial");
  for (int order=0; order<coeffs_.size(); order++)
    {
      for (Map<MultipleDeriv, RCP<ScalarExpr> >::const_iterator
             i = coeffs_[order].begin(); i != coeffs_[order].end(); i++)
        {
          const MultipleDeriv& key = i->first;
          Expr e = Expr::handle(i->second);
          XMLObject term("Term");
          XMLObject coeff("Coeff");
          coeff.addChild(e.toXML());
          term.addChild(coeff);
          for (MultipleDeriv::const_iterator
                 j = key.begin(); j != key.end(); j++)
            {
              XMLObject f("FunctionalDeriv");
              f.addAttribute("mu", j->toString());
              term.addChild(f);
            }
          rtn.addChild(term);
        }
    }
  return rtn;
}

Set<Deriv> FunctionalPolynomial
::findFuncsForSummation(const Set<MultipleDeriv>& prevSet,
                        const MultipleDeriv& currentDeriv) const
{
  Set<Deriv> rtn;

  for (Set<MultipleDeriv>::const_iterator
         i = prevSet.begin(); i != prevSet.end(); i++)
    {
      const MultipleDeriv& mdPrev = *i;
      TEUCHOS_TEST_FOR_EXCEPTION(currentDeriv.size()+1 != mdPrev.size(),
                         std::logic_error,
                         "deriv orders must differ by 1. Found "
                         "currentDeriv.size()=" << currentDeriv.size()
                         << " while mdPrev.size()=" << mdPrev.size());

      /* We are looking for cases where the previous multiple derivative
       * is equal to the current plus one *greater* element. In such cases, the
       * set difference will contain exactly one element, and that element
       * will be greater than or equal to the the upper bound of the current 
       * set */
      Set<Deriv> tmp;
      set_difference(mdPrev.begin(), mdPrev.end(),
                     currentDeriv.begin(), currentDeriv.end(),
                     inserter(tmp, tmp.begin()));
      if (tmp.size()==1)
        {
          const Deriv& d = *(tmp.begin());
          if (currentDeriv.upper_bound(d) == currentDeriv.end()) rtn.put(d);
        }
    }
  return rtn;
}


MultipleDeriv FunctionalPolynomial::successorTerm(const MultipleDeriv& md) const
{
  MultipleDeriv rtn;

  unsigned int k = 0;
  for (MultipleDeriv::const_iterator i=md.begin(); i!=md.end(); i++, k++)
    {
      if (k < md.size()-1) rtn.put(*i);
    }
  return rtn;
}

void FunctionalPolynomial
::stepRecurrence(int level, const Map<MultipleDeriv, std::string>& sPrev,
                 Map<MultipleDeriv, std::string>& sCurr) const 
{
  Set<MultipleDeriv> allKeys;
  Set<MultipleDeriv> inducedKeys;
  Set<MultipleDeriv> prevKeys;
  for (Map<MultipleDeriv, std::string>::const_iterator 
         i = sPrev.begin(); i != sPrev.end(); i++)
    {
      inducedKeys.put(successorTerm(i->first));
    }
  for (Map<MultipleDeriv, std::string>::const_iterator 
         i = sPrev.begin(); i != sPrev.end(); i++)
    {
      prevKeys.put(i->first);
    }

  Out::os() << "keys from prev level are " << prevKeys << std::endl;
  Out::os() << "induced keys are " << inducedKeys << std::endl;
  Out::os() << "natural keys are " << keys_[level] << std::endl;

  allKeys.merge(inducedKeys);
  allKeys.merge(keys_[level]);

  Out::os() << "keys active at this level are " << allKeys << std::endl;

  for (Set<MultipleDeriv>::const_iterator 
         i = allKeys.begin(); i != allKeys.end(); i++)
    {
      const MultipleDeriv& key = *i;
      Out::os() << "working on key " << key << std::endl;
      std::string str;
      if (coeffs_[level].containsKey(key))
        {
          str = coeffs_[level].get(key)->toString();
        }

      Set<Deriv> funcs = findFuncsForSummation(prevKeys, key);
      Out::os() << "funcs to sum are " << funcs << std::endl;
      for (Set<Deriv>::const_iterator 
             j = funcs.begin(); j != funcs.end(); j++)
        {
          MultipleDeriv prevKey = key;
          Out::os() << "prev key = " << prevKey << std::endl;
          Out::os() << "func = " << *j << std::endl;
          // if (key.size() > 0 && key.upper_bound(*j) == key.end()) 
          //  {
          //    Out::os() << "skipping" << std::endl;
          //    continue;
          // }
          prevKey.put(*j);
          TEUCHOS_TEST_FOR_EXCEPTION(!sPrev.containsKey(prevKey), std::logic_error,
                             "inconsisent key lookup");
          const std::string& prevStr = sPrev.get(prevKey);
          std::string funcStr = (*j).toString();
          if (str.length() > 0) str += " + ";
          str += funcStr + "*(" + prevStr + ")";
        }
      if (str.length() > 0) sCurr.put(key, str);
    }
}

string FunctionalPolynomial::evalString() const
{
  int maxLevel = coeffs_.size()-1;

  Map<MultipleDeriv, std::string> sPrev;
  Map<MultipleDeriv, std::string> sCurr;

  for (int level=maxLevel; level>=0; level--)
    {
      Out::os() << "********* Recurrence level = " << level << " ***************"
           << std::endl;      
      sCurr = Map<MultipleDeriv, std::string>();
      stepRecurrence(level, sPrev, sCurr);
      sPrev = sCurr;

      for (Map<MultipleDeriv, std::string>::const_iterator 
             j = sCurr.begin(); j != sCurr.end(); j++)
        {
          Out::os() << "key=" << j->first << std::endl;
          Out::os() << "value=" << j->second << std::endl;
        }
    }

  //  TEUCHOS_TEST_FOR_EXCEPTION(sCurr.size() != 1, std::logic_error,
  //                     "final value should have only one element");

  return sCurr.begin()->second;
}



bool FunctionalPolynomial::lessThan(const ScalarExpr* other) const
{
  const FunctionalPolynomial* f = dynamic_cast<const FunctionalPolynomial*>(other);
  TEUCHOS_TEST_FOR_EXCEPTION(f==0, std::logic_error, "cast should never fail at this point");
  
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "FunctionalPolynomial::lessThan() not "
                     "implemented");
}

  
