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

#include "SundanceDiffOpEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceDiffOp.hpp"


#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFuncEvaluator.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceZeroExpr.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;





DiffOpEvaluator
::DiffOpEvaluator(const DiffOp* expr,
  const EvalContext& context)
  : UnaryEvaluator<DiffOp>(expr, context),
    isConstant_(this->sparsity()->numDerivs()),
    resultIndices_(this->sparsity()->numDerivs()),
    constantMonomials_(this->sparsity()->numDerivs()),
    vectorMonomials_(this->sparsity()->numDerivs()),
    constantFuncCoeffs_(this->sparsity()->numDerivs()),
    vectorFuncCoeffs_(this->sparsity()->numDerivs()),
    funcEvaluators_(),
    constantCoeffFuncIndices_(this->sparsity()->numDerivs()),
    constantCoeffFuncMi_(this->sparsity()->numDerivs()),
    vectorCoeffFuncIndices_(this->sparsity()->numDerivs()),
    vectorCoeffFuncMi_(this->sparsity()->numDerivs())
{
  int verb = context.setupVerbosity();
  Tabs tabs;
  SUNDANCE_MSG1(verb, tabs << "initializing diff op evaluator for " 
    << expr->toString());

  {
    Tabs tab0;
  
    SUNDANCE_MSG2(verb, tab0 << "return sparsity " << std::endl << *(this->sparsity)());

    SUNDANCE_MSG2(verb, tab0 << "argument sparsity subset " << std::endl 
      << *(argSparsitySuperset()));

    Map<const DiscreteFuncElementEvaluator*, int> funcToIndexMap;

    int vecResultIndex = 0;
    int constResultIndex = 0;
  
    for (int i=0; i<this->sparsity()->numDerivs(); i++)
    {
      Tabs tab1;
      const MultipleDeriv& resultDeriv = this->sparsity()->deriv(i);
      SUNDANCE_MSG3(verb, tab0 << "working out procedure for computing " 
        << resultDeriv);

      if (this->sparsity()->state(i)==ConstantDeriv)
      {
        Tabs tab2;
        addConstantIndex(i, constResultIndex);
        resultIndices_[i] = constResultIndex++;
        isConstant_[i] = true;
        SUNDANCE_MSG3(verb, tab2 << "deriv is constant, will be stored at index "
          << resultIndices_[i] << " in the const result array");
      }
      else
      {
        Tabs tab2;
        addVectorIndex(i, vecResultIndex);
        resultIndices_[i] = vecResultIndex++;
        isConstant_[i] = false;
        SUNDANCE_MSG3(verb, tab2 << "deriv is variable, will be stored at index "
          << resultIndices_[i] << " in the var result array");
      }

      int order = resultDeriv.order();
      const Set<MultipleDeriv>& RArg 
        = argExpr()->findR(order, context);
      const Set<MultipleDeriv>& RArgPlus
        = argExpr()->findR(order+1, context);
      const Set<MultipleDeriv>& W1Arg 
        = argExpr()->findW(1, context);

        
      SUNDANCE_MSG3(verb, tab1 << "RArg = " << RArg);
      SUNDANCE_MSG3(verb, tab1 << "RArgPlus = " << RArgPlus);
      SUNDANCE_MSG3(verb, tab1 << "W1Arg = " << W1Arg);

      Set<MultipleDeriv> funcTermCoeffs 
        = RArgPlus.intersection(increasedDerivs(resultDeriv, W1Arg, verb));
      SUNDANCE_MSG3(verb, tab1 << "function term coeffs = " << funcTermCoeffs);

      
      if (funcTermCoeffs.size()==0)
      {
        SUNDANCE_MSG3(verb, tab1 << "no direct chain rule terms");
      }
      else
      {
        SUNDANCE_MSG3(verb, tab1 << "getting direct chain rule terms");
      }


      for (Set<MultipleDeriv>::const_iterator 
             j=funcTermCoeffs.begin(); j != funcTermCoeffs.end(); j++)
      {
        Tabs tab2;
        SUNDANCE_MSG3(verb, tab2 << "getting coefficient of " << *j);

        int argIndex = argSparsitySuperset()->getIndex(*j);
        TEUCHOS_TEST_FOR_EXCEPTION(argIndex==-1, std::runtime_error,
          "Derivative " << *j << " expected in argument "
          "but not found");

        Deriv lambda = remainder(*j, resultDeriv, verb);

        if (lambda.isCoordDeriv())
        {
          Tabs tab3;
          SUNDANCE_MSG3(verb, tab2 << "detected coordinate deriv");
          if (lambda.coordDerivDir()!=expr->mi().firstOrderDirection())
          {
            SUNDANCE_MSG3(verb, tab2 << "direction mismatch, skipping");
            continue;
          }
          const DerivState& argState = argSparsitySuperset()->state(argIndex);
          if (argState==ConstantDeriv)
          {
            int constArgIndex = argEval()->constantIndexMap().get(argIndex);
            constantMonomials_[i].append(constArgIndex);
          }
          else
          {
            int vectorArgIndex = argEval()->vectorIndexMap().get(argIndex);
            vectorMonomials_[i].append(vectorArgIndex);
          }
        }
        else if (lambda.opOnFunc().isPartial() || lambda.opOnFunc().isIdentity())
        {
          Tabs tab3;
          SUNDANCE_MSG3(verb, tab3 << "detected functional deriv " << lambda);
          const SymbolicFuncElement* f = lambda.symbFuncElem();
          const MultiIndex& mi = expr->mi() + lambda.opOnFunc().mi(); 
          SUNDANCE_MSG3(verb, tab3 << "modified multiIndex is " << mi);

          const TestFuncElement* t 
            = dynamic_cast<const TestFuncElement*>(f);
          if (t != 0) continue;

          const UnknownFuncElement* u 
            = dynamic_cast<const UnknownFuncElement*>(f);
          TEUCHOS_TEST_FOR_EXCEPTION(u==0, std::logic_error,
            "Non-unknown function detected where an unknown "
            "function was expected in "
            "DiffOpEvaluator ctor");


          const EvaluatableExpr* evalPt = u->evalPt();
          const ZeroExpr* z = dynamic_cast<const ZeroExpr*>(evalPt);
          if (z != 0) continue;
          TEUCHOS_TEST_FOR_EXCEPTION(z != 0, std::logic_error,
            "DiffOpEvaluator detected identically zero "
            "function");

          const DiscreteFuncElement* df 
            = dynamic_cast<const DiscreteFuncElement*>(evalPt);
          
          TEUCHOS_TEST_FOR_EXCEPTION(df==0, std::logic_error,
            "DiffOpEvaluator ctor: evaluation point of "
            "unknown function " << u->toString() 
            << " is not a discrete function");

          const SymbolicFuncElementEvaluator* uEval 
            = dynamic_cast<const SymbolicFuncElementEvaluator*>(u->evaluator(context).get());

          const DiscreteFuncElementEvaluator* dfEval = uEval->dfEval();


          TEUCHOS_TEST_FOR_EXCEPTION(dfEval==0, std::logic_error,
            "DiffOpEvaluator ctor: evaluator for "
            "evaluation point is not a "
            "DiscreteFuncElementEvaluator");

          TEUCHOS_TEST_FOR_EXCEPTION(!dfEval->hasMultiIndex(mi), std::logic_error,
            "DiffOpEvaluator ctor: evaluator for "
            "discrete function " << df->toString()
            << " does not know about multiindex "
            << mi.toString());
          
          int fIndex;
          int miIndex = dfEval->miIndex(mi);
          
          if (funcToIndexMap.containsKey(dfEval))
          {
            fIndex = funcToIndexMap.get(dfEval);
          }
          else
          {
            fIndex = funcEvaluators_.size();
            funcEvaluators_.append(dfEval);
            funcToIndexMap.put(dfEval, fIndex);
          }

            
          const DerivState& argState = argSparsitySuperset()->state(argIndex);
          if (argState==ConstantDeriv)
          {
            int constArgIndex = argEval()->constantIndexMap().get(argIndex);
            constantCoeffFuncIndices_[i].append(fIndex);
            constantCoeffFuncMi_[i].append(miIndex);
            constantFuncCoeffs_[i].append(constArgIndex);
          }
          else
          {
            int vectorArgIndex = argEval()->vectorIndexMap().get(argIndex);
            vectorCoeffFuncIndices_[i].append(fIndex);
            vectorCoeffFuncMi_[i].append(miIndex);
            vectorFuncCoeffs_[i].append(vectorArgIndex);
          }
        }
        else
        {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
            "DiffOpEvaluator has been asked to preprocess a Deriv that "
            "is not a simple partial derivative. The problem child is: "
            << lambda);
        }
      }
      
      
      Set<MultipleDeriv> isolatedTerms 
        = RArg.intersection(backedDerivs(resultDeriv, W1Arg, verb));
      
      if (isolatedTerms.size()==0)
      {
        SUNDANCE_MSG3(verb, tab1 << "no indirect chain rule terms");
      }
      else
      {
        SUNDANCE_MSG3(verb, tab1 << "getting indirect chain rule terms");
        SUNDANCE_MSG3(verb, tab1 << "isolated terms = " << isolatedTerms);
      }

      for (Set<MultipleDeriv>::const_iterator 
             j=isolatedTerms.begin(); j != isolatedTerms.end(); j++)
      {
        int argIndex = argSparsitySuperset()->getIndex(*j);
        TEUCHOS_TEST_FOR_EXCEPTION(argIndex==-1, std::runtime_error,
          "Derivative " << *j << " expected in argument "
          "but not found");
        const DerivState& argState = argSparsitySuperset()->state(argIndex);
        if (argState==ConstantDeriv)
        {
          int constArgIndex = argEval()->constantIndexMap().get(argIndex);
          constantMonomials_[i].append(constArgIndex);
        }
        else
        {
          int vectorArgIndex = argEval()->vectorIndexMap().get(argIndex);
          vectorMonomials_[i].append(vectorArgIndex);
        }
      }
    }
  }

  if (verb > 2)
  {
    Out::os() << tabs << "instruction tables for summing spatial/functional chain rule" << std::endl;
    for (int i=0; i<this->sparsity()->numDerivs(); i++)
    {
      Tabs tab1;
      Out::os() << tab1 << "deriv " << sparsity()->deriv(i) << std::endl;
      {
        Tabs tab2;
        Out::os() << tab2 << "constant monomials: " << constantMonomials_[i]
                  << std::endl;
        Out::os() << tab2 << "vector monomials: " << vectorMonomials_[i]
                  << std::endl;
            
        Out::os() << tab2 << "constant coeff functions: " << std::endl;
        for (int j=0; j<constantFuncCoeffs_[i].size(); j++)
        {
          Tabs tab3;
          Out::os() << tab3 << "func=" << constantCoeffFuncIndices_[i][j]
                    << " mi=" << constantCoeffFuncMi_[i][j] << std::endl;
        } 
        Out::os() << tab2 << "vector coeff functions: " << std::endl;
        for (int j=0; j<vectorFuncCoeffs_[i].size(); j++)
        {
          Tabs tab3;
          Out::os() << tab3 << "func=" << vectorCoeffFuncIndices_[i][j]
                    << " mi=" << vectorCoeffFuncMi_[i][j] << std::endl;
        }
            
      }
    }
  }
}


Deriv DiffOpEvaluator::remainder(const MultipleDeriv& big, 
  const MultipleDeriv& little, int verb) const 
{
  Tabs tab;
  SUNDANCE_MSG5(verb, tab << "computing remainder: big=" << big << ", little="
    << little);
  TEUCHOS_TEST_FOR_EXCEPT(big.order()-little.order() != 1);

  MultipleDeriv r;
  if (little.order()==0) r = big;
  else r = big.factorOutDeriv(little);

  SUNDANCE_MSG5(verb, tab << "remainder = " << r);

  TEUCHOS_TEST_FOR_EXCEPT(r.order() != 1);

  return *(r.begin());
}

Set<MultipleDeriv> DiffOpEvaluator
::increasedDerivs(const MultipleDeriv& mu,
  const Set<MultipleDeriv>& W1, int verb) const
{
  Tabs tabs;
  SUNDANCE_MSG3(verb, tabs << "computing increased derivs");
  Set<MultipleDeriv> rtn;
  for (Set<MultipleDeriv>::const_iterator i=W1.begin(); i!=W1.end(); i++)
  {
    MultipleDeriv md = *i;
    TEUCHOS_TEST_FOR_EXCEPT(md.order() != 1);
    Deriv lambda = *(md.begin());
    MultipleDeriv lambdaMu = mu;
    lambdaMu.put(lambda);
    rtn.put(lambdaMu);
  }
  SUNDANCE_MSG3(verb, tabs << "increased derivs = " << rtn);
  return rtn;
}

Set<MultipleDeriv> DiffOpEvaluator
::backedDerivs(const MultipleDeriv& mu,
  const Set<MultipleDeriv>& W1, int verb) const
{
  Tabs tabs;
  SUNDANCE_MSG3(verb, tabs << "computing backed-out derivs for mu= " << mu
    << ", W1=" << W1);
  Set<MultipleDeriv> rtn;
  if (mu.order() != 0) 
  {
    const MultiIndex& alpha = expr()->mi();

    for (Set<MultipleDeriv>::const_iterator i=W1.begin(); i!=W1.end(); i++)
    {
      const MultipleDeriv& md = *i;
      TEUCHOS_TEST_FOR_EXCEPT(md.order() != 1);
      Deriv lambda = *(md.begin());
      if (lambda.isCoordDeriv()) continue;
      TEUCHOS_TEST_FOR_EXCEPT(!lambda.isFunctionalDeriv());
      FunctionIdentifier lambda_fid = lambda.fid();
      const MultiIndex& lambda_mi = lambda.opOnFunc().mi(); 
      for (MultipleDeriv::const_iterator j=mu.begin(); j!=mu.end(); j++)
      {
        const Deriv& d = *j;
        if (d.isCoordDeriv()) continue;
        FunctionIdentifier d_fid = d.fid();
        const MultiIndex& d_mi = d.opOnFunc().mi(); 
        if (d_fid != lambda_fid) continue;
        if (!(alpha + lambda_mi == d_mi)) continue;
        MultipleDeriv z = mu.factorOutDeriv(d);
        z.put(lambda);
        rtn.put(z);
      }
    }
  }
  SUNDANCE_MSG3(verb, tabs << "backed-out derivs = " << rtn);
  return rtn;
}



void DiffOpEvaluator::internalEval(const EvalManager& mgr,
  Array<double>& constantResults,
  Array<RCP<EvalVector> >& vectorResults)  const
{
  Tabs tabs;
  SUNDANCE_MSG1(mgr.verb(), tabs << "DiffOpEvaluator::eval() expr=" 
    << expr()->toString());

  /* evaluate the argument */
  Array<RCP<EvalVector> > argVectorResults;
  Array<double> argConstantResults;

  SUNDANCE_MSG2(mgr.verb(), tabs << "evaluating operand");
  evalOperand(mgr, argConstantResults, argVectorResults);


  if (mgr.verb() > 2)
  {
    Tabs tab1;
    Out::os() << tabs << "DiffOp operand results" << std::endl;
    mgr.showResults(Out::os(), argSparsitySuperset(), argVectorResults,
      argConstantResults);
  }



  /* evaluate the required discrete functions */
  SUNDANCE_MSG2(mgr.verb(), tabs << "evaluating discrete functions, num funcs= " << funcEvaluators_.size());

  Array<Array<RCP<EvalVector> > > funcVectorResults(funcEvaluators_.size());
  Array<double> funcConstantResults;
  for (int i=0; i<funcEvaluators_.size(); i++)
  {
    funcEvaluators_[i]->eval(mgr, funcConstantResults, funcVectorResults[i]);
  }
  
  constantResults.resize(this->sparsity()->numConstantDerivs());
  vectorResults.resize(this->sparsity()->numVectorDerivs());
  
  SUNDANCE_MSG3(mgr.verb(), tabs << "summing spatial/functional chain rule");

  for (int i=0; i<this->sparsity()->numDerivs(); i++)
  {
    Tabs tab1;
    SUNDANCE_MSG4(mgr.verb(), tab1 << "working on deriv " 
      << this->sparsity()->deriv(i));

    /* add constant monomials */
    SUNDANCE_MSG4(mgr.verb(), tab1 << "have " <<  constantMonomials_[i].size()
      << " constant monomials");
    double constantVal = 0.0;
    for (int j=0; j<constantMonomials_[i].size(); j++)
    {
      SUNDANCE_MSG4(mgr.verb(), tab1 << "adding in constant monomial (index "
        << constantMonomials_[i][j] 
        << " in arg results)");
      constantVal += argConstantResults[constantMonomials_[i][j]];
    }
    if (isConstant_[i])
    {
      constantResults[resultIndices_[i]] = constantVal;
      SUNDANCE_MSG4(mgr.verb(), tab1 << "result is constant: value=" 
        << constantVal);
      continue;
    }

    RCP<EvalVector> result;
    bool vecHasBeenAllocated = false;

    /* add in the vector monomials */
    const Array<int>& vm = vectorMonomials_[i];
    SUNDANCE_MSG4(mgr.verb(), tab1 << "have " << vm.size() 
      << " vector monomials");
    for (int j=0; j<vm.size(); j++)
    {
      Tabs tab2;

      const RCP<EvalVector>& v = argVectorResults[vm[j]];

      SUNDANCE_MSG4(mgr.verb(), tab2 << "found vec monomial term " << v->str());

      /* if we've not yet allocated a vector for the results, 
       * do so now, and set it to the initial value */ 
      if (!vecHasBeenAllocated)
      {
        SUNDANCE_MSG4(mgr.verb(), tab2 << "allocated result vector");
        result = mgr.popVector();
        vecHasBeenAllocated = true;
        if (isZero(constantVal))
        {
          result->setTo_V(v.get());
        }
        else
        {
          result->setTo_S_add_V(constantVal, v.get());
        }
      }
      else
      {
        result->add_V(v.get());
      }
      SUNDANCE_MSG4(mgr.verb(), tab2 << "result is " << result->str());
    }
      
    /* add in the function terms with constant coeffs */
    const Array<int>& cf = constantFuncCoeffs_[i];
    SUNDANCE_MSG4(mgr.verb(), tab1 << "adding " << cf.size()
      << " func terms with constant coeffs");
    for (int j=0; j<cf.size(); j++)
    {
      Tabs tab2;
      const double& coeff = argConstantResults[cf[j]];
      int fIndex = constantCoeffFuncIndices_[i][j];
      int miIndex = constantCoeffFuncMi_[i][j];
      const RCP<EvalVector>& fValue 
        = funcVectorResults[fIndex][miIndex];

      SUNDANCE_MSG4(mgr.verb(), tab2 << "found term: coeff= " 
        << coeff << ", func value=" << fValue->str());
          
      /* if we've not yet allocated a vector for the results, 
       * do so now, and set it to the initial value */ 
      if (!vecHasBeenAllocated)
      {
        SUNDANCE_MSG4(mgr.verb(), tab2 << "allocated result vector");
        result = mgr.popVector();
        vecHasBeenAllocated = true;
        if (isOne(coeff))
        {
          if (isZero(constantVal))
          {
            result->setTo_V(fValue.get());
          }
          else
          {
            result->setTo_S_add_V(constantVal, fValue.get());
          }
        }
        else
        {
          if (isZero(constantVal))
          {
            result->setTo_SV(coeff, fValue.get());
          }
          else
          {
            result->setTo_S_add_SV(constantVal, coeff, fValue.get());
          }
        }
      }
      else
      {
        if (isOne(coeff))
        {
          result->add_V(fValue.get());
        }
        else
        {
          result->add_SV(coeff, fValue.get());
        }
      }
      SUNDANCE_MSG4(mgr.verb(), tab2 << "result is " << result->str());
    }

      
    /* add in the function terms with vector coeffs */
    const Array<int>& vf = vectorFuncCoeffs_[i];
    SUNDANCE_MSG4(mgr.verb(), tab1 << "adding " << vf.size()
      << " func terms with vector coeffs");
    for (int j=0; j<vf.size(); j++)
    {
      Tabs tab2;

      const RCP<EvalVector>& coeff = argVectorResults[vf[j]];
      int fIndex = vectorCoeffFuncIndices_[i][j];
      int miIndex = vectorCoeffFuncMi_[i][j];
      const RCP<EvalVector>& fValue 
        = funcVectorResults[fIndex][miIndex];

      SUNDANCE_MSG4(mgr.verb(), tab2 << "found term: coeff= " 
        << coeff->str() << ", func value=" 
        << fValue->str());
          
      /* if we've not yet allocated a vector for the results, 
       * do so now, and set it to the initial value */ 
      if (!vecHasBeenAllocated)
      {
        SUNDANCE_MSG4(mgr.verb(), tab2 << "allocated result vector");
        result = mgr.popVector();
        vecHasBeenAllocated = true;
        result->setTo_VV(coeff.get(), fValue.get());
      }
      else
      {
        result->add_VV(coeff.get(), fValue.get());
      }
      SUNDANCE_MSG4(mgr.verb(), tab2 << "result is " << result->str());
    }

    TEUCHOS_TEST_FOR_EXCEPTION(!vecHasBeenAllocated, std::logic_error,
      "created empty vector in DiffOpEvaluator::internalEval");
    vectorResults[resultIndices_[i]] = result;
  }

  if (mgr.verb() > 1)
  {
    Out::os() << tabs << "diff op results" << std::endl;
    mgr.showResults(Out::os(), sparsity(), vectorResults,
      constantResults);
  }
  SUNDANCE_MSG1(mgr.verb(), tabs << "done spatial/functional chain rule");
}



void DiffOpEvaluator::resetNumCalls() const 
{
  argEval()->resetNumCalls();
  for (int i=0; i<funcEvaluators_.size(); i++) 
  {
    funcEvaluators_[i]->resetNumCalls();
  }
  Evaluator::resetNumCalls();
}
