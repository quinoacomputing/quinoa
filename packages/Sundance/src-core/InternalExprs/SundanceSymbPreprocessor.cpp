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

#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceExpr.hpp"
#include "PlayaTabs.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceUnknownParameterElement.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"

#include "SundanceOut.hpp"
#include "Teuchos_Utils.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;





DerivSet SymbPreprocessor::setupFwdProblem(const Expr& expr, 
  const Expr& tests,
  const Expr& unks,
  const Expr& unkEvalPts, 
  const Expr& unkParams,
  const Expr& unkParamEvalPts,
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,
  const Expr& fixedFields,
  const Expr& fixedFieldEvalPts,
  const EvalContext& context,
  const ComputationType& compType)
{
  Expr zero;
  Expr v = tests.flatten();
  Array<Expr> z(v.size());
  for (int i=0; i<v.size(); i++) z[i] = new ZeroExpr();
  zero = new ListExpr(z);

  return setupVariations(expr, 
    tests, zero,
    unks, unkEvalPts,
    unkParams, unkParamEvalPts,
    fixedFields, fixedFieldEvalPts,
    fixedParams, fixedParamEvalPts,
    context, compType);
}




DerivSet SymbPreprocessor::setupSensitivities(const Expr& expr, 
  const Expr& tests,
  const Expr& unks,
  const Expr& unkEvalPts, 
  const Expr& unkParams,
  const Expr& unkParamEvalPts,
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,
  const Expr& fixedFields,
  const Expr& fixedFieldEvalPts,
  const EvalContext& context,
  const ComputationType& compType)
{
  Expr zero;
  Expr v = tests.flatten();
  Array<Expr> z(v.size());
  for (int i=0; i<v.size(); i++) z[i] = new ZeroExpr();
  zero = new ListExpr(z);

  return setupVariations(expr, 
    tests, zero,
    unks, unkEvalPts,
    unkParams, unkParamEvalPts,
    fixedFields, fixedFieldEvalPts,
    fixedParams, fixedParamEvalPts,
    context, compType);
}


DerivSet SymbPreprocessor::setupFunctional(const Expr& expr, 
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,
  const Expr& fixedFields,
  const Expr& fixedFieldEvalPts,
  const EvalContext& context,
  const ComputationType& compType)
{
  Expr vars;
  Expr varEvalPts;
  Expr unks;
  Expr unkEvalPts;
  Expr unkParams;
  Expr unkParamEvalPts;

  return setupVariations(expr, 
    vars, varEvalPts,
    unks, unkEvalPts,
    unkParams, unkParamEvalPts,
    fixedFields, fixedFieldEvalPts,
    fixedParams, fixedParamEvalPts,
    context,compType);
}





DerivSet SymbPreprocessor::setupGradient(const Expr& expr, 
  const Expr& vars,
  const Expr& varEvalPts,
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,
  const Expr& fixedFields,
  const Expr& fixedFieldEvalPts, 
  const EvalContext& context,
  const ComputationType& compType)
{
  Expr unks;
  Expr unkEvalPts;
  Expr unkParams;
  Expr unkParamEvalPts;

  return setupVariations(expr, 
    vars, varEvalPts,
    unks, unkEvalPts,
    unkParams, unkParamEvalPts,
    fixedFields, fixedFieldEvalPts,
    fixedParams, fixedParamEvalPts,
    context, compType);
}




DerivSet SymbPreprocessor::setupVariations(const Expr& expr, 
  const Expr& vars,
  const Expr& varEvalPts,
  const Expr& unks,
  const Expr& unkEvalPts,
  const Expr& unkParams,
  const Expr& unkParamEvalPts,
  const Expr& fixedFields,
  const Expr& fixedFieldEvalPts, 
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts, 
  const EvalContext& context,
  const ComputationType& compType)
{
  TimeMonitor t(preprocTimer());
  Tabs tab;

  const EvaluatableExpr* e 
    = dynamic_cast<const EvaluatableExpr*>(expr.ptr().get());

  Array<Set<MultiSet<int> > > funcDerivs(3);
  Array<Set<MultiIndex> > spatialDerivs(3);

  int verb=context.setupVerbosity();
  SUNDANCE_BANNER1(verb, tab, "in setupVariations()");
  verbosity<EvaluatableExpr>() = verb;
  SUNDANCE_MSG1(verb,
    tab << "************ setting up variations of expr: " 
    << expr 
    << std::endl << tab << "context is " << context 
    << std::endl << tab << "conp type is " << compType
    << std::endl << tab << "vars are " << vars
    << std::endl << tab << "unks are " << unks
    << std::endl << tab << "unk parameters " << unkParams
    << std::endl << tab << "fixed parameters " << fixedParams
    << std::endl << tab << "the eval points for the vars are " 
    << varEvalPts
    << std::endl << tab << "the eval points for the unks are " 
    << unkEvalPts
    << std::endl << tab 
    << "the eval points for the unknown parameters are " 
    << unkParamEvalPts 
    << std::endl << tab 
    << "the eval points for the fixed parameters are " 
    << fixedParamEvalPts 
    << tab << std::endl);

  TEUCHOS_TEST_FOR_EXCEPTION(e==0, std::logic_error,
    "Non-evaluatable expr " << expr.toString()
    << " given to SymbPreprocessor::setupExpr()");

  /* make flat lists of variations, unknowns, parameters, and fixed fields */
  Expr v = vars.flatten();
  Expr v0 = varEvalPts.flatten();
  Expr u = unks.flatten();
  Expr u0 = unkEvalPts.flatten();
  Expr alpha = unkParams.flatten();
  Expr alpha0 = unkParamEvalPts.flatten();
  Expr beta = fixedParams.flatten();
  Expr beta0 = fixedParamEvalPts.flatten();
  Expr f = fixedFields.flatten();
  Expr f0 = fixedFieldEvalPts.flatten();
  


  Set<int> varID = processInputFuncs<SymbolicFuncElement>(v, v0);

  Set<int> unkID = processInputFuncs<UnknownFuncElement>(u, u0);

  Set<int> fixedID = processInputFuncs<UnknownFuncElement>(f, f0);

  Set<int> unkParamID 
    = processInputParams<UnknownParameterElement>(alpha, alpha0);

  Set<int> fixedParamID 
    = processInputParams<UnknownParameterElement>(beta, beta0);


  

  /* put together the set of functions that are active differentiation
   * variables */

  SUNDANCE_MSG2(verb, tab << "forming active set");
  Array<Sundance::Set<MultiSet<int> > > activeFuncIDs(3);
  if (context.needsDerivOrder(0)) activeFuncIDs[0].put(MultiSet<int>());
  if (context.topLevelDiffOrder() >= 1)
  {
    for (Set<int>::const_iterator i=varID.begin(); i != varID.end(); i++)
    {
      if (context.needsDerivOrder(1)) activeFuncIDs[1].put(makeMultiSet<int>(*i));
      if (context.topLevelDiffOrder()==2)
      {
        for (Set<int>::const_iterator j=unkID.begin(); j != unkID.end(); j++)
        {
          activeFuncIDs[2].put(makeMultiSet<int>(*i, *j));
        }
        if (compType==MatrixAndVector)
        {
          for (Set<int>::const_iterator 
                 j=unkParamID.begin(); j != unkParamID.end(); j++)
          {
            activeFuncIDs[2].put(makeMultiSet<int>(*i, *j));
          }
        }
        else if (compType==Sensitivities)
        {
          for (Set<int>::const_iterator 
                 j=fixedParamID.begin(); j != fixedParamID.end(); j++)
          {
            activeFuncIDs[2].put(makeMultiSet<int>(*i, *j));
          }
        }
      }
    }
  }
  SUNDANCE_MSG1(verb, tab << std::endl << tab 
    << " ************* Finding nonzeros for expr " << std::endl << tab);
  for (int i=0; i<=context.topLevelDiffOrder(); i++)
  {
    Tabs tab2;
    SUNDANCE_MSG4(verb, tab2 << "diff order=" << i << ", active funcs="
      << activeFuncIDs[i]);
  }

  Set<MultiIndex> miSet;
  miSet.put(MultiIndex());
  e->registerSpatialDerivs(context, miSet);
  
  SUNDANCE_MSG2(verb,
    tab << std::endl << tab 
    << " ************* finding required functions" << std::endl << tab);
  SUNDANCE_MSG2(verb,
    tab << "activeFuncIDs are = " << activeFuncIDs);
  SUNDANCE_MSG2(verb,
    tab << "spatial derivs are = " << spatialDerivs);
  Array<Set<MultipleDeriv> > RInput 
    = e->computeInputR(context, activeFuncIDs, spatialDerivs);

  SUNDANCE_MSG3(verb,
    tab << std::endl << tab 
    << " ************* Top-level required funcs are " << RInput << std::endl << tab);


  SUNDANCE_MSG2(verb,
    tab << std::endl << tab 
    << " ************* Calling determineR()" << std::endl << tab);
  
  e->determineR(context, RInput);

  
  SUNDANCE_MSG1(verb,
    tab << std::endl << tab 
    << " ************* Finding sparsity structure " << std::endl << tab);
  DerivSet derivs = e->sparsitySuperset(context)->derivSet();
  SUNDANCE_MSG1(verb,
    tab << std::endl << tab 
    << "Nonzero deriv set = " << derivs);

  SUNDANCE_MSG1(verb,
    tab << std::endl << tab 
    << " ************* Setting up evaluators for expr " << std::endl << tab);

  int saveVerb = context.setupVerbosity();
  context.setSetupVerbosity(0);
  e->setupEval(context);

  if (verb>1)
  { 
    SUNDANCE_MSG1(verb,
      tab << std::endl << tab 
      << " ************* Nonzeros are:");
    e->displayNonzeros(Out::os(), context);

  }


  context.setSetupVerbosity(saveVerb);
  return derivs;
}


namespace Sundance
{
Expr makeZeros(const Expr& v)
{
  Array<Expr> z(v.size());
  for (int i=0; i<v.size(); i++) 
  {
    if (v[i].size()==1) z[i] = new ZeroExpr();
    else z[i] = makeZeros(v[i]);
  }
  return new ListExpr(z);
}
}
