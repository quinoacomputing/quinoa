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

#include "SundanceEquationSet.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceFunctionSupportResolver.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceSpectralExpr.hpp"
#include "SundanceUnknownParameterElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceSumOfIntegrals.hpp"
#include "SundanceSumOfBCs.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"

 



using namespace Sundance;
using namespace Teuchos;

EquationSet::EquationSet(const Expr& eqns, 
  const Expr& bcs, 
  const Expr& params,
  const Expr& paramEvalPts,
  const Array<Expr>& fields,
  const Array<Expr>& fieldValues)
  : fsr_(),
    varUnkPairsOnRegions_(),
    bcVarUnkPairsOnRegions_(),
    regionQuadCombos_(),
    bcRegionQuadCombos_(),
    regionQuadComboExprs_(),
    bcRegionQuadComboExprs_(),
    regionQuadComboNonzeroDerivs_(),
    bcRegionQuadComboNonzeroDerivs_(),
    rqcToContext_(),
    bcRqcToContext_(),
    rqcToSkip_(),
    bcRqcToSkip_(),
    unkLinearizationPts_(),
    unkParamEvalPts_(),
    fixedParamEvalPts_(paramEvalPts),
    compTypes_(),
    isNonlinear_(false),
    isVariationalProblem_(true),
    isFunctionalCalculator_(true),
    isSensitivityProblem_(false)
{
  Array<Expr> unks;
  Array<Expr> unkEvalPt;
  Array<Expr> vars;
  Array<Expr> varEvalPt;
  Expr unkParams;

  fsr_ = rcp(new FunctionSupportResolver(eqns, bcs, vars, unks, unkParams,
      params, fields, true));
  
  compTypes_.put(FunctionalOnly);

  rqcToContext_.put(FunctionalOnly, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(FunctionalOnly, Map<RegionQuadCombo, EvalContext>());

  regionQuadComboNonzeroDerivs_.put(FunctionalOnly,
    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(FunctionalOnly,
    Map<RegionQuadCombo, DerivSet>());

  

  init(varEvalPt, unkEvalPt, unkParamEvalPts_, paramEvalPts, fieldValues);
}

EquationSet::EquationSet(const Expr& eqns, 
  const Expr& bcs,
  const Array<Expr>& vars, 
  const Array<Expr>& unks,
  const Array<Expr>& unkLinearizationPts,
  const Expr& unkParams,
  const Expr& unkParamEvalPts, 
  const Expr& params,
  const Expr& paramEvalPts,
  const Array<Expr>& fixedFields,
  const Array<Expr>& fixedFieldValues)
  : fsr_(),
    varUnkPairsOnRegions_(),
    bcVarUnkPairsOnRegions_(),
    regionQuadCombos_(),
    bcRegionQuadCombos_(),
    regionQuadComboExprs_(),
    bcRegionQuadComboExprs_(),
    regionQuadComboNonzeroDerivs_(),
    bcRegionQuadComboNonzeroDerivs_(),
    rqcToContext_(),
    bcRqcToContext_(),
    rqcToSkip_(),
    bcRqcToSkip_(),
    unkLinearizationPts_(flattenSpectral(unkLinearizationPts)),
    unkParamEvalPts_(unkParamEvalPts),
    fixedParamEvalPts_(paramEvalPts),
    compTypes_(),
    isNonlinear_(false),
    isVariationalProblem_(false),
    isFunctionalCalculator_(false),
    isSensitivityProblem_(unkParams.size() > 0)
{
  compTypes_.put(MatrixAndVector);
  compTypes_.put(VectorOnly);

  fsr_ = rcp(new FunctionSupportResolver(eqns, bcs, vars, 
      unks, unkParams,
      params, fixedFields, false));

  rqcToContext_.put(MatrixAndVector, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(MatrixAndVector, Map<RegionQuadCombo, EvalContext>());

  rqcToContext_.put(VectorOnly, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(VectorOnly, Map<RegionQuadCombo, EvalContext>());

  regionQuadComboNonzeroDerivs_.put(MatrixAndVector, 
    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(MatrixAndVector, 
    Map<RegionQuadCombo, DerivSet>());

  regionQuadComboNonzeroDerivs_.put(VectorOnly, 
    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(VectorOnly, 
    Map<RegionQuadCombo, DerivSet>());

  if (params.size() > 0) 
  {
    compTypes_.put(Sensitivities);
    rqcToContext_.put(Sensitivities, Map<RegionQuadCombo, EvalContext>());
    bcRqcToContext_.put(Sensitivities, Map<RegionQuadCombo, EvalContext>());
    regionQuadComboNonzeroDerivs_.put(Sensitivities, 
      Map<RegionQuadCombo, DerivSet>());
    bcRegionQuadComboNonzeroDerivs_.put(Sensitivities, 
      Map<RegionQuadCombo, DerivSet>());
  }

  Array<Expr> zero;

  init(zero, flattenSpectral(unkLinearizationPts), unkParamEvalPts_, 
    paramEvalPts, fixedFieldValues);
}


EquationSet::EquationSet(const Expr& eqns, 
  const Expr& bcs,
  const Array<Expr>& vars,
  const Array<Expr>& varLinearizationPts, 
  const Array<Expr>& unks,
  const Array<Expr>& unkLinearizationPts, 
  const Expr& params,
  const Expr& paramEvalPts,
  const Array<Expr>& fixedFields,
  const Array<Expr>& fixedFieldValues)
  : fsr_(),
    varUnkPairsOnRegions_(),
    bcVarUnkPairsOnRegions_(),
    regionQuadCombos_(),
    bcRegionQuadCombos_(),
    regionQuadComboExprs_(),
    bcRegionQuadComboExprs_(),
    regionQuadComboNonzeroDerivs_(),
    bcRegionQuadComboNonzeroDerivs_(),
    rqcToContext_(),
    bcRqcToContext_(),
    rqcToSkip_(),
    bcRqcToSkip_(),
    unkLinearizationPts_(flattenSpectral(unkLinearizationPts)),
    unkParamEvalPts_(),
    fixedParamEvalPts_(paramEvalPts),
    compTypes_(),
    isNonlinear_(false),
    isVariationalProblem_(true),
    isFunctionalCalculator_(false),
    isSensitivityProblem_(false)
{
  Expr unkParams;
  fsr_ = rcp(new FunctionSupportResolver(eqns, bcs, vars, 
      unks, unkParams, params, fixedFields, 
      isVariationalProblem_));

  compTypes_.put(MatrixAndVector);
  compTypes_.put(VectorOnly);

  rqcToContext_.put(MatrixAndVector, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(MatrixAndVector, Map<RegionQuadCombo, EvalContext>());

  rqcToContext_.put(VectorOnly, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(VectorOnly, Map<RegionQuadCombo, EvalContext>());

  regionQuadComboNonzeroDerivs_.put(MatrixAndVector, 
    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(MatrixAndVector, 
    Map<RegionQuadCombo, DerivSet>());

  regionQuadComboNonzeroDerivs_.put(VectorOnly, 
    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(VectorOnly, 
    Map<RegionQuadCombo, DerivSet>());

  init(flattenSpectral(varLinearizationPts), 
    flattenSpectral(unkLinearizationPts), 
    unkParamEvalPts_, 
    paramEvalPts, fixedFieldValues);
}

EquationSet::EquationSet(const Expr& eqns, 
  const Expr& bcs,
  const Array<Expr>& vars,
  const Array<Expr>& varLinearizationPts, 
  const Expr& params,
  const Expr& paramEvalPts,
  const Array<Expr>& fixedFields,
  const Array<Expr>& fixedFieldValues)
  : fsr_(),
    varUnkPairsOnRegions_(),
    bcVarUnkPairsOnRegions_(),
    regionQuadCombos_(),
    bcRegionQuadCombos_(),
    regionQuadComboExprs_(),
    bcRegionQuadComboExprs_(),
    regionQuadComboNonzeroDerivs_(),
    bcRegionQuadComboNonzeroDerivs_(),
    rqcToContext_(),
    bcRqcToContext_(),
    rqcToSkip_(),
    bcRqcToSkip_(),
    unkLinearizationPts_(),
    unkParamEvalPts_(),
    fixedParamEvalPts_(paramEvalPts),
    compTypes_(),
    isNonlinear_(false),
    isVariationalProblem_(true),
    isFunctionalCalculator_(true),
    isSensitivityProblem_(false)
{
  compTypes_.put(FunctionalOnly);
  compTypes_.put(FunctionalAndGradient);

  Expr unkParams;
  Array<Expr> unks;
  fsr_ = rcp(new FunctionSupportResolver(eqns, bcs, vars, 
      unks, unkParams,
      params, fixedFields, isVariationalProblem_));

  rqcToContext_.put(FunctionalAndGradient, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(FunctionalAndGradient, Map<RegionQuadCombo, EvalContext>());

  rqcToContext_.put(FunctionalOnly, Map<RegionQuadCombo, EvalContext>());
  bcRqcToContext_.put(FunctionalOnly, Map<RegionQuadCombo, EvalContext>());

  regionQuadComboNonzeroDerivs_.put(FunctionalAndGradient, 
    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(FunctionalAndGradient, 
    Map<RegionQuadCombo, DerivSet>());

  regionQuadComboNonzeroDerivs_.put(FunctionalOnly, 
    Map<RegionQuadCombo, DerivSet>());
  bcRegionQuadComboNonzeroDerivs_.put(FunctionalOnly, 
    Map<RegionQuadCombo, DerivSet>());


  init(flattenSpectral(varLinearizationPts), 
    flattenSpectral(unkLinearizationPts_), 
    unkParamEvalPts_, 
    paramEvalPts, fixedFieldValues);
}




void EquationSet::init(
  const Array<Expr>& varLinearizationPts,
  const Array<Expr>& unkLinearizationPts,
  const Expr& unkParamEvalPts, 
  const Expr& fixedParamEvalPts,
  const Array<Expr>& fixedFieldValues)
{
  Tabs tab0(0);
  Tabs tab1;
  bool hasBCs = fsr_->hasBCs();
  const SumOfBCs* bcSum = fsr_->bcSum();
  const SumOfIntegrals* integralSum = fsr_->integralSum();
 
  int verb = 0;

  /* upgrade base verbosity level if one of the terms is being watched */
  if (integralSum->hasWatchedTerm() || (hasBCs && bcSum->hasWatchedTerm()))
  {
    int v1 = integralSum->eqnSetSetupVerb();
    int v2 = 0;
    if (hasBCs) v2 = bcSum->eqnSetSetupVerb();
    verb = std::max(v1, v2);
  }
  SUNDANCE_BANNER1(verb, tab0, "EquationSet setup");


  /* get symbolic functions from the function support resolver */
  const Array<Expr>& unks = fsr_->unks();
  const Array<Expr>& vars = fsr_->vars();
  const Expr& unkParams = fsr_->unkParams();
  const Expr& fixedParams = fsr_->fixedParams();
  const Array<Expr>& fixedFields = fsr_->fixedFields();

  SUNDANCE_MSG1(verb, tab0 << "fixed params = " << fixedParams);

  const Set<int>& varFuncSet = fsr_->varFuncSet();
  const Set<int>& unkFuncSet = fsr_->unkFuncSet();

  Set<RegionQuadCombo> rqcSet;
  Set<RegionQuadCombo> rqcBCSet;



  Array<int> contextID = tuple(EvalContext::nextID(),
    EvalContext::nextID(),
    EvalContext::nextID(),
    EvalContext::nextID(),
    EvalContext::nextID());


  /* for each computation type, initialize list of regions to skip */
  for (Set<ComputationType>::const_iterator 
         i=compTypes_.begin(); i!=compTypes_.end(); i++)
  {
    rqcToSkip_[*i] = Set<RegionQuadCombo>();
    bcRqcToSkip_[*i] = Set<RegionQuadCombo>();
  }

  SUNDANCE_MSG1(verb, tab0 << "computation types = " << compTypes_);



  /* Now compile a list of all regions appearing in either the eqns or
   * the BCs */

  /* Do the non-bc eqns first */
  SUNDANCE_MSG1(verb, tab1 << "processing integral terms");
  for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
         r=integralSum->rqcToExprMap().begin(); 
       r!=integralSum->rqcToExprMap().end(); r++)
  {
    Tabs tab15;
    Tabs tab2;
    RegionQuadCombo rqc = r->first;
    int rqcVerb = verb;
    int symbVerb = 0;
    int evalSetupVerb = 0;
    if (rqc.watch().isActive()) 
    {
      symbVerb = rqc.watch().param("symbolic preprocessing");
      rqcVerb=rqc.watch().param("equation set setup");
      evalSetupVerb=rqc.watch().param("evaluator setup");
    }
    SUNDANCE_MSG1(std::max(verb,rqcVerb), tab15 << "processing RQC = " << rqc);


    rqcSet.put(rqc);
    Expr term = r->second;
    OrderedHandle<CellFilterStub> reg = rqc.domain();
    OrderedHandle<QuadratureFamilyStub> quad = rqc.quad();

    regionQuadComboExprs_.put(rqc, term);

    /* prepare calculation of both stiffness matrix and load vector */
    if (compTypes_.contains(MatrixAndVector))
    {
      SUNDANCE_MSG2(rqcVerb, tab2 << "preparing matrix/vector calculation");
      Tabs tab3; 
      EvalContext context(rqc, makeSet(1,2), contextID[0]);
      context.setSetupVerbosity(symbVerb);
      context.setEvalSetupVerbosity(evalSetupVerb);
      DerivSet nonzeros;
      
      if (isVariationalProblem_)
      {
        nonzeros = SymbPreprocessor
          ::setupVariations(term, 
            toList(vars), 
            toList(varLinearizationPts),
            toList(unks), 
            toList(unkLinearizationPts),
            unkParams, 
            unkParamEvalPts,
            toList(fixedFields), 
            toList(fixedFieldValues),
            fixedParams, 
            fixedParamEvalPts,
            context, 
            MatrixAndVector);
      }
      else
      {
        nonzeros = SymbPreprocessor
          ::setupFwdProblem(term, toList(vars), 
            toList(unks), 
            toList(unkLinearizationPts),
            unkParams, 
            unkParamEvalPts,
            fixedParams, 
            fixedParamEvalPts,
            toList(fixedFields), 
            toList(fixedFieldValues),
            context, 
            MatrixAndVector);
      }
      SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
      if (nonzeros.size()==0) 
      {
        rqcToSkip_[MatrixAndVector].put(rqc);
      }
      else
      {
        addToVarUnkPairs(rqc.domain(), varFuncSet, unkFuncSet,
          nonzeros, false, rqcVerb);
        rqcToContext_[MatrixAndVector].put(rqc, context);
        regionQuadComboNonzeroDerivs_[MatrixAndVector].put(rqc, 
          nonzeros);
      }
    }



    /* prepare calculation of load vector only */
    if (compTypes_.contains(VectorOnly))
    {
      SUNDANCE_MSG2(rqcVerb, tab2 << "preparing vector-only calculation");
      Tabs tab3; 
      EvalContext context(rqc, makeSet(1), contextID[1]);
      context.setSetupVerbosity(symbVerb);
      context.setEvalSetupVerbosity(evalSetupVerb);
      DerivSet nonzeros;
      if (isVariationalProblem_)
      {
        nonzeros = SymbPreprocessor
          ::setupVariations(term, toList(vars), 
            toList(varLinearizationPts),
            toList(unks), 
            toList(unkLinearizationPts),
            unkParams, 
            unkParamEvalPts,
            toList(fixedFields),
            toList(fixedFieldValues),
            fixedParams, 
            fixedParamEvalPts,
            context,
            VectorOnly);
      }
      else
      {
        nonzeros = SymbPreprocessor
          ::setupFwdProblem(term, toList(vars), 
            toList(unks), 
            toList(unkLinearizationPts),
            unkParams, 
            unkParamEvalPts,
            fixedParams, 
            fixedParamEvalPts,
            toList(fixedFields), 
            toList(fixedFieldValues),
            context,
            VectorOnly);
      }
      SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
      if (nonzeros.size()==0) 
      {
        rqcToSkip_[VectorOnly].put(rqc);
      }
      else
      {
        rqcToContext_[VectorOnly].put(rqc, context);
        regionQuadComboNonzeroDerivs_[VectorOnly].put(rqc, nonzeros);
      }
    }


    /* prepare calculation of sensitivities */
    if (compTypes_.contains(Sensitivities))
    {
      SUNDANCE_MSG2(rqcVerb, tab2 << "preparing sensitivity calculation");
      Tabs tab3;
      EvalContext context(rqc, makeSet(2), contextID[4]);
      context.setSetupVerbosity(symbVerb);
      context.setEvalSetupVerbosity(evalSetupVerb);
      DerivSet nonzeros;
      nonzeros = SymbPreprocessor
        ::setupSensitivities(term, toList(vars), 
          toList(unks), 
          toList(unkLinearizationPts),
          unkParams, 
          unkParamEvalPts,
          fixedParams, 
          fixedParamEvalPts,
          toList(fixedFields), 
          toList(fixedFieldValues),
          context,
          Sensitivities);
      SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
      if (nonzeros.size()==0) 
      {
        rqcToSkip_[Sensitivities].put(rqc);
      }
      else
      {
        rqcToContext_[Sensitivities].put(rqc, context);
        regionQuadComboNonzeroDerivs_[Sensitivities].put(rqc, nonzeros);
      }
    }



    /* prepare calculation of functional value only */
    if (compTypes_.contains(FunctionalOnly))
    {
      SUNDANCE_MSG2(rqcVerb, tab2 << "preparing functional calculation");
      Tabs tab3;

      EvalContext context(rqc, makeSet(0), contextID[2]);
      context.setSetupVerbosity(symbVerb);
      context.setEvalSetupVerbosity(evalSetupVerb);
      DerivSet nonzeros;
      Expr fields;
      Expr fieldValues;
      if (fixedFields.size() > 0)
      {
        if (vars.size() > 0)
        {
          fields = List(toList(fixedFields), toList(vars));
          fields = fields.flatten();
        }
        else
        {
          fields = toList(fixedFields);
        }
      }
      else
      {
        if (vars.size() > 0)
        {
          fields = toList(vars);
        }
      }
      if (fixedFieldValues.size() > 0)
      {
        if (varLinearizationPts.size() > 0)
        {
          fieldValues = List(toList(fixedFieldValues), 
            toList(varLinearizationPts));
          fieldValues = fieldValues.flatten();
        }
        else
        {
          fieldValues = toList(fixedFieldValues);
        }
      }
      else
      {
        if (varLinearizationPts.size() > 0)
        {
          fieldValues = toList(varLinearizationPts);
        }
      }

      nonzeros = SymbPreprocessor
        ::setupFunctional(term, 
          fixedParams, 
          fixedParamEvalPts,
          fields,
          fieldValues,
          context,
          FunctionalOnly);
      SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);

      if (nonzeros.size()==0) 
      {
        rqcToSkip_[FunctionalOnly].put(rqc);
      }
      else
      {
        rqcToContext_[FunctionalOnly].put(rqc, context);
        regionQuadComboNonzeroDerivs_[FunctionalOnly].put(rqc, nonzeros);
      }
    }
    /* prepare calculation of functional value and gradient */
    if (compTypes_.contains(FunctionalAndGradient))
    {
      SUNDANCE_MSG2(rqcVerb, tab2 << "preparing functional/gradient calculation");
      Tabs tab3;
      EvalContext context(rqc, makeSet(0,1), contextID[3]);
      context.setSetupVerbosity(symbVerb);
      context.setEvalSetupVerbosity(evalSetupVerb);
      DerivSet nonzeros;
      nonzeros = SymbPreprocessor
        ::setupGradient(term, 
          toList(vars), toList(varLinearizationPts),
          fixedParams, 
          fixedParamEvalPts,
          toList(fixedFields), toList(fixedFieldValues),
          context,
          FunctionalAndGradient);

      SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);

      if (nonzeros.size()==0) 
      {
        rqcToSkip_[FunctionalAndGradient].put(rqc);
      }
      else
      {
        rqcToContext_[FunctionalAndGradient].put(rqc, context);
        regionQuadComboNonzeroDerivs_[FunctionalAndGradient].put(rqc, nonzeros);
      }
    }
  }
  
  /* now do the BCs */
  if (hasBCs)
  {
    /* functions found in the BCs both in the overall lists and 
     * also in the bc-specific lists */
    SUNDANCE_MSG1(verb, tab1 << "processing essential BC terms");
    for (Sundance::Map<RegionQuadCombo, Expr>::const_iterator 
           r=bcSum->rqcToExprMap().begin(); 
         r!=bcSum->rqcToExprMap().end(); r++)
    {
      Tabs tab15;
      RegionQuadCombo rqc = r->first;
      int rqcVerb = verb;
      int symbVerb = 0;
      int evalSetupVerb = 0;
      if (rqc.watch().isActive()) 
      {
        symbVerb = rqc.watch().param("symbolic preprocessing");
        rqcVerb=rqc.watch().param("equation set setup");
        evalSetupVerb=rqc.watch().param("evaluator setup");
      }
      SUNDANCE_MSG1(verb, tab15 << "processing BC RQC = " << rqc);

      rqcBCSet.put(rqc);
      Expr term = r->second;
      OrderedHandle<CellFilterStub> reg = rqc.domain();
      OrderedHandle<QuadratureFamilyStub> quad = rqc.quad();

      bcRegionQuadComboExprs_.put(rqc, term); 


              
      /* prepare calculation of both stiffness matrix and load vector */
      if (compTypes_.contains(MatrixAndVector))
      {
        Tabs tab3;
        SUNDANCE_MSG2(rqcVerb, tab3 << "preparing matrix/vector calculation");
        EvalContext context(rqc, makeSet(1,2), contextID[0]);
        context.setSetupVerbosity(symbVerb);
        context.setEvalSetupVerbosity(evalSetupVerb);
        DerivSet nonzeros;
              
        if (isVariationalProblem_)
        {
          nonzeros = SymbPreprocessor
            ::setupVariations(term, 
              toList(vars), 
              toList(varLinearizationPts),
              toList(unks), 
              toList(unkLinearizationPts),
              unkParams, 
              unkParamEvalPts,
              toList(fixedFields), 
              toList(fixedFieldValues),
              fixedParams, 
              fixedParamEvalPts,
              context,
              MatrixAndVector);
        }
        else
        {
          nonzeros = SymbPreprocessor
            ::setupFwdProblem(term, toList(vars), toList(unks), 
              toList(unkLinearizationPts),
              unkParams, 
              unkParamEvalPts,
              fixedParams, 
              fixedParamEvalPts,
              toList(fixedFields), 
              toList(fixedFieldValues),
              context,
              MatrixAndVector);
        }
        SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
        if (nonzeros.size()==0) 
        {
          bcRqcToSkip_[MatrixAndVector].put(rqc);
        }
        else
        {
          addToVarUnkPairs(rqc.domain(), varFuncSet, unkFuncSet,
            nonzeros, true, rqcVerb);
          bcRqcToContext_[MatrixAndVector].put(rqc, context);
          bcRegionQuadComboNonzeroDerivs_[MatrixAndVector].put(rqc, 
            nonzeros);
        }
      }





      /* prepare calculation of load vector only */
      if (compTypes_.contains(VectorOnly))
      {
        Tabs tab3;
        SUNDANCE_MSG2(rqcVerb, tab3 << "preparing vector-only calculation");
        EvalContext context(rqc, makeSet(1), contextID[1]);
        context.setSetupVerbosity(symbVerb);
        context.setEvalSetupVerbosity(evalSetupVerb);
        DerivSet nonzeros;
        if (isVariationalProblem_)
        {
          nonzeros = SymbPreprocessor
            ::setupVariations(term, toList(vars), 
              toList(varLinearizationPts),
              toList(unks), 
              toList(unkLinearizationPts),
              unkParams, 
              unkParamEvalPts,
              toList(fixedFields), 
              toList(fixedFieldValues),
              fixedParams, 
              fixedParamEvalPts,
              context,
              VectorOnly);
        }
        else
        {
          nonzeros = SymbPreprocessor
            ::setupFwdProblem(term, toList(vars), toList(unks), 
              toList(unkLinearizationPts),
              unkParams, 
              unkParamEvalPts,
              fixedParams, 
              fixedParamEvalPts,
              toList(fixedFields), 
              toList(fixedFieldValues),
              context,
              VectorOnly);
        }
        SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
        if (nonzeros.size()==0) 
        {
          bcRqcToSkip_[VectorOnly].put(rqc);
        }
        else
        {
          bcRqcToContext_[VectorOnly].put(rqc, context);
          bcRegionQuadComboNonzeroDerivs_[VectorOnly].put(rqc, nonzeros);
        }
      }








      /* prepare calculation of sensitivities */
      if (compTypes_.contains(Sensitivities))
      {
        Tabs tab3;
        SUNDANCE_MSG2(rqcVerb, tab3 << "preparing sensitivity calculation");
        EvalContext context(rqc, makeSet(2), contextID[4]);
        context.setSetupVerbosity(symbVerb);
        context.setEvalSetupVerbosity(evalSetupVerb);
        DerivSet nonzeros;
        nonzeros = SymbPreprocessor
          ::setupSensitivities(term, toList(vars), toList(unks), 
            toList(unkLinearizationPts),
            unkParams, 
            unkParamEvalPts,
            fixedParams, 
            fixedParamEvalPts,
            toList(fixedFields), 
            toList(fixedFieldValues),
            context,
            Sensitivities);

        SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
        if (nonzeros.size()==0) 
        {
          bcRqcToSkip_[Sensitivities].put(rqc);
        }
        else
        {
          bcRqcToContext_[Sensitivities].put(rqc, context);
          bcRegionQuadComboNonzeroDerivs_[Sensitivities].put(rqc, nonzeros);
        }
      }







      /* prepare calculation of functional value only */
      if (compTypes_.contains(FunctionalOnly))
      {
        Tabs tab3;
        SUNDANCE_MSG2(rqcVerb, tab3 << "preparing functional-only calculation");
        EvalContext context(rqc, makeSet(0), contextID[2]);
        context.setSetupVerbosity(symbVerb);
        context.setEvalSetupVerbosity(evalSetupVerb);
        DerivSet nonzeros;
        Expr fields;
        Expr fieldValues;
        if (fixedFields.size() > 0)
        {
          if (vars.size() > 0)
          {
            fields = List(toList(fixedFields), toList(vars));
            fields = fields.flatten();
          }
          else
          {
            fields = toList(fixedFields);
          }
        }
        else
        {
          if (vars.size() > 0)
          {
            fields = toList(vars);
          }
        }
        if (fixedFieldValues.size() > 0)
        {
          if (varLinearizationPts.size() > 0)
          {
            fieldValues = List(toList(fixedFieldValues), 
              toList(varLinearizationPts));
            fieldValues = fieldValues.flatten();
          }
          else
          {
            fieldValues = toList(fixedFieldValues);
          }
        }
        else
        {
          if (varLinearizationPts.size() > 0)
          {
            fieldValues = toList(varLinearizationPts);
          }
        }
        nonzeros = SymbPreprocessor
          ::setupFunctional(term, 
            fixedParams, 
            fixedParamEvalPts,
            fields, fieldValues,
            context,
            FunctionalOnly);

        SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);                  
        if (nonzeros.size()==0) 
        {
          bcRqcToSkip_[FunctionalOnly].put(rqc);
        }
        else
        {
          bcRqcToContext_[FunctionalOnly].put(rqc, context);
          bcRegionQuadComboNonzeroDerivs_[FunctionalOnly].put(rqc, nonzeros);
        }
      }




      /* prepare calculation of functional value and gradient */
      if (compTypes_.contains(FunctionalAndGradient))
      {
        Tabs tab3;
        SUNDANCE_MSG2(rqcVerb, tab3 << "preparing functional and gradient calculation");
        EvalContext context(rqc, makeSet(0,1), contextID[3]);
        context.setSetupVerbosity(symbVerb);
        context.setEvalSetupVerbosity(evalSetupVerb);
        DerivSet nonzeros;
        nonzeros = SymbPreprocessor
          ::setupGradient(term, 
            toList(vars), 
            toList(varLinearizationPts),
            fixedParams, 
            fixedParamEvalPts,
            toList(fixedFields), 
            toList(fixedFieldValues),
            context,
            FunctionalAndGradient);

        SUNDANCE_MSG2(rqcVerb, tab3 << "nonzeros are " << nonzeros);
        if (nonzeros.size()==0) 
        {
          bcRqcToSkip_[FunctionalAndGradient].put(rqc);
        }
        else
        {
          bcRqcToContext_[FunctionalAndGradient].put(rqc, context);
          bcRegionQuadComboNonzeroDerivs_[FunctionalAndGradient].put(rqc, nonzeros);
        }
      }
    }
  }

  


  /* convert sets to arrays */
  regionQuadCombos_ = rqcSet.elements();
  bcRegionQuadCombos_ = rqcBCSet.elements();

  


}


void EquationSet
::addToVarUnkPairs(const OrderedHandle<CellFilterStub>& domain,
  const Set<int>& vars,
  const Set<int>& unks,
  const DerivSet& nonzeros, 
  bool isBC,
  int verb)
{
  Tabs tab;
  SUNDANCE_MSG2(verb, tab << "finding var-unk pairs "
    "for domain " << domain);
  SUNDANCE_MSG2(verb, tab << "isBC=" << isBC);
  
  RCP<Set<OrderedPair<int, int> > > funcPairs;
  Map<OrderedHandle<CellFilterStub>, RCP<Set<OrderedPair<int, int> > > >* theMap;

  if (isBC) 
  {
    theMap = &(bcVarUnkPairsOnRegions_);
  }
  else 
  {
    theMap = &(varUnkPairsOnRegions_);
  } 

  if (theMap->containsKey(domain))
  {
    funcPairs = theMap->get(domain);
  }
  else
  {
    funcPairs = rcp(new Set<OrderedPair<int, int> >());
    theMap->put(domain, funcPairs);
  }

  for (DerivSet::const_iterator i=nonzeros.begin(); i!=nonzeros.end(); i++)
  {
    Tabs tab1;
    const MultipleDeriv& md = *i;
    if (md.order() != 2) continue;
      
    Array<Deriv> f;
    for (MultipleDeriv::const_iterator j=md.begin(); j != md.end(); j++)
    {
      const Deriv& d = *j;
      TEUCHOS_TEST_FOR_EXCEPTION(!d.isFunctionalDeriv(), 
        std::logic_error, "non-functional deriv "
        << d << " detected in EquationSet::"
        "addToVarUnkPairs()");
      f.append(d);
    }

    SUNDANCE_MSG2(verb, tab1 << "f1=" << f[0].dofID()
      << ", f2=" << f[1].dofID() << ", vars=" << vars 
      << ", unks=" << unks);
    
    bool gotIt=false;
    if (unks.contains(f[0].dofID())
      && vars.contains(f[1].dofID()))
    {
      int unkID = f[0].dofID();
      int varID = f[1].dofID();
      funcPairs->put(OrderedPair<int, int>(varID, unkID));
      gotIt=true;
    }
    if (unks.contains(f[1].dofID())
      && vars.contains(f[0].dofID()))
    {
      int unkID = f[1].dofID();
      int varID = f[0].dofID();
      funcPairs->put(OrderedPair<int, int>(varID, unkID));
      gotIt=true;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!gotIt, std::logic_error,
      "no valid (var,unk) pair could be extracted from "
      "derivative " << md);
  }

  SUNDANCE_MSG2(verb, tab << "found " << *funcPairs);
  
}

bool EquationSet::hasActiveWatchFlag() const 
{
  for (int i=0; i<regionQuadCombos().size(); i++)
  {
    if (regionQuadCombos()[i].watch().isActive()) return true;
  }
  for (int i=0; i<bcRegionQuadCombos().size(); i++)
  {
    if (bcRegionQuadCombos()[i].watch().isActive()) return true;
  }
  return false;
}


int EquationSet::maxWatchFlagSetting(const std::string& param) const 
{
  int rtn = 0;
  if (!hasActiveWatchFlag()) return 0;
  for (int i=0; i<regionQuadCombos().size(); i++)
  {
    int v = regionQuadCombos()[i].watch().param(param);
    if (v > rtn) rtn = v;
  }
  for (int i=0; i<bcRegionQuadCombos().size(); i++)
  {
    int v = bcRegionQuadCombos()[i].watch().param(param);
    if (v > rtn) rtn = v;
  }
  return rtn;
}

Array<Expr> EquationSet::flattenSpectral(const Array<Expr>& expr) const
{
  Array<Expr> rtn(expr.size());
  for (int i=0; i<expr.size(); i++)
  {
    const Expr& e = expr[i];
    rtn[i] = flattenSpectral(e);
  }
  return rtn;
}

Expr EquationSet::flattenSpectral(const Expr& expr) const
{
  Array<Expr> rtn(expr.size());
  for (int i=0; i<expr.size(); i++)
  {
    if (expr[i].size() == 1)
    {
      const SpectralExpr* se 
        = dynamic_cast<const SpectralExpr*>(expr[i][0].ptr().get());
      if (se != 0)
      {
        int nt = se->getSpectralBasis().nterms();
        Array<Expr> e(nt);
        for (int j=0; j<nt; j++)
        {
          e[j] = se->getCoeff(j);
        }
        rtn[i] = new ListExpr(e);
      }
      else
      {
        rtn[i] = expr[i];
      }
    }
    else
    {
      rtn[i] = flattenSpectral(expr[i]);
    }
  }
  Expr r = new ListExpr(rtn);
  return r.flatten();
                  
}

const RCP<Set<OrderedPair<int, int> > >& EquationSet::
bcVarUnkPairs(const OrderedHandle<CellFilterStub>& domain) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!bcVarUnkPairsOnRegions_.containsKey(domain),
    std::logic_error,
    "equation set does not have a var-unk pair list for "
    "bc region " << domain);
  const RCP<Set<OrderedPair<int, int> > >& rtn 
    = bcVarUnkPairsOnRegions_.get(domain);

  TEUCHOS_TEST_FOR_EXCEPTION(rtn.get()==0, std::logic_error, 
    "null var-unk pair list for BC region " << domain);
  return rtn;
}

bool EquationSet::isBCRegion(int d) const
{
  return fsr_->isBCRegion(d);
}


EvalContext EquationSet::rqcToContext(ComputationType compType, 
  const RegionQuadCombo& r) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!rqcToContext_.containsKey(compType),
    std::logic_error,
    "EquationSet::rqcToContext() did not find key " 
    << compType);
  TEUCHOS_TEST_FOR_EXCEPTION(!rqcToContext_.get(compType).containsKey(r),
    std::logic_error, 
    "EquationSet::rqcToContext(" << compType 
    << ") did not find expected key " 
    << r);
    
  return rqcToContext_.get(compType).get(r);
}

EvalContext EquationSet::bcRqcToContext(ComputationType compType, 
  const RegionQuadCombo& r) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!bcRqcToContext_.containsKey(compType),
    std::logic_error,
    "EquationSet::bcRqcToContext() did not find key " 
    << compType);
  TEUCHOS_TEST_FOR_EXCEPTION(!bcRqcToContext_.get(compType).containsKey(r),
    std::logic_error, 
    "EquationSet::bcRqcToContext(" << compType 
    << ") did not find expected key " 
    << r);
  return bcRqcToContext_.get(compType).get(r);
}


bool EquationSet::skipRqc(ComputationType compType, 
  const RegionQuadCombo& r) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!rqcToSkip_.containsKey(compType),
    std::logic_error,
    "EquationSet::skipRqc() did not find expected key " 
    << compType);
    
  return rqcToSkip_.get(compType).contains(r);
}

bool EquationSet::skipBCRqc(ComputationType compType, 
  const RegionQuadCombo& r) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!bcRqcToSkip_.containsKey(compType),
    std::logic_error,
    "EquationSet::skipBCRqc() did not find expected key " 
    << compType);
    
  return bcRqcToSkip_.get(compType).contains(r);
}

const DerivSet& EquationSet::nonzeroFunctionalDerivs(ComputationType compType,
  const RegionQuadCombo& r) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!regionQuadComboNonzeroDerivs_.containsKey(compType),
    std::logic_error,
    "EquationSet:nonzeroFunctionalDerivs() did not find key " 
    << compType);
  return regionQuadComboNonzeroDerivs_.get(compType).get(r);
}

const DerivSet& EquationSet::nonzeroBCFunctionalDerivs(ComputationType compType,
  const RegionQuadCombo& r) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!bcRegionQuadComboNonzeroDerivs_.containsKey(compType),
    std::logic_error,
    "EquationSet:nonzeroBCFunctionalDerivs() did not find key " 
    << compType);
  return bcRegionQuadComboNonzeroDerivs_.get(compType).get(r);
}



int EquationSet::reducedVarID(int varID) const 
{
  return fsr_->reducedVarID(varID);
}

int EquationSet::reducedUnkID(int unkID) const 
{
  return fsr_->reducedUnkID(unkID);
}


int EquationSet::reducedUnkParamID(int unkParamID) const 
{
  return fsr_->reducedUnkParamID(unkParamID);
}

int EquationSet::reducedFixedParamID(int paramID) const 
{
  return fsr_->reducedFixedParamID(paramID);
}


Expr EquationSet::toList(const Array<Expr>& e)
{
  return new ListExpr(e);
}


int EquationSet::blockForVarID(int varID) const 
{
  return fsr_->blockForVarID(varID);
}

int EquationSet::blockForUnkID(int unkID) const 
{
  return fsr_->blockForUnkID(unkID);
}

const Set<OrderedHandle<CellFilterStub> >&  EquationSet::regionsForTestFunc(int testID) const
{
  return fsr_->regionsForTestFunc(testID);
}

const Set<OrderedHandle<CellFilterStub> >&  EquationSet::regionsForUnkFunc(int unkID) const
{
  return fsr_->regionsForUnkFunc(unkID);
}

int EquationSet::indexForRegion(const OrderedHandle<CellFilterStub>& region) const
{
  return fsr_->indexForRegion(region);
}

int EquationSet::numRegions() const {return fsr_->numRegions();}

const RCP<CellFilterStub>& EquationSet::region(int d) const 
{return fsr_->region(d);}



/* Returns the number of variational function blocks */
int EquationSet::numVarBlocks() const 
{return fsr_->numVarBlocks();}

/* Returns the number of unknown function blocks */
int EquationSet::numUnkBlocks() const 
{return fsr_->numUnkBlocks();}

/* Returns the number of unknown parameters */
int EquationSet::numUnkParams() const 
{return fsr_->numUnkParams();}

/* Returns the number of fixed parameters */
int EquationSet::numFixedParams() const 
{return fsr_->numFixedParams();}

/* Returns the number of variational functions in this block */
int EquationSet::numVars(int block) const 
{return fsr_->numVars(block);}

/* Returns the number of unk functions in this block */
int EquationSet::numUnks(int block) const 
{return fsr_->numUnks(block);}

/* Returns the number of variational function IDs in this block */
int EquationSet::numVarIDs(int block) const 
{return fsr_->numVarIDs(block);}

/* Returns the number of unk function IDs in this block */
int EquationSet::numUnkIDs(int block) const 
{return fsr_->numUnkIDs(block);}

/* Returns the i-th variational function in block b */
RCP<const CommonFuncDataStub>  EquationSet::varFuncData(int b, int i) const 
{return fsr_->varFuncData(b,i);}
    

/* Returns the i-th unknown function in block b */
RCP<const CommonFuncDataStub>  EquationSet::unkFuncData(int b, int i) const 
{return fsr_->unkFuncData(b,i);}

/* Returns the i-th unknown parameter */
const Expr& EquationSet::unkParam(int i) const 
{return fsr_->unkParam(i);}

/* Determine whether a given func ID is listed as a 
 * variational function in this equation set */
bool EquationSet::hasVarID(int fid) const 
{return fsr_->hasVarID(fid);}

/* Determine whether a given func ID is listed as a unk function 
 * in this equation set */
bool EquationSet::hasUnkID(int fid) const 
{return fsr_->hasUnkID(fid);}

/* Determine whether a given func ID is listed as a unk parameter 
 * in this equation set */
bool EquationSet::hasUnkParamID(int fid) const 
{return fsr_->hasUnkParamID(fid);}

/* Determine whether a given func ID is listed as a fixed parameter 
 * in this equation set */
bool EquationSet::hasFixedParamID(int fid) const 
{return fsr_->hasFixedParamID(fid);}

 
/* Returns the variational functions that appear explicitly
 * on the d-th region */
const Set<int>& EquationSet::varsOnRegion(int d) const 
{return fsr_->varsOnRegion(d);}

/* Returns the unknown functions that appear explicitly on the
 * d-th region. */
const Set<int>& EquationSet::unksOnRegion(int d) const 
{return fsr_->unksOnRegion(d);}

/* Returns the variational functions that 
 * appear in BCs on the d-th region.
 * We can use this information to tag certain rows as BC rows */
const Set<int>& EquationSet::bcVarsOnRegion(int d) const 
{return fsr_->bcVarsOnRegion(d);}

/* Returns the unknown functions that appear in BCs on the d-th region.
 * We can use this information to tag certain columns as BC
 * columns in the event we're doing symmetrized BC application */
const Set<int>& EquationSet::bcUnksOnRegion(int d) const 
{return fsr_->bcUnksOnRegion(d);}

/* Returns the reduced variational functions that appear explicitly
 * on the d-th region */
const Array<Set<int> >& EquationSet::reducedVarsOnRegion(const OrderedHandle<CellFilterStub>& r) const 
{return fsr_->reducedVarsOnRegion(r);}

/* Returns the reduced unknown functions that appear explicitly on the
 * d-th region. */
const Array<Set<int> >& EquationSet::reducedUnksOnRegion(const OrderedHandle<CellFilterStub>& r) const 
{return fsr_->reducedVarsOnRegion(r);}


/** get the unreduced funcID for a variational function
 * as specified by a reduced ID and block index */
int EquationSet::unreducedVarID(int block, int reducedVarID) const 
{return fsr_->unreducedVarID(block, reducedVarID);}

/** get the unreduced funcID for an unknown function
 * as specified by a reduced ID and block index */
int EquationSet::unreducedUnkID(int block, int reducedUnkID) const 
{return fsr_->unreducedUnkID(block, reducedUnkID);}


/** get the unreduced funcID for an unknown parameter
 * as specified by a reduced ID */
int EquationSet::unreducedUnkParamID(int reducedUnkParamID) const 
{return fsr_->unreducedUnkParamID(reducedUnkParamID);}


/** get the unreduced funcID for a fixed parameter
 * as specified by a reduced ID */
int EquationSet::unreducedFixedParamID(int reducedParamID) const 
{return fsr_->unreducedFixedParamID(reducedParamID);}

