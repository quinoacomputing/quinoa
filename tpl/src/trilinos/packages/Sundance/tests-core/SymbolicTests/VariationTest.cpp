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


#include "SundanceExpr.hpp"
#include "SundanceStdMathOps.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceProductTransformation.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceParameter.hpp"
#include "SundanceUnknownParameter.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceStringEvalMediator.hpp"

using Sundance::List;
using namespace Sundance;
using namespace Teuchos;
using std::endl;
using std::cerr;
using std::exception;



static Time& totalTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}

static Time& doitTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("doit"); 
  return *rtn;
}



void doVariations(const Expr& e, 
  const Expr& vars,
  const Expr& varEvalPt,
  const Expr& unks,
  const Expr& unkEvalPt, 
  const Expr& unkParams,
  const Expr& unkParamEvalPts,
  const Expr& fixed,
  const Expr& fixedEvalPt, 
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts, 
  const EvalContext& region)
{
  TimeMonitor t0(doitTimer());
  EvalManager mgr;
  mgr.setRegion(region);

  static RCP<AbstractEvalMediator> mediator 
    = rcp(new StringEvalMediator());

  mgr.setMediator(mediator);

  const EvaluatableExpr* ev 
    = dynamic_cast<const EvaluatableExpr*>(e[0].ptr().get());

  DerivSet d = SymbPreprocessor::setupVariations(e[0], 
    vars,
    varEvalPt,
    unks,
    unkEvalPt,
    unkParams,
    unkParamEvalPts,
    fixed,
    fixedEvalPt,
    fixedParams,
    fixedParamEvalPts,
    region,
    MatrixAndVector);

  Tabs tab;

  Array<double> constantResults;
  Array<RCP<EvalVector> > vectorResults;

  ev->evaluate(mgr, constantResults, vectorResults);

  ev->sparsitySuperset(region)->print(cerr, vectorResults, constantResults);
}



void doGradient(const Expr& e, 
  const Expr& vars,
  const Expr& varEvalPt,
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,
  const Expr& fixed,
  const Expr& fixedEvalPts, 
  const EvalContext& region)
{
  TimeMonitor t0(doitTimer());
  EvalManager mgr;
  mgr.setRegion(region);

  static RCP<AbstractEvalMediator> mediator 
    = rcp(new StringEvalMediator());

  mgr.setMediator(mediator);

  const EvaluatableExpr* ev 
    = dynamic_cast<const EvaluatableExpr*>(e[0].ptr().get());

  DerivSet d = SymbPreprocessor::setupGradient(e[0], 
    vars,
    varEvalPt,
    fixedParams,
    fixedParamEvalPts,
    fixed,
    fixedEvalPts,
    region,
    FunctionalAndGradient);

  Tabs tab;
  //  std::cerr << tab << *ev->sparsitySuperset(region) << std::endl;
  //  ev->showSparsity(cerr, region);

  // RCP<EvalVectorArray> results;

  Array<double> constantResults;
  Array<RCP<EvalVector> > vectorResults;

  ev->evaluate(mgr, constantResults, vectorResults);

  ev->sparsitySuperset(region)->print(cerr, vectorResults, constantResults);

  
  // results->print(cerr, ev->sparsitySuperset(region).get());
}



void doFunctional(const Expr& e, 
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,
  const Expr& fixed,
  const Expr& fixedEvalPt, 
  const EvalContext& region)
{
  TimeMonitor t0(doitTimer());
  EvalManager mgr;
  mgr.setRegion(region);

  static RCP<AbstractEvalMediator> mediator 
    = rcp(new StringEvalMediator());

  mgr.setMediator(mediator);

  const EvaluatableExpr* ev 
    = dynamic_cast<const EvaluatableExpr*>(e[0].ptr().get());

  DerivSet d = SymbPreprocessor::setupFunctional(e[0], 
    fixedParams,
    fixedParamEvalPts,
    fixed,
    fixedEvalPt,
    region,
    FunctionalOnly);

  Tabs tab;
  //  std::cerr << tab << *ev->sparsitySuperset(region) << std::endl;
  //  ev->showSparsity(cerr, region);

  // RCP<EvalVectorArray> results;

  Array<double> constantResults;
  Array<RCP<EvalVector> > vectorResults;

  ev->evaluate(mgr, constantResults, vectorResults);

  ev->sparsitySuperset(region)->print(cerr, vectorResults, constantResults);

  
  // results->print(cerr, ev->sparsitySuperset(region).get());
}



void testVariations(const Expr& e,  
  const Expr& vars,
  const Expr& varEvalPt,
  const Expr& unks,
  const Expr& unkEvalPt, 
  const Expr& unkParams,
  const Expr& unkParamEvalPts,
  const Expr& fixed,
  const Expr& fixedEvalPt, 
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,  
  const EvalContext& region)
{
  std::cerr << std::endl 
       << "------------------------------------------------------------- " << std::endl;
  std::cerr  << "-------- testing " << e.toString() << " -------- " << std::endl;
  std::cerr << std::endl 
       << "------------------------------------------------------------- " << std::endl;

  try
  {
    doVariations(e, vars, varEvalPt, 
      unks, unkEvalPt, 
      unkParams,
      unkParamEvalPts,
      fixed, fixedEvalPt, 
      fixedParams,
      fixedParamEvalPts,
      region);
  }
  catch(std::exception& ex)
  {
    std::cerr << "EXCEPTION DETECTED!" << std::endl;
    std::cerr << ex.what() << std::endl;
    // std::cerr << "repeating with increased verbosity..." << std::endl;
    //       std::cerr << "-------- testing " << e.toString() << " -------- " << std::endl;
    //       Evaluator::verbosity() = 2;
    //       EvalVector::verbosity() = 2;
    //       EvaluatableExpr::verbosity() = 2;
    //       Expr::showAllParens() = true;
    //       doit(e, region);
    exit(1);
  }
}

void testGradient(const Expr& e,  
  const Expr& vars,
  const Expr& varEvalPt,
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,
  const Expr& fixed,
  const Expr& fixedEvalPt,  
  const EvalContext& region)
{
  std::cerr << std::endl 
       << "------------------------------------------------------------- " << std::endl;
  std::cerr  << "-------- testing " << e.toString() << " -------- " << std::endl;
  std::cerr << std::endl 
       << "------------------------------------------------------------- " << std::endl;

  try
  {
    doGradient(e, vars, varEvalPt,  
      fixedParams,
      fixedParamEvalPts,
      fixed, 
      fixedEvalPt, 
      region);
  }
  catch(std::exception& ex)
  {
    std::cerr << "EXCEPTION DETECTED!" << std::endl;
    std::cerr << ex.what() << std::endl;
    // std::cerr << "repeating with increased verbosity..." << std::endl;
    //       std::cerr << "-------- testing " << e.toString() << " -------- " << std::endl;
    //       Evaluator::verbosity() = 2;
    //       EvalVector::verbosity() = 2;
    //       EvaluatableExpr::verbosity() = 2;
    //       Expr::showAllParens() = true;
    //       doit(e, region);
    exit(1);
  }
}


void testFunctional(const Expr& e, 
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts, 
  const Expr& fixed,
  const Expr& fixedEvalPt,  
  const EvalContext& region)
{
  std::cerr << std::endl 
       << "------------------------------------------------------------- " << std::endl;
  std::cerr  << "-------- testing " << e.toString() << " -------- " << std::endl;
  std::cerr << std::endl 
       << "------------------------------------------------------------- " << std::endl;

  try
  {
    doFunctional(e,   
      fixedParams,
      fixedParamEvalPts,fixed, fixedEvalPt, region);
  }
  catch(std::exception& ex)
  {
    std::cerr << "EXCEPTION DETECTED!" << std::endl;
    std::cerr << ex.what() << std::endl;
    // std::cerr << "repeating with increased verbosity..." << std::endl;
    //       std::cerr << "-------- testing " << e.toString() << " -------- " << std::endl;
    //       Evaluator::verbosity() = 2;
    //       EvalVector::verbosity() = 2;
    //       EvaluatableExpr::verbosity() = 2;
    //       Expr::showAllParens() = true;
    //       doit(e, region);
    exit(1);
  }
}


int main(int argc, char** argv)
{
  
  try
  {
    GlobalMPISession session(&argc, &argv);

    TimeMonitor t(totalTimer());

    verbosity<SymbolicTransformation>() = 0;
    verbosity<Evaluator>() = 0;
    verbosity<EvalVector>() = 0;
    verbosity<EvaluatableExpr>() = 0;
    Expr::showAllParens() = true;
    ProductTransformation::optimizeFunctionDiffOps() = false;

    EvalVector::shadowOps() = true;

    Expr dx = new Derivative(0);
    Expr dy = new Derivative(1);
    Expr grad = List(dx, dy);

    Expr u = new UnknownFunctionStub("u");
    Expr lambda_u = new UnknownFunctionStub("lambda_u");
    Expr T = new UnknownFunctionStub("T");
    Expr lambda_T = new UnknownFunctionStub("lambda_T");
    Expr alpha = new UnknownParameter("alpha");

    Expr u0 = new DiscreteFunctionStub("u0");
    Expr lambda_u0 = new DiscreteFunctionStub("lambda_u0");
    Expr T0 = new DiscreteFunctionStub("T0");
    Expr lambda_T0 = new DiscreteFunctionStub("lambda_T0");
    Expr zero = new ZeroExpr();
    Expr alpha0 = new Parameter(3.14, "alpha0");

    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);

    Expr empty;

    Array<Expr> tests;

//#define BLAHBLAH 1
#ifdef BLAHBLAH
    verbosity<Evaluator>() = 5;
    verbosity<SparsitySuperset>() = 5;
    verbosity<EvaluatableExpr>() = 5;
#endif

    tests.append( 0.5*(u-1.404)*(u-1.404)
      + 0.5*(grad*u)*(grad*u)
      + 0.5*alpha*alpha
      + (grad*lambda_u)*(grad*u)
      + lambda_u*alpha);


    std::cerr << std::endl << "============== u STATE EQUATIONS =================" << std::endl;

    for (int i=0; i<tests.length(); i++)
    {
      RegionQuadCombo rqc(rcp(new CellFilterStub()), 
        rcp(new QuadratureFamilyStub(1)));
      EvalContext context(rqc, makeSet(1,2), EvalContext::nextID());
      testVariations(tests[i], 
        List(lambda_u),
        List(zero),
        List(u),
        List(u0),
        empty,
        empty,
        empty,
        empty,
        alpha,
        alpha0,
        context);
    }

    std::cerr << std::endl << "=============== u ADJOINT EQUATIONS =================" << std::endl;
    for (int i=0; i<tests.length(); i++)
    {
      RegionQuadCombo rqc(rcp(new CellFilterStub()), 
        rcp(new QuadratureFamilyStub(1)));
      EvalContext context(rqc, makeSet(1,2), EvalContext::nextID());
      testVariations(tests[i], 
        List(u),
        List(u0),
        List(lambda_u),
        List(zero),
        empty,
        empty,
        empty,
        empty,
        alpha,
        alpha0,
        context);
    }


    std::cerr << std::endl << "================ REDUCED GRADIENT ====================" << std::endl;
    for (int i=0; i<tests.length(); i++)
    {
      RegionQuadCombo rqc(rcp(new CellFilterStub()), 
        rcp(new QuadratureFamilyStub(1)));
      EvalContext context(rqc, makeSet(0,1), EvalContext::nextID());
      testGradient(tests[i], 
        alpha, 
        alpha0,
        empty,
        empty,
        List(u, lambda_u),
        List(u0, zero),
        context);
    }

    std::cerr << std::endl << "=================== FUNCTIONAL ====================" << std::endl;
    for (int i=0; i<tests.length(); i++)
    {
      RegionQuadCombo rqc(rcp(new CellFilterStub()), 
        rcp(new QuadratureFamilyStub(1)));
      EvalContext context(rqc, makeSet(0), EvalContext::nextID());
      testFunctional(tests[i], 
        alpha,
        alpha0,
        List(u, lambda_u),
        List(u0, zero),
        context);
    }

    std::cerr << std::endl << "=================== ALL-AT-ONCE ====================" << std::endl;
    for (int i=0; i<tests.length(); i++)
    {
      RegionQuadCombo rqc(rcp(new CellFilterStub()), 
        rcp(new QuadratureFamilyStub(1)));
      EvalContext context(rqc, makeSet(1,2), EvalContext::nextID());
      testVariations(tests[i], 
        List(lambda_u, u, alpha), 
        List(zero, zero, zero), 
        List(lambda_u, u, alpha), 
        List(lambda_u0, u0, alpha0), 
        empty,
        empty,
        empty,
        empty,
        empty,
        empty,
        context);
      
    }
    TimeMonitor::summarize();
  }
	catch(std::exception& e)
  {
    Out::println(e.what());
  }


  
}
