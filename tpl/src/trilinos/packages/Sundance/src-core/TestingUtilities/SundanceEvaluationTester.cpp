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


#include "SundanceEvaluationTester.hpp"
#include "SundanceOut.hpp"
#include "SundanceExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "PlayaTabs.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceMultiIndex.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace SundanceTesting;

using namespace Teuchos;
using namespace std;



EvaluationTester::EvaluationTester(const Expr& e, int maxDiffOrder)
  : e_(e),
    rqc_(),
    context_(),
    mgr_(),
    mediator_(),
    tem_(),
    ev_(dynamic_cast<const EvaluatableExpr*>(e[0].ptr().get())),
    sparsity_(),
    unkIDToDiscreteIDMap_(),
    maxDiffOrder_(maxDiffOrder)
{
  setVerb(classVerbosity());
  Tabs tabs;

  SUNDANCE_VERB_LOW(tabs << "creating tester for expression " << e.toString());

  /* ------------------------------------------------------ */
  /* create a dummy cell filter and quadrature rule         */
  /* ------------------------------------------------------ */
  rqc_ = RegionQuadCombo(rcp(new CellFilterStub()), 
    rcp(new QuadratureFamilyStub(1)));


  /* ------------------------------------------------------ */
  /* find the unknown functions in the expression           */
  /* ------------------------------------------------------ */
  Array<Expr> unks;
  RCP<Set<int> > unkID = rcp(new Set<int>());

  SUNDANCE_VERB_LOW(tabs << "finding unknowns...");

  const UnknownFuncElement* ue 
    = dynamic_cast<const UnknownFuncElement*>(e[0].ptr().get());
  if (ue!=0)
  {
    unkID->put(ue->fid().dofID());
    unks.append(e[0]);
  }
  else
  {
    ev_->getUnknowns(*unkID, unks);
  }


 
  /* ------------------------------------------------------ */
  /* create a context object                                */
  /* ------------------------------------------------------ */ 
  Set<int> d;
  for (int i=0; i<=maxDiffOrder; i++) d.put(i);
  context_ = EvalContext(rqc_, d, EvalContext::nextID()); 
 


  /* ------------------------------------------------------ 
   * tell the evaluation manager which context we're using  
   * ------------------------------------------------------ */  
  mgr_.setRegion(context_);


  /* ------------------------------------------------------ *
   * We next create a "TestEvalMediator" that knows how to 
   * compute phony "discrete functions" defined in terms of
   * AD objects. This evaluation mediator is a stand-in for
   * a simulation framework with real discrete function.
   * 
   * In this loop we are creating phony discrete functions
   * and associating them with the unknown functions in the
   * problem. 
   * ------------------------------------------------------ */  
  Array<Expr> uDisc;
  Array<Expr> uEval;

  for (int i=0; i<unks.size(); i++)
  {
    Tabs tabs1;
    const UnknownFuncElement* fe 
      = dynamic_cast<const UnknownFuncElement*>(unks[i].ptr().get());
    TEUCHOS_TEST_FOR_EXCEPTION(fe==0, std::logic_error,
      "unk " << unks[i] << " is not an UnknownFunction");

    RCP<const UnknownFuncDataStub> data 
      = rcp_dynamic_cast<const UnknownFuncDataStub>(fe->commonData());

    RCP<const TestUnknownFuncData> tufd  
      = rcp_dynamic_cast<const TestUnknownFuncData>(data);

    TEUCHOS_TEST_FOR_EXCEPTION(tufd.get()==0, std::logic_error,
      "unk " << unks[i] << " is not a TestUnknownFunction");

    Expr discFunc = tufd->createDiscreteFunction(fe->name());

    const DiscreteFuncElement* df 
      = dynamic_cast<const DiscreteFuncElement*>(discFunc[0].ptr().get());
    TEUCHOS_TEST_FOR_EXCEPTION(df==0, std::logic_error,
      "df " << discFunc 
      << " is not a DiscreteFuncElement");

    SUNDANCE_VERB_LOW(tabs1 << i << " unk=" << unks[i] 
      << " disc=" << discFunc);

    uDisc.append(discFunc);
    unkIDToDiscreteIDMap_.put(fe->fid().dofID(), df->fid().dofID());
    if (tufd->coeffIsZero())
    {
      Expr Z = new ZeroExpr();
      uEval.append(Z);
    }
    else
    {
      uEval.append(discFunc);
    }
  }

  Expr unkList = new ListExpr(unks);
  Expr discList = new ListExpr(uDisc);
  Expr evalList = new ListExpr(uEval);

  SUNDANCE_VERB_LOW(tabs << "creating test eval mediator...");
  
  /* ------------------------------------------------------------- 
   * register the evaluator mediator with the evaluation manager
   * ------------------------------------------------------------- */  
  tem_ = new TestEvalMediator(discList);
  mediator_ = rcp(tem_);
  mgr_.setMediator(mediator_);

  SUNDANCE_VERB_LOW(tabs << "setting up evaluation...");
  Expr dummy;
  
  /* ------------------------------------------------------------- 
   * do preprocessing for this expression
   * - compute sparsity sets
   * - create evaluators 
   * ------------------------------------------------------------- */  

  if (maxDiffOrder_ > 0)
  {
    ComputationType compType=VectorOnly;    
    if (maxDiffOrder_==2) compType=MatrixAndVector;

    SymbPreprocessor::setupVariations(e[0], 
      unkList,
      evalList,
      unkList,
      evalList,
      dummy,
      dummy,
      dummy,
      dummy,
      dummy,
      dummy,
      context_,
      compType);
  }
  else
  {
    SymbPreprocessor::setupFunctional(e[0], 
      dummy,
      dummy,
      unkList,
      evalList,
      context_,
      FunctionalOnly);
  }

  Out::os() << tabs << "sparsity pattern: " << std::endl;
  sparsity_ = ev_->sparsitySuperset(context_);
  sparsity_->print(Out::os());

}



double EvaluationTester::evaluate() const 
{
  Tabs tabs;

  Array<double> constantResults;
  Array<RCP<EvalVector> > vectorResults;

  
  SUNDANCE_VERB_MEDIUM(tabs << std::endl << tabs << "evaluating...");

  tem_->setEvalPoint(ADField::evalPoint());

  
  /* ------------------------------------------------------------- 
   * Reset the number of calls to zero; this will flush all
   * cached evaluation vectors. This needs to be done before
   * any top-level evaluation call.
   * ------------------------------------------------------------- */ 
  
  ev_->evaluator(context_)->resetNumCalls();


  /* ------------------------------------------------------------- 
   * Do the evaluation!
   * Return constant results and vector results in different arrays
   * ------------------------------------------------------------- */ 

  ev_->evaluate(mgr_, constantResults, vectorResults);



  /* ------------------------------------------------------------- 
   * Print results if asked to
   * ------------------------------------------------------------- */ 

  if (verb() > 1)
  {
    ev_->sparsitySuperset(context_)->print(Out::os(), vectorResults,
      constantResults);
  }

  /* ------------------------------------------------------------- 
   * Sum results to compute functional value
   * ------------------------------------------------------------- */ 
  int vectorCount=0;
  int constantCount=0;
  double rtn = 0.0;

  for (int i=0; i<sparsity_->numDerivs(); i++)
  {
    const MultipleDeriv& md = sparsity_->deriv(i);
      
    Array<int> fieldIndex;
    Array<MultiIndex> mi;

    SUNDANCE_VERB_MEDIUM("md=" << md);

    for (MultipleDeriv::const_iterator 
           iter=md.begin(); iter != md.end(); iter++)
    {
      const Deriv& d = *iter;
      TEUCHOS_TEST_FOR_EXCEPTION(d.isCoordDeriv(), std::logic_error,
        "coordinate deriv found in TestEvalMediator::"
        "sumFunctionalChainRule");
      int uid = d.fid().dofID();
      SUNDANCE_VERB_EXTREME("deriv=" << d << " uid=" << uid);
      TEUCHOS_TEST_FOR_EXCEPTION(!unkIDToDiscreteIDMap_.containsKey(uid),
        std::logic_error,
        "uid " << uid << " not found in map " 
        <<unkIDToDiscreteIDMap_ );
      int fid = unkIDToDiscreteIDMap_.get(uid);
      int m = tem_->funcIdToFieldNumberMap().get(fid);
      fieldIndex.append(m);
      mi.append(d.opOnFunc().mi());
    }

    SUNDANCE_VERB_MEDIUM("field indices are " << fieldIndex
      << ", multiindices are=" <<mi.toString() );


    int resultIndex;
    if (sparsity_->state(i)==ConstantDeriv)
    {
      resultIndex = constantCount++;
    }
    else
    {
      resultIndex = vectorCount++;
    }
      
    double fieldDerivs = 1.0;
    for (int k=0; k<fieldIndex.size(); k++)
    {
      fieldDerivs *= tem_->evalDummyBasis(fieldIndex[k], mi[k]);
    }

    double coeff;
    if (sparsity_->state(i)==ConstantDeriv)
    {
      coeff = constantResults[resultIndex];
    }
    else
    {
      coeff = vectorResults[resultIndex]->start()[0];
    }
    SUNDANCE_VERB_HIGH("field deriv " << md.toString() 
      << " value=" << fieldDerivs << " coeff=" << coeff);

    if (fieldIndex.size() == 0)
    {
      SUNDANCE_VERB_HIGH("adding " << coeff * fieldDerivs << " to result");
      rtn += coeff * fieldDerivs;
    }
  }

  return rtn;
}





double EvaluationTester::evaluate(Array<double>& firstDerivs) const 
{
  Tabs tabs;

  firstDerivs.resize(tem_->numFields());
  for (int i=0; i<tem_->numFields(); i++) 
  {
    firstDerivs[i] = 0.0;
  }

  Array<double> constantResults;
  Array<RCP<EvalVector> > vectorResults;

  
  SUNDANCE_VERB_MEDIUM(tabs << std::endl << tabs << "evaluating...");

  tem_->setEvalPoint(ADField::evalPoint());
  ev_->evaluator(context_)->resetNumCalls();
  ev_->evaluate(mgr_, constantResults, vectorResults);

  if (verb() > 1)
  {
    ev_->sparsitySuperset(context_)->print(Out::os(), vectorResults,
      constantResults);
  }

  int vectorCount=0;
  int constantCount=0;
  double rtn = 0.0;

  for (int i=0; i<sparsity_->numDerivs(); i++)
  {
    const MultipleDeriv& md = sparsity_->deriv(i);
      
    Array<int> fieldIndex;
    Array<MultiIndex> mi;

    SUNDANCE_VERB_MEDIUM("md=" << md);

    for (MultipleDeriv::const_iterator 
           iter=md.begin(); iter != md.end(); iter++)
    {
      const Deriv& d = *iter;
      TEUCHOS_TEST_FOR_EXCEPTION(d.isCoordDeriv(), std::logic_error,
        "coordinate deriv found in TestEvalMediator::"
        "sumFunctionalChainRule");
      int uid = d.fid().dofID();
      SUNDANCE_VERB_EXTREME("deriv=" << d << " uid=" << uid);
      TEUCHOS_TEST_FOR_EXCEPTION(!unkIDToDiscreteIDMap_.containsKey(uid),
        std::logic_error,
        "uid " << uid << " not found in map " 
        <<unkIDToDiscreteIDMap_ );
      int fid = unkIDToDiscreteIDMap_.get(uid);
      int m = tem_->funcIdToFieldNumberMap().get(fid);
      fieldIndex.append(m);
      mi.append(d.opOnFunc().mi());
    }

    SUNDANCE_VERB_MEDIUM("field indices are " << fieldIndex
      << ", multiindices are=" <<mi.toString() );


    int resultIndex;
    if (sparsity_->state(i)==ConstantDeriv)
    {
      resultIndex = constantCount++;
    }
    else
    {
      resultIndex = vectorCount++;
    }
      
    double fieldDerivs = 1.0;
    for (int k=0; k<fieldIndex.size(); k++)
    {
      fieldDerivs *= tem_->evalDummyBasis(fieldIndex[k], mi[k]);
    }

    double coeff;
    if (sparsity_->state(i)==ConstantDeriv)
    {
      coeff = constantResults[resultIndex];
    }
    else
    {
      coeff = vectorResults[resultIndex]->start()[0];
    }
    SUNDANCE_VERB_HIGH("field deriv " << md.toString() 
      << " value=" << fieldDerivs << " coeff=" << coeff);

    if (fieldIndex.size() == 0)
    {
      SUNDANCE_VERB_HIGH("adding " << coeff * fieldDerivs << " to result");
      rtn += coeff * fieldDerivs;
    }
    else if (fieldIndex.size() == 1)
    {
      SUNDANCE_VERB_HIGH("adding " << coeff * fieldDerivs << " to first deriv");
      firstDerivs[fieldIndex[0]] += coeff * fieldDerivs;
    }
  }

  return rtn;
}








double EvaluationTester::evaluate(Array<double>& firstDerivs, 
  Array<Array<double> >& secondDerivs) const 
{
  Tabs tabs;

  firstDerivs.resize(tem_->numFields());
  secondDerivs.resize(tem_->numFields());
  for (int i=0; i<tem_->numFields(); i++) 
  {
    firstDerivs[i] = 0.0;
    secondDerivs[i].resize(tem_->numFields());
    for (int j=0; j<tem_->numFields(); j++)
    {
      secondDerivs[i][j] = 0.0;
    }
  }

  Array<double> constantResults;
  Array<RCP<EvalVector> > vectorResults;

  
  SUNDANCE_VERB_MEDIUM(tabs << std::endl << tabs << "evaluating...");

  tem_->setEvalPoint(ADField::evalPoint());
  ev_->evaluator(context_)->resetNumCalls();
  ev_->evaluate(mgr_, constantResults, vectorResults);

  if (verb() > 1)
  {
    ev_->sparsitySuperset(context_)->print(Out::os(), vectorResults,
      constantResults);
  }

  int vectorCount=0;
  int constantCount=0;
  double rtn = 0.0;

  for (int i=0; i<sparsity_->numDerivs(); i++)
  {
    const MultipleDeriv& md = sparsity_->deriv(i);
      
    Array<int> fieldIndex;
    Array<MultiIndex> mi;

    SUNDANCE_VERB_MEDIUM("md=" << md);

    for (MultipleDeriv::const_iterator 
           iter=md.begin(); iter != md.end(); iter++)
    {
      const Deriv& d = *iter;
      TEUCHOS_TEST_FOR_EXCEPTION(d.isCoordDeriv(), std::logic_error,
        "coordinate deriv found in TestEvalMediator::"
        "sumFunctionalChainRule");
      int uid = d.fid().dofID();
      SUNDANCE_VERB_EXTREME("deriv=" << d << " uid=" << uid);
      TEUCHOS_TEST_FOR_EXCEPTION(!unkIDToDiscreteIDMap_.containsKey(uid),
        std::logic_error,
        "uid " << uid << " not found in map " 
        <<unkIDToDiscreteIDMap_ );
      int fid = unkIDToDiscreteIDMap_.get(uid);
      int m = tem_->funcIdToFieldNumberMap().get(fid);
      fieldIndex.append(m);
      mi.append(d.opOnFunc().mi());
    }

    SUNDANCE_VERB_MEDIUM("field indices are " << fieldIndex
      << ", multiindices are=" <<mi.toString() );


    int resultIndex;
    if (sparsity_->state(i)==ConstantDeriv)
    {
      resultIndex = constantCount++;
    }
    else
    {
      resultIndex = vectorCount++;
    }
      
      
    double fieldDerivs = 1.0;
    for (int k=0; k<fieldIndex.size(); k++)
    {
      fieldDerivs *= tem_->evalDummyBasis(fieldIndex[k], mi[k]);
    }

    double coeff;
    if (sparsity_->state(i)==ConstantDeriv)
    {
      coeff = constantResults[resultIndex];
    }
    else
    {
      coeff = vectorResults[resultIndex]->start()[0];
    }
    SUNDANCE_VERB_HIGH("field deriv " << md.toString() 
      << " value=" << fieldDerivs << " coeff=" << coeff);

    if (fieldIndex.size() == 0)
    {
      SUNDANCE_VERB_HIGH("adding " << coeff * fieldDerivs << " to result");
      rtn += coeff * fieldDerivs;
    }
    else if (fieldIndex.size() == 1)
    {
      SUNDANCE_VERB_HIGH("adding " << coeff * fieldDerivs << " to first deriv");
      firstDerivs[fieldIndex[0]] += coeff * fieldDerivs;
    }
    else
    {
      int multiplicity = 1;
      if (fieldIndex[0] != fieldIndex[1] || !(mi[0] == mi[1]))
      {
        multiplicity = 2;
      }
      SUNDANCE_VERB_HIGH("adding " << coeff * fieldDerivs << " to second deriv with multiplicity " << multiplicity);
      secondDerivs[fieldIndex[0]][fieldIndex[1]] += multiplicity * coeff * fieldDerivs;
      //     if (fieldIndex[0] != fieldIndex[1])
//             {
//               secondDerivs[fieldIndex[1]][fieldIndex[0]] += multiplicity * coeff * fieldDerivs;
//             }

    }
  }

  Array<Array<double> > tmp = secondDerivs;
  /* symmetrize the second derivs */
  for (int i=0; i<tem_->numFields(); i++) 
  {
    for (int j=0; j<tem_->numFields(); j++)
    {
      secondDerivs[i][j] = (tmp[i][j] + tmp[j][i])/2.0;
    }
  }
  return rtn;
}


double EvaluationTester
::fdEvaluate(const double& step, const double& tol, const double& tol2,
  bool& isOK)
{
  Array<double> afdFirst(tem_->numFields());
  Array<Array<double> > afdSecond(tem_->numFields()); 

  Array<double> fdFirst(tem_->numFields());
  Array<Array<double> > fdSecond(tem_->numFields());

  Array<double> firstPlus(tem_->numFields());
  Array<double> firstMinus(tem_->numFields());
  Array<Array<double> > tmpSecond(tem_->numFields());
  
  for (int i=0; i<tem_->numFields(); i++) 
  {
    afdSecond[i].resize(tem_->numFields());
    fdSecond[i].resize(tem_->numFields());
    tmpSecond[i].resize(tem_->numFields());
  }

  Out::os() << setprecision(14);
  double f0 = evaluate(afdFirst, afdSecond);

  for (int i=0; i<tem_->numFields(); i++)
  {
    Tabs tab2;
    double A0 = tem_->fieldCoeff(i);
    tem_->setFieldCoeff(i, A0 - step);
    double fMinus = evaluate(firstMinus, tmpSecond);
    tem_->setFieldCoeff(i, A0 + step);
    double fPlus = evaluate(firstPlus, tmpSecond);
    fdFirst[i] = (fPlus - fMinus)/2.0/step;
    tem_->setFieldCoeff(i, A0);
    double error1 = fabs(fdFirst[i] - afdFirst[i])/(step + fabs(afdFirst[i]));
    Out::os() << tab2 << std::endl;
    SUNDANCE_VERB_MEDIUM(tab2 << "f(A_i+h)=" << fPlus << "      f(A_i-h)=" << fMinus);

    Out::os() << tab2 << "field " << tem_->fieldName(i) << " exact first deriv=" << afdFirst[i]
              << "    fd=" << fdFirst[i] << "    |exact - fd|=" 
              << error1 << std::endl;
    if (error1 > tol) 
    {
      isOK = false;
      Out::os() << tab2 << "first deriv wrt field "
                << i << " calculation FAILED" << std::endl;
    }
    Out::os() << tab2 << std::endl;

    if (maxDiffOrder_ < 2) continue;
    /* second deriv wrt this field by finite differences 
     * on the AFD first derivs*/
    fdSecond[i][i] = (firstPlus[i] - firstMinus[i])/2.0/step;
    double error2 = fabs(fdSecond[i][i] - afdSecond[i][i])/(step + fabs(afdSecond[i][i]));
    Out::os() << tab2 << "field " << tem_->fieldName(i) << " exact second deriv=" << afdSecond[i][i]
              << "    fd=" << fdSecond[i][i] << "    |exact - fd|=" 
              << error2 << std::endl;
    if (error2 > tol2) 
    {
      isOK = false;
      Out::os() << tab2 << "second deriv calculation wrt field "
                << tem_->fieldName(i) << " FAILED" << std::endl;
      Tabs tab3;
      Out::os() << tab3 << "f'(A_i+h) = " << firstPlus[i] << std::endl;
      Out::os() << tab3 << "f'(A_i-h) = " << firstMinus[i] << std::endl;
      Out::os() << tab3 << "step=" << step << std::endl;
    }
      
    /* mixed partials by finite differences on the AFD first derivs*/
    for (int j=0; j<i; j++)
    {
      //  /* -1,-1 node */
//           double B0 = tem_->fieldCoeff(j);
//           tem_->setFieldCoeff(i, A0 - step);
//           tem_->setFieldCoeff(j, B0 - step);
//           double fMM = evaluate(tmpFirst, tmpSecond);
//           /* -1,+1 node */
//           tem_->setFieldCoeff(i, A0 - step);
//           tem_->setFieldCoeff(j, B0 + step);
//           double fMP = evaluate(tmpFirst, tmpSecond);
//           /* +1,+1 node */
//           tem_->setFieldCoeff(i, A0 + step);
//           tem_->setFieldCoeff(j, B0 + step);
//           double fPP = evaluate(tmpFirst, tmpSecond);
//           /* +1,-1 node */
//           tem_->setFieldCoeff(i, A0 + step);
//           tem_->setFieldCoeff(j, B0 - step);
//           double fPM = evaluate(tmpFirst, tmpSecond);
//           fdSecond[i][j] = (fPP + fMM - fPM - fMP)/4.0/step/step;
//           fdSecond[j][i] = fdSecond[i][j];
      fdSecond[i][j] = (firstPlus[j] - firstMinus[j])/2.0/step;
      error2 = fabs(fdSecond[i][j] - afdSecond[i][j])/(step + fabs(afdSecond[i][j]));
      Out::os() << tab2 << "(" << tem_->fieldName(i) << ", " << tem_->fieldName(j) 
                << ") exact mixed deriv=" << afdSecond[i][j]
                << "    fd=" << fdSecond[i][j] << "    |exact - fd|=" 
                << error2 << std::endl;
      if (error2 > tol2) 
      {
        isOK = false;
        Out::os() << tab2 << "mixed partial deriv calculation wrt fields "
                  << i << ", " << j << " FAILED" << std::endl;
        Tabs tab3;
        Out::os() << tab3 << "f'(A_i+h) = " << firstPlus[j] << std::endl;
        Out::os() << tab3 << "f'(A_i-h) = " << firstMinus[j] << std::endl;
        Out::os() << tab3 << "step=" << step << std::endl;
      }
      tem_->setFieldCoeff(i, A0);
      //   tem_->setFieldCoeff(j, B0);
    }
  }

  

  return f0;
}


