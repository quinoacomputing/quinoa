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


#include "PDEOptPDEConstrainedObjBase.hpp"
#include "Sundance.hpp"


using namespace Teuchos;
using namespace Playa;
using namespace Sundance;

using std::endl;
using std::runtime_error;


PDEConstrainedObjBase::PDEConstrainedObjBase(
  const Functional& lagrangian,
  const Array<Expr>& stateVarVals,
  const Array<Expr>& adjointVarVals,
  const Expr& designVarVal,
  int verb)
  : ObjectiveBase(verb),
    Lagrangian_(lagrangian),
    designVarVal_(designVarVal),
    stateVarVals_(stateVarVals),
    adjointVarVals_(adjointVarVals),
    fEval_(),
    numFuncEvals_(0),
    numGradEvals_(0),
    invHScale_(1.0),
    iterCallback_()
{}


PDEConstrainedObjBase::PDEConstrainedObjBase(
  const Functional& lagrangian,
  const Array<Expr>& stateVarVals,
  const Array<Expr>& adjointVarVals,
  const Expr& designVarVal,
  const RCP<IterCallbackBase>& iterCallback,
  int verb)
  : ObjectiveBase(verb),
    Lagrangian_(lagrangian),
    designVarVal_(designVarVal),
    stateVarVals_(stateVarVals),
    adjointVarVals_(adjointVarVals),
    fEval_(),
    numFuncEvals_(0),
    numGradEvals_(0),
    invHScale_(1.0),
    iterCallback_(iterCallback)
{}


void PDEConstrainedObjBase::init(
  const Array<Expr>& stateVars,
  const Array<Expr>& adjointVars,
  const Expr& designVar)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stateVars.size() != adjointVars.size(), 
    RuntimeError,
    "number of state and adjoint variables should be identical");

  TEUCHOS_TEST_FOR_EXCEPTION(stateVars.size() != stateVarVals_.size(), 
    RuntimeError,
    "number of state variables and state values "
    "should be identical");

  TEUCHOS_TEST_FOR_EXCEPTION(adjointVars.size() != adjointVarVals_.size(), 
    RuntimeError,
    "number of adjoint variables and adjoint values "
    "should be identical");

  Array<Array<Expr> > fixedVarsInStateEqns(stateVars.size());
  Array<Array<Expr> > fixedVarsInStateEqnsVals(stateVars.size());

  Array<Array<Expr> > fixedVarsInAdjointEqns(stateVars.size());
  Array<Array<Expr> > fixedVarsInAdjointEqnsVals(stateVars.size());

  for (int i=0; i<stateVars.size(); i++)
  {
    /* hold all the other variables fixed */
    fixedVarsInStateEqns[i].append(designVar);
    fixedVarsInStateEqnsVals[i].append(designVarVal_);
    for (int j=0; j<stateVars.size(); j++)
    {
      if (i != j)
      {
        fixedVarsInStateEqns[i].append(stateVars[j]);
        fixedVarsInStateEqnsVals[i].append(stateVarVals_[j]);
        fixedVarsInAdjointEqns[i].append(adjointVars[j]);
        fixedVarsInAdjointEqnsVals[i].append(adjointVarVals_[j]);
        fixedVarsInStateEqns[i].append(adjointVars[j]);
        fixedVarsInStateEqnsVals[i].append(adjointVarVals_[j]);
        fixedVarsInAdjointEqns[i].append(stateVars[j]);
        fixedVarsInAdjointEqnsVals[i].append(stateVarVals_[j]);   
      }
    }
  }

  initEquations(stateVars, adjointVars, 
    fixedVarsInStateEqns, fixedVarsInStateEqnsVals,
    fixedVarsInAdjointEqns, fixedVarsInAdjointEqnsVals);

  Array<Expr> allVarArray = stateVars;
  Array<Expr> allVarValArray = stateVarVals_;
  for (int i=0; i<adjointVars.size(); i++)
  {
    allVarArray.append(adjointVars[i]);
    allVarValArray.append(adjointVarVals_[i]);
  }
  Expr allVars = new ListExpr(allVarArray);
  Expr allVarVals = new ListExpr(allVarValArray);

  fEval_ = Lagrangian_.evaluator(designVar, designVarVal_,
    allVars, allVarVals);
}


void PDEConstrainedObjBase::evalGrad(const Vector<double>& x, double& f, 
  Vector<double>& grad) const
{
  Tabs tabs(0);
  PLAYA_MSG2(verb(), tabs << "in evalGrad()"); 

  /* Solve for the state and adjoint variables */
  solveStateAndAdjoint(x);

  /* Compute the gradient using the adjoint method */
  Expr gradExpr = fEval_.evalGradient(f);
  grad = getDiscreteFunctionVector(gradExpr).copy();

  PLAYA_MSG2(verb(), tabs << "function value = " << f);
  PLAYA_MSG5(verb(), tabs << "|gradient| is " << grad.norm2());
  PLAYA_MSG5(verb(), tabs << "gradient is " << grad);
  numFuncEvals_++;
  numGradEvals_++;
}

void PDEConstrainedObjBase::eval(const Vector<double>& x, double& f) const  
{
  Tabs tabs(0);
  PLAYA_MSG2(verb(), tabs << "in eval()");

  /* Solve for the state variables */
  solveState(x);

  /* Set all the adjoints to zero. We do this because the residual of the state
   * equation might not be exactly zero. This way, we can evaluate the Lagrangian
   * an have zero contribution from the constraint term. */
  for (int i=0; i<adjointVarVals_.size(); i++)
  {

    Vector<double> adjVec = getDiscreteFunctionVector(adjointVarVals_[i]);
    adjVec.zero();
    setDiscreteFunctionVector(adjointVarVals_[i], adjVec);
  }

  /* Evaluate the functional */
  f = fEval_.evaluate();
  PLAYA_MSG2(verb(), tabs << "function value = " << f);
  numFuncEvals_++;
}


Vector<double> PDEConstrainedObjBase::getInit() const
{
  return getDiscreteFunctionVector(designVarVal());
}

void PDEConstrainedObjBase:: iterationCallback(const Vector<double>& x, 
  int iter) const
{
  if (iterCallback_.ptr().get() != 0) 
  {
    iterCallback_->call(this, iter);
  }
}











