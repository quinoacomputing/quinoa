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


#include "PDEOptLinearPDEConstrainedObj.hpp"
#include "Sundance.hpp"

using namespace Teuchos;
using namespace Playa;

namespace Sundance
{

LinearPDEConstrainedObj::LinearPDEConstrainedObj(
  const Functional& lagrangian,
  const Expr& stateVars,
  const Expr& stateVarVals,
  const Expr& adjointVars,
  const Expr& adjointVarVals,
  const Expr& designVar,
  const Expr& designVarVal,
  const LinearSolver<double>& solver,
  const RCP<IterCallbackBase>& iterCallback,
  int verb)
  : PDEConstrainedObjBase(lagrangian, tuple(stateVarVals),
    tuple(adjointVarVals), designVarVal, iterCallback, verb),
    stateProbs_(),
    adjointProbs_(),
    solvers_(tuple(solver))
{
  init(tuple(stateVars), tuple(adjointVars), designVar);
}


LinearPDEConstrainedObj::LinearPDEConstrainedObj(
  const Functional& lagrangian,
  const Array<Expr>& stateVars,
  const Array<Expr>& stateVarVals,
  const Array<Expr>& adjointVars,
  const Array<Expr>& adjointVarVals,
  const Expr& designVar,
  const Expr& designVarVal,
  const Array<LinearSolver<double> >& solvers,
  const RCP<IterCallbackBase>& iterCallback,
  int verb)
  : PDEConstrainedObjBase(lagrangian, stateVarVals,
    adjointVarVals, designVarVal, iterCallback, verb),
    stateProbs_(),
    adjointProbs_(),
    solvers_(solvers)
{
  init(stateVars, adjointVars, designVar);
}


LinearPDEConstrainedObj::LinearPDEConstrainedObj(
  const Functional& lagrangian,
  const Expr& stateVars,
  const Expr& stateVarVals,
  const Expr& adjointVars,
  const Expr& adjointVarVals,
  const Expr& designVar,
  const Expr& designVarVal,
  const LinearSolver<double>& solver,
  int verb)
  : PDEConstrainedObjBase(lagrangian, tuple(stateVarVals),
    tuple(adjointVarVals), designVarVal, verb),
    stateProbs_(),
    adjointProbs_(),
    solvers_(tuple(solver))
{
  init(tuple(stateVars), tuple(adjointVars), designVar);
}


LinearPDEConstrainedObj::LinearPDEConstrainedObj(
  const Functional& lagrangian,
  const Array<Expr>& stateVars,
  const Array<Expr>& stateVarVals,
  const Array<Expr>& adjointVars,
  const Array<Expr>& adjointVarVals,
  const Expr& designVar,
  const Expr& designVarVal,
  const Array<LinearSolver<double> >& solvers,
  int verb)
  : PDEConstrainedObjBase(lagrangian, stateVarVals,
    adjointVarVals, designVarVal, verb),
    stateProbs_(),
    adjointProbs_(),
    solvers_(solvers)
{
  init(stateVars, adjointVars, designVar);
}


void LinearPDEConstrainedObj::initEquations(
  const Array<Expr>& stateVars,
  const Array<Expr>& adjointVars,
  const Array<Array<Expr> >& fixedVarsInStateEqns,
  const Array<Array<Expr> >& fixedVarsInStateEqnsVals,
  const Array<Array<Expr> >& fixedVarsInAdjointEqns,
  const Array<Array<Expr> >& fixedVarsInAdjointEqnsVals
  )
{
  Tabs tab(0);
  PLAYA_MSG2(verb(), tab << "setting up linear equations");
  
  for (int i=0; i<stateVars.size(); i++)
  {
    Tabs tab1;
    PLAYA_MSG3(verb(), tab1 << "setting up state equation #" << i);
    Expr fixedVars = new ListExpr(fixedVarsInStateEqns[i]);
    Expr fixedVarVals = new ListExpr(fixedVarsInStateEqnsVals[i]);
    LinearProblem stateProb 
      = Lagrangian().linearVariationalProb(adjointVars[i], adjointVarVals(i),
        stateVars[i],
        fixedVars, fixedVarVals);
                                   
    stateProbs_.append(stateProb);
  }

  for (int i=0; i<adjointVars.size(); i++)
  {
    Tabs tab1;
    PLAYA_MSG3(verb(), tab1 << "setting up adjoint equation #" << i);
    Expr fixedVars = new ListExpr(fixedVarsInAdjointEqns[i]);
    Expr fixedVarVals = new ListExpr(fixedVarsInAdjointEqnsVals[i]);
    LinearProblem adjointProb 
      = Lagrangian().linearVariationalProb(stateVars[i], stateVarVals(i),
        adjointVars[i],
        fixedVars, fixedVarVals);
                                   
    adjointProbs_.append(adjointProb);
  }

  PLAYA_MSG2(verb(), tab << "done setting up linear equations");
}




void LinearPDEConstrainedObj::solveState(const Vector<double>& x) const
{
  Tabs tab(0);
  PLAYA_MSG2(verb(), tab << "solving state"); 
  PLAYA_MSG3(verb(), tab << "|x|=" << x.norm2()); 
  PLAYA_MSG5(verb(), tab << "x=" << endl << tab << x.norm2());
  setDiscreteFunctionVector(designVarVal(), x);

  /* solve the state equations in order */
  for (int i=0; i<stateProbs_.size(); i++)
  {
    SolverState<double> status 
      = stateProbs_[i].solve(solvers_[i], stateVarVals(i));
    TEUCHOS_TEST_FOR_EXCEPTION(status.finalState() != SolveConverged,
      std::runtime_error,
      "state equation could not be solved: status="
      << status.stateDescription());
  }

  PLAYA_MSG2(verb(), tab << "done state solve"); 
  /* do postprocessing */
  statePostprocCallback();
}



void LinearPDEConstrainedObj
::solveStateAndAdjoint(const Vector<double>& x) const
{
  Tabs tab(0);
  PLAYA_MSG2(verb(), tab << "solving state and adjoint"); 
  PLAYA_MSG3(verb(), tab << "|x|=" << x.norm2()); 
  PLAYA_MSG5(verb(), tab << "x=" << endl << tab << x.norm2()); 

  Tabs tab1;
  setDiscreteFunctionVector(designVarVal(), x);

  PLAYA_MSG3(verb(), tab1 << "solving state eqns");
  /* solve the state equations in order */
  for (int i=0; i<stateProbs_.size(); i++)
  {
    SolverState<double> status 
      = stateProbs_[i].solve(solvers_[i], stateVarVals(i));

    /* if the solve failed, write out the design var and known state
     * variables */
    if (status.finalState() != SolveConverged)
    {
      FieldWriter w = new VTKWriter("badSolve");
      w.addMesh(Lagrangian().mesh());
      w.addField("designVar", new ExprFieldWrapper(designVarVal()));
      for (int j=0; j<i; j++)
      {
        Expr tmp = stateVarVals(j).flatten();
        for (int k=0; k<tmp.size(); k++)
        {
          w.addField("stateVar-"+Teuchos::toString(j)+"-"+Teuchos::toString(k),
            new ExprFieldWrapper(tmp[k]));
        }
      }
      w.write();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(status.finalState() != SolveConverged,
      std::runtime_error,
      "state equation " << i 
      << " could not be solved: status="
      << status.stateDescription());
  }

  PLAYA_MSG3(verb(), tab1 << "done solving state eqns");

  /* do postprocessing */
  statePostprocCallback();

  PLAYA_MSG3(verb(), tab1 << "solving adjoint eqns");

  /* solve the adjoint equations in reverse order */
  for (int i=adjointProbs_.size()-1; i>=0; i--)
  {
    SolverState<double> status 
      = adjointProbs_[i].solve(solvers_[i], adjointVarVals(i));

    /* if the solve failed, write out the design var and known state
     * and adjoint variables */
    if (status.finalState() != SolveConverged)
    {
      FieldWriter w = new VTKWriter("badSolve");
      w.addMesh(Lagrangian().mesh());
      w.addField("designVar", new ExprFieldWrapper(designVarVal()));
      for (int j=0; j<stateProbs_.size(); j++)
      {
        Expr tmp = stateVarVals(j).flatten();
        for (int k=0; k<tmp.size(); k++)
        {
          w.addField("stateVar-"+Teuchos::toString(j)+"-"+Teuchos::toString(k),
            new ExprFieldWrapper(tmp[k]));
        }
      }
      for (int j=adjointProbs_.size()-1; j>i; j--)
      {
        Expr tmp = adjointVarVals(j).flatten();
        for (int k=0; k<tmp.size(); k++)
        {
          w.addField("adjointVar-"+Teuchos::toString(j)+"-"+Teuchos::toString(k),
            new ExprFieldWrapper(tmp[k]));
        }

      }
      w.write();

    }
    TEUCHOS_TEST_FOR_EXCEPTION(status.finalState() != SolveConverged,
      std::runtime_error,
      "adjoint equation " << i 
      << " could not be solved: status="
      << status.stateDescription());
  }
  PLAYA_MSG3(verb(), tab1 << "done solving adjoint eqns");
  PLAYA_MSG2(verb(), tab1 << "done solving state and adjoint eqns");
}

  


}










