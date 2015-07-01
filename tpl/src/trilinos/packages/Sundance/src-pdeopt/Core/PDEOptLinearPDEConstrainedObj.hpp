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


#ifndef PDEOPT_LINEARPDECONSTRAINEDOBJ_H
#define PDEOPT_LINEARPDECONSTRAINEDOBJ_H

#include "PDEOptPDEConstrainedObjBase.hpp"
#include "SundanceExpr.hpp"
#include "SundanceFunctional.hpp"
#include "SundanceLinearProblem.hpp"

namespace Sundance
{
using namespace Playa;


/**
 * LinearPDEConstrainedObj is a base class for objective functions of the
 * reduced-space variable where the constraint is a linear PDE in the
 * state variables. 
 */
class LinearPDEConstrainedObj : public PDEConstrainedObjBase
{
public:
  /**  */
  LinearPDEConstrainedObj(
    const Functional& lagrangian,
    const Expr& stateVars,
    const Expr& stateVarVals,
    const Expr& adjointVars,
    const Expr& adjointVarVals,
    const Expr& designVars,
    const Expr& designVarVals,
    const LinearSolver<double>& solver,
    int verb=0);

  /**  */
  LinearPDEConstrainedObj(
    const Functional& lagrangian,
    const Array<Expr>& stateVars,
    const Array<Expr>& stateVarVals,
    const Array<Expr>& adjointVars,
    const Array<Expr>& adjointVarVals,
    const Expr& designVars,
    const Expr& designVarVals,
    const Array<LinearSolver<double> >& solvers,
    int verb=0);

  /**  */
  LinearPDEConstrainedObj(
    const Functional& lagrangian,
    const Expr& stateVars,
    const Expr& stateVarVals,
    const Expr& adjointVars,
    const Expr& adjointVarVals,
    const Expr& designVars,
    const Expr& designVarVals,
    const LinearSolver<double>& solver,
    const RCP<IterCallbackBase>& iterCallback,
    int verb=0);

  /**  */
  LinearPDEConstrainedObj(
    const Functional& lagrangian,
    const Array<Expr>& stateVars,
    const Array<Expr>& stateVarVals,
    const Array<Expr>& adjointVars,
    const Array<Expr>& adjointVarVals,
    const Expr& designVars,
    const Expr& designVarVals,
    const Array<LinearSolver<double> >& solvers,
    const RCP<IterCallbackBase>& iterCallback,
    int verb=0);

  /** virtual dtor */
  virtual ~LinearPDEConstrainedObj(){;}



  /** Solve the sequence of state equations, followed by postprocessing.
   * At the end of this call, the system is ready for evaluation of
   * the objective function or solution of the adjoint equations. */
  void solveState(const Vector<double>& x) const;

  /** Solve the sequence of state equations, then do postprocessing,
   * then finally the adjoint equations in <b>reverse</b> order. 
   * At the end of this call, the system is ready for evaluation
   * of the objective function and its gradient. */
  void solveStateAndAdjoint(const Vector<double>& x) const;

  /** Set up the linear equations */
  void initEquations(
    const Array<Expr>& stateVars,
    const Array<Expr>& adjointVars,
    const Array<Array<Expr> >& fixedVarsInStateEqns,
    const Array<Array<Expr> >& fixedVarsInStateEqnsVals,
    const Array<Array<Expr> >& fixedVarsInAdjointEqns,
    const Array<Array<Expr> >& fixedVarsInAdjointEqnsVals
    );


private:

  Array<LinearProblem> stateProbs_;

  Array<LinearProblem> adjointProbs_;

  Array<LinearSolver<double> > solvers_;

};

}

#endif
