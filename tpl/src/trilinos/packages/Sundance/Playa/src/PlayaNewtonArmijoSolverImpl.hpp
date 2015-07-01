/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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

#ifndef PLAYA_NEWTON_ARMIJO_SOLVER_IMPL_HPP
#define PLAYA_NEWTON_ARMIJO_SOLVER_IMPL_HPP

#include "PlayaNewtonArmijoSolverDecl.hpp"
#include "PlayaNonlinearOperator.hpp"
#include "PlayaTabs.hpp"
#include "PlayaOut.hpp"
#include "Teuchos_ParameterList.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif


using std::setw;

namespace Playa
{



template <class Scalar> inline
NewtonArmijoSolver<Scalar>::NewtonArmijoSolver(
  const ParameterList& params, 
  const LinearSolver<Scalar>& linSolver)
    : NonlinearSolverBase<Scalar>(params),
      linSolver_(linSolver),
      tauR_(10.0*Teuchos::ScalarTraits<Scalar>::eps()),
      tauA_(10.0*Teuchos::ScalarTraits<Scalar>::eps()),
      alpha_(1.0e-4),
      stepReduction_(0.5),
      maxIters_(20),
      maxLineSearch_(20),
      verb_(0)
  {
    if (params.isParameter("Tau Relative")) tauR_ = params.get<Scalar>("Tau Relative");
    if (params.isParameter("Tau Absolute")) tauA_ = params.get<Scalar>("Tau Absolute");
    if (params.isParameter("Alpha")) alpha_ = params.get<Scalar>("Alpha");
    if (params.isParameter("Step Reduction")) stepReduction_ = params.get<Scalar>("Step Reduction");
    if (params.isParameter("Max Iterations")) maxIters_ = params.get<int>("Max Iterations");
    if (params.isParameter("Max Backtracks")) maxLineSearch_ = params.get<int>("Max Backtracks");
    if (params.isParameter("Verbosity")) verb_ = params.get<int>("Verbosity");
  }

template <class Scalar> inline
SolverState<Scalar> NewtonArmijoSolver<Scalar>::solve(const NonlinearOperator<Scalar>& F,
  Vector<Scalar>& soln) const  
{
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  typedef typename Teuchos::ScalarTraits<Scalar> ST;
  
  Tabs tab0(0);
  PLAYA_MSG1(verb_, tab0 << " begin Playa::NewtonArmijoSolver::solve()");

  soln = F.getInitialGuess().copy();
  Vector<Scalar> newtonStep = soln.copy();
  
  F.setEvalPt(soln);
  Vector<Scalar> resid = F.getFunctionValue();

  ScalarMag r0 = resid.norm2();
  ScalarMag normF0 = r0;

  for (int i=0; i<maxIters_; i++)
  {
    Tabs tab1;
    PLAYA_MSG2(verb_, tab1 << "Newton iter #" << setw(6) << i << " |F|=" << setw(12) << normF0 << " |F|/|F0|="
      << setw(12) << normF0/r0);
    
    if (normF0 < r0*tauR_ + tauA_)
    {
      PLAYA_MSG3(verb_, tab1 << "|F|=" << setw(12) << normF0);
      PLAYA_MSG3(verb_, tab1 << "Relative tolerance tauR=" << setw(12) << tauR_);
      PLAYA_MSG3(verb_, tab1 << "Absolute tolerance tauA=" << setw(12) << tauA_);
      PLAYA_MSG3(verb_, tab1 << "  F0*tauR+tauA=" << setw(12) << r0*tauR_ + tauA_);
      PLAYA_MSG2(verb_, tab1 << "converged!");
      PLAYA_MSG1(verb_, tab0 << " done Playa::NewtonArmijoSolver::solve()");
      soln = F.currentEvalPt().copy();
      return SolverState<Scalar>(SolveConverged, "NewtonArmijoSolver::solve converged",
        i, normF0);
    }
    LinearOperator<Scalar> J = F.getJacobian();


    SolverState<Scalar> linSolverState = linSolver_.solve(J, resid, newtonStep);
    if (linSolverState.finalState() != SolveConverged)
    {
      PLAYA_MSG1(verb_, tab0 << " done Playa::NewtonArmijoSolver::solve()");
      return SolverState<Scalar>(SolveCrashed, 
        "NewtonArmijoSolver::solve: linear solve failed with message [" 
        + linSolverState.finalMsg() + "]", i, normF0);
    }
      
    
    Scalar t = ST::one();

    bool stepAccepted = false;
    soln = F.currentEvalPt().copy();
    
    for (int j=0; j<maxLineSearch_; j++)
    {
      Tabs tab2;
      Vector<Scalar> tmp = soln - t*newtonStep;
      F.setEvalPt( tmp );
      resid = F.getFunctionValue();
      ScalarMag normF1 = resid.norm2();
      PLAYA_MSG2(verb_, tab2 << "step t=" << setw(12) << t << " |F|=" << setw(12) << normF1);
      if (normF1 < (ST::one() - alpha_*t)*normF0)
      {
        stepAccepted = true;
        normF0 = normF1;
        break;
      }
      t = stepReduction_*t;
    }
    
    if (!stepAccepted)
    {
      PLAYA_MSG1(verb_, tab0 << " done Playa::NewtonArmijoSolver::solve()");
      return SolverState<Scalar>(SolveCrashed, 
        "NewtonArmijoSolver: line search failed",i, normF0);
    }
  }
  
  PLAYA_MSG1(verb_, tab0 << " done Playa::NewtonArmijoSolver::solve()");
  return SolverState<Scalar>(SolveFailedToConverge, "NewtonArmijoSolver: convergence failure after "
    + Teuchos::toString(maxIters_) + " steps.", maxIters_, normF0); 
}

}


#endif
