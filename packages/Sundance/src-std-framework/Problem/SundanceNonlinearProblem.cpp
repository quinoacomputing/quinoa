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

#include "SundanceNonlinearProblem.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceLinearSolveDriver.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#endif


using namespace Sundance;
using namespace Teuchos;
using namespace std;
using namespace Playa;


static Time& nlpCtorTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("NonlinearProblem ctor"); 
  return *rtn;
}


NonlinearProblem::NonlinearProblem() 
  : op_()
{
  TimeMonitor timer(nlpCtorTimer());
}


NonlinearProblem::NonlinearProblem(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const Expr& test, 
  const Expr& unk, 
  const Expr& u0, 
  const VectorType<double>& vecType)
  : op_(rcp(new NLOp(mesh, eqn, bc, test, unk, u0, vecType)))
{}

NonlinearProblem::NonlinearProblem(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const Expr& test, 
  const Expr& unk, 
  const Expr& u0, 
  const Expr& params, 
  const Expr& paramValues,  
  const VectorType<double>& vecType)
  : op_(rcp(new NLOp(mesh, eqn, bc, test, unk, u0, 
        params, paramValues, vecType)))
{}

NonlinearProblem::NonlinearProblem(
  const Mesh& mesh, const Expr& eqn, const Expr& bc,
  const BlockArray& test, const BlockArray& unk, const Expr& u0)
  : op_(rcp(new NLOp(mesh, eqn, bc, test, unk, u0)))
{}


NonlinearProblem::NonlinearProblem(const RCP<Assembler>& assembler, 
                                   const Expr& u0)
  : op_(rcp(new NLOp(assembler, u0)))
{}



SolverState<double>
NonlinearProblem::solve(const NOXSolver& solver) const
{
  RCP<NonlinearOperatorBase<double> > op = op_;
  NonlinearOperator<double> F = op;
  Vector<double> soln;
  SolverState<double> rtn = solver.solve(F, soln);
  F.setEvalPt(soln);
  return rtn;
}



SolverState<double>
NonlinearProblem::solve(const NonlinearSolver<double>& solver) const
{
  RCP<NonlinearOperatorBase<double> > op = op_;
  NonlinearOperator<double> F = op;
  Vector<double> soln;
  SolverState<double> rtn = solver.solve(F, soln);
  F.setEvalPt(soln);
  return rtn;
}

