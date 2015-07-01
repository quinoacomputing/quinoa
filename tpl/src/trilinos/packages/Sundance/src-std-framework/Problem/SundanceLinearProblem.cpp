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

#include "SundanceLinearProblem.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceExpr.hpp"
#include "SundanceListExpr.hpp"
#include "PlayaSolverState.hpp"
#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#endif


using namespace Sundance;
using namespace Teuchos;
using namespace Playa;
using namespace std;


static Time& lpCtorTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("LinearProblem ctor"); 
  return *rtn;
}


LinearProblem::LinearProblem() 
  : assembler_(),
    A_(),
    rhs_()
{
  TimeMonitor timer(lpCtorTimer());
}



LinearProblem::LinearProblem(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const Expr& test, 
  const Expr& unk, 
  const VectorType<double>& vecType)
  : assembler_(),
    A_(),
    rhs_(1),
    names_(1),
    solveDriver_(),
    params_()
{
  bool partitionBCs = false;
  TimeMonitor timer(lpCtorTimer());
  Expr u = unk.flattenSpectral();
  Expr v = test.flattenSpectral();

  Array<Expr> zero(u.size());
  for (int i=0; i<u.size(); i++) 
  {
    Expr z = new ZeroExpr();
    zero[i] = z;
    names_[0].append(u[i].toString());
  }

  Expr u0 = new ListExpr(zero);

  Expr unkParams;
  Expr fixedParams;
  Array<Expr> fixedFields;
  Expr unkParamValues;
  Expr fixedParamValues;
  Array<Expr> fixedFieldValues;


  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, tuple(v), tuple(u), tuple(u0),
        unkParams, unkParamValues,
        fixedParams, fixedParamValues,
        fixedFields, fixedFieldValues));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vecType), tuple(vecType), partitionBCs));
}


LinearProblem::LinearProblem(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const Expr& test, 
  const Expr& unk, 
  const Expr& params, 
  const Expr& paramVals, 
  const VectorType<double>& vecType)
  : assembler_(),
    A_(),
    rhs_(1),
    names_(1),
    solveDriver_(),
    params_(params)
{
  bool partitionBCs = false;
  TimeMonitor timer(lpCtorTimer());
  Expr u = unk.flattenSpectral();
  Expr v = test.flattenSpectral();

  Expr dumParams;
  Expr alpha = params.flattenSpectral();
  Expr alpha0 = paramVals.flattenSpectral();
  Array<Expr> zero(u.size());
  for (int i=0; i<u.size(); i++) 
  {
    Expr z = new ZeroExpr();
    zero[i] = z;
    names_[0].append(u[i].toString());
  }

  Expr u0 = new ListExpr(zero);

  Array<Expr> fixedFields;

  
  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, tuple(v), tuple(u), tuple(u0),
        dumParams, dumParams, 
        alpha, alpha0,
        fixedFields, fixedFields));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vecType), tuple(vecType), partitionBCs));
}



LinearProblem::LinearProblem(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const BlockArray& test, 
  const BlockArray& unk)
  : assembler_(),
    A_(),
    rhs_(1),
    names_(unk.size()),
    solveDriver_(),
    params_()
{
  bool partitionBCs = false;
  TimeMonitor timer(lpCtorTimer());
  Array<Expr> v(test.size());  
  Array<Expr> u(unk.size());
  Array<Expr> u0(unk.size());

  Array<VectorType<double> > testVecType(test.size());
  Array<VectorType<double> > unkVecType(unk.size());

  for (int i=0; i<test.size(); i++)
  {
    v[i] = test[i].expr().flattenSpectral();
    testVecType[i] = test[i].vecType();
  }

  for (int i=0; i<unk.size(); i++)
  {
    u[i] = unk[i].expr().flattenSpectral();
    unkVecType[i] = unk[i].vecType();
    Array<Expr> zero(u[i].size());
    for (int j=0; j<u[i].size(); j++) 
    {
      Expr z = new ZeroExpr();
      zero[j] = z;
      names_[i].append(u[i][j].toString());
    }
    u0[i] = new ListExpr(zero);

  }

  Expr unkParams;
  Expr fixedParams;
  Array<Expr> fixedFields;
  Expr unkParamValues;
  Expr fixedParamValues;
  Array<Expr> fixedFieldValues;
  
  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, v, u, u0,
        unkParams, unkParamValues,
        fixedParams, fixedParamValues,
        fixedFields, fixedFieldValues));

  assembler_ = rcp(new Assembler(mesh, eqnSet, testVecType, unkVecType, partitionBCs));
}


LinearProblem::LinearProblem(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const BlockArray& test, 
  const BlockArray& unk,
  const Expr& unkParams, 
  const Expr& unkParamVals)
  : assembler_(),
    A_(),
    rhs_(1),
    names_(unk.size()),
    solveDriver_(),
    params_(unkParams)
{
  bool partitionBCs = false;
  TimeMonitor timer(lpCtorTimer());
  Array<Expr> v(test.size());  
  Array<Expr> u(unk.size());
  Array<Expr> u0(unk.size());
  Array<VectorType<double> > testVecType(test.size());
  Array<VectorType<double> > unkVecType(unk.size());

  for (int i=0; i<test.size(); i++)
  {
    v[i] = test[i].expr().flattenSpectral();
    testVecType[i] = test[i].vecType();
  }

  for (int i=0; i<unk.size(); i++)
  {
    u[i] = unk[i].expr().flattenSpectral();
    unkVecType[i] = unk[i].vecType();
    Array<Expr> zero(u[i].size());
    for (int j=0; j<u[i].size(); j++) 
    {
      Expr z = new ZeroExpr();
      zero[j] = z;
      names_[i].append(u[i][j].toString());
    }
    u0[i] = new ListExpr(zero);
  }

  Expr fixedParams;
  Array<Expr> fixedFields;
  Expr fixedParamValues;
  Array<Expr> fixedFieldValues;
  
  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, v, u, u0,
        unkParams.flattenSpectral(), 
        unkParamVals.flattenSpectral(),
        fixedParams, fixedParamValues,
        fixedFields, fixedFieldValues));

  assembler_ = rcp(new Assembler(mesh, eqnSet, testVecType, unkVecType, partitionBCs));
}

LinearProblem::LinearProblem(const RCP<Assembler>& assembler)
  : assembler_(assembler),
    A_(),
    rhs_(1),
    names_(),
    params_()
{  
  TimeMonitor timer(lpCtorTimer());
  const RCP<EquationSet>& eqn = assembler->eqnSet();
  names_.resize(eqn->numUnkBlocks());
  for (int i=0; i<eqn->numUnkBlocks(); i++)
  {
    for (int j=0; j<eqn->numUnks(i); j++) 
    {
      names_[i].append("u(" + Teuchos::toString(i) + ", "
        + Teuchos::toString(j) + ")");
    }
  }
}

/* Return the map from cells and functions to row indices */
const RCP<DOFMapBase>& LinearProblem::rowMap(int blockRow) const 
{return assembler_->rowMap()[blockRow];}
    
/* Return the map from cells and functions to column indices */
const RCP<DOFMapBase>& LinearProblem::colMap(int blockCol) const 
{return assembler_->colMap()[blockCol];}

/* Return the discrete space in which solutions live */
const Array<RCP<DiscreteSpace> >& LinearProblem::solnSpace() const 
{return assembler_->solutionSpace();}
    
/* Return the set of row indices marked as 
 * essential boundary conditions */
const RCP<Set<int> >& LinearProblem::bcRows(int blockRow) const 
{return assembler_->bcRows()[blockRow];}

/* Return the number of block rows in the problem  */
int LinearProblem::numBlockRows() const {return assembler_->rowMap().size();}

/* Return the number of block cols in the problem  */
int LinearProblem::numBlockCols() const {return assembler_->colMap().size();}

Array<Vector<double> > LinearProblem::getRHS() const 
{
  Tabs tab;
  SUNDANCE_MSG1(assembler_->maxWatchFlagSetting("solve control"),
    tab << "LinearProblem::solve() building vector");
  assembler_->assemble(rhs_);
  return rhs_;
}


Playa::LinearOperator<double> LinearProblem::getOperator() const 
{
  Tabs tab;
  SUNDANCE_MSG1(assembler_->maxWatchFlagSetting("solve control"),
    tab << "LinearProblem::solve() building matrix and vector");
  assembler_->assemble(A_, rhs_);
  return A_;
}

Expr LinearProblem::solve(const LinearSolver<double>& solver) const 
{
  Tabs tab;
  int verb = assembler_->maxWatchFlagSetting("solve control");
  Array<Vector<double> > solnVec(rhs_.size());
  
  SUNDANCE_MSG1(verb, tab << "LinearProblem::solve() building system");

  assembler_->assemble(A_, rhs_);
  for (int i=0; i<rhs_.size(); i++)
    rhs_[i].scale(-1.0);

  SUNDANCE_MSG1(verb, tab << "LinearProblem::solve() solving system");

  Expr rtn;

  /* we're not checking the status of the solve, so failures should
   * be considered fatal */
  bool save = LinearSolveDriver::solveFailureIsFatal();
  LinearSolveDriver::solveFailureIsFatal() = true;

  solveDriver_.solve(solver, A_, rhs_, solnSpace(), names_, verb, rtn);

  SUNDANCE_MSG1(verb, tab << "LinearProblem::solve() done solving system");

  /* restore original failure-handling setting */
  LinearSolveDriver::solveFailureIsFatal()=save; 

  return rtn;
}

SolverState<double> LinearProblem
::solve(const LinearSolver<double>& solver,
  Expr& soln) const 
{
  Tabs tab;
  int verb = assembler_->maxWatchFlagSetting("solve control");

  Array<Vector<double> > solnVec(rhs_.size());
  
  SUNDANCE_MSG1(verb, tab << "LinearProblem::solve() building system");

  assembler_->assemble(A_, rhs_);
  for (int i=0; i<rhs_.size(); i++)
  {
    rhs_[i].scale(-1.0);
  }

  SUNDANCE_MSG1(verb, tab << "solving LinearProblem");
  
  return solveDriver_.solve(solver, A_, rhs_, solnSpace(), names_, verb, soln);
}


Expr LinearProblem::formSolutionExpr(const Array<Vector<double> >& vec) const
{
  int verb = assembler_->maxWatchFlagSetting("solve control");
  return solveDriver_.formSolutionExpr(vec, solnSpace(), names_, verb);
}


void LinearProblem::reAssembleProblem() const {
	assembler_->flushConfiguration();
}

