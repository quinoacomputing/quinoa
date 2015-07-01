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

#include "SundanceNLOp.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceLinearSolveDriver.hpp"
#include "PlayaLinearSolverDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#endif

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace std;
using namespace Playa;


static Time& nlpCtorTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("NLOp ctor"); 
  return *rtn;
}


NLOp::NLOp() 
  : NonlinearOperatorBase<double>(),
    assembler_(),
    u0_(),
    params_()

{
  TimeMonitor timer(nlpCtorTimer());
}


NLOp::NLOp(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const Expr& test, 
  const Expr& unk, 
  const Expr& u0, 
  const VectorType<double>& vecType,
  bool partitionBCs)
  : NonlinearOperatorBase<double>(),
    assembler_(),
    u0_(u0),
    params_()
{
  TimeMonitor timer(nlpCtorTimer());

  Expr unkParams;
  Expr fixedParams;
  Expr fixedFields;
  Expr unkParamValues;
  Expr fixedParamValues;
  Expr fixedFieldValues;

  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, tuple(test.flattenSpectral()), tuple(unk.flattenSpectral()), tuple(u0),
        unkParams, unkParamValues,
        fixedParams, fixedParamValues,
        tuple(fixedFields), tuple(fixedFieldValues)));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vecType), tuple(vecType), partitionBCs));

  VectorSpace<double> domain = assembler_->solnVecSpace();
  VectorSpace<double> range = assembler_->rowVecSpace();

  setDomainAndRange(domain, range);
}



NLOp::NLOp(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const BlockArray& test, 
  const BlockArray& unk, 
  const Expr& u0)
  : NonlinearOperatorBase<double>(),
    assembler_(),
    u0_(u0),
    params_()
{
  TimeMonitor timer(nlpCtorTimer());
  bool partitionBCs = false;

  Array<Expr> v(test.size());  
  Array<Expr> u(unk.size());
  Array<Expr> u0Array(unk.size());

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
    u0Array[i] = u0[i].flattenSpectral();
    unkVecType[i] = unk[i].vecType(); 
  }

  Expr unkParams;
  Expr fixedParams;
  Expr fixedFields;
  Expr unkParamValues;
  Expr fixedParamValues;
  Expr fixedFieldValues;

  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(eqn, bc, 
        v, u, u0Array,
        unkParams, unkParamValues,
        fixedParams, fixedParamValues,
        tuple(fixedFields), tuple(fixedFieldValues)));

  assembler_ = rcp(new Assembler(mesh, eqnSet, testVecType, unkVecType, partitionBCs));

  VectorSpace<double> domain = assembler_->solnVecSpace();
  VectorSpace<double> range = assembler_->rowVecSpace();

  setDomainAndRange(domain, range);
}



NLOp::NLOp(const Mesh& mesh, 
  const Expr& eqn, 
  const Expr& bc,
  const Expr& test, 
  const Expr& unk, 
  const Expr& u0, 
  const Expr& params, 
  const Expr& paramValues,  
  const VectorType<double>& vecType,
  bool partitionBCs)
  : NonlinearOperatorBase<double>(),
    assembler_(),
    J_(),
    u0_(u0),
    params_(params)
{
  TimeMonitor timer(nlpCtorTimer());

  Expr unkParams;
  Expr fixedFields;
  Expr unkParamValues;
  Expr fixedFieldValues;

  RCP<EquationSet> eqnSet 
    = rcp(new EquationSet(
            eqn, bc, tuple(test.flattenSpectral()), 
            tuple(unk.flattenSpectral()), tuple(u0), 
            unkParams, unkParamValues,
            params, paramValues,
            tuple(fixedFields), tuple(fixedFieldValues)));

  assembler_ = rcp(new Assembler(mesh, eqnSet, tuple(vecType), tuple(vecType), partitionBCs));

  VectorSpace<double> domain = assembler_->solnVecSpace();
  VectorSpace<double> range = assembler_->rowVecSpace();

  setDomainAndRange(domain, range);
}


NLOp::NLOp(const RCP<Assembler>& assembler, 
  const Expr& u0)
  : NonlinearOperatorBase<double>(),
    assembler_(assembler),
    J_(),
    u0_(u0),
    params_()
{
  TimeMonitor timer(nlpCtorTimer());

  VectorSpace<double> domain = assembler_->solnVecSpace();
  VectorSpace<double> range = assembler_->rowVecSpace();

  setDomainAndRange(domain, range);
}


void NLOp::updateDiscreteFunctionValue(const Vector<double>& vec) const
{
  setDiscreteFunctionVector(u0_, vec.copy());
}

Vector<double> NLOp::getInitialGuess() const 
{
  return getDiscreteFunctionVector(u0_);
}

void NLOp::setInitialGuess(const Expr& u0New) 
{
  Vector<double> vec = getDiscreteFunctionVector(u0New);
  setEvalPt(vec);
  updateDiscreteFunctionValue(vec);
}


LinearOperator<double> NLOp::allocateJacobian() const
{
  return assembler_->allocateMatrix();
}


LinearOperator<double> NLOp::
computeJacobianAndFunction(Vector<double>& functionValue) const
{
  updateDiscreteFunctionValue(currentEvalPt());

  Array<Vector<double> > mv(1);
  mv[0].acceptCopyOf(functionValue);
  assembler_->assemble(J_, mv);
  functionValue.acceptCopyOf(mv[0]);

  return J_;
}

void NLOp::
computeJacobianAndFunction(LinearOperator<double>& J,
  Vector<double>& resid) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, std::runtime_error,
    "null evaluation point in "
    "NLOp::jacobian()");

  TEUCHOS_TEST_FOR_EXCEPTION(J.ptr().get()==0, std::runtime_error,
    "null Jacobian pointer in "
    "NLOp::jacobian()");

  TEUCHOS_TEST_FOR_EXCEPTION(resid.ptr().get()==0, std::runtime_error,
    "null residual pointer in "
    "NLOp::jacobian()");

  Array<Vector<double> > mv(1);
  mv[0] = resid;

  updateDiscreteFunctionValue(currentEvalPt());
  assembler_->assemble(J, mv);

  resid.acceptCopyOf(mv[0]);

  J_ = J;
}




Vector<double> NLOp::computeFunctionValue() const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, std::runtime_error,
    "null evaluation point in "
    "NLOp::computeFunctionValue()");

  VectorSpace<double> spc = range();
  Vector<double> rtn = spc.createMember();

  Array<Vector<double> > mv(1);
  mv[0] = rtn;

  updateDiscreteFunctionValue(currentEvalPt());
  assembler_->assemble(mv);

  rtn.acceptCopyOf(mv[0]);

  return rtn;
}



void NLOp::computeFunctionValue(Vector<double>& resid) const 
{
  /* Set the vector underlying the discrete 
   * function to the evaluation point*/

  TEUCHOS_TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, std::runtime_error,
    "null evaluation point in "
    "NLOp::computeFunctionValue()");

  TEUCHOS_TEST_FOR_EXCEPTION(resid.ptr().get()==0, std::runtime_error,
    "null residual pointer in "
    "NLOp::computeFunctionValue()");

  updateDiscreteFunctionValue(currentEvalPt());


  Array<Vector<double> > mv(1);
  mv[0] = resid;

  assembler_->assemble(mv);

  resid.acceptCopyOf(mv[0]);
}



Expr NLOp::
computeSensitivities(const LinearSolver<double>& solver) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(currentEvalPt().ptr().get()==0, std::runtime_error,
    "null evaluation point in "
    "NLOp::computeSensitivities()");

  TEUCHOS_TEST_FOR_EXCEPTION(J_.ptr().get()==0, std::runtime_error,
    "null Jacobian pointer in "
    "NLOp::computeSensitivities()");

  TEUCHOS_TEST_FOR_EXCEPTION(params_.ptr().get()==0, std::runtime_error,
    "null parameters in NLOp::computeSensitivities()");

  updateDiscreteFunctionValue(currentEvalPt());

  Array<Vector<double> > mv(params_.size());

  assembler_->assembleSensitivities(J_, mv);

  LinearSolveDriver driver;

  Expr sens;
  int vrb = 0;
  Array<Array<string> > names(params_.size());
  for (int i=0; i<params_.size(); i++)
  {
    names[i].resize(u0_.size());
    for (int j=0; j<u0_.size(); j++) 
      names[i][j]="sens(" + u0_[j].toString() + ", " + params_[i].toString() + ")";
    mv[i].scale(-1.0);
  }

  driver.solve(solver, J_, mv, assembler_->solutionSpace(), names, vrb, sens);
  return sens;
}


void NLOp::reAssembleProblem() const
{
	assembler_->flushConfiguration();
}
