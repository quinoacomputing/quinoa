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

#include "SundanceLinearEigenproblem.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceGaussianQuadrature.hpp"

#include "PlayaLinearCombinationDecl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaSimpleDiagonalOpDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaSimpleDiagonalOpImpl.hpp"
#endif

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;

static Time& normalizationTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("Eigenfunction normalization"); 
  return *rtn;
}

static Time& makeEigensystemTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("Building eigensystem stiffness matrix"); 
  return *rtn;
}

LinearEigenproblem::LinearEigenproblem(
  const Mesh& mesh, const Expr& eqn,
  const Expr& v, const Expr& u,
  const VectorType<double>& vecType)
  : lumpMass_(false),
    kProb_(),
    mProb_(),
    M_(),
    MUnlumped_(),
    discSpace_()
{
  Expr empty;
  
  kProb_ = LinearProblem(mesh, eqn, empty, v, u, vecType);
  discSpace_ = *(kProb_.solnSpace()[0]);
}    

LinearEigenproblem::LinearEigenproblem(
  const Mesh& mesh, const Expr& eqn,
  const Expr& v, const Expr& u,
  const VectorType<double>& vecType,
  bool lumpedMass)
  : lumpMass_(lumpedMass),
    kProb_(),
    mProb_(),
    M_(),
    MUnlumped_(),
    discSpace_()
{
  Expr empty;
  
  kProb_ = LinearProblem(mesh, eqn, empty, v, u, vecType);
  mProb_ = makeMassProb(mesh, empty, v, u, vecType);
  discSpace_ = *(kProb_.solnSpace()[0]);
  MUnlumped_ = mProb_.getOperator();
  if (lumpMass_)
  {
    M_ = lumpedOperator(MUnlumped_);
  }
  else
  {
    M_ = MUnlumped_;
  }
}    


LinearEigenproblem::LinearEigenproblem(
  const Mesh& mesh, const Expr& eqn,
  const Expr& massExpr,
  const Expr& v, const Expr& u,
  const VectorType<double>& vecType,
  bool lumpedMass)
  : lumpMass_(lumpedMass),
    kProb_(),
    mProb_(),
    M_(),
    MUnlumped_(),
    discSpace_()
{
  Expr bc;
  kProb_ = LinearProblem(mesh, eqn, bc, v, u, vecType);
  mProb_ = makeMassProb(mesh, massExpr, v, u, vecType);
  discSpace_ = *(kProb_.solnSpace()[0]);

  MUnlumped_ = mProb_.getOperator();
  if (lumpMass_)
  {
    M_ = lumpedOperator(MUnlumped_);
  }
  else
  {
    M_ = MUnlumped_;
  }
  
}    

LinearProblem LinearEigenproblem::makeMassProb(
  const Mesh& mesh,
  const Expr& massExpr,
  const Expr& v, const Expr& u,
  const VectorType<double>& vecType) const
{
  Expr eqn;
  
  CellFilter interior = new MaximalCellFilter();
  QuadratureFamily quad = new GaussianQuadrature( 4 );
  if (massExpr.ptr().get()==0)
  {
    eqn = Integral(interior, v*u, quad);
  }
  else
  {
    eqn = Integral(interior, massExpr, quad);
  }
  Expr bc;
  LinearProblem rtn(mesh, eqn, bc, v, u, vecType);
  return rtn;
}


Array<Expr> LinearEigenproblem::makeEigenfunctions(
  Array<Vector<double> >& ev) const 
{
  TimeMonitor timer(normalizationTimer());

  Array<Expr> x(ev.size());
  CellFilter interior = new MaximalCellFilter();
  QuadratureFamily q = new GaussianQuadrature(2);
  for (int i=0; i<ev.size(); i++) 
  {
    x[i] = new DiscreteFunction(discSpace_, ev[i], "ev[" + Teuchos::toString(i)+"]");
    double N = 1.0;
    if (MUnlumped_.ptr().get())
    {
      N = ev[i] * (MUnlumped_ * ev[i]);
    }
    else
    {
      N = evaluateIntegral(discSpace_.mesh(), 
        Integral(interior, x[i]*x[i], q));
    }
    ev[i].scale(1.0/sqrt(N));
  }

  return x;
}


LinearOperator<double> 
LinearEigenproblem::lumpedOperator(const LinearOperator<double>& M) const 
{
  Vector<double> ones = M.domain().createMember();
  ones.setToConstant(1.0);
  Vector<double> m = M * ones;
  LinearOperator<double> rtn = diagonalOperator(m);

  return rtn;
}


Eigensolution LinearEigenproblem::solve(const Eigensolver<double>& solver) const 
{
  Array<std::complex<double> > ew;
  Array<Vector<double> > ev;

  LinearOperator<double> K;
  {
    TimeMonitor timer(makeEigensystemTimer());
    K = kProb_.getOperator();
  }
  
  solver.solve(K, M_, ev, ew);

  Array<Expr> eigenfuncs = makeEigenfunctions(ev);

  return Eigensolution(eigenfuncs, ew);
}

