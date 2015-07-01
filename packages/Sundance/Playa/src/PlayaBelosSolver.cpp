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


#include "PlayaBelosSolver.hpp"
#include "PlayaPreconditioner.hpp"
#include "PlayaPreconditionerFactory.hpp"
#include "PlayaParameterListPreconditionerFactory.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;


BelosSolver::BelosSolver(const ParameterList& params)
  : LinearSolverBase<double>(params), pf_(), hasSolver_(false)
{
  if (params.isSublist("Preconditioner"))
  {
    ParameterList precParams = params.sublist("Preconditioner");
    pf_ = new ParameterListPreconditionerFactory(precParams);
  }
}



SolverState<double> BelosSolver::solve(const LinearOperator<double>& A, 
  const Vector<double>& rhs, 
  Vector<double>& soln) const
{
  typedef Anasazi::SimpleMV                      MV;
  typedef LinearOperator<double>                 OP;
  typedef Belos::LinearProblem<double, MV, OP>   LP;

  TEUCHOS_TEST_FOR_EXCEPT(!A.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPT(!rhs.ptr().get());

  if (!soln.ptr().get()) 
  {
    soln = rhs.copy();
    /* KRL 8 Jun 2012: set x0 to zero to workaround bug in Belos */
    soln.zero();
  }

  if (rhs.norm2()==0.0)
  {
    soln.zero();
    SolverStatusCode code = SolveConverged;
    SolverState<double> state(code, "Detected trivial solution", 0, 0.0);
    
    return state;
  }


  RCP<OP> APtr = rcp(new LinearOperator<double>(A));
  RCP<MV> bPtr = rcp(new MV(1));
  (*bPtr)[0] = rhs;
  RCP<MV> ansPtr = rcp(new MV(1));
  (*ansPtr)[0] = soln;
  
  
  RCP<LP> prob = rcp(new LP(APtr, ansPtr, bPtr));

  TEUCHOS_TEST_FOR_EXCEPT(!prob->setProblem());

  
  if (pf_.ptr().get())
  {
    Preconditioner<double> P = pf_.createPreconditioner(A);
    if (P.hasLeft())
    {
      prob->setLeftPrec(rcp(new OP(P.left())));
    }
  
    if (P.hasRight())
    {
      prob->setRightPrec(rcp(new OP(P.right())));
    }
  }

  if (!hasSolver_)
  {

    ParameterList plist = parameters();

    RCP<ParameterList> belosList = rcp(&plist, false);

    std::string solverType = parameters().get<string>("Method");
      
    if (solverType=="GMRES")
    {
      solver_=rcp(new Belos::BlockGmresSolMgr<double, MV, OP>(prob, belosList));
    }
    else if (solverType=="CG")
    {
      solver_=rcp(new Belos::BlockCGSolMgr<double, MV, OP>(prob, belosList));
    }
    else if (solverType=="TFQMR")
    {
      solver_=rcp(new Belos::TFQMRSolMgr<double, MV, OP>(prob, belosList));
    }
    else if (solverType=="GCRODR")
    {
      solver_=rcp(new Belos::GCRODRSolMgr<double, MV, OP>(prob, belosList));
      hasSolver_ = true; // only cache recycling solvers
    }
    else if (solverType=="RCG")
    {
      solver_=rcp(new Belos::RCGSolMgr<double, MV, OP>(prob, belosList));
      hasSolver_ = true; // only cache recycling solvers
    }
    else
    {
      TEUCHOS_TEST_FOR_EXCEPT(!(solverType=="GMRES" || solverType=="CG"));
    }
  }
  else // reset problem
  {
    solver_->setProblem( prob );
  }
  
  Belos::ReturnType rtn = solver_->solve();

  int numIters = solver_->getNumIters();
  double resid = solver_->achievedTol();
  
  SolverStatusCode code = SolveFailedToConverge;
  if (rtn==Belos::Converged) code = SolveConverged;
  SolverState<double> state(code, "Belos solver completed", numIters, resid);
  
  return state;
}



