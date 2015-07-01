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


#include "PlayaAmesosSolver.hpp"
#include "PlayaEpetraVector.hpp"
#include "PlayaEpetraMatrix.hpp"



#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#endif

#include "Amesos.h"
#include "Amesos_BaseSolver.h"


using namespace Teuchos;

namespace Playa
{

AmesosSolver::AmesosSolver(const ParameterList& params)
  : LinearSolverBase<double>(params),
    kernel_()
{
  if (parameters().isParameter("Kernel"))
  {
    kernel_ = getParameter<string>(parameters(), "Kernel");
  }
  else
  {
    kernel_ = "Klu";
  }
}



SolverState<double> AmesosSolver::solve(const LinearOperator<double>& op, 
  const Vector<double>& rhs, 
  Vector<double>& soln) const
{
	Playa::Vector<double> bCopy = rhs.copy();
	Playa::Vector<double> xCopy = rhs.copy();

  Epetra_Vector* b = EpetraVector::getConcretePtr(bCopy);
  Epetra_Vector* x = EpetraVector::getConcretePtr(xCopy);

	Epetra_CrsMatrix& A = EpetraMatrix::getConcrete(op);

  Epetra_LinearProblem prob(&A, x, b);

  Amesos amFactory;
  RCP<Amesos_BaseSolver> solver 
    = rcp(amFactory.Create("Amesos_" + kernel_, prob));
  TEUCHOS_TEST_FOR_EXCEPTION(solver.get()==0, std::runtime_error, 
    "AmesosSolver::solve() failed to instantiate "
    << kernel_ << "solver kernel");

  int ierr = solver->Solve();
  
  soln = xCopy;

  SolverStatusCode state;
  std::string msg;

  switch(ierr)
  {
    case 0:
      state = SolveConverged;
      msg = "converged";
      break;
    default:
      state = SolveCrashed;
      msg = "amesos failed: ierr=" + Teuchos::toString(ierr);
  }

  SolverState<double> rtn(state, "Amesos solver " + msg, 0, 0);
  return rtn;
}

}
