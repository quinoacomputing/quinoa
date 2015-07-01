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

#include "SundanceLinearSolveDriver.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearSolverImpl.hpp"
#endif




using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace Playa;
using namespace std;


SolverState<double> 
LinearSolveDriver::solve(const LinearSolver<double>& solver,
  const LinearOperator<double>& A,
  const Array<Vector<double> >& rhs,
  const Array<RCP<DiscreteSpace> >& solutionSpace,
  const Array<Array<string> >& names,
  int verb,
  Expr& soln) const
{
  Tabs tab(0);
  Array<Vector<double> > solnVec(rhs.size());
  SolverState<double> state;

  for (int i=0; i<rhs.size(); i++)
  {
    Tabs tab1;

    solnVec[i] = rhs[i].copy();
    
    SUNDANCE_MSG2(verb, tab1 << "solving with RHS #" << i 
      << " of " << rhs.size());

    state = solver.solve(A, rhs[i], solnVec[i]);
    
    SUNDANCE_MSG2(verb, tab1 << "solve completed with status="
      << state.stateDescription());

    /* deal with a failure to converge */
    if (state.finalState() != SolveConverged)
    {
      TeuchosOStringStream ss;
      ss << "Solve failed! state = "
         << state.stateDescription()
         << "\nmessage=" << state.finalMsg()
         << "\niters taken = " << state.finalIters()
         << "\nfinal residual = " << state.finalResid();

      /* If requested, write the bad matrix and vector */
      if (dumpBadMatrix())
      {
        if (A.ptr().get() != 0)
        {
          ofstream osA(badMatrixFilename().c_str());
          A.print(osA);
          ss << "\nmatrix written to " << badMatrixFilename();
        }
        else
        {
          ss << "\nthe matrix is null! Evil is afoot in your code...";
        }
        if (rhs[i].ptr().get() != 0)
        {
          ofstream osb(badVectorFilename().c_str());
          rhs[i].print(osb);
          ss << "\nRHS vector written to " << badVectorFilename();
        }
        else
        {
          ss << "\nthe RHS vector is null! Evil is afoot in your code...";
        }
      }
      
      /* If solve errors are fatal, throw an exception */
      TEUCHOS_TEST_FOR_EXCEPTION(solveFailureIsFatal(),
        std::runtime_error, TEUCHOS_OSTRINGSTREAM_GET_C_STR(ss));

      /* otherwise, return the state information */
      return state;
    }
  }
   
  /* Put the solution vector into a discrete function */
  if (soln.ptr().get()==0)
  {
    soln = formSolutionExpr(solnVec, solutionSpace, names, verb);
  }
  else
  {
    writeIntoSolutionExpr(solnVec, soln, verb);
  }

  return state;
      
}



Expr LinearSolveDriver::formSolutionExpr(
  const Array<Vector<double> >& solnVector,
  const Array<RCP<DiscreteSpace> >& solutionSpace,
  const Array<Array<string> >& names,
  int verb) const
{
  Array<Expr> cols(solnVector.size());

  for (int m=0; m<solnVector.size(); m++)
  {
    Array<Expr> col(solutionSpace.size());

    for (int i=0; i<col.size(); i++)
    {
      std::string name = names[m][i];
      if (col.size() > 1) name += "[" + Teuchos::toString(i) + "]";
      col[i] = new DiscreteFunction(*(solutionSpace[i]),
        solnVector[m].getBlock(i), name);
    }
    if (col.size() > 1)
    {
      cols[m] = new ListExpr(col);
    }
    else
    {
      cols[m] = col[0];
    }
  }

  if (cols.size() > 1)
  {
    return new ListExpr(cols);;
  }
  else
  {
    return cols[0];
  }
}


void LinearSolveDriver::writeIntoSolutionExpr(
  const Array<Vector<double> >& solnVector,
  Expr soln,
  int verb) const 
{
#ifdef BLAH
  TEUCHOS_TEST_FOR_EXCEPTION(soln.size() != solnVector.size(),
    std::runtime_error,
    "soln=" << soln << " soln.size()=" << soln.size() << " while solnVector.size()=" << solnVector.size());
#endif
  for (int i=0; i<solnVector.size(); i++)
  {
    Expr u = soln[i];
    setDiscreteFunctionVector(u, solnVector[i]);
  }
}
