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

#ifndef PLAYA_BLOCKTRIANGULARSOLVER_IMPL_HPP
#define PLAYA_BLOCKTRIANGULARSOLVER_IMPL_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearSolverDecl.hpp" 
#include "PlayaLinearCombinationImpl.hpp" 
#include "PlayaSimpleZeroOpDecl.hpp" 
#include "PlayaBlockTriangularSolverDecl.hpp" 


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearSolverImpl.hpp" 
#include "PlayaSimpleZeroOpDecl.hpp" 
#endif

namespace Playa
{
using namespace PlayaExprTemplates;

template <class Scalar> inline
BlockTriangularSolver<Scalar>
::BlockTriangularSolver(const LinearSolver<Scalar>& solver)
  : LinearSolverBase<Scalar>(ParameterList()), solvers_(tuple(solver)) {;}

template <class Scalar> inline
BlockTriangularSolver<Scalar>
::BlockTriangularSolver(const Array<LinearSolver<Scalar> >& solvers)
  : LinearSolverBase<Scalar>(ParameterList()), solvers_(solvers) {;}

template <class Scalar> inline
SolverState<Scalar> BlockTriangularSolver<Scalar>
::solve(const LinearOperator<Scalar>& op,
  const Vector<Scalar>& rhs,
  Vector<Scalar>& soln) const
{
  int nRows = op.numBlockRows();
  int nCols = op.numBlockCols();

  soln = op.domain().createMember();
  //    bool converged = false;

  TEUCHOS_TEST_FOR_EXCEPTION(nRows != rhs.space().numBlocks(), std::runtime_error,
    "number of rows in operator " << op
    << " not equal to number of blocks on RHS "
    << rhs);

  TEUCHOS_TEST_FOR_EXCEPTION(nRows != nCols, std::runtime_error,
    "nonsquare block structure in block triangular "
    "solver: nRows=" << nRows << " nCols=" << nCols);

  bool isUpper = false;
  bool isLower = false;

  for (int r=0; r<nRows; r++)
  {
    for (int c=0; c<nCols; c++)
    {
      if (op.getBlock(r,c).ptr().get() == 0 ||
        dynamic_cast<const SimpleZeroOp<Scalar>* >(op.getBlock(r,c).ptr().get()))
      {
        TEUCHOS_TEST_FOR_EXCEPTION(r==c, std::runtime_error,
          "zero diagonal block (" << r << ", " << c 
          << " detected in block "
          "triangular solver. Operator is " << op);
        continue;
      }
      else
      {
        if (r < c) isUpper = true;
        if (c < r) isLower = true;
      }
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION(isUpper && isLower, std::runtime_error, 
    "block triangular solver detected non-triangular operator "
    << op);

  bool oneSolverFitsAll = false;
  if ((int) solvers_.size() == 1 && nRows != 1) 
  {
    oneSolverFitsAll = true;
  }

  for (int i=0; i<nRows; i++)
  {
    int r = i;
    if (isUpper) r = nRows - 1 - i;
    Vector<Scalar> rhs_r = rhs.getBlock(r);
    for (int j=0; j<i; j++)
    {
      int c = j;
      if (isUpper) c = nCols - 1 - j;
      if (op.getBlock(r,c).ptr().get() != 0)
      {
        rhs_r = rhs_r - op.getBlock(r,c) * soln.getBlock(c);
      }
    }

    SolverState<Scalar> state;
    Vector<Scalar> soln_r;
    if (oneSolverFitsAll)
    {
      state = solvers_[0].solve(op.getBlock(r,r), rhs_r, soln_r);
    }
    else
    {
      state = solvers_[r].solve(op.getBlock(r,r), rhs_r, soln_r);
    }
    if (nRows > 1) soln.setBlock(r, soln_r);
    else soln = soln_r;
    if (state.finalState() != SolveConverged)
    {
      return state;
    }
  }

  return SolverState<Scalar>(SolveConverged, "block solves converged",
    0, ScalarTraits<Scalar>::zero());
}
  
}

#endif
