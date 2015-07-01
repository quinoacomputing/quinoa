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

#ifndef PLAYA_KRYLOVSOLVER_HPP
#define PLAYA_KRYLOVSOLVER_HPP

#include "PlayaDefs.hpp"
#include "PlayaIterativeSolver.hpp"
#include "PlayaPreconditionerFactory.hpp"
#include "PlayaILUKPreconditionerFactory.hpp"
#include "PlayaSimpleComposedOpDecl.hpp"

namespace Playa
{
using namespace Teuchos;

/**
 *
 */
template <class Scalar>
class KrylovSolver : public IterativeSolver<Scalar>
{
public:
  /** */
  KrylovSolver(const ParameterList& params);
  /** */
  KrylovSolver(const ParameterList& params,
    const PreconditionerFactory<Scalar>& precond);

  /** */
  virtual ~KrylovSolver(){;}

  /** */
  virtual SolverState<Scalar> solve(const LinearOperator<Scalar>& op,
    const Vector<Scalar>& rhs,
    Vector<Scalar>& soln) const ;
protected:
  virtual SolverState<Scalar> solveUnprec(const LinearOperator<Scalar>& op,
    const Vector<Scalar>& rhs,
    Vector<Scalar>& soln) const = 0 ;

  const PreconditionerFactory<Scalar>& precond() const {return precond_;}

private:
  PreconditionerFactory<Scalar> precond_;
};

  
template <class Scalar> inline
KrylovSolver<Scalar>::KrylovSolver(const ParameterList& params)
  : IterativeSolver<Scalar>(params), precond_()
{
  if (!params.isParameter("Precond")) return;

  const std::string& precondType = params.template get<string>("Precond");

  if (precondType=="ILUK")
  {
    precond_ = new ILUKPreconditionerFactory<Scalar>(params);
  }
}

template <class Scalar> inline
KrylovSolver<Scalar>::KrylovSolver(const ParameterList& params,
  const PreconditionerFactory<Scalar>& precond)
  : IterativeSolver<Scalar>(params), precond_(precond)
{
  TEUCHOS_TEST_FOR_EXCEPTION(params.isParameter("Precond"), std::runtime_error,
    "ambiguous preconditioner specification in "
    "KrylovSolver ctor: parameters specify "
    << params.template get<string>("Precond") 
    << " but preconditioner argument is " 
    << precond);
}

template <class Scalar> inline
SolverState<Scalar> KrylovSolver<Scalar>
::solve(const LinearOperator<Scalar>& op,
  const Vector<Scalar>& rhs,
  Vector<Scalar>& soln) const
{
  if (precond_.ptr().get()==0) 
  {
    return solveUnprec(op, rhs, soln);
  }


  Preconditioner<Scalar> p = precond_.createPreconditioner(op);
    
  if (!p.hasRight())
  {
    LinearOperator<Scalar> A = p.left()*op;
    Vector<Scalar> newRHS = rhs.space().createMember();
    p.left().apply(rhs, newRHS);
    return solveUnprec(A, newRHS, soln);
  }
  else if (!p.hasLeft())
  {
    LinearOperator<Scalar> A = op * p.right();
    Vector<Scalar> intermediateSoln;
    SolverState<Scalar> rtn 
      = solveUnprec(A, rhs, intermediateSoln);
    if (rtn.finalState()==SolveConverged) 
    {
      p.right().apply(intermediateSoln, soln);
    }
    return rtn;
  }
  else
  {
    LinearOperator<Scalar> A = p.left() * op * p.right();
    Vector<Scalar> newRHS;
    p.left().apply(rhs, newRHS);
    Vector<Scalar> intermediateSoln;
    SolverState<Scalar> rtn 
      = solveUnprec(A, newRHS, intermediateSoln);
    if (rtn.finalState()==SolveConverged) 
    {
      p.right().apply(intermediateSoln, soln);
    }
    return rtn;
  }
}
  
}

#endif

