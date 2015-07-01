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


#ifndef PLAYA_LINEARSOLVERDECL_HPP
#define PLAYA_LINEARSOLVERDECL_HPP

#include "PlayaTabs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"
#include "Teuchos_TimeMonitor.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearSolverBaseImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

inline static Teuchos::Time& solveTimer() 
{
  static Teuchos::RCP<Teuchos::Time> rtn 
    = Teuchos::TimeMonitor::getNewTimer("linear solve"); 
  return *rtn;
}

namespace Playa
{
using namespace Teuchos;

  
/**
 * \brief User-level linear solver object
 */
template <class Scalar>
class LinearSolver : public Playa::Handle<LinearSolverBase<Scalar> >
{
public:
  /** */
  LinearSolver() : Playa::Handle<LinearSolverBase<Scalar> >() {;}
  /** */
  LinearSolver( Playa::Handleable<LinearSolverBase<Scalar> >* rawPtr) 
    : Playa::Handle<LinearSolverBase<Scalar> >(rawPtr) {;}
  /** */
  LinearSolver(const RCP<LinearSolverBase<Scalar> >& smartPtr)
    : Playa::Handle<LinearSolverBase<Scalar> >(smartPtr) {;}


  /** Change the convergence tolerance. Default does nothing. */
  void updateTolerance(const double& tol) {this->ptr()->updateTolerance(tol);}

  /** Set a user-defined preconditioner */
  void setUserPrec(const LinearOperator<Scalar>& op,
    const LinearSolver<Scalar>& pSolver) ;

  /** Set a user-defined preconditioner */
  void setUserPrec(const PreconditionerFactory<Scalar>& pf);


  /** */
  SolverState<Scalar> solve(const LinearOperator<Scalar>& op,
    const Vector<Scalar>& rhs,
    Vector<Scalar>& soln) const ;
    
    

  /** */
  const ParameterList& parameters() const ;

  /** */
  ParameterList& parameters() ;
};

  
template <class Scalar> inline 
SolverState<Scalar> LinearSolver<Scalar>
::solve(const LinearOperator<Scalar>& op,
  const Vector<Scalar>& rhs,
  Vector<Scalar>& soln) const
{
  Tabs tab;
  TEUCHOS_TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
    "null pointer in LinearSolver<Scalar>::solve()");

  TEUCHOS_TEST_FOR_EXCEPTION(rhs.ptr().get()==0, std::runtime_error,
    "null rhs pointer in LinearSolver<Scalar>::solve()");

  TEUCHOS_TEST_FOR_EXCEPTION(op.ptr().get()==0, std::runtime_error,
    "null op pointer in LinearSolver<Scalar>::solve()");

  TimeMonitor timer(solveTimer());

  PLAYA_MSG1(this->ptr()->verb() * (MPIComm::world().getRank()==0), 
    tab << "Solver(" << this->description() << ") starting solve");

  SolverState<Scalar> rtn = this->ptr()->solve(op, rhs, soln);

  PLAYA_MSG1(this->ptr()->verb() * (MPIComm::world().getRank()==0), 
    tab << "Solver(" << this->description() << ") done solve:");
  Tabs tab1;
  PLAYA_MSG2(this->ptr()->verb() * (MPIComm::world().getRank()==0), 
    tab << "state=" << rtn);

  return rtn;    
}

template <class Scalar> inline 
const ParameterList& LinearSolver<Scalar>::parameters() const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
    "null pointer in LinearSolver<Scalar>::parameters()");
  return this->ptr()->parameters();
}

template <class Scalar> inline 
ParameterList& LinearSolver<Scalar>::parameters() 
{
  TEUCHOS_TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
    "null pointer in LinearSolver<Scalar>::parameters()");
  return this->ptr()->parameters();
}

  

  

}

#endif
