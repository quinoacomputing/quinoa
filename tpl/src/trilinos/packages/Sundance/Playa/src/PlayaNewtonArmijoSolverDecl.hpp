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

#ifndef PLAYA_NEWTON_ARMIJO_SOLVER_DECL_HPP
#define PLAYA_NEWTON_ARMIJO_SOLVER_DECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaNonlinearSolverBase.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Playa
{
using namespace Teuchos;

/**
 * Playa implementation of Newton's method with Armijo line search.
 *
 * The solver's behavior is controlled by parameters in a ParameterList.
 * <ul>
 * <li> Scalar "Tau Relative" tolerance for relative error. Default value: 10 times machine epsilon
 * <li> Scalar "Tau Absolute" tolerance for absolute error. Default value 10 times machine epsilon
 * <li> Scalar "Alpha" constant in Armijo sufficient decrease condition. Default value 1.0e-4
 * <li> double "Step Reduction" factor by which to reduce step during line search. Default value: 0.5 
 * <li> int "Max Iterations" number of iterations to allow before failure. Default value 20.
 * <li> int "Max Backtracks" number of step reductions to allow before failure. Default value 20.
 * <li> int "Verbosity" amount of diagnostic output. Default value 0.
 * </ul>
 */
template <class Scalar>
class NewtonArmijoSolver : public NonlinearSolverBase<Scalar> 
{
public:
  /** */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** */
  NewtonArmijoSolver(const ParameterList& params, 
    const LinearSolver<Scalar>& linSolver);

  /** */
  virtual ~NewtonArmijoSolver(){;}

  /** */
  SolverState<Scalar> solve(const NonlinearOperator<Scalar>& F,
    Vector<Scalar>& soln) const ;

  /* */
  GET_RCP(NonlinearSolverBase<Scalar>);

private:
  LinearSolver<Scalar> linSolver_;
  ScalarMag tauR_;
  ScalarMag tauA_;
  ScalarMag alpha_;
  ScalarMag stepReduction_;
  int maxIters_;
  int maxLineSearch_;
  int verb_;
    
};

  
}

#endif
