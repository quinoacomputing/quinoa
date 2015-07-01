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

#include "SundanceTransientStepProblem.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"

namespace Sundance
{

TransientStepProblem::TransientStepProblem(
    const NonlinearProblem& stepProb, 
    const Expr& tPrev, const Expr& uPrev, 
    const Expr& tNext, const Expr& uNext,
    const Expr& dt,
    int verbosity
    )
  : prob_(stepProb),
    tPrev_(tPrev),
    tNext_(tNext),
    uPrev_(uPrev),
    uNext_(uNext),
    dt_(dt),
    verb_(verbosity)
{
}

bool TransientStepProblem::step(
  double tCur, const Expr& uCur,
  double tNext, Expr uNext,
  const NonlinearSolver<double>& solver) const
{
  Tabs tab;
  
  PLAYA_MSG1(verb_, tab << "step from t=" << tCur << " to " << tNext);
  tPrev_.setParameterValue(tCur);
  tNext_.setParameterValue(tNext);
  dt_.setParameterValue(tNext-tCur);

  PLAYA_MSG2(verb_, tab << "updating uPrev");
  updateDiscreteFunction(uCur, uPrev_);

  PLAYA_MSG2(verb_, tab << "updating uNext");
  updateDiscreteFunction(uNext, uNext_);

  PLAYA_MSG2(verb_, tab << "Solving NLP");
  SolverState<double> state = prob_.solve(solver);

  PLAYA_MSG2(verb_, tab << "Updating solution");
  updateDiscreteFunction(uNext_, uNext);

  return state.finalState() == SolveConverged;
}

}
