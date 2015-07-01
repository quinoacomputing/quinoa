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


#include "PlayaDefaultOptConvergenceTest.hpp"
#include "PlayaOptState.hpp"
#include "PlayaTabs.hpp"
#include "PlayaLinearCombinationImpl.hpp"

namespace Playa
{
using std::endl;

DefaultOptConvergenceTest::DefaultOptConvergenceTest(
  const ParameterList& params)
  : OptConvergenceTestBase(params),
    minIters_(params.get<int>("Min Iterations")), 
    maxIters_(params.get<int>("Max Iterations")),
    requiredPasses_(params.get<int>("Num Required Passes")),
    objTol_(params.get<double>("Objective Tolerance")), 
    gradTol_(params.get<double>("Gradient Tolerance")),
    stepTol_(params.get<double>("Step Tolerance")),
    xTyp_(1.0),
    fTyp_(1.0)
{
  if (params.isParameter("Typical X Scale")) xTyp_ = params.get<double>("Typical X Scale");
  if (params.isParameter("Typical F Scale")) fTyp_ = params.get<double>("Typical F Scale");
}


OptStatus DefaultOptConvergenceTest::test(const OptState& state) const
{
  Tabs tab(0);
  int i = state.iter();
  int conv = 0;

  PLAYA_MSG1(verb(), tab << "DefaultOptConvergenceTest testing iter #" << i);
  Tabs tab1;

  if (i < std::max(1, minIters_)) 
  {
    PLAYA_MSG2(verb(), tab1 << "iter #" << i 
      << " below minimum, skipping test");
    return Opt_Continue;
  }

  /* Get problem scale */
  double fCur = state.fCur();
  double fPrev = state.fPrev();
  double fScale = std::max( std::fabs(fPrev), std::fabs(fTyp_) );
  double xScale = std::max( state.xCur().normInf(), std::fabs(xTyp_) );

  /* check function value change */
  double objConv = std::fabs(fCur - fPrev)/fScale;
  if (objConv <= objTol_) conv++;
  PLAYA_MSG2(verb(), tab1 << "obj test: " << objConv << " tol=" << objTol_);

  /* check gradient */
  double gradConv = state.gradCur().normInf() * xScale / fScale ;
  PLAYA_MSG2(verb(), tab1 << "|grad|: " << gradConv << " tol=" << gradTol_);
  if (gradConv <= gradTol_) conv++;

  /* compute |xPrev_k - xCur_k| / xScale */
  double stepConv = (state.xCur() - state.xPrev()).normInf()/xScale;
  if (stepConv <= stepTol_) conv++;
  PLAYA_MSG2(verb(), tab1 << "step test " << stepConv << " tol=" << stepTol_);
  
  PLAYA_MSG2(verb(), tab1 << conv << " of " << requiredPasses_ << " criteria "
    "passed");

  if (conv >= requiredPasses_) 
  {
    PLAYA_MSG2(verb(), tab1 << "convergence detected!");
    return Opt_Converged;
  }

  if (i >= maxIters_) 
  {
    PLAYA_MSG2(verb(), "iter #" << i << " above maxiters, giving up");
    return Opt_ExceededMaxiters;
  }

  PLAYA_MSG2(verb(), tab1 << "not yet converged");
  return Opt_Continue;
}

void DefaultOptConvergenceTest::print(std::ostream& os) const
{
  Tabs tab(0);
  os << tab << "DefaultOptConvergenceTest[" << endl;
  {
    Tabs tab1;
    os << tab1 << "min iterations " << minIters_ << endl;
    os << tab1 << "max iterations " << maxIters_ << endl;
    os << tab1 << "objective tol  " << objTol_ << endl;
    os << tab1 << "gradient tol   " << gradTol_ << endl;
    os << tab1 << "step tol       " << stepTol_ << endl;
  }
  os << tab << "]" << endl;
}


}
