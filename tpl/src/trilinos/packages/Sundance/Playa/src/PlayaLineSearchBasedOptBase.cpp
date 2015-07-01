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


#include "PlayaLineSearchBasedOptBase.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaLineSearchBuilder.hpp"
#include "PlayaOptConvergenceTestBuilder.hpp"


namespace Playa
{
using std::endl;
using std::setw;

LineSearchBasedOptBase::LineSearchBasedOptBase(
  const ParameterList& params
  )
  : lineSearch_(),
    convTest_()
{
  const ParameterList& lsParams = params.sublist("Line Search");
  lineSearch_ = LineSearchBuilder::createLineSearch(lsParams);

  const ParameterList& ctParams = params.sublist("Convergence Test");
  convTest_ = OptConvergenceTestBuilder::createConvTest(ctParams);
}


OptState LineSearchBasedOptBase::run(const RCP<ObjectiveBase>& obj,
  const Vector<double>& xInit,
  const RCP<ConvergenceMonitor>& convMonitor) const
{
  Tabs tab0(0);
  PLAYA_MSG1(verb(), tab0 << "in " << description() << "::run()");

  RCP<DirectionGeneratorBase> dirGen = makeDirectionGenerator();

  VectorSpace<double> vs = xInit.space();

  Vector<double> xCur = xInit.copy();
  Vector<double> gradCur =vs.createMember();
  double fCur;

  obj->evalGrad(xCur, fCur, gradCur);

  OptState state(xCur, fCur, gradCur);

  Vector<double> xNew = vs.createMember();
  Vector<double> gradNew = vs.createMember();
  Vector<double> p = vs.createMember();
  double fNew;

  /* -------- Main optimization loop ------------- */
  while (state.status()==Opt_Continue)
  {
    Tabs tab1;
    PLAYA_MSG2(verb(), 
      tab1 << "--------------------------------------------------------------");
    PLAYA_MSG2(verb(), tab1 << description() << " iter " 
      << setw(5) << state.iter() 
      << " f=" << setw(14) << state.fCur() << " |grad|=" << setw(14)
      << state.gradCur().norm2());
      

    Tabs tab2;
    Tabs tab3;
    PLAYA_MSG4(verb(), tab2 << "current x " << endl << tab3 << state.xCur());

    obj->iterationCallback(state.xCur(), state.iter());

    if (convMonitor.ptr()!=null) 
    {
      convMonitor->addRecord(state.iter(), tuple(state.fCur()));
    }


    /* Find a direction for the line search */
    bool dirOK = dirGen->generateDirection(obj, state.xCur(), 
      state.gradCur(), state.fCur(), p);

    if (!dirOK)
    {
      state.setStatus(Opt_DirectionFailure);
      return state;
    }

    /* Perform the line search */
    LineSearchStatus info = lineSearch_->search(obj, 
      state.xCur(), state.fCur(), p,
      1.0, xNew, gradNew, fNew);

    if (info != LS_Success)
    {
      state.setStatus(Opt_LineSearchFailed);
      return state;
    }

    PLAYA_MSG4(verb(), tab2 << "new x " << endl << tab3 << xNew);

    /* Update the state to know about the new approximation to the min */
    state.update(xNew, gradNew, fNew);

    /* Check for convergence */
    OptStatus iterStatus = convTest_->test(state);
    state.setStatus(iterStatus);
  }

  /* -------- End main optimization loop ------------- */
  if (convMonitor.ptr()!=null) 
  {
    convMonitor->addRecord(state.iter(), tuple(state.fCur()));
  }

  return state;
}

void LineSearchBasedOptBase::print(std::ostream& os) const 
{
  Tabs tab(0);
  os << tab << description() << "[" << endl;
  Tabs tab1;
  os << tab1 << "Line search=" << endl;
  lineSearch_->print(os);
  os << tab1 << "Convergence test=" << endl;
  convTest_->print(os);
  os << tab << "]";
  
}


}
