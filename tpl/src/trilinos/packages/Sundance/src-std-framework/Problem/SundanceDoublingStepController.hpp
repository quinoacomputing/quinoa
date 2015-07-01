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

#ifndef SUNDANCE_DOUBLING_STEP_CONTROLLER_H
#define SUNDANCE_DOUBLING_STEP_CONTROLLER_H

#include "SundanceFieldWriter.hpp"
#include "SundanceFieldWriterFactory.hpp"
#include "SundanceTransientStepProblem.hpp"
#include "SundanceExprComparison.hpp"
#include "SundanceEventDetector.hpp"



namespace Sundance
{

/** */
class StepControlParameters
{
public:
  /** */
  StepControlParameters()
    {initDefaults();}

  double tStart_;

  double tStop_;

  double tau_;

  double initialStepsize_;

  double minStepsizeFactor_;

  double maxStepsizeFactor_;

  double stepsizeReductionSafetyFactor_;

  int maxSteps_;

  int verbosity_;

  int stepOrder_;

private:
  void initDefaults()
    {
      tStart_ = 0.0;
      tStop_ = 0.0;
      tau_ = 1.0e-6;
      initialStepsize_ = 0.01;
      minStepsizeFactor_ = 0.01;
      maxStepsizeFactor_ = 10.0;
      stepsizeReductionSafetyFactor_ = 0.9;
      maxSteps_ = 100000;
      verbosity_ = 0;
      stepOrder_ = 2;
    }
};

/** */
class StepHookBase
{
public:
  virtual void call(const double& tCur, const Expr& uCur) const = 0 ;
};

/** */
class OutputControlParameters
{
public:
  /** */
  OutputControlParameters(
    const FieldWriterFactory& wf,
    const string& filename,
    const double& writeInterval,
    int verb=0)
    : 
    writeInterval_(writeInterval),
    wf_(wf),
    filename_(filename),
    verbosity_(verb)
    {}
  
  double writeInterval_;
  FieldWriterFactory wf_;
  string filename_;
  int verbosity_;
};


/** */
class DoublingStepController
{
public:
  /** */
  DoublingStepController(
    const TransientStepProblem& prob,
    const NonlinearSolver<double>& solver,
    const StepControlParameters& stepControl,
    const OutputControlParameters& outputControl,
    const RCP<ExprComparisonBase>& compare)
    : prob_(prob), 
      stepControl_(stepControl), 
      outputControl_(outputControl),
      solver_(solver),
      compare_(compare),
      eventHandler_()
    {}

  /** */
  void setEventHandler(RCP<EventDetectorBase> e)
    {eventHandler_ = e;}

  /** */
  void setStepHook(RCP<StepHookBase> h)
    {stepHook_ = h;}
      

  /** */
  bool run() const ;

  /** */
  void write(int index, double t, const Expr& u) const ;

private:
  TransientStepProblem prob_;
  StepControlParameters stepControl_;
  OutputControlParameters outputControl_;
  NonlinearSolver<double> solver_;
  RCP<ExprComparisonBase> compare_;
  RCP<EventDetectorBase> eventHandler_;
  RCP<StepHookBase> stepHook_;
};


}


#endif
