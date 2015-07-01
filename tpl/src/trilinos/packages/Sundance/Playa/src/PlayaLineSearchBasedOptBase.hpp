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


#ifndef PLAYA_LINE_SEARCH_BASED_OPT_BASE_H
#define PLAYA_LINE_SEARCH_BASED_OPT_BASE_H

#include "PlayaUnconstrainedOptimizerBase.hpp"
#include "PlayaLineSearchBase.hpp"
#include "PlayaOptConvergenceTestBase.hpp"
#include "PlayaObjectiveBase.hpp"

namespace Playa
{

/** 
 * 
 */
class DirectionGeneratorBase : public ObjectWithVerbosity
{
public:
  /** */
  DirectionGeneratorBase() {}

  /** */
  ~DirectionGeneratorBase() {}

  /** */
  virtual bool generateDirection(const RCP<ObjectiveBase>& obj,
    const Vector<double>& xCur,
    const Vector<double>& gradCur,
    const double& fCur,
    Vector<double>& p) = 0 ;
};


/** 
 * Base class for optimizers based on line search methods.
 *
 * @author Kevin Long
 */
class LineSearchBasedOptBase : public UnconstrainedOptimizerBase 
{
public:
  /** */
  LineSearchBasedOptBase(const ParameterList& params);
  /** */
  virtual ~LineSearchBasedOptBase(){}
        
  /** Main method to apply the algorithm starting with x and
      returning the result in x */
  OptState run(const RCP<ObjectiveBase>& obj,
    const Vector<double>& xInit,
    const RCP<ConvergenceMonitor>& convMonitor = null) const ;

  /** */
  virtual RCP<DirectionGeneratorBase> makeDirectionGenerator() const = 0 ; 

  /** */
  virtual void print(std::ostream& os) const ;

protected:
  /** */
  const RCP<LineSearchBase>& lineSearch() const {return lineSearch_;}

  /** */
  const RCP<OptConvergenceTestBase>& convTest() const {return convTest_;}

private:
  RCP<LineSearchBase> lineSearch_;
  RCP<OptConvergenceTestBase> convTest_;
  
};


}

#endif
