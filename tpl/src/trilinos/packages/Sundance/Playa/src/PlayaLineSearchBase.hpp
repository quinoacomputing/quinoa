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


#ifndef PLAYA_LINESEARCHBASE_H
#define PLAYA_LINESEARCHBASE_H


#include "PlayaObjectiveBase.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Playa
{
  
using Teuchos::RCP;
using Teuchos::ParameterList;

/** */
enum LineSearchStatus {LS_Success, LS_ExceededMaxiters, LS_StepTooSmall, 
                       LS_Crashed};


/**
 * Base class for line search methods. 
 * \author Paul Boggs and Kevin Long
 */
class LineSearchBase : public ObjectWithVerbosity,
                       public Describable, 
                       public Printable
  
{
public:
  /** */
  LineSearchBase(const ParameterList& params); 


  /** */
  virtual LineSearchStatus search(const RCP<ObjectiveBase>& obj,
    const Vector<double>& x0,
    const double& f0,
    const Vector<double>& direction,
    const double& alphaMax,
    Vector<double>& xn, 
    Vector<double>& gradF,
    double& fVal) const = 0 ;
    
  /** */
  virtual double minStepSize() const {return minStepSize_;}

  /** */
  virtual int maxSteps() const {return maxSteps_;}

  /** */
  const ParameterList& params() const {return params_;}

private:
  ParameterList params_;

  int maxSteps_;

  double minStepSize_;
};
}

#endif
