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


#include "PlayaBasicLMBFGS.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaLineSearchBuilder.hpp"
#include "PlayaOptConvergenceTestBuilder.hpp"

namespace Playa
{
using std::endl;

BasicLMBFGS::BasicLMBFGS(
  const ParameterList& params
  )
  : LineSearchBasedOptBase(params),
    memSize_(getParameter<int>(params, "Max Memory Size"))
{}


RCP<DirectionGeneratorBase> 
BasicLMBFGS::makeDirectionGenerator() const 
{
  return rcp(new BasicLMBFGSDirection(memSize_));
}


BasicLMBFGSDirection::BasicLMBFGSDirection(int memSize)
  : memSize_(memSize),
    xPrev_(),    
    gradPrev_(),
    sMem_(),
    yMem_()
{}

bool BasicLMBFGSDirection::generateDirection(
  const RCP<ObjectiveBase>& obj,
  const Vector<double>& xCur,
  const Vector<double>& gradCur,
  const double& fCur,
  Vector<double>& p)
{
  Vector<double> q = gradCur.copy();
  int numStored = sMem_.size();
  Array<double> alpha(numStored);
  Array<double> rho(numStored);

  for (int i=numStored-1; i>=0; i--)
  {
    rho[i] = 1.0/(sMem_[i]*yMem_[i]);
    alpha[i] = rho[i] * (sMem_[i]*q);
    q = q - alpha[i]*yMem_[i];
  }
  
  double gamma;
  if (numStored > 0)
  {
    int j = numStored-1;
    gamma = (sMem_[j]*yMem_[j])/(yMem_[j]*yMem_[j]);
  }
  else
  {
    gamma = obj->getInvHScale();
  }

  Vector<double> r = gamma*q;

  for (int i=0; i<numStored; i++)
  {
    double beta = rho[i]*(yMem_[i]*r);
    r = r + (alpha[i]-beta)*sMem_[i];
  }

  p = -1.0*r;

  if (xPrev_.ptr().get() != 0)
  {
    Vector<double> s = xCur - xPrev_;
    Vector<double> y = gradCur - gradPrev_;
    sMem_.push_back(s);
    yMem_.push_back(y);
    if ((int) sMem_.size() > memSize_)
    {
      sMem_.pop_front();
      yMem_.pop_front();
    }
  }
  
  xPrev_.acceptCopyOf(xCur);
  gradPrev_.acceptCopyOf(gradCur);

  return true;
}


}
