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


#include "PlayaSimpleBacktracking.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaLinearCombinationImpl.hpp"


namespace Playa
{
using std::endl;

SimpleBacktracking::SimpleBacktracking(const ParameterList& params)
  : LineSearchBase(params)
{;}

std::string SimpleBacktracking::description() const
{
  std::ostringstream oss;
  oss << "SimpleBacktracking(maxSteps=" << maxSteps()
      << ", minStepSize=" << minStepSize() << ")";
  return oss.str();
}

LineSearchStatus SimpleBacktracking::search(const RCP<ObjectiveBase>& obj,
  const Vector<double>& x0,
  const double& f0,
  const Vector<double>& direction,
  const double& alphaMax,
  Vector<double>& xn, 
  Vector<double>& gradF,
  double& fVal) const
{
  Tabs tab(0);
  try
  {
    if (verb() > 0)
    {
      Out::root() << tab << "backtracking line search" << endl;
    }
    double alpha = alphaMax;
    if (verb() > 3)
    {
      Out::root() << tab << "line search initial vector " << endl;
      Out::os() << x0 << endl;
      Out::root() << tab << "line search direction " << endl;
      Out::os() << direction << endl;
      Out::root() << tab << "line search initial value " << f0 << endl;
    }

    for (int i = 0; i < maxSteps(); i++)
    {
      Tabs tab1;
      if (verb() > 1)
      {
        Out::root() << tab1 << "line search step " 
                    << i << " alpha=" << alpha << endl;
      }

      xn = x0 + alpha * direction;

      if (verb() > 3)
      {
        Out::root() << tab1 << "line search trial vector " << endl;
        Out::os() << xn << endl;
      }
          
      obj->evalGrad(xn, fVal, gradF);

      if (verb() > 1)
      {
        Out::root() << tab1 << "f=" << fVal << " f0=" << f0 << endl;
      }

      if (alpha < minStepSize()) return LS_StepTooSmall;
          
      if (fVal < f0)
      {
        if (verb() > 0)
        {
          Out::root() << tab1 << "Line search successful: steps = " << i << endl;
        }
        return LS_Success;
      }
      alpha = alpha/2.0;
    }
    return LS_ExceededMaxiters;
  }

  catch(std::exception& e)
  {
    Out::root() << "Exception detected in SimpleBacktracking::search(): " 
                << e.what() << endl;
    return LS_Crashed;
  }
}



}
