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


#ifndef PLAYA_BASIC_LMBFGS_OPT_H
#define PLAYA_BASIC_LMBFGS_OPT_H


#include "PlayaLineSearchBasedOptBase.hpp"
#include <deque>

namespace Playa
{
/** 
 * Implements the L-BFGS algorithm
 *
 * @author Kevin Long
 */
class BasicLMBFGS : public LineSearchBasedOptBase
{
public:
  /** */
  BasicLMBFGS(const ParameterList& params);
        
  /** */
  std::string description() const {return "LM-BFGS";}

  /** */
  RCP<DirectionGeneratorBase> makeDirectionGenerator() const ; 

private:
  int memSize_;
};

/** */
class BasicLMBFGSDirection : public DirectionGeneratorBase
{
public:
  /** */
  BasicLMBFGSDirection(int memSize);

  /** */
  ~BasicLMBFGSDirection() {}

  /** */
  bool generateDirection(const RCP<ObjectiveBase>& obj,
    const Vector<double>& xCur,
    const Vector<double>& gradCur,
    const double& fCur,
    Vector<double>& p) ;

private:
  int memSize_;
  Vector<double> xPrev_;
  Vector<double> gradPrev_;
  std::deque<Vector<double> > sMem_;
  std::deque<Vector<double> > yMem_;
};

}

#endif
