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

#ifndef PLAYA_SOLVERSTATE_HPP
#define PLAYA_SOLVERSTATE_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Playa
{
using namespace Teuchos;

/** 
 * SolverStatusCode is an enum that encapsulates whether a solver succeeded or failed
 */
enum SolverStatusCode {SolveCrashed, SolveFailedToConverge, SolveConverged};


/**
 * SolverState provides information about the result of a linear or nonlinear solve
 */
template <class Scalar>
class SolverState
{
public:
  /** */
  SolverState(SolverStatusCode finalState, const std::string& msg, 
    int finalIters, const Scalar& finalResid)
    : finalState_(finalState),
      finalResid_(finalResid),
      finalIters_(finalIters),
      msg_(msg)
    {;}

  /** */
  SolverState() {;}

  /** */
  const Scalar& finalResid() const {return finalResid_;}

  /** */
  int finalIters() const {return finalIters_;}

  /** */
  const SolverStatusCode& finalState() const {return finalState_;}

  /** */
  const std::string& finalMsg() const {return msg_;}

  /** */
  std::string stateDescription() const 
    {
      switch (finalState_)
      {
        case SolveCrashed:
          return "Crashed";
        case SolveFailedToConverge:
          return "Failed to converge";
        case SolveConverged:
          return "Converged";
      }
      return "Crashed";
    }

private:
    
  SolverStatusCode finalState_;

  Scalar finalResid_;

  int finalIters_;

  std::string msg_;
};


template <class Scalar> inline
std::ostream& operator<<(std::ostream& os, 
  const Playa::SolverState<Scalar>& state)
{
  os << "Solver final state: " << state.stateDescription() << std::endl;
  os << "message: " << state.finalMsg() << std::endl;
  os << "iters taken: " << state.finalIters() << std::endl;
  os << "final residual: " << state.finalResid() << std::endl;
  return os;
}

}


#endif
