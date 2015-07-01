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


#ifndef PLAYA_OPT_STATE_H
#define PLAYA_OPT_STATE_H


#include "PlayaVectorDecl.hpp"

namespace Playa
{

/** 
 * OptStatus provides diagnostic information on the current state
 * of an optimization run.
 */
enum OptStatus
{
  Opt_Continue,
  Opt_Converged,
  Opt_DirectionFailure,
  Opt_ExceededMaxiters,
  Opt_LineSearchFailed,
  Opt_Crashed
};

/** \relates OptStatus */
inline std::ostream& operator<<(std::ostream& os, const OptStatus& s)
{
  switch (s)
  {
    case Opt_Continue:
      os << "Opt_Continue"; break;
    case Opt_Converged:
      os << "Opt_Converged"; break;
    case Opt_DirectionFailure:
      os << "Opt_DirectionFailure"; break;
    case Opt_ExceededMaxiters:
      os << "Opt_ExceededMaxiters"; break;
    case Opt_LineSearchFailed:
      os << "Opt_LineSearchFailed"; break;
    default:
      os << "Opt_Crashed";
  }
  return os;
}


/** 
 * OptState encapsulates the current state of an optimization run, for
 * use in convergence testing.
 */
class OptState
{
public:
  /** */
  OptState(const Vector<double>& xCur,
    const double& fCur,
    const Vector<double>& gradCur);

  /** */
  OptStatus status() const {return status_;}

  /** */
  void setStatus(const OptStatus status) {status_ = status;} 

  /** Return the current iteration count */
  int iter() const {return iter_;}

  /** Return the current objective function value */
  double fCur() const {return fCur_;}

  /** Return the previous objective function value */
  double fPrev() const {return fPrev_;}

  /** Return the current evaluation point */
  Vector<double> xCur() const {return xCur_;}

  /** Return the previous evaluation point */
  Vector<double> xPrev() const {return xPrev_;}

  /** Return the current gradient */
  Vector<double> gradCur() const {return gradCur_;}

  /** Return the previous gradientx */
  Vector<double> gradPrev() const {return gradPrev_;}

  /** */
  void update(const Vector<double>& xNew, const Vector<double>& gradNew, 
    const double& fNew);

private:
  OptStatus status_;
  int iter_;

  Vector<double> xCur_;
  Vector<double> xPrev_;

  Vector<double> gradCur_;
  Vector<double> gradPrev_;

  double fCur_;
  double fPrev_;

};

}

#endif

