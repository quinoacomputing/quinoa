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

#ifndef PLAYA_NONLINEAROPERATOR_HPP
#define PLAYA_NONLINEAROPERATOR_HPP

#include "PlayaDefs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaNonlinearOperatorBase.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Playa
{
using namespace Teuchos;

/** 
 * User-level nonlinear operator class
 */
template <class Scalar>
class NonlinearOperator : public Handle<NonlinearOperatorBase<Scalar> >
{
public:
  /* boilerplate ctors */
  HANDLE_CTORS(NonlinearOperator<Scalar>, NonlinearOperatorBase<Scalar>);

  /** */
  VectorSpace<Scalar> domain() const 
    {return this->ptr()->domain();}

  /** */
  VectorSpace<Scalar>  range() const 
    {return this->ptr()->range();}

  /** */
  void setEvalPt(const Vector<double>& evalPt) const 
    {
      this->ptr()->setEvalPt(evalPt);
    }
      
  /** */
  LinearOperator<Scalar> getJacobian() const 
    {
      return this->ptr()->getJacobian();
    }

  /** */
  Vector<double> getFunctionValue() const 
    {
      return this->ptr()->getFunctionValue();
    }

      

  /** */
  Vector<double> getInitialGuess() const 
    {
      return this->ptr()->getInitialGuess();
    }

  /** */
  Vector<double> currentEvalPt() const 
    {
      return this->ptr()->currentEvalPt();
    }

private:
};
}


#endif
