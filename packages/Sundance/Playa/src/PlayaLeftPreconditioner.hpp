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

#ifndef PLAYA_GENERICLEFTPRECONDITIONER_HPP
#define PLAYA_GENERICLEFTPRECONDITIONER_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_ParameterList.hpp"


namespace Playa
{
  using namespace Teuchos;

  /**
   * A one-size-fits-most left preconditioner that can be constructed by
   * accepting an operator for the left op of the preconditioner. 
   */
  template <class Scalar>
  class GenericLeftPreconditioner : public PreconditionerBase<Scalar>
  {
  public:
    /** construct with  */
    GenericLeftPreconditioner(const LinearOperator<Scalar>& left) 
    : PreconditionerBase<Scalar>(), left_(left) {;}

    /** virtual dtor */
    virtual ~GenericLeftPreconditioner(){;}

    
    /** Return the left operator */
    virtual LinearOperator<Scalar> left() const {return left_;}

    /** A call to right() results in an error for a left precond. */
    virtual LinearOperator<Scalar> right() const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, logic_error, "right() called for a "
                         "preconditioner known to be a left precond");
      return LinearOperator<Scalar>();
    }

    /** return true because 
     * this preconditioner has a nontrivial left component. */
    virtual bool hasLeft() const {return true;}

    /** return false, because this preconditioner has
     * no nontrivial right component */
    virtual bool hasRight() const {return false;}

    /* Handleable boilerplate */
    GET_RCP(PreconditionerBase<Scalar>);

  private:
    LinearOperator<Scalar> left_;
  };

}
