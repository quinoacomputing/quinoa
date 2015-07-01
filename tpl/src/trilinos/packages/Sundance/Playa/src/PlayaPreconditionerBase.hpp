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

#ifndef PLAYA_PRECONDITIONERBASE_HPP
#define PLAYA_PRECONDITIONERBASE_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_ParameterList.hpp"


namespace Playa
{
  using namespace Teuchos;

  /**
   * Base class for preconditioners. A general preconditioner object
   * is split into a left preconditioner M1^-1 and a right
   * preconditioner M2^-1. To solve A x = b, we define the auxiliary
   * system M2^-1 y = x, and solve M1^-1 A M2^-1 y = M1^-1 b to obtain y.
   * Having y, we can quickly recover x by applying M2^-1 to y.
   *
   * The base class implements neither a left nor a right preconditioner.
   */
  template <class Scalar>
  class PreconditionerBase : public Playa::Handleable<PreconditionerBase<Scalar> >
  {
  public:
    /** empty ctor */
    PreconditionerBase() {;}

    /** virtual dtor */
    virtual ~PreconditionerBase(){;}

    
    /** */
    virtual LinearOperator<Scalar> left() const = 0 ;

    /** */
    virtual LinearOperator<Scalar> right() const = 0 ;

    /** return true if this preconditioner has a nontrivial left component */
    virtual bool hasLeft() const = 0 ;

    /** return true if this preconditioner has
     * a nontrivial right component */
    virtual bool hasRight() const = 0 ;

  private:
  };

}

#endif
