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

#ifndef PLAYA_PRECONDITIONER_HPP
#define PLAYA_PRECONDITIONER_HPP

#include "PlayaDefs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaPreconditionerBase.hpp"

namespace Playa
{
  /**
   * \brief Preconditioner stores left and/or right operators for
   * use in preconditioning.
   */
  template <class Scalar> 
  class Preconditioner : public Playa::Handle<PreconditionerBase<Scalar> >
  {
  public:
    /* Boilerplate ctors */
    HANDLE_CTORS(Preconditioner, PreconditionerBase<Scalar>);

    /** Change the value of a double parameter */
    void changeParameter(const std::string& name, const double& value);

    /** Change the value of an integer parameter */
    void changeParameter(const std::string& name, int value);

    
    
    /** Left preconditioner */
    LinearOperator<Scalar> left() const ;
    
    /** Right preconditioner */
    LinearOperator<Scalar> right() const ;
    
    /** return true if this preconditioner has both left and
     * right components. */
    bool isTwoSided() const {return hasLeft() && hasRight();}
    
    /** return true if this preconditioner has a nontrivial left component */
    bool hasLeft() const ;
    
    /** return true if this preconditioner has
     * a nontrivial right component */
    bool hasRight() const ;
    
    /** return true if this preconditioner has neither left nor
     * right operators defined */
    bool isIdentity() const {return !hasLeft() && !hasRight();}
  };

  

  template <class Scalar> inline 
  LinearOperator<Scalar> Preconditioner<Scalar>::left() const 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
                       "null pointer in Preconditioner<Scalar>::left()");
    return this->ptr()->left();
  }

  template <class Scalar> inline 
  LinearOperator<Scalar> Preconditioner<Scalar>::right() const 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
                       "null pointer in Preconditioner<Scalar>::right()");
    return this->ptr()->right();
  }

  template <class Scalar> inline
  bool Preconditioner<Scalar>::hasLeft() const 
  {
    return (this->ptr().get()!=0 && this->ptr()->hasLeft());
  }

  template <class Scalar> inline
  bool Preconditioner<Scalar>::hasRight() const 
  {
    return (this->ptr().get()!=0 && this->ptr()->hasRight());
  }

  
}

#endif
