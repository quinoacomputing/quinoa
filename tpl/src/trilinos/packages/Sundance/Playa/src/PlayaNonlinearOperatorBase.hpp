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

#ifndef PLAYA_NONLINEAROPERATORBASE_HPP
#define PLAYA_NONLINEAROPERATORBASE_HPP

#include "PlayaDefs.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaOut.hpp"
#include "PlayaObjectWithVerbosity.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaLinearCombinationDecl.hpp"

namespace Playa
{
using namespace Teuchos;

/** 
 * Base class for nonlinear operators
 */
template <class Scalar>
class NonlinearOperatorBase 
  : public Handleable<NonlinearOperatorBase<Scalar> >,
    public ObjectWithVerbosity
{
public:
  /** Empty ctor, for contexts in which we don't know the
   * domain and range spaces at the beginning of construction time */
  NonlinearOperatorBase() 
    : domain_(), range_(), 
      jacobianIsValid_(false),
      residualIsValid_(false),
      currentEvalPt_(),
      currentFunctionValue_(),
      currentJ_()
    {;}

  /** Construct a nonlinear operator with a domain and range */
  NonlinearOperatorBase(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range) 
    : domain_(domain.ptr()), range_(range.ptr()), 
      jacobianIsValid_(false),
      residualIsValid_(false),
      currentEvalPt_(),
      currentFunctionValue_(),
      currentJ_()
    {;}
                            
  /** Return the domain space */
  const RCP<const VectorSpaceBase<Scalar> >& domain() const 
    {return domain_;}

  /** Return the range space */
  const RCP<const VectorSpaceBase<Scalar> >& range() const 
    {return range_;}

  /** Set the evaluation point */
  void setEvalPt(const Vector<Scalar>& x) const 
    {
      if (this->verb() >= 1)
      {
        Out::os() << "NonlinearOperatorBase Setting new eval pt";
        if (this->verb() > 3)
        {
          Out::os() << " to " << std::endl ;
          x.print(Out::os());
        }
        Out::os() << std::endl;
      }
      jacobianIsValid_ = false;
      residualIsValid_ = false;

      currentEvalPt_.acceptCopyOf(x);
    }

  /** Get the current point at which the function is to be 
   * evaluated */
  const Vector<double>& currentEvalPt() const {return currentEvalPt_;}

  /** Return the Jacobian at the current evaluation point */
  LinearOperator<double> getJacobian() const 
    {
      if (this->verb() > 1)
      {
        Out::os() << "NonlinearOperatorBase getting Jacobian" << std::endl;
      }
      if (!jacobianIsValid_)
      {
        if (this->verb() > 3)
        {
          Out::os() << "...computing new J and F" << std::endl;
        }
        currentJ_ 
          = computeJacobianAndFunction(currentFunctionValue_);
        jacobianIsValid_ = true;
        residualIsValid_ = true;
      }
      else
      {
        if (this->verb() > 1)
        {
          Out::os() << "...reusing valid J" << std::endl;
        }
      }
      if (this->verb() > 3)
      {
        Out::os() << "J is " << std::endl;
        currentJ_.print(Out::os());
        Out::os() << std::endl;
      }
      return currentJ_;
    }

      

  /** Return the function value at the current evaluation point */
  Vector<double> getFunctionValue() const 
    {
      if (this->verb() > 1)
      {
        Out::os() << "NonlinearOperatorBase getting function value" << std::endl;
      }
      if (!residualIsValid_)
      {
        if (this->verb() > 1)
        {
          Out::os() << "...computing new F" << std::endl;
        }
        currentFunctionValue_ = computeFunctionValue();
        residualIsValid_ = true;
      }
      else
      {
        if (this->verb() > 1)
        {
          Out::os() << "...reusing valid F" << std::endl;
        }
      }

      if (this->verb() > 3)
      {
        Out::os() << "F is " << std::endl;
        currentFunctionValue_.print(Out::os());
        Out::os() << std::endl;
      }
      return currentFunctionValue_;
    }


  /** Return an initial guess appropriate to this problem */
  virtual Vector<double> getInitialGuess() const = 0 ;


protected:

  /** Compute the Jacobian at the current eval point */
  virtual LinearOperator<Scalar> computeJacobianAndFunction(Vector<double>& functionValue) const = 0 ;

  /** Compute the function value at the current eval point */
  virtual Vector<Scalar> computeFunctionValue() const 
    {
      computeJacobianAndFunction(currentFunctionValue_);
      return currentFunctionValue_;
    }

      
  /** Set the domain and range. This is protected so that solver
   * developers don't try to change the spaces on the fly */
  void setDomainAndRange(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range)
    {
      domain_ = domain.ptr();
      range_ = range.ptr();
    }


private:
  /** */
  RCP<const VectorSpaceBase<Scalar> > domain_;

  /** */
  RCP<const VectorSpaceBase<Scalar> > range_;

  /** */
  mutable bool jacobianIsValid_;

  /** */
  mutable bool residualIsValid_;

  /** */
  mutable Vector<double> currentEvalPt_;

  /** */
  mutable Vector<double> currentFunctionValue_;

  /** */
  mutable LinearOperator<double> currentJ_;
};



 
}


#endif
