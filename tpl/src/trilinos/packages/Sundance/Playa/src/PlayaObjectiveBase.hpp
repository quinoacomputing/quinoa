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


#ifndef PLAYA_OBJECTIVE_BASE_H
#define PLAYA_OBJECTIVE_BASE_H

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaObjectWithVerbosity.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#endif




namespace Playa
{
/**
 * Base class for differentiable objective functions.
 * @author Paul Boggs and Kevin Long
 *
 */
class ObjectiveBase : public ObjectWithVerbosity,
                      public Describable
{
public:
  /** */
  ObjectiveBase(int verb=0) : ObjectWithVerbosity(verb), contextString_() {;}

  /** virtual dtor */
  virtual ~ObjectiveBase(){;}

  /** evaluate objective function and gradient */
  virtual void  evalGrad(const Vector<double>& x, double& f, 
    Vector<double>& grad) const = 0;

  /** evaluate objective function without gradient. */
  virtual void  eval(const Vector<double>& x, double& f) const = 0 ;

  /** return an initial guess for the design vector  */
  virtual Vector<double> getInit() const = 0;

  /** return an initial approximation to the scale for the 
   * inverse of the Hessian */
  virtual double getInvHScale() const {return 1.0;}

  /** User-overrideable hook for any callbacks to be done at the 
   * end of each iteration. Default is a no-op.  */
  virtual void iterationCallback(const Vector<double>& x, int iter) const {;}

  /** User-overrideable hook for any callbacks to be done at the 
   * end of an optimization loop. Default is a no-op.  */
  virtual void finalCallback(const Vector<double>& x) const {;}

  /** 
   * Set a string describing the context in which the function is being called.
   * This is intended to aid in reading diagnostic output. 
   */
  void setContextString(const string& str) const {contextString_ = str;}

  /** 
   * Set a string describing the context in which the function is being called.
   * This is intended to aid in reading diagnostic output. 
   */
  const string& contextString() const {return contextString_;}

  /** 
   * Return the number of evaluations
   */
  virtual int numFuncEvals() const {return -1;}

  /** 
   *
   */
  virtual string description() const {return "ObjectiveBase";}

  /** Debugging utility to check the gradient 
   * by comparing to a finite difference 
   * gradient calculation. This is will be expensive. */
  bool fdCheck(const Vector<double>& x, double tol, int verbosity=0) const ;

protected:
  /** */
  virtual double fdStep() const {return 1.0e-4;}


private:
  mutable string contextString_;
};

}

#endif
