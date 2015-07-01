// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef MERIT_FUNC_NLP_H
#define MERIT_FUNC_NLP_H

#include <iosfwd>

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** \brief Base class for all merit functions for NonLinear Programs (NLP) {abstract}.
  */
class MeritFuncNLP {
public:

  /** \brief . */
  class InvalidInitialization : public std::logic_error
  {public: InvalidInitialization(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /** \brief . */
  virtual ~MeritFuncNLP() {}

  /** @name To be overridden by subclasses */
  //@{

  /** \brief Assign the state of one Merit functions.
   *
   * The default implementation throws an <tt>std::logic_error</tt> exception
   * unless it is assignment to self.
   */
  virtual MeritFuncNLP& operator=(const MeritFuncNLP&);

  /** \brief Return the value of the merit function at f(x), c(x), h(x).
   * This interface requires the client to compute f(x)
   * c(x) and h(x) and pass it to this function to have
   * the value of phi(f,c,h,hl,hu) calculated.
   *
   * If the merit function has not been initialized properly
   * then a <tt>InvalidInitialization</tt> exception will be thrown.
   */
  virtual value_type value(
    value_type             f
    ,const Vector    *c
    ,const Vector    *h
    ,const Vector    *hl
    ,const Vector    *hu
    ) const = 0;

  /** \brief Return the value of the directional derivative of the 
   * merit function w.r.t. alpha at alpha = 0.  In other words
   * compute return d( phi(f(x),c(x),h(x)) ) / d(alpha_k) at alpha_k = 0
   * where x = x_k + alpha_k * d_k.
   *
   * If the merit function has not been initialized properly
   * then a <tt>InvalidInitialization</tt> exception will be thrown.
   */
  virtual value_type deriv() const = 0;

  /** \brief Print the merit funciton
   */
  virtual void print_merit_func(
    std::ostream& out, const std::string& leading_str ) const = 0;

  //@}

};	// end class MeritFuncNLP

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_NLP_H
