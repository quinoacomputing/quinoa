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

#ifndef MERIT_FUNC_CALC_1D_H
#define MERIT_FUNC_CALC_1D_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** \brief Abstracts a 1D merit function {abstract}.
  *
  * This is the interface that line search algorithms use to compute
  * the value of the merit function at alpha (phi(alpha)) and
  * to retrieve the initial descent derivative of the merit function
  * (using \c deriv()).
  */
class MeritFuncCalc1D {
public:

  /** \brief . */
  virtual ~MeritFuncCalc1D() {}

  /// Return the value of the merit function at alpha.
  virtual value_type operator()( value_type alpha ) const = 0;

  /// Return the derivative of the merit function at alpha = 0
  virtual value_type deriv() const = 0;

  /// Print the particular merit function
  virtual void print_merit_func(
    std::ostream& out, const std::string& leading_str ) const = 0;

};	// end class MeritFuncCalc1D

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_CALC_1D_H
