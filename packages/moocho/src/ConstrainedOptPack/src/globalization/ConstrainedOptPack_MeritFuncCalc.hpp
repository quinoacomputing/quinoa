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

#ifndef MERIT_FUNC_CALC_H
#define MERIT_FUNC_CALC_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** \brief Abstract iterface for n-D merit functions {abstract}.
  *
  * Used to compute the value of the merit at a point \a x (\c phi(x)) and
  * to retrieve the derivative (\c phi.deriv()) along some direction \a d from
  * some base point \a xo.
  */
class MeritFuncCalc  {
public:

  /** \brief . */
  virtual ~MeritFuncCalc() {}

  /** \brief Return the value of the merit function at \a x.
    */
  virtual value_type operator()(const Vector& x) const= 0;

  /// Calls value(d_k) on aggregate merit_func.
  virtual value_type deriv() const = 0;

  /// Print what this merit function is
  virtual void print_merit_func(
    std::ostream& out, const std::string& leading_str ) const = 0;

};	// end class MeritFuncCalc

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_CALC_H
