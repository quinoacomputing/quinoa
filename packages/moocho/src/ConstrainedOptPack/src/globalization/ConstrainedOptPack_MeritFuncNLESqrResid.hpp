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

#ifndef MERIT_FUNC_NLE_SQR_RESID_H
#define MERIT_FUNC_NLE_SQR_RESID_H

#include "ConstrainedOptPack_MeritFuncNLE.hpp"

namespace ConstrainedOptPack {

/** \brief A merit function for the square of the constriant values.
  *
  * phi(x) = 1/2 * c(x)'*c(x)
  *
  * Dphi(x_k,d_k) = - c(x)'*c(x)
  *
  * Note that the definition of Dphi(x_k,d_k) assumes
  * that Gc_k'*d_k + c_k = 0.  In otherwords, d_k must
  * satisfiy the linearized equality constraints at
  * at x_k.
  *
  * Implicit copy constructor and assignment operators
  * are allowed.
  */
class MeritFuncNLESqrResid : public MeritFuncNLE {
public:

  /// Initializes deriv() = 0
  MeritFuncNLESqrResid();

  /** \brief . */
  value_type calc_deriv( const Vector& c_k );

  // ////////////////////////////////
  // Overridden from MeritFuncNLE

  /** \brief . */
  value_type value(const Vector& c) const;

  /** \brief . */
  value_type deriv() const;

  /** \brief . */
  void print_merit_func(std::ostream& out
    , const std::string& leading_str) const;

private:
  value_type deriv_;

};	// end class MeritFuncNLESqrResid

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_NLE_SQR_RESID_H
