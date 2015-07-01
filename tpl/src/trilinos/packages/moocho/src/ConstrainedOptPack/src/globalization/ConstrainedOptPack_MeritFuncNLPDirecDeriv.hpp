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

#ifndef MERIT_FUNC_NLP_DIREC_DERIV_H
#define MERIT_FUNC_NLP_DIREC_DERIV_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

/** \brief This class provides a mix-in interface for allowing subclass merit
 * functions to compute the directional 1D derivative at a base point.
 *
 * The quantities Gf(xo) (gradient of f(xo))
 * c(xo), h(xo) and d are used by several
 * types of merit functions to calculate the derivative of:<br>
 * d(phi(x_k + alpha_k*d_k))/d(alpha_k) at alpha_k = 0.
 *
 * It is generally assumed that d satisfies Gc_k'*d_k + c_k = 0 otherwise the
 * merit function would need Gc_k to compute this directional derivative
 * properly.
 */
class MeritFuncNLPDirecDeriv {
public:

  /** \brief . */
  virtual ~MeritFuncNLPDirecDeriv() {}

  /** @name To be overridden by subclasses */
  //@{

  /** \brief Calculate d(phi(x_k + alpha_k*d_k))/d(alpha_k) at alpha_k = 0.
    *
    * The value is stored internally by the subclass are returned by its
    * deriv() member usually.  The value is also returned from this
    * function.
    *
    * If the sizes of the vectors input do not aggree then
    * #std::length_error# exception will be thrown.
    */
  virtual value_type calc_deriv(
    const Vector    &Gf_k
    ,const Vector   *c_k
    ,const Vector   *h_k
    ,const Vector   *hl
    ,const Vector   *hu
    ,const Vector   &d_k
    ) = 0;

  //@}

};	// end class MeritFuncNLPDirecDeriv

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_NLP_DIREC_DERIV_H
