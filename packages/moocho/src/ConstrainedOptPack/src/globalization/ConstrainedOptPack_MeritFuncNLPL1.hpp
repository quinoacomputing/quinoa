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

#ifndef MERIT_FUNC_NLP_L1_H
#define MERIT_FUNC_NLP_L1_H

#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "ConstrainedOptPack_MeritFuncNLPDirecDeriv.hpp"
#include "ConstrainedOptPack_MeritFuncPenaltyParam.hpp"

namespace ConstrainedOptPack {

/** \brief The L1 merit function.
 *
 * phi(x) = f(x) + mu * norm(c(x),1)
 *
 * Dphi(x_k,d_k) = Gf_k' * d_k - mu * norm(c_k,1)
 *
 * Note that the definition of Dphi(x_k,d_k) assumes
 * that Gc_k'*d_k + c_k = 0.  In otherwords, d_k must
 * satisfiy the linearized equality constraints at
 * at x_k.
 *
 * ToDo: Add a term for general inequalities hl <= h <= hu
 *
 * Implicit copy constructor and assignment operators
 * are allowed.
 */
class MeritFuncNLPL1
  : public MeritFuncNLP
  , public MeritFuncNLPDirecDeriv
  , public MeritFuncPenaltyParam
{
public:

  /// Initializes deriv() = 0 and mu() = 0
  MeritFuncNLPL1();

  /** @name Overridden from MeritFuncNLP */
  //@{

  /** \brief . */
  MeritFuncNLP& operator=(const MeritFuncNLP&);
  /** \brief . */
  value_type value(
    value_type             f
    ,const Vector    *c
    ,const Vector    *h
    ,const Vector    *hl
    ,const Vector    *hu
    ) const;
  /** \brief . */
  value_type deriv() const;
  /** \brief . */
  void print_merit_func(
    std::ostream& out, const std::string& leading_str ) const;

  //@}

  /** @name Overridden from MeritFuncNLPDirecDeriv */
  //@{

  /** \brief . */
  value_type calc_deriv(
    const Vector    &Gf_k
    ,const Vector   *c_k
    ,const Vector   *h_k
    ,const Vector   *hl
    ,const Vector   *hu
    ,const Vector   &d_k
    );

  //@}

  /** @name Overridden from MeritFuncPenaltyParam */
  //@{

  /** \brief . */
  void mu(value_type mu);

  /** \brief . */
  value_type mu() const;

  //@}

private:
  value_type deriv_;
  value_type mu_;

};	// end class MeritFuncNLPL1

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_NLP_L1_H
