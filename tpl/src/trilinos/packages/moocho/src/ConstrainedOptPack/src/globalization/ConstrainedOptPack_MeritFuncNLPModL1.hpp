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

#ifndef MERIT_FUNC_NLP_MOD_L1_H
#define MERIT_FUNC_NLP_MOD_L1_H

#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "ConstrainedOptPack_MeritFuncNLPDirecDeriv.hpp"
#include "ConstrainedOptPack_MeritFuncPenaltyParams.hpp"

namespace ConstrainedOptPack {

/** \brief The modified L1 merit function using different penatly parameters for each constriant.
  *
  * phi(x) = f) + sum( mu(j) * abs(c(j)), j = 1,...,m )
  *
  * Dphi(x_k,d_k) = Gf_k' * d_k - sum( mu(j) * abs(c(j)), j = 1,...,m )
  *
  * Note that the definition of Dphi(x_k,d_k) assumes
  * that Gc_k'*d_k + c_k = 0.  In otherwords, d_k must
  * satisfiy the linearized equality constraints at
  * at x_k.
  *
  * Implicit copy constructor and assignment operators
  * are allowed.
  */
class MeritFuncNLPModL1
  : public MeritFuncNLP
  , public MeritFuncNLPDirecDeriv
  , public MeritFuncPenaltyParams
{
public:

  /// Initializes deriv() = 0 and mu() = 0
  MeritFuncNLPModL1();

  /** @name Overridden from MeritFuncNLP */
  //@{

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

  /** \brief If the value n passed to resize(n) does not
    * equal the size of the vector parameters then
    * an exception #MeritFuncNLP::InvalidInitialization#
    * will be thrown.
    */
  value_type calc_deriv(
    const Vector    &Gf_k
    ,const Vector   *c_k
    ,const Vector   *h_k
    ,const Vector   *hl
    ,const Vector   *hu
    ,const Vector   &d_k
    );
  
  //@}

  /** @name Overridden from MeritFuncPenaltyParams */
  //@{

  /** \brief . */
  void set_space_c( const VectorSpace::space_ptr_t& space_c );

  /** \brief . */
  VectorMutable& set_mu();

  /** \brief . */
  const Vector& get_mu() const;

  //@}

private:
  value_type                   deriv_;
  VectorSpace::vec_mut_ptr_t   mu_;

};	// end class MeritFuncNLPModL1

}	// end namespace ConstrainedOptPack

#endif	// MERIT_FUNC_NLP_MOD_L1_H
