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

#ifndef DAMPEN_CROSS_TERM_STD_STEP_H
#define DAMPEN_CROSS_TERM_STD_STEP_H

#include "rSQPAlgo_Step.h"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Compute a dampening term zeta_k for the cross term w_k such that
 * Gf'*Z*pz <= 0.
 *
 * This condition Gf'*Z*pz <= 0 is needed to ensure descent of many
 * merit functions.
 *
 * This implementation ensures that Gf'*Z*pz <= 0 only if
 * there will not be any active constraints when the reduced QP subproblem
 * is solved (nu_k = 0) and there are no undecomposed equality constraints
 * or if there they are linearly dependent (lambda_k(undecomp_con) = 0).
 * 
 * In particular this implementation computes zeta_k such that:
 * 
 * Gf'*Z*pz <= frac_descent * rGf'inv(B)*rGf
 * 
 * where: 0 < frac_descent < 1
 * 
 * To ensure strong descent (and hopefully deal with the cases where
 * nu_k != 0 and lambda_k(undecomp_con) != 0) the parameter frac_descent
 * is set to frac_descent = 0.9 by default.
 * 
 * The basis derivation goes like this:
 * 
 * ToDo: Finish documentation!
 */
class DampenCrossTermStd_Step : public rSQPAlgo_Step {
public:

  /// «std comp» members for frac_descent
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, frac_descent );

  /** \brief . */
  DampenCrossTermStd_Step(const value_type& frac_descent = 0.9);

  // ////////////////////
  // Overridden

  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);

  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;

};	// end class DampenCrossTermStd_Step

}	// end namespace MoochoPack 

#endif	// DAMPEN_CROSS_TERM_STD_STEP_H
