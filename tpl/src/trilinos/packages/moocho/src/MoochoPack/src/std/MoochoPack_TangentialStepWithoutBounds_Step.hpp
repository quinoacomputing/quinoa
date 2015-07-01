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

#ifndef NULL_SPACE_STEP_WITHOUT_BOUNDS_STEP_H
#define NULL_SPACE_STEP_WITHOUT_BOUNDS_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Solves the unconstrained QP subproblem: <tt>min  qp_grad' * pz + (1/2) * pz' * rHL * pz</tt>.
  *
  * The solution to this system is just:<br>
  * <tt>pz = inv(rHL) *qp_grad</tt>.
  *
  * If use_qp_correc is false then:<br>
  *   <tt>qp_grad = rGf</tt>
  * else<br>
  *   <tt>qp_grad = rGf + zeta * ZtHLYpy.<br>
  *
  * Then <tt>Zpz = Z * pz</tt>
  */
class TangentialStepWithoutBounds_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /** \brief Set the maximum size for ||pz|| dampening.
   *
   * A value of <tt>max_pz_norm <= 0.0</tt> means not to dampen pz!
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_pz_norm );

  /** \brief Set the number of iterations to dampen pz for.
   *
   * A value of <tt>num_pz_damp_iters <= 0</tt> means not to dampen pz for any
   * iterations!
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, num_pz_damp_iters );

  /** \brief . */
  TangentialStepWithoutBounds_Step();

  /** @name Overridden from AlgorithmStep */
  //@{
  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);
  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
  //@}

};	// end class TangentialStepWithoutBounds_Step

}	// end namespace MoochoPack 

#endif	// NULL_SPACE_STEP_WITHOUT_BOUNDS_STEP_H
