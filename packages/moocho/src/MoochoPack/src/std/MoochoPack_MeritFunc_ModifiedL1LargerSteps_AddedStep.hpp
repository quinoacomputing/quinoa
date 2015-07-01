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

#ifndef MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_H
#define MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_H

#include "rSQPAlgo_Step.h"
#include "ConstrainedOptPack_MeritFuncNLP.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief This function increases the penalty parameters of the modifed L1
  * merit function to allow for larger steps by taking advantage
  * of constraints that are reduced for a full step.
  */
class MeritFunc_ModifiedL1LargerSteps_AddedStep
  : public rSQPAlgo_Step
{
public:

  /// <<std comp>> members for merit_func
  STANDARD_COMPOSITION_MEMBERS(MeritFuncNLP,merit_func);

  /** \brief . */
  MeritFunc_ModifiedL1LargerSteps_AddedStep(
      const merit_func_ptr_t& merit_func
    , value_type	eta
    , int			after_k_iter				= 3
    , value_type	obj_increase_threshold		= 1e-4
    , value_type	max_pos_penalty_increase	= 1.0
    , value_type	pos_to_neg_penalty_increase	= 1.0
    , value_type	incr_mult_factor			= 1e-4 );

  /// eta.
  STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,eta);

  /// after_k_iter.
  STANDARD_MEMBER_COMPOSITION_MEMBERS(int,after_k_iter);

  /// obj_increase_threshold.
  STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,obj_increase_threshold);

  /// max_pos_penalty_increase
  STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,max_pos_penalty_increase);

  /// pos_to_neg_penalty_increase
  STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,pos_to_neg_penalty_increase);

  /// incr_mult_factor
  STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,incr_mult_factor);

  // ///////////////////////////////
  // Overridden from AlgorithmStep

  /** \brief . */
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);

  /** \brief . */
  void print_step( const Algorithm& algo, poss_type step_poss
    , IterationPack::EDoStepType type, poss_type assoc_step_poss
    , std::ostream& out, const std::string& leading_str ) const;

private:
  // not defined and not to be called.
  MeritFunc_ModifiedL1LargerSteps_AddedStep();
  MeritFunc_ModifiedL1LargerSteps_AddedStep(const MeritFunc_ModifiedL1LargerSteps_AddedStep&);
  MeritFunc_ModifiedL1LargerSteps_AddedStep& operator=(const MeritFunc_ModifiedL1LargerSteps_AddedStep&);
  
};	// end class MeritFunc_ModifiedL1LargerSteps_AddedStep

}	// end namespace MoochoPack 

#endif	// MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_H
