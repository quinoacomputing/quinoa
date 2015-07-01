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

#ifndef PRE_EVAL_NEW_POINT_BARRIER_STEP_H
#define PRE_EVAL_NEW_POINT_BARRIER_STEP_H


#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"

#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Standard evaluation step class for extra parameters in primal/dual barrier method.
 *
 * This class calculates \c invXu, \c invXl \c invXu_m_invXl
 *
 */

class PreEvalNewPointBarrier_Step
  : public IterationPack::AlgorithmStep // doxygen needs full path
  {
  public:

    /** \brief relative fraction for initializing x within
     *   bounds.
     *   xl_sb = min(xl+relative_bound_push*(xu-xl),
     *               xl + absolute_bound_push)
     *   xu_sb = max(xu-relative_bound_push*(xu-xl),
     *               xu - absolute_bound_push)
     *   if (xl_sb > xu_sb) then
     *      x = (xl + (xu-xl)/2
     *   else if (x < xl_sb) then 
     *      x = xl_sb
     *   else if (x > xu_sb) then
     *      x = xu_sb
     */
    STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, relative_bound_push );

    /** \brief absolute measure for initializing x within
     *   bounds.
     *   xl_sb = min(xl+relative_bound_push*(xu-xl),
     *               xl + absolute_bound_push)
     *   xu_sb = max(xu-relative_bound_push*(xu-xl),
     *               xu - absolute_bound_push)
     *   if (xl_sb > xu_sb) then
     *      x = (xl + (xu-xl)/2
     *   else if (x < xl_sb) then 
     *      x = xl_sb
     *   else if (x > xu_sb) then
     *      x = xu_sb
     */
    STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, absolute_bound_push );

    /** @name Overridden from AlgorithmStep */
    //@{
    /** \brief . */
    bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
           , poss_type assoc_step_poss);
    
    
    void print_step( const IterationPack::Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
             , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
    //@}

    /** Constructor.
     */
    PreEvalNewPointBarrier_Step(
      const value_type relative_bound_push = 0.01,
      const value_type absolute_bound_push = 0.001
      );
    //@}

  }; // end class PreEvalNewPointBarrier_Step

class PreEvalNewPointBarrier_StepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode,
    public OptionsFromStreamPack::SetOptionsToTargetBase< PreEvalNewPointBarrier_Step >
  {
  public:
    PreEvalNewPointBarrier_StepSetOptions(
      PreEvalNewPointBarrier_Step* target = 0,
      const char opt_grp_name[] = "PreEvalNewPointBarrier" );

  protected:

    /// Overridden from SetOptionsFromStreamNode
    void setOption( int option_num, const std::string& option_value );
  
  };	// end class PreEvalNewPointBarrier_StepSetOptions


}  // end namespace MoochoPack

#endif // PRE_EVAL_NEW_POINT_BARRIER_STEP_H
