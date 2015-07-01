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

#ifndef MERIT_FUNC_PENALTY_PARAM_UPDATE_WITH_MULTBASE_ADDED_STEP_SET_OPTIONS_H
#define MERIT_FUNC_PENALTY_PARAM_UPDATE_WITH_MULTBASE_ADDED_STEP_SET_OPTIONS_H

#include "MoochoPack_MeritFunc_PenaltyParamUpdate_AddedStep.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for MeritFunc_PenaltyParamUpdate_AddedStep from a
 * OptionsFromStream object.
 *
 * The options group is:
 *
 \verbatim
  options_group MeritFuncPenaltyParamUpdate {
      small_mu     = 1e-6;
      min_mu_ratio = 1e-8;
      mult_factor  = 1e-4;
      kkt_near_sol = 1.0;
  }
 \endverbatim
 *
 * <ul>
 *	<li>[small_mu] The smallest mu allows when away from the soltion.<br>
 *		Example: small_mu = 1e-6;
 *	<li>[min_mu_ratio] This bounds the smallest mu(i) as:<br>
 *		min(mu(i))/max(mu(i)) >= min_mu_ratio.<br>
 *		Example: min_mu_ratio = 1e-4;
 *	<li>[mult_factor] Multiplicative factor for mu(j) = (1.0+mult_factor) * abs(lambda(j)).<br>
 *		Example: mult_factor = 1e-4;
 *	<li>[kkt_near_sol] When the total kkt_error is below kkt_near_sol a safer
 *		penalty update will be used.<br>
 *		Example: kkt_near_sol = 1.0;
 * </ul>
 */
class MeritFunc_PenaltyParamUpdate_AddedStepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      MeritFunc_PenaltyParamUpdate_AddedStep >
{
public:

  /** \brief . */
  MeritFunc_PenaltyParamUpdate_AddedStepSetOptions(
    MeritFunc_PenaltyParamUpdate_AddedStep* target = 0 );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class MeritFunc_PenaltyParamUpdate_AddedStepSetOptions

}	// end namespace MoochoPack

#endif	// MERIT_FUNC_PENALTY_PARAM_UPDATE_WITH_MULTBASE_ADDED_STEP_SET_OPTIONS_H
