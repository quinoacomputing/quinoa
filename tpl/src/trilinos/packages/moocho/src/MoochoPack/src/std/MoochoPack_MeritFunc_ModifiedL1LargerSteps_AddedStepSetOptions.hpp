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

#ifndef MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_SET_OPTIONS_H
#define MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_SET_OPTIONS_H

#include "MoochoPack_MeritFunc_ModifiedL1LargerSteps_AddedStep.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for MeritFunc_ModifiedL1LargerSteps_AddedStep from a
  * OptionsFromStream object.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group MeritFuncModifiedL1LargerSteps {
    after_k_iter                = 3;
    obj_increase_threshold      = 1e-3;
    max_pos_penalty_increase    = 1.0;
    pos_to_neg_penalty_increase = 1.0;
    incr_mult_factor            = 1e-4;
  }
  \end{verbatim}
  *
  * \begin{description}
  *	\item[after_k_iter] ToDo : Finish.
  *		Example: after_k_iter = 4;
  *	\item[obj_increase_threshold] ToDo : Finish.
  *		Example: obj_increase_threshold = 1e-4;
  *	\item[max_pos_penalty_increase] ToDo : Finish.
  *		Example: max_pos_penalty_increase = 1.0;
  *	\item[pos_to_neg_penalty_increase] ToDo : Finish.
  *		Example: pos_to_neg_penalty_increase = 1.0;
  *	\item[incr_mult_factor] ToDo : Finish.
  *		Example: incr_mult_factor = 1e-4;
  *	\end{description}
  */
class MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      MeritFunc_ModifiedL1LargerSteps_AddedStep >
{
public:

  /** \brief . */
  MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions(
    MeritFunc_ModifiedL1LargerSteps_AddedStep* target = 0 );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions

}	// end namespace MoochoPack

#endif	// MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_SET_OPTIONS_H
