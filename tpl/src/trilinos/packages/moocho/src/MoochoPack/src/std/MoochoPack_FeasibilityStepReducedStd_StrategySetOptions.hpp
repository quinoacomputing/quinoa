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

#ifndef FEASIBILITY_STEP_REDUCED_STD_STRATEGY_SET_OPTIONS_H
#define FEASIBILITY_STEP_REDUCED_STD_STRATEGY_SET_OPTIONS_H

#include "MoochoPack_FeasibilityStepReducedStd_Strategy.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for FeasibilityStepReducedStd_Strategy from an
  * OptionsFromStream object.
  *
  * The default options group name is IndepDirecWithBoundsStd.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group FeasibilityStepReducedStd_Strategy {
  *    qp_objective = OBJ_MIN_FULL_STEP;
  *    qp_objective = OBJ_MIN_NULL_SPACE_STEP;
  *    qp_objective = OBJ_RSQP;
  *    qp_testing   = QP_TEST_DEFAULT;
  *    qp_testing   = QP_TEST;
  *    qp_testing   = QP_NO_TEST;
  }
  \end{verbatim}
  */
class FeasibilityStepReducedStd_StrategySetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      FeasibilityStepReducedStd_Strategy >
{
public:

  /** \brief . */
  FeasibilityStepReducedStd_StrategySetOptions(
      FeasibilityStepReducedStd_Strategy* target = 0
    , const char opt_grp_name[] = "FeasibilityStepReducedStd" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class FeasibilityStepReducedStd_StrategySetOptions

}	// end namespace MoochoPack

#endif	// FEASIBILITY_STEP_REDUCED_STD_STRATEGY_SET_OPTIONS_H
