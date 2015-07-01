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

#ifndef DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_SET_OPTIONS_H
#define DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_SET_OPTIONS_H

#include "ConstrainedOptPack_DirectLineSearchArmQuad_Strategy.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace ConstrainedOptPack {

/** \brief Set options for DirectLineSearchArmQuad_Strategy from a
  * OptionsFromStream object.
  *
  * The default options group name is DirectLineSearchArmQuad.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group DirectLineSearchArmQuad {
    slope_frac = ?;
    min_frac_step = ?:
    max_frac_step = ?;
    max_ls_iter = ?;
    max_out_ls_iter = ?;
  }
  \end{verbatim}
  *
  * \begin{description}
  *	\item[slope_frac] Fraction of intial decent slope of merit
  *		function required by armijo test.  Mapps to eta in linesearch
  *		algorithm.  Must be in the range [0, 1].
  *		Example: slope_frac = 1.0e-4;
  *	\item[min_frac_step] Minimum fractional change in the step length
  *		per linesearch iteration.  Mapps to min_frac in the linesearch
  *		algorithm.  Must be in the range (0, 1].
  *		Example: min_frac_step = 0.1;
  *	\item[max_frac_step] Maximum fractional change in the step length
  *		per linesearch iteration.  Mapps to max_frac in the linesearch
  *		algorithm.  Must be in the range (0, 1].  Note that
  *		max_frac_step must be greater than min_frac_step.
  *		Example: min_frac_step = 0.5;
  *	\item[max_ls_iter] The maximum number of linesearch iterations
  *		to take before giving up and declaring a line search failure.
  *		Mapps to max_ls_iter in linesearch algorithm.
  *		Example: max_ls_iter = 20.
  *	\item[max_out_ls_iter] A flag to max out on line search iterations.
  *   Mostly just used for debugging, not very useful in general.
  *	\end{description}
  */
class DirectLineSearchArmQuad_StrategySetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      DirectLineSearchArmQuad_Strategy >
{
public:

  /** \brief . */
  DirectLineSearchArmQuad_StrategySetOptions(
      DirectLineSearchArmQuad_Strategy* qp_solver = 0
    , const char opt_grp_name[] = "DirectLineSearchArmQuad" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class DirectLineSearchArmQuad_StrategySetOptions

}	// end namespace ConstrainedOptPack

#endif	// DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_SET_OPTIONS_H
