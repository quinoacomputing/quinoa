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

#ifndef LINE_SEARCH_2ND_ORDER_CORRECT_STEP_SET_OPTIONS_H
#define LINE_SEARCH_2ND_ORDER_CORRECT_STEP_SET_OPTIONS_H

#include "MoochoPack_LineSearch2ndOrderCorrect_Step.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for LineSearch2ndOrderCorrect_Step from an
  * OptionsFromStream object.
  *
  * The default options group name is LineSearch2ndOrderCorrect.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group LineSearch2ndOrderCorrect {
      newton_olevel           = PRINT_NOTHING;
      constr_norm_threshold   = 1e-3;
      constr_incr_ratio       = 5.0;
      after_k_iter            = 3;
      forced_constr_reduction = LESS_X_D;
      forced_reduct_ratio     = 1.0;
      max_step_ratio          = 0.7;
      max_newton_iter         = 3;
  }
  \end{verbatim}
  *
  * \begin{description}
  *	\item[newton_olevel] This is the output level for the internal newton iterations
  *		that compute the second order correction.  The value of this is usually
  *		determined by default to be compatable with the output level for the
  *		whole algorithm.  This option does not effect the rest of the output
  *		in the LinSearch2ndOrderCorrect step.
  *		The output level can be set specifically to the values of:
  *		\begin{description}
  *		\item[PRINT_USE_DEFAULT] Let it be determined by overall algorithm print level.
  *		\item[PRINT_NOTHING] No output about the newton iterations is performed.
  *		\item[PRINT_SUMMARY_INFO] A compact summary table is created.
  *		\item[PRINT_STEPS] Print more detailed info about the steps.
  *		\item[PRINT_VECTORS] Also print out calculated vectors.  Careful, this
  *			could produced a lot of output for large problems.
  *		\end{description}
  *	\item[constr_norm_threshold] See after_k_iter.
  *	\item[after_k_iter] When ||c_k||inf < |constr_norm_threshold| and
  *		k >= |after_k_iter| then the second order correction will be considered
  *		if full steps are not taken.  Having a dual test allows an initail
  *		nonoptimal ||c(x)|| = 0 to force a second order correction too early.
  *	\item[forced_constr_reduction] Determines how much the constraint error
  *		should be reduced for each SQP iteration.  The legal values are:
  *		\begin{description}
  *		\item[LESS_X_D] Just find a point where ||c(x_k+d_k+w)|| < ||c(x_k+d_k)||.
  *			The algorithm should be able to compute this point in one iteration
  *			unless it can't compute a descent direction for ||c(x)||.
  *		\item[LESS_X] Find a point where ||c(x_k+d_k+w)|| < ||c(x_k)||.
  *			This is a much more difficult test requirement.  If the algorithm
  *			can not find an acceptable point in the alotted iterations
  *			(see max_newton_iter) then a point ||c(x_k+d_k+w)|| < ||c(x_k+d_k)||
  *			will be accepted if possible.
  *		\end{description}
  *	\item[max_step_ratio] This limits the size of the correction term to be:
  *		w = min( 1, |max_step_ratio|/(||w||inf/||d_k||) ) * w.  This option
  *		keeps the x_k+d_k+w from getting too far away from x_k+d_k in case
  *		bad corrections are computed.  A value less than one is probably best.
  *	\item[max_newton_iter] Limits the number of newton iterations performed.
  *		It may be best to keep this a small number.
  *	\end{description}
  */
class LineSearch2ndOrderCorrect_StepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      LineSearch2ndOrderCorrect_Step >
{
public:

  /** \brief . */
  LineSearch2ndOrderCorrect_StepSetOptions(
      LineSearch2ndOrderCorrect_Step* target = 0
    , const char opt_grp_name[] = "LineSearch2ndOrderCorrect" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class LineSearch2ndOrderCorrect_StepSetOptions

}	// end namespace MoochoPack

#endif	// LINE_SEARCH_2ND_ORDER_CORRECT_STEP_SET_OPTIONS_H
