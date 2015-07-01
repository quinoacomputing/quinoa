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

#ifndef INIT_FIN_DIFF_REDUCED_HESSIAN_STEP_SET_OPTIONS_H
#define INIT_FIN_DIFF_REDUCED_HESSIAN_STEP_SET_OPTIONS_H

#include "MoochoPack_InitFinDiffReducedHessian_Step.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for InitFinDiffReducedHessian_Step from an
  * OptionsFromStream object.
  *
  * The default options group name is InitFinDiffReducedHessian.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group InitFinDiffReducedHessian {
    initialization_method   = SCALE_DIAGONAL_ABS;
    max_cond                = 1e+1;
    min_diag                = 1e-8;
    step_scale              = 1e-1;
  }
  \end{verbatim}
  *
  * \begin{description}
  *	\item[initialization_method] Determines how the diagonal is initialized.
  *		from the finite difference taken.
  *		\begin{description}
  *		\item[SCALE_IDENTITY]      diag(i) = max( ||rGf_fd||inf , smallest_ele )
  *		\item[SCALE_DIAGONAL]      diag(i) = max( rGf_fd(i)     , smallest_ele )
  *		\item[SCALE_DIAGONAL_ABS]  diag(i) = max( abs(rGf_fd(i)), smallest_ele )
  *		\end{description}
  *		where: smallest_ele = max( ||rGf_fd||inf / max_cond , min_diag )
  *	\item[max_cond] The maximum condition of the initialized matrix.
  *		See initialization_method.\\
  *		Example: max_cond = 1e+1.
  *	\item[min_diag] The smallest absolute diagonal element.\\
  *		Example: min_diag = 1e-8.
  *	\item[step_scale] scales the step for the finite difference by
  *		#u = scale_step / ||Z*e||inf#.
  *		The finite difference is then taken as:\\
  *		#rGf_fd = ( Z_k * g(x_k + u * Z*e - rGf_k ) / u#\\
  *		Example: step_scale = 1.0.
  */
class InitFinDiffReducedHessian_StepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      InitFinDiffReducedHessian_Step >
{
public:

  /** \brief . */
  InitFinDiffReducedHessian_StepSetOptions(
      InitFinDiffReducedHessian_Step* target = 0
    , const char opt_grp_name[] = "InitFinDiffReducedHessian" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class InitFinDiffReducedHessian_StepSetOptions

}	// end namespace MoochoPack

#endif	// INIT_FIN_DIFF_REDUCED_HESSIAN_STEP_SET_OPTIONS_H
