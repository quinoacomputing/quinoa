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

#ifndef CHECK_CONVERGENCE_STD_ADDED_STEP_SET_OPTIONS_H
#define CHECK_CONVERGENCE_STD_ADDED_STEP_SET_OPTIONS_H

#include "MoochoPack_CheckConvergenceStd_AddedStep.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for CheckConvergenceStd_AddedStep from an
  * OptionsFromStream object.
  *
  * The default options group name is CheckConvergenceStd.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group CheckConvergenceStd {
    scale_kkt_error_by   = SCALE_BY_ONE;
    scale_opt_error_by_Gf = true;
  }
  \end{verbatim}
  *
  * \begin{description}
  *	\item[scale_kkt_error_by] Determines if and how the optimality (opt_kkt_err)
  *		and feasiblity (feas_kkt_err)
  *		errors for the convergence check are scaled by for the unkowns x before
  *		comparing it to the set tolerances of opt_tol and feas_tol (see the
  *		class \Ref{CheckConvergenceStd_AddedStep} and its printed algorithm
  *		for more details).
  *		\begin{description}
  *		\item[SCALE_BY_ONE]			no scaling by x
  *		\item[SCALE_BY_NORM_2_X]    scale opt_kkt_err and feas_kkt_err by 1/||x||2
  *		\item[SCALE_BY_NORM_INF_X]  scale opt_kkt_err and feas_kkt_err by 1/||x||inf
  *		\end{description}
  *	\item[scale_opt_error_by_Gf] Determines if opt_kkt_err is scaled by
  *		||Gf_k||inf or not.
  *	\end{description}
  */
class CheckConvergenceStd_AddedStepSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      CheckConvergenceStd_AddedStep >
{
public:

  /** \brief . */
  CheckConvergenceStd_AddedStepSetOptions(
      CheckConvergenceStd_AddedStep* target = 0
    , const char opt_grp_name[] = "CheckConvergenceStd" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class CheckConvergenceStd_AddedStepSetOptions

}	// end namespace MoochoPack

#endif	// CHECK_CONVERGENCE_STD_ADDED_STEP_SET_OPTIONS_H
