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

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_SET_OPTIONS_H
#define REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_SET_OPTIONS_H

#include "MoochoPack_ReducedHessianSecantUpdateBFGSProjected_Strategy.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace MoochoPack {

/** \brief Set options for ReducedHessianSecantUpdateBFGSProjected_Strategy
  * from a OptionsFromStream object.
  *
  * The options group is (with the default name):
  *
  \begin{verbatim}
  options_group ReducedHessianSecantUpdateBFGSProjected {
    act_set_frac_proj_start   = 0.8;   *** (+dbl)
    project_error_tol         = 1e-5;  *** (+dbl) [0.0, 1.0]
        super_basic_mult_drop_tol = 1e-5;  *** (+dbl)
  }
  \end{verbatim}
  *
  * \begin{description}
  *	\item[act_set_frac_proj_start] ToDo : Finish.
  *	\item[project_error_tol] ToDo : Finish.
  * \item[super_basic_mult_drop_tol] ToDo: Finish.
  *	\end{description}
  */
class ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      ReducedHessianSecantUpdateBFGSProjected_Strategy >
{
public:

  /** \brief . */
  ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions(
    ReducedHessianSecantUpdateBFGSProjected_Strategy* target = 0
    , const char opt_grp_name[] = "ReducedHessianSecantUpdatePBFGS" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions

}	// end namespace MoochoPack

#endif	// REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_SET_OPTIONS_H
